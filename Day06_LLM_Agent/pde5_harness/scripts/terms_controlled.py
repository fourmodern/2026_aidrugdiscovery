#!/usr/bin/env python3
"""골격 통제 세트의 스코어 항 분석 — 어떤 항이 실제로 기여하는가.

n=163 은 이전 n=30 보다 항 5개를 적합하기에 훨씬 낫다. 그래도 다음을 전부 함께 낸다.
  - LOO 교차검증 Q2 (학습 R2 만으로 성능을 주장하지 않는다)
  - 라벨섞기 귀무분포 (우연 수준을 명시한다)
  - 표준화 계수 (비표준화 계수는 단위가 달라 서로 비교 불가)
  - 순열 중요도 (그 항만 섞었을 때 Q2 하락)
  - VIF (항끼리 겹치면 개별 계수를 해석할 수 없다)
  - **Tanimoto 를 넣은 모델과의 비교** (항이 역가를 보는지 골격을 보는지)
"""
from __future__ import annotations
import argparse, json, math, random, re, subprocess, sys, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"
WORK2 = SR / "structures" / "work_controlled"
SMINA = os.environ.get("SMINA", "/home/hjpark/vina/smina")
TERMS = ["gauss(o=0,_w=0.5,_c=8)", "gauss(o=3,_w=2,_c=8)", "repulsion(o=0,_c=8)",
         "hydrophobic(g=0.5,_b=1.5,_c=8)", "non_dir_h_bond(g=-0.7,_b=0,_c=8)"]
SHORT = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]


def extract(job):
    """한 화합물의 점수 1위 자세에서 항 값을 뽑는다."""
    cid, mode = job
    pose = WORK2 / f"{cid}_poses.sdf"
    if not pose.exists():
        return cid, None
    # 1위 모드만 잘라낸다
    txt = pose.read_text().split("$$$$\n")
    if len(txt) < mode or not txt[mode - 1].strip():
        return cid, None
    single = WORK2 / f"{cid}_m{mode}.sdf"
    single.write_text(txt[mode - 1] + "$$$$\n")
    cmd = [SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(single),
           "--score_only", "--custom_scoring", str(WORK2 / "terms.txt"), "--cpu", "1"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return cid, None
    vals = {}
    for line in p.stdout.splitlines():
        m = re.match(r"^##\s+(\S+)\s+(.*)$", line.strip())
        if m and m.group(1) not in ("Name",):
            nums = m.group(2).split()
            if len(nums) >= len(TERMS):
                try:
                    vals = dict(zip(SHORT, [float(x) for x in nums[:len(TERMS)]]))
                except ValueError:
                    pass
    if not vals:
        # 헤더/값 두 줄 형식 대응
        lines = [l for l in p.stdout.splitlines() if l.startswith("##")]
        for l in lines:
            parts = l.replace("##", "").split()
            nums = [x for x in parts if re.match(r"^-?\d+\.?\d*(e-?\d+)?$", x)]
            if len(nums) >= len(TERMS):
                vals = dict(zip(SHORT, [float(x) for x in nums[-len(TERMS):]]))
                break
    return cid, (vals or None)


def ols(X, y):
    k = len(X[0]); n = len(y)
    A = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] +
         [sum(X[i][a] * y[i] for i in range(n))] for a in range(k)]
    for c in range(k):
        piv = max(range(c, k), key=lambda r: abs(A[r][c]))
        if abs(A[piv][c]) < 1e-12: return None
        A[c], A[piv] = A[piv], A[c]
        for r in range(k):
            if r == c: continue
            f = A[r][c] / A[c][c]
            for j in range(c, k + 1): A[r][j] -= f * A[c][j]
    return [A[i][k] / A[i][i] for i in range(k)]


def r2(y, pred):
    m = sum(y) / len(y); ss = sum((a - m) ** 2 for a in y)
    return 1 - sum((a - b) ** 2 for a, b in zip(y, pred)) / ss if ss else float("nan")


def q2_loo(X, y):
    n = len(y); pr = []
    for i in range(n):
        b = ols([X[j] for j in range(n) if j != i], [y[j] for j in range(n) if j != i])
        if b is None: return float("nan")
        pr.append(sum(c * v for c, v in zip(b, X[i])))
    return r2(y, pr)


def zs(v):
    n = len(v); m = sum(v) / n
    sd = (sum((a - m) ** 2 for a in v) / (n - 1)) ** 0.5
    return [(a - m) / sd for a in v] if sd else [0.0] * n


def pearson(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    nu = sum((a - mx) * (b - my) for a, b in zip(x, y))
    de = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
    return nu / de if de else float("nan")


def fit(X, y, null_iters, seed):
    b = ols(X, y)
    if b is None: return None
    fitv = [sum(c * v for c, v in zip(b, x)) for x in X]
    rng = random.Random(seed); null = []
    for _ in range(null_iters):
        ys = y[:]; rng.shuffle(ys)
        bb = ols(X, ys)
        if bb: null.append(r2(ys, [sum(c * v for c, v in zip(bb, x)) for x in X]))
    null.sort()
    return {"coefficients": [round(v, 4) for v in b], "n_params": len(b),
            "R2_fit": round(r2(y, fitv), 3), "Q2_loo": round(q2_loo(X, y), 3),
            "null_R2_median": round(null[len(null) // 2], 3) if null else None,
            "null_R2_p95": round(null[int(len(null) * .95)], 3) if null else None,
            "null_iters": len(null)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dock", default=str(SR / "docking_controlled.json"))
    ap.add_argument("--out", default=str(SR / "terms_controlled.json"))
    ap.add_argument("--null-iters", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=32)
    a = ap.parse_args()

    WORK2.mkdir(parents=True, exist_ok=True)
    (WORK2 / "terms.txt").write_text("\n".join(f"1.0 {t}" for t in TERMS) + "\n")
    D = json.loads(Path(a.dock).read_text())["result"]
    rows = D["rows"]
    jobs = [(r["chembl_id"], 1) for r in rows]      # 점수 1위 자세
    print(f"항 추출 {len(jobs)}건…", file=sys.stderr)
    got = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(extract, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            cid, v = f.result()
            if v: got[cid] = v
            if i % 40 == 0: print(f"  {i}/{len(jobs)}", file=sys.stderr)
    print(f"  성공 {len(got)}/{len(jobs)}", file=sys.stderr)

    use = [r for r in rows if r["chembl_id"] in got]
    for r in use:
        r.update(got[r["chembl_id"]])
    n = len(use)
    y = [r["pchembl_value"] for r in use]
    tan = [r["tanimoto_to_sildenafil"] for r in use]
    sc = [r["top_pose_score"] for r in use]
    z = {t: zs([r[t] for r in use]) for t in SHORT}
    yz = zs(y)

    models = {
        "single_vina_score": fit([[1.0, v] for v in zs(sc)], yz, a.null_iters, 42),
        "five_terms": fit([[1.0] + [z[t][i] for t in SHORT] for i in range(n)],
                          yz, a.null_iters, 42),
        "tanimoto_only": fit([[1.0, v] for v in zs(tan)], yz, a.null_iters, 42),
        "five_terms_plus_tanimoto": fit(
            [[1.0] + [z[t][i] for t in SHORT] + [zs(tan)[i]] for i in range(n)],
            yz, a.null_iters, 42),
    }
    full = models["five_terms"]
    fq = full["Q2_loo"]
    Xfull = [[1.0] + [z[t][i] for t in SHORT] for i in range(n)]

    rng = random.Random(42); perm = {}
    for t in SHORT:
        drops = []
        for _ in range(30):
            col = z[t][:]; rng.shuffle(col)
            Xp = [[1.0] + [(col[i] if s == t else z[s][i]) for s in SHORT] for i in range(n)]
            qq = q2_loo(Xp, yz)
            if qq == qq: drops.append(fq - qq)
        perm[t] = {"mean_Q2_drop": round(sum(drops) / len(drops), 4) if drops else None,
                   "n_shuffles": len(drops)}

    vif = {}
    for t in SHORT:
        others = [s for s in SHORT if s != t]
        Xo = [[1.0] + [z[s][i] for s in others] for i in range(n)]
        b = ols(Xo, z[t])
        rr = r2(z[t], [sum(c * v for c, v in zip(b, x)) for x in Xo]) if b else float("nan")
        vif[t] = round(1 / (1 - rr), 1) if rr == rr and rr < .999 else None

    uni = {t: {"pearson_r": round(pearson([r[t] for r in use], y), 3),
               "pearson_r_vs_tanimoto": round(pearson([r[t] for r in use], tan), 3)}
           for t in SHORT}
    std = dict(zip(SHORT, [round(v, 3) for v in full["coefficients"][1:]]))

    out = {"n": n, "terms": TERMS, "short": SHORT, "models": models,
           "standardized_coefficients": std, "univariate": uni,
           "permutation_importance": perm, "vif": vif,
           "term_correlation": {t: {s: round(pearson([r[t] for r in use],
                                                     [r[s] for r in use]), 2)
                                    for s in SHORT} for t in SHORT},
           "rows": [{k: r[k] for k in ("chembl_id", "pchembl_value", "top_pose_score",
                                       "tanimoto_to_sildenafil", "potency_bin",
                                       "similarity_bin", *SHORT)} for r in use]}
    checks = [("항 추출 성공률 >= 95%", len(got) >= .95 * len(jobs)),
              ("네 모델 모두 적합", all(v is not None for v in models.values())),
              ("전 항 VIF 산출", all(v is not None for v in vif.values())),
              ("라벨섞기 >= 1000회", full["null_iters"] >= 1000)]
    env = make_result(out, "smina --score_only 항 추출 + 표준화 OLS",
                      f"n={n}, 항 {len(TERMS)}개", checks,
                      notes=(f"5항 Q2 {full['Q2_loo']}, 단일점수 Q2 "
                             f"{models['single_vina_score']['Q2_loo']}, "
                             f"Tanimoto 단독 Q2 {models['tanimoto_only']['Q2_loo']}."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\nn={n}")
    print(f"{'모델':<26}{'R2':>8}{'Q2':>8}{'섞기중앙':>10}{'섞기95%':>10}")
    for k, m in models.items():
        print(f"{k:<26}{m['R2_fit']:>+8.3f}{m['Q2_loo']:>+8.3f}"
              f"{m['null_R2_median']:>+10.3f}{m['null_R2_p95']:>+10.3f}")
    print(f"\n{'항':<12}{'표준화계수':>11}{'단변량r':>10}{'r vs Tanimoto':>15}"
          f"{'순열Q2하락':>12}{'VIF':>7}")
    for t in SHORT:
        print(f"{t:<12}{std[t]:>+11.3f}{uni[t]['pearson_r']:>+10.3f}"
              f"{uni[t]['pearson_r_vs_tanimoto']:>+15.3f}"
              f"{perm[t]['mean_Q2_drop']:>+12.4f}{str(vif[t]):>7}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
