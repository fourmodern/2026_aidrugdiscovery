#!/usr/bin/env python3
"""적합 가중치를 smina 커스텀 스코어로 넣어 held-out 세트를 다시 도킹한다.

이전 n=10 시험은 검정력이 없었다 (차이의 95% CI 가 [-0.93, +0.73]). 여기서는 통제 세트를
**층화 분할**해 훈련/시험을 나눈다 — 9칸 격자의 각 칸에서 같은 비율로 뽑아 두 집합이 같은
역가·골격 분포를 갖게 한다. 무작위 분할은 칸이 얇은 곳에서 쏠린다.

순환 방지: 가중치는 **훈련 집합에서만** 적합하고, 시험 집합은 가중치 적합에 일절 쓰이지 않는다.
비교는 같은 시험 화합물을 두 채점 함수로 각각 도킹해서 한다 (화합물·상자·시드 동일).
"""
from __future__ import annotations
import argparse, json, os, random, subprocess, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result                                    # noqa: E402
from analyze_controlled import spearman, fisher_ci                # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"
WORK2 = SR / "structures" / "work_controlled"
OUT = SR / "structures" / "custom"
SMINA = os.environ.get("SMINA", "/home/hjpark/vina/smina")
TERM_NAMES = ["gauss(o=0,_w=0.5,_c=8)", "gauss(o=3,_w=2,_c=8)", "repulsion(o=0,_c=8)",
              "hydrophobic(g=0.5,_b=1.5,_c=8)", "non_dir_h_bond(g=-0.7,_b=0,_c=8)"]
SHORT = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]


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


def dock(job):
    cid, tag, scoring, exh, seed = job
    OUT.mkdir(parents=True, exist_ok=True)
    lig = WORK2 / f"{cid}.sdf"
    if not lig.exists():
        return {"chembl_id": cid, "arm": tag, "error": "리간드 SDF 없음"}
    o = OUT / f"{cid}_{tag}.sdf"
    cmd = [SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(lig),
           "--autobox_ligand", str(WORK / "ligand_ref.sdf"), "--autobox_add", "3",
           "-o", str(o), "--exhaustiveness", str(exh), "--seed", str(seed),
           "--num_modes", "9", "--cpu", "4"]
    if scoring:
        cmd += ["--custom_scoring", scoring]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=5400)
    except subprocess.TimeoutExpired:
        return {"chembl_id": cid, "arm": tag, "error": "timeout"}
    sc = [float(x.split()[1]) for x in p.stdout.splitlines()
          if len(x.split()) >= 2 and x.split()[0].isdigit()]
    if not sc:
        return {"chembl_id": cid, "arm": tag, "error": f"점수 없음 rc={p.returncode}"}
    return {"chembl_id": cid, "arm": tag, "top_score": sc[0], "n_modes": len(sc)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-frac", type=float, default=0.30)
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--out", default=str(SR / "custom_scoring_controlled.json"))
    a = ap.parse_args()

    tc = json.loads((SR / "terms_controlled.json").read_text())["result"]
    rows = tc["rows"]
    # 9칸 층화 분할 — 무작위 분할은 얇은 칸에서 쏠린다
    rng = random.Random(a.seed)
    cells = defaultdict(list)
    for r in rows:
        cells[(r["potency_bin"], r["similarity_bin"])].append(r)
    train, test = [], []
    for key, g in sorted(cells.items()):
        g = sorted(g, key=lambda x: x["chembl_id"]); rng.shuffle(g)
        k = max(1, int(round(len(g) * a.test_frac)))
        test += g[:k]; train += g[k:]
    print(f"층화 분할: 훈련 {len(train)} / 시험 {len(test)}", file=sys.stderr)

    def zs(v):
        n = len(v); m = sum(v) / n
        sd = (sum((x - m) ** 2 for x in v) / (n - 1)) ** 0.5
        return [(x - m) / sd for x in v] if sd else [0.0] * n, m, sd

    y = [r["pchembl_value"] for r in train]
    cols, stats_ = [], []
    for t in SHORT:
        z, m, sd = zs([r[t] for r in train]); cols.append(z); stats_.append((m, sd))
    X = [[1.0] + [cols[j][i] for j in range(len(SHORT))] for i in range(len(train))]
    b = ols(X, y)
    if b is None:
        raise SystemExit("훈련 회귀 실패")
    # 표준화 계수를 원 척도로 되돌린 뒤 부호를 뒤집는다.
    # 회귀는 pIC50(클수록 좋음)을 예측하고 smina 점수는 에너지(작을수록 좋음)라 부호가 반대다.
    raw_w = [b[j + 1] / stats_[j][1] for j in range(len(SHORT))]
    scoring = OUT / "fitted_scoring.txt"
    OUT.mkdir(parents=True, exist_ok=True)
    scoring.write_text("\n".join(f"{-w:.8f} {t}" for w, t in zip(raw_w, TERM_NAMES)) + "\n")
    print("적합 가중치 (부호 반전, 원 척도):", file=sys.stderr)
    for w, t in zip(raw_w, TERM_NAMES):
        print(f"   {-w:+.6f}  {t}", file=sys.stderr)

    jobs = ([(r["chembl_id"], "base", None, a.exhaustiveness, a.seed) for r in test] +
            [(r["chembl_id"], "custom", str(scoring), a.exhaustiveness, a.seed) for r in test])
    got = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(dock, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if "error" not in r:
                got[(r["chembl_id"], r["arm"])] = r["top_score"]
            if i % 20 == 0:
                print(f"  {i}/{len(jobs)}", file=sys.stderr)

    ids = [r["chembl_id"] for r in test
           if (r["chembl_id"], "base") in got and (r["chembl_id"], "custom") in got]
    pv = {r["chembl_id"]: r["pchembl_value"] for r in test}
    yy = [pv[i] for i in ids]
    bb = [got[(i, "base")] for i in ids]
    cc = [got[(i, "custom")] for i in ids]
    rb, rc = spearman(yy, bb), spearman(yy, cc)

    rng2 = random.Random(a.seed); diffs = []
    for _ in range(10000):
        idx = [rng2.randrange(len(ids)) for _ in ids]
        try:
            d = (spearman([yy[i] for i in idx], [cc[i] for i in idx])
                 - spearman([yy[i] for i in idx], [bb[i] for i in idx]))
            if d == d: diffs.append(d)
        except Exception:
            pass
    diffs.sort()
    ci = [round(diffs[int(.025 * len(diffs))], 3), round(diffs[int(.975 * len(diffs))], 3)]
    favor = round(sum(1 for d in diffs if d < 0) / len(diffs), 3)

    out = {"design": ("가중치는 훈련 집합에서만 적합. 시험 집합은 9칸 층화로 분리했고 "
                      "가중치 적합에 쓰이지 않았다."),
           "n_train": len(train), "n_test": len(test), "n_scored": len(ids),
           "test_frac": a.test_frac, "split_seed": a.seed,
           "fitted_weights_raw_scale": {t: round(-w, 6) for t, w in zip(TERM_NAMES, raw_w)},
           "arm_default": {"spearman": round(rb, 3), "ci95": fisher_ci(rb, len(ids))},
           "arm_custom": {"spearman": round(rc, 3), "ci95": fisher_ci(rc, len(ids))},
           "difference": {"delta_spearman": round(rc - rb, 3), "bootstrap_ci95": ci,
                          "frac_resamples_favoring_custom": favor,
                          "bootstrap_iters": len(diffs), "bootstrap_seed": a.seed},
           "rows": [{"chembl_id": i, "pchembl_value": pv[i],
                     "base_score": got[(i, "base")], "custom_score": got[(i, "custom")]}
                    for i in ids],
           "reading": ("차이의 신뢰구간이 0 을 포함하면 개선을 주장할 수 없다. "
                       "점 추정값만 보고하면 우연한 방향 일치를 성능으로 오인한다.")}
    checks = [("두 팔 모두 전 시험 화합물 채점", len(ids) == len(test)),
              ("훈련/시험 분리", not ({r["chembl_id"] for r in train}
                                      & {r["chembl_id"] for r in test})),
              ("시험 집합 20건 이상", len(ids) >= 20),
              ("부트스트랩 5000회 이상", len(diffs) >= 5000)]
    env = make_result(out, "smina --custom_scoring, 9칸 층화 분할",
                      f"훈련 {len(train)} / 시험 {len(ids)}", checks,
                      notes=(f"기본 {rb:+.3f} vs 커스텀 {rc:+.3f}, 차이 {rc-rb:+.3f}, "
                             f"부트스트랩 95% CI {ci}."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n기본 Spearman {rb:+.3f}  |  커스텀 {rc:+.3f}  |  차이 {rc-rb:+.3f}")
    print(f"부트스트랩 95% CI {ci}  ·  커스텀 우세 재표본 비율 {favor:.3f}")
    print(f"게이트 {'PASS' if env['verification']['passed'] else 'FAIL'}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
