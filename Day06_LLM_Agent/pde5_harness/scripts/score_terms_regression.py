#!/usr/bin/env python3
"""도킹 스코어 항으로 pIC50 회귀 — 작은 표본의 한계를 정면으로 다룬다.

단일 Vina 점수는 congeneric series 의 역가 순위를 잘 재현하지 못하는 것으로 알려져 있다.
여기서는 점수를 구성하는 개별 항을 꺼내 회귀에 넣고, 단일 점수보다 나은지 본다.

**표본이 30건이다.** 항이 5개면 자유도가 넉넉하지 않다. 따라서
  - 항상 leave-one-out 교차검증(Q2)을 함께 보고하고,
  - 학습 R2 만으로 성능을 주장하지 않으며,
  - 무작위 라벨 대조(y-scrambling)로 우연 수준을 함께 낸다.
이 셋이 없으면 소표본 회귀는 과적합을 성능으로 오인하게 만든다.
"""
from __future__ import annotations
import argparse, json, re, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402
WORK = ROOT / "sample_run" / "structures" / "work"
SMINA = "/home/hjpark/vina/smina"

# Vina 기본 스코어의 5개 항 (smina --score_only 가 이 순서로 낸다)
TERMS = ["gauss(o=0,_w=0.5,_c=8)", "gauss(o=3,_w=2,_c=8)", "repulsion(o=0,_c=8)",
         "hydrophobic(g=0.5,_b=1.5,_c=8)", "non_dir_h_bond(g=-0.7,_b=0,_c=8)"]
SHORT = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]


def score_terms(receptor: Path, pose: Path) -> dict | None:
    """한 자세의 항별 값을 뽑는다. smina --score_only 출력에서 파싱."""
    p = subprocess.run([SMINA, "-r", str(receptor), "-l", str(pose), "--score_only"],
                       capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        return None
    vals, grab = {}, False
    for line in p.stdout.splitlines():
        if line.startswith("Term values, before weighting:"):
            grab = True
            nums = re.findall(r"-?\d+\.?\d*(?:e[-+]?\d+)?", line.split(":", 1)[1])
            if nums:
                for k, v in zip(SHORT, nums):
                    vals[k] = float(v)
            continue
        if grab and not vals:
            nums = re.findall(r"-?\d+\.?\d*(?:e[-+]?\d+)?", line)
            if len(nums) >= len(SHORT):
                for k, v in zip(SHORT, nums):
                    vals[k] = float(v)
            grab = False
    return vals or None


def ols(X, y):
    """절편 포함 최소제곱. 의존성을 늘리지 않으려고 정규방정식을 직접 푼다."""
    n, k = len(X), len(X[0])
    A = [[1.0] + list(row) for row in X]
    AT = list(zip(*A))
    ATA = [[sum(AT[i][m] * AT[j][m] for m in range(n)) for j in range(k + 1)]
           for i in range(k + 1)]
    ATy = [sum(AT[i][m] * y[m] for m in range(n)) for i in range(k + 1)]
    # 가우스 소거 (부분 피벗) + 미세 리지로 특이행렬 회피
    for i in range(k + 1):
        ATA[i][i] += 1e-8
    M = [row[:] + [ATy[i]] for i, row in enumerate(ATA)]
    for c in range(k + 1):
        piv = max(range(c, k + 1), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-12:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(k + 1):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for cc in range(c, k + 2):
                M[r][cc] -= f * M[c][cc]
    return [M[i][k + 1] / M[i][i] for i in range(k + 1)]


def predict(b, row):
    return b[0] + sum(bi * xi for bi, xi in zip(b[1:], row))


def r2(y, yhat):
    m = sum(y) / len(y)
    ss_res = sum((a - b) ** 2 for a, b in zip(y, yhat))
    ss_tot = sum((a - m) ** 2 for a in y)
    return 1 - ss_res / ss_tot if ss_tot else float("nan")


def loo_q2(X, y):
    preds = []
    for i in range(len(X)):
        Xt = X[:i] + X[i + 1:]; yt = y[:i] + y[i + 1:]
        b = ols(Xt, yt)
        if b is None:
            return float("nan"), []
        preds.append(predict(b, X[i]))
    return r2(y, preds), preds


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--docking", default=str(ROOT / "sample_run" / "docking.json"))
    ap.add_argument("--out", default=str(ROOT / "sample_run" / "regression.json"))
    ap.add_argument("--scramble", type=int, default=200)
    a = ap.parse_args()

    env = json.loads(Path(a.docking).read_text())
    if not env["verification"]["passed"]:
        sys.exit("도킹 게이트가 통과하지 않았다 — 회귀를 수행하지 않는다.")
    docked = env["result"]["docked"]
    receptor = WORK / "receptor.pdbqt"

    rows = []
    for d in docked:
        pose = WORK / f"{d['chembl_id']}_selected.sdf"
        if not pose.exists():
            pose = WORK / f"{d['chembl_id']}_pose.sdf"
        if not pose.exists():
            continue
        t = score_terms(receptor, pose)
        if t is None:
            continue
        rows.append({**d, **t})
    if len(rows) < 20:
        sys.exit(f"항을 뽑은 화합물이 {len(rows)}건뿐 — 20건 미만이면 회귀를 수행하지 않는다.")

    y = [r["pIC50"] for r in rows]
    single = [[r["dock_score"]] for r in rows]
    multi = [[r[k] for k in SHORT] for r in rows]

    def _ranks(v):
        """동점은 평균 순위를 준다. 이것을 빼먹으면 Spearman 값이 틀린다.

        [수정] 이전 판은 정렬 위치를 그대로 순위로 썼다. dock_score 에 -9.5 가
        다섯 번 나오는 등 동점이 많아 보고값이 실제와 달랐다 (-0.498 vs -0.474).
        """
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v); i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1                      # 1-based 평균 순위
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    def spearman(x, y):
        """동점 보정 Spearman. d^2 공식은 동점에서 성립하지 않으므로 Pearson-on-ranks 로 계산."""
        rx, ry = _ranks(x), _ranks(y); n = len(x)
        mx, my = sum(rx) / n, sum(ry) / n
        num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
        return round(num / den, 3) if den else float("nan")

    def spearman_perm_p(x, y, iters=10000, seed=42):
        """순열 검정으로 양측 p 값. n=30 에서 정규근사 대신 직접 센다."""
        import random
        obs = abs(spearman(x, y)); rnd = random.Random(seed)
        ys = list(y); hit = 0
        for _ in range(iters):
            rnd.shuffle(ys)
            if abs(spearman(x, ys)) >= obs:
                hit += 1
        return round((hit + 1) / (iters + 1), 4)

    score = [r["dock_score"] for r in rows]
    out = {"n": len(rows), "terms": SHORT, "rows": rows, "models": {},
           "correlation": {
               "method": "tie-corrected Spearman (average ranks, Pearson on ranks)",
               "spearman_pIC50_vs_dock_score": spearman(y, score),
               "perm_p_dock_score": spearman_perm_p(y, score),
               "spearman_pIC50_vs_top_pose_score":
                   spearman(y, [r.get("top_pose_score", r["dock_score"]) for r in rows]),
               "perm_p_top_pose_score":
                   spearman_perm_p(y, [r.get("top_pose_score", r["dock_score"]) for r in rows]),
               "perm_iters": 10000, "perm_seed": 42,
           },
           "strata": {s: sum(1 for r in rows if r.get("stratum") == s)
                      for s in ("strong", "medium", "weak")}}
    for name, X in (("single_vina_score", single), ("five_terms", multi)):
        b = ols(X, y)
        if b is None:
            out["models"][name] = {"error": "특이행렬"}
            continue
        fit = [predict(b, x) for x in X]
        q2, _ = loo_q2(X, y)
        # y-scrambling: 라벨을 섞었을 때 얻어지는 R2 분포
        import random
        rnd = random.Random(42); null = []   # 시드 42, 셔플 200회
        for _ in range(a.scramble):
            ys = y[:]; rnd.shuffle(ys)
            bb = ols(X, ys)
            if bb: null.append(r2(ys, [predict(bb, x) for x in X]))
        null.sort()
        out["models"][name] = {
            "coefficients": [round(v, 4) for v in b],
            "R2_fit": round(r2(y, fit), 3),
            "Q2_loo": round(q2, 3),
            "null_R2_median": round(null[len(null) // 2], 3) if null else None,
            "null_R2_p95": round(null[int(len(null) * 0.95)], 3) if null else None,
            "n_params": len(b),
        }
    checks = [
        ("전 화합물에서 항 추출 성공", all(all(t in r for t in SHORT) for r in out["rows"])),
        ("두 모델 모두 LOO 수행", all(m.get("Q2_loo") is not None
                                      for m in out["models"].values())),
        ("라벨섞기 귀무분포 산출", all(m.get("null_R2_median") is not None
                                       for m in out["models"].values())),
        ("동점 보정 Spearman 사용", "tie-corrected" in out["correlation"]["method"]),
    ]
    env = make_result(out, "smina --score_only 항 추출 + OLS",
                      f"n={out['n']}, 항 {len(TERMS)}개", checks,
                      notes=(f"Spearman {out['correlation']['spearman_pIC50_vs_dock_score']} "
                             f"(p={out['correlation']['perm_p_dock_score']}). "
                             "Q2 가 음수면 과적합이다."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")

    print(f"n = {out['n']}  (pIC50 {min(y):.2f}–{max(y):.2f})  층 {out['strata']}")
    c = out["correlation"]
    print(f"  Spearman(동점보정)  선택자세 {c['spearman_pIC50_vs_dock_score']:+.3f} "
          f"(순열 p={c['perm_p_dock_score']})  |  1위자세 "
          f"{c['spearman_pIC50_vs_top_pose_score']:+.3f} (p={c['perm_p_top_pose_score']})")
    for name, m in out["models"].items():
        if "error" in m:
            print(f"  {name}: {m['error']}"); continue
        print(f"  {name:18s} R2_fit {m['R2_fit']:+.3f}  Q2_loo {m['Q2_loo']:+.3f}  "
              f"파라미터 {m['n_params']}  라벨섞기 R2 중앙 {m['null_R2_median']:+.3f} "
              f"/ 95% {m['null_R2_p95']:+.3f}")
    print(f"\n{Path(a.out)} 저장")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
