#!/usr/bin/env python3
"""어떤 스코어 항이 실제로 기여하는가 — 네 각도에서 본다.

Table 6 의 비표준화 계수는 서로 비교할 수 없다. 가우스 항은 수백 단위이고 수소결합 항은
1 미만이라, 계수의 절댓값이 큰 것이 중요한 것이 아니라 단위가 작은 것일 뿐이다.
따라서 다음 네 가지를 함께 낸다.

  1. 표준화 계수   — 예측변수를 z 점수로 바꾼 뒤의 계수. 서로 비교 가능해진다.
  2. 단변량 상관   — 그 항 하나만으로 pIC50 을 얼마나 설명하는가.
  3. 순열 중요도   — 그 항의 값만 섞었을 때 예측이 얼마나 나빠지는가.
  4. 공선성(VIF)   — 항끼리 겹치면 개별 계수는 해석할 수 없다.

네 지표가 서로 다른 답을 준다면, 그것이 결론이다 — n=30 에서 개별 항의 기여를 특정할 수 없다.
"""
from __future__ import annotations
import argparse, json, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
SHORT = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]


def zscore(col):
    n = len(col); m = sum(col) / n
    sd = (sum((v - m) ** 2 for v in col) / (n - 1)) ** 0.5
    return [(v - m) / sd for v in col] if sd else [0.0] * n, m, sd


def ols(X, y):
    """정규방정식 + 가우스 소거. 외부 의존 없이 돌아야 한다."""
    k = len(X[0]); n = len(y)
    A = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] + 
         [sum(X[i][a] * y[i] for i in range(n))] for a in range(k)]
    for c in range(k):
        piv = max(range(c, k), key=lambda r: abs(A[r][c]))
        if abs(A[piv][c]) < 1e-12:
            return None
        A[c], A[piv] = A[piv], A[c]
        for r in range(k):
            if r == c: continue
            f = A[r][c] / A[c][c]
            for j in range(c, k + 1):
                A[r][j] -= f * A[c][j]
    return [A[i][k] / A[i][i] for i in range(k)]


def r2(y, pred):
    m = sum(y) / len(y)
    ss = sum((a - m) ** 2 for a in y)
    return 1 - sum((a - b) ** 2 for a, b in zip(y, pred)) / ss if ss else float("nan")


def pearson(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    nu = sum((a - mx) * (b - my) for a, b in zip(x, y))
    de = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
    return nu / de if de else float("nan")


def q2_loo(X, y):
    n = len(y); preds = []
    for i in range(n):
        Xt = [X[j] for j in range(n) if j != i]; yt = [y[j] for j in range(n) if j != i]
        b = ols(Xt, yt)
        if b is None: return float("nan")
        preds.append(sum(c * v for c, v in zip(b, X[i])))
    return r2(y, preds)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(SR / "term_importance.json"))
    ap.add_argument("--perm", type=int, default=2000)
    a = ap.parse_args()

    reg = json.loads((SR / "regression.json").read_text()); reg = reg.get("result", reg)
    rows = reg["rows"]; n = len(rows)
    y = [r["pIC50"] for r in rows]
    raw = {t: [r[t] for r in rows] for t in SHORT}
    z = {t: zscore(raw[t])[0] for t in SHORT}

    # 1) 표준화 계수 (y 도 표준화해 완전 비교 가능하게)
    yz, ymean, ysd = zscore(y)
    Xz = [[1.0] + [z[t][i] for t in SHORT] for i in range(n)]
    bz = ols(Xz, yz)
    std_coef = dict(zip(SHORT, [round(v, 3) for v in bz[1:]])) if bz else {}
    full_r2 = r2(yz, [sum(c * v for c, v in zip(bz, x)) for x in Xz]) if bz else float("nan")
    full_q2 = q2_loo(Xz, yz)

    # 2) 단변량
    uni = {t: {"pearson_r": round(pearson(raw[t], y), 3),
               "r2_alone": round(pearson(raw[t], y) ** 2, 3)} for t in SHORT}

    # 3) 순열 중요도 — 해당 항만 섞고 Q2 하락을 본다
    rng = random.Random(42); perm_imp = {}
    for ti, t in enumerate(SHORT):
        drops = []
        for _ in range(max(20, a.perm // 100)):
            col = z[t][:]; rng.shuffle(col)
            Xp = [[1.0] + [(col[i] if s == t else z[s][i]) for s in SHORT] for i in range(n)]
            qq = q2_loo(Xp, yz)
            if qq == qq: drops.append(full_q2 - qq)
        perm_imp[t] = {"mean_Q2_drop": round(sum(drops) / len(drops), 3) if drops else None,
                       "n_shuffles": len(drops)}

    # 4) 공선성 — 각 항을 나머지로 회귀했을 때의 VIF
    vif = {}
    for t in SHORT:
        others = [s for s in SHORT if s != t]
        Xo = [[1.0] + [z[s][i] for s in others] for i in range(n)]
        b = ols(Xo, z[t])
        rr = r2(z[t], [sum(c * v for c, v in zip(b, x)) for x in Xo]) if b else float("nan")
        vif[t] = round(1 / (1 - rr), 1) if rr == rr and rr < 0.999 else None

    # 항끼리 상관 행렬
    cmat = {t: {s: round(pearson(raw[t], raw[s]), 2) for s in SHORT} for t in SHORT}

    out = {"n": n, "standardized_coefficients": std_coef,
           "full_model_R2_standardized": round(full_r2, 3),
           "full_model_Q2_standardized": round(full_q2, 3),
           "univariate": uni, "permutation_importance": perm_imp,
           "vif": vif, "term_correlation": cmat,
           "note": ("표준화 계수는 서로 비교 가능하지만, VIF 가 높으면 개별 계수의 해석 자체가 "
                    "불안정하다. 네 지표가 일치하지 않으면 기여를 특정할 수 없다는 뜻이다.")}
    ranks = {t: sorted(SHORT, key=lambda s: -abs(std_coef.get(s, 0))).index(t) + 1 for t in SHORT}
    uranks = {t: sorted(SHORT, key=lambda s: -abs(uni[s]["pearson_r"])).index(t) + 1 for t in SHORT}
    out["rank_agreement"] = {"standardized": ranks, "univariate": uranks,
                             "same_top_term": (min(ranks, key=ranks.get)
                                               == min(uranks, key=uranks.get))}
    checks = [("전 항 표준화 계수 산출", len(std_coef) == len(SHORT)),
              ("전 항 VIF 산출", all(v is not None for v in vif.values())),
              ("순열 중요도 산출", all(v["mean_Q2_drop"] is not None
                                       for v in perm_imp.values()))]
    env = make_result(out, "표준화 OLS + 순열 중요도 + VIF", f"n={n}, 항 {len(SHORT)}개", checks,
                      notes=(f"최대 VIF {max(v for v in vif.values() if v)}. "
                             f"표준화·단변량 1위 항 일치: {out['rank_agreement']['same_top_term']}."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"n={n}  표준화 R2 {full_r2:+.3f}  Q2 {full_q2:+.3f}")
    print(f"{'항':<12}{'표준화계수':>10}{'단변량r':>10}{'순열Q2하락':>12}{'VIF':>8}")
    for t in SHORT:
        print(f"{t:<12}{std_coef.get(t,0):>+10.3f}{uni[t]['pearson_r']:>+10.3f}"
              f"{perm_imp[t]['mean_Q2_drop']:>+12.3f}{str(vif[t]):>8}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
