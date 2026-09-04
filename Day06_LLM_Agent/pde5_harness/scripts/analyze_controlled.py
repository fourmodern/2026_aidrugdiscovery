#!/usr/bin/env python3
"""골격 통제 설계의 결정적 분석 — 도킹 점수가 역가를 보는가 골격을 보는가.

이전 결과의 문제는 두 축이 섞여 있어 "선별 성능"의 정체를 알 수 없었다는 것이다.
이 설계에서는 세 가지로 판정한다.

  1. 유사도 구간 내 상관 — 골격이 고정된 상태에서 역가 신호가 남는가. **이것이 결정적이다.**
  2. 편상관            — Tanimoto 와 **분자 크기**를 통제한 뒤의 역가-점수 상관.
                        Vina 계열 점수는 크기와 상관되므로 골격만 통제하면 부족하다.
  3. 이중 회귀         — pIC50 ~ 점수 + Tanimoto. 점수 계수가 살아남는가.

세 답이 일치하면 결론을 신뢰할 수 있다. 갈리면 그 사실이 결론이다.
"""
from __future__ import annotations
import argparse, json, math, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"


def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3: return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    nu = sum((a - mx) * (b - my) for a, b in zip(x, y))
    de = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
    return nu / de if de else float("nan")


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def perm_p(x, y, iters=20000, seed=42):
    obs = abs(spearman(x, y))
    if obs != obs: return None
    rng = random.Random(seed); yy = list(y); c = 0
    for _ in range(iters):
        rng.shuffle(yy)
        if abs(spearman(x, yy)) >= obs - 1e-12:
            c += 1
    return round((c + 1) / (iters + 1), 5)


def fisher_ci(r, n):
    if r != r or abs(r) >= 1 or n < 4: return None
    z = 0.5 * math.log((1 + r) / (1 - r)); se = 1 / math.sqrt(n - 3)
    return [round(math.tanh(z - 1.96 * se), 3), round(math.tanh(z + 1.96 * se), 3)]


def partial_spearman(x, y, z):
    """z 를 통제한 x-y 순위 편상관. 순위 공간에서 계산한다."""
    rxy, rxz, ryz = spearman(x, y), spearman(x, z), spearman(y, z)
    den = ((1 - rxz ** 2) * (1 - ryz ** 2)) ** 0.5
    return (rxy - rxz * ryz) / den if den else float("nan")


def partial_perm_p(x, y, z, iters=20000, seed=42):
    """편상관의 순열 검정 — Freedman-Lane 절차.

    y 를 통째로 섞으면 y~z 관계까지 깨진다. 그러면 귀무가설이 "y 는 x 와도 z 와도 무관"이
    되어, 실제로 검정하고 싶은 "x 가 z 너머의 정보를 주지 않는다"와 달라진다. z 와 y 가 강하게
    엮여 있을수록 그 차이가 커지고 p 값이 낙관적으로 나온다.

    Freedman-Lane 은 이렇게 한다. (1) y 를 z 로 회귀해 적합값과 잔차로 나누고,
    (2) **잔차만** 섞은 뒤, (3) 적합값에 도로 더해 y* 를 만든다. y*~z 관계는 보존되고
    x 와의 추가 관계만 파괴된다. Spearman 을 쓰므로 순위 공간에서 수행한다.
    """
    obs = abs(partial_spearman(x, y, z))
    if obs != obs: return None
    ry, rz = ranks(y), ranks(z)
    n = len(ry)
    mz, my = sum(rz) / n, sum(ry) / n
    den = sum((v - mz) ** 2 for v in rz)
    b = (sum((a - mz) * (c - my) for a, c in zip(rz, ry)) / den) if den else 0.0
    a0 = my - b * mz
    fit = [a0 + b * v for v in rz]
    resid = [c - f for c, f in zip(ry, fit)]
    rng = random.Random(seed); idx = list(range(n)); cnt = 0
    for _ in range(iters):
        rng.shuffle(idx)
        ystar = [f + resid[i] for f, i in zip(fit, idx)]
        if abs(partial_spearman(x, ystar, z)) >= obs - 1e-12:
            cnt += 1
    return round((cnt + 1) / (iters + 1), 5)


def partial_spearman_multi(x, y, covs):
    """여러 공변량을 동시에 통제한 순위 편상관.

    Vina 계열 점수는 분자 크기와 상관되는 것으로 알려져 있다 [R5]. 골격 유사성만 통제하고
    크기를 놓치면 같은 종류의 오독이 반복된다. 순위 공간에서 x 와 y 를 각각 공변량들로
    회귀한 뒤, 그 잔차끼리의 상관을 낸다.
    """
    rx, ry = ranks(x), ranks(y)
    rc = [ranks(c) for c in covs]
    n = len(rx)
    X = [[1.0] + [rc[j][i] for j in range(len(rc))] for i in range(n)]
    def resid(v):
        b = ols(X, v)
        if b is None: return None
        return [vi - sum(bk * xk for bk, xk in zip(b, xi)) for vi, xi in zip(v, X)]
    ex, ey = resid(rx), resid(ry)
    if ex is None or ey is None: return float("nan")
    return pearson(ex, ey)


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


def zs(v):
    n = len(v); m = sum(v) / n
    sd = (sum((a - m) ** 2 for a in v) / (n - 1)) ** 0.5
    return [(a - m) / sd for a in v] if sd else [0.0] * n


def auc(scores, labels):
    pos = [-s for s, l in zip(scores, labels) if l]
    neg = [-s for s, l in zip(scores, labels) if not l]
    if not pos or not neg: return None
    return round(sum(1.0 if a > b else 0.5 if a == b else 0.0
                     for a in pos for b in neg) / (len(pos) * len(neg)), 3)


def auc_ci(a, n_pos, n_neg):
    """Hanley-McNeil 근사 신뢰구간. AUC 를 CI 없이 보고하면 표본이 작을 때 오독된다."""
    if a is None or n_pos < 2 or n_neg < 2:
        return None
    q1 = a / (2 - a); q2 = 2 * a * a / (1 + a)
    se = ((a * (1 - a) + (n_pos - 1) * (q1 - a * a)
           + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg)) ** 0.5
    return [round(max(0.0, a - 1.96 * se), 3), round(min(1.0, a + 1.96 * se), 3)]


MCS_FLOOR = 8          # dock_controlled.py 의 MCS 최소 원자수
MCS_INTERPRETABLE = 15  # 이보다 적으면 RMSD 를 자세 정확도로 읽을 수 없다


def pose_validity(rows, key="top_pose_mcs_rmsd", thr=2.0):
    """점수를 매긴 자세가 결정 자세의 틀 안에 있는가.

    **주의: 이 지표는 공유 부분구조가 적으면 성립하지 않는다.** MCS 원자수가 하한(8) 근처면
    RMSD 는 "자세가 맞는가" 가 아니라 "우연히 겹치는 원자 몇 개가 어디 있나" 를 재게 된다.
    그래서 해석 가능 여부를 함께 돌려주고, 통합값은 해석 가능한 대역만으로도 따로 낸다.
    """
    v = [r[key] for r in rows if r.get(key) is not None]
    if not v:
        return None
    a = [r["mcs_atoms"] for r in rows if r.get(key) is not None and r.get("mcs_atoms")]
    v2 = sorted(v)
    mean_atoms = round(sum(a) / len(a), 1) if a else None
    at_floor = sum(1 for x in a if x <= MCS_FLOOR)
    return {"n": len(v), "median_rmsd": round(v2[len(v2) // 2], 3),
            "frac_under_threshold": round(sum(1 for x in v if x < thr) / len(v), 3),
            "threshold_angstrom": thr,
            "mean_mcs_atoms": mean_atoms, "n_at_mcs_floor": at_floor,
            "interpretable": bool(mean_atoms and mean_atoms >= MCS_INTERPRETABLE
                                  and at_floor == 0),
            "note": ("MCS 원자수가 적어 RMSD 를 자세 정확도로 읽을 수 없다"
                     if not (mean_atoms and mean_atoms >= MCS_INTERPRETABLE and at_floor == 0)
                     else "MCS 원자수가 충분해 해석 가능")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dock", default=str(SR / "docking_controlled.json"))
    ap.add_argument("--out", default=str(SR / "analysis_controlled.json"))
    a = ap.parse_args()

    D = json.loads(Path(a.dock).read_text())["result"]
    rows = D["rows"]; n = len(rows)
    y = [r["pchembl_value"] for r in rows]
    tan = [r["tanimoto_to_sildenafil"] for r in rows]
    hac = [r["heavy_atoms"] for r in rows]
    arms = {"top_pose": [r["top_pose_score"] for r in rows]}
    ref_ok = [r for r in rows if r.get("ref_pose_score") is not None]

    res = {"n": n, "design": D.get("pose_selection"), "arms": {}}

    for nm, sc in arms.items():
        rho = spearman(y, sc)
        res["arms"][nm] = {
            "spearman": round(rho, 3), "ci95": fisher_ci(rho, n), "perm_p": perm_p(y, sc),
            "spearman_score_vs_tanimoto": round(spearman(sc, tan), 3),
            "perm_p_score_vs_tanimoto": perm_p(sc, tan),
            "spearman_pIC50_vs_tanimoto": round(spearman(y, tan), 3),
            "perm_p_pIC50_vs_tanimoto": perm_p(y, tan),
            "partial_spearman_controlling_tanimoto": round(partial_spearman(y, sc, tan), 3),
            "partial_perm_p": partial_perm_p(y, sc, tan),
            "spearman_score_vs_heavy_atoms": round(spearman(sc, hac), 3),
            "spearman_pIC50_vs_heavy_atoms": round(spearman(y, hac), 3),
            "partial_spearman_controlling_size": round(partial_spearman(y, sc, hac), 3),
            "partial_spearman_controlling_both": round(
                partial_spearman_multi(sc, y, [tan, hac]), 3),
        }
        # 유사도 구간 내 상관 — 결정적 검정
        within = {}
        for b in ("far", "mid", "near"):
            idx = [i for i in range(n) if rows[i]["similarity_bin"] == b]
            if len(idx) < 6: continue
            yy = [y[i] for i in idx]; ss = [sc[i] for i in idx]
            r_ = spearman(yy, ss)
            within[b] = {"n": len(idx), "spearman": round(r_, 3),
                         "ci95": fisher_ci(r_, len(idx)), "perm_p": perm_p(yy, ss),
                         "pIC50_range": [round(min(yy), 2), round(max(yy), 2)]}
        res["arms"][nm]["within_similarity_bin"] = within
        vals = [v["spearman"] for v in within.values()]
        res["arms"][nm]["within_bin_mean_spearman"] = round(sum(vals) / len(vals), 3) if vals else None
        # 이중 회귀 (표준화)
        Z = [[1.0, a_, b_] for a_, b_ in zip(zs(sc), zs(tan))]
        b = ols(Z, zs(y))
        res["arms"][nm]["dual_regression_standardized"] = (
            {"intercept": round(b[0], 3), "beta_score": round(b[1], 3),
             "beta_tanimoto": round(b[2], 3)} if b else None)
        # 선별 지표
        sw = [(s, 1 if r["potency_bin"] == "strong" else 0)
              for s, r in zip(sc, rows) if r["potency_bin"] in ("strong", "weak")]
        res["arms"][nm]["auc_strong_vs_weak"] = auc([s for s, _ in sw], [l for _, l in sw])
        aucb = {}
        for bnm in ("far", "mid", "near"):
            sub = [(s, 1 if r["potency_bin"] == "strong" else 0)
                   for s, r in zip(sc, rows)
                   if r["similarity_bin"] == bnm and r["potency_bin"] in ("strong", "weak")]
            if len(sub) >= 8:
                aucb[bnm] = auc([s for s, _ in sub], [l for _, l in sub])
        res["arms"][nm]["auc_within_similarity_bin"] = aucb
        k = max(1, n // 10)
        topk = sorted(range(n), key=lambda i: sc[i])[:k]
        res["arms"][nm][f"top{k}_potency_bins"] = {
            b: sum(1 for i in topk if rows[i]["potency_bin"] == b)
            for b in ("strong", "medium", "weak")}
        res["arms"][nm][f"top{k}_similarity_bins"] = {
            b: sum(1 for i in topk if rows[i]["similarity_bin"] == b)
            for b in ("far", "mid", "near")}

        # 대역별 편상관 (크기·골격 각각/동시) — 전체에만 적용하면 이질성이 숨는다
        wctl = {}
        for b in ("far", "mid", "near"):
            idx = [i for i in range(n) if rows[i]["similarity_bin"] == b]
            if len(idx) < 8: continue
            yy = [y[i] for i in idx]; ss = [sc[i] for i in idx]
            tt = [tan[i] for i in idx]; hh = [hac[i] for i in idx]
            wctl[b] = {"n": len(idx), "raw": round(spearman(yy, ss), 3),
                       "controlling_size": round(partial_spearman(yy, ss, hh), 3),
                       "controlling_tanimoto": round(partial_spearman(yy, ss, tt), 3),
                       "controlling_both": round(partial_spearman_multi(ss, yy, [tt, hh]), 3)}
        res["arms"][nm]["within_bin_partial"] = wctl

        # AUC 신뢰구간
        cis = {}
        for b, au in list(aucb.items()) + [("all", res["arms"][nm]["auc_strong_vs_weak"])]:
            sub = [r for r in rows if r["potency_bin"] in ("strong", "weak")
                   and (b == "all" or r["similarity_bin"] == b)]
            npos = sum(1 for r in sub if r["potency_bin"] == "strong")
            # 음성군이 전수(census) 칸으로만 이루어져 있으면 그 AUC 는 그 칸 없이는
            # 재계산 자체가 불가능하다. 민감도 분석의 한계를 명시하기 위해 표시한다.
            neg_rows = [r for r in sub if r["potency_bin"] == "weak"]
            cis[b] = {"auc": au, "n_pos": npos, "n_neg": len(sub) - npos,
                      "ci95": auc_ci(au, npos, len(sub) - npos),
                      "neg_class_is_single_cell": bool(neg_rows) and len(
                          {r["similarity_bin"] for r in neg_rows}) == 1 and b != "all"}
        res["arms"][nm]["auc_ci95"] = cis

        # 3개 대역 다중비교 보정 — 계산해 보고한다 (미보정 고지만으로는 부족)
        ps = [(b, v["perm_p"]) for b, v in within.items() if v.get("perm_p") is not None]
        ps_sorted = sorted(ps, key=lambda t: t[1])
        m = len(ps_sorted)
        res["arms"][nm]["multiplicity"] = {
            "family": "유사도 3대역 내 상관 검정",
            "bonferroni": {b: round(min(1.0, pv * m), 4) for b, pv in ps},
            "benjamini_hochberg": {b: round(min(1.0, pv * m / (i + 1)), 4)
                                   for i, (b, pv) in enumerate(ps_sorted)}}

        # 자세 타당도 — 점수를 매긴 자세가 결정 자세 틀 안에 있는가
        pv = {"all": pose_validity(rows),
              **{b: pose_validity([r for r in rows if r["similarity_bin"] == b])
                 for b in ("far", "mid", "near")}}
        # 해석 불가한 대역을 뺀 통합값 — 헤드라인은 이 쪽을 써야 한다.
        # 해석할 수 없다고 적은 대역을 통합값에는 넣어 두면 그 통합값도 오염된다.
        bad = [b for b in ("far", "mid", "near") if pv[b] and not pv[b]["interpretable"]]
        keep = [r for r in rows if r["similarity_bin"] not in bad]
        pv["interpretable_only"] = pose_validity(keep)
        pv["excluded_bands"] = bad
        res["arms"][nm]["pose_validity"] = pv

        # 대역별 평균 점수 — 어느 대역이 좋은 점수를 받는가 (기전 서술의 근거).
        # 차이에 p 를 붙이지 않고 서술하면, 다른 곳에서 p=0.3 을 "신호 없음"이라 한 것과
        # 기준이 어긋난다. 순열로 검정한다.
        by_bin = {b: [r["top_pose_score"] for r in rows if r["similarity_bin"] == b]
                  for b in ("far", "mid", "near")}
        res["arms"][nm]["mean_score_by_bin"] = {
            b: round(sum(v) / len(v), 3) for b, v in by_bin.items() if v}
        _rng = random.Random(42)
        obs_gap = (sum(by_bin["far"]) / len(by_bin["far"])
                   - sum(by_bin["near"]) / len(by_bin["near"]))
        pool = by_bin["far"] + by_bin["near"]; nf = len(by_bin["far"]); cnt = 0
        for _ in range(20000):
            _rng.shuffle(pool)
            g = sum(pool[:nf]) / nf - sum(pool[nf:]) / (len(pool) - nf)
            if abs(g) >= abs(obs_gap) - 1e-12: cnt += 1
        res["arms"][nm]["mean_score_gap_far_minus_near"] = {
            "gap_kcal_mol": round(obs_gap, 3), "perm_p": round((cnt + 1) / 20001, 4)}

        # 상위 k 구성의 귀무 기대값 — 대역 크기가 다르므로 n/3 이 아니다
        kk = max(1, n // 10)
        res["arms"][nm][f"top{kk}_expected_similarity"] = {
            b: round(len(v) * kk / n, 2) for b, v in by_bin.items()}

        # near 대역 민감도 — 어세이 바닥값(pChEMBL 5.00)이 결과를 만드는가
        near_idx = [i for i in range(n) if rows[i]["similarity_bin"] == "near"]
        drop_floor = [i for i in near_idx if abs(y[i] - 5.0) > 1e-9]
        drop_weak = [i for i in near_idx if rows[i]["potency_bin"] != "weak"]
        res["arms"][nm]["near_sensitivity"] = {
            "full": {"n": len(near_idx),
                     "spearman": round(spearman([y[i] for i in near_idx],
                                                [sc[i] for i in near_idx]), 3)},
            "drop_pchembl_5.00": {"n": len(drop_floor), "n_dropped": len(near_idx) - len(drop_floor),
                                  "spearman": round(spearman([y[i] for i in drop_floor],
                                                             [sc[i] for i in drop_floor]), 3)},
            "drop_weak_cell": {"n": len(drop_weak), "n_dropped": len(near_idx) - len(drop_weak),
                               "spearman": round(spearman([y[i] for i in drop_weak],
                                                          [sc[i] for i in drop_weak]), 3)},
            "note": ("weak x near 칸은 풀 전수(census)라 표본이 아니다. "
                     "그 칸의 pChEMBL 이 어세이 바닥값에 몰려 있으면 상관이 그 한 점에 얹힌다.")}

    # 참조 기준 자세 민감도
    if len(ref_ok) >= 20:
        yy = [r["pchembl_value"] for r in ref_ok]
        ss = [r["ref_pose_score"] for r in ref_ok]
        tt = [r["tanimoto_to_sildenafil"] for r in ref_ok]
        rr = spearman(yy, ss)
        res["sensitivity_reference_pose"] = {
            "n": len(ref_ok), "spearman": round(rr, 3), "ci95": fisher_ci(rr, len(ref_ok)),
            "partial_controlling_tanimoto": round(partial_spearman(yy, ss, tt), 3)}

    # 대조
    # 표본 구성 — 본문에서 "중앙값으로 합쳐 잡음을 줄였다"고 쓰려면 몇 건이 실제로
    # 합쳐졌는지 적어야 한다.
    from collections import Counter as _C
    nmeas = _C(r.get("n_measurements", 1) for r in rows)
    res["sample_composition"] = {
        "n_with_single_measurement": nmeas.get(1, 0),
        "n_with_replicates": sum(v for k, v in nmeas.items() if k > 1),
        "assay_types": dict(_C(r.get("standard_type", "?") for r in rows)),
        "near_bin_assay_types": dict(_C(r.get("standard_type", "?") for r in rows
                                        if r["similarity_bin"] == "near"))}
    res["confound_metrics"] = {
        "pearson_r_potency_vs_similarity": round(pearson(y, tan), 3),
        "spearman_rho_potency_vs_similarity": round(spearman(y, tan), 3),
        "spearman_tanimoto_vs_heavy_atoms": round(spearman(tan, hac), 3),
        "note": ("두 지표를 모두 적는다. 작은 쪽만 골라 보고하면 게이트를 통과시키려고 "
                 "고른 것으로 읽힌다. Tanimoto 와 중원자수의 상관이 높으면 '두 교란'이 "
                 "사실상 한 축이므로 동시 통제의 의미가 줄어든다.")}
    res["control_redock"] = D["control_redock"]
    T = res["arms"]["top_pose"]
    res["verdict"] = {
        "raw_signal": T["spearman"], "partial_signal": T["partial_spearman_controlling_tanimoto"],
        "partial_controlling_size": T["partial_spearman_controlling_size"],
        "partial_controlling_both": T["partial_spearman_controlling_both"],
        "within_bin_mean": T["within_bin_mean_spearman"],
        "confound_share": (round(1 - abs(T["partial_spearman_controlling_tanimoto"])
                                 / abs(T["spearman"]), 3)
                           if T["spearman"] else None),
        "reading": ("편상관과 구간내 상관이 원상관과 비슷하면 신호는 역가에서 온 것이다. "
                    "크게 줄면 골격 유사성이 원인이었다는 뜻이다.")}

    checks = [
        ("전 화합물 상관 산출", T["spearman"] == T["spearman"]),
        ("크기·골격 동시 통제 편상관 산출",
         T["partial_spearman_controlling_both"] == T["partial_spearman_controlling_both"]),
        ("세 유사도 구간 모두 분석됨", len(T["within_similarity_bin"]) == 3),
        ("편상관 산출", T["partial_spearman_controlling_tanimoto"] ==
                        T["partial_spearman_controlling_tanimoto"]),
        ("순열 p 산출", T["perm_p"] is not None and T["partial_perm_p"] is not None),
        ("자세 타당도 산출", T["pose_validity"]["all"] is not None),
        ("AUC 신뢰구간 산출", all(v.get("ci95") for v in T["auc_ci95"].values())),
        ("다중비교 보정 산출", bool(T["multiplicity"]["bonferroni"])),
    ]
    env = make_result(res, "골격 통제 설계 분석", f"n={n}", checks,
                      notes=(f"원상관 {T['spearman']}, 편상관 "
                             f"{T['partial_spearman_controlling_tanimoto']}, "
                             f"구간내 평균 {T['within_bin_mean_spearman']}."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")

    print(f"n={n}")
    print(f"  원상관   pIC50 vs 점수      {T['spearman']:+.3f}  CI{T['ci95']}  p={T['perm_p']}")
    print(f"  교란     점수 vs Tanimoto   {T['spearman_score_vs_tanimoto']:+.3f}")
    print(f"           pIC50 vs Tanimoto  {T['spearman_pIC50_vs_tanimoto']:+.3f}")
    print(f"  교란     점수 vs 중원자수    {T['spearman_score_vs_heavy_atoms']:+.3f}")
    print(f"           pIC50 vs 중원자수   {T['spearman_pIC50_vs_heavy_atoms']:+.3f}")
    print(f"  편상관   Tanimoto 통제       {T['partial_spearman_controlling_tanimoto']:+.3f}  "
          f"p={T['partial_perm_p']}")
    print(f"           크기 통제           {T['partial_spearman_controlling_size']:+.3f}")
    print(f"           둘 다 통제          {T['partial_spearman_controlling_both']:+.3f}")
    print(f"  구간내 상관:")
    for b, v in T["within_similarity_bin"].items():
        print(f"      {b:5s} n={v['n']:3d}  rho {v['spearman']:+.3f}  CI{v['ci95']}  p={v['perm_p']}")
    print(f"  AUC 전체 {T['auc_strong_vs_weak']}   구간별 {T['auc_within_similarity_bin']}")
    if T["dual_regression_standardized"]:
        d = T["dual_regression_standardized"]
        print(f"  이중회귀 (표준화)  beta_score {d['beta_score']:+.3f}   "
              f"beta_tanimoto {d['beta_tanimoto']:+.3f}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
