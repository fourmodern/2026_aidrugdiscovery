#!/usr/bin/env python3
"""우선순위(triage) 지표를 계산한다 — 서론이 던진 질문에 실제로 답하는 지표.

전역 상관과 R2 는 "점수가 pIC50 을 정량 예측하는가"를 묻는다. 그런데 도킹의 실무 용도는
"상위 N 건만 합성할 때 헛수고를 줄이는가"이다. 이 둘은 다른 질문이고 답도 다를 수 있다.
ROC-AUC / EF / BEDROC 없이 부정 결론을 내리면 자기 데이터를 읽지 않은 것이다.
"""
from __future__ import annotations
import argparse, json, math, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402
ROOT = HERE.parent


def auc_and_p(scores, labels, iters=20000, seed=42):
    """Mann-Whitney AUC. 점수는 낮을수록 좋으므로 부호를 뒤집어 '높을수록 활성'으로 맞춘다."""
    pos = [-s for s, l in zip(scores, labels) if l]
    neg = [-s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return None, None
    def _a(p, n):
        return sum(1.0 if a > b else 0.5 if a == b else 0.0
                   for a in p for b in n) / (len(p) * len(n))
    obs = _a(pos, neg)
    rng = random.Random(seed)
    allv = pos + neg; k = len(pos); c = 0
    for _ in range(iters):
        rng.shuffle(allv)
        if abs(_a(allv[:k], allv[k:]) - 0.5) >= abs(obs - 0.5) - 1e-12:
            c += 1
    return round(obs, 3), round((c + 1) / (iters + 1), 4)


def enrichment(scores, labels, frac):
    """상위 frac 비율에서의 농축 계수 EF. 1.0 이면 무작위와 같다."""
    n = len(scores); k = max(1, int(round(n * frac)))
    order = sorted(range(n), key=lambda i: scores[i])[:k]
    hits = sum(labels[i] for i in order)
    base = sum(labels) / n
    return round((hits / k) / base, 2) if base else None, hits, k


def bedroc(scores, labels, alpha=20.0):
    """초기 인식에 가중치를 주는 지표. alpha=20 은 상위 8% 를 주로 본다."""
    n = len(scores); a = sum(labels)
    if a == 0 or a == n:
        return None
    order = sorted(range(n), key=lambda i: scores[i])
    ranks = [i + 1 for i, idx in enumerate(order) if labels[idx]]
    ra = a / n
    s = sum(math.exp(-alpha * r / n) for r in ranks)
    rie = s / (a * (1 - math.exp(-alpha)) / (n * (math.exp(alpha / n) - 1)))
    num = ra * math.sinh(alpha / 2) / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra))
    return round(rie * ra * (math.sinh(alpha / 2)) / (math.cosh(alpha / 2)
                 - math.cosh(alpha / 2 - alpha * ra)) / num if num else rie, 3)


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


def spearman(x, y):
    rx, ry = ranks(x), ranks(y); n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return round(num / den, 3) if den else float("nan")


def fisher_ci(rho, n, z=1.96):
    """Fisher z 변환 신뢰구간. n 이 작아 넓게 나오며, 그 사실 자체가 결과다."""
    if abs(rho) >= 1 or n < 4:
        return None
    zr = 0.5 * math.log((1 + rho) / (1 - rho)); se = 1 / math.sqrt(n - 3)
    lo, hi = zr - z * se, zr + z * se
    return [round(math.tanh(lo), 3), round(math.tanh(hi), 3)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reg", default=str(ROOT / "sample_run" / "regression.json"))
    ap.add_argument("--out", default=str(ROOT / "sample_run" / "enrichment.json"))
    a = ap.parse_args()

    reg = json.loads(Path(a.reg).read_text())
    reg = reg.get("result", reg) if "rows" not in reg else reg
    rows = reg["rows"]; n = len(rows)
    p50 = [r["pIC50"] for r in rows]; stra = [r["stratum"] for r in rows]
    arms = {"selected_pose": [r["dock_score"] for r in rows],
            "top_pose": [r["top_pose_score"] for r in rows]}

    out = {"n": n, "arms": {}, "within_stratum": {}, "mcs_by_stratum": {}}
    for name, sc in arms.items():
        # strong(pIC50>=7) vs weak(<6) 이분. medium 은 경계라 제외한다.
        sw = [(s, st) for s, st in zip(sc, stra) if st in ("strong", "weak")]
        lab = [1 if st == "strong" else 0 for _, st in sw]
        au, ap_ = auc_and_p([s for s, _ in sw], lab)
        top10 = sorted(range(n), key=lambda i: sc[i])[:10]
        comp = {k: sum(1 for i in top10 if stra[i] == k) for k in ("strong", "medium", "weak")}
        ef10, hits, k = enrichment([s for s, _ in sw], lab, 10 / len(sw))
        rho = spearman(p50, sc)
        out["arms"][name] = {
            "spearman": rho, "spearman_ci95": fisher_ci(rho, n),
            "auc_strong_vs_weak": au, "auc_perm_p": ap_,
            "bedroc_a20": bedroc([s for s, _ in sw], lab),
            "ef_top10of20": ef10, "top10_strata": comp,
        }
    for st in ("strong", "medium", "weak"):
        idx = [i for i in range(n) if stra[i] == st]
        out["within_stratum"][st] = {
            "n": len(idx),
            "spearman_selected": spearman([p50[i] for i in idx],
                                          [arms["selected_pose"][i] for i in idx]),
            "pIC50_range": [round(min(p50[i] for i in idx), 2),
                            round(max(p50[i] for i in idx), 2)]}
        at = [r["mcs_atoms"] for r in rows if r["stratum"] == st and r.get("mcs_atoms")]
        rm = [r["mcs_rmsd_to_ref"] for r in rows
              if r["stratum"] == st and r.get("mcs_rmsd_to_ref") is not None]
        out["mcs_by_stratum"][st] = {
            "mean_mcs_atoms": round(sum(at) / len(at), 1) if at else None,
            "mean_mcs_rmsd": round(sum(rm) / len(rm), 2) if rm else None,
            "at_8atom_floor": sum(1 for x in at if x <= 8)}

    sel = out["arms"]["selected_pose"]; top = out["arms"]["top_pose"]
    checks = [
        ("두 자세 팔 모두 AUC 산출", sel["auc_strong_vs_weak"] is not None
                                    and top["auc_strong_vs_weak"] is not None),
        ("상위10 층 구성 합 = 10", sum(top["top10_strata"].values()) == 10),
        ("층별 MCS 통계 산출", all(v["mean_mcs_atoms"] for v in out["mcs_by_stratum"].values())),
    ]
    env = make_result(out, "regression.json 재계산", f"n={n}", checks,
                      notes=(f"AUC(strong vs weak) 선택 {sel['auc_strong_vs_weak']} / "
                             f"1위 {top['auc_strong_vs_weak']}. "
                             f"1위 자세 점수 상위10 층구성 {top['top10_strata']}. "
                             "정량 예측과 순위 분리는 서로 다른 질문이며 답도 다르다."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=2))
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
