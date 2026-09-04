#!/usr/bin/env python3
"""항 기여도 시각화 — 네 지표가 서로 다른 답을 준다는 사실 자체를 보여준다.

Fig 13  네 지표 비교 (표준화 계수 / 단변량 r / 순열 Q2 하락 / VIF)
Fig 14  항끼리 상관 히트맵 — 공선성 때문에 개별 계수를 못 읽는 이유
Fig 15  항별 산점도 격자 — 눈으로 봐도 신호가 없다는 것
Fig 16  선별 지표 요약 (ROC 곡선 + 상위10 층 구성 + 농축)
Fig 17  커스텀 채점 비교 (독립 시험 세트, 부트스트랩 CI 포함)

라벨은 영문이다 — DejaVu Sans 에 한글 글리프가 없어 네모로 깨진다. 본문 서술이 한국어다.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SR = ROOT / "sample_run"
FIG = SR / "report" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
# Okabe-Ito — 색각 이상에서도 구분된다
OK = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight"})
plt.rcParams["font.family"] = "Pretendard"      # 사용자 지정 서체
plt.rcParams["axes.unicode_minus"] = False         # Pretendard 의 마이너스 글리프 사용

SHORT = ["gauss1", "gauss2", "repulsion", "hydrophobic", "hbond"]


def load(n):
    p = SR / n
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return d.get("result", d)


def fig13(ti):
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.4))
    x = np.arange(len(SHORT))
    panels = [
        ("Standardized coefficient", [ti["standardized_coefficients"][t] for t in SHORT],
         "Comparable across terms\n(predictors z-scored)"),
        ("Univariate Pearson r", [ti["univariate"][t]["pearson_r"] for t in SHORT],
         "Each term alone vs pIC50"),
        ("Permutation Q² drop", [ti["permutation_importance"][t]["mean_Q2_drop"] for t in SHORT],
         "Positive = term carries signal"),
        ("VIF (collinearity)", [ti["vif"][t] for t in SHORT],
         "VIF > 5 = coefficient unstable"),
    ]
    for k, (title, vals, sub) in enumerate(panels):
        cols = [OK[1] if v < 0 else OK[0] for v in vals] if k < 3 else \
               [OK[1] if v and v > 5 else OK[2] for v in vals]
        ax[k].bar(x, vals, color=cols, edgecolor="black", linewidth=.5)
        ax[k].set_xticks(x); ax[k].set_xticklabels(SHORT, rotation=40, ha="right")
        ax[k].set_title(title, fontsize=9.5, fontweight="bold")
        ax[k].text(.5, -.42, sub, transform=ax[k].transAxes, ha="center",
                   fontsize=7.5, color="#555")
        if k < 3:
            ax[k].axhline(0, color="black", lw=.8)
        else:
            ax[k].axhline(5, color=OK[1], lw=.9, ls="--")
    ax[2].text(.5, .92, "ALL NEGATIVE:\nshuffling any term improves Q²",
               transform=ax[2].transAxes, ha="center", va="top", fontsize=7.5,
               color=OK[1], fontweight="bold",
               bbox=dict(fc="white", ec=OK[1], lw=.8, alpha=.9, pad=2))
    fig.suptitle("Which scoring term matters? Four metrics disagree", fontweight="bold", y=1.06)
    p = FIG / "fig13_term_importance.png"; fig.savefig(p); plt.close(fig); return p


def fig14(ti):
    M = np.array([[ti["term_correlation"][a][b] for b in SHORT] for a in SHORT])
    fig, ax = plt.subplots(figsize=(4.6, 4.1))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(SHORT))); ax.set_xticklabels(SHORT, rotation=40, ha="right")
    ax.set_yticks(range(len(SHORT))); ax.set_yticklabels(SHORT)
    for i in range(len(SHORT)):
        for j in range(len(SHORT)):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7.5,
                    color="white" if abs(M[i, j]) > .55 else "black")
    fig.colorbar(im, ax=ax, shrink=.8, label="Pearson r")
    ax.set_title("Scoring terms are correlated with each other\n"
                 "(so individual coefficients cannot be read alone)",
                 fontsize=9, fontweight="bold")
    p = FIG / "fig14_term_correlation.png"; fig.savefig(p); plt.close(fig); return p


def fig15(reg, ti):
    rows = reg["rows"]; y = [r["pIC50"] for r in rows]
    strat = [r["stratum"] for r in rows]
    cmap = {"strong": OK[0], "medium": OK[4], "weak": OK[1]}
    fig, ax = plt.subplots(1, 5, figsize=(14, 2.9), sharey=True)
    for k, t in enumerate(SHORT):
        v = [r[t] for r in rows]
        for s in ("strong", "medium", "weak"):
            xs = [a for a, b in zip(v, strat) if b == s]
            ys = [a for a, b in zip(y, strat) if b == s]
            ax[k].scatter(xs, ys, s=22, c=cmap[s], edgecolor="black", linewidth=.3,
                          label=s if k == 0 else None, alpha=.85)
        r = ti["univariate"][t]["pearson_r"]
        ax[k].set_title(f"{t}\nr = {r:+.3f}", fontsize=9,
                        color=OK[1] if abs(r) < .25 else "black",
                        fontweight="bold" if abs(r) >= .25 else "normal")
        ax[k].set_xlabel(f"{t} term value", fontsize=8)
    ax[0].set_ylabel("measured pIC50")
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.suptitle("No single term separates potency (n=30)", fontweight="bold", y=1.09)
    p = FIG / "fig15_term_scatter.png"; fig.savefig(p); plt.close(fig); return p


def fig16(enr, reg):
    rows = reg["rows"]
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
    # A. ROC
    for nm, key, col in (("selected pose", "dock_score", OK[0]),
                         ("top-scored pose", "top_pose_score", OK[1])):
        sw = [(r[key], 1 if r["stratum"] == "strong" else 0)
              for r in rows if r["stratum"] in ("strong", "weak")]
        sw.sort(key=lambda t: t[0])
        P = sum(l for _, l in sw); N = len(sw) - P
        tp = fp = 0; xs, ys = [0], [0]
        for _, l in sw:
            tp += l; fp += (1 - l); xs.append(fp / N); ys.append(tp / P)
        au = enr["arms"]["selected_pose" if key == "dock_score" else "top_pose"]["auc_strong_vs_weak"]
        ax[0].plot(xs, ys, color=col, lw=2, label=f"{nm} (AUC {au})")
    ax[0].plot([0, 1], [0, 1], ls="--", color="grey", lw=1, label="random")
    ax[0].set_xlabel("false positive rate"); ax[0].set_ylabel("true positive rate")
    ax[0].set_title("A. Strong vs weak discrimination", fontsize=9.5, fontweight="bold")
    ax[0].legend(frameon=False, fontsize=7.5, loc="lower right")
    # B. 상위10 층 구성
    lbl = ["strong", "medium", "weak"]; w = .35
    for i, (nm, arm) in enumerate((("selected", "selected_pose"), ("top-scored", "top_pose"))):
        c = enr["arms"][arm]["top10_strata"]
        ax[1].bar(np.arange(3) + (i - .5) * w, [c[l] for l in lbl], w,
                  color=OK[i], edgecolor="black", lw=.5, label=f"{nm} pose")
    ax[1].axhline(10 / 3, color="grey", ls="--", lw=1)
    ax[1].text(2.35, 10 / 3 + .12, "random\nexpectation", fontsize=7, color="grey", ha="right")
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels(lbl)
    ax[1].set_ylabel("compounds in score top-10")
    ax[1].set_title("B. Composition of the top-10 by score", fontsize=9.5, fontweight="bold")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].annotate("zero weak\ncompounds", xy=(2.18, .12), xytext=(1.55, 3.4),
                   fontsize=7.5, color=OK[1], fontweight="bold",
                   arrowprops=dict(arrowstyle="->", color=OK[1], lw=1))
    # C. 두 질문의 답
    q = ["Quantitative\nprediction\n(Q²)", "Rank\nseparation\n(AUC)"]
    val = [reg["models"]["single_vina_score"]["Q2_loo"],
           enr["arms"]["top_pose"]["auc_strong_vs_weak"]]
    thr = [0.3, 0.7]
    for i, (vv, tt) in enumerate(zip(val, thr)):
        ok = vv >= tt
        ax[2].bar(i, vv, .5, color=OK[2] if ok else OK[1], edgecolor="black", lw=.6)
        ax[2].plot([i - .3, i + .3], [tt, tt], color="black", lw=1.4, ls="--")
        ax[2].text(i, max(vv, tt) + .05, "PASS" if ok else "FAIL", ha="center",
                   fontsize=9, fontweight="bold", color=OK[2] if ok else OK[1])
    ax[2].set_xticks(range(2)); ax[2].set_xticklabels(q, fontsize=8)
    ax[2].axhline(0, color="black", lw=.8); ax[2].set_ylim(-.15, 1.05)
    ax[2].set_ylabel("metric value")
    ax[2].set_title("C. Two questions, two answers", fontsize=9.5, fontweight="bold")
    ax[2].text(.5, -.32, "dashed line = pre-stated success threshold",
               transform=ax[2].transAxes, ha="center", fontsize=7.5, color="#555")
    fig.suptitle("Docking triages but does not quantify", fontweight="bold", y=1.04)
    p = FIG / "fig16_triage_summary.png"; fig.savefig(p); plt.close(fig); return p


def fig17(cust):
    if not cust:
        return None
    d = cust; fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    nm = ["default Vina", "fitted custom"]
    v = [d["arm_default"]["spearman"], d["arm_custom"]["spearman"]]
    ax[0].bar(range(2), v, .5, color=[OK[0], OK[3]], edgecolor="black", lw=.6)
    ax[0].set_xticks(range(2)); ax[0].set_xticklabels(nm)
    ax[0].set_ylabel("held-out Spearman ρ"); ax[0].axhline(0, color="black", lw=.8)
    ax[0].set_title(f"A. Independent test set (n={d['test_n']})", fontsize=9.5, fontweight="bold")
    for i, vv in enumerate(v):
        ax[0].text(i, vv - .04, f"{vv:+.3f}", ha="center", va="top", fontsize=9)
    ci = d["difference"]["bootstrap_ci95"]; dd = d["difference"]["delta_spearman"]
    ax[1].errorbar([dd], [0], xerr=[[dd - ci[0]], [ci[1] - dd]], fmt="o", color=OK[3],
                   capsize=6, markersize=9, lw=2)
    ax[1].axvline(0, color=OK[1], lw=1.6, ls="--")
    ax[1].set_yticks([]); ax[1].set_xlabel("Δ Spearman (custom − default)")
    ax[1].set_xlim(min(ci) - .12, max(ci) + .12)
    ax[1].set_title("B. The difference is not distinguishable from zero",
                    fontsize=9.5, fontweight="bold")
    ax[1].text(.5, .74, f"Δ = {dd:+.3f}   95% CI [{ci[0]:+.2f}, {ci[1]:+.2f}]\n"
                        f"{d['difference']['frac_resamples_favoring_custom']:.0%} of resamples "
                        f"favour custom — a coin flip",
               transform=ax[1].transAxes, ha="center", fontsize=8,
               bbox=dict(fc="white", ec=OK[1], lw=.8, alpha=.95, pad=4))
    fig.suptitle("Re-docking with fitted weights: the number moved, the evidence did not",
                 fontweight="bold", y=1.04)
    p = FIG / "fig17_custom_scoring.png"; fig.savefig(p); plt.close(fig); return p


def main() -> int:
    ti, reg, enr, cust = load("term_importance.json"), load("regression.json"), \
                         load("enrichment.json"), load("custom_score.json")
    if not (ti and reg and enr):
        raise SystemExit("필요한 산출 파일이 없다 — 그림을 만들지 않는다.")
    made = [f for f in (fig13(ti), fig14(ti), fig15(reg, ti), fig16(enr, reg), fig17(cust)) if f]
    for p in made:
        print(f"  {p.name}  ({p.stat().st_size // 1024} KB)")
    print(f"그림 {len(made)}장")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
