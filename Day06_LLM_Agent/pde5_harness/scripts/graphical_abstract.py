#!/usr/bin/env python3
"""그래픽 초록 — 문제 → 설계 → 결과 3막 구조.

패널을 늘어놓는 것과 이야기를 만드는 것은 다르다. 독자가 그림 하나만 보고도
"무엇이 문제였고, 어떻게 설계했고, 무엇이 나왔는지"를 알 수 있어야 한다.
모든 수치는 sample_run/ 산출 파일에서 읽는다.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SR = ROOT / "sample_run"
FIG = SR / "report" / "figures_controlled"
OK = {"blue": "#0072B2", "orange": "#D55E00", "green": "#009E73", "pink": "#CC79A7",
      "yellow": "#E69F00", "sky": "#56B4E9", "grey": "#6E6E6E", "ink": "#1A1A1A"}
plt.rcParams.update({"font.size": 9, "figure.dpi": 300, "savefig.bbox": "tight",
                     "savefig.facecolor": "white"})
plt.rcParams["font.family"] = "Pretendard"      # 사용자 지정 서체
plt.rcParams["axes.unicode_minus"] = False         # Pretendard 의 마이너스 글리프 사용



def load(n):
    p = SR / n
    if not p.exists(): return None
    d = json.loads(p.read_text()); return d.get("result", d)


def panel_frame(ax, title, num, color):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.005, 0.005), 0.99, 0.99, transform=ax.transAxes,
                                boxstyle="round,pad=0.008,rounding_size=0.02",
                                fc="white", ec=color, lw=1.6, zorder=0))
    ax.add_patch(Circle((0.055, 0.925), 0.032, transform=ax.transAxes,
                        fc=color, ec="none", zorder=3))
    ax.text(0.055, 0.925, num, transform=ax.transAxes, ha="center", va="center",
            fontsize=10, color="white", fontweight="bold", zorder=4)
    ax.text(0.105, 0.925, title, transform=ax.transAxes, ha="left", va="center",
            fontsize=11.5, color=color, fontweight="bold")


def main() -> int:
    ds, dk, an = load("dataset_controlled.json"), load("docking_controlled.json"), \
                 load("analysis_controlled.json")
    tc, sw = load("terms_controlled.json"), load("exhaustiveness_sweep.json")
    prot, col = load("protonation_test.json"), load("collapse_diagnosis.json")
    if not (ds and dk and an):
        raise SystemExit("산출 파일 없음")
    T = an["arms"]["top_pose"]; rows = dk["rows"]; PV = T["pose_validity"]

    fig = plt.figure(figsize=(15, 5.6))
    gs = fig.add_gridspec(1, 3, wspace=0.10, left=0.015, right=0.985,
                          top=0.805, bottom=0.055)

    # ── ① 문제 ───────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0]); panel_frame(ax, "The problem", "1", OK["orange"])
    ax.text(.5, .845, "Benchmarks measure two things at once", ha="center",
            fontsize=10, color=OK["ink"], style="italic")
    ax.add_patch(FancyBboxPatch((.06, .60), .26, .13, boxstyle="round,pad=0.012",
                                fc=OK["sky"] + "33", ec=OK["sky"], lw=1.2))
    ax.text(.19, .665, "looks like the\nco-crystal ligand", ha="center", va="center", fontsize=8.5)
    ax.add_patch(FancyBboxPatch((.68, .60), .26, .13, boxstyle="round,pad=0.012",
                                fc=OK["blue"] + "33", ec=OK["blue"], lw=1.2))
    ax.text(.81, .665, "docking score\nlooks good", ha="center", va="center", fontsize=8.5)
    ax.add_patch(FancyBboxPatch((.37, .34), .26, .13, boxstyle="round,pad=0.012",
                                fc=OK["green"] + "33", ec=OK["green"], lw=1.2))
    ax.text(.50, .405, "high measured\npotency", ha="center", va="center", fontsize=8.5)
    for a, b, lab in (((.33, .665), (.67, .665), "fits the induced-fit pocket"),
                      ((.24, .59), (.42, .48), "same optimised series"),
                      ((.60, .48), (.79, .59), "?")):
        ax.add_patch(FancyArrowPatch(a, b, transform=ax.transAxes, arrowstyle="-|>",
                                     mutation_scale=13, lw=1.5,
                                     color=OK["orange"] if lab == "?" else OK["grey"],
                                     connectionstyle="arc3,rad=0.12"))
    ax.text(.50, .715, "fits the induced-fit pocket", ha="center", fontsize=7.5, color=OK["grey"])
    ax.text(.245, .50, "same series", ha="center", fontsize=7.5, color=OK["grey"], rotation=-32)
    ax.text(.735, .50, "what we want\nto measure", ha="center", fontsize=7.5,
            color=OK["orange"], rotation=32, fontweight="bold")
    ax.text(.5, .20, "A single correlation cannot tell these apart.\n"
                     "Is the score reading potency — or scaffold?",
            ha="center", va="center", fontsize=9.5, color=OK["orange"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.5", fc=OK["orange"] + "18",
                      ec=OK["orange"], lw=1.0))

    # ── ② 설계 ───────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[1]); panel_frame(ax, "The design", "2", OK["green"])
    ax.text(.5, .845, f"Orthogonalise the two axes  (n = {an['n']})", ha="center",
            fontsize=10, color=OK["ink"], style="italic")
    cnt = {(c["potency"], c["similarity"]): c["taken"] for c in ds["cells"]}
    P = ["strong", "medium", "weak"]; S = ["far", "mid", "near"]
    x0, y0, w, h = .20, .30, .19, .155
    mx = max(cnt.values())
    for i, pn in enumerate(P):
        for j, sn in enumerate(S):
            v = cnt.get((pn, sn), 0)
            ax.add_patch(Rectangle((x0 + j * w, y0 + (2 - i) * h), w * .93, h * .88,
                                   transform=ax.transAxes,
                                   fc=plt.cm.Greens(0.18 + 0.62 * v / mx), ec="white", lw=1.6))
            ax.text(x0 + j * w + w * .465, y0 + (2 - i) * h + h * .44, str(v),
                    transform=ax.transAxes, ha="center", va="center", fontsize=11,
                    color="white" if v > mx * .55 else OK["ink"], fontweight="bold")
        ax.text(x0 - .022, y0 + (2 - i) * h + h * .44, pn, transform=ax.transAxes,
                ha="right", va="center", fontsize=8.5)
    for j, sn in enumerate(S):
        ax.text(x0 + j * w + w * .465, y0 - .035, sn, transform=ax.transAxes,
                ha="center", va="top", fontsize=8.5)
    ax.text(x0 + 1.5 * w, y0 - .105, "Tanimoto to co-crystal ligand  →",
            transform=ax.transAxes, ha="center", fontsize=8.5, color=OK["grey"])
    ax.text(x0 - .085, y0 + 1.5 * h, "← measured potency", transform=ax.transAxes,
            ha="center", va="center", fontsize=8.5, color=OK["grey"], rotation=90)
    ax.text(.5, .155, f"confound  r = "
                      f"{ds['confound_pearson_r_potency_vs_similarity']:+.3f}"
                      f"   (potency-only design: +0.281)",
            ha="center", fontsize=9.5, color=OK["green"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.42", fc=OK["green"] + "18",
                      ec=OK["green"], lw=1.0))
    ax.text(.5, .065, "ChEMBL PDE5A, all "
                      f"{ds['records_scanned']} records → median per compound → 9-cell grid",
            ha="center", fontsize=7.6, color=OK["grey"])

    # ── ③ 결과 ───────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[2]); panel_frame(ax, "The result", "3", OK["blue"])
    ax.text(.5, .845, "No signal — and we know what it is not", ha="center",
            fontsize=10, color=OK["ink"], style="italic")
    sub = fig.add_axes([0.700, 0.455, 0.132, 0.225])
    SIMC = {"near": OK["pink"], "mid": OK["green"], "far": OK["sky"]}
    for b in ("far", "mid", "near"):
        xs = [r["top_pose_score"] for r in rows if r["similarity_bin"] == b]
        ys = [r["pchembl_value"] for r in rows if r["similarity_bin"] == b]
        sub.scatter(xs, ys, s=7, c=SIMC[b], edgecolor="none", alpha=.75)
    sub.set_xlabel("docking score", fontsize=7); sub.set_ylabel("pChEMBL", fontsize=7)
    sub.tick_params(labelsize=6); sub.spines[["top", "right"]].set_visible(False)
    ax.text(.585, .445, f"ρ = {T['spearman']:+.3f}\n95% CI "
                        f"[{T['ci95'][0]:+.2f}, {T['ci95'][1]:+.2f}]\np = {T['perm_p']}",
            transform=ax.transAxes, ha="center", va="top", fontsize=9,
            fontweight="bold", color=OK["orange"])
    ax.text(.855, .755, "n.s.", transform=ax.transAxes, ha="center", fontsize=13,
            color=OK["orange"], fontweight="bold")
    elim = [("search depth 4→128", f"C2 flat, {sw['diagnosis']['c2_range_angstrom']} Å" if sw else "—"),
            ("protonation (2×2)", "score identical" if prot else "—"),
            ("scaffold confound", f"old-data attenuation "
                                  f"{col['old_design']['attenuation_from_scaffold_control']:.1%}" if col else "—")]
    ax.text(.07, .335, "Ruled out by experiment", transform=ax.transAxes, fontsize=9,
            fontweight="bold", color=OK["ink"])
    for k, (a_, b_) in enumerate(elim):
        yy = .275 - k * .062
        ax.text(.085, yy, "✕", transform=ax.transAxes, fontsize=9, color=OK["green"],
                fontweight="bold")
        ax.text(.135, yy, a_, transform=ax.transAxes, fontsize=8.2, va="center")
        ax.text(.60, yy, b_, transform=ax.transAxes, fontsize=7.6, va="center",
                color=OK["grey"])
    ax.text(.5, .065, f"Ceiling: {1 - PV['all']['frac_under_threshold']:.0%} of scored poses "
                      f"lie outside 2 Å of the crystal frame",
            transform=ax.transAxes, ha="center", fontsize=8.6, color=OK["orange"],
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.36", fc=OK["orange"] + "15",
                      ec=OK["orange"], lw=0.9))

    fig.suptitle("Orthogonalising scaffold similarity and potency in a PDE5A docking benchmark",
                 fontsize=15, fontweight="bold", y=0.975, color=OK["ink"])
    fig.text(0.5, 0.905, "A docking score that ranked potency at n = 30 did not replicate at "
                         "n = 163 once scaffold and potency were decoupled",
             ha="center", fontsize=10.5, color=OK["grey"])
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"fig01_graphical_abstract.{ext}")
    plt.close(fig)
    p = FIG / "fig01_graphical_abstract.png"
    print(f"  {p.name}  ({p.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
