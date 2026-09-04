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


def _problem_panel(ax, big=False):
    """문제 도식 — 그래픽 초록의 A 패널이자 발표 2번 슬라이드의 전면 그림.

    단독 렌더(big=True)에서는 상자를 넓히고 화살표 라벨을 상자 밖으로 빼서 겹침을 없앤다.
    같은 좌표를 두 크기에 그대로 쓰면 라벨이 상자를 침범한다.
    """
    k = 1.95 if big else 1.0
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    bw, bh = (.28, .15) if big else (.26, .13)
    ty, by = (.72, .40) if big else (.665, .405)

    def node(cx, cy, txt, fc, ec):
        ax.add_patch(FancyBboxPatch((cx - bw / 2, cy - bh / 2), bw, bh,
                                    boxstyle="round,pad=0.014,rounding_size=0.02",
                                    fc=fc + "30", ec=ec, lw=1.3 * k, zorder=3))
        ax.text(cx, cy, txt, ha="center", va="center", fontsize=8.6 * k, zorder=4)

    lx, rx = (.245, .755) if big else (.20, .80)
    node(lx, ty, "공결정 리간드를\n닮은 화합물", OK["sky"], OK["sky"])
    node(rx, ty, "도킹 점수가\n좋게 나온다", OK["blue"], OK["blue"])
    node(.50, by, "실측 역가가\n높다", OK["green"], OK["green"])

    def arrow(a, b, colr, rad=0.0, lw=1.5):
        ax.add_patch(FancyArrowPatch(a, b, transform=ax.transAxes, arrowstyle="-|>",
                                     mutation_scale=13 * k, lw=lw * k, color=colr,
                                     connectionstyle=f"arc3,rad={rad}", zorder=2))

    # 가로 화살표는 상자 사이 직선, 라벨은 그 위쪽에 둔다
    arrow((lx + bw / 2 + .01, ty), (rx - bw / 2 - .01, ty), OK["grey"])
    ax.text(.50, ty + bh / 2 + .045, "유도적합된 포켓에 잘 맞아서",
            ha="center", fontsize=7.6 * k, color=OK["grey"])
    # 대각 화살표는 라벨을 바깥쪽으로
    arrow((lx, ty - bh / 2 - .01), (.50 - bw / 2 - .005, by + bh / 2 - .01),
          OK["grey"], rad=0.14)
    ax.text(lx - bw / 2 - .015, (ty + by) / 2, "같은 계열의\n최적화 산물이라",
            ha="right", va="center", fontsize=7.6 * k, color=OK["grey"])
    arrow((.50 + bw / 2 + .005, by + bh / 2 - .01), (rx, ty - bh / 2 - .01),
          OK["orange"], rad=0.14, lw=2.0)
    ax.text(rx + bw / 2 + .015, (ty + by) / 2, "우리가\n재려던 것", ha="left",
            va="center", fontsize=7.8 * k, color=OK["orange"], fontweight="bold")

    ax.text(.5, .085 if big else .105, "상관 하나로는 이 둘을 구분할 수 없다.\n"
                      "점수는 역가를 읽은 것인가, 골격 닮음을 읽은 것인가?",
            ha="center", va="center", fontsize=9.2 * k, color=OK["orange"],
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.55", fc=OK["orange"] + "18",
                      ec=OK["orange"], lw=1.1 * k))


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
    ax = fig.add_subplot(gs[0]); panel_frame(ax, "The problem  문제", "A", OK["orange"])
    ax.text(.5, .855, "Benchmarks measure two things at once", ha="center",
            fontsize=10, color=OK["ink"], style="italic")
    ax.text(.5, .805, "벤치마크는 두 가지를 동시에 재고 있다", ha="center",
            fontsize=9, color=OK["grey"])
    _problem_panel(ax)

    # ── ② 설계 ───────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[1]); panel_frame(ax, "The design  설계", "B", OK["green"])
    ax.text(.5, .855, f"Orthogonalise the two axes  (n = {an['n']})", ha="center",
            fontsize=10, color=OK["ink"], style="italic")
    ax.text(.5, .805, "역가와 골격 유사도를 직교화한다", ha="center",
            fontsize=9, color=OK["grey"])
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
    ax = fig.add_subplot(gs[2]); panel_frame(ax, "The result  결과", "C", OK["blue"])
    ax.text(.5, .855, "No signal — and we know what it is not", ha="center",
            fontsize=10, color=OK["ink"], style="italic")
    ax.text(.5, .805, "신호 없음 — 다만 귀무를 채택한 것은 아니다", ha="center",
            fontsize=9, color=OK["grey"])
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
    ax.text(.07, .385, "Ruled out by experiment", transform=ax.transAxes, fontsize=9,
            fontweight="bold", color=OK["ink"])
    for k, (a_, b_) in enumerate(elim):
        yy = .325 - k * .062
        ax.text(.085, yy, "✕", transform=ax.transAxes, fontsize=9, color=OK["green"],
                fontweight="bold")
        ax.text(.135, yy, a_, transform=ax.transAxes, fontsize=8.2, va="center")
        ax.text(.60, yy, b_, transform=ax.transAxes, fontsize=7.6, va="center",
                color=OK["grey"])
    ax.text(.5, .115, f"Ceiling: {1 - PV['all']['frac_under_threshold']:.0%} of scored poses "
                      f"lie outside 2 Å of the crystal frame",
            transform=ax.transAxes, ha="center", fontsize=8.6, color=OK["orange"],
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.32", fc=OK["orange"] + "15",
                      ec=OK["orange"], lw=0.9))
    ax.text(.5, .042, f"Not an accepted null: the CI still admits |ρ| up to "
                      f"{abs(T['ci95'][0]):.2f}", transform=ax.transAxes, ha="center",
            fontsize=8, color=OK["grey"], style="italic")

    fig.suptitle("Orthogonalising scaffold similarity and potency in a PDE5A docking benchmark",
                 fontsize=15, fontweight="bold", y=0.975, color=OK["ink"])
    fig.text(0.5, 0.905, "A docking score that ranked potency at n = 30 did not replicate at "
                         "n = 163 once scaffold and potency were decoupled",
             ha="center", fontsize=10.5, color=OK["grey"])
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"fig01_graphical_abstract.{ext}")
    plt.close(fig)

    # 발표 슬라이드 2 용 — 문제 패널만 단독으로. 전체 초록을 2번에 띄우면 뒤 여섯 장을
    # 미리 다 보여주는 데다, 16:9 투사에서 본문이 6-7 pt 라 읽히지도 않는다.
    f2 = plt.figure(figsize=(10.0, 5.4))
    ax = f2.add_axes([0.02, 0.02, 0.96, 0.96])
    _problem_panel(ax, big=True)
    for ext in ("png", "svg"):
        f2.savefig(FIG / f"fig00_problem.{ext}")
    plt.close(f2)
    print(f"  fig00_problem.png  ({(FIG / 'fig00_problem.png').stat().st_size // 1024} KB)")
    p = FIG / "fig01_graphical_abstract.png"
    print(f"  {p.name}  ({p.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
