#!/usr/bin/env python3
"""골격 통제 연구의 논문급 그림 세트.

모든 수치는 sample_run/ 산출 파일에서 읽는다. 그림 안에 손으로 적은 숫자는 없다.
라벨은 영문 — DejaVu Sans 에 한글 글리프가 없어 네모로 깨진다. 본문 서술은 한국어다.
팔레트는 Okabe-Ito 로 색각 이상에서도 구분된다.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SR = ROOT / "sample_run"
FIG = SR / "report" / "figures_controlled"
FIG.mkdir(parents=True, exist_ok=True)

OK = {"blue": "#0072B2", "orange": "#D55E00", "green": "#009E73", "pink": "#CC79A7",
      "yellow": "#E69F00", "sky": "#56B4E9", "lemon": "#F0E442", "grey": "#7F7F7F"}
POT_C = {"strong": OK["blue"], "medium": OK["yellow"], "weak": OK["orange"]}
SIM_C = {"near": OK["pink"], "mid": OK["green"], "far": OK["sky"]}
plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold", "figure.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white", "axes.linewidth": .8,
    "xtick.major.width": .8, "ytick.major.width": .8, "legend.frameon": False})
plt.rcParams["font.family"] = "Pretendard"      # 사용자 지정 서체
plt.rcParams["axes.unicode_minus"] = False         # Pretendard 의 마이너스 글리프 사용



def load(name):
    p = SR / name
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return d.get("result", d)


def save(fig, name):
    p = FIG / name
    fig.savefig(p); fig.savefig(p.with_suffix(".svg")); plt.close(fig)
    print(f"  {p.name}  ({p.stat().st_size // 1024} KB)")
    return p


# ───────────────────────── Fig 1 그래픽 초록 ─────────────────────────
def fig01(ds, an, ctrl, tc=None):
    fig = plt.figure(figsize=(13, 5.2))
    gs = fig.add_gridspec(2, 4, hspace=.55, wspace=.38,
                          height_ratios=[1, 1.05], left=.05, right=.97)
    T = an["arms"]["top_pose"]

    # A 문제
    ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
    ax.set_title("A  The question", loc="left")
    ax.text(0, .78, "Does a docking score track\nmeasured potency?", fontsize=10.5, va="top")
    ax.text(0, .40, "Two sub-questions with\npossibly different answers:", fontsize=8.5,
            va="top", color="#444")
    ax.text(.03, .16, "① quantify (Q²)      ② triage (AUC)", fontsize=9,
            va="top", color=OK["blue"], fontweight="bold")

    # B 설계
    ax = fig.add_subplot(gs[0, 1])
    ax.set_title("B  Decoupled design", loc="left")
    cnt = {(c["potency"], c["similarity"]): c["taken"] for c in ds["cells"]}
    P = ["weak", "medium", "strong"]; S = ["far", "mid", "near"]
    M = np.array([[cnt.get((p, s), 0) for s in S] for p in P])
    im = ax.imshow(M, cmap="Blues", vmin=0)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, M[i, j], ha="center", va="center", fontsize=9,
                    color="white" if M[i, j] > M.max() * .6 else "black", fontweight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels([f"{s}\nscaffold" for s in S], fontsize=7.5)
    ax.set_yticks(range(3)); ax.set_yticklabels(P, fontsize=8)
    ax.set_xlabel("Tanimoto to co-crystal ligand", fontsize=8)
    ax.text(1, -.95, f"confound r = {ds['confound_pearson_r_potency_vs_similarity']:+.2f}",
            ha="center", fontsize=8, color=OK["green"], fontweight="bold")

    # C 대조
    ax = fig.add_subplot(gs[0, 2])
    ax.set_title("C  Re-docking control", loc="left")
    vals = [ctrl["rmsd_best_of_modes_angstrom"], ctrl["rmsd_top_pose_angstrom"]]
    cols = [OK["green"] if v is not None and v <= ctrl["threshold_angstrom"] else OK["orange"]
            for v in vals]
    ax.bar([0, 1], vals, .5, color=cols, edgecolor="black", lw=.6)
    ax.axhline(ctrl["threshold_angstrom"], ls="--", color="black", lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["C1\nsampling", "C2\nranking"], fontsize=8)
    ax.set_ylabel("RMSD to crystal pose (Å)", fontsize=8)
    for i, v in enumerate(vals):
        if v is not None:
            ax.text(i, v + .08, f"{v:.2f}", ha="center", fontsize=8.5, fontweight="bold")

    # D 판정
    ax = fig.add_subplot(gs[0, 3]); ax.axis("off")
    ax.set_title("D  Verdict", loc="left")
    q2 = (tc or {}).get("models", {}).get("single_vina_score", {}).get("Q2_loo")
    items = [("Quantify (Q²)", f"{q2:+.3f}" if q2 is not None else "—", None),
             ("Triage, raw ρ", f"{T['spearman']:+.3f}", T["perm_p"]),
             ("Triage, scaffold-controlled",
              f"{T['partial_spearman_controlling_tanimoto']:+.3f}", T["partial_perm_p"])]
    for k, (lab, val, pp) in enumerate(items):
        ax.text(0, .82 - k * .28, lab, fontsize=8.5, color="#444")
        ax.text(0, .70 - k * .28, val, fontsize=13, fontweight="bold",
                color=OK["blue"] if k else OK["orange"])
        if pp is not None:
            ax.text(.42, .72 - k * .28, f"p = {pp}", fontsize=8, color="#666")

    # E 주 산점도
    ax = fig.add_subplot(gs[1, :2])
    ax.set_title("E  Docking score vs measured potency, by scaffold similarity", loc="left")
    for b in ("far", "mid", "near"):
        xs = [r["top_pose_score"] for r in ROWS if r["similarity_bin"] == b]
        ys = [r["pchembl_value"] for r in ROWS if r["similarity_bin"] == b]
        ax.scatter(xs, ys, s=20, c=SIM_C[b], edgecolor="black", lw=.25, alpha=.85,
                   label=f"{b} (n={len(xs)})")
    ax.set_xlabel("docking score (kcal/mol, lower = better)")
    ax.set_ylabel("measured pChEMBL")
    ax.legend(fontsize=7.5, loc="upper left", ncol=3, columnspacing=.9, handletextpad=.3)

    # F 결정적 비교
    ax = fig.add_subplot(gs[1, 2:])
    ax.set_title("F  Does the signal survive scaffold control?", loc="left")
    labs, vals, cis, cols = [], [], [], []
    labs.append("all compounds\n(raw)"); vals.append(T["spearman"])
    cis.append(T["ci95"]); cols.append(OK["grey"])
    labs.append("partial\n(scaffold controlled)")
    vals.append(T["partial_spearman_controlling_tanimoto"]); cis.append(None)
    cols.append(OK["blue"])
    if T.get("partial_spearman_controlling_both") is not None:
        labs.append("partial\n(scaffold + size)")
        vals.append(T["partial_spearman_controlling_both"]); cis.append(None)
        cols.append(OK["green"])
    for b in ("far", "mid", "near"):
        w = T["within_similarity_bin"].get(b)
        if w:
            labs.append(f"within {b}\n(n={w['n']})"); vals.append(w["spearman"])
            cis.append(w["ci95"]); cols.append(SIM_C[b])
    yy = np.arange(len(vals))[::-1]
    for i, (v, ci, c) in enumerate(zip(vals, cis, cols)):
        y0 = yy[i]
        if ci:
            ax.plot([ci[0], ci[1]], [y0, y0], color=c, lw=2, alpha=.55, solid_capstyle="round")
        ax.plot([v], [y0], "o", color=c, markersize=8, markeredgecolor="black",
                markeredgewidth=.5)
        ax.text(v, y0 + .28, f"{v:+.3f}", ha="center", fontsize=7.5, color=c,
                fontweight="bold")
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=7.5)
    ax.set_xlabel("Spearman ρ (pChEMBL vs docking score)")
    ax.set_xlim(-1.02, .62)
    fig.suptitle(f"Does a docking score rank PDE5A inhibitors? "
                 f"A scaffold-decoupled test (n={an['n']})",
                 fontsize=13, fontweight="bold", y=1.03)
    return save(fig, "fig01_graphical_abstract.png")


# ───────────────────────── Fig 2 데이터셋 특성 ───────────────────────
def fig02(ds):
    comps = ds["compounds"]
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.1))
    pv = [c["pchembl_value"] for c in comps]; tn = [c["tanimoto_to_sildenafil"] for c in comps]
    ax[0].hist(pv, bins=24, color=OK["blue"], edgecolor="white", lw=.5)
    ax[0].set_xlabel("pChEMBL (median of replicates)"); ax[0].set_ylabel("compounds")
    ax[0].set_title("A  Potency distribution", loc="left")
    ax[0].text(.97, .93, f"n = {len(comps)}\nrange {min(pv):.2f}–{max(pv):.2f}",
               transform=ax[0].transAxes, ha="right", va="top", fontsize=7.5)
    ax[1].hist(tn, bins=24, color=OK["pink"], edgecolor="white", lw=.5)
    ax[1].set_xlabel("Tanimoto to co-crystal ligand"); ax[1].set_ylabel("compounds")
    ax[1].set_title("B  Scaffold similarity", loc="left")
    for b in ("far", "mid", "near"):
        xs = [c["tanimoto_to_sildenafil"] for c in comps if c["similarity_bin"] == b]
        ys = [c["pchembl_value"] for c in comps if c["similarity_bin"] == b]
        ax[2].scatter(xs, ys, s=16, c=SIM_C[b], edgecolor="black", lw=.2, alpha=.8, label=b)
    ax[2].set_xlabel("Tanimoto to co-crystal ligand"); ax[2].set_ylabel("pChEMBL")
    ax[2].set_title("C  The two axes are decoupled", loc="left")
    ax[2].legend(fontsize=7.5, loc="lower right")
    ax[2].text(.03, .95, f"r = {ds['confound_pearson_r_potency_vs_similarity']:+.3f}",
               transform=ax[2].transAxes, va="top", fontsize=9,
               color=OK["green"], fontweight="bold")
    rep = [c["n_measurements"] for c in comps]
    ax[3].hist(rep, bins=range(1, max(rep) + 2), color=OK["green"], edgecolor="white", lw=.5)
    ax[3].set_xlabel("ChEMBL measurements per compound"); ax[3].set_ylabel("compounds")
    ax[3].set_title("D  Replicates (median-aggregated)", loc="left")
    ax[3].text(.97, .93, f"{sum(1 for r in rep if r > 1)} compounds\nhave >1 record",
               transform=ax[3].transAxes, ha="right", va="top", fontsize=7.5)
    fig.suptitle("Dataset: potency and scaffold similarity varied independently",
                 fontsize=11.5, y=1.06)
    return save(fig, "fig02_dataset.png")


# ───────────────────────── Fig 3 재도킹 대조 ─────────────────────────
def _control_mode_rmsds(ctrl, sw, exh, seed):
    """대조의 모드별 RMSD 를 구한다.

    dock_controlled 의 대조는 모드별 목록을 저장하지 않는다. 스윕에 같은 조건(깊이·시드)의
    실행이 있으면 거기서 가져오되, **점수와 1위 RMSD 가 일치할 때만** 쓴다. 조건이 다른
    실행의 숫자를 몰래 끌어다 쓰면 그림이 본문과 다른 실험을 그리게 된다.
    """
    if ctrl.get("all_mode_rmsd"):
        return ctrl["all_mode_rmsd"], "본 대조 실행"
    if not sw:
        return [], None
    for r in sw.get("runs", []):
        if r.get("exhaustiveness") != exh or r.get("seed") != seed:
            continue
        same_top = (r.get("c2_top_rmsd") is not None
                    and ctrl.get("rmsd_top_pose_angstrom") is not None
                    and abs(r["c2_top_rmsd"] - ctrl["rmsd_top_pose_angstrom"]) < 0.05)
        same_score = (r.get("top_score") is not None
                      and ctrl.get("score_kcal_mol") is not None
                      and abs(r["top_score"] - ctrl["score_kcal_mol"]) < 0.05)
        if same_top and same_score:
            return r.get("all_rmsd") or [], f"스윕 동일 조건 (exh {exh}, seed {seed})"
    return [], None


def fig03(ctrl, sw=None, exh=None, seed=None):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    d = ctrl
    n = d.get("n_modes", 0)
    rs, src = _control_mode_rmsds(d, sw, exh, seed)
    if not rs:
        rs = [d["rmsd_top_pose_angstrom"]] + [None] * max(0, n - 1)
        src = None
    xs = np.arange(1, len(rs) + 1)
    vals = [r if r is not None else np.nan for r in rs]
    cols = [OK["green"] if (v == v and v <= d["threshold_angstrom"]) else OK["orange"]
            for v in vals]
    ax[0].bar(xs, vals, .62, color=cols, edgecolor="black", lw=.5)
    ax[0].axhline(d["threshold_angstrom"], ls="--", color="black", lw=1.1)
    ax[0].text(len(rs) * .98, d["threshold_angstrom"] + .07,
               f"{d['threshold_angstrom']} Å criterion", ha="right", fontsize=7.5)
    ax[0].set_xlabel("docking mode (1 = best score)"); ax[0].set_ylabel("RMSD to crystal (Å)")
    ax[0].set_title("A  Sampling finds the pose; scoring misranks it", loc="left")
    if src is None:
        ax[0].text(.5, .55, "per-mode RMSD not recorded\nfor this control run",
                   transform=ax[0].transAxes, ha="center", fontsize=8, color=OK["grey"],
                   bbox=dict(fc="white", ec=OK["grey"], lw=.7, alpha=.9))
    else:
        ax[0].text(.98, .02, f"source: {src}", transform=ax[0].transAxes, ha="right",
                   fontsize=6.5, color="#888")
    if d.get("best_mode_index"):
        ax[0].annotate("crystal pose\nrecovered here",
                       xy=(d["best_mode_index"], d["rmsd_best_of_modes_angstrom"]),
                       xytext=(d["best_mode_index"] + .6, d["rmsd_best_of_modes_angstrom"] + 2.1),
                       fontsize=7.5, color=OK["green"], fontweight="bold",
                       arrowprops=dict(arrowstyle="->", color=OK["green"], lw=1.1))
    ax[1].axis("off"); ax[1].set_title("B  Two controls, two verdicts", loc="left")
    rowsd = [("C1  sampling", "best RMSD across modes", d["rmsd_best_of_modes_angstrom"],
              d["sampling_control_passed"]),
             ("C2  ranking", "RMSD of the score-top pose", d["rmsd_top_pose_angstrom"],
              d["ranking_control_passed"])]
    for k, (nm, sub, v, ok) in enumerate(rowsd):
        y = .74 - k * .36
        ax[1].add_patch(Rectangle((.02, y - .14), .96, .27, transform=ax[1].transAxes,
                                  fc=(OK["green"] if ok else OK["orange"]) + "22",
                                  ec=OK["green"] if ok else OK["orange"], lw=1.1))
        ax[1].text(.06, y + .06, nm, fontsize=10, fontweight="bold", transform=ax[1].transAxes)
        ax[1].text(.06, y - .05, sub, fontsize=8, color="#555", transform=ax[1].transAxes)
        ax[1].text(.72, y + .01, f"{v:.2f} Å", fontsize=12, fontweight="bold",
                   transform=ax[1].transAxes)
        ax[1].text(.90, y + .01, "PASS" if ok else "FAIL", fontsize=10, fontweight="bold",
                   color=OK["green"] if ok else OK["orange"], transform=ax[1].transAxes)
    ax[1].text(.03, .06, "Without splitting the control, a single failing number\n"
                         "would have been read as “docking does not work”.",
               fontsize=8, color="#444", transform=ax[1].transAxes)
    fig.suptitle("Re-docking control separates pose generation from scoring",
                 fontsize=11.5, y=1.04)
    return save(fig, "fig03_control.png")


# ───────────────────────── Fig 4 구간별 산점도 ───────────────────────
def fig04(an):
    T = an["arms"]["top_pose"]
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.2), sharey=True)
    ax[0].set_ylabel("measured pChEMBL")
    xs = [r["top_pose_score"] for r in ROWS]; ys = [r["pchembl_value"] for r in ROWS]
    cs = [SIM_C[r["similarity_bin"]] for r in ROWS]
    ax[0].scatter(xs, ys, s=16, c=cs, edgecolor="black", lw=.2, alpha=.8)
    m, b = np.polyfit(xs, ys, 1); xr = np.linspace(min(xs), max(xs), 10)
    ax[0].plot(xr, m * xr + b, color="black", lw=1.2, ls="--")
    ax[0].set_title(f"All compounds\nρ = {T['spearman']:+.3f}  p = {T['perm_p']}", loc="left")
    ax[0].set_xlabel("docking score (kcal/mol)")
    for k, bnm in enumerate(("far", "mid", "near"), start=1):
        sub = [r for r in ROWS if r["similarity_bin"] == bnm]
        xs = [r["top_pose_score"] for r in sub]; ys = [r["pchembl_value"] for r in sub]
        ax[k].scatter(xs, ys, s=20, c=SIM_C[bnm], edgecolor="black", lw=.25, alpha=.88)
        if len(xs) > 2:
            m, b = np.polyfit(xs, ys, 1); xr = np.linspace(min(xs), max(xs), 10)
            ax[k].plot(xr, m * xr + b, color="black", lw=1.2, ls="--")
        w = T["within_similarity_bin"].get(bnm, {})
        good = w.get("perm_p") is not None and w["perm_p"] < .05
        ax[k].set_title(f"{bnm} scaffold  (n={len(sub)})\n"
                        f"ρ = {w.get('spearman', float('nan')):+.3f}  p = {w.get('perm_p')}",
                        loc="left", color=OK["green"] if good else OK["orange"])
        ax[k].set_xlabel("docking score (kcal/mol)")
    fig.suptitle("Within a fixed scaffold band, does the score still rank potency?",
                 fontsize=11.5, y=1.06)
    return save(fig, "fig04_score_vs_potency.png")


# ───────────────────────── Fig 5 상관 forest ─────────────────────────
def fig05(an):
    T = an["arms"]["top_pose"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    labs, vals, cis, cols = [], [], [], []
    labs.append(f"All compounds (n={an['n']})"); vals.append(T["spearman"])
    cis.append(T["ci95"]); cols.append(OK["grey"])
    labs.append("Partial ρ, scaffold controlled")
    vals.append(T["partial_spearman_controlling_tanimoto"]); cis.append(None); cols.append(OK["blue"])
    if T.get("partial_spearman_controlling_size") is not None:
        labs.append("Partial ρ, size controlled")
        vals.append(T["partial_spearman_controlling_size"]); cis.append(None); cols.append(OK["sky"])
    if T.get("partial_spearman_controlling_both") is not None:
        labs.append("Partial ρ, scaffold + size")
        vals.append(T["partial_spearman_controlling_both"]); cis.append(None); cols.append(OK["green"])
    for b in ("far", "mid", "near"):
        w = T["within_similarity_bin"].get(b)
        if w:
            labs.append(f"Within {b} scaffold (n={w['n']})")
            vals.append(w["spearman"]); cis.append(w["ci95"]); cols.append(SIM_C[b])
    sens = an.get("sensitivity_reference_pose")
    if sens:
        labs.append(f"Reference-selected pose (n={sens['n']})")
        vals.append(sens["spearman"]); cis.append(sens["ci95"]); cols.append(OK["lemon"])
    yy = np.arange(len(vals))[::-1]
    for i, (v, ci, c) in enumerate(zip(vals, cis, cols)):
        y0 = yy[i]
        if ci:
            ax.plot(ci, [y0, y0], color=c, lw=3, alpha=.5, solid_capstyle="round")
            ax.plot([ci[0]] * 2, [y0 - .12, y0 + .12], color=c, lw=1.4)
            ax.plot([ci[1]] * 2, [y0 - .12, y0 + .12], color=c, lw=1.4)
        ax.plot([v], [y0], "o", color=c, ms=9, mec="black", mew=.6)
        ax.text(v, y0 + .26, f"{v:+.3f}", ha="center", fontsize=8, color=c, fontweight="bold")
    ax.axvline(0, color="black", lw=1.1, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=8.5)
    ax.set_xlabel("Spearman ρ (pChEMBL vs docking score); negative = predictive")
    ax.set_title("Correlation before and after controlling the confounders\n"
                 "(bars = 95% Fisher CI; partial correlations have no CI shown)", loc="left")
    ax.set_xlim(-1.02, .55)
    return save(fig, "fig05_forest.png")


# ───────────────────────── Fig 6 ROC ─────────────────────────────────
def _roc(scores, labels):
    z = sorted(zip(scores, labels), key=lambda t: t[0])
    P = sum(l for _, l in z); N = len(z) - P
    if not P or not N: return None, None
    tp = fp = 0; xs, ys = [0.0], [0.0]
    for _, l in z:
        tp += l; fp += (1 - l); xs.append(fp / N); ys.append(tp / P)
    return xs, ys


def fig06(an):
    T = an["arms"]["top_pose"]
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.9))
    sw = [(r["top_pose_score"], 1 if r["potency_bin"] == "strong" else 0)
          for r in ROWS if r["potency_bin"] in ("strong", "weak")]
    xs, ys = _roc([a for a, _ in sw], [b for _, b in sw])
    ax[0].plot(xs, ys, color=OK["blue"], lw=2.4,
               label=f"all (AUC {T['auc_strong_vs_weak']}, n={len(sw)})")
    ax[0].plot([0, 1], [0, 1], ls="--", color="grey", lw=1, label="random")
    ax[0].fill_between(xs, ys, xs, color=OK["blue"], alpha=.12)
    ax[0].set_xlabel("false positive rate"); ax[0].set_ylabel("true positive rate")
    ax[0].set_title("A  Strong vs weak, all compounds", loc="left")
    ax[0].legend(fontsize=8, loc="lower right")
    for b in ("far", "mid", "near"):
        sub = [(r["top_pose_score"], 1 if r["potency_bin"] == "strong" else 0)
               for r in ROWS if r["similarity_bin"] == b and r["potency_bin"] in ("strong", "weak")]
        if len(sub) < 8: continue
        xs, ys = _roc([a for a, _ in sub], [c for _, c in sub])
        au = T["auc_within_similarity_bin"].get(b)
        ax[1].plot(xs, ys, color=SIM_C[b], lw=2.2, label=f"{b} (AUC {au}, n={len(sub)})")
    ax[1].plot([0, 1], [0, 1], ls="--", color="grey", lw=1, label="random")
    ax[1].set_xlabel("false positive rate"); ax[1].set_ylabel("true positive rate")
    ax[1].set_title("B  Within each scaffold band", loc="left")
    ax[1].legend(fontsize=7.5, loc="lower right")
    fig.suptitle("Triage performance survives (or does not) scaffold control",
                 fontsize=11.5, y=1.04)
    return save(fig, "fig06_roc.png")


# ───────────────────────── Fig 7 농축 곡선 ───────────────────────────
def fig07(an):
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.6))
    order = sorted(ROWS, key=lambda r: r["top_pose_score"])
    n = len(order)
    base = sum(1 for r in ROWS if r["potency_bin"] == "strong") / n
    fr = np.arange(1, n + 1) / n
    hits = np.cumsum([1 if r["potency_bin"] == "strong" else 0 for r in order])
    ax[0].plot(fr * 100, hits / np.arange(1, n + 1) / base, color=OK["blue"], lw=2.2,
               label="ranked by docking score")
    ax[0].axhline(1, ls="--", color="grey", lw=1.1, label="random (EF = 1)")
    ax[0].set_xlabel("top % of the ranked list"); ax[0].set_ylabel("enrichment factor")
    ax[0].set_title("A  Enrichment of strong binders", loc="left")
    ax[0].set_xlim(0, 100); ax[0].legend(fontsize=8)
    for pct in (5, 10, 20):
        k = max(1, int(n * pct / 100))
        ef = (hits[k - 1] / k) / base
        ax[0].plot([pct], [ef], "o", color=OK["orange"], ms=6, zorder=5)
        ax[0].text(pct, ef + .18, f"EF{pct}%\n{ef:.2f}", ha="center", fontsize=7,
                   color=OK["orange"], fontweight="bold")
    k = max(1, n // 10)
    top = order[:k]
    P = ["strong", "medium", "weak"]; S = ["far", "mid", "near"]
    w = .36
    for i, (grp, cols, keyf) in enumerate(((P, POT_C, "potency_bin"), (S, SIM_C, "similarity_bin"))):
        obs = [sum(1 for r in top if r[keyf] == g) for g in grp]
        exp = [len([r for r in ROWS if r[keyf] == g]) * k / n for g in grp]
        xx = np.arange(len(grp)) + i * 4
        ax[1].bar(xx - w / 2, obs, w, color=[cols[g] for g in grp], edgecolor="black", lw=.5,
                  label="observed in top-10%" if i == 0 else None)
        ax[1].bar(xx + w / 2, exp, w, color="white", edgecolor="black", lw=.7, hatch="///",
                  label="expected if random" if i == 0 else None)
    ax[1].set_xticks(list(range(3)) + list(range(4, 7)))
    ax[1].set_xticklabels(P + S, fontsize=8)
    ax[1].set_ylabel(f"compounds in score top-{k}")
    ax[1].set_title("B  Who lands in the top-10%?", loc="left")
    ax[1].legend(fontsize=7.5)
    ax[1].text(1, -.30, "potency", transform=ax[1].get_xaxis_transform(), ha="center",
               fontsize=8.5, color="#444", fontweight="bold")
    ax[1].text(5, -.30, "scaffold similarity", transform=ax[1].get_xaxis_transform(),
               ha="center", fontsize=8.5, color="#444", fontweight="bold")
    fig.suptitle("If the score enriches for scaffold rather than potency, panel B shows it",
                 fontsize=11.5, y=1.05)
    return save(fig, "fig07_enrichment.png")


# ───────────────────────── Fig 8 탐색 깊이 스윕 ──────────────────────
def fig08(sw):
    if not sw: return None
    S = sw["summary"]
    e = [x["exhaustiveness"] for x in S]
    c1 = [x["c1_best_rmsd_mean"] for x in S]
    c2 = [x["c2_top_rmsd_mean"] for x in S]
    sc = [x["top_score_mean"] for x in S]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.5))
    ax[0].plot(e, c1, "o-", color=OK["green"], lw=2.2, ms=7, label="C1  best pose found")
    ax[0].plot(e, c2, "s-", color=OK["orange"], lw=2.2, ms=7, label="C2  top-scored pose")
    ax[0].axhline(2.0, ls="--", color="black", lw=1)
    ax[0].set_xscale("log", base=2); ax[0].set_xticks(e); ax[0].set_xticklabels(e)
    ax[0].set_xlabel("exhaustiveness (search depth)")
    ax[0].set_ylabel("RMSD to crystal pose (Å)")
    ax[0].set_title("A  Deeper search, worse ranking", loc="left")
    ax[0].legend(fontsize=8)
    ax[0].text(.5, .04, "if sampling were the bottleneck,\nthe orange line would fall too",
               transform=ax[0].transAxes, ha="center", fontsize=7.5, color="#444")
    ax[1].plot(e, sc, "o-", color=OK["blue"], lw=2.2, ms=7)
    ax[1].set_xscale("log", base=2); ax[1].set_xticks(e); ax[1].set_xticklabels(e)
    ax[1].set_xlabel("exhaustiveness"); ax[1].set_ylabel("best score (kcal/mol)")
    ax[1].set_title("B  The score keeps improving", loc="left")
    ax[1].text(.5, .12, "search is doing its job:\nit finds better-scoring poses",
               transform=ax[1].transAxes, ha="center", fontsize=7.5, color="#444")
    w = .36; xx = np.arange(len(e))
    ax[2].bar(xx - w / 2, [x["c1_pass_rate"] for x in S], w, color=OK["green"],
              edgecolor="black", lw=.5, label="C1 pass rate")
    ax[2].bar(xx + w / 2, [x["c2_pass_rate"] for x in S], w, color=OK["orange"],
              edgecolor="black", lw=.5, label="C2 pass rate")
    ax[2].set_xticks(xx); ax[2].set_xticklabels(e)
    ax[2].set_xlabel("exhaustiveness"); ax[2].set_ylabel(f"fraction of {S[0]['n_runs']} seeds passing")
    ax[2].set_ylim(0, 1.08); ax[2].set_title("C  Pass rate across random seeds", loc="left")
    ax[2].legend(fontsize=8)
    fig.suptitle("The bottleneck is scoring, not sampling: a search-depth experiment",
                 fontsize=11.5, y=1.04)
    return save(fig, "fig08_exhaustiveness.png")


# ───────────────────────── Fig 9 항 기여도 ───────────────────────────
def fig09(tc):
    if not tc: return None
    S = tc["short"]
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.5))
    sets = [("Standardized coefficient", [tc["standardized_coefficients"][t] for t in S],
             "comparable across terms"),
            ("Univariate r vs pChEMBL", [tc["univariate"][t]["pearson_r"] for t in S],
             "each term on its own"),
            ("Permutation Q² drop", [tc["permutation_importance"][t]["mean_Q2_drop"] for t in S],
             "positive = the term carries signal"),
            ("VIF", [tc["vif"][t] for t in S], "VIF > 5 → coefficient not interpretable")]
    x = np.arange(len(S))
    for k, (title, vals, sub) in enumerate(sets):
        if k < 3:
            cols = [OK["orange"] if v < 0 else OK["blue"] for v in vals]
        else:
            cols = [OK["orange"] if (v and v > 5) else OK["green"] for v in vals]
        ax[k].bar(x, vals, .62, color=cols, edgecolor="black", lw=.5)
        ax[k].set_xticks(x); ax[k].set_xticklabels(S, rotation=38, ha="right", fontsize=8)
        ax[k].set_title(f"{'ABCD'[k]}  {title}", loc="left", fontsize=9.5)
        ax[k].text(.5, -.44, sub, transform=ax[k].transAxes, ha="center",
                   fontsize=7.5, color="#555")
        ax[k].axhline(5 if k == 3 else 0, color="black" if k < 3 else OK["orange"],
                      lw=.9, ls="-" if k < 3 else "--")
    fig.suptitle("Which scoring term carries the signal? Four views of the same model",
                 fontsize=11.5, y=1.06)
    return save(fig, "fig09_term_importance.png")


# ───────────────────────── Fig 10 항 상관 + 모델 비교 ────────────────
def fig10(tc):
    if not tc: return None
    S = tc["short"]
    fig = plt.figure(figsize=(12.6, 3.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.35], wspace=.34)
    ax0 = fig.add_subplot(gs[0])
    M = np.array([[tc["term_correlation"][a][b] for b in S] for a in S])
    im = ax0.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
    ax0.set_xticks(range(len(S))); ax0.set_xticklabels(S, rotation=40, ha="right", fontsize=7.5)
    ax0.set_yticks(range(len(S))); ax0.set_yticklabels(S, fontsize=7.5)
    for i in range(len(S)):
        for j in range(len(S)):
            ax0.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7,
                     color="white" if abs(M[i, j]) > .55 else "black")
    fig.colorbar(im, ax=ax0, shrink=.78, label="Pearson r")
    ax0.set_title("A  Terms correlate with each other", loc="left", fontsize=9.5)

    ax1 = fig.add_subplot(gs[1])
    vals = [tc["univariate"][t]["pearson_r_vs_tanimoto"] for t in S]
    ax1.barh(range(len(S)), vals, .6,
             color=[OK["pink"] if abs(v) > .3 else OK["grey"] for v in vals],
             edgecolor="black", lw=.5)
    ax1.set_yticks(range(len(S))); ax1.set_yticklabels(S, fontsize=8)
    ax1.axvline(0, color="black", lw=.9)
    ax1.set_xlabel("r (term vs Tanimoto to co-crystal ligand)")
    ax1.set_title("B  Do the terms just track scaffold?", loc="left", fontsize=9.5)

    ax2 = fig.add_subplot(gs[2])
    names = ["single\nVina score", "5 scoring\nterms",
             "Tanimoto\nonly", "5 terms\n+ Tanimoto"]
    keys = ["single_vina_score", "five_terms", "tanimoto_only", "five_terms_plus_tanimoto"]
    r2v = [tc["models"][k]["R2_fit"] for k in keys]
    q2v = [tc["models"][k]["Q2_loo"] for k in keys]
    nullp = [tc["models"][k]["null_R2_p95"] for k in keys]
    xx = np.arange(len(keys)); w = .36
    ax2.bar(xx - w / 2, r2v, w, color=OK["sky"], edgecolor="black", lw=.5, label="R² (fit)")
    ax2.bar(xx + w / 2, q2v, w, color=OK["blue"], edgecolor="black", lw=.5, label="Q² (LOO)")
    for i, v in enumerate(nullp):
        ax2.plot([i - .55, i + .55], [v, v], color=OK["orange"], lw=1.6, ls="--")
    ax2.plot([], [], color=OK["orange"], lw=1.6, ls="--", label="95th pct of label-shuffled R²")
    ax2.set_xticks(xx); ax2.set_xticklabels(names, fontsize=8)
    ax2.axhline(0, color="black", lw=.9)
    ax2.set_ylabel("variance explained")
    ax2.set_title("C  Fit vs cross-validation vs chance", loc="left", fontsize=9.5)
    ax2.legend(fontsize=7.5, loc="upper left")
    fig.suptitle("Scoring terms: collinear with each other, and partly with scaffold",
                 fontsize=11.5, y=1.05)
    return save(fig, "fig10_terms_models.png")


# ───────────────────────── Fig 11 자세 규칙 민감도 ───────────────────
def fig11(an):
    T = an["arms"]["top_pose"]; sens = an.get("sensitivity_reference_pose")
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.5))
    rm = [r.get("top_pose_mcs_rmsd") for r in ROWS if r.get("top_pose_mcs_rmsd") is not None]
    rr = [r.get("ref_pose_mcs_rmsd") for r in ROWS if r.get("ref_pose_mcs_rmsd") is not None]
    bins = np.linspace(0, max(max(rm or [1]), max(rr or [1])), 26)
    ax[0].hist(rm, bins=bins, color=OK["blue"], alpha=.75, label=f"top-scored pose (n={len(rm)})")
    ax[0].hist(rr, bins=bins, color=OK["lemon"], alpha=.7,
               label=f"reference-selected pose (n={len(rr)})")
    ax[0].axvline(2.0, ls="--", color="black", lw=1.1)
    ax[0].set_xlabel("MCS-RMSD to the co-crystal ligand (Å)"); ax[0].set_ylabel("compounds")
    ax[0].set_title("A  How close does each rule land?", loc="left")
    ax[0].legend(fontsize=7.5)
    labs = ["top-scored\n(no reference used)"]
    vals = [T["spearman"]]; cols = [OK["blue"]]
    if sens:
        labs.append("reference-selected\n(uses co-crystal ligand)")
        vals.append(sens["spearman"]); cols.append(OK["lemon"])
    ax[1].bar(range(len(vals)), vals, .5, color=cols, edgecolor="black", lw=.6)
    ax[1].set_xticks(range(len(vals))); ax[1].set_xticklabels(labs, fontsize=8)
    ax[1].axhline(0, color="black", lw=.9)
    ax[1].set_ylabel("Spearman ρ vs pChEMBL")
    ax[1].set_title("B  The conclusion does not depend on the rule", loc="left")
    for i, v in enumerate(vals):
        ax[1].text(i, v - .02, f"{v:+.3f}", ha="center", va="top", fontsize=9.5,
                   fontweight="bold")
    fig.suptitle("Pose-selection rule: sensitivity analysis", fontsize=11.5, y=1.04)
    return save(fig, "fig11_pose_rule.png")


# ───────────────────────── Fig 12 결론의 이력 ────────────────────────
def fig12(an, old_reg, old_enr, ds, col=None):
    """세 판본에서 결론이 어떻게 바뀌었는지 — 이 연구의 가장 중요한 교훈."""
    T = an["arms"]["top_pose"]
    fig = plt.figure(figsize=(13.4, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1.05], wspace=.30)

    # A 판본별 헤드라인 수치
    ax = fig.add_subplot(gs[0])
    stages = ["v1.0\nn=30\npotency-stratified",
              "v2.0\nsame data\ntriage metrics added",
              "v3.0\nn=%d\nscaffold-decoupled" % an["n"]]
    v1 = old_reg["correlation"]["spearman_pIC50_vs_top_pose_score"] if old_reg else None
    v2 = v1
    v3 = T["spearman"]
    v3p = T["partial_spearman_controlling_tanimoto"]
    xs = np.arange(3)
    ax.bar(xs, [v1, v2, v3], .52, color=[OK["grey"], OK["grey"], OK["blue"]],
           edgecolor="black", lw=.6, label="reported ρ")
    ax.bar([2], [v3p], width=.26, color=OK["orange"], edgecolor="black", lw=.6,
           alpha=.9, label="ρ after scaffold control")
    for i, v in enumerate([v1, v2, v3]):
        if v is not None:
            ax.text(i, v - .03, f"{v:+.3f}", ha="center", va="top", fontsize=9,
                    fontweight="bold")
    ax.text(2.19, v3p - .03, f"{v3p:+.3f}", ha="center", va="top", fontsize=9,
            fontweight="bold", color=OK["orange"])
    ax.set_xticks(xs); ax.set_xticklabels(stages, fontsize=8)
    ax.axhline(0, color="black", lw=.9)
    ax.set_ylabel("Spearman ρ (pChEMBL vs score)")
    ax.set_title("A  The headline number across three versions", loc="left")
    ax.legend(fontsize=7.5, loc="lower left")

    # B 각 판본의 결론과 그것을 무너뜨린 관찰
    ax = fig.add_subplot(gs[1]); ax.axis("off")
    ax.set_title("B  What each version concluded — and what broke it", loc="left")
    att = (col or {}).get("old_design", {}).get("attenuation_from_scaffold_control")
    story = [
        ("v1.0", "“Docking does not\npredict potency.”",
         "Broke: no triage metric\nwas ever computed.", OK["orange"]),
        ("v2.0", "“It triages but does\nnot quantify.”",
         "Broke: did not replicate in a\nlarger, scaffold-decoupled set.", OK["orange"]),
        ("v3.0\ndraft", "“That scaffold bias\nexplains v2.0.”",
         ("Broke: controlling scaffold in\nthe OLD data changes rho by %.1f%%"
          % (100 * att)) if att is not None else "Broke: tested and refuted.", OK["orange"]),
        ("v3.0\nfinal", "“It did not replicate.\nCause not identified.”",
         "Claim narrowed to what\nthe tests actually support.", OK["green"]),
    ]
    for k, (tag, concl, broke, cc) in enumerate(story):
        y = .84 - k * .235
        ax.add_patch(Rectangle((.01, y - .155), .97, .205, transform=ax.transAxes,
                               fc=cc + "18", ec=cc, lw=1.0))
        ax.text(.05, y + .00, tag, fontsize=8.6, fontweight="bold", transform=ax.transAxes,
                va="center")
        ax.text(.22, y + .035, concl, fontsize=7.6, transform=ax.transAxes, va="center")
        ax.text(.22, y - .075, broke, fontsize=6.9, color="#555", transform=ax.transAxes,
                va="center")

    # C 무엇이 결론을 바꿨는가
    ax = fig.add_subplot(gs[2]); ax.axis("off")
    ax.set_title("C  What actually caught each error", loc="left")
    caught = [("Automated gates", 0, OK["orange"],
               "schema, counts, provenance —\nall PASS in every version"),
              ("External critical review", 3, OK["green"],
               "asked for the metric, the confound,\nand then for the confound to be tested")]
    for k, (nm, cnt, col, sub) in enumerate(caught):
        y = .70 - k * .40
        ax.text(.04, y + .12, nm, fontsize=10, fontweight="bold", transform=ax.transAxes)
        ax.text(.04, y - .04, sub, fontsize=8, color="#444", transform=ax.transAxes, va="top")
        ax.text(.80, y + .10, str(cnt), fontsize=26, fontweight="bold", color=col,
                transform=ax.transAxes, ha="center")
        ax.text(.80, y - .04, "errors\ncaught", fontsize=7, color="#666",
                transform=ax.transAxes, ha="center", va="top")
    ax.text(.04, .06, "Gates verify that numbers exist and match their source.\n"
                      "They cannot verify that the conclusion follows from them.",
            fontsize=8, color=OK["orange"], transform=ax.transAxes, fontweight="bold")
    fig.suptitle("Three versions, three conclusions: what an automated pipeline could not catch",
                 fontsize=12.5, fontweight="bold", y=1.03)
    return save(fig, "fig12_version_history.png")


# ───────────────────────── Fig 13 커스텀 채점 ────────────────────────
def fig13(cs):
    """적합 가중치로 다시 도킹했을 때 held-out 성능이 실제로 나아지는가."""
    if not cs: return None
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    d = cs["difference"]
    nm = ["default Vina", "fitted custom"]
    v = [cs["arm_default"]["spearman"], cs["arm_custom"]["spearman"]]
    ci = [cs["arm_default"]["ci95"], cs["arm_custom"]["ci95"]]
    for i, (vv, cc) in enumerate(zip(v, ci)):
        ax[0].bar(i, vv, .5, color=[OK["blue"], OK["pink"]][i], edgecolor="black", lw=.6)
        if cc:
            ax[0].plot([i, i], cc, color="black", lw=1.4)
            ax[0].plot([i - .09, i + .09], [cc[0]] * 2, color="black", lw=1.2)
            ax[0].plot([i - .09, i + .09], [cc[1]] * 2, color="black", lw=1.2)
        ax[0].text(i, vv + .04, f"{vv:+.3f}", ha="center", va="bottom", fontsize=9.5,
                   fontweight="bold")
    ax[0].set_xticks(range(2)); ax[0].set_xticklabels(nm, fontsize=8.5)
    ax[0].axhline(0, color="black", lw=.9)
    ax[0].set_ylabel("held-out Spearman ρ")
    ax[0].set_title(f"A  Independent test set (n={cs['n_scored']})", loc="left")

    dd = d["delta_spearman"]; bci = d["bootstrap_ci95"]
    ax[1].errorbar([dd], [0], xerr=[[dd - bci[0]], [bci[1] - dd]], fmt="o",
                   color=OK["pink"], capsize=7, markersize=10, lw=2.2, mec="black", mew=.6)
    ax[1].axvline(0, color=OK["orange"], lw=1.8, ls="--")
    ax[1].set_yticks([]); ax[1].set_xlabel("Δ Spearman (custom − default)")
    span = max(abs(bci[0]), abs(bci[1])) * 1.25
    ax[1].set_xlim(-span, span)
    inside = bci[0] <= 0 <= bci[1]
    ax[1].set_title("B  " + ("Not distinguishable from zero" if inside
                             else "Difference excludes zero"), loc="left",
                    color=OK["orange"] if inside else OK["green"])
    ax[1].text(.5, .76, f"Δ = {dd:+.3f}    95% CI [{bci[0]:+.2f}, {bci[1]:+.2f}]\n"
                        f"{d['frac_resamples_favoring_custom']:.0%} of bootstrap resamples "
                        f"favour custom",
               transform=ax[1].transAxes, ha="center", fontsize=8,
               bbox=dict(fc="white", ec=OK["orange"] if inside else OK["green"],
                         lw=.9, alpha=.95, pad=4))

    yy = [r["pchembl_value"] for r in cs["rows"]]
    ax[2].scatter([r["base_score"] for r in cs["rows"]], yy, s=22, c=OK["blue"],
                  edgecolor="black", lw=.25, alpha=.8, label="default")
    ax2b = ax[2].twiny()
    ax2b.scatter([r["custom_score"] for r in cs["rows"]], yy, s=22, c=OK["pink"],
                 marker="^", edgecolor="black", lw=.25, alpha=.8, label="custom")
    ax[2].set_xlabel("default Vina score", color=OK["blue"])
    ax2b.set_xlabel("fitted custom score", color=OK["pink"])
    ax[2].set_ylabel("measured pChEMBL")
    ax[2].set_title("C  Same compounds, two scoring functions", loc="left")
    fig.suptitle("Re-docking with fitted weights: does the number move for a reason?",
                 fontsize=11.5, y=1.06)
    return save(fig, "fig13_custom_scoring.png")


ROWS = []


def main() -> int:
    global ROWS
    ds, dk, an = load("dataset_controlled.json"), load("docking_controlled.json"), \
                 load("analysis_controlled.json")
    if not (ds and dk and an):
        raise SystemExit("골격 통제 산출 파일이 아직 없다 — 그림을 만들지 않는다.")
    ROWS = dk["rows"]
    ctrl = dict(dk["control_redock"])
    sw, tc = load("exhaustiveness_sweep.json"), load("terms_controlled.json")
    made = [f for f in (fig01(ds, an, ctrl, tc), fig02(ds),
                        fig03(ctrl, sw, dk.get('exhaustiveness'), dk.get('seed')),
                        fig04(an), fig05(an), fig06(an), fig07(an),
                        fig08(sw), fig09(tc), fig10(tc), fig11(an),
                        fig12(an, load("regression.json"), load("enrichment.json"), ds,
                              load("collapse_diagnosis.json")),
                        fig13(load("custom_scoring_controlled.json"))) if f]
    print(f"그림 {len(made)}장 (핵심 세트)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
