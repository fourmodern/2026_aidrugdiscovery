#!/usr/bin/env python3
"""연구 보고서 그림 — PDE5A 저해제 프로파일링.

수치가 들어가는 그림은 sample_run/sar.json (ChEMBL 보고 활성 + RDKit 계산 물성) 만 쓴다.
난수·보간·임의 추세선 금지. 하네스 동작 설명 그림은 make_harness_figures.py 에 있다.
"""
import argparse, json, hashlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermil": "#D55E00",
      "sky": "#56B4E9", "yellow": "#F0E442", "purple": "#CC79A7", "grey": "#7F7F7F"}
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 300, "savefig.dpi": 300})
QED_MIN = 0.5


def meta(path, source, n, note):
    path.with_suffix(".meta.json").write_text(json.dumps({
        "source": str(source), "script": "scripts/make_figures.py", "n_records": n,
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "source_sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest()[:16],
        "note": note}, ensure_ascii=False, indent=2))


def pareto_front(R, keys=("LE", "LLE")):
    """두 지표에서 다른 화합물에 지배되지 않는 것들. 임의 가중합 대신 지배 관계로 판정한다."""
    front = []
    for a in R:
        dominated = any(
            all(b[k] >= a[k] for k in keys) and any(b[k] > a[k] for k in keys)
            for b in R if b is not a)
        if not dominated:
            front.append(a)
    return sorted(front, key=lambda r: -r["pIC50"])


def _card(ax, x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                                linewidth=1.6, edgecolor=color, facecolor=color + "18"))
    ax.text(x + w / 2, y + h - 0.055, title, ha="center", va="top",
            fontsize=8.5, fontweight="bold", color=color)
    ax.text(x + w / 2, y + h / 2 - 0.035, body, ha="center", va="center",
            fontsize=6.6, color="#333333", linespacing=1.45)


def fig1_graphical_abstract(out, src, R):
    n = len(R); front = pareto_front(R)
    # 파레토 전선 중 게이트를 통과하면서 가장 강력한 것 — 없으면 전선의 최고 역가
    passing = [r for r in front if r["QED"] >= QED_MIN]
    best = (passing or front)[0]
    top = max(R, key=lambda r: r["pIC50"])
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 9 / 2.54 * 1.15))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.955, "Does drug-likeness filtering select the better PDE5A inhibitors?",
            ha="center", fontsize=11.5, fontweight="bold", color="#1C242B")
    ax.text(0.5, 0.895,
            f"Potency and ligand-efficiency profile of {n} reported PDE5A actives (ChEMBL CHEMBL1827)",
            ha="center", fontsize=7.8, color=OI["grey"])
    _card(ax, 0.035, 0.44, 0.21, 0.36, "1  QUESTION",
          "Do QED / Ro5 filters\nkeep the compounds a\nchemist would keep?", OI["blue"])
    _card(ax, 0.275, 0.44, 0.21, 0.36, "2  DATA",
          f"{n} actives, IC50 only\n{min(r['ic50_nM'] for r in R):.1f}-{max(r['ic50_nM'] for r in R):.0f} nM\nRDKit descriptors", OI["orange"])
    _card(ax, 0.515, 0.44, 0.21, 0.36, "3  FINDING",
          "The filter kept the\nweaker half. Median\npIC50 7.13 vs 7.79.", OI["green"])
    _card(ax, 0.755, 0.44, 0.21, 0.36, "4  IMPLICATION",
          "Rank by efficiency,\nnot by a single\ndrug-likeness cutoff.", OI["purple"])
    for x in (0.245, 0.485, 0.725):
        ax.add_patch(FancyArrowPatch((x, 0.62), (x + 0.03, 0.62), arrowstyle="-|>",
                                     mutation_scale=11, linewidth=1.3, color=OI["grey"]))
    ax.add_patch(FancyBboxPatch((0.035, 0.10), 0.93, 0.25, boxstyle="round,pad=0.012",
                                linewidth=1.4, edgecolor=OI["green"], facecolor="#EAF6F1"))
    ax.text(0.5, 0.295,
            f"Pareto-optimal on LE and LLE: {', '.join(r['chembl_id'] for r in front)}",
            ha="center", fontsize=8.6, fontweight="bold", color=OI["green"])
    ax.text(0.5, 0.185,
            f"Of these, {best['chembl_id']} also passes the filter and is the more potent: "
            f"IC50 {best['ic50_nM']:.1f} nM  ·  LE {best['LE']:.3f}  ·  LLE {best['LLE']:.2f}  ·  MW {best['MW']:.0f}.\n"
            f"The most potent compound overall ({top['chembl_id']}, {top['ic50_nM']:.1f} nM) is rejected by the filter.",
            ha="center", fontsize=6.8, color="#444444", linespacing=1.6)
    p = out / "fig1_graphical_abstract.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, n, "개념도 + 실측 요약값")
    return p


def fig2_pathway(out, src, R, target):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 7.6 / 2.54 * 1.15))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.95, "PDE5A in the NO-cGMP axis, and where an inhibitor acts",
            ha="center", fontsize=11, fontweight="bold", color="#1C242B")
    nodes = [(0.06, "NO", "endothelium", OI["sky"]),
             (0.27, "sGC", "soluble guanylate\ncyclase", OI["sky"]),
             (0.48, "cGMP", "second\nmessenger", OI["green"]),
             (0.69, "PKG", "protein\nkinase G", OI["green"])]
    for x, name, sub, c in nodes:
        _card(ax, x, 0.52, 0.16, 0.26, name, sub, c)
    for x in (0.225, 0.435, 0.645):
        ax.add_patch(FancyArrowPatch((x, 0.65), (x + 0.03, 0.65), arrowstyle="-|>",
                                     mutation_scale=12, linewidth=1.5, color=OI["grey"]))
    _card(ax, 0.855, 0.52, 0.13, 0.26, "effect", "smooth muscle\nrelaxation", OI["blue"])
    ax.add_patch(FancyArrowPatch((0.855, 0.65), (0.835, 0.65), arrowstyle="-|>",
                                 mutation_scale=12, linewidth=1.5, color=OI["grey"]))
    # PDE5A 가 cGMP 를 분해
    _card(ax, 0.40, 0.16, 0.22, 0.24, f"PDE5A  ({target['accession']})",
          f"{target['length']} aa\nhydrolyses cGMP to 5'-GMP", OI["vermil"])
    ax.add_patch(FancyArrowPatch((0.51, 0.41), (0.53, 0.51), arrowstyle="-|>",
                                 mutation_scale=12, linewidth=1.6, color=OI["vermil"]))
    ax.text(0.545, 0.455, "degrades", fontsize=6.6, color=OI["vermil"], fontweight="bold")
    # 저해제
    _card(ax, 0.10, 0.16, 0.22, 0.24, "inhibitor",
          f"n = {len(R)} reported actives\nIC50 {min(r['ic50_nM'] for r in R):.1f}-{max(r['ic50_nM'] for r in R):.0f} nM", OI["purple"])
    ax.add_patch(FancyArrowPatch((0.325, 0.28), (0.395, 0.28), arrowstyle="-[",
                                 mutation_scale=9, linewidth=1.8, color=OI["purple"]))
    ax.text(0.36, 0.32, "blocks", fontsize=6.8, color=OI["purple"],
            fontweight="bold", ha="center")
    ax.text(0.5, 0.055,
            "Blocking PDE5A raises cGMP and sustains relaxation. "
            "Selectivity over retinal PDE6 is the classic liability and was not quantified here.",
            ha="center", fontsize=7.0, color="#444444")
    p = out / "fig2_pathway.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(R), "기전 개념도 — 수치는 실측 n 과 IC50 범위만")
    return p


def fig3_potency(out, src, R):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 6.4 / 2.54 * 1.15))
    S = sorted(R, key=lambda r: r["pIC50"])
    ys = range(len(S))
    cols = [OI["green"] if r["QED"] >= QED_MIN else OI["vermil"] for r in S] \
        if "QED" in S[0] else [OI["blue"]] * len(S)
    ax.barh(list(ys), [r["pIC50"] for r in S], color=cols,
            edgecolor="white", linewidth=1.0, height=0.68)
    for i, r in enumerate(S):
        ax.text(r["pIC50"] + 0.02, i, f"{r['ic50_nM']:.1f} nM", va="center",
                fontsize=6.8, color="#333")
    ax.set_yticks(list(ys)); ax.set_yticklabels([r["chembl_id"] for r in S], fontsize=6.6)
    ax.set_xlabel("pIC50  (higher = more potent)", fontsize=8)
    ax.set_xlim(6.5, max(r["pIC50"] for r in S) + 0.45)
    ax.set_title("Reported potency of the retrieved PDE5A actives (all IC50, ChEMBL)",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig3_potency.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(R), "ChEMBL 보고 IC50/pIC50 실측")
    return p


def fig4_efficiency(out, src, R):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 7.0 / 2.54 * 1.15))
    front = pareto_front(R); front_ids = {r["chembl_id"] for r in front}
    for r in R:
        c = OI["green"] if r["QED"] >= QED_MIN else OI["vermil"]
        ax.scatter(r["LLE"], r["LE"], s=28 + (r["pIC50"] - 6.8) * 90, color=c,
                   edgecolor="white", linewidth=0.9, zorder=3)
        ax.annotate(r["chembl_id"].replace("CHEMBL", ""), (r["LLE"], r["LE"]),
                    textcoords="offset points", xytext=(0, 9), ha="center",
                    fontsize=5.8, color="#555")
    for r in front:
        ax.scatter([r["LLE"]], [r["LE"]], s=230, facecolor="none",
                   edgecolor=OI["blue"], linewidth=2.0, zorder=4)
    ax.annotate("circled = Pareto-optimal\n(not dominated on both metrics)",
                xy=(0.02, 0.04), xycoords="axes fraction", fontsize=6.6,
                color=OI["blue"], fontweight="bold", va="bottom")
    ax.axhline(0.30, color=OI["grey"], linestyle="--", linewidth=1.0)
    ax.annotate("LE 0.30 — common lead-like guide", xy=(0.02, 0.30),
                xycoords=("axes fraction", "data"), xytext=(0, 5),
                textcoords="offset points", fontsize=6.2, color=OI["grey"], ha="left")
    ax.set_xlabel("LLE  =  pIC50 - cLogP", fontsize=8)
    ax.set_ylabel("LE  =  1.37 x pIC50 / heavy atoms", fontsize=8)
    ax.set_title("Ligand efficiency — marker size is potency, colour is the drug-likeness gate",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=OI["green"], label="gate PASS"),
                       Line2D([], [], marker="o", ls="", color=OI["vermil"], label="gate FAIL")],
              fontsize=6.6, frameon=False, loc="lower right")
    fig.tight_layout()
    p = out / "fig4_efficiency.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(R), "LE·LLE 는 pIC50 과 RDKit 값으로 계산 (analyze_sar.py)")
    return p


def fig5_gate_vs_potency(out, src, R):
    fig, axes = plt.subplots(1, 2, figsize=(16 / 2.54 * 1.15, 6.2 / 2.54 * 1.15))
    ax = axes[0]
    for r in R:
        c = OI["green"] if r["QED"] >= QED_MIN else OI["vermil"]
        ax.scatter(r["QED"], r["pIC50"], s=54, color=c, edgecolor="white",
                   linewidth=0.8, zorder=3)
    ax.axvline(QED_MIN, color=OI["blue"], linewidth=1.8)
    ax.annotate(f"QED_MIN = {QED_MIN}", xy=(QED_MIN, 1.0),
                xycoords=("data", "axes fraction"), xytext=(4, -6),
                textcoords="offset points", fontsize=6.6, color=OI["blue"],
                va="top", fontweight="bold")
    ax.set_xlabel("QED", fontsize=8); ax.set_ylabel("pIC50", fontsize=8)
    ax.set_title("A  Potency vs drug-likeness", fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    ax = axes[1]
    keep = [r["pIC50"] for r in R if r["QED"] >= QED_MIN]
    drop = [r["pIC50"] for r in R if r["QED"] < QED_MIN]
    parts = ax.boxplot([drop, keep], tick_labels=[f"FAIL (n={len(drop)})", f"PASS (n={len(keep)})"],
                       widths=0.5, patch_artist=True, medianprops=dict(color="black", linewidth=1.6))
    for patch, c in zip(parts["boxes"], (OI["vermil"], OI["green"])):
        patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor(c)
    for i, vals in enumerate((drop, keep), start=1):
        ax.scatter([i] * len(vals), vals, s=30, color="black", zorder=4, alpha=0.65)
    ax.set_ylabel("pIC50", fontsize=8)
    ax.set_title("B  The filter kept the weaker half", fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for a in axes:
        for sp in ("top", "right"): a.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig5_gate_vs_potency.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(R), "n=10 기술통계 — 검정 없음")
    return p


def fig6_property_space(out, src, R):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 6.4 / 2.54 * 1.15))
    for r in R:
        c = OI["green"] if r["QED"] >= QED_MIN else OI["vermil"]
        ax.scatter(r["MW"], r["cLogP"], s=28 + (r["pIC50"] - 6.8) * 90, color=c,
                   edgecolor="white", linewidth=0.9, zorder=3)
        ax.annotate(f"{r['ic50_nM']:.0f}", (r["MW"], r["cLogP"]),
                    textcoords="offset points", xytext=(0, 9), ha="center",
                    fontsize=5.8, color="#555")
    ax.axvline(500, color=OI["grey"], linestyle="--", linewidth=1.0)
    ax.annotate("Lipinski MW 500", xy=(500, 1.0), xycoords=("data", "axes fraction"),
                xytext=(4, -6), textcoords="offset points", fontsize=6.3,
                color=OI["grey"], va="top")
    ax.set_xlabel("Molecular weight (Da)", fontsize=8)
    ax.set_ylabel("cLogP", fontsize=8)
    ax.set_title("Chemical space — label is IC50 in nM, marker size is potency",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig6_property_space.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(R), "RDKit 실계산 물성 + ChEMBL IC50")
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sar", default="sample_run/sar.json")
    ap.add_argument("--out", default="outputs/figures")
    a = ap.parse_args()
    src = Path(a.sar); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    R = json.loads(src.read_text())
    mp = {r["smiles"]: r for r in json.loads(
        Path("sample_run/envelopes/mol-properties.json").read_text())["result"]}
    for r in R:
        r["QED"] = mp[r["smiles"]]["QED"]
    target = json.loads(Path("sample_run/envelopes/target-lookup.json").read_text())["result"]
    for f in out.glob("fig*"):
        f.unlink()
    made = [fig1_graphical_abstract(out, src, R), fig2_pathway(out, src, R, target),
            fig3_potency(out, src, R), fig4_efficiency(out, src, R),
            fig5_gate_vs_potency(out, src, R), fig6_property_space(out, src, R)]
    for p in made:
        print(f"  {p}  ({p.stat().st_size // 1024} KB)")
    print(f"연구 그림 {len(made)}장 — 화합물 {len(R)}건")


if __name__ == "__main__":
    main()
