#!/usr/bin/env python3
"""도킹 연구의 핵심 결과 그림 — 점수 대 실측 역가, 회귀 비교, 재도킹 대조."""
import json, hashlib
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
SR = ROOT / "sample_run"
OUT = ROOT / "outputs" / "figures"
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermil": "#D55E00",
      "sky": "#56B4E9", "purple": "#CC79A7", "grey": "#7F7F7F"}
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 300, "savefig.dpi": 300})
COL = {"strong": OI["green"], "medium": OI["orange"], "weak": OI["vermil"]}


def meta(p, src, n, note):
    p.with_suffix(".meta.json").write_text(json.dumps({
        "source": str(src), "script": "scripts/make_docking_figures.py", "n_records": n,
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "source_sha256": hashlib.sha256(Path(src).read_bytes()).hexdigest()[:16],
        "note": note}, ensure_ascii=False, indent=2))


def fig10_score_vs_potency(reg, src):
    R = reg["rows"]; rho = reg["correlation"]["spearman_pIC50_vs_dock_score"]
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.05, 6.8 / 2.54 * 1.05))
    for s in ("strong", "medium", "weak"):
        xs = [r["dock_score"] for r in R if r.get("stratum") == s]
        ys = [r["pIC50"] for r in R if r.get("stratum") == s]
        ax.scatter(xs, ys, s=52, color=COL[s], edgecolor="white", linewidth=0.9,
                   label=f"{s} (n={len(xs)})", zorder=3)
    ax.set_xlabel("docking score (kcal/mol, lower = better predicted binding)", fontsize=8)
    ax.set_ylabel("pIC50 (measured)", fontsize=8)
    ax.set_title(f"Docking score does not rank measured potency  —  Spearman {rho:+.3f}  (n={reg['n']})",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.legend(fontsize=6.8, frameon=False, loc="upper left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.annotate("a perfect predictor would fall on a tight downward line",
                xy=(0.98, 0.04), xycoords="axes fraction", ha="right",
                fontsize=6.4, color=OI["grey"])
    fig.tight_layout()
    p = OUT / "fig10_score_vs_potency.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, reg["n"], "선택 자세 점수 대 실측 pIC50")
    return p


def fig11_regression(reg, src):
    m = reg["models"]
    names = ["single Vina score", f"{len(reg['terms'])} terms"]
    keys = ["single_vina_score", "five_terms"]
    fit = [m[k]["R2_fit"] for k in keys]
    q2 = [m[k]["Q2_loo"] for k in keys]
    nul = [m[k]["null_R2_median"] for k in keys]
    p95 = [m[k]["null_R2_p95"] for k in keys]
    x = range(len(names)); w = 0.34
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.05, 6.4 / 2.54 * 1.05))
    ax.bar([i - w / 2 for i in x], fit, w, label="R² (fit)", color=OI["sky"],
           edgecolor="white", linewidth=1.0)
    ax.bar([i + w / 2 for i in x], q2, w, label="Q² (leave-one-out)", color=OI["blue"],
           edgecolor="white", linewidth=1.0)
    for i, (f, q) in enumerate(zip(fit, q2)):
        ax.text(i - w / 2, f + (0.02 if f >= 0 else -0.05), f"{f:+.3f}", ha="center",
                fontsize=7.5, fontweight="bold", color="#1C242B")
        ax.text(i + w / 2, q + (0.02 if q >= 0 else -0.05), f"{q:+.3f}", ha="center",
                fontsize=7.5, fontweight="bold",
                color=OI["vermil"] if q < 0 else "#1C242B")
    for i, (nm, p9) in enumerate(zip(nul, p95)):
        ax.hlines(nm, i - 0.45, i + 0.45, color=OI["grey"], linestyle="--", linewidth=1.2)
        ax.hlines(p9, i - 0.45, i + 0.45, color=OI["vermil"], linestyle=":", linewidth=1.2)
    ax.axhline(0, color="black", linewidth=1.0)
    ax.set_xticks(list(x)); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("explained variance", fontsize=8)
    ax.set_title("Adding terms raises nothing and breaks cross-validation",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.annotate("dashed = median R² from label shuffling\ndotted = its 95th percentile",
                xy=(0.98, 0.96), xycoords="axes fraction", ha="right", va="top",
                fontsize=6.4, color=OI["grey"])
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = OUT / "fig11_regression.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, reg["n"], "R2/Q2 및 라벨섞기 귀무분포")
    return p


def fig12_control(dock, src):
    c = dock["result"]["control_redock"]
    modes = [r for r in c.get("rmsd_all_modes", []) if r is not None]
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.05, 5.6 / 2.54 * 1.05))
    cols = [OI["green"] if r <= c["threshold_angstrom"] else OI["vermil"] for r in modes]
    ax.bar(range(1, len(modes) + 1), modes, color=cols, edgecolor="white", linewidth=1.0)
    for i, r in enumerate(modes, start=1):
        ax.text(i, r + 0.06, f"{r:.2f}", ha="center", fontsize=7, color="#333")
    ax.axhline(c["threshold_angstrom"], color=OI["blue"], linewidth=1.8)
    ax.annotate(f"{c['threshold_angstrom']} Å criterion", xy=(0.99, c["threshold_angstrom"]),
                xycoords=("axes fraction", "data"), xytext=(0, 5),
                textcoords="offset points", ha="right", fontsize=6.8,
                color=OI["blue"], fontweight="bold")
    ax.set_xlabel("docking mode, ranked by score (1 = best score)", fontsize=8)
    ax.set_ylabel("heavy-atom RMSD to crystal pose (Å)", fontsize=8)
    ax.set_title("Re-docking control: the crystal pose is found, but not ranked first",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = OUT / "fig12_redock_control.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(modes), "모드별 결정 자세 대비 RMSD")
    return p


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    reg = json.loads((SR / "regression.json").read_text())
    dock = json.loads((SR / "docking.json").read_text())
    made = [fig10_score_vs_potency(reg, SR / "regression.json"),
            fig11_regression(reg, SR / "regression.json"),
            fig12_control(dock, SR / "docking.json")]
    for p in made:
        print(f"  {p.name}  ({p.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
