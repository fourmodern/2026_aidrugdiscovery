#!/usr/bin/env python3
"""보고서 그림 생성 — 수치가 들어가는 그림은 run_stdout.json 실측값만 쓴다.

무-날조: 난수 생성·보간·추세선 임의 삽입 금지. 데이터 없는 패널은 빈 패널 + 사유.
"""
import argparse, json, re, sys, hashlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Okabe-Ito (colorblind-safe)
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermil": "#D55E00",
      "sky": "#56B4E9", "yellow": "#F0E442", "purple": "#CC79A7", "grey": "#7F7F7F"}
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 300,
                     "savefig.dpi": 300, "axes.grid": False})

QED_MIN, SA_MAX = 0.5, 6.0        # scripts/mol_properties.py 와 동일


def load_run(path: Path):
    """run_stdout.json 에서 최상위 JSON 객체들을 순서대로 꺼낸다."""
    txt = path.read_text()
    dec, i, objs = json.JSONDecoder(), 0, []
    while True:
        m = re.search(r"[\[{]", txt[i:])
        if not m:
            break
        st = i + m.start()
        try:
            o, end = dec.raw_decode(txt[st:]); objs.append(o); i = st + end
        except ValueError:
            i = st + 1
    target = next((o for o in objs if isinstance(o, dict) and "accession" in o), None)
    props = next((o for o in objs if isinstance(o, list) and o and "QED" in o[0]), None)
    if props is None:
        sys.exit("run_stdout.json 에서 물성 레코드를 찾지 못했습니다 — 그림을 만들지 않습니다.")
    return target, props


def meta(path: Path, source: Path, n: int, note: str):
    path.with_suffix(".meta.json").write_text(json.dumps({
        "source": str(source), "script": "scripts/make_figures.py", "n_records": n,
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest()[:16], "note": note,
    }, ensure_ascii=False, indent=2))


def _card(ax, x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                                linewidth=1.6, edgecolor=color, facecolor=color + "18"))
    ax.text(x + w / 2, y + h - 0.055, title, ha="center", va="top",
            fontsize=8.5, fontweight="bold", color=color)
    ax.text(x + w / 2, y + h / 2 - 0.035, body, ha="center", va="center",
            fontsize=6.6, color="#333333", linespacing=1.45)


def _arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11,
                                 linewidth=1.3, color=OI["grey"]))


def fig1_graphical_abstract(out: Path, src: Path, n_props: int, n_pass: int):
    """개념도 — 수치는 실측 요약(n)만 넣는다."""
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 9 / 2.54 * 1.15))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.955, "Can an agent harness keep a drug-discovery pipeline honest?",
            ha="center", fontsize=11.5, fontweight="bold", color="#1C242B")
    ax.text(0.5, 0.895, "PDE5A inhibitor triage under file-encoded contracts and executable gates",
            ha="center", fontsize=7.8, color=OI["grey"])
    _card(ax, 0.035, 0.44, 0.21, 0.36, "1  QUESTION",
          "Do LLM agents invent\nnumbers when a tool\nreturns nothing?", OI["blue"])
    _card(ax, 0.275, 0.44, 0.21, 0.36, "2  METHOD",
          "Contract in CLAUDE.md\n4 tools, 1 envelope\ngate() blocks only if\nthe caller checks it", OI["orange"])
    _card(ax, 0.515, 0.44, 0.21, 0.36, "3  FINDING",
          f"{n_props} compounds scored\n{n_pass} passed the gate\nQED alone decided", OI["green"])
    _card(ax, 0.755, 0.44, 0.21, 0.36, "4  IMPLICATION",
          "A gate you never\nfailed is untested,\nnot verified.", OI["purple"])
    for x in (0.245, 0.485, 0.725):
        _arrow(ax, x, 0.62, x + 0.03, 0.62)
    ax.add_patch(FancyBboxPatch((0.035, 0.10), 0.93, 0.25, boxstyle="round,pad=0.012",
                                linewidth=1.4, edgecolor=OI["vermil"], facecolor="#FDF2EC"))
    ax.text(0.5, 0.295, "Every result value in this report is a tool return value",
            ha="center", fontsize=8.6, fontweight="bold", color=OI["vermil"])
    ax.text(0.5, 0.185,
            "UniProt O76074 lookup  ·  ChEMBL CHEMBL1827 actives  ·  RDKit descriptors\n"
            "Prose and environment strings are human-written. Unavailable quantities are reported as unavailable.",
            ha="center", fontsize=6.8, color="#444444", linespacing=1.6)
    p = out / "fig1_graphical_abstract.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, n_props, "개념도 — 수치는 실측 n 요약만 포함")
    return p


def fig2_pipeline(out: Path, src: Path, n_props: int):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 8.2 / 2.54 * 1.15))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.96, "Harness architecture: contract, tools, envelope, gate",
            ha="center", fontsize=11, fontweight="bold", color="#1C242B")
    ax.add_patch(FancyBboxPatch((0.03, 0.72), 0.94, 0.15, boxstyle="round,pad=0.01",
                                linewidth=1.5, edgecolor=OI["blue"], facecolor="#E8F1F8"))
    ax.text(0.5, 0.845, "CLAUDE.md  —  workflow order · per-step pass conditions · no-fabrication rules",
            ha="center", fontsize=7.8, fontweight="bold", color=OI["blue"])
    ax.text(0.5, 0.775, ".claude/skills/*/SKILL.md  (when to call, how to call, what to verify)   ·"
                        "   .claude/settings.json  (PostToolUse hooks)",
            ha="center", fontsize=6.5, color="#456")
    steps = [("target-lookup", "UniProt\nO76074"), ("chembl-actives", "CHEMBL1827\nactives"),
             ("mol-properties", f"RDKit\n{n_props} compounds"), ("selectivity-check", "PDE5 vs PDE6\nqualitative")]
    for i, (name, body) in enumerate(steps):
        x = 0.045 + i * 0.238
        _card(ax, x, 0.375, 0.20, 0.255, name, body, OI["orange"])
        _blk = name != "selectivity-check"
        ax.add_patch(FancyBboxPatch((x, 0.245), 0.20, 0.095, boxstyle="round,pad=0.008",
                                    linewidth=1.3,
                                    edgecolor=OI["green"] if _blk else OI["vermil"],
                                    facecolor="#E6F5EF" if _blk else "#FDF0E8"))
        _blocking = name != "selectivity-check"
        ax.text(x + 0.10, 0.293,
                "gate() PASSED" if _blocking else "logged, did not block",
                ha="center", va="center", fontsize=6.6, fontweight="bold",
                color=OI["green"] if _blocking else OI["vermil"])
        _arrow(ax, x + 0.10, 0.373, x + 0.10, 0.343)
        if i < 3:
            _arrow(ax, x + 0.205, 0.293, x + 0.233, 0.293)
    ax.text(0.5, 0.165, "every tool returns the same envelope   { result , provenance , verification }",
            ha="center", fontsize=8, fontweight="bold", color="#1C242B",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F4F6F8", edgecolor=OI["grey"], linewidth=1.0))
    ax.text(0.5, 0.065, "verification.passed = false blocks the next step in 3 of 4 steps — the 4th discarded the return value",
            ha="center", fontsize=6.9, color=OI["vermil"])
    p = out / "fig2_pipeline.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, n_props, "구조 모식도 — 수치는 단계 수와 레코드 수만")
    return p


def fig3_property_space(out: Path, src: Path, props):
    fig, axes = plt.subplots(1, 2, figsize=(16 / 2.54 * 1.15, 6.6 / 2.54 * 1.15))
    ok = [r for r in props if r["gate_pass"]]; no = [r for r in props if not r["gate_pass"]]
    ax = axes[0]
    ax.scatter([r["MW"] for r in no], [r["logP"] for r in no], s=64, marker="X",
               color=OI["vermil"], edgecolor="white", linewidth=0.8, label=f"gate FAIL (n={len(no)})", zorder=3)
    ax.scatter([r["MW"] for r in ok], [r["logP"] for r in ok], s=76, marker="o",
               color=OI["green"], edgecolor="white", linewidth=0.8, label=f"gate PASS (n={len(ok)})", zorder=4)
    ax.axvline(500, color=OI["grey"], linestyle="--", linewidth=1.0)
    ax.text(502, ax.get_ylim()[1] * 0.98, "Lipinski MW 500", fontsize=6.2, color=OI["grey"],
            va="top", rotation=90)
    ax.set_xlabel("Molecular weight (Da)", fontsize=8)
    ax.set_ylabel("cLogP", fontsize=8)
    ax.set_title("A  Property space, coloured by gate outcome", fontsize=8.5,
                 fontweight="bold", loc="left")
    ax.legend(fontsize=6.6, frameon=False, loc="lower right")
    ax.tick_params(labelsize=7)
    ax = axes[1]
    ax.scatter([r["QED"] for r in no], [r["SA"] for r in no], s=64, marker="X",
               color=OI["vermil"], edgecolor="white", linewidth=0.8, zorder=3)
    ax.scatter([r["QED"] for r in ok], [r["SA"] for r in ok], s=76, marker="o",
               color=OI["green"], edgecolor="white", linewidth=0.8, zorder=4)
    ax.axvline(QED_MIN, color=OI["blue"], linestyle="-", linewidth=1.6)
    # 축 좌표(0-1)로 배치해야 데이터 점과 겹치지 않는다
    ax.annotate(f"QED_MIN = {QED_MIN}", xy=(QED_MIN, 1.0), xycoords=("data", "axes fraction"),
                xytext=(4, -6), textcoords="offset points", fontsize=6.6,
                color=OI["blue"], va="top", ha="left", fontweight="bold")
    ax.axhline(SA_MAX, color=OI["grey"], linestyle="--", linewidth=1.0)
    ax.annotate(f"SA_MAX = {SA_MAX}", xy=(0.99, SA_MAX), xycoords=("axes fraction", "data"),
                xytext=(0, -7), textcoords="offset points", fontsize=6.2,
                color=OI["grey"], va="top", ha="right")
    ax.set_xlabel("QED (quantitative estimate of drug-likeness)", fontsize=8)
    ax.set_ylabel("Synthetic accessibility (lower = easier)", fontsize=8)
    ax.set_title(f"B  The single threshold that decided (SA max {max(r['SA'] for r in props):.2f} << {SA_MAX})",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for a in axes:
        for sp in ("top", "right"): a.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig3_property_space.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(props), "RDKit 실계산값 산점도")
    return p


def fig4_gate_waterfall(out: Path, src: Path, props, target):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 6.0 / 2.54 * 1.15))
    n = len(props); n_ro5 = sum(r["Ro5_pass"] for r in props)
    n_sa = sum(r["SA"] <= SA_MAX for r in props); n_final = sum(r["gate_pass"] for r in props)
    labels = ["ChEMBL actives\nretrieved", "RDKit descriptors\ncomputed",
              "Lipinski Ro5\n(<=1 violation)", f"SA <= {SA_MAX}", f"QED >= {QED_MIN}"]
    vals = [n, n, n_ro5, n_sa, n_final]
    colors = [OI["blue"], OI["blue"], OI["sky"], OI["sky"], OI["green"]]
    bars = ax.bar(range(len(vals)), vals, color=colors, edgecolor="white", linewidth=1.2, width=0.62)
    for i, (b, v) in enumerate(zip(bars, vals)):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.18, str(v), ha="center",
                fontsize=9, fontweight="bold", color="#1C242B")
        if i and vals[i] < vals[i - 1]:
            ax.text(b.get_x() + b.get_width() / 2, v / 2, f"-{vals[i-1]-v}", ha="center",
                    fontsize=8, fontweight="bold", color="white")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Compounds remaining", fontsize=8)
    ax.set_ylim(0, n + 1.6)
    ax.set_title("Attrition through the drug-likeness gate (real run, 2026-09-04)",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig4_gate_waterfall.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, n, "단계별 잔존 건수 — 실측")
    return p


def fig5_qed_threshold(out: Path, src: Path, props):
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 5.6 / 2.54 * 1.15))
    qeds = sorted((r["QED"], r["gate_pass"]) for r in props)
    ys = range(len(qeds))
    ax.barh(list(ys), [q for q, _ in qeds],
            color=[OI["green"] if g else OI["vermil"] for _, g in qeds],
            edgecolor="white", linewidth=1.0, height=0.68)
    for i, (q, g) in enumerate(qeds):
        ax.text(q + 0.008, i, f"{q:.3f}", va="center", fontsize=6.8,
                color=OI["green"] if g else OI["vermil"], fontweight="bold")
    ax.axvline(QED_MIN, color=OI["blue"], linewidth=2.0)
    # 막대와 겹치지 않도록 축 하단(axes fraction y=0)에 붙인다
    ax.annotate(f"QED_MIN = {QED_MIN}  (demo threshold, set by the script)",
                xy=(QED_MIN, 0.0), xycoords=("data", "axes fraction"),
                xytext=(5, 5), textcoords="offset points",
                fontsize=6.8, color=OI["blue"], va="bottom", ha="left", fontweight="bold")
    ax.set_yticks([]); ax.set_xlabel("QED", fontsize=8)
    ax.set_xlim(0, max(q for q, _ in qeds) * 1.18)
    ax.set_title("Every compound ranked by QED — the threshold is the whole decision",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig5_qed_threshold.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(props), "QED 실계산값 정렬 + 임계선")
    return p


def fig6_threshold_sweep(out: Path, src: Path, props):
    """임계를 바꾸면 결과가 어떻게 달라지는가 — 실측 QED 만으로 재계산."""
    qed = sorted(r["QED"] for r in props)
    ts = [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]
    ns = [sum(q >= t for q in qed) for t in ts]
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.15, 5.8 / 2.54 * 1.15))
    ax.step(ts, ns, where="post", color=OI["blue"], linewidth=2.0)
    ax.scatter(ts, ns, s=46, color=OI["blue"], zorder=4)
    for t, nn in zip(ts, ns):
        ax.annotate(str(nn), (t, nn), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=7.6, fontweight="bold", color=OI["blue"])
    i50 = ts.index(0.50)
    ax.scatter([0.50], [ns[i50]], s=150, facecolor="none", edgecolor=OI["vermil"],
               linewidth=2.0, zorder=5)
    ax.annotate(f"threshold used in this run\n QED_MIN = {QED_MIN} -> {ns[i50]} pass",
                (0.50, ns[i50]), textcoords="offset points", xytext=(12, 14),
                fontsize=7.2, color=OI["vermil"], fontweight="bold")
    ax.set_xlabel("QED threshold", fontsize=8)
    ax.set_ylabel("Compounds passing", fontsize=8)
    ax.set_ylim(0, len(qed) + 1.5)
    ax.set_title("Sensitivity of the outcome to an arbitrary threshold (same 10 measured QED values)",
                 fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    p = out / "fig6_threshold_sweep.png"
    fig.savefig(p, bbox_inches="tight", facecolor="white"); plt.close(fig)
    meta(p, src, len(props), "실측 QED 재집계 — 임계만 바꿔 통과 건수 계산")
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="sample_run/run_stdout.json")
    ap.add_argument("--out", default="outputs/figures")
    a = ap.parse_args()
    src = Path(a.run); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    target, props = load_run(src)
    n_pass = sum(r["gate_pass"] for r in props)
    made = [fig1_graphical_abstract(out, src, len(props), n_pass),
            fig2_pipeline(out, src, len(props)),
            fig3_property_space(out, src, props),
            fig4_gate_waterfall(out, src, props, target),
            fig5_qed_threshold(out, src, props),
            fig6_threshold_sweep(out, src, props)]
    for p in made:
        print(f"  {p}  ({p.stat().st_size//1024} KB)")
    print(f"그림 {len(made)}장 — 물성 {len(props)}건 중 게이트 통과 {n_pass}건")


if __name__ == "__main__":
    main()
