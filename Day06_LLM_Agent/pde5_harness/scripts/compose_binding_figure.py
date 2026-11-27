"""PyMOL 패널 4장을 보고서용 2x2 그림으로 합성한다.

· 수치는 전부 sample_run/docking_controlled.json 에서 읽는다 (하드코딩 금지)
· 잔기/금속 라벨은 PyMOL 마커 렌더로 얻은 anchors.json 좌표에 그린다
· 해상도는 보고서의 다른 그림과 같은 300 DPI, 폭 약 3,430 px
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

PANELS = os.environ.get("PANEL_DIR",
                        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "..", "sample_run", "report", "_pymol_panels"))
SR = os.environ.get("SAMPLE_RUN",
                    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "sample_run"))
OUTDIR = f"{SR}/report/figures_controlled"

for w in ("Regular", "SemiBold", "Bold"):
    p = f"/home/hjpark/.fonts/pretendard/Pretendard-{w}.ttf"
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
plt.rcParams["font.family"] = "Pretendard"
plt.rcParams["axes.unicode_minus"] = False

C_REF, C_TOP, C_BEST = "#009E73", "#D55E00", "#0072B2"
INK, SUB = "#1A2733", "#55606B"

d = json.load(open(f"{SR}/docking_controlled.json"))["result"]
cr = d["control_redock"]
rows = {r["chembl_id"]: r for r in d["rows"]}
TOP_ID, NEAR_ID = "CHEMBL536568", "CHEMBL47108"
rt, rn = rows[TOP_ID], rows[NEAR_ID]
anchors = json.load(open(f"{PANELS}/anchors.json"))

n_modes_top = rt["all_mode_mcs_rmsd"]
best_mode_of_top = min(range(len(n_modes_top)), key=lambda i: n_modes_top[i])

panels = [
    ("panelA_1udt.png", "A",
     "결합 부위 — 실데나필이 PDE5A 포켓에 붙는 방식",
     "공결정 리간드 VIA(녹색) · 4.5 Å 이내 잔기 17개(회색)\n"
     "빨간 점선 = Gln817 수소결합 · 구 = 촉매 금속 Zn·Mg · 물 제거"),
    ("panelB_redock.png", "B",
     "공결정 리간드 재도킹 — 샘플링은 맞고 채점이 틀렸다",
     f"공결정 자세(녹색) · 점수 1위 자세(주황, {cr['rmsd_top_pose_angstrom']} Å) · "
     f"최소 RMSD 모드 {cr['best_mode_index']}번(파랑, {cr['rmsd_best_of_modes_angstrom']} Å)\n"
     f"판정 기준 {cr['threshold_angstrom']} Å"),
    ("panelC_top_scorer.png", "C",
     f"점수 1위 화합물 {TOP_ID}",
     f"점수 {rt['top_pose_score']} kcal/mol (163종 중 1위) · "
     f"공결정 기준 MCS-RMSD {rt['top_pose_mcs_rmsd']} Å · "
     f"공유 부분구조 {rt['mcs_atoms']}원자 · {rt['similarity_bin']}/{rt['potency_bin']} 칸"),
    ("panelD_closest.png", "D",
     f"공결정에 가장 가까운 자세 {NEAR_ID}",
     f"점수 {rn['top_pose_score']} kcal/mol · "
     f"MCS-RMSD {rn['top_pose_mcs_rmsd']} Å (163종 최소) · "
     f"공유 부분구조 {rn['mcs_atoms']}원자 · {rn['similarity_bin']}/{rn['potency_bin']} 칸"),
]

# 네 패널의 내용 영역 합집합으로 공통 크롭 상자를 만든다 — 모든 패널에 같은
# 상자를 적용하므로 배율과 정렬은 그대로 유지된다.
import numpy as _np
from PIL import Image as _Image

_boxes = []
for fn, *_ in panels:
    _a = _np.asarray(_Image.open(f"{PANELS}/{fn}").convert("RGB")).astype(int)
    _nonwhite = (_a.sum(axis=2) < 735)        # 255*3 - 여유
    _ys, _xs = _np.nonzero(_nonwhite)
    if len(_xs):
        _boxes.append((_xs.min(), _ys.min(), _xs.max(), _ys.max()))
_H_img, _W_img = _np.asarray(
    _Image.open(f"{PANELS}/{panels[0][0]}")).shape[:2]
_x0 = max(min(b[0] for b in _boxes) - 40, 0)
_y0 = max(min(b[1] for b in _boxes) - 40, 0)
_x1 = min(max(b[2] for b in _boxes) + 40, _W_img)
_y1 = min(max(b[3] for b in _boxes) + 40, _H_img)
CROP = (_x0, _y0, _x1, _y1)
print("공통 크롭:", CROP, "->", (_x1 - _x0, _y1 - _y0))


def panel_image(fn):
    return _np.asarray(_Image.open(f"{PANELS}/{fn}").convert("RGB").crop(CROP))


def to_crop(nx, ny):
    """원본 정규화 좌표 -> 크롭된 이미지 기준 정규화 좌표."""
    px, py = nx * _W_img, ny * _H_img
    return ((px - CROP[0]) / (CROP[2] - CROP[0]),
            (py - CROP[1]) / (CROP[3] - CROP[1]))


FIG_W = 11.43          # × 300 dpi = 3,429 px — 보고서의 다른 그림과 같은 폭
_ar = (CROP[3] - CROP[1]) / (CROP[2] - CROP[0])
_fig_h = 2 * (FIG_W * 0.49 * _ar) + 2.05
fig, axes = plt.subplots(2, 2, figsize=(FIG_W, _fig_h), dpi=300)
fig.patch.set_facecolor("white")

for idx, (ax, (fn, letter, title, sub)) in enumerate(zip(axes.ravel(), panels)):
    ax.imshow(panel_image(fn))
    ax.axis("off")
    ax.text(-0.015, 1.105, letter, transform=ax.transAxes,
            fontsize=18, fontweight="bold", color=INK, ha="left", va="top")
    ax.text(0.040, 1.105, title, transform=ax.transAxes,
            fontsize=11.4, fontweight="bold", color=INK, ha="left", va="top")
    ax.text(0.040, 1.040, sub, transform=ax.transAxes,
            fontsize=8.4, color=SUB, ha="left", va="top", linespacing=1.45)

    if letter == "A":
        import math
        placed = []
        for name, (nx0, ny0) in sorted(anchors.items()):
            nx, ny = to_crop(nx0, ny0)
            ax_x, ax_y = nx, 1 - ny
            dx, dy = ax_x - 0.5, ax_y - 0.5
            r = math.hypot(dx, dy) or 1e-6
            tx = min(max(ax_x + dx / r * 0.19, 0.05), 0.95)
            ty = min(max(ax_y + dy / r * 0.19, 0.05), 0.95)
            # 이미 놓인 라벨과 너무 가까우면 바깥으로 더 민다
            for _ in range(14):
                if all(math.hypot(tx - px, ty - py) > 0.115 for px, py in placed):
                    break
                tx = min(max(tx + dx / r * 0.055, 0.05), 0.95)
                ty = min(max(ty + dy / r * 0.055, 0.05), 0.95)
            placed.append((tx, ty))
            metal = name in ("Zn", "Mg")
            ax.annotate(
                name, xy=(ax_x, ax_y), xycoords=ax.transAxes,
                xytext=(tx, ty), textcoords=ax.transAxes,
                fontsize=8.6, color="#8A5A00" if metal else INK,
                fontweight="semibold", ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.7,
                                color="#8893A0", shrinkA=1, shrinkB=2),
                bbox=dict(boxstyle="round,pad=0.22", fc="white",
                          ec="#C9D1D9", lw=0.5, alpha=0.94))

handles = [
    Line2D([0], [0], color=C_REF, lw=6, label="공결정 리간드 (1UDT VIA, 실험 결정 자세)"),
    Line2D([0], [0], color=C_TOP, lw=6, label="점수 1위 자세 (채점 함수가 고른 자세)"),
    Line2D([0], [0], color=C_BEST, lw=6, label="공결정에 가장 가까운 자세"),
]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
           fontsize=9.4, bbox_to_anchor=(0.5, 0.030),
           handlelength=1.5, columnspacing=2.0)

fig.text(0.5, 0.008,
         "네 패널은 같은 방향·같은 배율이다 — 촉매 도메인 N-말단(537번)이 위, 리간드 포켓이 정면. "
         "오픈소스 PyMOL 3.2.0a 렌더링 · 수용체 PDB 1UDT · "
         "도킹 smina (AutoDock Vina 1.1.2 채점), exhaustiveness 64, seed 42",
         ha="center", va="bottom", fontsize=7.6, color="#78838E")

fig.subplots_adjust(left=0.008, right=0.992, top=0.935, bottom=0.075,
                    wspace=0.02, hspace=0.155)

os.makedirs(OUTDIR, exist_ok=True)
for ext in ("png", "svg"):
    out = f"{OUTDIR}/fig14_pymol_poses.{ext}"
    fig.savefig(out, dpi=300, facecolor="white")
    print("saved", out, os.path.getsize(out) // 1024, "KB")

from PIL import Image
im = Image.open(f"{OUTDIR}/fig14_pymol_poses.png")
print("크기:", im.size, "dpi:", im.info.get("dpi"))

meta = {
    "figure": "fig14_pymol_poses",
    "renderer": "PyMOL open-source 3.2.0a (pymol-open-source==3.2.0a0.post5859, PyPI)",
    "receptor": "PDB 1UDT chain A, waters removed, ZN/MG retained",
    "reference_ligand": "VIA (sildenafil, co-crystal)",
    "orientation": "N-terminus (residue 537) up, ligand pocket facing viewer; "
                   "camera axes computed from coordinates, identical in all panels",
    "pocket_definition": "byres (receptor within 4.5 A of VIA) — 17 residues",
    "labels": "anchors measured by re-rendering each residue as a uniquely coloured "
              "sphere and taking that colour's pixel centroid (no projection math)",
    "panels": {
        "A": "co-crystal binding mode with Gln817 hydrogen bonds",
        "B": f"reference redock: top pose {cr['rmsd_top_pose_angstrom']} A, "
             f"best mode {cr['best_mode_index']} {cr['rmsd_best_of_modes_angstrom']} A",
        "C": f"{TOP_ID} top-scoring pose, score {rt['top_pose_score']}, "
             f"MCS-RMSD {rt['top_pose_mcs_rmsd']} A",
        "D": f"{NEAR_ID} top pose, score {rn['top_pose_score']}, "
             f"MCS-RMSD {rn['top_pose_mcs_rmsd']} A (minimum of 163)",
    },
    "source_data": "sample_run/docking_controlled.json",
    "anchors_normalised": anchors,
}
with open(f"{OUTDIR}/fig14_pymol_poses.meta.json", "w", encoding="utf-8") as f:
    json.dump(meta, f, ensure_ascii=False, indent=2)
print("saved meta.json")
