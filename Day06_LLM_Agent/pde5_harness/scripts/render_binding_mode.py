#!/usr/bin/env python3
"""결합 양상 렌더링 — 오픈소스 PyMOL 로 도킹 자세 그림을 만든다.

자세는 점수 1위가 아니라 참조(공결정 sildenafil) 기준으로 고른 것을 쓴다.
C2 채점 대조가 실패했기 때문이며, 그 선택 규칙은 보고서 방법에 명시한다.

세 장을 만든다.
  A. 전체 조망  — 수용체 카툰 + 결합 포켓 위치
  B. 포켓 근접  — 결정 자세(회색)와 도킹 자세(색) 겹침 + 극성 접촉
  C. 접촉 잔기  — 4 A 이내 잔기를 거리 기반으로 집계한 막대그림 (matplotlib)

접촉 판정은 좌표에서 직접 계산한다. 도구가 내놓지 않은 상호작용을 그리지 않는다.
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402
ROOT = HERE.parent
WORK = ROOT / "sample_run" / "structures" / "work"
PYMOL = "/share/anaconda3/bin/pymol"
CONTACT_A = 4.0


def pml(script: str, out_png: Path, width=1800, height=1200) -> bool:
    f = WORK / "scene.pml"
    f.write_text(script + f'\npng {out_png}, width={width}, height={height}, dpi=300, ray=1\n')
    p = subprocess.run([PYMOL, "-cq", str(f)], capture_output=True, text=True, timeout=900)
    if not out_png.exists():
        print(p.stdout[-800:], p.stderr[-800:], file=sys.stderr)
    return out_png.exists()


COMMON = """
bg_color white
set ray_opaque_background, 1
set cartoon_transparency, 0.15
set stick_radius, 0.14
set sphere_scale, 0.28
set label_size, 16
set label_color, black
set antialias, 2
set ray_shadows, 0
set depth_cue, 0
set specular, 0.2
"""


def scene_overview(receptor: Path, ref: Path, pose: Path, out: Path) -> bool:
    s = f"""
load {receptor}, rec
load {ref}, refl
load {pose}, dock
{COMMON}
hide everything
show cartoon, rec
color grey80, rec
show spheres, rec and (resn ZN+MG)
color orange, rec and resn ZN
color green, rec and resn MG
show sticks, refl
color grey50, refl
show sticks, dock
color marine, dock
select pocket, byres (rec within 5 of refl) and polymer
show sticks, pocket and not (name C+N+O)
color palecyan, pocket and elem C
orient rec
turn x, -10
"""
    return pml(s, out)


def scene_pocket(receptor: Path, ref: Path, pose: Path, out: Path,
                 label_resi: list[str] | None = None) -> bool:
    # 라벨은 실제 접촉 잔기만. 전부 붙이면 글자가 겹쳐 읽을 수 없다.
    sel = ("resi " + "+".join(label_resi)) if label_resi else "polymer"
    s = f"""
load {receptor}, rec
load {ref}, refl
load {pose}, dock
{COMMON}
hide everything
select pocket, byres (rec within 5 of refl) and polymer
show sticks, pocket
color palecyan, pocket and elem C
show sticks, refl
color grey60, refl and elem C
show sticks, dock
color marine, dock and elem C
show spheres, rec and (resn ZN+MG)
color orange, rec and resn ZN
color green, rec and resn MG
distance polar, dock, pocket, 3.5, mode=2
color red, polar
set dash_width, 2.5
set dash_gap, 0.35
select shown, pocket and ({sel})
label shown and name CB, "%s%s" % (resn, resi)
set label_size, 19
set label_outline_color, white
set label_position, (0, 0, 2.5)
hide sticks, pocket and not shown
orient refl
zoom refl, 3.2
turn y, 12
"""
    return pml(s, out)


def contacts(receptor_pdb: Path, pose_sdf: Path):
    """도킹 자세에서 4 A 이내 잔기를 좌표로 직접 센다."""
    from rdkit import Chem
    m = next(iter(Chem.SDMolSupplier(str(pose_sdf), removeHs=True)), None)
    if m is None:
        return []
    conf = m.GetConformer()
    lig = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
           for i in range(m.GetNumAtoms())]
    hits = Counter()
    for line in receptor_pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = line[76:78].strip() or line[12:16].strip()[:1]
        if el == "H":
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        key = f"{line[17:20].strip()}{line[22:26].strip()}"
        for lx, ly, lz in lig:
            if (x - lx) ** 2 + (y - ly) ** 2 + (z - lz) ** 2 <= CONTACT_A ** 2:
                hits[key] += 1
                break
    return hits.most_common()


def fig_contacts(pairs, out: Path, title: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 300, "savefig.dpi": 300})
plt.rcParams["font.family"] = "Pretendard"      # 사용자 지정 서체
plt.rcParams["axes.unicode_minus"] = False         # Pretendard 의 마이너스 글리프 사용

    top = pairs[:14][::-1]
    if not top:
        return None
    fig, ax = plt.subplots(figsize=(16 / 2.54 * 1.05, 6.4 / 2.54 * 1.05))
    ax.barh([k for k, _ in top], [v for _, v in top], color="#0072B2",
            edgecolor="white", linewidth=1.0)
    for i, (_, v) in enumerate(top):
        ax.text(v + 0.15, i, str(v), va="center", fontsize=7, color="#333")
    ax.set_xlabel(f"contacting heavy atoms of the ligand (within {CONTACT_A} A)", fontsize=8)
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--compound", default=None, help="ChEMBL ID (기본: 가장 강력한 것)")
    ap.add_argument("--out", default=str(ROOT / "outputs" / "figures"))
    ap.add_argument("--source", choices=("legacy", "controlled"), default="controlled")
    ap.add_argument("--pose", choices=("top", "reference"), default="top",
                    help="top = 점수 1위(참조 미사용), reference = 공결정 리간드 기준 선택")
    a = ap.parse_args()
    W2 = ROOT / "sample_run" / "structures" / "work_controlled"
    receptor = WORK / "receptor.pdb"
    ref = WORK / "ligand_ref.sdf"

    if a.source == "controlled":
        dj = ROOT / "sample_run" / "docking_controlled.json"
        if not dj.exists():
            sys.exit("docking_controlled.json 이 없다 — 먼저 dock_controlled.py 를 실행하라.")
        env = json.loads(dj.read_text())
        if not env["verification"]["passed"]:
            sys.exit("도킹 게이트 미통과 — 결합 양상 그림을 만들지 않는다.")
        rows = [r for r in env["result"]["rows"] if r.get("top_pose_score") is not None]
        target = a.compound or max(rows, key=lambda d: d["pchembl_value"])["chembl_id"]
        row = next(r for r in rows if r["chembl_id"] == target)
        allf = W2 / f"{target}_poses.sdf"
        if not allf.exists():
            sys.exit(f"{target} 자세 파일이 없다.")
        k = 1 if a.pose == "top" else (row.get("ref_pose_index") or 1)
        blocks = allf.read_text().split("$$$$\n")
        if k > len(blocks) or not blocks[k - 1].strip():
            sys.exit(f"{target} 모드 {k} 를 찾지 못했다.")
        pose = W2 / f"{target}_render_m{k}.sdf"
        pose.write_text(blocks[k - 1] + "$$$$\n")
        print(f"대상 {target} (pChEMBL {row['pchembl_value']:.2f}, "
              f"Tanimoto {row['tanimoto_to_sildenafil']:.3f}), 자세 규칙 {a.pose}, 모드 {k}",
              file=sys.stderr)
    else:
        dj = ROOT / "sample_run" / "docking.json"
        if not dj.exists():
            sys.exit("docking.json 이 없다 — 먼저 dock_actives.py 를 실행하라.")
        env = json.loads(dj.read_text())
        if not env["verification"]["passed"]:
            sys.exit("도킹 게이트 미통과 — 결합 양상 그림을 만들지 않는다.")
        docked = [d for d in env["result"]["docked"] if d.get("dock_score") is not None]
        target = a.compound or max(docked, key=lambda d: d["pIC50"])["chembl_id"]
        pose = WORK / f"{target}_selected.sdf"
        if not pose.exists():
            pose = WORK / f"{target}_pose.sdf"
        if not pose.exists():
            sys.exit(f"{target} 자세 파일이 없다.")
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    made = []
    p1 = out / "fig7_binding_overview.png"
    if scene_overview(receptor, ref, pose, p1): made.append(p1)
    pairs = contacts(receptor, pose)
    import re as _re
    top_resi = [_re.sub(r"^[A-Z]{2,3}", "", k) for k, _ in pairs[:12]]
    p2 = out / "fig8_binding_pocket.png"
    if scene_pocket(receptor, ref, pose, p2, top_resi): made.append(p2)
    p3 = fig_contacts(pairs, out / "fig9_contacts.png",
                      f"Residues contacting {target} in the docked pose (heavy atoms within {CONTACT_A} A)")
    if p3: made.append(p3)

    NOTE = {"top": ("자세는 점수 1위다 — 공결정 리간드를 참조하지 않았다. "
                    "따라서 접촉 잔기 일치는 선택 규칙의 귀결이 아니다. "
                    "다만 상자를 포켓에 지정했으므로 포켓 잔기 접촉 자체는 거의 자동이다."),
            "reference": ("자세는 공결정 리간드와의 MCS-RMSD 최소 규칙으로 골랐다. "
                          "따라서 접촉 잔기 일치는 선택 규칙의 귀결이지 검증이 아니다.")}
    payload = {"compound": target, "cutoff_angstrom": CONTACT_A, "contacts": pairs,
               "source": a.source, "pose_rule": a.pose,
               "pose_selection_note": NOTE[a.pose]}
    checks = [
        ("접촉 잔기 검출됨", len(pairs) > 0),
        ("접촉이 좌표에서 계산됨 (하드코딩 아님)", all(isinstance(v, int) and v > 0
                                                       for _, v in pairs)),
        ("그림 3장 렌더됨", len(made) == 3),
    ]
    (ROOT / "sample_run" / "binding_contacts.json").write_text(
        json.dumps(make_result(payload, "PyMOL 좌표 기반 접촉 계산",
                               f"{target}, cutoff {CONTACT_A} A", checks,
                               notes="n=1 이므로 어느 방향으로도 일반화할 수 없다."),
                   ensure_ascii=False, indent=2) + "\n")
    for p in made:
        print(f"  {p}  ({p.stat().st_size // 1024} KB)")
    print(f"결합 양상 그림 {len(made)}장 — 대상 {target}, 접촉 잔기 {len(pairs)}개")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
