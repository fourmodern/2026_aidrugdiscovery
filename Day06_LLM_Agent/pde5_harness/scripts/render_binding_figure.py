"""오픈소스 PyMOL 로 PDE5A 결합 구조를 설명하는 그림을 렌더링한다.

방향 규약 (좌표에서 계산, 네 패널 공통)
  · 촉매 도메인 N-말단(537번)이 화면 위
  · 리간드 결합 포켓이 정면

가독성 규약
  · 포켓 표면은 리간드보다 '뒤'에 있는 잔기에만 씌운다 → 앞벽이 리간드를 가리지 않는다
  · 카툰은 포켓 주변 12 Å 로 제한 (전체 접힘은 설명에 기여하지 않아 생략)
  · 깊이(slab)로 앞쪽 잡음을 잘라낸다

패널 A  결합 부위 — 실데나필과 포켓, 핵심 잔기와 수소결합
패널 B  공결정 리간드 재도킹 — 샘플링 vs 채점
패널 C  점수 1위 화합물 CHEMBL536568
패널 D  공결정에 가장 가까운 자세 CHEMBL47108

실행: SAMPLE_RUN=<...>/sample_run python -m pymol -cq render_binding_figure.py
      (오픈소스 PyMOL: pip install pymol-open-source)
"""
import os

import numpy as np
from pymol import cmd

SR = os.environ.get(
    "SAMPLE_RUN",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
        os.environ.get("HARNESS_DIR", __file__)))), "sample_run"))
# PyMOL 은 스크립트 실행 시 __file__ 을 자기 패키지로 덮으므로 경로는 환경변수로 받는다.
if not os.path.isdir(SR):
    SR = "/home/hjpark/2026_aidrugdiscovery/Day06_LLM_Agent/pde5_harness/sample_run"
OUT = os.environ.get("PANEL_DIR",
                     os.path.join(SR, "report", "_pymol_panels"))
os.makedirs(OUT, exist_ok=True)

W, H = 1680, 1290
SLAB = float(os.environ.get("NB_SLAB", "34"))

# Okabe-Ito 계열 — 색맹 안전
C_REF = "0x009E73"   # 공결정 리간드
C_TOP = "0xD55E00"   # 점수 1위 자세
C_BEST = "0x0072B2"  # 공결정에 가장 가까운 자세
C_POCK = "0x51606E"  # 포켓 잔기
C_SURF = "0xB9C6D2"  # 포켓 표면

# 설명용으로 라벨을 다는 잔기 (PDE5 억제제 결합의 교과서적 기준점)
LABELLED = {"GLN817": "Gln817", "PHE820": "Phe820",
            "VAL782": "Val782", "TYR612": "Tyr612"}


def camera_axes():
    """N-말단 위 · 포켓 정면이 되는 카메라 축 (x, up, z) 를 좌표에서 계산한다."""
    prot = np.array(cmd.centerofmass("rec"), dtype=float)
    lig = np.array(cmd.centerofmass("via"), dtype=float)

    cas = []
    cmd.iterate_state(1, "rec and name CA",
                      "cas.append((int(resi), x, y, z))", space={"cas": cas})
    n_resi, nx, ny, nz = min(cas, key=lambda t: t[0])

    z = lig - prot
    z /= np.linalg.norm(z)
    up = np.array([nx, ny, nz], dtype=float) - prot
    up -= z * float(np.dot(up, z))
    up /= np.linalg.norm(up)
    x = np.cross(up, z)
    return x, up, z, n_resi, lig


def base_scene(label=False, contacts=False, surface=True):
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.26)
    cmd.set("specular", 0.18)
    cmd.set("auto_zoom", 0)
    cmd.set("two_sided_lighting", 1)
    cmd.set("ray_interior_color", "0x58636E")
    cmd.set("cartoon_transparency", 0.86)
    cmd.set("stick_radius", 0.155)
    cmd.set("sphere_scale", 0.30)
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.08)
    cmd.set("surface_color", C_SURF)
    cmd.set("label_size", 20)
    cmd.set("label_color", "0x1A2733")
    cmd.set("label_outline_color", "white")
    cmd.set("float_labels", 1)   # 표면/막대 뒤에 가리지 않게
    cmd.set("label_position", (0, 0, 2.2))
    cmd.set("dash_color", "0xCC2200")
    cmd.set("dash_width", 2.6)
    cmd.set("dash_gap", 0.32)
    cmd.set("dash_radius", 0.055)   # ray 렌더링에서 점선이 보이려면 반지름이 있어야 한다
    cmd.set("dash_length", 0.45)
    cmd.set("valence", 0)

    cmd.load(f"{SR}/structures/1UDT.pdb", "rec")
    cmd.remove("rec and solvent")
    cmd.remove("rec and hydro")
    cmd.create("via", "rec and resn VIA")
    cmd.create("metals", "rec and resn ZN+MG")
    cmd.remove("rec and resn VIA+ZN+MG")

    cmd.hide("everything")

    # 포켓 잔기 — 좌표에서 계산 (VIA 4.5 Å 이내)
    cmd.select("pocket", "byres (rec within 4.5 of via)")
    cmd.show("sticks", "pocket and sidechain")
    cmd.set("stick_radius", 0.135, "pocket")
    cmd.color(C_POCK, "pocket")
    cmd.util.cnc("pocket")

    # 맥락용 카툰은 포켓 주변으로만 — 전체 접힘은 설명에 기여하지 않는다

    cmd.show("spheres", "metals")
    cmd.color("0xE69F00", "metals and elem ZN")
    cmd.color("0xCC79A7", "metals and elem MG")

    cmd.show("sticks", "via")
    cmd.color(C_REF, "via and elem C")
    cmd.util.cnc("via and not elem C")

    x, up, z, n_resi, lig = camera_axes()

    # 표면은 쓰지 않는다 — 어떤 투명도에서도 리간드를 가려 결합 구조 설명을 방해한다.
    # 포켓은 리간드를 둘러싼 잔기 17개의 배치로 보여준다.

    if contacts:
        # 리간드–포켓 극성 접촉 (PyMOL 기본 기하 기준)
        # 그림에 그리는 점선과 캡션에 적는 거리가 같은 집합이어야 한다.
        # PyMOL mode=2 는 자체 기준으로 일부만 그려서 둘이 어긋날 수 있으므로,
        # 아래에서 좌표로 잰 접촉을 그대로 하나씩 그린다.
        contacts_found = polar_contacts()
        for k, (dist, a, b) in enumerate(contacts_found):
            la = a.split()[1]
            rn, ri = b.split()[0][:3], b.split()[0][3:]
            rb = b.split()[1]
            cmd.distance(f"polar{k}",
                         f"via and name {la}",
                         f"pocket and resn {rn} and resi {ri} and name {rb}")
            cmd.hide("labels", f"polar{k}")
            cmd.show("dashes", f"polar{k}")
            cmd.enable(f"polar{k}")
        # [주의] cmd.distance() 의 반환값은 생성된 측정들의 **평균 거리**다.
        # 그것을 단일 수소결합 거리로 적으면 존재하지 않는 값을 보고하게 된다
        # (실제로 한 번 그렇게 적었다). 그래서 접촉을 하나씩 좌표에서 다시 잰다.
        for dist, a, b in contacts_found:
            print(f"  극성 접촉 {dist:.3f} Å — {a} ··· {b}")
        import json as _j
        with open(os.path.join(OUT, "polar_contacts.json"), "w",
                  encoding="utf-8") as _f:
            _j.dump([{"distance_angstrom": round(d, 3), "ligand_atom": a,
                      "residue_atom": b} for d, a, b in contacts_found],
                    _f, ensure_ascii=False, indent=1)

    if label:
        for resn_resi, text in LABELLED.items():
            resn, resi = resn_resi[:3], resn_resi[3:]
            sel = f"pocket and resn {resn} and resi {resi} and name CB"
            if cmd.count_atoms(sel) == 0:
                sel = f"pocket and resn {resn} and resi {resi} and name CA"
            if cmd.count_atoms(sel):
                pass  # 라벨은 합성 단계에서 그린다

    v = list(cmd.get_view())
    v[0:9] = [float(t) for t in np.vstack([x, up, z]).flatten()]
    cmd.set_view(v)
    cmd.zoom("pocket", 1.25)
    cmd.center("via")
    cmd.clip("slab", SLAB)
    print(f"  방향: N-말단 {n_resi}번 위 · 포켓 정면 · slab {SLAB}")
    return cmd.get_view()



def polar_contacts(cutoff=3.3):
    """리간드와 포켓 잔기 사이의 극성 접촉을 좌표에서 하나씩 잰다.

    평균이 아니라 개별 거리를 돌려준다. 반환: [(거리, 리간드원자, 잔기원자), ...]
    """
    lig, poc = [], []
    cmd.iterate_state(1, "via and elem N+O",
                      "lig.append((name, x, y, z))", space={"lig": lig})
    cmd.iterate_state(1, "pocket and elem N+O",
                      "poc.append((resn, resi, name, x, y, z))",
                      space={"poc": poc})
    out = []
    for ln, lx, ly, lz in lig:
        for rn, ri, pn, px, py, pz in poc:
            d = float(np.linalg.norm(np.array([lx - px, ly - py, lz - pz])))
            if d <= cutoff:
                out.append((d, f"VIA {ln}", f"{rn}{ri} {pn}"))
    return sorted(out)


def load_pose(path, name, state, color):
    cmd.load(path, "_tmp", discrete=1)
    cmd.create(name, "_tmp", state, 1)
    cmd.delete("_tmp")
    cmd.remove(f"{name} and hydro")
    cmd.show("sticks", name)
    cmd.set("stick_radius", 0.165, name)
    cmd.color(color, f"{name} and elem C")
    cmd.util.cnc(f"{name} and not elem C")


def png(path):
    cmd.ray(W, H)
    cmd.png(path, dpi=300)
    print(f"  saved {os.path.basename(path)}")



def label_anchors():
    """잔기마다 고유 색 구를 한 번 더 렌더링해, 그 색의 픽셀 중심으로 위치를 얻는다.

    투영 수식을 쓰지 않으므로 카메라 설정과 무관하게 항상 맞는다.
    """
    targets = [
        ("Gln817", "pocket and resn GLN and resi 817 and sidechain"),
        ("Phe820", "pocket and resn PHE and resi 820 and sidechain"),
        ("Val782", "pocket and resn VAL and resi 782 and sidechain"),
        ("Tyr612", "pocket and resn TYR and resi 612 and sidechain"),
        ("Zn", "metals and elem ZN"),
        ("Mg", "metals and elem MG"),
    ]
    palette = [(255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 0, 255),
               (0, 255, 255), (255, 255, 0)]
    assert len(palette) >= 6

    keep = cmd.get_view()
    cmd.hide("everything")
    # 마커는 색을 그대로 읽어야 하므로 조명·반사를 끈다
    _light = (cmd.get("ambient"), cmd.get("specular"), cmd.get("reflect"))
    cmd.set("ambient", 1.0)
    cmd.set("specular", 0.0)
    cmd.set("reflect", 0.0)
    marks = []
    for i, (text, sel) in enumerate(targets):
        if cmd.count_atoms(sel) == 0:
            continue
        name = f"mk{i}"
        cmd.create(name, sel)
        cmd.show("spheres", name)
        cmd.set("sphere_scale", 0.9, name)
        r, g, b = palette[i]
        cmd.set_color(f"mkc{i}", [r / 255, g / 255, b / 255])
        cmd.color(f"mkc{i}", name)
        marks.append((text, (r, g, b), name))
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.set_view(keep)
    mk_png = f"{OUT}/_marker.png"
    cmd.ray(W, H)
    cmd.png(mk_png, dpi=300)
    cmd.refresh()

    import time
    for _ in range(100):
        if os.path.exists(mk_png) and os.path.getsize(mk_png) > 0:
            break
        time.sleep(0.1)

    from PIL import Image
    import numpy as _np
    a = _np.asarray(Image.open(mk_png).convert("RGB")).astype(int)
    hh, ww, _ = a.shape
    out = {}
    for text, (r, g, b), name in marks:
        d = _np.abs(a - _np.array([r, g, b])).sum(axis=2)
        m = d < 120
        if m.sum() >= 20:
            ys, xs = _np.nonzero(m)
            out[text] = [round(float(xs.mean() / ww), 4),
                         round(float(ys.mean() / hh), 4)]
        cmd.delete(name)
    cmd.set("ambient", _light[0])
    cmd.set("specular", _light[1])
    cmd.set("reflect", _light[2])
    cmd.set_view(keep)
    return out



# ------------------------------------------------------------------ A
view = base_scene(label=True, contacts=True)
res = []
cmd.iterate("pocket and name CA", "res.append((int(resi), resn))",
            space={"res": res})
print("포켓 잔기 (VIA 4.5 Å 이내):",
      ", ".join(f"{n}{i}" for i, n in sorted(set(res))))
png(f"{OUT}/panelA_1udt.png")
# 앵커 계산은 장면을 숨기므로 패널 렌더 뒤에 한다
anchors = label_anchors()
import json as _json
with open(f"{OUT}/anchors.json", "w", encoding="utf-8") as _f:
    _json.dump(anchors, _f, ensure_ascii=False, indent=1)
print("  앵커:", anchors)

# ------------------------------------------------------------------ B
base_scene(surface=False)
cmd.set_view(view)
cp = f"{SR}/structures/work_controlled/control_poses.sdf"
load_pose(cp, "redock_top", 1, C_TOP)
load_pose(cp, "redock_best", 5, C_BEST)
png(f"{OUT}/panelB_redock.png")

# ------------------------------------------------------------------ C
# 해석 가능 대역(mid·near)의 점수 1위. far 대역은 공유 부분구조가 작아 RMSD 가
# 자세 정확도를 재지 못하므로(§3.6) 이 패널에 쓰지 않는다.
base_scene(surface=False)
cmd.set_view(view)
_c = f"{SR}/structures/work_controlled/CHEMBL3944571_poses.sdf"
load_pose(_c, "top_scorer", 1, C_TOP)    # 점수 1위 자세
load_pose(_c, "top_best", 5, C_BEST)     # 같은 화합물에서 RMSD 가 가장 작은 모드
png(f"{OUT}/panelC_top_scorer.png")

# ------------------------------------------------------------------ D
base_scene(surface=False)
cmd.set_view(view)
load_pose(f"{SR}/structures/work_controlled/CHEMBL47108_poses.sdf",
          "closest", 1, C_BEST)
png(f"{OUT}/panelD_closest.png")

print("PYMOL-DONE")
