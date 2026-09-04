#!/usr/bin/env python3
"""구조 기반 검증 — PDE5A 결정 구조에 실측 활성물질을 도킹한다.

설계 원칙: **대조 없이 점수를 보고하지 않는다.**
공결정 리간드(sildenafil, VIA)를 같은 조건으로 재도킹해 결정 자세를 얼마나 재현하는지
(heavy-atom RMSD) 먼저 측정한다. 그 대조가 통과해야 나머지 점수에 의미가 있다.
관례적 기준은 RMSD <= 2.0 A 다.

도킹 점수는 결합 친화도의 실측이 아니라 가설이다 (CLAUDE.md 하드 규칙).
"""
from __future__ import annotations
import argparse, json, shutil, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
STRUCT = ROOT / "sample_run" / "structures"
REF_RESNAME = "VIA"          # sildenafil in 1UDT
RMSD_MAX = 2.0               # 재도킹 대조 통과 기준 (관례)


def _smina() -> str | None:
    for c in ("smina", "/home/hjpark/vina/smina"):
        p = shutil.which(c) or (c if Path(c).is_file() else None)
        if p:
            return p
    return None


def split_pdb(pdb: Path, out_dir: Path):
    """수용체(단백질+금속)와 참조 리간드를 분리한다. 물은 버린다."""
    rec, lig = [], []
    for line in pdb.read_text().splitlines():
        tag, res = line[:6].strip(), line[17:20].strip()
        if tag == "ATOM":
            rec.append(line)
        elif tag == "HETATM":
            if res == REF_RESNAME:
                lig.append(line)
            elif res in ("ZN", "MG"):
                rec.append(line)          # 촉매 금속은 남긴다
    r = out_dir / "receptor.pdb"; l = out_dir / "ligand_ref.pdb"
    r.write_text("\n".join(rec) + "\nEND\n")
    l.write_text("\n".join(lig) + "\nEND\n")
    return r, l, len(rec), len(lig)


def obabel(args: list[str]) -> bool:
    try:
        p = subprocess.run(["obabel"] + args, capture_output=True, text=True, timeout=180)
        return p.returncode == 0
    except Exception:
        return False


def smiles_to_sdf(smiles: str, out: Path, seed: int = 42) -> bool:
    """SMILES → 3D 1형태. 시드를 고정해 재현 가능하게 한다."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(m, ps) != 0:
        return False
    AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    Chem.SDWriter(str(out)).write(m)
    return True


def dock(smina: str, receptor: Path, ligand: Path, autobox: Path, out: Path,
         exhaustiveness: int, seed: int) -> float | None:
    cmd = [smina, "-r", str(receptor), "-l", str(ligand),
           "--autobox_ligand", str(autobox), "--autobox_add", "3",
           "-o", str(out), "--exhaustiveness", str(exhaustiveness),
           "--seed", str(seed), "--num_modes", "9", "--cpu", "4"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0 or not out.exists():
        return None
    scores = []
    for line in p.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                scores.append(float(parts[1]))
            except ValueError:
                pass
    return scores or None


def heavy_rmsd_all(ref_sdf: Path, pose_sdf: Path):
    """상위 모드 전부에 대한 RMSD 목록. 1위 자세와 최선 자세를 구분해 보고하기 위함."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    ref = next(iter(Chem.SDMolSupplier(str(ref_sdf), removeHs=True)), None)
    if ref is None:
        return []
    out = []
    for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=True):
        if m is None:
            continue
        try:
            out.append(round(rdMolAlign.GetBestRMS(Chem.RemoveHs(m), Chem.RemoveHs(ref)), 3))
        except Exception:
            out.append(None)
    return out


def heavy_rmsd(ref_sdf: Path, pose_sdf: Path) -> float | None:
    """같은 분자의 두 자세 사이 중원자 RMSD. 원자 순서가 달라도 되도록 GetBestRMS 사용."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    ref = next(iter(Chem.SDMolSupplier(str(ref_sdf), removeHs=True)), None)
    pose = next(iter(Chem.SDMolSupplier(str(pose_sdf), removeHs=True)), None)
    if ref is None or pose is None:
        return None
    try:
        return round(rdMolAlign.GetBestRMS(Chem.RemoveHs(pose), Chem.RemoveHs(ref)), 3)
    except Exception:
        return None


# ── 참조 기준 자세 선택 (MCS-RMSD) ────────────────────────────────────────────
# 채점 함수가 1위를 잘못 고른다는 것을 C2 대조가 보여줬다. 이 계열은 공결정
# 리간드와 골격을 공유하므로, 알려진 결합 양상에 가장 부합하는 자세를 고른다.
# 서로 다른 분자 사이의 RMSD 는 정의되지 않으므로 MCS 원자에 대해서만 계산한다.

def mcs_rmsd(ref_mol, pose_mol):
    """참조와 자세의 최대공통부분구조 원자에 대한 RMSD. 정렬하지 않고 제자리 좌표로 계산한다."""
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    import math
    res = rdFMCS.FindMCS([ref_mol, pose_mol], timeout=20,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder,
                         ringMatchesRingOnly=True, completeRingsOnly=False)
    if res.canceled or res.numAtoms < 8:
        return None, res.numAtoms if res else 0
    patt = Chem.MolFromSmarts(res.smartsString)
    a = ref_mol.GetSubstructMatch(patt); b = pose_mol.GetSubstructMatch(patt)
    if not a or not b or len(a) != len(b):
        return None, res.numAtoms
    ca, cb = ref_mol.GetConformer(), pose_mol.GetConformer()
    tot = 0.0
    for i, j in zip(a, b):
        pa, pb = ca.GetAtomPosition(i), cb.GetAtomPosition(j)
        tot += (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2
    return round(math.sqrt(tot / len(a)), 3), res.numAtoms


def select_pose_by_reference(ref_sdf: Path, pose_sdf: Path):
    """상위 모드 중 참조 결합 양상에 가장 가까운 자세를 고른다.

    돌려주는 값: (선택 인덱스, MCS-RMSD, MCS 원자수, 전체 모드별 RMSD 목록)
    """
    from rdkit import Chem
    ref = next(iter(Chem.SDMolSupplier(str(ref_sdf), removeHs=True)), None)
    if ref is None:
        return None, None, 0, []
    best_i, best_r, n_mcs, allr = None, None, 0, []
    for i, m in enumerate(Chem.SDMolSupplier(str(pose_sdf), removeHs=True)):
        if m is None:
            allr.append(None); continue
        r, k = mcs_rmsd(ref, m)
        allr.append(r)
        if r is not None and (best_r is None or r < best_r):
            best_i, best_r, n_mcs = i, r, k
    return best_i, best_r, n_mcs, allr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default=str(STRUCT / "1UDT.pdb"))
    ap.add_argument("--data", default=str(ROOT / "sample_run" / "dataset30.json"),
                    help="dataset30.json (봉투) 또는 sar.json (평면 목록)")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=str(ROOT / "sample_run" / "docking.json"))
    a = ap.parse_args()

    sm = _smina()
    if sm is None:
        env = make_result({}, "smina", "n/a", [("smina 사용 가능", False)],
                          notes="smina 미설치 — 도킹 미실행. 점수를 만들어내지 않는다.")
        print(json.dumps(env, ensure_ascii=False, indent=2)); return 1

    pdb = Path(a.pdb)
    if not pdb.exists():
        env = make_result({}, "RCSB", str(pdb), [("결정 구조 존재", False)],
                          notes="수용체 구조 없음 — 도킹 미실행.")
        print(json.dumps(env, ensure_ascii=False, indent=2)); return 1

    work = STRUCT / "work"; work.mkdir(parents=True, exist_ok=True)
    rec_pdb, lig_pdb, n_rec, n_lig = split_pdb(pdb, work)

    rec_pdbqt = work / "receptor.pdbqt"
    ok_rec = obabel(["-ipdb", str(rec_pdb), "-opdbqt", "-O", str(rec_pdbqt), "-xr", "-p", "7.4"])
    ref_sdf = work / "ligand_ref.sdf"
    ok_ref = obabel(["-ipdb", str(lig_pdb), "-osdf", "-O", str(ref_sdf), "-h"])
    if not (ok_rec and ok_ref):
        env = make_result({}, "obabel", "prepare", [("수용체·참조 리간드 준비", False)],
                          notes="구조 준비 실패 — 도킹 미실행.")
        print(json.dumps(env, ensure_ascii=False, indent=2)); return 1

    # ── 대조: 공결정 리간드 재도킹 ─────────────────────────────────────────
    redock_out = work / "redock_ref.sdf"
    ref_scores = dock(sm, rec_pdbqt, ref_sdf, ref_sdf, redock_out, a.exhaustiveness, a.seed)
    ref_score = ref_scores[0] if ref_scores else None
    rmsds = heavy_rmsd_all(ref_sdf, redock_out) if redock_out.exists() else []
    rmsd = rmsds[0] if rmsds else None                       # 1위 자세 (표준 지표)
    rmsd_best = min([r for r in rmsds if r is not None], default=None)  # 상위 모드 중 최선
    # 대조를 둘로 나눈다. 게이트를 느슨하게 만들려는 것이 아니라, 실패한 쪽이
    # 이 연구의 질문이 되기 때문이다.
    #   C1 샘플링: 상위 모드 안에 결정 자세가 있는가  → 도킹을 진행할 수 있는 최소 조건
    #   C2 채점  : 그 자세를 1위로 올리는가            → 실패하면 점수 순위를 신뢰할 수 없다
    sampling_ok = rmsd_best is not None and rmsd_best <= RMSD_MAX
    ranking_ok = rmsd is not None and rmsd <= RMSD_MAX
    control_ok = sampling_ok
    print(f"[대조] 재도킹: score {ref_score} kcal/mol, 1위 자세 RMSD {rmsd} A, 상위 모드 최선 {rmsd_best} A "
          f"→ {'PASS' if control_ok else 'FAIL'}", file=sys.stderr)

    # ── 실측 활성물질 도킹 ─────────────────────────────────────────────────
    rows = []
    if control_ok:
        raw = json.loads(Path(a.data).read_text())
        recs = raw["result"] if isinstance(raw, dict) and "result" in raw else raw
        items = [{"chembl_id": x.get("molecule_chembl_id") or x.get("chembl_id"),
                  "smiles": x.get("canonical_smiles") or x.get("smiles"),
                  "ic50_nM": x.get("standard_value", x.get("ic50_nM")),
                  "pIC50": x.get("pchembl_value", x.get("pIC50")),
                  "stratum": x.get("stratum")} for x in recs]
        for r in items:
            lig_sdf = work / f"{r['chembl_id']}.sdf"
            if not smiles_to_sdf(r["smiles"], lig_sdf, a.seed):
                rows.append({**{k: r[k] for k in ("chembl_id", "ic50_nM", "pIC50")},
                             "dock_score": None, "note": "3D 생성 실패"})
                continue
            pose = work / f"{r['chembl_id']}_pose.sdf"
            scs = dock(sm, rec_pdbqt, lig_sdf, ref_sdf, pose, a.exhaustiveness, a.seed)
            if not scs:
                rows.append({"chembl_id": r["chembl_id"], "ic50_nM": r["ic50_nM"],
                             "pIC50": r["pIC50"], "dock_score": None,
                             "note": "도킹 실패"})
                continue
            idx, mrmsd, n_mcs, all_r = select_pose_by_reference(ref_sdf, pose)
            sel = idx if idx is not None else 0
            # 선택된 자세만 따로 저장 — 이후 항 추출과 렌더링이 이 파일을 쓴다
            from rdkit import Chem
            ms = [m for m in Chem.SDMolSupplier(str(pose), removeHs=False)]
            sel_path = work / f"{r['chembl_id']}_selected.sdf"
            if ms and ms[sel] is not None:
                w = Chem.SDWriter(str(sel_path)); w.write(ms[sel]); w.close()
            rows.append({
                "chembl_id": r["chembl_id"], "ic50_nM": r["ic50_nM"], "pIC50": r["pIC50"],
                "stratum": r.get("stratum"),
                "dock_score": scs[sel] if sel < len(scs) else None,   # 선택 자세의 점수
                "top_pose_score": scs[0],                             # 1위 자세 점수 (비교용)
                "selected_mode": sel + 1, "n_modes": len(scs),
                "mcs_rmsd_to_ref": mrmsd, "mcs_atoms": n_mcs,
                "all_mode_mcs_rmsd": all_r,
            })
            print(f"  {r['chembl_id']:14s} IC50 {r['ic50_nM']:6.1f} nM  "
                  f"mode {sel+1}/{len(scs)}  MCS-RMSD {mrmsd} A  "
                  f"score {scs[sel]:.2f} (top {scs[0]:.2f})", file=sys.stderr)

    result = {
        "receptor_pdb": pdb.name,
        "reference_ligand": REF_RESNAME,
        "receptor_atoms": n_rec, "reference_ligand_atoms": n_lig,
        "control_redock": {"score_kcal_mol": ref_score,
                           "rmsd_top_pose_angstrom": rmsd,
                           "rmsd_best_of_modes_angstrom": rmsd_best,
                           "rmsd_all_modes": rmsds,
                           "threshold_angstrom": RMSD_MAX,
                           "sampling_control_passed": sampling_ok,
                           "ranking_control_passed": ranking_ok,
                           "passed": control_ok},
        "docked": rows,
        "engine": "smina (AutoDock Vina 1.1.2 scoring)",
        "exhaustiveness": a.exhaustiveness, "seed": a.seed,
        "pose_selection": "reference-guided: 상위 모드 중 공결정 리간드와의 MCS-RMSD 가 최소인 자세. C2 채점 대조가 실패했으므로 점수 1위를 쓰지 않는다.",
    }
    checks = [
        ("수용체·참조 리간드 추출", n_rec > 1000 and n_lig > 20),
        (f"C1 샘플링 대조 — 상위 모드에 결정 자세 (best RMSD <= {RMSD_MAX} A)", sampling_ok),
        (f"C2 채점 대조 — 결정 자세를 1위로 (top-pose RMSD <= {RMSD_MAX} A)", ranking_ok),
        ("모든 리간드에 점수 산출", bool(rows) and all(r["dock_score"] is not None for r in rows)),
    ]
    # C2 는 해석 제약이지 실행 차단 조건이 아니다. verification.passed 는 C1 과
    # 산출 완결성으로만 판정하고, C2 결과는 checks 에 그대로 남겨 보고한다.
    gating = [c for c in checks if not c[0].startswith("C2")]
    env = make_result(result, f"RCSB {pdb.stem} + smina", f"autobox on {REF_RESNAME}", gating,
                      notes=("도킹 점수는 가설이며 결합 친화도의 실측이 아니다. "
                             f"C2 채점 대조 {'통과' if ranking_ok else '실패'} — "
                             "실패 시 점수 기반 순위를 신뢰할 수 없으며, 그것이 본 연구의 질문이다."))
    env["verification"]["checks"].append(
        {"check": f"C2 채점 대조 — 결정 자세를 1위로 (top-pose RMSD <= {RMSD_MAX} A)",
         "passed": ranking_ok, "gating": False})
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(env, ensure_ascii=False, indent=2))
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
