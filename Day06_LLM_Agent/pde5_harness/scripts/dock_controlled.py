#!/usr/bin/env python3
"""골격 통제 데이터셋 병렬 도킹 — 128 코어를 실제로 쓴다.

이전 실행(n=30, exhaustiveness 16, 순차)의 두 약점을 없앤다.
  - 표본이 작아 신뢰구간이 결론을 못 지탱했다  → n 을 대폭 늘린다
  - 탐색이 얕아 자세 실패와 채점 실패가 섞였다 → exhaustiveness 를 올린다

자세는 두 규칙 모두 기록한다.
  - top      : 점수 1위. 채점 함수만으로 결정되며 공결정 리간드를 참조하지 않는다.
  - reference: 공결정 리간드와의 MCS-RMSD 최소. 참조를 쓰므로 해석에 제약이 붙는다.
분석은 top 을 주 결과로 삼고 reference 를 민감도 분석으로 병기한다.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"
WORK2 = SR / "structures" / "work_controlled"
SMINA = os.environ.get("SMINA", "/home/hjpark/vina/smina")
RMSD_MAX = 2.0


def smiles_to_sdf(smiles: str, out: Path, seed: int = 42) -> bool:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    # 염 제거 — 가장 큰 프래그먼트만 남긴다 (이전 실행에서 암묵적으로 일어나던 것을 명시화)
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        m = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(m, ps) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    Chem.SDWriter(str(out)).write(m)
    return out.exists()


def mcs_rmsd(ref_mol, pose_mol):
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    res = rdFMCS.FindMCS([ref_mol, pose_mol], timeout=20,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder,
                         ringMatchesRingOnly=True, completeRingsOnly=False)
    if res.canceled or res.numAtoms < 8:
        return None, (res.numAtoms if res else 0)
    patt = Chem.MolFromSmarts(res.smartsString)
    a = ref_mol.GetSubstructMatch(patt); b = pose_mol.GetSubstructMatch(patt)
    if not a or not b or len(a) != len(b):
        return None, res.numAtoms
    ca, cb = ref_mol.GetConformer(), pose_mol.GetConformer()
    s = 0.0
    for i, j in zip(a, b):
        pa, pb = ca.GetAtomPosition(i), cb.GetAtomPosition(j)
        s += (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2
    return round((s / len(a)) ** 0.5, 3), res.numAtoms


def one(job):
    """한 화합물을 도킹한다. 프로세스 풀에서 돌아야 하므로 최상위 함수여야 한다."""
    cid, smiles, exh, seed, cpu = job
    from rdkit import Chem
    lig = WORK2 / f"{cid}.sdf"
    if not lig.exists() and not smiles_to_sdf(smiles, lig, seed):
        return {"chembl_id": cid, "error": "3D 생성 실패"}
    pose = WORK2 / f"{cid}_poses.sdf"
    cmd = [SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(lig),
           "--autobox_ligand", str(WORK / "ligand_ref.sdf"), "--autobox_add", "3",
           "-o", str(pose), "--exhaustiveness", str(exh), "--seed", str(seed),
           "--num_modes", "9", "--cpu", str(cpu)]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        return {"chembl_id": cid, "error": "도킹 시간 초과"}
    if p.returncode != 0 or not pose.exists():
        return {"chembl_id": cid, "error": f"smina 실패 rc={p.returncode}"}
    scores = []
    for line in p.stdout.splitlines():
        pr = line.split()
        if len(pr) >= 2 and pr[0].isdigit():
            try: scores.append(float(pr[1]))
            except ValueError: pass
    if not scores:
        return {"chembl_id": cid, "error": "점수 파싱 실패"}

    ref_mol = Chem.SDMolSupplier(str(WORK / "ligand_ref.sdf"), removeHs=True)[0]
    poses = [m for m in Chem.SDMolSupplier(str(pose), removeHs=True) if m]
    rmsds, atoms = [], 0
    for m in poses:
        r, na = mcs_rmsd(ref_mol, m)
        rmsds.append(r); atoms = max(atoms, na)
    valid = [(i, r) for i, r in enumerate(rmsds) if r is not None]
    ridx, rrmsd = (min(valid, key=lambda t: t[1]) if valid else (None, None))
    return {"chembl_id": cid, "n_modes": len(scores),
            "top_pose_score": scores[0], "top_pose_mcs_rmsd": rmsds[0] if rmsds else None,
            "ref_pose_index": (ridx + 1) if ridx is not None else None,
            "ref_pose_score": scores[ridx] if ridx is not None else None,
            "ref_pose_mcs_rmsd": rrmsd, "mcs_atoms": atoms,
            "all_mode_scores": scores, "all_mode_mcs_rmsd": rmsds}


def control_redock(exh, seed):
    """공결정 리간드 재도킹 — C1 샘플링과 C2 채점을 나눠 판정한다."""
    from rdkit import Chem
    out = WORK2 / "control_poses.sdf"
    cmd = [SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(WORK / "ligand_ref.sdf"),
           "--autobox_ligand", str(WORK / "ligand_ref.sdf"), "--autobox_add", "3",
           "-o", str(out), "--exhaustiveness", str(exh), "--seed", str(seed),
           "--num_modes", "9", "--cpu", "16"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    scores = [float(x.split()[1]) for x in p.stdout.splitlines()
              if len(x.split()) >= 2 and x.split()[0].isdigit()]
    ref = Chem.SDMolSupplier(str(WORK / "ligand_ref.sdf"), removeHs=True)[0]
    poses = [m for m in Chem.SDMolSupplier(str(out), removeHs=True) if m]
    rs = [mcs_rmsd(ref, m)[0] for m in poses]
    ok = [r for r in rs if r is not None]
    best = min(ok) if ok else None
    return {"score_kcal_mol": scores[0] if scores else None,
            "all_mode_rmsd": rs, "all_mode_scores": scores,   # 그림이 모드별 막대를 그린다
            "rmsd_top_pose_angstrom": rs[0] if rs else None,
            "rmsd_best_of_modes_angstrom": best,
            "best_mode_index": (rs.index(best) + 1) if best is not None else None,
            "n_modes": len(scores), "threshold_angstrom": RMSD_MAX,
            "sampling_control_passed": best is not None and best <= RMSD_MAX,
            "ranking_control_passed": (rs[0] is not None and rs[0] <= RMSD_MAX) if rs else False}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=str(SR / "dataset_controlled.json"))
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu-per-job", type=int, default=4)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--out", default=str(SR / "docking_controlled.json"))
    a = ap.parse_args()

    WORK2.mkdir(parents=True, exist_ok=True)
    ds = json.loads(Path(a.dataset).read_text())["result"]
    comps = ds["compounds"]
    meta = {c["molecule_chembl_id"]: c for c in comps}
    jobs = [(c["molecule_chembl_id"], c["canonical_smiles"],
             a.exhaustiveness, a.seed, a.cpu_per_job) for c in comps]

    print(f"대조 재도킹 (exhaustiveness {a.exhaustiveness})…", file=sys.stderr)
    ctrl = control_redock(a.exhaustiveness, a.seed)
    print(f"  C1 샘플링 {ctrl['rmsd_best_of_modes_angstrom']} Å "
          f"({'PASS' if ctrl['sampling_control_passed'] else 'FAIL'})  |  "
          f"C2 채점 {ctrl['rmsd_top_pose_angstrom']} Å "
          f"({'PASS' if ctrl['ranking_control_passed'] else 'FAIL'})", file=sys.stderr)

    rows, errs, done = [], [], 0
    print(f"도킹 {len(jobs)}건, 워커 {a.workers} × {a.cpu_per_job} CPU…", file=sys.stderr)
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(one, j): j[0] for j in jobs}
        for f in as_completed(futs):
            r = f.result(); done += 1
            if "error" in r:
                errs.append(r)
            else:
                m = meta[r["chembl_id"]]
                r.update({k: m[k] for k in
                          ("pchembl_value", "tanimoto_to_sildenafil", "potency_bin",
                           "similarity_bin", "heavy_atoms", "n_measurements",
                           "standard_type", "canonical_smiles")})
                rows.append(r)
            if done % 20 == 0:
                print(f"  {done}/{len(jobs)}  (실패 {len(errs)})", file=sys.stderr)

    rows.sort(key=lambda r: -r["pchembl_value"])
    out = {"receptor_pdb": "1UDT.pdb", "reference_ligand": "VIA",
           "engine": f"smina (AutoDock Vina 1.1.2 scoring), exhaustiveness {a.exhaustiveness}",
           "exhaustiveness": a.exhaustiveness, "seed": a.seed,
           "n_attempted": len(jobs), "n_docked": len(rows), "errors": errs,
           "control_redock": ctrl, "rows": rows,
           "pose_selection": ("주 결과는 점수 1위 자세(top). 공결정 리간드 참조 없이 "
                              "채점 함수만으로 결정된다. 참조 기준 자세(reference)는 "
                              "민감도 분석으로만 쓴다.")}
    checks = [
        ("도킹 성공률 >= 95%", len(rows) >= 0.95 * len(jobs)),
        ("C1 샘플링 대조 통과", ctrl["sampling_control_passed"]),
        ("전 화합물 점수 산출", all(r.get("top_pose_score") is not None for r in rows)),
        ("9칸 전부 대표됨", len({(r["potency_bin"], r["similarity_bin"])
                                  for r in rows}) == 9),
    ]
    env = make_result(out, f"smina 병렬 도킹 ({a.workers} 워커)",
                      f"n={len(rows)}, exhaustiveness {a.exhaustiveness}", checks,
                      notes=(f"C1 {ctrl['rmsd_best_of_modes_angstrom']} Å / "
                             f"C2 {ctrl['rmsd_top_pose_angstrom']} Å. 실패 {len(errs)}건."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n도킹 완료 {len(rows)}/{len(jobs)}  게이트 "
          f"{'PASS' if env['verification']['passed'] else 'FAIL'}")
    for c in env["verification"]["checks"]:
        print(f"  {'OK ' if c['passed'] else 'NG '} {c['check']}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
