#!/usr/bin/env python3
"""탐색 깊이를 바꿔가며 두 대조를 따로 본다 — 자세 실패와 채점 실패를 갈라내는 결정적 실험.

관찰: exhaustiveness 16 에서 C2(채점) RMSD 가 2.98 Å 였는데 64 에서 8.37 Å 로 **나빠졌다.**
탐색이 부족해서 생긴 문제라면 반대여야 한다. 탐색을 늘리면 채점 함수가 결정 자세보다
좋게 평가하는 가짜 최소점을 더 많이 찾아내기 때문이라는 설명이 이 방향과 맞는다.

이 스크립트는 그 관찰을 통제 실험으로 만든다. 같은 리간드·수용체·상자에서
exhaustiveness 만 바꾸고, 시드도 여러 개 써서 우연을 배제한다.
  C1 = 상위 모드 중 결정 자세와 가장 가까운 것의 RMSD  (샘플링 능력)
  C2 = 점수 1위 모드의 RMSD                             (채점 능력)
C1 이 개선/유지되는데 C2 가 악화되면, 병목은 탐색이 아니라 채점이다.
"""
from __future__ import annotations
import argparse, json, subprocess, sys, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"
OUT = SR / "structures" / "sweep"
SMINA = os.environ.get("SMINA", "/home/hjpark/vina/smina")
RMSD_MAX = 2.0


def mcs_rmsd(ref_mol, pose_mol):
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    res = rdFMCS.FindMCS([ref_mol, pose_mol], timeout=20,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder,
                         ringMatchesRingOnly=True, completeRingsOnly=False)
    if res.canceled or res.numAtoms < 8:
        return None
    from rdkit import Chem as C
    patt = C.MolFromSmarts(res.smartsString)
    a = ref_mol.GetSubstructMatch(patt); b = pose_mol.GetSubstructMatch(patt)
    if not a or not b or len(a) != len(b):
        return None
    ca, cb = ref_mol.GetConformer(), pose_mol.GetConformer()
    s = 0.0
    for i, j in zip(a, b):
        pa, pb = ca.GetAtomPosition(i), cb.GetAtomPosition(j)
        s += (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2
    return round((s / len(a)) ** 0.5, 3)


def run(job):
    exh, seed, nmodes = job
    from rdkit import Chem
    OUT.mkdir(parents=True, exist_ok=True)
    o = OUT / f"ctrl_e{exh}_s{seed}.sdf"
    cmd = [SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(WORK / "ligand_ref.sdf"),
           "--autobox_ligand", str(WORK / "ligand_ref.sdf"), "--autobox_add", "3",
           "-o", str(o), "--exhaustiveness", str(exh), "--seed", str(seed),
           "--num_modes", str(nmodes), "--cpu", "4"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        return {"exhaustiveness": exh, "seed": seed, "error": "timeout"}
    scores = [float(x.split()[1]) for x in p.stdout.splitlines()
              if len(x.split()) >= 2 and x.split()[0].isdigit()]
    ref = Chem.SDMolSupplier(str(WORK / "ligand_ref.sdf"), removeHs=True)[0]
    poses = [m for m in Chem.SDMolSupplier(str(o), removeHs=True) if m]
    rs = [mcs_rmsd(ref, m) for m in poses]
    ok = [r for r in rs if r is not None]
    best = min(ok) if ok else None
    return {"exhaustiveness": exh, "seed": seed, "n_modes": len(scores),
            "top_score": scores[0] if scores else None,
            "best_score_of_modes": min(scores) if scores else None,
            "c1_best_rmsd": best,
            "c1_best_mode": (rs.index(best) + 1) if best is not None else None,
            "c2_top_rmsd": rs[0] if rs else None,
            "score_of_nearest_pose": (scores[rs.index(best)]
                                      if best is not None and rs.index(best) < len(scores) else None),
            "all_rmsd": rs, "all_scores": scores,
            "c1_pass": best is not None and best <= RMSD_MAX,
            "c2_pass": bool(rs and rs[0] is not None and rs[0] <= RMSD_MAX)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", default="4,8,16,32,64,128")
    ap.add_argument("--seeds", default="42,7,2024")
    ap.add_argument("--num-modes", type=int, default=20)
    ap.add_argument("--workers", type=int, default=18)
    ap.add_argument("--out", default=str(SR / "exhaustiveness_sweep.json"))
    a = ap.parse_args()

    levels = [int(x) for x in a.levels.split(",")]
    seeds = [int(x) for x in a.seeds.split(",")]
    jobs = [(e, s, a.num_modes) for e in levels for s in seeds]
    print(f"탐색 깊이 {levels} × 시드 {seeds} = {len(jobs)}회", file=sys.stderr)

    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % 3 == 0:
                print(f"  {i}/{len(jobs)}", file=sys.stderr)
    rows.sort(key=lambda r: (r["exhaustiveness"], r["seed"]))

    summary = []
    for e in levels:
        g = [r for r in rows if r["exhaustiveness"] == e and "error" not in r]
        if not g: continue
        c1 = [r["c1_best_rmsd"] for r in g if r["c1_best_rmsd"] is not None]
        c2 = [r["c2_top_rmsd"] for r in g if r["c2_top_rmsd"] is not None]
        ts = [r["top_score"] for r in g if r["top_score"] is not None]
        summary.append({
            "exhaustiveness": e, "n_runs": len(g),
            "c1_best_rmsd_mean": round(sum(c1) / len(c1), 3) if c1 else None,
            "c1_best_rmsd_min": round(min(c1), 3) if c1 else None,
            "c2_top_rmsd_mean": round(sum(c2) / len(c2), 3) if c2 else None,
            "c2_top_rmsd_max": round(max(c2), 3) if c2 else None,
            "top_score_mean": round(sum(ts) / len(ts), 3) if ts else None,
            "c1_pass_rate": round(sum(r["c1_pass"] for r in g) / len(g), 2),
            "c2_pass_rate": round(sum(r["c2_pass"] for r in g) / len(g), 2)})

    first, last = summary[0], summary[-1]
    c1 = [x["c1_best_rmsd_mean"] for x in summary if x["c1_best_rmsd_mean"] is not None]
    c2 = [x["c2_top_rmsd_mean"] for x in summary if x["c2_top_rmsd_mean"] is not None]
    sc = [x["top_score_mean"] for x in summary if x["top_score_mean"] is not None]
    c1_better = bool(c1) and (first["c1_best_rmsd_mean"] - last["c1_best_rmsd_mean"]) > 0.3
    c2_worse = bool(c2) and (last["c2_top_rmsd_mean"] - first["c2_top_rmsd_mean"]) > 0.5
    c2_flat = bool(c2) and (max(c2) - min(c2)) < 0.5
    score_better = bool(sc) and (first["top_score_mean"] - last["top_score_mean"]) > 0.2
    score_flat = bool(sc) and (max(sc) - min(sc)) < 0.2
    c1_fixed = (first["c1_pass_rate"] < 0.5) and (last["c1_pass_rate"] >= 0.99)
    c2_never = all(x["c2_pass_rate"] == 0 for x in summary)

    out = {"levels": levels, "seeds": seeds, "num_modes": a.num_modes,
           "runs": rows, "summary": summary,
           "diagnosis": {
               "c1_improved_with_depth": c1_better,
               "c1_fixed_by_depth": c1_fixed,
               "c2_degraded_with_depth": c2_worse,
               "c2_flat_across_depth": c2_flat,
               "c2_never_passes": c2_never,
               "score_improved_with_depth": score_better,
               "score_plateaued": score_flat,
               "c2_range_angstrom": round(max(c2) - min(c2), 3) if c2 else None,
               "score_range_kcal": round(max(sc) - min(sc), 3) if sc else None,
               "verdict": (
                   "탐색이 병목: 깊이를 올리면 C1·C2 가 함께 개선된다"
                   if (c1_better and not c2_flat and not c2_worse) else
                   "채점이 병목: 깊이가 샘플링(C1)을 고치지만 채점(C2)은 어느 깊이에서도 실패한다. "
                   "점수가 얕은 깊이에서 이미 포화한다면 탐색은 전역 최소를 이미 찾은 것이고, "
                   "그 최소가 결정 자세와 다른 위치에 있다는 뜻이다 — 탐색으로 고칠 수 없다."
                   if (c1_fixed and c2_never) else
                   "판정 보류: C1·C2 패턴이 두 전형에 맞지 않는다. 상자·리간드 준비를 재검토하라.")}}
    checks = [("모든 조건 실행됨", len([r for r in rows if "error" not in r]) == len(jobs)),
              ("최소 4개 깊이 비교", len(summary) >= 4),
              ("시드 복수 사용", len(seeds) >= 2)]
    env = make_result(out, f"smina 대조 재도킹 스윕 ({len(jobs)}회)",
                      f"levels={levels}, seeds={seeds}", checks,
                      notes=(f"C1 {first['c1_best_rmsd_mean']}→{last['c1_best_rmsd_mean']} Å, "
                             f"C2 {first['c2_top_rmsd_mean']}→{last['c2_top_rmsd_mean']} Å, "
                             f"점수 {first['top_score_mean']}→{last['top_score_mean']} kcal/mol."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n{'exh':>5}{'C1 최선':>10}{'C2 1위':>10}{'점수':>10}{'C1통과':>8}{'C2통과':>8}")
    for s in summary:
        print(f"{s['exhaustiveness']:>5}{s['c1_best_rmsd_mean']:>10.2f}"
              f"{s['c2_top_rmsd_mean']:>10.2f}{s['top_score_mean']:>10.2f}"
              f"{s['c1_pass_rate']:>8.2f}{s['c2_pass_rate']:>8.2f}")
    print(f"\n진단: {out['diagnosis']['verdict']}")
    print(f"  C1 개선 {c1_better} (통과율 {first['c1_pass_rate']}→{last['c1_pass_rate']}) · "
          f"C2 변동폭 {out['diagnosis']['c2_range_angstrom']} Å · "
          f"점수 변동폭 {out['diagnosis']['score_range_kcal']} kcal/mol")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
