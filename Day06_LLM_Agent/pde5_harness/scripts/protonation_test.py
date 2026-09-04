#!/usr/bin/env python3
"""프로토네이션 비대칭이 C2 실패의 원인인가 — 2×2 통제 실험.

현행 파이프라인은 비대칭이다. 수용체는 `obabel -xr -p 7.4` 로 pH 7.4 상태를 만들지만
리간드는 `Chem.AddHs` 중성 그대로다. 실데나필의 피페라진 질소는 생리적 pH 에서
양성자화되므로(obabel 확인: 형식전하 0 → +1), **중성 리간드를 pH 7.4 수용체에 넣은 것이
C2 실패의 직접 원인일 수 있다.** 이전 판은 이를 한계로만 적고 검정하지 않았다.

설계
    수용체 {pH 7.4, 중성} × 리간드 {중성, pH 7.4} = 4 조건 × 시드 3개
  - (pH7.4, pH7.4)  : 물리적으로 옳은 짝. C2 가 여기서 좋아지면 비대칭이 원인이다.
  - (pH7.4, 중성)   : 현행. 기준선.
  - (중성, 중성)    : 짝은 맞지만 단백질로서는 비현실적. 비대칭 효과를 격리하는 대조.
  - (중성, pH7.4)   : 반대 방향 비대칭. 대칭성 확인용.

RMSD 는 중원자 기준(removeHs=True)이라 프로토네이션이 지표 자체를 바꾸지 않는다.
따라서 조건 간 비교가 깨끗하다.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result   # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"
OUT = SR / "structures" / "protonation"
SMINA = os.environ.get("SMINA", "/home/hjpark/vina/smina")
RMSD_MAX = 2.0


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, timeout=1800, **kw)


def prepare():
    """수용체 2종·리간드 2종을 만든다. 이미 있으면 다시 만들지 않는다."""
    OUT.mkdir(parents=True, exist_ok=True)
    made = {}
    # 수용체 — pH 7.4 (현행) 는 기존 파일 재사용, 중성은 새로 만든다
    made["rec_pH74"] = WORK / "receptor.pdbqt"
    rn = OUT / "receptor_neutral.pdbqt"
    if not rn.exists():
        sh([os.environ.get("OBABEL", "obabel"), "-ipdb", str(WORK / "receptor.pdb"),
            "-opdbqt", "-O", str(rn), "-xr", "-h"])
    made["rec_neutral"] = rn
    # 리간드 — 중성(현행) 재사용, pH 7.4 는 obabel 로
    made["lig_neutral"] = WORK / "ligand_ref.sdf"
    lp = OUT / "ligand_ref_pH74.sdf"
    if not lp.exists():
        sh([os.environ.get("OBABEL", "obabel"), str(WORK / "ligand_ref.sdf"),
            "-osdf", "-p", "7.4", "-O", str(lp)])
    made["lig_pH74"] = lp
    return made


def charges(path: Path):
    """리간드의 총 형식전하 — 조건이 실제로 달라졌음을 산출물에 남긴다."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    return None if m is None else Chem.GetFormalCharge(m)


def mcs_rmsd(ref_mol, pose_mol):
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    res = rdFMCS.FindMCS([ref_mol, pose_mol], timeout=20,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder,
                         ringMatchesRingOnly=True, completeRingsOnly=False)
    if res.canceled or res.numAtoms < 8:
        return None
    patt = Chem.MolFromSmarts(res.smartsString)
    a = ref_mol.GetSubstructMatch(patt); b = pose_mol.GetSubstructMatch(patt)
    if not a or not b or len(a) != len(b):
        return None
    ca, cb = ref_mol.GetConformer(), pose_mol.GetConformer()
    s = sum((ca.GetAtomPosition(i).x - cb.GetAtomPosition(j).x) ** 2 +
            (ca.GetAtomPosition(i).y - cb.GetAtomPosition(j).y) ** 2 +
            (ca.GetAtomPosition(i).z - cb.GetAtomPosition(j).z) ** 2
            for i, j in zip(a, b))
    return round((s / len(a)) ** 0.5, 3)


def run(job):
    tag, rec, lig, seed, exh, nmodes = job
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    o = OUT / f"{tag}_s{seed}.sdf"
    cmd = [SMINA, "-r", str(rec), "-l", str(lig),
           "--autobox_ligand", str(WORK / "ligand_ref.sdf"), "--autobox_add", "3",
           "-o", str(o), "--exhaustiveness", str(exh), "--seed", str(seed),
           "--num_modes", str(nmodes), "--cpu", "4"]
    p = sh(cmd)
    if p.returncode != 0 or not o.exists():
        return {"condition": tag, "seed": seed, "error": f"rc={p.returncode}"}
    scores = [float(x.split()[1]) for x in p.stdout.splitlines()
              if len(x.split()) >= 2 and x.split()[0].isdigit()]
    # 기준 자세는 언제나 결정 구조의 공결정 리간드다 (조건과 무관하게 고정)
    ref = Chem.SDMolSupplier(str(WORK / "ligand_ref.sdf"), removeHs=True)[0]
    poses = [m for m in Chem.SDMolSupplier(str(o), removeHs=True) if m]
    rs = [mcs_rmsd(ref, m) for m in poses]
    ok = [r for r in rs if r is not None]
    best = min(ok) if ok else None
    return {"condition": tag, "seed": seed, "n_modes": len(scores),
            "top_score": scores[0] if scores else None,
            "c1_best_rmsd": best,
            "c1_best_mode": (rs.index(best) + 1) if best is not None else None,
            "c2_top_rmsd": rs[0] if rs else None,
            "c1_pass": best is not None and best <= RMSD_MAX,
            "c2_pass": bool(rs and rs[0] is not None and rs[0] <= RMSD_MAX),
            "all_rmsd": rs}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=64)
    ap.add_argument("--seeds", default="42,7,2024")
    ap.add_argument("--num-modes", type=int, default=20)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out", default=str(SR / "protonation_test.json"))
    a = ap.parse_args()

    f = prepare()
    seeds = [int(x) for x in a.seeds.split(",")]
    conds = [("recpH74_ligpH74", f["rec_pH74"], f["lig_pH74"]),
             ("recpH74_ligneutral", f["rec_pH74"], f["lig_neutral"]),
             ("recneutral_ligneutral", f["rec_neutral"], f["lig_neutral"]),
             ("recneutral_ligpH74", f["rec_neutral"], f["lig_pH74"])]
    lig_charge = {"lig_neutral": charges(f["lig_neutral"]),
                  "lig_pH74": charges(f["lig_pH74"])}
    print(f"리간드 총 형식전하 — 중성 {lig_charge['lig_neutral']} / "
          f"pH7.4 {lig_charge['lig_pH74']}", file=sys.stderr)
    if lig_charge["lig_neutral"] == lig_charge["lig_pH74"]:
        print("경고: 두 리간드의 전하가 같다 — 조건이 실제로 다르지 않다", file=sys.stderr)

    jobs = [(t, r, l, s, a.exhaustiveness, a.num_modes)
            for t, r, l in conds for s in seeds]
    print(f"{len(conds)} 조건 × 시드 {seeds} = {len(jobs)}회", file=sys.stderr)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            rows.append(fu.result())
            if i % 3 == 0: print(f"  {i}/{len(jobs)}", file=sys.stderr)

    summary = []
    for t, _, _ in conds:
        g = [r for r in rows if r["condition"] == t and "error" not in r]
        if not g: continue
        c1 = [r["c1_best_rmsd"] for r in g if r["c1_best_rmsd"] is not None]
        c2 = [r["c2_top_rmsd"] for r in g if r["c2_top_rmsd"] is not None]
        ts = [r["top_score"] for r in g if r["top_score"] is not None]
        summary.append({
            "condition": t, "n_runs": len(g),
            "c1_best_rmsd_mean": round(sum(c1) / len(c1), 3) if c1 else None,
            "c2_top_rmsd_mean": round(sum(c2) / len(c2), 3) if c2 else None,
            "c2_top_rmsd_min": round(min(c2), 3) if c2 else None,
            "top_score_mean": round(sum(ts) / len(ts), 3) if ts else None,
            "c1_pass_rate": round(sum(r["c1_pass"] for r in g) / len(g), 2),
            "c2_pass_rate": round(sum(r["c2_pass"] for r in g) / len(g), 2)})

    by = {s["condition"]: s for s in summary}
    base = by.get("recpH74_ligneutral", {})       # 현행 (비대칭)
    matched = by.get("recpH74_ligpH74", {})        # 물리적으로 옳은 짝
    improved = (matched.get("c2_top_rmsd_mean") is not None
                and base.get("c2_top_rmsd_mean") is not None
                and matched["c2_top_rmsd_mean"] < base["c2_top_rmsd_mean"] - 0.5)
    fixed = matched.get("c2_pass_rate", 0) > 0

    if fixed:
        verdict = ("프로토네이션 비대칭이 원인이다: 리간드를 pH 7.4 로 맞추면 C2 대조가 통과한다")
    elif improved:
        verdict = ("프로토네이션이 일부 기여한다: C2 가 개선되나 기준을 넘지는 못한다")
    else:
        verdict = ("프로토네이션 비대칭은 C2 실패의 원인이 아니다: 리간드를 pH 7.4 로 맞춰도 "
                   "C2 가 개선되지 않는다")

    out = {"design": "수용체 {pH7.4, 중성} × 리간드 {중성, pH7.4}, 시드 3개",
           "seeds": seeds, "exhaustiveness": a.exhaustiveness, "num_modes": a.num_modes,
           "ligand_formal_charge": lig_charge,
           "protonation_site": "피페라진 N-메틸 질소 (obabel -p 7.4 에서 +1)",
           "rmsd_note": "RMSD 는 중원자 기준이라 프로토네이션이 지표를 바꾸지 않는다",
           "runs": rows, "summary": summary,
           "baseline_current_pipeline": "recpH74_ligneutral",
           "c2_change_matched_minus_baseline": (
               round(matched["c2_top_rmsd_mean"] - base["c2_top_rmsd_mean"], 3)
               if matched.get("c2_top_rmsd_mean") is not None
               and base.get("c2_top_rmsd_mean") is not None else None),
           "verdict": verdict}
    checks = [("4 조건 모두 실행", len(summary) == 4),
              ("두 리간드 전하가 실제로 다름",
               lig_charge["lig_neutral"] != lig_charge["lig_pH74"]),
              ("시드 복수", len(seeds) >= 2),
              ("전 실행 성공", not [r for r in rows if "error" in r])]
    env = make_result(out, f"smina 프로토네이션 2×2 ({len(jobs)}회)",
                      f"exhaustiveness {a.exhaustiveness}", checks, notes=verdict)
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")

    print(f"\n{'조건':<26}{'C1 최선':>9}{'C2 1위':>9}{'점수':>9}{'C1통과':>8}{'C2통과':>8}")
    for s in summary:
        print(f"{s['condition']:<26}{s['c1_best_rmsd_mean']:>9.2f}"
              f"{s['c2_top_rmsd_mean']:>9.2f}{s['top_score_mean']:>9.2f}"
              f"{s['c1_pass_rate']:>8.2f}{s['c2_pass_rate']:>8.2f}")
    print(f"\n판정: {verdict}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
