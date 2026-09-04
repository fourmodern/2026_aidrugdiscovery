#!/usr/bin/env python3
"""골격 유사성과 역가를 설계로 분리한 데이터셋 — 이전 30건 설계의 교란을 없앤다.

**이전 설계의 결함.** 역가로만 층화해 뽑았더니 강한 화합물은 대부분 실데나필 유사체였고
약한 화합물은 무관한 골격이었다. 수용체가 실데나필 결합 구조이므로, 도킹 점수가 역가를
반영하는지 골격 유사성을 반영하는지 구분할 수 없었다. 관찰된 "선별 성능"이 전부 골격
편향일 가능성을 배제하지 못한다.

**이 설계.** 역가 3구간 × 공결정 리간드와의 Tanimoto 유사도 3구간 = 9칸 격자에서 같은
수씩 뽑는다. 그러면 각 역가 구간이 동일한 유사도 분포를 갖게 되어 두 축이 직교한다.
유사도 구간 안에서 역가 신호가 남는지를 보면 교란 없이 판정할 수 있다.

반복 측정은 화합물당 **중앙값**으로 집계한다 (이전에는 첫 레코드만 썼다).
"""
from __future__ import annotations
import argparse, json, statistics, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
TARGET = "CHEMBL1827"          # PDE5A
# 공결정 리간드 = 실데나필 (PDB 1UDT 의 VIA)
SILDENAFIL = ("CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1")
POT_BINS = [("weak", None, 6.0), ("medium", 6.0, 7.5), ("strong", 7.5, None)]
SIM_BINS = [("far", None, 0.25), ("mid", 0.25, 0.45), ("near", 0.45, None)]


def evenly(seq, k):
    """정렬된 목록에서 균등 간격으로 k 개. 결정적이라 시드가 필요 없다."""
    if len(seq) <= k:
        return list(seq)
    step = (len(seq) - 1) / (k - 1) if k > 1 else 0
    return [seq[int(round(i * step))] for i in range(k)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cell", type=int, default=20)
    ap.add_argument("--out", default=str(ROOT / "sample_run" / "dataset_controlled.json"))
    a = ap.parse_args()

    from chembl_webresource_client.new_client import new_client
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator, DataStructs
    RDLogger.DisableLog("rdApp.*")

    ref = Chem.MolFromSmiles(SILDENAFIL)
    if ref is None:
        raise SystemExit("참조 리간드 SMILES 파싱 실패")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ref_fp = gen.GetFingerprint(ref)

    print("ChEMBL 활성 데이터 수집 중 (전수)…", file=sys.stderr)
    act = new_client.activity.filter(
        target_chembl_id=TARGET, standard_type__in=["IC50", "Ki"],
        standard_relation="=", pchembl_value__isnull=False).only(
        ["molecule_chembl_id", "canonical_smiles", "standard_type",
         "standard_value", "standard_units", "pchembl_value"])

    # 화합물별로 모든 측정을 모은다 — 중앙값 집계를 위해
    by_mol = defaultdict(lambda: {"pv": [], "vals": [], "types": set(), "smi": None})
    seen_records = 0
    for r in act:
        cid, smi, pv = r.get("molecule_chembl_id"), r.get("canonical_smiles"), r.get("pchembl_value")
        if not (cid and smi and pv):
            continue
        seen_records += 1
        m = by_mol[cid]
        m["pv"].append(float(pv)); m["smi"] = m["smi"] or smi
        m["types"].add(r.get("standard_type"))
        if r.get("standard_value"):
            try: m["vals"].append(float(r["standard_value"]))
            except (TypeError, ValueError): pass
        if seen_records % 2000 == 0:
            print(f"  … {seen_records} 레코드 / {len(by_mol)} 화합물", file=sys.stderr)
    print(f"  수집 완료: {seen_records} 레코드 / {len(by_mol)} 화합물", file=sys.stderr)

    pool = []
    for cid, m in by_mol.items():
        mol = Chem.MolFromSmiles(m["smi"])
        if mol is None or mol.GetNumHeavyAtoms() < 12 or mol.GetNumHeavyAtoms() > 70:
            continue
        sim = DataStructs.TanimotoSimilarity(ref_fp, gen.GetFingerprint(mol))
        pool.append({"molecule_chembl_id": cid, "canonical_smiles": m["smi"],
                     "pchembl_value": round(statistics.median(m["pv"]), 3),
                     "n_measurements": len(m["pv"]),
                     "pchembl_spread": round(max(m["pv"]) - min(m["pv"]), 2),
                     "standard_type": "/".join(sorted(m["types"])),
                     "standard_value": round(statistics.median(m["vals"]), 4) if m["vals"] else None,
                     "tanimoto_to_sildenafil": round(sim, 3),
                     "heavy_atoms": mol.GetNumHeavyAtoms()})
    print(f"  유효 화합물 {len(pool)}건", file=sys.stderr)

    def pbin(v):
        for nm, lo, hi in POT_BINS:
            if (lo is None or v >= lo) and (hi is None or v < hi):
                return nm
    def sbin(v):
        for nm, lo, hi in SIM_BINS:
            if (lo is None or v >= lo) and (hi is None or v < hi):
                return nm

    grid = defaultdict(list)
    for r in pool:
        r["potency_bin"] = pbin(r["pchembl_value"])
        r["similarity_bin"] = sbin(r["tanimoto_to_sildenafil"])
        grid[(r["potency_bin"], r["similarity_bin"])].append(r)

    picked, cells = [], []
    for pnm, _, _ in POT_BINS:
        for snm, _, _ in SIM_BINS:
            avail = sorted(grid[(pnm, snm)], key=lambda x: x["pchembl_value"])
            take = evenly(avail, a.per_cell)
            picked += take
            cells.append({"potency": pnm, "similarity": snm,
                          "available": len(avail), "taken": len(take)})
            print(f"  {pnm:7s} × {snm:5s}  후보 {len(avail):5d} → {len(take):3d}", file=sys.stderr)

    pv = [r["pchembl_value"] for r in picked]
    sims = [r["tanimoto_to_sildenafil"] for r in picked]
    # 교란이 실제로 제거되었는지 — 역가와 유사도의 상관이 낮아야 한다
    n = len(picked); mp = sum(pv) / n; ms = sum(sims) / n
    num = sum((a_ - mp) * (b - ms) for a_, b in zip(pv, sims))
    den = (sum((a_ - mp) ** 2 for a_ in pv) * sum((b - ms) ** 2 for b in sims)) ** 0.5
    confound_r = round(num / den, 3) if den else 0.0

    out = {"compounds": picked, "cells": cells, "n": n,
           "confound_pearson_r_potency_vs_similarity": confound_r,
           "reference_ligand_smiles": SILDENAFIL,
           "pool_size": len(pool), "records_scanned": seen_records,
           "aggregation": "화합물당 pChEMBL 중앙값 (반복 측정 병합)",
           "design": "역가 3구간 × 실데나필 Tanimoto 3구간 = 9칸, 칸당 균등 추출"}
    checks = [
        ("9칸 모두 채워짐", all(c["taken"] > 0 for c in cells)),
        ("칸별 최소 10건", all(c["taken"] >= 10 for c in cells)),
        ("역가 범위 >= 3 로그", (max(pv) - min(pv)) >= 3.0),
        ("교란 상관 |r| < 0.35 (이전 설계 대비 감소)", abs(confound_r) < 0.35),
        ("중복 없음", len({r["molecule_chembl_id"] for r in picked}) == n),
    ]
    env = make_result(out, "ChEMBL webresource_client + RDKit Morgan Tanimoto",
                      f"target={TARGET}, 9칸 격자, 칸당 {a.per_cell}", checks,
                      notes=(f"n={n}, pChEMBL {min(pv):.2f}–{max(pv):.2f}, "
                             f"Tanimoto {min(sims):.3f}–{max(sims):.3f}. "
                             f"역가-유사도 교란 상관 r={confound_r}."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\nn={n}  pChEMBL {min(pv):.2f}–{max(pv):.2f}  Tanimoto {min(sims):.3f}–{max(sims):.3f}")
    print(f"교란 상관 (역가 vs 유사도) r = {confound_r}   게이트 "
          f"{'PASS' if env['verification']['passed'] else 'FAIL'}")
    for c in env["verification"]["checks"]:
        print(f"  {'OK ' if c['passed'] else 'NG '} {c['check']}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
