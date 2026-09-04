#!/usr/bin/env python3
"""활성-물성 통합 분석 — 리간드 효율 지표를 실측값으로 계산한다.

pChEMBL 과 RDKit 기술자를 결합해 의약화학에서 표준으로 쓰는 두 효율 지표를 낸다.
  LE  = 1.37 x pIC50 / heavy atom count      (원자당 결합 효율)
  LLE = pIC50 - cLogP                        (지용성 대비 효율)
어느 값도 추정하지 않는다. pIC50 은 ChEMBL 보고값, cLogP·중원자수는 RDKit 계산값이다.
"""
import json, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load_envelope(name):
    return json.loads((ROOT / "sample_run" / "envelopes" / name).read_text())


def main():
    ch = load_envelope("chembl-actives.json")["result"]
    mp = load_envelope("mol-properties.json")["result"]
    by_smiles = {r["smiles"]: r for r in mp}

    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        sys.exit("RDKit 이 필요합니다 — 중원자수를 계산할 수 없어 중단합니다.")

    rows = []
    for a in ch:
        smi = a.get("canonical_smiles")
        prop = by_smiles.get(smi)
        if prop is None:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        hac = mol.GetNumHeavyAtoms()
        p = float(a["pchembl_value"])
        rows.append({
            "chembl_id": a["molecule_chembl_id"],
            "smiles": smi,
            "assay_type": a.get("standard_type"),
            "ic50_nM": float(a["standard_value"]),
            "pIC50": p,
            "heavy_atoms": hac,
            "MW": prop["MW"], "cLogP": prop["logP"], "TPSA": prop["TPSA"],
            "HBD": prop["HBD"], "HBA": prop["HBA"],
            "QED": prop["QED"], "SA": prop["SA"],
            "Ro5_pass": prop["Ro5_pass"],
            # 표준 정의 그대로. 반올림은 표시용으로만 한다.
            "LE": round(1.37 * p / hac, 3),
            "LLE": round(p - prop["logP"], 2),
        })
    rows.sort(key=lambda r: -r["pIC50"])

    out = ROOT / "sample_run" / "sar.json"
    out.write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n")

    print(f"{out}  ({len(rows)}건)")
    print(f"{'ChEMBL ID':14s} {'IC50(nM)':>9s} {'pIC50':>6s} {'HAC':>4s} "
          f"{'MW':>7s} {'cLogP':>6s} {'LE':>6s} {'LLE':>6s}")
    for r in rows:
        print(f"{r['chembl_id']:14s} {r['ic50_nM']:9.1f} {r['pIC50']:6.2f} "
              f"{r['heavy_atoms']:4d} {r['MW']:7.1f} {r['cLogP']:6.2f} "
              f"{r['LE']:6.3f} {r['LLE']:6.2f}")
    pic = [r["pIC50"] for r in rows]
    le = [r["LE"] for r in rows]; lle = [r["LLE"] for r in rows]
    print(f"\npIC50 {min(pic):.2f}–{max(pic):.2f}  |  "
          f"LE {min(le):.3f}–{max(le):.3f}  |  LLE {min(lle):.2f}–{max(lle):.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
