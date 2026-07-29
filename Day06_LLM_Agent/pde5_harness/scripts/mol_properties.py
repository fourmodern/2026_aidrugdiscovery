"""RDKit 물성/약물성 계산 + Lipinski·QED·SA 게이트. 무-날조: 실계산값만.

사용:
  python scripts/mol_properties.py "SMILES" ["SMILES2" ...]
  echo '<chembl_actives.py 출력 JSON>' | python scripts/mol_properties.py --stdin
출력: 표준 결과 봉투(JSON) — 후보별 MW/logP/HBD/HBA/TPSA/QED/SA + 게이트 통과여부.
"""
from __future__ import annotations
import sys, json
from verify import make_result, valid_smiles

QED_MIN = 0.5   # 데모 임계(교육용) — 실제 프로젝트는 목적에 맞게 조정
SA_MAX = 6.0


def _sa_score(mol):
    try:
        from rdkit.Chem import RDConfig
        import os
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer
        return round(sascorer.calculateScore(mol), 2)
    except Exception:
        return None  # SA 스크립트 없으면 None(날조 금지)


def props(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, Lipinski, Crippen, rdMolDescriptors
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    mw = round(Descriptors.MolWt(m), 1)
    logp = round(Crippen.MolLogP(m), 2)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)
    tpsa = round(rdMolDescriptors.CalcTPSA(m), 1)
    qed = round(QED.qed(m), 3)
    sa = _sa_score(m)
    ro5 = sum([mw <= 500, logp <= 5, hbd <= 5, hba <= 10])  # Lipinski 위반 허용 1
    passed = (ro5 >= 3) and (qed >= QED_MIN) and (sa is None or sa <= SA_MAX)
    return {"smiles": smiles, "MW": mw, "logP": logp, "HBD": hbd, "HBA": hba,
            "TPSA": tpsa, "QED": qed, "SA": sa, "Ro5_pass": ro5 >= 3, "gate_pass": passed}


def collect_smiles():
    if len(sys.argv) > 1 and sys.argv[1] == "--stdin":
        env = json.load(sys.stdin)
        rows = env.get("result", env) if isinstance(env, dict) else env
        return [r.get("canonical_smiles") or r.get("smiles") for r in rows]
    return sys.argv[1:]


def main():
    try:
        import rdkit  # noqa
    except Exception:
        env = make_result([], "RDKit", "n/a", [("RDKit 설치", False)],
                          notes="RDKit 미설치 — pip install rdkit 후 재실행(데모 수치 날조 금지).")
        print(json.dumps(env, ensure_ascii=False, indent=2)); return
    smis = [s for s in collect_smiles() if s]
    out, sane = [], True
    for s in smis:
        if not valid_smiles(s):
            sane = False; continue
        p = props(s)
        if p is None:
            sane = False; continue
        # sanity 범위 점검
        if not (0 < p["MW"] < 2000 and -10 < p["logP"] < 15 and 0 <= p["QED"] <= 1):
            sane = False
        out.append(p)
    checks = [
        ("후보 존재", len(out) > 0),
        ("값 sanity 범위 내", sane and len(out) > 0),
    ]
    npass = sum(1 for p in out if p["gate_pass"])
    env = make_result(out, "RDKit Descriptors/QED", f"n_smiles={len(smis)}", checks,
                      notes=f"물성 계산 {len(out)}건, 게이트(Ro5≥3·QED≥{QED_MIN}·SA≤{SA_MAX}) 통과 {npass}건.")
    print(json.dumps(env, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
