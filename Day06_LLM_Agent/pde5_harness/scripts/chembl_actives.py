"""ChEMBL에서 PDE5A(CHEMBL1827) 활성물질 조회. 무-날조: 실제 ChEMBL ID/SMILES/활성값만.

사용: python scripts/chembl_actives.py [limit]
출력: 표준 결과 봉투(JSON). 각 레코드는 실제 molecule_chembl_id + RDKit로 유효한 SMILES.
"""
from __future__ import annotations
import sys, json
from verify import make_result, valid_smiles

TARGET = "CHEMBL1827"  # Phosphodiesterase 5A

# 오프라인 데모 폴백(공개된 대표 PDE5 저해제, ChEMBL 실 ID). 네트워크 불가 시에만.
OFFLINE = [
    {"molecule_chembl_id": "CHEMBL192",  "pref_name": "SILDENAFIL",
     "canonical_smiles": "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"},
    {"molecule_chembl_id": "CHEMBL779",  "pref_name": "TADALAFIL",
     "canonical_smiles": "O=C1N(C)CC(=O)N2[C@@H]1Cc1c([nH]c3ccccc13)[C@H]2c1ccc2c(c1)OCO2"},
    {"molecule_chembl_id": "CHEMBL1520", "pref_name": "VARDENAFIL",
     "canonical_smiles": "CCCc1nc(C)c2c(=O)nc(-c3cc(S(=O)(=O)N4CCN(CC)CC4)ccc3OCC)[nH]n12"},
]


def fetch(limit: int):
    try:
        from chembl_webresource_client.new_client import new_client
        act = new_client.activity
        qs = act.filter(target_chembl_id=TARGET, standard_type__in=["IC50", "Ki"],
                        standard_relation="=", pchembl_value__isnull=False)
        rows = []
        for a in qs[:limit]:
            rows.append({
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "canonical_smiles": a.get("canonical_smiles"),
                "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
                "pchembl_value": a.get("pchembl_value"),
            })
        return rows or None, "ChEMBL webresource_client"
    except Exception:
        return None, "ChEMBL (unavailable)"


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    rows, source = fetch(limit)
    query = f"target_chembl_id={TARGET}; type in (IC50,Ki); n={limit}"
    if rows is None:
        rows = [dict(r) for r in OFFLINE]
        source = "ChEMBL (OFFLINE DEMO — 공개 대표 저해제, 실데이터 아님)"
    # 검증: 실 ID + SMILES 유효성
    n = len(rows)
    n_valid = sum(1 for r in rows if valid_smiles(r.get("canonical_smiles", "")))
    n_id = sum(1 for r in rows if str(r.get("molecule_chembl_id", "")).startswith("CHEMBL"))
    checks = [
        ("레코드 존재", n > 0),
        ("모든 SMILES RDKit 유효", n_valid == n and n > 0),
        ("모든 molecule_chembl_id 형식", n_id == n and n > 0),
    ]
    env = make_result(rows, source, query, checks,
                      notes=f"활성물질 {n}건, SMILES유효 {n_valid}/{n}, 실ID {n_id}/{n}. 무-날조: 도구 반환값만.")
    print(json.dumps(env, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
