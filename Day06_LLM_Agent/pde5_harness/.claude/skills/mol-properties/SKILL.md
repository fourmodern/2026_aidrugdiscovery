---
name: mol-properties
description: RDKit로 후보 화합물의 물성·약물성(MW/logP/HBD/HBA/TPSA/QED/SA)과 Lipinski 게이트를 계산한다.
---
# mol-properties

## 언제 사용하는가
활성물질 조회(chembl-actives) 후, 실행 단계 (c). 후보의 물리화학적 성질과 약물유사성을 정량 평가해
선별 게이트(Ro5·QED·SA)를 적용한다.

## 어떻게 (호출 명령)
```bash
# chembl-actives 출력을 파이프(권장)
python scripts/chembl_actives.py | python scripts/mol_properties.py --stdin

# 개별 SMILES 직접 평가
python scripts/mol_properties.py "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"
```
- RDKit `Descriptors/QED/Lipinski/Crippen/rdMolDescriptors` 로 실계산. SA는 RDKit contrib(sascorer) 있으면 계산, 없으면 `null`(날조 금지).
- 게이트: `Ro5 위반 ≤1 (≥3 충족)` AND `QED ≥ 0.5` AND `SA ≤ 6.0`(있을 때).

## 반환 검증
표준 봉투를 반환. 다음을 확인:
- `verification.passed == true` (checks: `후보 존재`, `값 sanity 범위 내`).
- 각 후보의 `gate_pass`(선별 통과 여부)와 `Ro5_pass` 확인.
- RDKit 미설치 시: `RDKit 설치=false` 로 정직하게 실패 보고(수치 날조 금지) → 안내대로 `pip install rdkit`.
- 물리적 범위 sanity: 0<MW<2000, -10<logP<15, 0≤QED≤1.
