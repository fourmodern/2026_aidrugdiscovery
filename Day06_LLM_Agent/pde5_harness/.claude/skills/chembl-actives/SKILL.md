---
name: chembl-actives
description: ChEMBL(CHEMBL1827)에서 PDE5A 활성물질(실 ID·SMILES·활성값)을 조회한다. 물성 평가의 입력.
---
# chembl-actives

## 언제 사용하는가
표적 검증(target-lookup) 통과 후, 실행 단계 (b). PDE5A에 대해 보고된 활성 화합물을 확보해
다음 물성 평가(mol-properties)의 입력으로 넘긴다.

## 어떻게 (호출 명령)
```bash
python scripts/chembl_actives.py           # 기본 25건
python scripts/chembl_actives.py 50        # limit 지정
```
- `chembl_webresource_client` 로 `target_chembl_id=CHEMBL1827`, type in (IC50, Ki), pchembl 존재 레코드 조회.
- 각 레코드: `molecule_chembl_id`, `canonical_smiles`, `standard_type/value/units`, `pchembl_value`.
- 네트워크·라이브러리 불가 시 오프라인 폴백(sildenafil CHEMBL192 / tadalafil CHEMBL779 / vardenafil CHEMBL1520 — 실 ID, "실데이터 아님" 표기).

## 반환 검증
표준 봉투를 반환. 다음을 확인:
- `verification.passed == true` (아니면 재시도 2회 후 플래그).
- checks: `레코드 존재`, `모든 SMILES RDKit 유효`, `모든 molecule_chembl_id 형식(CHEMBL 접두)`.
- 무-날조: SMILES는 RDKit로 파싱되는 것만, ChEMBL ID는 도구 반환값만. LLM이 ID·활성값을 지어내지 말 것.
- 다음 단계로 파이프: `python scripts/chembl_actives.py | python scripts/mol_properties.py --stdin`.
