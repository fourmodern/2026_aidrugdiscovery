---
name: antigen-lookup
description: 항원(표적) 단백질의 실재·기능·세포외 도메인을 UniProt REST 로 확인한다. 모든 후속 단계의 첫 게이트.
---
# antigen-lookup

## 언제 사용하는가
파이프라인 첫 단계. 항체를 찾거나 설계하기 전에 **항원이 실재하고 어떤 단백질인지** 확인한다.
기본 표적은 HER2/ERBB2 (UniProt P04626).

## 어떻게 (호출 명령)
```bash
python scripts/antigen_lookup.py            # 기본 P04626 (HER2)
python scripts/antigen_lookup.py P04626     # 명시적
python scripts/antigen_lookup.py P01375     # 다른 표적 (TNF-alpha)
```
`https://rest.uniprot.org/uniprotkb/<acc>.json` 을 조회해 accession / gene / protein name /
organism / length / sequence / FUNCTION / subcellular location / topology feature 를 추출한다.

## 반환 검증
표준 봉투 `{result, provenance, verification}`. 다음을 확인:
- `verification.passed == true` — 아니면 다음 단계 진행 금지 (최대 2회 재시도 후 "미확인" 플래그).
- checks: HTTP 200 / accession 일치 / gene name / protein name / 서열 아미노산 알파벳 유효 /
  `length` 와 실제 서열 길이 일치 / FUNCTION 서술 존재.
- `result.topology_features` 에서 `Topological domain: Extracellular` 구간을 확인 —
  **항체가 접근할 수 있는 영역**이며, 설계 시 에피토프는 이 범위 안에서 골라야 한다.
  (HER2 는 signal 1-22, extracellular 23-652, TM 653-675, cytoplasmic 676-1255.)

## 무-날조
네트워크/HTTP 실패 시 **오프라인 폴백이 없다.** `result: null` + `passed: false` 로 정직 보고한다.
길이·서열·기능을 기억에서 채워 넣지 말 것.
