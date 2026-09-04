---
name: antibody-search
description: 해당 항원에 결합하는 알려진 항체 구조를 RCSB PDB 에서 검색하고 실제 heavy/light 사슬 서열을 가져온다.
---
# antibody-search

## 언제 사용하는가
antigen-lookup 통과 후 2단계. 새로 설계하기 전에 **이미 알려진 항체를 벤치마크로 확보**한다.
설계 결과를 "트라스투주맙 대비" 평가하려면 이 단계가 선행되어야 한다.

## 어떻게 (호출 명령)
```bash
python scripts/antibody_search.py                    # 기본 P04626, 최대 8 entry
python scripts/antibody_search.py P04626 6
python scripts/antibody_search.py P04626 6 --offline # 오프라인 캐시 강제 (교육용)
```
1. **RCSB Search API** (`search.rcsb.org/rcsbsearch/v2/query`) 로
   `reference_sequence_identifiers.database_accession == <acc>` AND `database_name == UniProt`
   AND `polymer_entity_count_protein >= 2` 인 experimental entry 검색.
2. 각 entry 를 **Data API** (`data.rcsb.org/rest/v1/core/entry|polymer_entity`) 로 조회해
   `entity_poly.pdbx_seq_one_letter_code_can` 원값을 가져온다.
3. 사슬을 `antigen` / `heavy` / `light` / `scfv` / `unknown` 으로 분류
   (UniProt 참조 일치 → antigen, 그 외는 보존 Ig 모티프 휴리스틱).
4. 항원 + 항체 사슬을 모두 가진 entry 만 남긴다.

## 반환 검증
- checks: 검색 성공 / 복합체 ≥1 / PDB ID 형식 유효 / 모든 사슬 서열 유효 / 중쇄·경쇄 각 ≥1.
- 각 사슬의 `classification_evidence` 로 왜 heavy/light 로 판정했는지 확인 (FR4 모티프 위치 등).
- `chain_type == "scfv"` 는 단일 사슬에 VH·VL 이 함께 있다는 뜻 — cdr-analysis 가 자동 분할한다.

## 무-날조
- 검색·조회 실패 → `result: []` + `passed: false`. 서열을 만들지 않는다.
- 오프라인 폴백은 **P04626 에 한해** PDB 1N8Z(트라스투주맙 Fab) 실서열 캐시를 쓴다.
  `provenance.source` 가 `offline cache (PDB 1N8Z, ...)` 이고 `passed=false` 이므로
  **게이트 통과가 아니다.** 보고서에 "실시간 조회 아님"을 반드시 표기하고 네트워크 복구 후 재실행.
