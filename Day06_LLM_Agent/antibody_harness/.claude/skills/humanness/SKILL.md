---
name: humanness
description: 인간 germline V 유전자와의 서열 동일성(%)을 실계산한다. 발명된 "휴먼성 점수"를 만들지 않는다.
---
# humanness

## 언제 사용하는가
항체 사슬 서열이 확보된 뒤. **알려진 항체와 설계 항체 모두**에 적용해 비교한다.

## 어떻게 (호출 명령)
```bash
python scripts/antibody_search.py | python scripts/humanness.py --stdin
python scripts/humanness.py "EVQLVESGG..." --top 5
```

## 방법
1. 입력 사슬에서 가변영역(V-domain) 경계를 추정 (`seq_utils.variable_domain`, 휴리스틱).
2. **UniProt REST** 에서 human germline V 서열을 조회:
   `(gene:IGHV*|IGKV*|IGLV*) AND (organism_id:9606) AND (reviewed:true) AND (fragment:false)`.
   결과는 `outputs/cache/human_germline_v.json` 에 캐시된다.
3. BioPython `PairwiseAligner` (BLOSUM62, **local**, gap open -11 / extend -1) 로 정렬.
4. `germline_identity_percent = 100 × 동일 잔기 수 / 정렬 길이`.
5. 중쇄는 IGHV pool, 경쇄는 IGKV+IGLV pool 과만 비교한다.

## 반환 검증
- checks: germline pool ≥50 / 도메인 ≥1 / 모든 도메인에 nearest germline / identity 0-100 / 스킵 0건.
- `nearest_germline` 에 germline 유전자명 + UniProt accession + identity + 정렬 길이가 들어 있다.

## 무-날조 (중요)
- **"휴먼성 점수 8.5/10" 같은 발명 지표를 만들지 말 것.** 보고 가능한 값은
  `germline_identity_percent` 하나뿐이며, 정의(`metric_definition`)를 항상 함께 인용한다.
- **local alignment 정의상** germline 의 signal peptide 와 항체의 CDR3/FR4 는 정렬 밖이다.
  → 이 값은 프레임워크 유사도에 가깝게 편향된다. 이 사실을 보고서에 쓸 것.
- germline 조회 실패 + 캐시 없음 → `result: []` + `passed: false`. 값을 만들지 않는다.
- **이 지표는 면역원성(ADA) 예측이 아니다.** 상관 지표일 뿐이며 T-cell epitope / MHC-II 결합 /
  응집 / 투여경로·용량은 반영되지 않는다. `limitations` 필드를 그대로 인용하라.
