---
name: compare-designs
description: 경로 A(ESMFold2 inversion)와 경로 B(RFantibody) 산출물을 동일한 서열 기반 지표로 비교한다. 통합 점수를 만들지 않는다.
---
# compare-designs

## 언제 사용하는가
두 설계 경로 중 하나 이상이 결과를 낸 뒤. 한쪽만 있어도 동작한다(없는 쪽은 null).

## 어떻게 (호출 명령)
```bash
python scripts/compare_designs.py \
    --track-a outputs/design_a.json \
    --track-b outputs/design_b.json \
    --out outputs/06_comparison.json

python scripts/compare_designs.py --track-a outputs/design_a.json   # 한쪽만
```

## 무엇을 하는가
1. 두 봉투에서 설계 서열을 추출한다.
2. 서열이 있는 설계에 대해 **기존 평가 3종을 그대로 재사용**한다
   (`cdr_analysis.py`, `developability.py`, `humanness.py` 를 `--stdin` 으로 호출).
3. 공통 축과 경로별 축을 **분리해서** 표로 낸다.

## 공통 비교축 (동일 방법으로 계산 — 직접 비교 정당)
`chain_type`, `cdrs` / `cdr_lengths`, `molecular_weight_Da`, `isoelectric_point_pI`,
`gravy_kyte_doolittle`, `instability_index`, `n_liability_hits`, `n_liabilities_in_cdr`,
`n_hydrophobic_patches`, `germline_identity_percent`, `nearest_germline`.

## 경로 고유 지표 (원본 이름 유지 — 교차 비교 금지)
- 경로 A: `iptm`, `distogram_iptm_proxy`, `cdr_distogram_iptm_proxy`, `final_loss`, `critic_name`
- 경로 B: RF2 score 파일 컬럼 그대로(pAE 등), `ca_rmsd_self_consistency_A`

## 절대 하면 안 되는 것 (무-날조)
- **경로 A 의 ipTM 과 경로 B 의 pAE 를 하나의 "종합 점수"로 합치지 말 것.**
  서로 다른 모델의 서로 다른 척도이며, 통합 점수는 근거 없는 발명이다.
  스크립트의 `interpretation_rules` 필드에 이 규칙이 함께 출력되니 그대로 인용하라.
- 한쪽 경로 결과가 없는 칸을 채워 넣지 말 것 (null 유지).
- RFantibody score 파일에 서열 컬럼이 없으면 그 설계는 공통 지표를 계산할 수 없다
  (`sequence_available: false`). 서열을 만들어내지 말 것.

## 반환 검증
checks: 최소 한 경로 로드 / 서열 있는 설계 ≥1 / 공통 평가 3종 성공.
`evaluation_errors` 가 비어 있는지 확인한다.

## 결론 서술 규약
- germline identity 와 liability 플래그는 **개발 가능성 위험 신호**이지 결합력·면역원성 예측이 아니다.
- 어느 경로도 wet-lab 을 대체하지 않는다. 실측(발현·SPR/BLI·DSF·SEC)만이 hit 여부를 결정한다.
- RFantibody 저자가 명시한 "신뢰할 만한 필터의 부재"와 "10k 규모 캠페인" 전망을 함께 인용한다.
