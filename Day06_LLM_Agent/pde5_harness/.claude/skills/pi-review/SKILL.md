---
name: pi-review
description: 보고서를 PI 수준으로 채점한다. 9차원 가중평균 9.7 이상 AND 전 차원 9.0 이상이어야 통과. 자가 통과 금지.
---
# pi-review

## 언제 사용하는가
REPORT 산출 직후 첫 번째 게이트.

## 통과 알고리즘 (둘 다 만족)
1. 가중평균 **≥ 9.7**
2. **어떤 차원도 9.0 미만이 아닐 것** — 하나라도 미만이면 평균과 무관하게 FAIL

## 채점 차원 (가중치 합 100%)
| # | 차원 | 가중 |
|---|------|------|
| 1 | Significance | 20% |
| 2 | Innovation | 15% |
| 3 | Approach — Rigor | 15% |
| 4 | Approach — Feasibility | 10% |
| 5 | Pitfalls / Alternatives | 10% |
| 6 | Preliminary Data | 5% |
| 7 | Writing / Clarity | 10% |
| 8 | Ethics / Reproducibility | 10% |
| 9 | Stats (보고서) | 5% |

## 외부 평가 의무
채점은 **작성자가 아닌 외부 리뷰어**가 한다. 자동화 환경에서는
`Agent(subagent_type="oh-my-claudecode:critic")` 과 `Agent(subagent_type="oh-my-claudecode:scientist")`
를 **매 iteration 새로 호출**하고, 두 결과의 평균으로 판정한다. 이전 iteration 결과 재사용 금지.

## 반환 형식
```json
{"weighted_score": 9.72, "per_dim": {"Significance": 9.8, "...": 9.5},
 "dims_below_9_0": [], "required_changes": ["..."], "pass": true}
```

## 기록
`checklists/pi_review_<iter>.md` 에 점수·근거·요구 수정사항을 남긴다.
