---
name: devils-advocate
description: "왜 이 결과가 틀렸을 수 있는가"만 본다. 10항목 전부 CLEAR 여야 통과. 1개라도 UNCLEAR 면 FAIL.
---
# devils-advocate

## 언제 사용하는가
`pi-review` 통과 직후. 긍정 평가는 이미 끝났으므로 **반증만** 본다.

## 체크리스트 (전부 CLEAR 여야 10/10)
1. **확증편향** — 반박 증거를 찾았는가. null 가능성이 논의에 있는가.
2. **체리피킹** — 대상 선정이 체계적인가. 제외 사유가 타당한가.
3. **대안 가설** — 경쟁 설명 2개 이상과 구별 방법.
4. **인과 vs 상관** — 관찰 데이터로 "유발한다"고 쓰지 않았는가.
5. **표본/선택 편향** — 코호트가 대표성 있는가.
6. **통계 무결성** — 다중검정 보정, 효과크기, 사전등록, HARKing 아님.
7. **재현성** — 코드·데이터 공개, seed 고정, 버전 명시.
8. **일반화 한계** — 적용 범위를 넘어선 주장이 없는가.
9. **임상 번역 격차** — in silico 결과를 임상 효능처럼 쓰지 않았는가.
10. **윤리** — 데이터 라이선스, 이해상충, 오남용 가능성.

## 외부 평가 의무
`Agent(subagent_type="oh-my-claudecode:critic")` 를 **매 iteration 새로** 호출한다.
자가 채점 금지. UNCLEAR 항목은 근거 문장을 인용해 지적한다.

## 반환 형식
```json
{"score": 10, "items": {"1_confirmation_bias": "CLEAR", "...": "UNCLEAR"},
 "unclear_reasons": ["..."], "required_fixes": ["..."], "pass": true}
```
