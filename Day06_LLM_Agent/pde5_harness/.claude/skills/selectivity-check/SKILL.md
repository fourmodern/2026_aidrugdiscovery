---
name: selectivity-check
description: PDE5 vs PDE6 선택성을 정직하게 평가한다. 정량 예측 없으면 정성·확인 필요로 보고(수치 날조 금지).
---
# selectivity-check

## 언제 사용하는가
물성 평가(mol-properties) 후, 실행 단계 (d). PDE5 저해제의 표적-외(off-target) 위험,
특히 망막 PDE6 교차 억제에 따른 시각 부작용 맥락을 평가한다.

## 어떻게 (호출 명령)
```bash
python scripts/selectivity.py
python scripts/selectivity.py --smiles "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"
```
- 이 도구는 **정량 fold-selectivity를 예측하지 않는다.** PDE5/PDE6 아이소자임 기능 맥락,
  교차 억제 기전(망막 PDE6 억제 → 일시적 청녹 시각이상/cyanopsia), 실제 선택성 확인에
  필요한 후속 근거(병렬 IC50 어세이 / 인용 문헌)를 반환한다.

## 반환 검증
표준 봉투를 반환. 다음을 확인:
- checks: `정량 fold-selectivity 미날조(정성 보고)`, `교차 억제 기전(PDE6 망막) 서술 존재`,
  `후속 확인 근거 목록 존재`, `참고문헌 실재(Boolell 1996 / Ghofrani 2006)`.
- `result.quantitative_selectivity_predicted == false` — 정량 수치를 만들지 않았음을 확인.
- 무-날조: 화합물별 선택성 비(예: "PDE5/PDE6 = 10배")는 병렬 어세이 또는 인용 문헌 없이 기재 금지.
  불확실은 "확인 권장"으로 표기.
