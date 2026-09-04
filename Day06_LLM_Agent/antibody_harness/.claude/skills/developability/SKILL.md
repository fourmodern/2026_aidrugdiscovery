---
name: developability
description: BioPython ProtParam 실계산 물성 + 규칙 기반 liability 모티프 스캔으로 항체 사슬의 개발 가능성 위험을 플래그한다.
---
# developability

## 언제 사용하는가
CDR 이 확보된 뒤. **알려진 항체와 설계 항체 모두**에 같은 자로 적용한다.

## 어떻게 (호출 명령)
```bash
python scripts/antibody_search.py | python scripts/developability.py --stdin
python scripts/developability.py "EVQLVESGG..."
python scripts/developability.py chains.fasta
```

## 무엇을 계산하는가

### (a) 물성 — BioPython `ProteinAnalysis` 실계산
`molecular_weight_Da`, `isoelectric_point_pI`, `gravy_kyte_doolittle`,
`instability_index` (+ Guruprasad 1990 기준 >40 = unstable 판정),
`aromaticity`, `extinction_coefficient` (환원/시스틴 2값), 2차구조 분율.
비표준 잔기(X/B/Z/U/O)는 제외하고 `n_nonstandard_residues_excluded` 에 기록한다.

### (b) Liability 모티프 — 정규식 규칙 (예측 모델 아님)
| 이름 | 규칙 | 심각도 |
|------|------|--------|
| N-glycosylation sequon | `N[^P][ST]` | high |
| deamidation NG | `NG` | high |
| deamidation moderate | `N[STNH]` | moderate |
| isomerization DG | `DG` | high |
| isomerization moderate | `D[STD]` | moderate |
| acid hydrolysis / clip | `DP` | moderate |
| oxidation | `M`, `W` | moderate |

구조적 liability: 홀수 개 Cys(유리 Cys 가능), N-말단 Q/E(pyroglutamate).
응집 프록시: Kyte-Doolittle 5-잔기 이동평균 ≥ 1.5 인 소수성 패치.

각 히트는 `position_1based` 와 **`in_cdr`**(어느 CDR 내부인지)를 함께 보고한다.
**CDR 내부 liability 가 가장 중요하다** — 항원 결합에 직접 영향을 준다.

## 반환 검증
- checks: 사슬 ≥1 / 물성 물리적 범위(0<MW, 0<pI<14, -5<GRAVY<5) / 스킵 0건.
- `summary.n_in_cdr` 로 CDR 내부 위험 개수를 먼저 본다.

## 무-날조
- BioPython 미설치 시 수치를 만들지 않고 `passed: false` 로 보고한다.
- liability 는 **규칙 기반 플래그**다. "이 항체는 응집한다"가 아니라
  "NG 탈아미드화 모티프가 CDR-H2 55번에 있다(규칙 기반 플래그)"로 서술하라.
  실제 위험도는 pH·온도·제형·구조 노출도에 따라 달라지며, 실험(가속 안정성·HIC·SEC)으로만 확인된다.
