---
name: figure-builder
description: 보고서용 그림을 만든다 — graphical abstract, 전체 모식도, 데이터 그림. 수치가 들어가는 그림은 도구 산출값만 사용한다.
---
# figure-builder

## 언제 사용하는가
REPORT 직전. `report-writer` 가 본문을 쓰기 전에 그림을 먼저 확정한다
(본문이 그림을 참조하므로 순서를 지킨다).

## 무엇을 만드는가

| 파일 | 종류 | 데이터 출처 |
|------|------|-------------|
| `fig1_graphical_abstract.png` | 개념도 | 수치 미포함 — 질문 → 방법 → 발견 → 함의 |
| `fig2_pipeline.png` | 전체 모식도 | 하네스 단계와 게이트 구조 (수치 미포함) |
| `fig3_property_space.png` | 산점도 | `sample_run/run_stdout.json` 물성 실계산값 |
| `fig4_gate_waterfall.png` | 폭포도 | 단계별 잔존 건수 (실측) |
| `fig5_qed_threshold.png` | 분포+임계선 | QED 실계산값과 `QED_MIN` |

## 어떻게 (호출 명령)
```bash
python scripts/make_figures.py --run sample_run/run_stdout.json --out outputs/figures
```

## 하드 규칙 (무-날조)
- **수치가 들어간 그림은 `run_stdout.json` 의 값만** 쓴다. `np.random.*` 금지.
- 데이터가 없는 패널은 **빈 패널 + 사유 텍스트**로 남긴다. 그럴듯한 곡선을 그리지 않는다.
- graphical abstract 와 모식도는 개념도이므로 수치를 넣지 않는다 (넣으면 근거 필요).
- 색은 colorblind-safe (Okabe-Ito) 팔레트, 300 DPI, 축·단위·n 표기 필수.

## 반환 검증
- 파일 5개 존재 + 각 PNG 옆에 `<name>.meta.json` (source, script, n, generated_at).
- 수치 그림의 n 이 `run_stdout.json` 의 레코드 수와 일치.
