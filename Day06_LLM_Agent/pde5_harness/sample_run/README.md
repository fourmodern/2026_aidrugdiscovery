# sample_run — 실제 실행 결과와 심사 기록 (읽기용)

`outputs/` 는 실행 산출물이라 저장소에서 제외된다. 강의에서 **실습 없이 파일과 결과만
읽을** 수 있도록, 한 번의 실제 실행과 그에 대한 심사 기록을 여기에 고정한다.

- **실행**: `python run_harness.py run` · 2026-09-04 · Python 3.12.3

| 파일 | 내용 |
|------|------|
| `run_gatelog.txt` | 단계별 `[VERIFY-GATE:<step>] PASSED/FAILED` 4줄 |
| `run_stdout.json` | 표적 봉투와 RDKit 물성 10건 |
| `envelopes/*.json` | 단계별 봉투 전문 (result + provenance + verification) |
| `negative_control_gatelog.txt` | **실패 봉투를 주입해 차단을 관측한 로그** |
| `negative_control.txt` | 음성 대조 실행 출력 |
| `references.md` | 보고서 참고문헌 20건 (Crossref·NCBI API 로 실재 검증) |
| `REVIEW_LOG.md` | **1차 심사에서 드러난 하네스 결함과 조치** |

## 읽는 순서

1. `CLAUDE.md` 의 단계별 검증 기준 표 — 무엇을 통과해야 하는가.
2. `run_gatelog.txt` — 4단계 전부 PASSED.
3. **여기서 멈추지 말 것.** `REVIEW_LOG.md` 를 읽는다.
   통과 로그만으로는 게이트가 강제되는지 알 수 없다는 것이 이 자료의 요점이다.
4. `negative_control_gatelog.txt` — 실패를 주입했을 때 비로소 차단이 관측된다.
5. `envelopes/*.json` — 값의 출처(provenance)와 판정 근거(verification)를 확인한다.

## 이 실행에서 관찰된 것

물성 게이트는 `Ro5 위반 1개까지 허용 AND QED >= 0.5 AND SA <= 6.0` 이다.
10건 중 4건이 통과했는데, **10건 전부가 Ro5 를 통과**하고 SA 도 최대 2.97 이라
세 조건 중 QED 하나가 결과를 전부 갈랐다.

주의: Ro5 가 걸러내지 못한 이유는 위반이 드물어서가 **아니다**.
실제로는 10건 중 7건이 MW 500 을 넘고, 위반이 0개인 화합물은 3건뿐이다.
게이트가 위반 1개를 허용하기 때문에 통과한 것이다. 초판 보고서는 이를 반대로
해석했다가 심사에서 지적받아 정정했다.

## 이 자료의 요점

게이트 통과 로그는 그 자체로 강제의 증거가 아니다. 1차 심사에서 selectivity 게이트가
어떤 입력으로도 실패할 수 없는 항진명제였음이 드러났다. 실패가 관측 가능해야
통과가 정보를 갖는다. 자세한 내용은 `REVIEW_LOG.md`.
