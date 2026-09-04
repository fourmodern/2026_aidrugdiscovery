# sample_run — 실제 실행 결과 (읽기용)

`outputs/` 는 실행 산출물이라 저장소에서 제외된다. 강의에서 **실습 없이 결과만 읽을** 수 있도록,
한 번의 실제 실행 결과를 여기에 고정해 둔다.

- **실행 명령**: `python run_harness.py run`
- **실행 일시**: 2026-09-04 (UTC+9)
- **환경**: Python 3.12.3 · rdkit · chembl-webresource-client · requests

| 파일 | 내용 |
|------|------|
| `run_gatelog.txt` | stderr — 단계별 `[VERIFY-GATE:<step>] PASSED/FAILED` 판정 4줄 |
| `run_stdout.json` | stdout — 각 단계가 돌려준 값 (UniProt 표적 봉투, RDKit 물성 10건) |

## 읽는 순서

1. `CLAUDE.md` 의 **단계별 검증 기준** 표를 먼저 읽는다 (무엇을 통과해야 하는가).
2. `run_gatelog.txt` 로 그 기준이 실제로 판정되었는지 본다 (4단계 모두 PASSED).
3. `run_stdout.json` 에서 값의 출처를 확인한다 — `accession: O76074` 는 UniProt 실조회,
   `MW/logP/QED` 는 RDKit 실계산이다. 모델이 적은 수치는 하나도 없다.

## 이 실행에서 관찰된 것

물성 게이트는 `Ro5 위반 1개까지 허용 AND QED >= 0.5 AND SA <= 6.0` 이다 (`scripts/mol_properties.py`).
10건 중 4건이 통과했는데, **10건 전부 Ro5 를 통과했고 SA 도 최대 2.97** 이었다.
즉 세 조건 중 **QED 하나만으로** 통과·탈락이 갈렸다. 임계 0.5 는 스크립트가 스스로
"데모 임계(교육용)" 라고 적어 둔 값이며, 목적에 맞게 정해야 하는 숫자다.
