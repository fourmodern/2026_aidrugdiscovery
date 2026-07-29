# PDE5 저해제 탐색 에이전트 하네스

[![Claude Code](https://img.shields.io/badge/Claude%20Code-harness-8A2BE2)](https://docs.claude.com/en/docs/claude-code)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-cheminformatics-green.svg)](https://www.rdkit.org/)
[![No-Fabrication](https://img.shields.io/badge/policy-no--fabrication-red.svg)](#하드-규칙-무-날조)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)

> 6일차(LLM·에이전트) 종합 실습 — **강사 데모(B 방식)**용 Claude Code 하네스.
> 자연어 목표를 받아 **계획 → 단계별 검증 게이트 → 근거 인용 보고서**까지, 과학자의
> 엄밀함으로 PDE5 저해제 후보를 조사·평가한다. 수치·식별자는 **도구 실계산값만** 사용한다.

---

## 개요

이 하네스는 "에이전트가 어떻게 과학적 엄밀함(무-날조·provenance·단계별 검증)을 유지하며
연구 파이프라인을 자동 수행하는가"를 시연한다. Claude Code가 [`CLAUDE.md`](CLAUDE.md)의 규약을
따라 스킬과 스크립트를 순서대로 호출하고, 각 단계는 표준 봉투 `{result, provenance, verification}`
를 남기며 검증 게이트를 통과해야 다음 단계로 진행한다.

- **표적**: PDE5A (phosphodiesterase 5A, UniProt **O76074**, ChEMBL **CHEMBL1827**)
- **기전**: cGMP 분해 억제 → NO-cGMP 경로 평활근 이완·혈관확장
- **사례 약물**: sildenafil (CHEMBL192) · tadalafil (CHEMBL779) · vardenafil · avanafil

> **주의**: 이 하네스는 알려진 화학공간의 후보를 선별·정리하는 **시연**이며 신약 발견이 아니다.
> 모든 **결과는 가설**이고, 실검증(합성·효소 어세이)은 별도다.

---

## 워크플로

```
                ┌─────────────┐
                │    PLAN      │  plan-research → outputs/plan.json
                │ (JSON 계획)  │  objective/questions/tools/criteria/rules
                └──────┬──────┘
                       ▼
   ┌───────────── EXECUTE (단계별, 검증 게이트) ─────────────┐
   │                                                          │
   │  target-lookup  →  chembl-actives  →  mol-properties  →  selectivity-check
   │  (PDE5A 실재)      (CHEMBL1827)       (QED/SA/Lipinski)   (PDE5 vs PDE6, 정직)
   │       │                 │                   │                   │
   │   [VERIFY]          [VERIFY]            [VERIFY]            [VERIFY]
   │   passed?→next      passed?→next        passed?→next        passed?→next
   └──────────────────────────────┬───────────────────────────────┘
                                   ▼
                          ┌─────────────────┐
                          │     REPORT       │  report-writer
                          │  IMRAD + 검증로그 │  → outputs/report_pde5.md
                          │  + 한계 + 참고문헌 │  (근거 없는 수치 0)
                          └─────────────────┘

     (옵션) mol_utils.py: RCSB 구조 확보  ─┐
     (옵션) docking.py:   smina 도킹        ─┴─→ 있을 때만, 없으면 graceful 스킵
```

각 단계는 `verification.passed=false` 이면 **다음 단계 진행 금지**. 최대 2회 재시도 후에도
실패하면 그 항목을 "미확인/플래그"로 표시하고 보고서 한계 섹션에 명시한다.

---

## 스킬 (`.claude/skills/<name>/SKILL.md`)

| 스킬 | 역할 | 호출 스크립트 |
|------|------|---------------|
| `plan-research` | 연구 계획을 구조화 JSON으로 작성(실행 전 필수) | — (JSON 작성) |
| `target-lookup` | PDE5A 실재·기능 확인 (UniProt O76074) | `scripts/target_lookup.py` |
| `chembl-actives` | CHEMBL1827 활성물질 조회 (실 ID·SMILES) | `scripts/chembl_actives.py` |
| `mol-properties` | 물성·QED·SA·Lipinski 게이트 | `scripts/mol_properties.py` |
| `selectivity-check` | PDE5 vs PDE6 선택성 (정직한 한계) | `scripts/selectivity.py` |
| `report-writer` | IMRAD 보고서 종합 → `outputs/report_pde5.md` | (봉투 종합 + `verify.py`) |

## 스크립트 (`scripts/`)

| 스크립트 | 설명 | 반환 |
|----------|------|------|
| `verify.py` | 검증 공용 유틸(make_result/gate/valid_smiles/numbers_backed) | — |
| `target_lookup.py` | UniProt O76074 조회 | 표준 봉투 |
| `chembl_actives.py` | ChEMBL CHEMBL1827 활성물질 | 표준 봉투 |
| `mol_properties.py` | RDKit 물성/약물성 | 표준 봉투 |
| `selectivity.py` | PDE5/PDE6 정성 선택성 | 표준 봉투 |
| `mol_utils.py` (옵션) | RCSB PDB 다운로드 + 구조 정보 | 표준 봉투 |
| `docking.py` (옵션) | smina 도킹(없으면 graceful 스킵) | 표준 봉투 |

모든 과학 스크립트는 `{result, provenance, verification}` 표준 봉투(JSON)를 반환한다.

## 훅 (`.claude/hooks/`, `.claude/settings.json`)

| 훅 | 이벤트 | 동작 |
|----|--------|------|
| `verify_provenance.py` | PostToolUse(Bash) | 출력에 provenance/verification 필드 존재 점검 → 없으면 경고(비차단) |
| `no_fabrication_guard.py` | PostToolUse(Bash) | 근거 없는 수치 표현·가짜 ID 패턴 리마인더(비차단) |

훅은 데모용 **경고형**(차단하지 않음)이며, `additionalContext`로 에이전트에 리마인더를 주입한다.

---

## 강사 데모 (B 방식) 실행법

1. Claude Code에서 이 폴더를 연다:
   ```bash
   cd /home/hjpark/2026_aidrugdiscovery/Day06_LLM_Agent/pde5_harness
   claude
   ```
2. 자연어로 지시한다:
   > **"PDE5 저해제 후보 조사해 보고서 써줘"**
3. 에이전트가 `CLAUDE.md` 규약대로 PLAN → 단계별 EXECUTE(검증 게이트) → REPORT를 수행하고
   `outputs/report_pde5.md` 를 생성한다. 각 Bash 호출 후 훅이 provenance/무-날조를 점검한다.

### 스크립트 직접 실행 (노트북 없이도 동작)

```bash
python scripts/target_lookup.py
python scripts/chembl_actives.py | python scripts/mol_properties.py --stdin
python scripts/selectivity.py
python scripts/mol_utils.py 1UDT      # (옵션) PDE5A 구조 확보
```

---

## 오프라인 / 네트워크·API 안내

- **네트워크(UniProt/ChEMBL) 불가 시**: `target_lookup.py`·`chembl_actives.py` 는 공개된 참조값
  기반 **오프라인 데모 모드**로 폴백하되 `"OFFLINE DEMO ... 실데이터 아님"` 을 반드시 표기한다.
- **RDKit 미설치 시**: `mol_properties.py` 는 수치를 만들지 않고 정직하게 실패를 보고한다
  (`pip install rdkit` 안내). SA 점수는 RDKit contrib(sascorer) 없으면 `null`.
- **smina 미설치 시**: `docking.py` 는 도킹을 실행하지 않고 "도킹 미실행(구조·smina 필요)"을
  표준 봉투로 반환한다(스코어 날조 금지).
- 모든 스크립트는 import/네트워크 실패에서 **죽지 않고** 무-날조 안내로 graceful 폴백한다.

---

## 하드 규칙 (무-날조)

- **수치는 도구 실계산값만.** LLM이 IC50·물성·ChEMBL ID·PMID 등 수치·식별자를 지어내는 것 금지.
- **모든 사실 주장에 출처**(ChEMBL ID / UniProt / PMID / DOI). 없으면 "확인 필요".
- **단계별 검증 게이트** 통과 후에만 다음 단계.
- **불확실은 불확실로**: 도킹 스코어·LLM 추론은 "가설". 결과=가설, 신약 발견 아님.

주요 서지: Boolell M et al. *Int J Impot Res.* 1996;8(2):47-52 · Ghofrani HA et al.
*Nat Rev Drug Discov.* 2006;5(8):689-702 (doi:10.1038/nrd2030).

---

## Docker

```bash
docker build -t pde5-harness .
docker run --rm -it -v "$PWD/outputs:/workspace/outputs" pde5-harness \
    python scripts/target_lookup.py
```

`requirements.txt` 로 로컬 설치도 가능:
```bash
pip install -r requirements.txt
```

---

## License

MIT
