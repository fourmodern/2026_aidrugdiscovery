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
| `dataset-builder` | 역가 × 골격유사도 직교 데이터셋 (9칸 격자) | `scripts/build_controlled_dataset.py` |
| `structure-fetch` | RCSB 결정 구조 + 공결정 리간드 확보 | `scripts/dock_controlled.py` (내장) |
| `docking-validation` | 도킹 + 재도킹 대조 C1(샘플링)/C2(채점) 분리 | `scripts/dock_controlled.py` |
| `exhaustiveness-sweep` | 탐색 깊이를 바꿔 자세 실패와 채점 실패를 가름 | `scripts/exhaustiveness_sweep.py` |
| `confound-control` | 골격 유사성 통제 — 편상관·구간내 상관 필수 | `scripts/analyze_controlled.py` |
| `score-regression` | 항 회귀 — 표준화 계수·VIF·Q²·라벨섞기 | `scripts/terms_controlled.py` |
| `statistics-validation` | 직접 구현한 통계 함수를 scipy 와 대조 | `scripts/test_statistics.py` |
| `binding-mode` | PyMOL 결합 양상 + 좌표 기반 접촉 집계 | `scripts/render_binding_mode.py` |
| `figure-builder` | 논문급 그림 세트 (PNG + SVG, 색각 안전) | `scripts/figures_controlled.py` |
| `report-writer` | IMRAD 보고서 종합 → `sample_run/report/` | `scripts/report_controlled.py` |
| `report-docs` | 본문 → docx · pptx (신선도 검사 포함) | `scripts/make_docs.py` |
| `pi-review` / `devils-advocate` / `galley-proof` | 외부 critic 리뷰 게이트 | (Agent tool) |

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
| `build_controlled_dataset.py` | ChEMBL 전수 → 중앙값 집계 → 9칸 격자 추출 | 표준 봉투 |
| `dock_controlled.py` | 병렬 도킹 + 재도킹 대조 (C1/C2) | 표준 봉투 |
| `exhaustiveness_sweep.py` | 깊이 × 시드 스윕, 자세/채점 진단 | 표준 봉투 |
| `analyze_controlled.py` | 편상관·구간내 상관·이중회귀·ROC | 표준 봉투 |
| `terms_controlled.py` | 항 추출 + 표준화 OLS + VIF + 순열 중요도 | 표준 봉투 |
| `test_statistics.py` | 통계 함수 scipy 대조 (동점 주입) | 표준 봉투 |
| `contact_concordance.py` | 접촉 잔기 — 두 자세 규칙 비교 | 표준 봉투 |
| `figures_controlled.py` | 논문급 그림 12종 (PNG + SVG) | 파일 |
| `report_controlled.py` | v3.0 보고서 생성 (수치 전량 주입) | 파일 |
| `make_docs.py` | 마크다운 → docx · pptx + 신선도 검사 | 파일 |

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

## 연구 산출물 (`sample_run/`)

이 하네스가 실제로 수행한 연구 결과가 들어 있다. **세 판본이 모두 보존되어 있으며, 그 자체가
이 하네스의 가장 중요한 교육 자료다** — 자동 게이트를 전부 통과한 결론이 두 번 틀렸다.

| 판본 | 설계 | 결론 | 무엇이 무너뜨렸나 |
|------|------|------|-------------------|
| v1.0 | n=30, 역가로만 층화 | "도킹은 역가를 예측하지 못한다" | 선별 지표를 하나도 계산하지 않았다 |
| v2.0 | 같은 데이터 + 선별 지표 | "선별은 되고 정량은 안 된다" | 강한 화합물이 전부 공결정 리간드 유사체였다 |
| **v3.0** | n=163, 역가 × 골격 직교 설계 | `report_controlled.md` 참조 | — (현행) |

```
sample_run/
├── dataset_controlled.json     역가 3구간 × Tanimoto 3구간 9칸 격자 (v3.0)
├── docking_controlled.json     도킹 결과 + 재도킹 대조 (C1/C2 분리)
├── exhaustiveness_sweep.json   탐색 깊이 스윕 — 자세 실패와 채점 실패를 가른다
├── analysis_controlled.json    편상관 · 구간내 상관 · 선별 지표
├── terms_controlled.json       스코어 항 회귀 (표준화 계수 · VIF · 순열 중요도)
├── statistics_validation.json  직접 구현한 통계 함수의 scipy 대조 결과
├── contact_concordance.json    접촉 잔기 — 참조 선택 자세 vs 점수 1위 자세
├── structures/work_controlled/ 도킹 자세 SDF 전수 (재현 검증용)
└── report/
    ├── report_controlled.md    v3.0 본문 (현행)
    ├── report_pde5.md          v2.0 본문 (대체됨 — 감사 추적용 보존)
    ├── figures_controlled/     v3.0 그림 (PNG + SVG)
    └── docs/                   docx · pptx
```

모든 산출 파일은 표준 봉투 `{result, provenance, verification}` 형식이다. 보고서의 모든
수치는 이 파일들에서 스크립트가 주입한 값이며, 사람이 타이핑한 것은 서술문과 절 제목뿐이다.

### 재현

```bash
python scripts/build_controlled_dataset.py --per-cell 20      # ChEMBL 전수 → 9칸 추출
python scripts/dock_controlled.py --exhaustiveness 64 --workers 24
python scripts/exhaustiveness_sweep.py                        # 탐색 깊이 진단
python scripts/analyze_controlled.py                          # 교란 통제 분석
python scripts/terms_controlled.py                            # 항 기여도
python scripts/test_statistics.py                             # 통계 함수 검증
python scripts/figures_controlled.py                          # 그림 전체
python scripts/report_controlled.py                           # 보고서
```

---

## 하네스 없이 재현하기 — `prompts/claude_for_science.md`

로컬 파일(`CLAUDE.md`·스킬·스크립트)이 없는 환경(Claude for Science 등)에서 같은 연구를
수행하기 위한 프롬프트다. **하네스가 파일로 강제하던 계약을 프롬프트 본문이 대신 진다.**

- **A. 마스터 프롬프트** (약 2,000단어) — 그대로 붙여넣는다
- **B. 각 조항이 무엇을 막는가** — 21개 조항을 실제 사건과 짝지은 표
- **C. 다른 표적으로 바꿀 때** — 치환할 자리, 이 연구에 맞춰진 값, 바꾸면 안 되는 다섯

조항은 전부 이 연구가 한 번씩 틀린 뒤에 추가된 것이다. 프롬프트 자체도 외부 리뷰를
받았고, 그 리뷰가 **프롬프트가 이 연구의 철회된 인과 주장을 되풀이하고 있다**는 것과
**틀린 결론 횟수를 세 번에서 두 번으로 축소해 적었다**는 것을 잡아 고쳤다.

파일이 없어 사라지는 보호(게이트의 차단 동작, 별도 리뷰어 에이전트, 음성 대조)는 B 절
말미에 대체 가능 여부와 함께 명시했다. **별도 리뷰어와 음성 대조는 대체 불가**다.

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
