# 항체 설계 에이전트 하네스

[![Claude Code](https://img.shields.io/badge/Claude%20Code-harness-8A2BE2)](https://docs.claude.com/en/docs/claude-code)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![BioPython](https://img.shields.io/badge/BioPython-ProtParam%20%2B%20Align-green.svg)](https://biopython.org/)
[![ESMFold2](https://img.shields.io/badge/design-ESMFold2%20inversion-orange.svg)](https://github.com/Biohub/esm)
[![RFantibody](https://img.shields.io/badge/baseline-RFantibody-lightgrey.svg)](https://github.com/RosettaCommons/RFantibody)
[![No-Fabrication](https://img.shields.io/badge/policy-no--fabrication-red.svg)](#하드-규칙-무-날조)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)

> 9/5 「대규모 언어모델과 신약개발」 실습 — 저분자(PDE5) 하네스와 짝을 이루는 **바이오로직스(항체) 하네스**.
> **실제 de novo 항체 설계 2경로**를 돌리고, 그 결과와 알려진 항체를 **동일한 서열 기반 지표**로 평가한다.
> 서열·수치·식별자·API 함수명은 **도구 실계산값 또는 공식 문서 확인값만** 사용한다.

---

## 개요

Claude Code 같은 자율 에이전트가 이 폴더를 열고 [`CLAUDE.md`](CLAUDE.md) 규약에 따라
스킬과 스크립트를 순서대로 호출한다. 각 단계는 표준 봉투
`{result, provenance, verification}` 를 남기며, 검증 게이트를 통과해야 다음 단계로 간다.

```
 ┌──────────────── 설계 (GPU 필요 — RunPod) ────────────────┐
 │  경로 A (주축)              경로 B (고전 비교군)          │
 │  ESMFold2 inversion         RFantibody 3단계              │
 │  design_esmfold2.py         design_rfantibody.py          │
 │      │                          │                         │
 │      └───────────┬──────────────┘                         │
 │                  ▼                                         │
 │           compare_designs.py                               │
 └──────────────────┬─────────────────────────────────────────┘
                    │ 설계 서열
                    ▼
 ┌──── 설계 후 평가 (GPU 불필요 — 노트북에서 즉시 동작) ────┐
 │  antigen_lookup → antibody_search → cdr_analysis →       │
 │  developability → humanness → report                     │
 └──────────────────────────────────────────────────────────┘
```

평가 5종은 **알려진 항체(PDB)에도, 새로 설계한 서열에도 똑같이** 적용된다.
그래서 "트라스투주맙 대비 우리 설계는 어떤가"를 같은 자로 잴 수 있다.

- **기본 항원**: HER2 / ERBB2 (UniProt **P04626**, 세포외 도메인 23-652)
- **참조 항체**: 트라스투주맙(PDB **1N8Z**), 퍼투주맙(PDB **1S78**) — 스크립트가 PDB 에서 실제로 가져온 것

> **주의**: 이 하네스의 in silico 결과는 **가설**이다. 실검증(발현·정제·SPR/BLI 친화도·
> DSF 열안정성·SEC 응집·in vivo PK)만이 hit 여부를 결정한다.

---

## 빠른 시작

```bash
cd /home/hjpark/lecture_drug1/0905_agent/antibody_harness

python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt

./.venv/bin/python run_harness.py check          # 환경 점검 (PASS 확인)
./.venv/bin/python run_harness.py design-check   # 설계 경로 dry-run (GPU 불필요)
./.venv/bin/python run_harness.py run P04626     # 평가 파이프라인 전체 실행
```

또는 Claude Code 로 폴더를 열고 자연어로:

> **"CLAUDE.md 읽고 HER2 항체 하네스 자율 실행해줘"**

---

## 설계 경로

### 경로 A — ESMFold2 inversion (주축, 2026 SOTA)

구조 예측 모델을 **역방향으로 미분**해 바인더 서열을 gradient 로 최적화한다.

| 항목 | 값 |
|------|-----|
| 저장소 | <https://github.com/Biohub/esm> (MIT) |
| 프로토콜 | `cookbook/tutorials/binder_design.py` (타깃 서열 → ranked binder 엔드투엔드) |
| 가중치 | `biohub/ESMFold2`, `biohub/ESMFold2-Fast`, `biohub/ESMC-6B` (HuggingFace, MIT) |
| 설치 | `pip install 'esm@git+https://github.com/Biohub/esm.git@main' abnumber modal` |
| 프리프린트 | <https://www.biorxiv.org/content/10.64898/2026.06.03.729735> |
| GPU | A100 40GB 이상 (80GB 권장) |

**wet-lab 검증**: 5개 타깃(EGFR·PDGFRβ·PD-L1·CTLA-4·CD45)에서 항체 포맷 hit rate **15–29%**,
미니바인더 **36–88%**, nM 친화도. FoldBench 항체-항원 DockQ pass-rate 에서 **AF3 상회**.

`design_esmfold2.py` 는 공식 `binder_design.py` 를 `vendor/` 로 내려받아 **그대로 import** 한다
(API 재구현 없음 → 함수명·인자를 지어낼 여지 없음). sha256 을 검증 시점 값과 대조해
upstream 변경 시 경고한다.

```bash
python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder --dry-run
```

### 경로 B — RFantibody (고전 비교군)

RFdiffusion(백본) → ProteinMPNN(서열) → RoseTTAFold2(구조 검증) 3단계.

| 항목 | 값 |
|------|-----|
| 저장소 | <https://github.com/RosettaCommons/RFantibody> (MIT) |
| 가중치 | `bash include/download_weights.sh` |
| 설치 | `uv sync` / Docker / Apptainer |
| 요구 | NVIDIA GPU, CUDA 11.8+, Ubuntu 22.04 권장 |
| 프리프린트 | <https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1> |

**저자가 명시한 한계 (보고서에 반드시 인용)**:
> "The lack of an effective filter is the main limitation of the RFantibody pipeline at the moment."

일부 타깃은 95 designs 로 VHH binder 를 찾았지만, 일반적으로 **10k 규모** 캠페인이 필요할 것으로
저자가 예상한다. → **소수 설계의 순위는 약한 증거다.**

권장 최소 필터(저자 문서 값): `RF2 pAE < 10`, `RMSD < 2 Å`, (선택) `Rosetta ddG < -20`.
`--self-consistency` 로 설계 백본 vs RF2 예측의 Cα RMSD 를 추가 계산한다.

```bash
python scripts/design_rfantibody.py --target antigen.pdb --framework framework_HLT.pdb \
    --loops "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11" --hotspots "T570,T593" \
    -n 10 --out outputs/design_rfab --dry-run
```

### 비교

```bash
python scripts/compare_designs.py --track-a outputs/design_a.json \
    --track-b outputs/design_b.json --out outputs/06_comparison.json
```

**두 경로의 신뢰도 지표를 하나의 종합 점수로 합치지 않는다.** ipTM 과 pAE 는 서로 다른 모델의
서로 다른 척도다. 공통 비교축은 **서열만으로 계산되는 지표**뿐이다.

---

## 스킬 (`.claude/skills/<name>/SKILL.md`)

| 스킬 | 역할 | 호출 스크립트 | GPU |
|------|------|---------------|-----|
| `design-esmfold2` | 경로 A — ESMFold2 inversion 설계 | `scripts/design_esmfold2.py` | 필요 |
| `design-rfantibody` | 경로 B — RFantibody 3단계 | `scripts/design_rfantibody.py` | 필요 |
| `compare-designs` | 두 경로 동일 지표 비교 | `scripts/compare_designs.py` | 불필요 |
| `antigen-lookup` | 항원 실재·기능·ECD 확인 | `scripts/antigen_lookup.py` | 불필요 |
| `antibody-search` | 알려진 항체 복합체 + 실서열 | `scripts/antibody_search.py` | 불필요 |
| `cdr-analysis` | CDR-H1/H2/H3 · L1/L2/L3 추출 | `scripts/cdr_analysis.py` | 불필요 |
| `developability` | ProtParam 물성 + liability 규칙 | `scripts/developability.py` | 불필요 |
| `humanness` | germline identity (%) | `scripts/humanness.py` | 불필요 |
| `report-writer` | IMRAD 보고서 → `outputs/` | (봉투 종합 + `verify.py`) | 불필요 |

## 스크립트 (`scripts/`)

| 스크립트 | 설명 | 반환 |
|----------|------|------|
| `verify.py` | 검증 공용 유틸 (`make_result`/`gate`/`valid_protein_seq`/`numbers_backed`) | — |
| `seq_utils.py` | 서열 파싱·사슬 분류·V-domain·CDR 휴리스틱·scFv 분할 | — |
| `antigen_lookup.py` | UniProt REST 항원 조회 | 표준 봉투 |
| `antibody_search.py` | RCSB Search+Data API 항체 복합체 + 실서열 | 표준 봉투 |
| `cdr_analysis.py` | CDR 추출 (anarci 정식 / 휴리스틱 근사) | 표준 봉투 |
| `developability.py` | BioPython ProtParam + liability 규칙 스캔 | 표준 봉투 |
| `humanness.py` | UniProt germline + BLOSUM62 local identity | 표준 봉투 |
| `design_esmfold2.py` | 경로 A (공식 프로토콜 vendor import) | 표준 봉투 |
| `design_rfantibody.py` | 경로 B (공식 CLI 오케스트레이션 + Cα RMSD) | 표준 봉투 |
| `compare_designs.py` | 두 경로 동일 지표 비교 | 표준 봉투 |

모든 과학 스크립트는 `{result, provenance, verification}` 표준 봉투(JSON)를 반환한다.

## 훅 (`.claude/hooks/`, `.claude/settings.json`)

| 훅 | 이벤트 | 동작 |
|----|--------|------|
| `verify_provenance.py` | PostToolUse(Bash) | provenance/verification 필드 점검, `passed=false`·휴리스틱·dry-run 리마인더 |
| `no_fabrication_guard.py` | PostToolUse(Bash) | 서열 창작·잘못된 PDB ID·"GPU 없는데 설계 결과" 모순 경고 |

경고형(비차단)이며 `additionalContext` 로 에이전트에 리마인더를 주입한다.

---

## GPU 실행 (RunPod)

**[`RUNPOD_가이드.md`](RUNPOD_가이드.md) 를 읽을 것.** 요약:

- RunPod 제어는 기존 검증 도구
  `/home/hjpark/foundation_model_research/projects/_shared_infra/runpod_ctl.py` 를 **재사용**한다.
  이 하네스에 RunPod 제어 코드를 새로 구현하지 않는다.
- 표준 흐름: `deploy` → `wait` → `put` → `ssh`(가중치+설계) → `get` → **`rm`**
- 환경변수는 **이름만** 사용: `RUNPOD_API_KEY`, `HF_TOKEN`. 값은 문서·로그·코드에 남기지 않는다
  (`.env` 는 `.gitignore` 에 포함).
- ⚠️ **작업 후 `rm` 하지 않으면 계속 과금된다.**

---

## 오프라인 / 네트워크 실패 시 동작

| 스크립트 | 네트워크 불가 시 |
|----------|------------------|
| `antigen_lookup.py` | `result: null` + `passed: false`. 값을 만들지 않는다 |
| `antibody_search.py` | PDB **1N8Z 오프라인 캐시(실서열)** 로 폴백. `passed: false` 유지 → 게이트 통과 아님 |
| `humanness.py` | 로컬 캐시가 있으면 사용(출처에 `[local cache]` 표기), 없으면 `result: []` + `passed: false` |
| `developability.py` / `cdr_analysis.py` | 네트워크 불필요 (BioPython 없으면 정직 실패) |
| `design_esmfold2.py` | 프로토콜 다운로드 실패 시 `mode: aborted`. API 를 추측 재구현하지 않는다 |

모든 스크립트는 import/네트워크 실패에서 **죽지 않고** 무-날조 안내로 graceful 폴백한다.

---

## 하드 규칙 (무-날조)

1. **서열을 지어내지 않는다** — UniProt/RCSB 또는 설계 모델 반환값만.
2. **식별자를 지어내지 않는다** — PDB ID / UniProt accession / DOI / PMID.
3. **수치를 지어내지 않는다** — MW·pI·GRAVY·germline identity·ipTM·pAE·RMSD.
4. **API 를 지어내지 않는다** — 경로 A 는 공식 스크립트를 그대로 import, 경로 B 는 공식 README CLI.
5. **GPU 없이 돌린 척하지 않는다** — dry-run 결과만 보고 + "RunPod GPU 에서 검증 필요" 명시.
6. **근사를 정확이라 말하지 않는다** — `heuristic (approximate)` 는 그대로 표기.
7. **지표를 발명하지 않는다** — "휴먼성 점수" 없음, 두 경로 통합 점수 없음.
8. **규칙 ≠ 예측** — liability 는 정규식 플래그이지 예측 모델이 아니다.
9. **키를 출력하지 않는다** — `RUNPOD_API_KEY` / `HF_TOKEN` 은 이름만.

---

## 이 하네스가 하지 않는 것

- 친화도(KD) 예측 — ipTM·pAE 는 구조 신뢰도다
- 면역원성(ADA) 예측 — germline identity 는 상관 지표다
- 에피토프 자동 발견 — hotspot 은 사용자가 지정한다
- 발현·정제 성공 예측 — developability 는 규칙 기반 위험 플래그다

---

## License

MIT

주요 출처: UniProt · RCSB PDB · [Biohub/esm](https://github.com/Biohub/esm) (MIT) ·
[RosettaCommons/RFantibody](https://github.com/RosettaCommons/RFantibody) (MIT) · BioPython
