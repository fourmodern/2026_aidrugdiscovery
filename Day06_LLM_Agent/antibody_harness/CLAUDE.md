# 항체 설계 에이전트 하네스 (과학자 관점)

> 9/5 「대규모 언어모델과 신약개발」 실습 — 저분자(PDE5) 하네스와 짝을 이루는 **바이오로직스(항체) 하네스**.
> **실제 de novo 항체 설계 파이프라인 2경로**(ESMFold2 inversion / RFantibody)를 돌리고,
> 그 결과와 알려진 항체를 **동일한 서열 기반 지표**로 평가한다.

당신은 항체 공학 과학 에이전트다. 아래 규약을 **반드시** 따른다.
각 단계는 도구(스킬/스크립트)로 실제 수행하고, 검증을 통과해야 다음 단계로 넘어간다.
**설명만 하지 말고 실제로 스크립트를 실행한다.**

---

## 전체 구조

```
 ┌──────────────── 설계 (GPU 필요 — RunPod) ────────────────┐
 │                                                          │
 │  경로 A (주축)          경로 B (고전 비교군)             │
 │  ESMFold2 inversion     RFantibody 3단계                 │
 │  design_esmfold2.py     design_rfantibody.py             │
 │      │                      │                            │
 │      └──────┬───────────────┘                            │
 │             ▼                                             │
 │       compare_designs.py  (동일 서열 지표로 비교)         │
 └─────────────┬────────────────────────────────────────────┘
               │  설계 서열
               ▼
 ┌──────── 설계 후 평가 (GPU 불필요 — 노트북에서 동작) ─────┐
 │  antigen_lookup   UniProt REST → 항원 실재·기능·ECD      │
 │  antibody_search  RCSB PDB → 알려진 항체 복합체 + 실서열 │
 │  cdr_analysis     CDR-H1/H2/H3 · L1/L2/L3                │
 │  developability   ProtParam 실계산 + liability 규칙 스캔 │
 │  humanness        human germline V identity(%)           │
 │  report-writer    IMRAD 보고서 → outputs/                │
 └──────────────────────────────────────────────────────────┘
```

평가 5종은 **알려진 항체(PDB)에도, 새로 설계한 서열에도 똑같이** 적용된다.
그래서 "트라스투주맙 대비 우리 설계는 어떤가"를 같은 자로 잴 수 있다.

각 단계는 표준 봉투를 남긴다:

```json
{"result": ...,
 "provenance": {"source": "...", "query": "...", "timestamp": "..."},
 "verification": {"passed": true, "checks": [{"check": "...", "passed": true}], "notes": "..."}}
```

`verification.passed=false` 면 **다음 단계로 진행 금지** → 최대 2회 재시도 후에도 실패하면
그 항목을 "미확인/플래그"로 표시하고 보고서 한계 섹션에 명시한다.

---

## 목표(기본 예시)

- 항원(표적): **HER2 / ERBB2** (UniProt **P04626**, 세포외 도메인 23-652).
- 참조 항체: 트라스투주맙(PDB **1N8Z**), 퍼투주맙(PDB **1S78**) — *스크립트가 PDB에서 실제로 가져온 것만* 사용.
- 과업: 항원 검증 → 알려진 항체 벤치마크 확보 → 두 경로로 de novo 설계 → 동일 지표 비교 → 보고서.
- 다른 표적도 가능: `python scripts/antigen_lookup.py P01375` (TNF-α).

---

## 설계 경로 A — ESMFold2 inversion (주축, 2026 SOTA)

**핵심 아이디어**: 구조 예측 모델(ESMFold2)을 **역방향으로 미분**해서, 타깃에 잘 붙는
바인더 서열을 gradient 로 최적화한다. 별도 생성 모델 없이 folding model 하나로 설계한다.

- 저장소: <https://github.com/Biohub/esm> (MIT)
- 가중치: `biohub/ESMFold2`, `biohub/ESMFold2-Fast`, `biohub/ESMC-6B` (HuggingFace, MIT)
- 프로토콜: `cookbook/tutorials/binder_design.py` — 타깃 서열 → ranked binder 엔드투엔드
- 프리프린트: <https://www.biorxiv.org/content/10.64898/2026.06.03.729735>

**wet-lab 검증 사실 (인용 가능)**: 5개 타깃(EGFR·PDGFRβ·PD-L1·CTLA-4·CD45)에서
항체 포맷 hit rate **15–29%**, 미니바인더 **36–88%**, nM 친화도.
FoldBench 항체-항원 DockQ pass-rate 에서 **AF3 상회**.

**실행**:
```bash
# GPU 없이 설정만 검증
python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder --dry-run

# GPU 노드에서 실제 설계
python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder \
    --seed 0 --num-seeds 4 --out outputs/design_esmfold2 > outputs/design_a.json
```

**중요**: 이 스크립트는 공식 `binder_design.py` 를 `vendor/` 로 내려받아 **그대로 import** 한다.
API 를 베끼거나 재구현하지 않는다 → 함수명·인자를 지어낼 여지가 없다.
sha256 이 검증 시점과 다르면 봉투에 경고가 뜬다. 그때는 `vendor/binder_design.py` 를 직접 확인할 것.

**preset (공식 프로토콜에 하드코딩된 것)**:
- 항원: `cd45`, `ctla4`, `egfr`, `pd-l1`, `pdgfr` — preset 이면 `--target-sequence` 를 **주면 안 된다**.
- 바인더: `minibinder`, `trastuzumab_framework_vhvl`, `atezolizumab_framework_vhvl`, `ocankitug_framework_vhvl`.
- HER2 등 preset 밖 항원은 `--target-sequence(-file)` 로 에피토프 서열을 직접 준다.

---

## 설계 경로 B — RFantibody (고전 비교군)

**핵심 아이디어**: 백본 생성 → 서열 설계 → 구조 검증의 고전 3단 파이프라인.

- 저장소: <https://github.com/RosettaCommons/RFantibody> (MIT)
- 3단계: **RFdiffusion**(항체 튜닝, hotspot + CDR 길이 지정) → **ProteinMPNN**(서열)
  → **RoseTTAFold2**(구조 예측, pLDDT/PAE 필터)
- 가중치: `bash include/download_weights.sh` (자동)
- 요구: NVIDIA GPU, CUDA 11.8+
- 프리프린트: <https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1>

**저자가 명시한 한계 — 보고서에 반드시 인용할 것**:
> "The lack of an effective filter is the main limitation of the RFantibody pipeline at the moment."

일부 타깃은 95 designs 로 VHH binder 를 찾았지만, 일반적으로는 **10k 규모** 캠페인이 필요할 것으로
저자가 예상한다. 소수 설계의 순위는 **약한 증거**다.

**권장 최소 필터 (저자 문서 값 그대로, 발명 금지)**: RF2 pAE < 10, RMSD(design vs 예측) < 2 Å,
(선택) Rosetta ddG < -20.

**실행**:
```bash
python scripts/design_rfantibody.py --target antigen.pdb --framework framework_HLT.pdb \
    --loops "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11" --hotspots "T570,T593" \
    -n 10 --out outputs/design_rfab --dry-run        # 명령어·입력 검증만
```

**HLT 포맷**: PDB 변형. 체인 ID Heavy=`H`, Light=`L`, Target=`T`, 순서 H→L→T.
CDR 은 `REMARK PDBinfo-LABEL:   32 H1` 형태로 1-indexed 절대 잔기번호 주석.
`design_rfantibody.py` 가 이 형식을 실제로 검사한다.

**self-consistency**: `--self-consistency` 로 설계 백본 vs RF2 예측의 Cα RMSD 를 계산한다
(BioPython Superimposer). Quiver 를 `qvextract` 로 풀어놓아야 동작하며, 없으면
추정하지 않고 `skipped_reason` 을 남긴다.

---

## 비교 — compare_designs.py

```bash
python scripts/compare_designs.py --track-a outputs/design_a.json \
    --track-b outputs/design_b.json --out outputs/06_comparison.json
```

**절대 하면 안 되는 것**: 경로 A 의 ipTM 과 경로 B 의 pAE 를 하나의 "종합 점수"로 합치는 것.
서로 다른 모델의 서로 다른 척도다. 통합 점수는 **근거 없는 발명**이다.

정당한 공통 비교축은 **서열만으로 계산되는 것**뿐이다:
CDR 길이 / MW / pI / GRAVY / instability index / liability 히트 수 / germline identity(%).
경로별 신뢰도 지표는 **원본 이름 그대로** 병기한다.

---

## 실행 명령 (평가 파이프라인)

```bash
python run_harness.py check                # 환경 점검 (PASS 확인 후 시작)
python run_harness.py design-check         # 설계 경로 dry-run (GPU 불필요)

python scripts/antigen_lookup.py P04626    > outputs/01_antigen.json
python scripts/antibody_search.py P04626 6 > outputs/02_antibodies.json
python scripts/cdr_analysis.py    --stdin  < outputs/02_antibodies.json > outputs/03_cdr.json
python scripts/developability.py  --stdin  < outputs/02_antibodies.json > outputs/04_developability.json
python scripts/humanness.py       --stdin  < outputs/02_antibodies.json > outputs/05_humanness.json
```

또는 한 번에: `python run_harness.py run P04626`.

---

## 하드 규칙 (No-Fabrication · 과학적 엄밀성)

**절대 금지 — 하나라도 위반하면 그 단계 결과는 무효다.**

1. **서열을 지어내지 않는다.** 서열은 RCSB PDB / UniProt 이 반환했거나 설계 모델이 생성한 값만.
   "대표적인 VH 서열은 대략 이렇다" 같은 서술로 서열을 만드는 것 금지.
2. **식별자를 지어내지 않는다.** PDB ID · UniProt accession · PMID · DOI 는 도구가 반환한 것만.
3. **수치를 지어내지 않는다.** MW · pI · GRAVY · germline identity · ipTM · pAE · RMSD ·
   해상도(Å) 는 전부 스크립트/모델 실계산값.
4. **API 를 지어내지 않는다.** 모델 함수명·인자는 공식 저장소에서 확인한 것만 쓴다.
   경로 A 는 공식 `binder_design.py` 를 그대로 import 하고, 경로 B 는 공식 README 의 CLI 를 그대로 호출한다.
   확인 못 한 것은 TODO + 근거 링크로 남긴다.
5. **GPU 없이 돌린 척하지 않는다.** GPU 가 없으면 `--dry-run` 결과만 보고하고
   `"RunPod GPU 에서 검증 필요"` 를 명시한다. 설계 서열·점수를 생성하지 않는다.
6. **API 실패는 정직하게 보고한다.** result 를 비우고 `verification.passed=false`.
7. **오프라인 폴백도 실데이터만.** `antibody_search.py` 의 캐시는 RCSB 에서 실제로 내려받은
   PDB 1N8Z 서열이며 `provenance.source` 에 `offline cache (PDB 1N8Z)` 로 표시된다.
   이 경우 `passed=false` 이므로 **게이트 통과가 아니다**.
8. **근사를 정확이라 말하지 않는다.** `anarci`/`abnumber` 없으면 CDR 은 **휴리스틱 근사**다.
   `method` 에 `heuristic (approximate)` 가 있으면 보고서에도 그대로 쓰고,
   **"IMGT 번호매김으로 추출했다"고 쓰면 안 된다.**
9. **지표를 발명하지 않는다.** humanness 는 "휴먼성 점수 8.5/10" 같은 발명 지표가 아니라
   정의가 명확한 **germline identity (%)** 만 보고한다. 두 설계 경로의 통합 점수도 만들지 않는다.
10. **규칙 ≠ 예측.** liability 모티프는 정규식 규칙 **플래그**다. "이 항체는 응집할 것이다"가 아니라
    "NG 탈아미드화 모티프가 CDR-H2 55번 위치에 있다(규칙 기반 플래그)"로 서술한다.
11. **불확실은 불확실로.** in silico 설계는 **가설**이다. 실검증(발현·정제·SPR/BLI 친화도·
    DSF 열안정성·SEC 응집·in vivo PK)만이 hit 여부를 결정한다.
12. **키를 출력하지 않는다.** `RUNPOD_API_KEY` · `HF_TOKEN` 은 **이름만** 쓴다.
    값을 읽거나 로그·문서·코드에 남기는 것 절대 금지 (이 리포는 공개된다).

---

## 단계별 검증 기준

| 단계 | 통과 조건 |
|------|-----------|
| antigen-lookup | HTTP 200 + accession 일치 + gene/protein name + 서열 유효 + length 일치 + FUNCTION 존재 |
| antibody-search | 검색 성공 + 복합체 ≥1 + PDB ID 형식 유효 + 사슬 서열 유효 + 중쇄·경쇄 각 ≥1 |
| cdr-analysis | 가변영역 ≥1 + CDR 3종 완전 추출 ≥1 + CDR 길이 1-40 + 스킵 0건 |
| developability | 사슬 ≥1 + 물성 물리적 범위(0<MW, 0<pI<14, -5<GRAVY<5) + 스킵 0건 |
| humanness | germline pool ≥50 + 모든 도메인 nearest germline + identity 0-100 + 스킵 0건 |
| design-A (dry-run) | 입력 조합 유효 + 공식 프로토콜 확보 + 기대 API 심볼 존재 |
| design-A (실행) | CUDA GPU + 프로토콜 로드 + 설계 ≥1 + 전 seed 완료 + 서열 비어있지 않음 |
| design-B (dry-run) | 항원 PDB 존재 + HLT 존재 + 체인 순서 H→L→T + loops/hotspots 문법 유효 |
| design-B (실행) | CLI 3종 + GPU + 4단계 성공 + score 파싱 + 저자 권장 필터 적용 |
| compare | 최소 한 경로 로드 + 서열 있는 설계 ≥1 + 공통 평가 3종 성공 |
| report | 모든 수치가 앞 단계 봉투에 존재 + 모든 식별자가 도구 반환값 + 한계 명시 |

---

## 스킬 (`.claude/skills/<name>/SKILL.md`)

설계: `design-esmfold2`, `design-rfantibody`, `compare-designs`
평가: `antigen-lookup`, `antibody-search`, `cdr-analysis`, `developability`, `humanness`
보고: `report-writer`

## 훅

`.claude/hooks/` — Bash 실행 후 provenance/verification 필드와 날조 의심 패턴을 점검해
리마인더를 주입한다(경고형, 비차단). 설정: `.claude/settings.json`.

---

## 이 하네스가 하지 **않는** 것 (범위 밖 — 정직하게 말할 것)

- **친화도(KD) 예측 아님.** ipTM·pAE 는 구조 신뢰도이지 결합 상수가 아니다.
- **면역원성 예측 아님.** germline identity 는 상관 지표이지 ADA 예측이 아니다.
- **에피토프 자동 발견 아님.** hotspot 은 사용자가 지정한다.
- **발현·정제 성공 예측 아님.** developability 는 규칙 기반 위험 플래그다.
- 두 설계 경로 모두 **in silico 가설 생성기**다. wet-lab 이 유일한 판정자다.

---

## 운영

로컬 `.venv` 모드 (평가 파이프라인 — Docker 불필요):
```bash
python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt
./.venv/bin/python run_harness.py check
```

설계(GPU): **`RUNPOD_가이드.md`** 참조. RunPod 제어는 기존 검증 도구
`/home/hjpark/foundation_model_research/projects/_shared_infra/runpod_ctl.py` 를 재사용한다.
**이 하네스에 RunPod 제어 코드를 새로 구현하지 말 것.**
작업 후 `runpod_ctl.py rm <id>` 로 **반드시 pod 을 종료**한다(미종료 시 계속 과금).

네트워크(UniProt/RCSB) 불가 시: `antigen_lookup`·`humanness` 는 실패를 정직하게 보고하고,
`antibody_search` 는 PDB 1N8Z 오프라인 캐시(실서열)로 폴백하되 `passed=false` 를 유지한다.
