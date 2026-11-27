# 2026 AI 신약개발 — 대규모 언어모델 & 에이전트 (6일차 실습)

재직자 과정 **6일차: 대규모 언어모델과 신약개발(LLM·에이전트)** 의 실습 자료입니다.
2025 워크숍(`fourmodern/2025_aidrugdiscovery`)에서 **LLM·에이전트 관련 실습만** 선별해 재구성했습니다.

이론(슬라이드) 5교시 + 실습 4H 구성입니다. **2026-11-28 회차부터 메인 실습은 Claude Science** 로 진행하며,
아래 노트북들은 무료 대안·보조 실습으로 남습니다. 모든 노트북은 Google Colab에서 바로 실행됩니다.

---

## ⭐ 메인 실습 — Claude Science 로 PDE5 저해제 선별

| 무엇 | 어디 |
|------|------|
| **⭐ 실습 핸드북 (PDF, 28쪽)** — 설치(Windows 상세)부터 프롬프트 전문까지 | [`1128_실습핸드북.pdf`](Day06_LLM_Agent/1128_실습핸드북.pdf) |
| 실습 교재 (Markdown, 단계별 프롬프트 + '확인할 것') | [`claude_science_pde5_practice.md`](Day06_LLM_Agent/claude_science_pde5_practice.md) |
| 강의 화면용 덱 (PDF, 40쪽) | [`1128_실습_클로드사이언스.pdf`](Day06_LLM_Agent/1128_실습_클로드사이언스.pdf) |
| 도구 사실 근거 (공식 문서 확인 2026-11-27) | [`claude_science_facts_202611.md`](Day06_LLM_Agent/claude_science_facts_202611.md) |

> **수강생은 핸드북 PDF 하나만 있으면 된다.** 0장 준비물 6개 → 1장 설치(Windows 단계별, macOS·Linux 요약)
> → 2장 수업 전 설정(커넥터 2 · 키 1 · 이메일 1 · 검토와 Auto mode) → 3장 여덟 번의 복사·붙여넣기
> → 4장 막히는 곳·더 해볼 것.

## 이론 강의 덱 (5교시, 292장)

이론 덱은 **미러 저장소**에 있다 —
[`fourmodern/2025_aidrugdiscovery` → `20260905/이론덱_1128/`](https://github.com/fourmodern/2025_aidrugdiscovery/tree/main/20260905/%EC%9D%B4%EB%A1%A0%EB%8D%B1_1128)

| 교시 | 장수 | 주제 |
|------|------|------|
| 1교시 | 66 | LLM 원리 · 단백질/오믹스 언어모델 · 분자 언어모델 · 범용과 전용의 경계 |
| 2교시 | 71 | 프롬프트·커스터마이징 4방법 · 환각과 RAG · 검색 대 기억 · 벤치마크 판별 |
| 3교시 | 53 | 에이전트 원리(ReAct·도구·계획·메모리) · 사례 · 폐루프와 한계 |
| 4교시 | 47 | 컨텍스트 엔지니어링 · 서브에이전트 · MCP · 평가 · Agentic 보안 |
| 5교시 | 55 | 자율 과학 · 자율실험실 · world model과 가상세포 · lab in the loop |

> 이 저장소가 아니라 미러에 둔 이유: 이 저장소는 수강생이 Claude Science 로 **가져오는**
> 대상이라 100 MB 제한을 지켜야 한다([`skills/README.md`](Day06_LLM_Agent/skills/README.md)).
> 이론 덱 5종이 약 16 MB 라서 여기 넣으면 제한을 넘어 스킬 가져오기가 막힌다.

동향 정리 → 선행 논문 → 활성·decoy 데이터셋 → 1UDT 도킹 → **점수 분포로 판정** →
보고서(docx·pptx) → **리뷰어 통과까지 반복**.
공개 DB(FDA·ChEMBL·PubChem·PDB·OpenAlex·PubMed)와 DUD-E 벤치마크만 쓰고, 도킹은 CPU 에서 도는
smina 를 쓴다(Windows·macOS·Linux 모두 conda-forge 빌드 있음).

실습에 쓰는 과학 스킬 5종(scientific-writing · scientific-visualization · citation-management ·
peer-review · scientific-critical-thinking)은 [`Day06_LLM_Agent/skills/`](Day06_LLM_Agent/skills)
에 담아 두었다 — 원본 저장소가 약 471 MB 라 Claude Science 의 가져오기 100 MB 제한을 넘기 때문이다.
`Settings > Skills > Add skill > Import from GitHub` 에서 이 저장소를 지정하면 된다.
직접 만드는 경로(Add skill > Chat with Claude)는 핸드북 ⑥-0 장에 있다.

수강생 준비물: **Pro 이상 Claude 플랜 + 앱 설치·로그인**(첫 실행에 약 5 GB 다운로드 — 수업 전에 끝낼 것).
수업 직전 체크 5개 — ChEMBL·PubMed 디렉터리 커넥터 추가 / OpenAlex 무료 키 / 연락 이메일 /
**Auto-review 켜기**(Pro 는 기본 off) / **Auto mode 확인**(꺼져 있으면 코드 셀마다 승인을 눌러야 한다).
둘은 다른 장치다 — Auto-review 는 주장이 실행 기록과 맞는지 사후 검사하고, Auto mode 는 매 단계
승인을 묻지 않는 실행 권한이다.

---

## 실습 목록 (`Day06_LLM_Agent/`)

| 교시 매핑 | 노트북 | 주제 | 설명 | Colab |
|------|------|------|------|-------|
| 1교시 p54–57 · 화학 LM | t042_molt5 | MolT5 | T5 기반 분자↔텍스트 상호 변환(분자 설명 생성 / 텍스트→분자) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/notebooks/t042_molt5.ipynb) |
| 1교시 p48–50 · 단백질 LM (자율 과제) | t110_esm2 | ESM-2 | ESM-2 단백질 언어모델로 펩타이드 결합 최적화(masked LM 스코어링) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/notebooks/t110_esm2_peptide_optimization_tutorial.ipynb) |
| 1교시 p48–50 · 단백질 LM | t111_esm3 | ESM3 | 생성형 멀티모달 단백질 모델 — 서열 설계(inpainting) + 구조 예측 + 신뢰도(pTM/pLDDT) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/notebooks/t111_esm3_protein_design.ipynb) |
| 2교시 p45–48 · RAG 구현·chunking·유사도 | t050_rag | RAG | PDF 문서 기반 로컬 RAG 파이프라인 구축(신약개발 질의응답) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/notebooks/t050_simple-local-rag.ipynb) |
| 3교시 p7–14 · ReAct·도구 (2교시 p67–70 재방문) | nb1_scientific_agent | 직접 만드는 과학 에이전트 | 프레임워크 없이 ReAct 루프·도구 호출·메모리를 손으로 구현 (PubMed·UniProt·ChEMBL·RDKit + 사내 위키 도구). 무-키 경로로 전체 실행 가능 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/notebooks/nb1_scientific_agent.ipynb) |
| 3교시 p28 · AI co-scientist | T004_llm_coscientist | LLM Co-Scientist | LangChain + LLM 기반 PDE5/PAH 표적 탐색 AI 에이전트 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/notebooks/T004_pde5_llm_coscientist.ipynb) |

> 이론 흐름(참고): ① LLM 원리 & 과학 언어모델 → ② 활용 & 환각 통제(RAG) → ③ 에이전트 → ④ 하네스·보안 → ⑤ 자율 과학.
> 쪽수는 이론 덱 PDF 기준이다(미러 저장소 `20260905/이론덱_1128/`). 4교시(하네스·보안)와 5교시(자율 과학)에는 끼우는 노트북이 없다 —
> 대신 nb1 의 도구 정의부를 다시 열어 "도구에 어떤 권한을 주는가"로 이어 간다.

---

## 에이전트 하네스 (`Day06_LLM_Agent/pde5_harness/`)

노트북과 별개로, **검증 게이트를 갖춘 연구 에이전트 하네스**와 그 하네스가 실제로 수행한
연구가 들어 있다. 실습 없이 파일과 결과만 읽는 용도다.

| 무엇 | 어디 |
|------|------|
| 하네스 계약·게이트·스킬 19종 | [`pde5_harness/CLAUDE.md`](Day06_LLM_Agent/pde5_harness/CLAUDE.md) · [`.claude/skills/`](Day06_LLM_Agent/pde5_harness/.claude/skills) |
| **연구 보고서** (논문 형식 28쪽, n=163) | [`sample_run/report/report_controlled.md`](Day06_LLM_Agent/pde5_harness/sample_run/report/report_controlled.md) |
| 학회 포스터 (A0) · 구두발표 덱 (16:9) | [`sample_run/report/docs/`](Day06_LLM_Agent/pde5_harness/sample_run/report/docs) |
| **⭐ Claude for Science 용 프롬프트** | [`prompts/claude_for_science.md`](Day06_LLM_Agent/pde5_harness/prompts/claude_for_science.md) |

### 프롬프트만 필요하다면

[**`prompts/claude_for_science.md`**](Day06_LLM_Agent/pde5_harness/prompts/claude_for_science.md)
하나만 열면 된다. 로컬 하네스 파일 없이(Claude for Science, 일반 Claude 대화 등) 같은
연구를 수행하기 위한 프롬프트로, 하네스가 파일로 강제하던 계약을 프롬프트 본문이 대신 진다.

- **A. 마스터 프롬프트** (약 2,000단어) — 그대로 붙여넣는다
- **B. 각 조항이 무엇을 막는가** — 21개 조항을 실제 실패 사건과 짝지은 표
- **C. 다른 표적으로 바꿀 때** — 치환할 자리 · 이 연구에 맞춰진 값 · 바꾸면 안 되는 다섯

> 이 연구는 같은 표적·같은 도구로 **세 판본**을 냈고 결론이 세 번 달랐다. 세 번 모두 자동
> 검증 게이트를 전부 통과했으며, 오류를 잡은 것은 매번 외부 비평이었다. 프롬프트의 조항은
> 전부 그 과정에서 한 번씩 틀린 뒤에 추가된 것이다.

---

## 실행 환경

- 모든 노트북은 **Google Colab**에서 바로 실행됩니다(각 노트북 첫 셀에서 필요한 패키지 자동 설치).
- Colab 메뉴 > 런타임 > 런타임 유형 변경 > T4 GPU 로 GPU를 켤 수 있습니다.

| 노트북 | GPU | API 키 / 로그인 | 내려받는 용량 | 비고 |
|---|---|---|---|---|
| `t042_molt5` | 선택 (CPU 가능) | 불필요 | 약 6GB (molt5-large 2종) | CPU에서 프롬프트 1건당 약 15초 |
| `t110_esm2` | 권장 (CPU 가능) | 불필요 | 약 2.6GB (PepMLM-650M) | CPU에서 binder 1개당 약 40초 |
| `t111_esm3` | 선택 (CPU 가능) | 불필요 | 약 2.8GB (esm3-sm-open-v1) | CPU에서 각 생성 단계 15초 내외 |
| `t050_rag` | **강력 권장** | 불필요 (기본 Qwen3) | 약 3.9GB + PDF 27MB | Gemma를 쓰려면 HF 라이선스 동의 + 토큰 필요 |
| `nb1_scientific_agent` | 불필요 | 불필요 (무-키 경로) · Gemini 무료 티어 키는 선택 | 없음 (`rdkit`·`requests` 만) | 외부 호출은 전부 무료 공개 API — 실시간 시연에 가장 안전 |
| `T004_llm_coscientist` | 불필요 | **Google Gemini API 키 필수** | — | https://aistudio.google.com/app/apikey |

- API 키는 노트북 안에서 `getpass` 로 입력하며 파일에 저장되지 않습니다.

### 수업 중 시연할 때 (Colab 기준)

1. **콜드 스타트를 미리 끝내 두세요.** 시연 시점에 노트북을 처음 열면 설치와 모델 내려받기만
   몇 분 걸립니다(t042 약 6GB · t050 약 3.9GB + PDF 27MB · t111 약 2.8GB · t110 약 2.6GB).
   1교시 시작할 때 쓸 노트북 탭을 모두 열어 첫 셀만 돌려 두면 됩니다. 무료 Colab 은 유휴 상태로
   약 90분이면 런타임이 끊기니 탭을 가끔 만져 주세요.
2. **파인튜닝(DPO/ORPO)은 노트북 없이 다룹니다.** 2교시 PEFT·비용 비교(p29–30)는 슬라이드로만
   진행합니다 — 수업 시간 안에 돌 수 있는 정렬 학습 실습이 없습니다.
3. **`T004` 는 수업 중에 키를 입력하게 됩니다.** `getpass` 로 Gemini 키를 받는 화면이 그대로
   투사되는데, 이것은 4교시에서 "하지 말 것"으로 가르치는 항목과 부딪칩니다. 수업 전에 미리
   입력해 런타임을 살려 두거나, 저장된 출력만 넘기세요.
4. **라이브로 돌릴 것은 `nb1` 하나로 좁히는 편이 안전합니다.** GPU·API 키가 모두 불필요하고
   외부 호출이 전부 무료 공개 API(PubMed·UniProt·ChEMBL·RDKit) 라 실패 지점이 거의 없습니다.
   무료 GPU 할당은 보장되지 않으므로 t050 은 CPU 로 떨어질 경우를 대비해 저장 출력(38셀 중 31)을
   그대로 넘기는 경로를 준비해 두세요.

## 출처

본 실습은 [`fourmodern/2025_aidrugdiscovery`](https://github.com/fourmodern/2025_aidrugdiscovery)에서 LLM·에이전트 관련 항목을 선별·재구성한 것입니다.
