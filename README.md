# 2026 AI 신약개발 — 대규모 언어모델 & 에이전트 (6일차 실습)

재직자 과정 **6일차: 대규모 언어모델과 신약개발(LLM·에이전트)** 의 실습 자료입니다.
2025 워크숍(`fourmodern/2025_aidrugdiscovery`)에서 **LLM·에이전트 관련 실습만** 선별해 재구성했습니다.

이론(슬라이드) 5교시 + 본 실습(3H) 구성이며, 모든 노트북은 Google Colab에서 바로 실행됩니다.

---

## 실습 목록 (`Day06_LLM_Agent/`)

| 교시 매핑 | 노트북 | 주제 | 설명 | Colab |
|------|------|------|------|-------|
| 1교시 · 과학 언어모델(화학) | t042_molt5 | MolT5 | T5 기반 분자↔텍스트 상호 변환(분자 설명 생성 / 텍스트→분자) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/t042_molt5.ipynb) |
| 1교시 · 과학 언어모델(단백질) | t110_esm2 | ESM-2 | ESM-2 단백질 언어모델로 펩타이드 결합 최적화(masked LM 스코어링) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/t110_esm2_peptide_optimization_tutorial.ipynb) |
| 1교시 · 과학 언어모델(단백질) | t111_esm3 | ESM3 | 생성형 멀티모달 단백질 모델 — 서열 설계(inpainting) + 구조 예측 + 신뢰도(pTM/pLDDT) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/t111_esm3_protein_design.ipynb) |
| 2교시 · 활용·RAG | t050_rag | RAG | PDF 문서 기반 로컬 RAG 파이프라인 구축(신약개발 질의응답) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/t050_simple-local-rag.ipynb) |
| 3교시 · 에이전트 | T004_llm_coscientist | LLM Co-Scientist | LangChain + LLM 기반 NSCLC 표적 탐색 AI 에이전트 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/T004_nsclc_llm_coscientist.ipynb) |
| (심화) LLM 정렬·파인튜닝 | qwen_tcga | Qwen 파인튜닝 | Qwen3를 TCGA 데이터에 선호도 최적화(DPO/ORPO)로 파인튜닝. GPU 필수 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fourmodern/2026_aidrugdiscovery/blob/main/Day06_LLM_Agent/qwen_tcga_advanced_2024.ipynb) |

> 이론 흐름(참고): ① LLM 원리 & 과학 언어모델 → ② 활용 & 환각 통제(RAG) → ③ 에이전트 → ④ 하네스·보안 → ⑤ 자율 과학.

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
| `T004_llm_coscientist` | 불필요 | **Google Gemini API 키 필수** | — | https://aistudio.google.com/app/apikey |
| `qwen_tcga` | **필수** | 불필요 | 대용량 | 심화 자료 |

- API 키는 노트북 안에서 `getpass` 로 입력하며 파일에 저장되지 않습니다.

## 출처

본 실습은 [`fourmodern/2025_aidrugdiscovery`](https://github.com/fourmodern/2025_aidrugdiscovery)에서 LLM·에이전트 관련 항목을 선별·재구성한 것입니다.
