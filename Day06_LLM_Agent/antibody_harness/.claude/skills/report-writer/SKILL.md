---
name: report-writer
description: 설계·평가 단계의 표준 봉투를 모아 근거 인용 IMRAD 보고서를 작성해 outputs/ 에 저장한다.
---
# report-writer

## 언제 사용하는가
마지막 단계. 앞 단계들의 봉투(`outputs/*.json`)를 종합한다.

## 무엇을 쓰는가 (IMRAD)
- **Introduction** — 항원(UniProt accession, gene, length, 세포외 도메인 범위), 왜 이 표적인가.
- **Methods** — 실제로 실행한 스킬/스크립트 파이프라인과 각 단계 검증 기준.
  설계를 했다면 경로 A/B 의 모델·버전·출처 URL·라이선스를 명시.
  **CDR 추출 방법(anarci 정식 vs 휴리스틱 근사)을 반드시 구분해 서술.**
- **Results** —
  - 알려진 항체 표 (PDB ID / 해상도 / 사슬 / CDR / 서열 길이)
  - developability 표 (MW / pI / GRAVY / instability / liability 히트 · CDR 내부)
  - humanness 표 (nearest germline / identity %)
  - 설계를 했다면: 경로별 설계 수, 경로 고유 신뢰도 지표(원본 이름), 공통 축 비교표
- **검증 로그** — 단계별 `verification` (passed / checks / provenance) 요약.
  `passed=false` 인 단계는 **반드시 표에 남긴다** (오프라인 캐시·dry-run 포함).
- **Limitations** — 아래 "한계 필수 항목" 전부.
- **References** — 실재하는 것만 (UniProt accession, PDB ID, GitHub URL, bioRxiv DOI).

## 한계 필수 항목 (하나도 빼지 말 것)
1. CDR 이 휴리스틱 근사로 추출됐다면 그 사실과 검증 범위(트라스투주맙·퍼투주맙에서만 확인).
2. liability 는 규칙 기반 플래그이지 예측 모델이 아니다.
3. germline identity 는 면역원성(ADA) 예측이 아니며, local alignment 특성상 프레임워크
   유사도로 편향된다.
4. ipTM / pAE 는 구조 신뢰도이지 결합상수(KD)가 아니다.
5. RFantibody 저자 인용: "The lack of an effective filter is the main limitation of the
   RFantibody pipeline at the moment." + 일반적으로 10k 규모 캠페인 필요.
6. GPU 미보유로 dry-run 만 수행했다면 **"설계를 실행하지 않았음"** 을 명시.
7. 두 경로의 신뢰도 지표를 통합 점수로 합치지 않았다는 사실과 그 이유.
8. 모든 in silico 결과는 **가설**이며 발현·정제·SPR/BLI·DSF·SEC·in vivo 로만 검증된다.

## 무-날조 검증 (필수)
- 보고서의 모든 수치가 앞 단계 봉투에 실제로 존재하는지
  `scripts/verify.py` 의 `collect_numbers()` + `numbers_backed(text, allowed)` 로 점검 —
  근거 없는 수치 목록이 **비어 있어야** 통과.
- 모든 식별자(PDB ID / UniProt accession / DOI / HuggingFace repo)가 도구 반환값 또는
  공식 문서에서 확인된 것인지 확인. 기억에 의존한 인용 금지.
- `passed=false` 였던 단계의 결과를 "성공"으로 서술하지 않았는지 확인.
- API 키(`RUNPOD_API_KEY`, `HF_TOKEN`) 값이 보고서에 들어가지 않았는지 확인 — **이름만 허용.**

## 저장
`outputs/report_antibody.md`
