---
name: report-writer
description: 앞 단계들의 표준 봉투(result+provenance+verification)를 모아 근거 인용 IMRAD 보고서를 작성해 outputs/에 저장한다.
---
# report-writer

## 언제 사용하는가
모든 실행 단계(target-lookup → chembl-actives → mol-properties → selectivity-check)가
끝난 뒤 마지막 REPORT 단계. 각 단계의 표준 봉투를 종합해 과학 보고서를 작성한다.

## 무엇을 하는가
IMRAD 구조의 Markdown 보고서를 작성한다:
- **Introduction** — 표적 PDE5A(UniProt O76074, ChEMBL CHEMBL1827), cGMP-NO 경로, 사례 약물.
- **Methods** — 사용한 도구·스킬 파이프라인과 각 단계 검증 기준.
- **Results** — 활성물질 표(ChEMBL ID·SMILES·활성값), 물성/약물성 표(MW/logP/QED/SA/게이트), 선택성 정성 평가.
- **검증 로그** — 단계별 `verification`(passed·checks·provenance) 요약.
- **한계(Discussion/Limitations)** — 결과=가설, 신약 발견 아님, 도킹/LLM 추론은 가설, 미확인 항목 플래그.
- **참고문헌** — 실재 인용만(Boolell 1996 IJIR; Ghofrani 2006 Nat Rev Drug Discov; UniProt; ChEMBL).

## 어떻게
앞 단계 JSON 봉투들을 입력으로 받아 종합한다. 모든 수치는 도구 결과(`result`)에서만 인용한다.
저장: `outputs/report_pde5.md`.

## 반환 검증 (무-날조 필수)
- 보고서의 모든 수치가 앞 단계 도구 결과에 실제로 존재하는지 `scripts/verify.py`의
  `numbers_backed(text, allowed_numbers)` 로 점검 — 근거 없는 수치 목록이 **비어 있어야** 통과.
  (`allowed_numbers` = 모든 단계 봉투에서 추출한 수치 집합.)
- 각 사실 주장에 출처(ChEMBL ID / UniProt / PMID / DOI)가 있는지 확인. 없으면 "확인 필요" 표기.
- 참고문헌이 실재하는지 확인(가짜 PMID/DOI 금지).
- 미확인/플래그 항목이 한계 섹션에 명시되었는지 확인.
