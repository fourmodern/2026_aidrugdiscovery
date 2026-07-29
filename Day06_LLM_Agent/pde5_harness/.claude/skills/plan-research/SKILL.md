---
name: plan-research
description: 연구 목표를 실행 전 구조화된 JSON 계획으로 변환한다. 모든 실행의 첫 단계 — 계획 없이 도구 실행 금지.
---
# plan-research

## 언제 사용하는가
워크플로의 **맨 처음**. 사용자의 자연어 목표(예: "PDE5 저해제 후보 조사해 보고서 써줘")를 받으면,
어떤 도구도 실행하기 전에 이 스킬로 연구 계획을 구조화한다. CLAUDE.md의 PLAN 단계에 해당한다.

## 무엇을 하는가
아래 필드를 가진 **JSON 계획**을 작성한다:
- `objective` — 한 문장 연구 목표 (표적 PDE5A, 후보 저해제 조사·평가·보고).
- `questions` — 답해야 할 구체적 질문 목록 (표적 실재? 활성물질? 물성·약물성? PDE6 선택성?).
- `tools` — 사용할 스킬/스크립트 순서 (`target-lookup` → `chembl-actives` → `mol-properties` → `selectivity-check` → `report-writer`).
- `success_criteria` — 성공 판정 기준 (각 단계 `verification.passed=true`, 보고서에 근거 없는 수치 0).
- `stopping_rules` — 중단 규칙 (단계 2회 재시도 후 실패 시 "미확인/플래그" 처리 후 진행).
- `verification_criteria` — 단계별 검증 규칙 (CLAUDE.md "단계별 검증 기준" 표 인용).

## 어떻게
전용 스크립트는 없다. 위 스키마대로 JSON을 직접 작성하고 `outputs/plan.json` 에 저장한다.
표적·사례 약물 사실은 CLAUDE.md의 "목표(고정 예시)"를 근거로 하되, 수치·ID는 이후 도구 단계에서 확정한다.

## 반환 검증
- JSON이 위 6개 키를 모두 포함하는지 확인.
- `tools` 순서가 CLAUDE.md 워크플로와 일치하는지 확인.
- 계획 단계에서는 수치·ChEMBL ID·PMID를 확정하지 않는다(도구 단계의 몫). 계획에 미리 수치를 적지 말 것(무-날조).
