#!/usr/bin/env python3
"""PostToolUse 훅 — 도구 출력에 provenance/verification 필드 존재를 점검(경고형, 비차단).

Bash 도구가 표준 봉투 JSON({result, provenance, verification})을 출력했는지 확인한다.
- provenance/verification 이 있으면: 통과(추가 안내 없음 또는 검증 상태 요약).
- 없으면: additionalContext 로 "검증 필드 누락 경고"를 에이전트에 주입(차단 아님).

데모용: stdin JSON이 없거나 파싱 실패해도 죽지 않고 exit 0. Claude Code 훅 스키마
(hookSpecificOutput.additionalContext)를 사용한다.
"""
import sys, json


def emit(context: str):
    """PostToolUse additionalContext 로 컨텍스트 주입."""
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": context,
        }
    }, ensure_ascii=False))


def main():
    try:
        raw = sys.stdin.read()
    except Exception:
        sys.exit(0)
    if not raw or not raw.strip():
        sys.exit(0)

    try:
        payload = json.loads(raw)
    except Exception:
        # 훅 입력이 JSON이 아니어도 데모에서는 조용히 통과.
        sys.exit(0)

    # 도구 출력 텍스트 추출(Claude Code PostToolUse 페이로드 구조에 관대하게 대응).
    tool_output = ""
    resp = payload.get("tool_response", payload.get("toolResponse", {}))
    if isinstance(resp, dict):
        tool_output = resp.get("stdout") or resp.get("output") or json.dumps(resp, ensure_ascii=False)
    elif isinstance(resp, str):
        tool_output = resp
    if not tool_output:
        tool_output = raw

    has_prov = '"provenance"' in tool_output
    has_ver = '"verification"' in tool_output

    if has_prov and has_ver:
        # 통과 — verification.passed 상태를 가볍게 요약(있으면).
        note = "[verify_provenance] OK — provenance/verification 필드 확인됨."
        if '"passed": false' in tool_output or '"passed":false' in tool_output:
            note += (" 주의: verification.passed=false 인 단계가 있음 → 다음 단계 진행 전 "
                     "재시도/플래그 검토. (오프라인 캐시·GPU 미보유 dry-run 도 여기에 해당)")
        if '"heuristic (approximate)"' in tool_output or "heuristic (approximate)" in tool_output:
            note += (" CDR 은 휴리스틱 근사로 추출되었습니다 — 보고서에 'IMGT 정확 번호'라고 "
                     "쓰지 말고 근사임을 명시하십시오.")
        if '"mode": "dry-run"' in tool_output:
            note += (" DRY-RUN 출력입니다 — 설계가 실제로 수행되지 않았습니다. "
                     "설계 서열/점수를 보고하지 마십시오.")
        emit(note)
    else:
        missing = []
        if not has_prov:
            missing.append("provenance")
        if not has_ver:
            missing.append("verification")
        # 표준 봉투가 아닌 일반 Bash 출력(ls 등)에는 과도 경고를 피하되, 힌트는 남긴다.
        emit(f"[verify_provenance] 경고: 도구 출력에 {', '.join(missing)} 필드가 없습니다. "
             f"과학 단계(antigen-lookup/antibody-search/cdr-analysis/developability/humanness/design-esmfold2/design-rfantibody/compare-designs)의 출력이라면 "
             f"표준 봉투 {{result, provenance, verification}} 형식을 갖췄는지 확인하세요(무-날조·검증 게이트).")
    sys.exit(0)


if __name__ == "__main__":
    main()
