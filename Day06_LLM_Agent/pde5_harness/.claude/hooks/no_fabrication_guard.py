#!/usr/bin/env python3
"""PostToolUse 훅 — 무-날조 리마인더(경고형, 비차단).

도구 출력에서 날조 의심 패턴(근거 없는 수치 표현, 가짜 식별자 힌트)을 가볍게 스캔해
에이전트에 리마인더를 주입한다. 차단하지 않으며(데모용), stdin이 없거나 파싱 실패해도 exit 0.
"""
import sys, json, re

# 날조 의심 표현(한/영). 실제 수치 대신 추정/근사를 암시하는 문구.
SUSPICIOUS_PHRASES = [
    "추정", "예상 effect", "예상 IC50", "대략", "약 ", "typical value", "plausible",
    "estimated", "predicted value", "representative sample", "assumed", "가정하면",
]


def emit(context: str):
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
        sys.exit(0)

    resp = payload.get("tool_response", payload.get("toolResponse", {}))
    if isinstance(resp, dict):
        text = resp.get("stdout") or resp.get("output") or json.dumps(resp, ensure_ascii=False)
    elif isinstance(resp, str):
        text = resp
    else:
        text = raw

    hits = [p for p in SUSPICIOUS_PHRASES if p in text]
    # 가짜 ChEMBL ID 힌트: CHEMBL 뒤 숫자가 아닌 경우(형식 오류) — 가벼운 경고.
    bad_id = bool(re.search(r"CHEMBL(?![0-9])", text))

    if hits or bad_id:
        parts = ["[no_fabrication_guard] 리마인더: 무-날조 규칙."]
        if hits:
            parts.append(f"의심 표현 감지: {hits}. 수치·효과크기는 도구 실계산값만 인용하고, "
                         f"근거 없으면 '확인 필요'로 표기하세요.")
        if bad_id:
            parts.append("형식이 이상한 식별자(CHEMBL...) 감지: 실 ID(도구 반환값)만 사용하세요.")
        emit(" ".join(parts))
    else:
        # 조용한 통과지만 상시 원칙 리마인더는 남긴다(데모 가시성).
        emit("[no_fabrication_guard] OK — 수치·ID·PMID는 도구 근거만. 불확실은 '확인 필요'로 표기(무-날조).")
    sys.exit(0)


if __name__ == "__main__":
    main()
