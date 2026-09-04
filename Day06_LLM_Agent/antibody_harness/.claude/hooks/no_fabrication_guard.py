#!/usr/bin/env python3
"""PostToolUse 훅 — 무-날조 리마인더(경고형, 비차단).

도구 출력에서 날조 의심 패턴(근거 없는 수치 표현, 가짜 식별자 힌트)을 가볍게 스캔해
에이전트에 리마인더를 주입한다. 차단하지 않으며(데모용), stdin이 없거나 파싱 실패해도 exit 0.
"""
import sys, json, re

# 날조 의심 표현(한/영). 실제 값 대신 추정/근사를 암시하는 문구.
SUSPICIOUS_PHRASES = [
    "추정", "대략", "typical value", "plausible", "estimated", "predicted value",
    "representative sample", "assumed", "가정하면", "예상 친화도", "예상 KD",
    "아마도", "일반적으로 알려진 값", "approximately the sequence", "likely sequence",
]

# 항체 하네스 특화: 근거 없이 서열/식별자/설계점수를 만들어내는 패턴.
SEQ_CLAIM_PHRASES = [
    "대표적인 서열", "전형적인 CDR", "일반적인 VH 서열", "일반적인 VL 서열",
    "typical CDR", "canonical sequence would be", "the sequence is roughly",
]
GPU_CLAIM_PHRASES = [
    "설계를 완료했", "designed sequences:", "ipTM =", "pAE =",
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
    seq_hits = [p for p in SEQ_CLAIM_PHRASES if p in text]

    # 형식이 어긋난 식별자 탐지 (실 ID 는 도구 반환값이어야 한다).
    #  - UniProt: 6 또는 10자 영숫자.  - PDB entry: 숫자1 + 영숫자3.
    bad_pdb = bool(re.search(r"\bPDB[: ]+(?![0-9][A-Za-z0-9]{3}\b)[A-Za-z0-9]{1,6}\b", text))

    # dry-run 없이 설계 수치를 말하면서 GPU 미사용을 함께 보고한 경우 → 모순 경고.
    contradiction = (("CUDA 사용 불가" in text or "GPU 미검출" in text)
                     and any(p in text for p in GPU_CLAIM_PHRASES))

    if hits or seq_hits or bad_pdb or contradiction:
        parts = ["[no_fabrication_guard] 리마인더: 무-날조 규칙."]
        if hits:
            parts.append(f"의심 표현 감지: {hits}. 수치는 도구/모델 실계산값만 인용하고, "
                         f"근거 없으면 '확인 필요'로 표기하세요.")
        if seq_hits:
            parts.append(f"서열 창작 의심 표현: {seq_hits}. 아미노산 서열은 UniProt/RCSB 또는 "
                         f"설계 모델이 반환한 값만 사용하세요.")
        if bad_pdb:
            parts.append("형식이 이상한 PDB ID 감지: 실 ID(도구 반환값)만 사용하세요 "
                         "(4자: 숫자1 + 영숫자3).")
        if contradiction:
            parts.append("모순 경고: GPU 미사용이 보고되었는데 설계 결과 수치가 함께 나타납니다. "
                         "GPU 없이 설계 서열/ipTM/pAE 를 생성하지 마세요 — dry-run 결과만 보고하고 "
                         "'RunPod GPU 에서 검증 필요'로 표기하십시오.")
        emit(" ".join(parts))
    else:
        emit("[no_fabrication_guard] OK — 서열·ID·수치·API 함수명은 도구/공식문서 근거만. "
             "GPU 미보유 시 설계는 dry-run 만. 불확실은 '확인 필요'로 표기(무-날조).")
    sys.exit(0)


if __name__ == "__main__":
    main()
