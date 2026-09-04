"""단계별 검증 공용 유틸 — 과학적 엄밀성/무-날조 게이트.

각 도구 결과는 {"result":..., "provenance":..., "verification":{"passed":bool,"checks":[...],"notes":...}}
형태를 갖춰야 한다. 이 모듈은 그 검증을 강제하는 헬퍼를 제공한다.
"""
from __future__ import annotations
import json, sys, datetime


def now() -> str:
    # 재현성: 호출 시각 기록(provenance)
    return datetime.datetime.utcnow().isoformat() + "Z"


def make_result(result, source: str, query: str, checks: list[tuple[str, bool]], notes: str = ""):
    """표준 결과 봉투 생성. checks = [(설명, 통과여부), ...]."""
    passed = all(ok for _, ok in checks)
    return {
        "result": result,
        "provenance": {"source": source, "query": query, "timestamp": now()},
        "verification": {
            "passed": passed,
            "checks": [{"check": c, "passed": ok} for c, ok in checks],
            "notes": notes,
        },
    }


def gate(envelope: dict, step: str):
    """검증 게이트: passed 를 판정해 bool 로 돌려준다.

    [정정 2026-09-04] 이전 docstring 은 "비정상 종료"라고 적었으나 이 함수는
    raise 도 sys.exit 도 하지 않는다. 차단은 전적으로 호출자가 반환값을 보고
    수행한다. 호출자는 반드시 `if not gate(...): return 1` 형태로 써야 한다.
    """
    v = envelope.get("verification", {})
    if not v.get("passed"):
        failed = [c["check"] for c in v.get("checks", []) if not c["passed"]]
        sys.stderr.write(f"[VERIFY-GATE:{step}] FAILED → 다음 단계 진행 금지. 실패 항목: {failed}\n")
        return False
    sys.stderr.write(f"[VERIFY-GATE:{step}] PASSED\n")
    return True


def valid_smiles(smiles: str) -> bool:
    """RDKit로 SMILES 유효성 확인(무-날조: 파싱 안 되면 가짜)."""
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        # RDKit 미설치 환경: 최소 형식 체크만(데모)
        return isinstance(smiles, str) and len(smiles) > 3


def numbers_backed(text: str, allowed_numbers: set[str]) -> list[str]:
    """보고서 텍스트의 숫자가 도구 결과(allowed_numbers)에 근거하는지 점검.
    근거 없는 수치 목록을 반환(빈 리스트면 통과). LLM 생성 수치 탐지용."""
    import re
    found = re.findall(r"\d+\.?\d*", text)
    unbacked = [n for n in found if n not in allowed_numbers and len(n) >= 2]
    return unbacked


if __name__ == "__main__":
    # 자가 점검
    demo = make_result({"x": 1}, "self-test", "n/a", [("형식", True)])
    print(json.dumps(demo, ensure_ascii=False, indent=2))
    print("gate:", gate(demo, "self-test"))
