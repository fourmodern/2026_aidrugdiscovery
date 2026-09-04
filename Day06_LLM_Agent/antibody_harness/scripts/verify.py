"""단계별 검증 공용 유틸 — 과학적 엄밀성/무-날조 게이트 (항체 하네스).

각 도구 결과는
    {"result": ...,
     "provenance": {"source": ..., "query": ..., "timestamp": ...},
     "verification": {"passed": bool, "checks": [{"check":..., "passed":...}], "notes": ...}}
형태를 갖춰야 한다. 이 모듈은 그 계약을 강제하는 헬퍼를 제공한다.

무-날조 원칙:
- API 실패 = 빈 결과 + verification.passed=false. 값을 지어내지 않는다.
- 오프라인 폴백은 "실제 공개 데이터 + 출처 명시" 일 때만 허용하며 provenance.source 에 표기한다.

사용: python scripts/verify.py            # 자가 점검
      python scripts/verify.py --gate     # stdin 봉투를 읽어 게이트 판정 (exit code 0/1)
"""
from __future__ import annotations

import datetime
import json
import re
import sys

# 표준 20 아미노산 + 모호 코드(X: unknown, B/Z/J: ambiguous, U: selenocysteine, O: pyrrolysine)
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
AA_EXTENDED = AA20 | set("XBZJUO")


def now() -> str:
    """재현성: 호출 시각 기록(provenance)."""
    return datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z")


def make_result(result, source: str, query: str, checks, notes: str = "") -> dict:
    """표준 결과 봉투 생성.

    checks: [(설명, 통과여부), ...]. 하나라도 False 면 verification.passed=False.
    """
    checks = list(checks)
    passed = bool(checks) and all(bool(ok) for _, ok in checks)
    return {
        "result": result,
        "provenance": {"source": source, "query": query, "timestamp": now()},
        "verification": {
            "passed": passed,
            "checks": [{"check": c, "passed": bool(ok)} for c, ok in checks],
            "notes": notes,
        },
    }


def gate(envelope: dict, step: str) -> bool:
    """검증 게이트: passed=False 면 False 반환(다음 단계 진행 차단)."""
    v = (envelope or {}).get("verification", {})
    if not v.get("passed"):
        failed = [c["check"] for c in v.get("checks", []) if not c.get("passed")]
        sys.stderr.write(
            f"[VERIFY-GATE:{step}] FAILED → 다음 단계 진행 금지. 실패 항목: {failed}\n"
        )
        return False
    sys.stderr.write(f"[VERIFY-GATE:{step}] PASSED\n")
    return True


def valid_protein_seq(seq, min_len: int = 10) -> bool:
    """단백질 서열 유효성(무-날조: 형식이 깨졌으면 신뢰하지 않는다).

    - 문자열이고 min_len 이상
    - 모든 문자가 확장 아미노산 알파벳 안에 있음
    - 표준 20 아미노산 비율 >= 90% (핵산 서열이나 placeholder 문자열 배제)
    """
    if not isinstance(seq, str):
        return False
    s = seq.strip().upper()
    if len(s) < min_len:
        return False
    if any(ch not in AA_EXTENDED for ch in s):
        return False
    std = sum(1 for ch in s if ch in AA20)
    if std / len(s) < 0.9:
        return False
    # ACGT 만으로 이루어진 문자열은 핵산일 가능성이 높다 → 단백질로 인정하지 않는다.
    if set(s) <= set("ACGTU"):
        return False
    return True


def valid_uniprot_acc(acc) -> bool:
    """UniProt accession 형식 점검(가짜 accession 조기 탐지)."""
    if not isinstance(acc, str):
        return False
    pat = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
    return bool(re.match(pat, acc.strip().upper()))


def valid_pdb_id(pdb) -> bool:
    """PDB entry ID 형식 점검(4자: 숫자 1개 + 영숫자 3개)."""
    if not isinstance(pdb, str):
        return False
    return bool(re.match(r"^[0-9][A-Za-z0-9]{3}$", pdb.strip()))


def numbers_backed(text: str, allowed_numbers) -> list:
    """보고서 텍스트의 숫자가 도구 결과(allowed_numbers)에 근거하는지 점검.

    근거 없는 수치 목록을 반환(빈 리스트면 통과). LLM 생성 수치 탐지용.
    """
    allowed = {str(a) for a in allowed_numbers}
    found = re.findall(r"\d+\.?\d*", text or "")
    return [n for n in found if n not in allowed and len(n) >= 2]


def collect_numbers(obj, acc=None) -> set:
    """봉투/딕셔너리에서 모든 수치를 문자열로 수집(report 단계 numbers_backed 용)."""
    acc = set() if acc is None else acc
    if isinstance(obj, dict):
        for v in obj.values():
            collect_numbers(v, acc)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            collect_numbers(v, acc)
    elif isinstance(obj, bool):
        pass
    elif isinstance(obj, (int, float)):
        acc.add(str(obj))
    elif isinstance(obj, str):
        for n in re.findall(r"\d+\.?\d*", obj):
            acc.add(n)
    return acc


def emit(envelope: dict) -> None:
    """표준 봉투를 stdout 으로 출력(모든 스크립트 공통 출구)."""
    print(json.dumps(envelope, ensure_ascii=False, indent=2))


def _self_test() -> int:
    demo = make_result({"x": 1}, "self-test", "n/a", [("형식", True)])
    print(json.dumps(demo, ensure_ascii=False, indent=2))
    print("gate:", gate(demo, "self-test"))
    print("valid_protein_seq(EVQLVESGGG...):", valid_protein_seq("EVQLVESGGGLVQPGGSLRLSCAAS"))
    print("valid_protein_seq(ACGTACGTACGT):", valid_protein_seq("ACGTACGTACGT"))
    print("valid_uniprot_acc(P04626):", valid_uniprot_acc("P04626"))
    print("valid_pdb_id(1N8Z):", valid_pdb_id("1N8Z"))
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--gate":
        try:
            env = json.load(sys.stdin)
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"[verify] stdin JSON 파싱 실패: {exc}\n")
            raise SystemExit(1)
        raise SystemExit(0 if gate(env, "stdin") else 1)
    raise SystemExit(_self_test())
