"""CDR 추출/번호매김 — anarci/abnumber 가 있으면 정식 번호매김, 없으면 문서화된 휴리스틱.

무-날조 정책 (중요):
- `anarci` / `abnumber` 미설치 시 결과는 **근사 휴리스틱**이며, `method` 필드에
  `"heuristic (approximate)"` 로 명시된다. **정확한 IMGT 번호라고 주장하지 않는다.**
- 경계를 찾지 못한 CDR 은 `null` 로 남긴다. 절대 서열을 만들어내지 않는다.
- 휴리스틱 규칙(보존 Cys / FR2 Trp / J-region 모티프)은 `seq_utils.py` 상단에 문서화되어 있다.

휴리스틱 검증(개발 시 확인한 사실): 트라스투주맙(PDB 1N8Z) · 퍼투주맙(PDB 1S78) 사슬에서
Kabat CDR 6종이 문헌 공개값과 일치하게 재현되었다. 그러나 모든 항체에 대해 보장되지 않는다.

사용:
    python scripts/cdr_analysis.py "EVQLVESGG..."                # 서열 직접
    python scripts/cdr_analysis.py chains.fasta                  # FASTA
    python scripts/antibody_search.py | python scripts/cdr_analysis.py --stdin
    python scripts/cdr_analysis.py --stdin --scheme kabat        # anarci 설치 시 scheme 지정

출력: 표준 결과 봉투(JSON).
"""
from __future__ import annotations

import json
import sys

from seq_utils import (cdr_positions, cdrs_heavy, cdrs_light, classify_chain,
                       clean_seq, numbering_backend, read_sequences, split_scfv,
                       try_anarci, variable_domain)
from verify import emit, make_result, valid_protein_seq

HEURISTIC_METHOD = "heuristic (approximate) — Kabat-like, conserved-motif regex (seq_utils)"
HEURISTIC_CAVEAT = (
    "anarci/abnumber 미설치 → 보존 모티프 기반 근사 추출. IMGT/Kabat 공식 번호가 아니며, "
    "비정형 프레임워크·삽입(insertion)·비인간 항체에서는 경계가 틀릴 수 있음. "
    "정확한 번호매김이 필요하면 `pip install anarci` 또는 `pip install abnumber` 후 재실행."
)


def analyze_domain(seq: str, seq_id: str, scheme: str) -> dict:
    """단일 가변영역(중쇄 또는 경쇄)의 CDR 추출."""
    s = clean_seq(seq)
    cls = classify_chain(s)
    chain_type = cls["chain_type"]

    # 1) 정식 번호매김 시도
    cdrs, method = try_anarci(s, scheme=scheme)
    warnings = []
    if cdrs is None:
        # 2) 문서화된 휴리스틱 폴백
        if chain_type == "heavy":
            raw = cdrs_heavy(s)
        elif chain_type == "light":
            raw = cdrs_light(s)
        else:
            return {
                "id": seq_id, "chain_type": chain_type, "length": len(s),
                "method": HEURISTIC_METHOD, "cdrs": None, "cdr_lengths": None,
                "cdr_positions_0based": None,
                "variable_domain": variable_domain(s),
                "classification_evidence": cls["evidence"],
                "warnings": ["항체 가변영역으로 판별되지 않아 CDR 추출을 수행하지 않음 (무-날조)"],
            }
        warnings = raw.pop("warnings", [])
        cdrs = raw
        method = HEURISTIC_METHOD
        warnings.append(HEURISTIC_CAVEAT)

    pos = cdr_positions(s, cdrs)
    return {
        "id": seq_id,
        "chain_type": chain_type,
        "length": len(s),
        "method": method,
        "scheme_requested": scheme,
        "cdrs": cdrs,
        "cdr_lengths": {k: (len(v) if v else None) for k, v in cdrs.items()},
        "cdr_positions_0based": pos or None,
        "variable_domain": variable_domain(s),
        "classification_evidence": cls["evidence"],
        "warnings": warnings,
    }


def analyze(seq: str, seq_id: str, scheme: str) -> list:
    """사슬 하나를 분석. scFv 면 VL/VH 로 분할해 각각 분석."""
    s = clean_seq(seq)
    cls = classify_chain(s)
    if cls["chain_type"] != "scfv":
        return [analyze_domain(s, seq_id, scheme)]

    segments, order = split_scfv(s)
    if not segments:
        return [analyze_domain(s, seq_id, scheme)]
    out = []
    for seg in segments:
        d = analyze_domain(seg["sequence"], f"{seq_id}:{seg['chain_type']}", scheme)
        d["parent_id"] = seq_id
        d["scfv_order"] = order
        d["offset_in_parent_0based"] = seg["start"]
        d.setdefault("warnings", []).append(
            f"scFv 로 판별되어 '{order}' 로 분할 후 분석 (분할 지점은 J-region 모티프 기반 근사)"
        )
        out.append(d)
    return out


def _from_stdin() -> list:
    """antibody_search 봉투에서 항체 사슬 서열 추출."""
    try:
        env = json.load(sys.stdin)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"[cdr_analysis] stdin JSON 파싱 실패: {exc}\n")
        return []
    rows = env.get("result", env) if isinstance(env, dict) else env
    if isinstance(rows, dict):
        rows = [rows]
    out = []
    for r in rows or []:
        if not isinstance(r, dict):
            continue
        if "chains" in r:
            for c in r["chains"]:
                if c.get("role") == "antigen" or c.get("chain_type") == "unknown":
                    continue
                cid = "".join(c.get("auth_asym_ids") or []) or str(c.get("entity_id", "?"))
                out.append({
                    "id": f"{r.get('pdb_id','?')}_{cid}",
                    "sequence": c.get("sequence", ""),
                    "source_pdb": r.get("pdb_id"),
                    "description": c.get("description", ""),
                })
        elif r.get("sequence"):
            out.append({"id": r.get("id", "input"), "sequence": r["sequence"]})
    return out


def main() -> int:
    scheme = "imgt"
    argv = sys.argv[1:]
    if "--scheme" in argv:
        i = argv.index("--scheme")
        if i + 1 < len(argv):
            scheme = argv[i + 1].lower()
            argv = argv[:i] + argv[i + 2:]

    if "--stdin" in sys.argv:
        records = _from_stdin()
        query = f"stdin envelope; scheme={scheme}"
    else:
        args = [a for a in argv if not a.startswith("--")]
        if not args:
            emit(make_result([], "input validation", "no input",
                             [("입력 서열 존재", False)],
                             notes="서열/FASTA 경로를 인자로 주거나 --stdin 을 사용하세요."))
            return 1
        records = []
        for a in args:
            records.extend(read_sequences(a))
        query = f"n_inputs={len(args)}; scheme={scheme}"

    backend = numbering_backend()
    results, errors = [], []
    for rec in records:
        s = clean_seq(rec.get("sequence", ""))
        if not valid_protein_seq(s, min_len=60):
            errors.append(f"{rec.get('id')}: 가변영역으로 보기엔 너무 짧거나 무효 (len={len(s)})")
            continue
        try:
            for d in analyze(s, rec.get("id", "input"), scheme):
                d["source_pdb"] = rec.get("source_pdb")
                d["source_description"] = rec.get("description", "")
                results.append(d)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{rec.get('id')}: 추출 실패 {type(exc).__name__}: {exc}")

    n_full = sum(1 for r in results
                 if r.get("cdrs") and all(v for k, v in r["cdrs"].items()))
    n_h3 = sum(1 for r in results if (r.get("cdrs") or {}).get("CDR-H3"))
    n_l3 = sum(1 for r in results if (r.get("cdrs") or {}).get("CDR-L3"))

    checks = [
        ("가변영역 분석 결과 ≥ 1", len(results) >= 1),
        ("CDR 3종 모두 추출된 도메인 ≥ 1", n_full >= 1),
        ("CDR 길이가 생물학적 범위(1-40) 내",
         all(1 <= (l or 1) <= 40
             for r in results for l in (r.get("cdr_lengths") or {}).values())),
        ("입력 사슬 전부 처리됨 (스킵 0건)", len(errors) == 0),
    ]
    method_note = ("정식 번호매김(anarci/abnumber) 사용" if backend != "heuristic"
                   else "휴리스틱 근사 사용 — IMGT 정확 번호 아님")
    notes = (f"도메인 {len(results)}건 분석 (CDR 3종 완전 추출 {n_full}건, "
             f"H3 {n_h3}건 · L3 {n_l3}건). 백엔드={backend} → {method_note}. "
             f"경계를 못 찾은 CDR 은 null 로 남긴다(무-날조). "
             + (f"처리 실패: {errors}" if errors else ""))
    emit(make_result(results, f"CDR extraction backend={backend}", query, checks, notes=notes))
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result([], "cdr_analysis (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
