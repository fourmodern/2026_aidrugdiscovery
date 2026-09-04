"""Developability 평가 — BioPython 실계산 + 규칙 기반 liability 모티프 스캔.

무-날조 정책:
- 모든 수치는 BioPython `ProteinAnalysis` 실계산값. 예측 모델이 아니다.
- Liability 는 **정규식 규칙**이며 실험적 검증값이 아니다. "위험 가능성 플래그"로만 해석한다.
- BioPython 미설치 시 수치를 만들지 않고 verification.passed=false 로 정직 보고한다.

사용:
    python scripts/developability.py "EVQLVESGG..."                   # 서열 직접
    python scripts/developability.py chains.fasta                     # FASTA 파일
    python scripts/antibody_search.py | python scripts/developability.py --stdin
    python scripts/cdr_analysis.py --stdin < ab.json | python scripts/developability.py --stdin

출력: 표준 결과 봉투(JSON).

계산 항목 (ProtParam, 실계산):
    molecular_weight_Da, isoelectric_point_pI, gravy (Kyte-Doolittle 평균),
    instability_index (Guruprasad), aromaticity (Lobry),
    extinction_coefficient (환원/산화 Cys 기준 2값), aa_composition

Liability 규칙 (항체 엔지니어링에서 통용되는 서열 모티프):
    N-glycosylation sequon   N[^P][ST]
    Deamidation              NG (high) / NS, NT, NN, NH (moderate)
    Isomerization            DG (high) / DS, DT, DD (moderate)
    Acid hydrolysis / clip   DP
    Oxidation                M (Met), W (Trp)
    Free (unpaired) cysteine Cys 개수가 홀수
    N-terminal pyroglutamate 1번 잔기가 Q 또는 E
    Hydrophobic patch        Kyte-Doolittle 5-잔기 이동평균 >= 1.5
"""
from __future__ import annotations

import json
import re
import sys

from seq_utils import (KD_HYDROPATHY, cdr_positions, cdrs_heavy, cdrs_light,
                       classify_chain, clean_seq, read_sequences)
from verify import emit, make_result, valid_protein_seq

# (이름, 정규식, 심각도, 설명) — 규칙 기반. 예측 모델 아님.
LIABILITY_RULES = [
    ("N-glycosylation_sequon", r"N[^P][ST]", "high",
     "N-결합 당쇄화 sequon. CDR 내부면 항원 결합·불균질성 위험."),
    ("deamidation_NG", r"NG", "high", "Asn 탈아미드화 고위험 모티프."),
    ("deamidation_moderate", r"N[STNH]", "moderate", "Asn 탈아미드화 중간 위험 모티프."),
    ("isomerization_DG", r"DG", "high", "Asp 이성질화 고위험 모티프."),
    ("isomerization_moderate", r"D[STD]", "moderate", "Asp 이성질화 중간 위험 모티프."),
    ("acid_hydrolysis_DP", r"DP", "moderate", "산성 조건 절단(clip) 위험 모티프."),
    ("oxidation_Met", r"M", "moderate", "Met 산화 위험 잔기."),
    ("oxidation_Trp", r"W", "moderate", "Trp 산화 위험 잔기."),
]

HYDROPHOBIC_WINDOW = 5
HYDROPHOBIC_CUTOFF = 1.5
INSTABILITY_CUTOFF = 40.0  # Guruprasad 1990: >40 이면 in vivo 불안정 예측


def _in_cdr(idx: int, cdr_pos: dict):
    """잔기 인덱스가 어느 CDR 안에 있는지. 아니면 None."""
    for name, (a, b) in (cdr_pos or {}).items():
        if a <= idx < b:
            return name
    return None


def scan_liabilities(seq: str, cdr_pos: dict) -> list:
    """규칙 기반 liability 모티프 스캔. 위치(0-based)와 CDR 내부 여부 표시."""
    hits = []
    for name, pattern, severity, desc in LIABILITY_RULES:
        for m in re.finditer(f"(?=({pattern}))", seq):  # overlapping 허용
            start = m.start()
            motif = m.group(1)
            hits.append({
                "liability": name,
                "motif": motif,
                "severity": severity,
                "position_0based": start,
                "position_1based": start + 1,
                "in_cdr": _in_cdr(start, cdr_pos),
                "rule": pattern,
                "description": desc,
            })
    hits.sort(key=lambda h: (h["position_0based"], h["liability"]))
    return hits


def hydrophobic_patches(seq: str) -> list:
    """Kyte-Doolittle 이동평균으로 소수성 패치(응집 경향 프록시) 검출."""
    n = len(seq)
    if n < HYDROPHOBIC_WINDOW:
        return []
    vals = [KD_HYDROPATHY.get(ch, 0.0) for ch in seq]
    patches, run = [], None
    for i in range(n - HYDROPHOBIC_WINDOW + 1):
        avg = sum(vals[i:i + HYDROPHOBIC_WINDOW]) / HYDROPHOBIC_WINDOW
        if avg >= HYDROPHOBIC_CUTOFF:
            if run is None:
                run = [i, i + HYDROPHOBIC_WINDOW, avg]
            else:
                run[1] = i + HYDROPHOBIC_WINDOW
                run[2] = max(run[2], avg)
        elif run is not None:
            patches.append(run)
            run = None
    if run is not None:
        patches.append(run)
    return [{"start_0based": a, "end_0based": b, "sequence": seq[a:b],
             "max_window_kd": round(c, 3),
             "window": HYDROPHOBIC_WINDOW, "cutoff": HYDROPHOBIC_CUTOFF}
            for a, b, c in patches]


def analyze(seq: str, seq_id: str = "input") -> dict:
    """단일 사슬 developability 분석. BioPython 실계산 + 규칙 스캔."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    s = clean_seq(seq)
    # ProtParam 은 비표준 잔기(X/B/Z/U/O)에서 실패하므로 제거하고 그 사실을 기록한다.
    s_std = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    removed = len(s) - len(s_std)

    cls = classify_chain(s)
    if cls["chain_type"] == "heavy":
        cdrs = cdrs_heavy(s)
    elif cls["chain_type"] == "light":
        cdrs = cdrs_light(s)
    else:
        cdrs = {}
    cdr_pos = cdr_positions(s, cdrs) if cdrs else {}

    pa = ProteinAnalysis(s_std)
    ext_red, ext_ox = pa.molar_extinction_coefficient()
    n_cys = s.count("C")

    props = {
        "length": len(s),
        "molecular_weight_Da": round(pa.molecular_weight(), 2),
        "isoelectric_point_pI": round(pa.isoelectric_point(), 2),
        "gravy_kyte_doolittle": round(pa.gravy(), 4),
        "instability_index": round(pa.instability_index(), 2),
        "instability_verdict": ("unstable (>40)" if pa.instability_index() > INSTABILITY_CUTOFF
                                else "stable (<=40)"),
        "aromaticity": round(pa.aromaticity(), 4),
        "extinction_coefficient_reduced_M1cm1": ext_red,
        "extinction_coefficient_cystines_M1cm1": ext_ox,
        "secondary_structure_fraction_helix_turn_sheet":
            [round(x, 4) for x in pa.secondary_structure_fraction()],
        "n_nonstandard_residues_excluded": removed,
    }

    liabilities = scan_liabilities(s, cdr_pos)
    patches = hydrophobic_patches(s)

    structural = []
    if n_cys % 2 == 1:
        structural.append({
            "liability": "free_unpaired_cysteine", "severity": "high",
            "count": n_cys,
            "positions_1based": [i + 1 for i, ch in enumerate(s) if ch == "C"],
            "description": "Cys 개수가 홀수 → 짝을 이루지 못한 유리 Cys 존재 가능 "
                           "(이량체화·불균질성 위험). 단, 사슬 간 이황화결합을 이루는 "
                           "정상 Cys 일 수 있으므로 구조 맥락 확인 필요.",
        })
    if s[:1] in ("Q", "E"):
        structural.append({
            "liability": "N_terminal_pyroglutamate", "severity": "moderate",
            "residue": s[0], "position_1based": 1,
            "description": "N-말단 Gln/Glu 는 pyroglutamate 로 고리화되어 전하 이질성을 만든다.",
        })

    by_sev = {"high": 0, "moderate": 0}
    for h in liabilities:
        by_sev[h["severity"]] = by_sev.get(h["severity"], 0) + 1
    for h in structural:
        by_sev[h["severity"]] = by_sev.get(h["severity"], 0) + 1

    cdr_liab = [h for h in liabilities if h["in_cdr"]]

    return {
        "id": seq_id,
        "chain_type": cls["chain_type"],
        "chain_classification_evidence": cls["evidence"],
        "sequence": s,
        "properties": props,
        "cdrs": {k: v for k, v in cdrs.items() if k != "warnings"} or None,
        "cdr_positions_0based": cdr_pos or None,
        "cdr_extraction_method": "heuristic (approximate, Kabat-like) — seq_utils",
        "liabilities": liabilities,
        "structural_liabilities": structural,
        "liabilities_in_cdr": cdr_liab,
        "hydrophobic_patches": patches,
        "summary": {
            "n_liability_hits": len(liabilities) + len(structural),
            "n_high_severity": by_sev.get("high", 0),
            "n_moderate_severity": by_sev.get("moderate", 0),
            "n_in_cdr": len(cdr_liab),
            "n_hydrophobic_patches": len(patches),
        },
    }


def _from_stdin() -> list:
    """상류 봉투(antibody_search / cdr_analysis)에서 서열 목록 추출."""
    try:
        env = json.load(sys.stdin)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"[developability] stdin JSON 파싱 실패: {exc}\n")
        return []
    rows = env.get("result", env) if isinstance(env, dict) else env
    out = []
    if isinstance(rows, dict):
        rows = [rows]
    for r in rows or []:
        if not isinstance(r, dict):
            continue
        # antibody_search 형식: {"pdb_id":..., "chains":[{...}]}
        if "chains" in r:
            for c in r["chains"]:
                if c.get("role") == "antigen" or c.get("chain_type") == "unknown":
                    continue
                out.append({"id": f"{r.get('pdb_id','?')}_{''.join(c.get('auth_asym_ids') or [c.get('entity_id','?')])}",
                            "sequence": c.get("sequence", "")})
        # cdr_analysis 형식: {"id":..., "sequence":...}
        elif r.get("sequence"):
            out.append({"id": r.get("id", "input"), "sequence": r["sequence"]})
    return out


def main() -> int:
    if "--stdin" in sys.argv:
        records = _from_stdin()
        query = "stdin envelope"
    else:
        args = [a for a in sys.argv[1:] if not a.startswith("--")]
        if not args:
            emit(make_result([], "input validation", "no input",
                             [("입력 서열 존재", False)],
                             notes="서열 또는 FASTA 파일 경로를 인자로 주거나 --stdin 을 사용하세요."))
            return 1
        records = []
        for a in args:
            records.extend(read_sequences(a))
        query = f"n_inputs={len(args)}"

    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        emit(make_result([], "BioPython ProtParam (MISSING)", query,
                         [("BioPython 설치", False)],
                         notes=f"BioPython 사용 불가 ({type(exc).__name__}: {exc}) — "
                               f"`pip install biopython` 후 재실행. "
                               f"무-날조 정책상 물성 수치를 생성하지 않습니다."))
        return 1

    results, errors = [], []
    for rec in records:
        seq = clean_seq(rec.get("sequence", ""))
        if not valid_protein_seq(seq, min_len=20):
            errors.append(f"{rec.get('id')}: 유효한 단백질 서열 아님 (len={len(seq)})")
            continue
        try:
            results.append(analyze(seq, rec.get("id", "input")))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{rec.get('id')}: 계산 실패 {type(exc).__name__}: {exc}")

    sane = all(
        0 < r["properties"]["molecular_weight_Da"] < 1_000_000
        and 0 < r["properties"]["isoelectric_point_pI"] < 14
        and -5 < r["properties"]["gravy_kyte_doolittle"] < 5
        for r in results
    )
    checks = [
        ("분석된 사슬 ≥ 1", len(results) >= 1),
        ("모든 물성값 물리적 범위 내 (0<MW, 0<pI<14, -5<GRAVY<5)", sane and len(results) >= 1),
        ("입력 서열 전부 처리됨 (스킵 0건)", len(errors) == 0),
    ]
    n_liab = sum(r["summary"]["n_liability_hits"] for r in results)
    n_cdr = sum(r["summary"]["n_in_cdr"] for r in results)
    notes = (f"사슬 {len(results)}건 분석. liability 규칙 히트 총 {n_liab}건 (CDR 내부 {n_cdr}건). "
             f"물성은 BioPython ProtParam 실계산, liability 는 정규식 규칙 기반 플래그이며 "
             f"예측 모델·실험값이 아니다. instability index 는 Guruprasad 1990 기준(>40 불안정). "
             + (f"처리 실패: {errors}" if errors else ""))
    emit(make_result(results, "BioPython ProtParam + rule-based liability motif scan",
                     query, checks, notes=notes))
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result([], "developability (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
