"""알려진 항체 수집 — RCSB PDB 실조회. 무-날조: PDB가 반환한 실제 사슬 서열만 사용.

주어진 항원(UniProt accession)이 포함된 구조를 RCSB Search API 로 찾고,
각 entry 의 polymer entity 를 Data API 로 조회해 **실제 사슬 서열**을 가져온다.
중쇄/경쇄 판별은 `seq_utils.classify_chain` (보존 Ig 모티프 휴리스틱).

사용:
    python scripts/antibody_search.py                       # 기본 P04626(HER2), 최대 8 entry
    python scripts/antibody_search.py P04626 8
    python scripts/antibody_search.py P04626 8 --offline    # 오프라인 캐시 강제(교육용)

출력: 표준 결과 봉투(JSON).

무-날조 정책:
- 검색/조회 실패 → result=[] + verification.passed=false. 서열을 만들지 않는다.
- 오프라인 폴백은 **RCSB 에서 실제로 내려받은 1N8Z(트라스투주맙 Fab–HER2 ECD) 서열**만
  상수로 보관하며 provenance.source 에 "offline cache (PDB 1N8Z)" 로 명시한다.
"""
from __future__ import annotations

import json
import sys

try:
    import requests
except Exception:  # noqa: BLE001
    requests = None

from seq_utils import classify_chain, clean_seq
from verify import emit, make_result, valid_pdb_id, valid_protein_seq, valid_uniprot_acc

DEFAULT_ACC = "P04626"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}"
ENTRY_PAGE = "https://www.rcsb.org/structure/{pdb}"

# --- 오프라인 캐시 -----------------------------------------------------------
# 아래 두 서열은 2026-08-31 에 RCSB Data API
#   https://data.rcsb.org/rest/v1/core/polymer_entity/1N8Z/{1,2}
# 에서 실제로 내려받은 `entity_poly.pdbx_seq_one_letter_code_can` 값이다.
# (지어낸 서열이 아니며, 재현하려면 위 URL 을 직접 호출해 대조할 수 있다.)
OFFLINE_1N8Z = {
    "pdb_id": "1N8Z",
    "title": "Crystal structure of extracellular domain of human HER2 complexed with Herceptin Fab",
    "experimental_method": "X-ray",
    "resolution_A": 2.52,
    "antigen_accession": "P04626",
    "url": ENTRY_PAGE.format(pdb="1N8Z"),
    "chains": [
        {
            "entity_id": "1",
            "auth_asym_ids": ["A"],
            "description": "Herceptin Fab (antibody) - light chain",
            "sequence": (
                "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISS"
                "LQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNAL"
                "QSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
            ),
        },
        {
            "entity_id": "2",
            "auth_asym_ids": ["B"],
            "description": "Herceptin Fab (antibody) - heavy chain",
            "sequence": (
                "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKN"
                "TAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFP"
                "EPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEP"
            ),
        },
    ],
}


def _get(url: str, timeout: int = 30):
    if requests is None:
        return None, "requests 미설치"
    try:
        r = requests.get(url, timeout=timeout, headers={"Accept": "application/json"})
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"
    try:
        return r.json(), None
    except Exception as exc:  # noqa: BLE001
        return None, f"JSON 파싱 실패: {exc}"


def search_entries(acc: str, rows: int, timeout: int = 40):
    """항원 accession 이 포함되고 단백질 entity 가 2개 이상인 PDB entry 검색."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers.database_accession",
                    "operator": "exact_match", "value": acc}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers.database_name",
                    "operator": "exact_match", "value": "UniProt"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                    "operator": "greater_or_equal", "value": 2}},
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max(rows * 3, 25)},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    if requests is None:
        return None, "requests 미설치", SEARCH_URL
    try:
        r = requests.post(SEARCH_URL, json=query, timeout=timeout)
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}", SEARCH_URL
    if r.status_code == 204:
        return [], "검색 결과 0건", SEARCH_URL
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}", SEARCH_URL
    try:
        d = r.json()
    except Exception as exc:  # noqa: BLE001
        return None, f"JSON 파싱 실패: {exc}", SEARCH_URL
    ids = [x["identifier"] for x in d.get("result_set", [])]
    return ids, None, SEARCH_URL


def fetch_entry(pdb: str, acc: str):
    """entry + polymer entity 를 조회해 사슬 서열/분류를 채운 dict 반환."""
    entry, err = _get(ENTRY_URL.format(pdb=pdb))
    if entry is None:
        return None, f"{pdb} entry 조회 실패: {err}"

    info = entry.get("rcsb_entry_info", {})
    res = info.get("resolution_combined") or []
    ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])

    chains, antigen_present = [], False
    for eid in ids:
        ent, err = _get(ENTITY_URL.format(pdb=pdb, eid=eid))
        if ent is None:
            continue
        seq = clean_seq(ent.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", ""))
        if not seq:
            continue
        cid = ent.get("rcsb_polymer_entity_container_identifiers", {})
        refs = [r.get("database_accession") for r in
                cid.get("reference_sequence_identifiers", []) or []]
        desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description", "") or ""
        cls = classify_chain(seq)
        role = "antigen" if acc in refs else cls["chain_type"]
        if role == "antigen":
            antigen_present = True
        chains.append({
            "entity_id": eid,
            "auth_asym_ids": cid.get("auth_asym_ids", []),
            "description": desc,
            "uniprot_refs": [r for r in refs if r],
            "role": role,
            "chain_type": cls["chain_type"],
            "classification_evidence": cls["evidence"],
            "length": len(seq),
            "sequence": seq,
        })

    return {
        "pdb_id": pdb,
        "title": entry.get("struct", {}).get("title", ""),
        "experimental_method": info.get("experimental_method"),
        "resolution_A": res[0] if res else None,
        "deposit_date": entry.get("rcsb_accession_info", {}).get("deposit_date"),
        "antigen_accession": acc if antigen_present else None,
        "antigen_present": antigen_present,
        "url": ENTRY_PAGE.format(pdb=pdb),
        "data_api_url": ENTRY_URL.format(pdb=pdb),
        "chains": chains,
    }, None


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    force_offline = "--offline" in sys.argv
    acc = (args[0] if args else DEFAULT_ACC).strip().upper()
    limit = int(args[1]) if len(args) > 1 else 8

    if not valid_uniprot_acc(acc):
        emit(make_result([], "input validation", f"accession={acc}",
                         [("accession 형식 유효", False)],
                         notes=f"'{acc}' 는 UniProt accession 형식이 아닙니다."))
        return 1

    query_desc = (f"RCSB Search: UniProt accession={acc} AND "
                  f"polymer_entity_count_protein>=2; max_entries={limit}")
    attempts = []

    if force_offline:
        ids, err = [], "사용자가 --offline 지정"
    else:
        ids, err, _ = search_entries(acc, limit)
        attempts.append(f"RCSB Search API: {'ok, ' + str(len(ids)) + ' hits' if err is None else err}")

    entries, fetch_errors = [], []
    if ids:
        for pdb in ids:
            if len(entries) >= limit:
                break
            if not valid_pdb_id(pdb):
                continue
            e, ferr = fetch_entry(pdb, acc)
            if e is None:
                fetch_errors.append(ferr)
                continue
            has_ab = any(c["chain_type"] in ("heavy", "light", "scfv") for c in e["chains"])
            if e["antigen_present"] and has_ab:
                entries.append(e)

    # --- 온라인 성공 -------------------------------------------------------
    if entries:
        n_h = sum(1 for e in entries for c in e["chains"] if c["chain_type"] == "heavy")
        n_l = sum(1 for e in entries for c in e["chains"] if c["chain_type"] == "light")
        all_seqs_ok = all(valid_protein_seq(c["sequence"])
                          for e in entries for c in e["chains"])
        checks = [
            ("RCSB 검색 성공", True),
            ("항체-항원 복합체 entry ≥ 1", len(entries) >= 1),
            ("모든 PDB ID 형식 유효", all(valid_pdb_id(e["pdb_id"]) for e in entries)),
            ("모든 사슬 서열이 유효한 아미노산 서열", all_seqs_ok),
            ("중쇄 후보 ≥ 1", n_h >= 1),
            ("경쇄 후보 ≥ 1", n_l >= 1),
        ]
        notes = (f"{acc} 항원을 포함한 항체 복합체 {len(entries)} entry 수집 "
                 f"(중쇄 후보 {n_h}, 경쇄 후보 {n_l}). "
                 f"사슬 서열은 RCSB Data API `entity_poly.pdbx_seq_one_letter_code_can` 원값. "
                 f"중쇄/경쇄 판별은 보존 Ig 모티프 휴리스틱 (실험적 확인 아님). "
                 + (f"일부 entry 조회 실패: {fetch_errors[:3]}" if fetch_errors else ""))
        emit(make_result(entries, "RCSB PDB Search + Data REST API", query_desc, checks, notes=notes))
        return 0

    # --- 실패 → 오프라인 캐시 (실제 PDB 서열) ------------------------------
    if acc == "P04626":
        cached = json.loads(json.dumps(OFFLINE_1N8Z))  # deep copy
        for c in cached["chains"]:
            cls = classify_chain(c["sequence"])
            c["chain_type"] = cls["chain_type"]
            c["classification_evidence"] = cls["evidence"]
            c["role"] = cls["chain_type"]
            c["length"] = len(c["sequence"])
            c["uniprot_refs"] = []
        cached["antigen_present"] = False
        cached["_note"] = ("OFFLINE CACHE — 실시간 조회 아님. 아래 서열은 RCSB Data API "
                           "polymer_entity/1N8Z/{1,2} 에서 실제로 내려받은 값이며, "
                           "항원(HER2 ECD) 사슬은 캐시에 포함하지 않았습니다.")
        checks = [
            ("RCSB 실시간 검색 성공", False),
            ("오프라인 캐시 서열 유효", all(valid_protein_seq(c["sequence"]) for c in cached["chains"])),
        ]
        notes = (f"RCSB 실시간 조회 실패({err or fetch_errors or '결과 0건'}) → PDB 1N8Z 오프라인 캐시 사용. "
                 f"**실데이터이지만 실시간 조회 결과가 아님** — 보고서에 반드시 표기하고, "
                 f"네트워크 복구 후 재실행하십시오. verification.passed=false 유지 (게이트 통과 아님). "
                 f"시도 로그: {attempts}")
        emit(make_result([cached], "offline cache (PDB 1N8Z, RCSB Data API 2026-08-31 취득)",
                         query_desc, checks, notes=notes))
        return 1

    # --- 실패 + 캐시 없음 → 빈 결과 정직 보고 -------------------------------
    checks = [
        ("RCSB 검색 성공", err is None),
        ("항체-항원 복합체 entry ≥ 1", False),
    ]
    notes = (f"{acc} 에 대한 항체-항원 복합체를 찾지 못했습니다 "
             f"(검색 오류: {err}; entry 조회 오류: {fetch_errors[:3]}). "
             f"무-날조 정책에 따라 서열을 생성하지 않습니다. "
             f"다른 accession 을 시도하거나 네트워크를 확인하십시오. 시도 로그: {attempts}")
    emit(make_result([], "RCSB PDB Search + Data REST API (no usable hits)",
                     query_desc, checks, notes=notes))
    return 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result([], "antibody_search (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
