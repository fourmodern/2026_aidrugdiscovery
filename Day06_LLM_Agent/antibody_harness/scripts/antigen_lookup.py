"""항원(표적) 조사 — UniProt REST 실조회. 무-날조: 실제 REST 응답만 사용.

기본 표적: HER2 / ERBB2 (UniProt P04626).

사용:
    python scripts/antigen_lookup.py                # 기본 P04626 (HER2)
    python scripts/antigen_lookup.py P04626         # accession 명시
    python scripts/antigen_lookup.py P01375         # 다른 표적(TNF-alpha 등)

출력: 표준 결과 봉투(JSON) — {result, provenance, verification}.

무-날조 정책:
- 네트워크/API 실패 시 서열·길이·기능을 **지어내지 않는다**. result=null,
  verification.passed=false 로 정직하게 보고한다.
- 오프라인 폴백은 제공하지 않는다 (UniProt 전체 레코드를 상수로 박는 것은
  캐시가 아니라 사실상 하드코딩이며, 다른 accession 요청에 대응할 수 없다).
"""
from __future__ import annotations

import json
import sys

try:
    import requests
except Exception:  # noqa: BLE001
    requests = None

from verify import emit, make_result, valid_protein_seq, valid_uniprot_acc

DEFAULT_ACC = "P04626"  # ERBB2 / HER2_HUMAN
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"


def _first(seq, default=None):
    return seq[0] if seq else default


def fetch(acc: str, timeout: int = 30):
    """UniProt REST 조회. (data|None, url, error|None) 반환."""
    url = UNIPROT_URL.format(acc=acc)
    if requests is None:
        return None, url, "requests 미설치 — pip install requests"
    try:
        r = requests.get(url, timeout=timeout, headers={"Accept": "application/json"})
    except Exception as exc:  # noqa: BLE001
        return None, url, f"네트워크 오류: {type(exc).__name__}: {exc}"
    if r.status_code != 200:
        return None, url, f"HTTP {r.status_code}"
    try:
        d = r.json()
    except Exception as exc:  # noqa: BLE001
        return None, url, f"JSON 파싱 실패: {exc}"

    desc = d.get("proteinDescription", {})
    protein = (desc.get("recommendedName", {}).get("fullName", {}).get("value")
               or _first([s.get("fullName", {}).get("value") for s in desc.get("submissionNames", [])], ""))
    alt_names = [a.get("fullName", {}).get("value", "") for a in desc.get("alternativeName", [])]
    short_names = [s.get("value", "") for s in
                   desc.get("recommendedName", {}).get("shortNames", [])]

    genes = d.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
    gene_synonyms = ([s.get("value", "") for s in genes[0].get("synonyms", [])]
                     if genes else [])

    function, subcellular = "", []
    for c in d.get("comments", []):
        ct = c.get("commentType")
        if ct == "FUNCTION" and not function:
            texts = c.get("texts", [])
            if texts:
                function = texts[0].get("value", "")
        elif ct == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations", []):
                v = loc.get("location", {}).get("value")
                if v:
                    subcellular.append(v)

    # 세포외 도메인(항체 에피토프가 위치하는 영역) 특징 추출
    topology = []
    for ft in d.get("features", []):
        if ft.get("type") in ("Topological domain", "Transmembrane", "Signal", "Domain"):
            loc = ft.get("location", {})
            topology.append({
                "type": ft.get("type"),
                "description": ft.get("description", ""),
                "start": loc.get("start", {}).get("value"),
                "end": loc.get("end", {}).get("value"),
            })

    seq = d.get("sequence", {}).get("value", "")
    organism = d.get("organism", {})

    return {
        "accession": d.get("primaryAccession", acc),
        "uniprot_id": d.get("uniProtkbId", ""),
        "gene": gene,
        "gene_synonyms": gene_synonyms,
        "protein_name": protein,
        "short_names": short_names,
        "alternative_names": alt_names,
        "organism": organism.get("scientificName", ""),
        "taxon_id": organism.get("taxonId"),
        "length": d.get("sequence", {}).get("length"),
        "sequence": seq,
        "function": function,
        "subcellular_location": subcellular,
        "topology_features": topology[:20],
        "reviewed": d.get("entryType", "").lower().startswith("uniprotkb reviewed"),
    }, url, None


def main() -> int:
    acc = (sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ACC).strip().upper()

    if not valid_uniprot_acc(acc):
        env = make_result(None, "input validation", f"accession={acc}",
                          [("accession 형식 유효", False)],
                          notes=f"'{acc}' 는 UniProt accession 형식이 아닙니다. 예: P04626 (HER2).")
        emit(env)
        return 1

    data, url, err = fetch(acc)

    if data is None:
        env = make_result(
            None, "UniProt REST (FAILED)", url,
            [("UniProt 조회 성공", False)],
            notes=(f"UniProt 조회 실패: {err}. 무-날조 정책에 따라 항원 정보를 생성하지 않습니다. "
                   f"네트워크 확인 후 재시도하세요. 재시도 실패 시 이 단계를 '미확인'으로 플래그하고 "
                   f"보고서 한계 섹션에 명시하십시오."),
        )
        emit(env)
        return 1

    seq_ok = valid_protein_seq(data.get("sequence", ""), min_len=50)
    checks = [
        ("UniProt 조회 성공 (HTTP 200)", True),
        (f"accession == {acc}", data.get("accession") == acc),
        ("gene name 존재", bool(data.get("gene"))),
        ("protein name 존재", bool(data.get("protein_name"))),
        ("서열 존재·아미노산 알파벳 유효", seq_ok),
        ("length 값과 실제 서열 길이 일치",
         data.get("length") == len(data.get("sequence", ""))),
        ("기능(FUNCTION) 서술 존재", bool(data.get("function"))),
    ]

    ecd = [t for t in data.get("topology_features", [])
           if t.get("type") == "Topological domain" and "extracellular" in (t.get("description") or "").lower()]
    notes = (f"항원 {data.get('gene')} ({data.get('protein_name')}), "
             f"{data.get('organism')}, length={data.get('length')}. "
             f"세포외 도메인(항체 접근 가능 영역) 특징 {len(ecd)}건 검출. "
             f"모든 값은 UniProt REST 응답에서만 추출 (무-날조).")

    emit(make_result(data, "UniProt REST", url, checks, notes=notes))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001  — 어떤 예외에도 표준 봉투로 정직 보고
        emit(make_result(None, "antigen_lookup (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
