"""Humanness 프록시 — 인간 germline V 유전자와의 서열 동일성(%) 실계산.

무-날조 정책 (중요):
- **"휴먼성 점수"라는 발명된 지표를 만들지 않는다.** 계산 가능하고 정의가 명확한
  **germline identity (%)** 하나만 보고한다.
- Germline 참조 서열은 **UniProt REST 에서 실시간 조회한 human IGHV/IGKV/IGLV 항목**이다
  (예: IGHV3-23 = P01764). 조회 실패 시 값을 만들지 않고 passed=false 로 보고한다.
- 이 지표는 임상 면역원성(ADA 발생률)을 예측하지 않는다. 상관은 있으나 인과가 아니며,
  T-cell epitope 예측·MHC 결합·응집·투여경로 등 다른 요인이 면역원성을 좌우한다.

방법:
1. 입력 사슬에서 가변영역(V-domain)을 추정 (seq_utils.variable_domain).
2. UniProt 에서 human germline V 유전자(reviewed) 서열 집합을 조회 (결과 캐시: outputs/cache/).
3. BioPython PairwiseAligner (BLOSUM62, local) 로 정렬.
4. identity(%) = (정렬 구간의 동일 잔기 수) / (정렬 구간 길이) × 100.
   → **local alignment 기준**이므로 germline 의 signal peptide 와 항체의 CDR3/FR4 는
     정렬에서 제외되는 것이 정상이다. 이 정의를 반드시 함께 보고한다.

사용:
    python scripts/humanness.py "EVQLVESGG..."
    python scripts/antibody_search.py | python scripts/humanness.py --stdin
    python scripts/humanness.py --stdin --top 5
"""
from __future__ import annotations

import json
import os
import sys

try:
    import requests
except Exception:  # noqa: BLE001
    requests = None

from seq_utils import (classify_chain, clean_seq, read_sequences, split_scfv,
                       variable_domain)
from verify import emit, make_result, valid_protein_seq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(ROOT, "outputs", "cache")
CACHE_FILE = os.path.join(CACHE_DIR, "human_germline_v.json")

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FAMILIES = {"heavy": ["IGHV"], "light": ["IGKV", "IGLV"]}

METRIC_DEFINITION = (
    "germline_identity_percent = 100 × (동일 잔기 수) / (BLOSUM62 local alignment 길이). "
    "local alignment 이므로 germline signal peptide 와 항체 CDR3/FR4 는 정렬 밖에 있다. "
    "이 값은 면역원성(ADA) 예측치가 아니라 서열 유사도 지표다."
)


def fetch_germlines(timeout: int = 60):
    """UniProt 에서 human germline V 서열 조회 (reviewed, non-fragment). 캐시 사용."""
    if os.path.isfile(CACHE_FILE):
        try:
            with open(CACHE_FILE, "r", encoding="utf-8") as fh:
                cached = json.load(fh)
            if cached.get("entries"):
                cached["source"] = cached.get("source", "UniProt REST") + " [local cache]"
                return cached, None
        except Exception:  # noqa: BLE001
            pass

    if requests is None:
        return None, "requests 미설치"

    entries, urls = [], []
    for fam in ["IGHV", "IGKV", "IGLV"]:
        params = {
            "query": f"(gene:{fam}*) AND (organism_id:9606) AND (reviewed:true) AND (fragment:false)",
            "fields": "accession,id,gene_names,protein_name,sequence",
            "format": "json",
            "size": 500,
        }
        try:
            r = requests.get(UNIPROT_SEARCH, params=params, timeout=timeout)
        except Exception as exc:  # noqa: BLE001
            return None, f"{fam} 조회 네트워크 오류: {type(exc).__name__}: {exc}"
        if r.status_code != 200:
            return None, f"{fam} 조회 HTTP {r.status_code}"
        urls.append(r.url)
        for e in r.json().get("results", []):
            genes = e.get("genes", [])
            gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
            seq = e.get("sequence", {}).get("value", "")
            if not gene.upper().startswith(fam) or not seq:
                continue
            entries.append({
                "accession": e.get("primaryAccession"),
                "uniprot_id": e.get("uniProtkbId"),
                "gene": gene,
                "family": fam,
                "sequence": seq,
                "length": len(seq),
            })

    if not entries:
        return None, "UniProt germline 조회 결과 0건"

    payload = {"entries": entries, "source": "UniProt REST", "urls": urls}
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(CACHE_FILE, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False)
    except Exception:  # noqa: BLE001
        pass
    return payload, None


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices

    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11.0
    a.extend_gap_score = -1.0
    a.mode = "local"
    return a


def identity(aligner, query: str, target: str):
    """BLOSUM62 local alignment identity(%). (percent, aln_len, n_identical, score)."""
    try:
        aln = aligner.align(query, target)[0]
    except Exception:  # noqa: BLE001
        return None
    qa, ta = str(aln[0]), str(aln[1])
    n = sum(1 for x, y in zip(qa, ta) if x == y and x != "-")
    L = len(qa)
    if L == 0:
        return None
    return round(100.0 * n / L, 2), L, n, float(aln.score)


def analyze(seq: str, seq_id: str, germlines: list, top: int) -> dict:
    s = clean_seq(seq)
    cls = classify_chain(s)
    vd = variable_domain(s)
    vseq = vd.get("sequence") or s
    fams = FAMILIES.get(cls["chain_type"], ["IGHV", "IGKV", "IGLV"])
    pool = [g for g in germlines if g["family"] in fams] or germlines

    aligner = _aligner()
    scored = []
    for g in pool:
        res = identity(aligner, vseq, g["sequence"])
        if res is None:
            continue
        pct, L, n, sc = res
        scored.append({
            "germline_gene": g["gene"],
            "germline_accession": g["accession"],
            "family": g["family"],
            "germline_identity_percent": pct,
            "alignment_length": L,
            "n_identical": n,
            "alignment_score_blosum62": sc,
        })
    scored.sort(key=lambda x: (-x["germline_identity_percent"], -x["alignment_length"]))
    best = scored[0] if scored else None

    return {
        "id": seq_id,
        "chain_type": cls["chain_type"],
        "variable_domain_used": {
            "start_0based": vd.get("start"), "end_0based": vd.get("end"),
            "length": len(vseq),
            "note": ("V-domain 경계 추정 실패 → 전체 사슬로 정렬"
                     if vd.get("sequence") is None else "seq_utils.variable_domain 근사 경계"),
        },
        "germline_pool_families": fams,
        "germline_pool_size": len(pool),
        "nearest_germline": best,
        "top_matches": scored[:top],
        "metric_definition": METRIC_DEFINITION,
        "limitations": [
            "germline identity 는 면역원성(ADA) 예측값이 아니다 — 상관 지표일 뿐이다.",
            "T-cell epitope / MHC-II 결합 / 응집 / 투여경로·용량은 이 지표에 반영되지 않는다.",
            "V-domain 경계는 보존 모티프 휴리스틱 추정치이며, 경계가 어긋나면 identity 도 변한다.",
            "germline pool 은 UniProt reviewed 항목에 한정되며 IMGT 전체 allele 집합이 아니다.",
            "local alignment 정의상 CDR3·FR4 는 정렬 밖 → 프레임워크 유사도에 가깝게 편향된다.",
        ],
    }


def _from_stdin() -> list:
    try:
        env = json.load(sys.stdin)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"[humanness] stdin JSON 파싱 실패: {exc}\n")
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
                out.append({"id": f"{r.get('pdb_id','?')}_{cid}", "sequence": c.get("sequence", "")})
        elif r.get("sequence"):
            out.append({"id": r.get("id", "input"), "sequence": r["sequence"]})
    return out


def main() -> int:
    argv = sys.argv[1:]
    top = 3
    if "--top" in argv:
        i = argv.index("--top")
        if i + 1 < len(argv):
            top = int(argv[i + 1])
            argv = argv[:i] + argv[i + 2:]

    if "--stdin" in sys.argv:
        records = _from_stdin()
        query = f"stdin envelope; top={top}"
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
        query = f"n_inputs={len(args)}; top={top}"

    try:
        from Bio import Align  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        emit(make_result([], "BioPython Align (MISSING)", query,
                         [("BioPython 설치", False)],
                         notes=f"BioPython 사용 불가 ({type(exc).__name__}: {exc}) — "
                               f"`pip install biopython` 후 재실행. 수치를 생성하지 않습니다."))
        return 1

    gl, err = fetch_germlines()
    if gl is None:
        emit(make_result([], "UniProt germline (FAILED)", query,
                         [("human germline V 서열 조회 성공", False)],
                         notes=f"germline 참조 서열 조회 실패: {err}. "
                               f"무-날조 정책상 germline identity 를 생성하지 않습니다. "
                               f"네트워크 확인 후 재시도하십시오."))
        return 1

    germlines = gl["entries"]
    results, errors = [], []
    for rec in records:
        s = clean_seq(rec.get("sequence", ""))
        if not valid_protein_seq(s, min_len=60):
            errors.append(f"{rec.get('id')}: 유효 서열 아님 (len={len(s)})")
            continue
        try:
            if classify_chain(s)["chain_type"] == "scfv":
                segs, order = split_scfv(s)
                for seg in segs or []:
                    d = analyze(seg["sequence"], f"{rec.get('id')}:{seg['chain_type']}",
                                germlines, top)
                    d["scfv_order"] = order
                    results.append(d)
            else:
                results.append(analyze(s, rec.get("id", "input"), germlines, top))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{rec.get('id')}: 계산 실패 {type(exc).__name__}: {exc}")

    n_hit = sum(1 for r in results if r.get("nearest_germline"))
    pcts = [r["nearest_germline"]["germline_identity_percent"]
            for r in results if r.get("nearest_germline")]
    checks = [
        (f"human germline V 참조 서열 로드 (n={len(germlines)})", len(germlines) >= 50),
        ("분석된 도메인 ≥ 1", len(results) >= 1),
        ("모든 도메인에 nearest germline 존재", n_hit == len(results) and len(results) >= 1),
        ("identity 값이 0-100 범위 내", all(0 <= p <= 100 for p in pcts) and bool(pcts)),
        ("입력 사슬 전부 처리됨 (스킵 0건)", len(errors) == 0),
    ]
    notes = (f"germline pool {len(germlines)}개(UniProt reviewed IGHV/IGKV/IGLV)와 "
             f"BLOSUM62 local alignment. 도메인 {len(results)}건 계산. "
             f"보고 지표는 germline identity(%) 하나뿐이며 '휴먼성 점수'를 발명하지 않는다. "
             f"{METRIC_DEFINITION} "
             + (f"처리 실패: {errors}" if errors else ""))
    emit(make_result(results, gl.get("source", "UniProt REST"), query, checks, notes=notes))
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result([], "humanness (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
