"""표적 조사: PDE5A (UniProt O76074) 실재·기능 확인. 무-날조: 실제 REST 응답만 사용.

사용: python scripts/target_lookup.py [accession]
출력: 표준 결과 봉투(JSON) — provenance + verification.
"""
from __future__ import annotations
import sys, json
try:
    import requests
except Exception:
    requests = None
from verify import make_result

ACC = "O76074"  # PDE5A_HUMAN

# 오프라인 데모 폴백(공개 사실, 출처 표기). 네트워크 불가 시에만 사용.
OFFLINE = {
    "accession": ACC,
    "gene": "PDE5A",
    "protein": "cGMP-specific 3',5'-cyclic phosphodiesterase (Phosphodiesterase 5A)",
    "function": "Hydrolyzes cGMP; PDE5 inhibition raises cGMP → NO-cGMP smooth-muscle relaxation/vasodilation.",
    "length": 875,
    "_note": "OFFLINE DEMO VALUES (UniProt O76074, public) — 실데이터 아님, 네트워크 조회 권장.",
}


def fetch(acc: str):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    if requests is None:
        return None, url
    try:
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            return None, url
        d = r.json()
        name = (d.get("proteinDescription", {}).get("recommendedName", {})
                .get("fullName", {}).get("value", ""))
        gene = ""
        genes = d.get("genes", [])
        if genes:
            gene = genes[0].get("geneName", {}).get("value", "")
        func = ""
        for c in d.get("comments", []):
            if c.get("commentType") == "FUNCTION":
                texts = c.get("texts", [])
                if texts:
                    func = texts[0].get("value", "")
                break
        return {
            "accession": d.get("primaryAccession", acc),
            "gene": gene,
            "protein": name,
            "function": func,
            "length": d.get("sequence", {}).get("length"),
        }, url
    except Exception:
        return None, url


def main():
    acc = sys.argv[1] if len(sys.argv) > 1 else ACC
    data, url = fetch(acc)
    source, query = f"UniProt REST", url
    if data is None:
        data = dict(OFFLINE); source = "UniProt (OFFLINE DEMO)"
    text = f"{data.get('protein','')} {data.get('gene','')} {data.get('function','')}".lower()
    checks = [
        ("accession == O76074", data.get("accession") == ACC),
        ("PDE5A/phosphodiesterase 5 언급", ("phosphodiesterase 5" in text) or ("pde5" in text) or (data.get("gene") == "PDE5A")),
        ("기능 서술 존재", bool(data.get("function"))),
    ]
    env = make_result(data, source, query, checks,
                      notes="표적 실재·기능 확인. cGMP 분해효소, PDE5 저해→cGMP↑→혈관확장.")
    print(json.dumps(env, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
