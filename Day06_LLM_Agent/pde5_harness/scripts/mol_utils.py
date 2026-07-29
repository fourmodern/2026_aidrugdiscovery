"""mol_utils.py — 경량 구조 유틸 (ddtut1 스타일, 무거운 의존성 없이).

RCSB에서 PDB 파일을 내려받고 간단한 구조 메타정보(HEADER/TITLE/해상도/HETATM 리간드)를
파싱한다. requests 또는 표준 urllib만 사용하며 RDKit/MDAnalysis/OpenBabel을 요구하지 않는다.
도킹 파이프라인의 전처리 보조용이며, 네트워크 실패 시 graceful하게 처리한다(무-날조).

PDE5A 대표 구조(공개 PDB ID):
  1UDT — PDE5A 촉매 도메인 + sildenafil 공결정
  1XP0 — PDE5A + vardenafil
  1XOZ — PDE5A + sildenafil (대체)
  1UDU — PDE5A apo/촉매 도메인
  2H42 — PDE5A + 억제제 공결정

사용:
  python scripts/mol_utils.py 1UDT          # 다운로드 + 구조 정보(표준 봉투 JSON)
  from mol_utils import fetch_pdb, get_structure_info
"""
from __future__ import annotations
import os, sys, json
try:
    import requests
except Exception:
    requests = None
import urllib.request

# PDE5A 대표 공개 PDB 코드(참조용 — 실재 코드만 기재).
PDE5A_PDB = {
    "1UDT": "PDE5A catalytic domain in complex with sildenafil",
    "1XP0": "PDE5A in complex with vardenafil",
    "1XOZ": "PDE5A in complex with sildenafil",
    "1UDU": "PDE5A catalytic domain",
    "2H42": "PDE5A in complex with an inhibitor",
}

_WATER = {"HOH", "WAT", "SOL", "DOD"}


def fetch_pdb(code: str, save_dir: str = "outputs/structures"):
    """RCSB에서 PDB 파일 다운로드. 반환: 로컬 경로(성공) 또는 None(실패).

    무-날조: 다운로드 실패 시 파일을 만들지 않고 None을 반환한다.
    """
    code = code.strip().lower()
    os.makedirs(save_dir, exist_ok=True)
    out_path = os.path.join(save_dir, f"{code}.pdb")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path
    url = f"https://files.rcsb.org/download/{code}.pdb"
    try:
        if requests is not None:
            r = requests.get(url, timeout=30)
            if r.status_code != 200 or not r.text.strip():
                return None
            with open(out_path, "w") as fh:
                fh.write(r.text)
        else:
            urllib.request.urlretrieve(url, out_path)
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            return out_path
        return None
    except Exception:
        return None


def _parse_pdb_text(text: str) -> dict:
    """PDB 텍스트에서 간단한 메타정보 파싱(외부 의존성 없이)."""
    info = {"header": "", "title": "", "resolution_A": None,
            "n_atoms": 0, "n_hetatm": 0, "chains": [], "hetero_ligands": []}
    title_parts, chains, hets = [], set(), {}
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "HEADER":
            info["header"] = line[10:50].strip()
        elif rec == "TITLE":
            title_parts.append(line[10:80].strip())
        elif rec == "REMARK" and "RESOLUTION." in line and "ANGSTROM" in line:
            for tok in line.replace(":", " ").split():
                try:
                    info["resolution_A"] = float(tok)
                    break
                except ValueError:
                    continue
        elif rec == "ATOM":
            info["n_atoms"] += 1
            ch = line[21:22].strip()
            if ch:
                chains.add(ch)
        elif rec == "HETATM":
            info["n_hetatm"] += 1
            resn = line[17:20].strip()
            if resn and resn not in _WATER:
                hets[resn] = hets.get(resn, 0) + 1
    info["title"] = " ".join(title_parts).strip()
    info["chains"] = sorted(chains)
    # 리간드 후보(물 제외 HETATM 그룹) — 원자 3개 이상만.
    info["hetero_ligands"] = sorted([k for k, v in hets.items() if v >= 3])
    return info


def get_structure_info(code: str, save_dir: str = "outputs/structures"):
    """PDB 다운로드 + 파싱된 구조 정보 반환(dict) 또는 None(실패)."""
    path = fetch_pdb(code, save_dir=save_dir)
    if path is None:
        return None
    try:
        with open(path) as fh:
            text = fh.read()
    except Exception:
        return None
    info = _parse_pdb_text(text)
    info["pdb_code"] = code.upper()
    info["path"] = path
    return info


def main():
    from verify import make_result
    code = sys.argv[1] if len(sys.argv) > 1 else "1UDT"
    known = code.upper() in PDE5A_PDB
    info = get_structure_info(code)
    url = f"https://files.rcsb.org/download/{code.lower()}.pdb"
    if info is None:
        result = {
            "pdb_code": code.upper(),
            "downloaded": False,
            "known_pde5a_reference": PDE5A_PDB.get(code.upper()),
            "note": "구조 다운로드 실패(네트워크·코드 확인 필요). 파일 미생성 — 날조 금지.",
            "pde5a_reference_codes": PDE5A_PDB,
        }
        checks = [
            ("PDB 다운로드 성공", False),
            ("PDE5A 참조 코드 목록 제공", len(PDE5A_PDB) > 0),
        ]
        env = make_result(result, "RCSB PDB", url, checks,
                          notes="다운로드 실패 시에도 참조 코드는 제공하되 구조 수치는 만들지 않음.")
    else:
        result = {
            "pdb_code": info["pdb_code"],
            "downloaded": True,
            "path": info["path"],
            "header": info["header"],
            "title": info["title"],
            "resolution_A": info["resolution_A"],
            "n_atoms": info["n_atoms"],
            "n_hetatm": info["n_hetatm"],
            "chains": info["chains"],
            "hetero_ligands": info["hetero_ligands"],
            "known_pde5a_reference": PDE5A_PDB.get(info["pdb_code"]),
        }
        checks = [
            ("PDB 다운로드 성공", True),
            ("원자 레코드 파싱됨", info["n_atoms"] > 0),
            ("체인 식별됨", len(info["chains"]) > 0),
            ("PDE5A 알려진 참조 코드 여부(참고)", known),
        ]
        env = make_result(result, "RCSB PDB (files.rcsb.org)", url, checks,
                          notes=f"{info['pdb_code']} 구조 정보 파싱(원자 {info['n_atoms']}, "
                                f"체인 {info['chains']}, 리간드 {info['hetero_ligands']}). 도킹 전처리 보조.")
    print(json.dumps(env, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
