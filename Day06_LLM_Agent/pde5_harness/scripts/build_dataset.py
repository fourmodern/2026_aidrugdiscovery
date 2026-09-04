#!/usr/bin/env python3
"""도킹·회귀용 데이터셋 구축 — 역가 층화 30건.

10건(강한 것만)으로는 회귀의 역가 범위가 1.6 로그밖에 안 돼 기울기를 추정할 수 없다.
ChEMBL 에 실측 약한 화합물이 있으므로 세 층에서 고르게 뽑아 범위를 넓힌다.

  강함   pChEMBL >= 7      10건
  중간   6 <= pChEMBL < 7  10건
  약함   pChEMBL < 6       10건

선택 규칙은 결정적이다: 각 층에서 pChEMBL 로 정렬한 뒤 균등 간격으로 집는다.
무작위 추출이 아니므로 시드가 필요 없고, 같은 ChEMBL 릴리스에서 같은 결과가 나온다.
음성은 '가정된 비활성'이 아니라 **실측된 약한 활성**이다.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
TARGET = "CHEMBL1827"
STRATA = [("strong", 7.0, None), ("medium", 6.0, 7.0), ("weak", None, 6.0)]


def evenly(seq, k):
    """정렬된 목록에서 균등 간격으로 k 개. 앞뒤 끝을 포함한다."""
    n = len(seq)
    if n <= k:
        return list(seq)
    step = (n - 1) / (k - 1)
    return [seq[round(i * step)] for i in range(k)]


def fetch(per_stratum: int):
    from chembl_webresource_client.new_client import new_client
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    act = new_client.activity
    picked, log = [], []
    for name, lo, hi in STRATA:
        q = act.filter(target_chembl_id=TARGET, standard_type__in=["IC50", "Ki"],
                       standard_relation="=", pchembl_value__isnull=False)
        if lo is not None:
            q = q.filter(pchembl_value__gte=lo)
        if hi is not None:
            q = q.filter(pchembl_value__lt=hi)
        rows = []
        seen = set()
        for a in q:
            cid = a.get("molecule_chembl_id"); smi = a.get("canonical_smiles")
            if not cid or not smi or cid in seen:
                continue
            if Chem.MolFromSmiles(smi) is None:
                continue
            seen.add(cid)
            rows.append({"molecule_chembl_id": cid, "canonical_smiles": smi,
                         "standard_type": a.get("standard_type"),
                         "standard_value": float(a["standard_value"]),
                         "standard_units": a.get("standard_units"),
                         "pchembl_value": float(a["pchembl_value"]),
                         "stratum": name})
            if len(rows) >= 400:
                break
        rows.sort(key=lambda r: r["pchembl_value"])
        take = evenly(rows, per_stratum)
        picked += take
        log.append({"stratum": name, "candidates": len(rows), "taken": len(take),
                    "pchembl_range": [rows[0]["pchembl_value"], rows[-1]["pchembl_value"]] if rows else None})
        print(f"  {name:7s} 후보 {len(rows):4d} → 선택 {len(take):2d}  "
              f"pChEMBL {take[0]['pchembl_value']:.2f}–{take[-1]['pchembl_value']:.2f}", file=sys.stderr)
    return picked, log


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-stratum", type=int, default=10)
    ap.add_argument("--out", default=str(ROOT / "sample_run" / "dataset30.json"))
    a = ap.parse_args()
    picked, log = fetch(a.per_stratum)
    pv = [r["pchembl_value"] for r in picked]
    checks = [
        ("세 층 모두 채워짐", all(l["taken"] == a.per_stratum for l in log)),
        ("전 화합물 SMILES 파싱", True),
        ("역가 범위 >= 3 로그", (max(pv) - min(pv)) >= 3.0 if pv else False),
        ("중복 ChEMBL ID 없음", len({r["molecule_chembl_id"] for r in picked}) == len(picked)),
    ]
    env = make_result(picked, "ChEMBL webresource_client",
                      f"target={TARGET}; type in (IC50,Ki); relation='='; pchembl not null; "
                      f"층화 {[s[0] for s in STRATA]} 각 {a.per_stratum}건; 정렬 후 균등 간격 선택",
                      checks,
                      notes=(f"n={len(picked)}, pChEMBL {min(pv):.2f}–{max(pv):.2f}. "
                             "음성은 가정된 비활성이 아니라 실측된 약한 활성이다. "
                             f"층별 후보 수: {log}"))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n{a.out}  n={len(picked)}  pChEMBL {min(pv):.2f}–{max(pv):.2f}  "
          f"게이트 {'PASS' if env['verification']['passed'] else 'FAIL'}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
