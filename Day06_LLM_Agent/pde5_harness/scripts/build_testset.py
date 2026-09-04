#!/usr/bin/env python3
"""held-out 시험 세트 10건 — 기존 30건과 겹치지 않게 뽑는다.

가중치를 적합한 화합물로 그 가중치를 시험하면 순환 논리다. 훈련 30건을 그대로 두고
별도 10건을 뽑아 독립 검증에 쓴다. 층화는 훈련과 같은 규칙을 따르되 이미 쓴 ID 는 제외한다.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402
from build_dataset import STRATA, TARGET, evenly  # noqa: E402

ROOT = HERE.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-stratum", type=int, default=4)   # 4+3+3 = 10
    ap.add_argument("--train", default=str(ROOT / "sample_run" / "dataset30.json"))
    ap.add_argument("--out", default=str(ROOT / "sample_run" / "testset10.json"))
    a = ap.parse_args()

    used = {r["molecule_chembl_id"]
            for r in json.loads(Path(a.train).read_text())["result"]}
    from chembl_webresource_client.new_client import new_client
    from rdkit import Chem
    import rdkit.RDLogger as RDLogger
    RDLogger.DisableLog("rdApp.*")
    act = new_client.activity

    want = [a.per_stratum, 3, 3]
    picked, log = [], []
    for (name, lo, hi), k in zip(STRATA, want):
        q = act.filter(target_chembl_id=TARGET, standard_type__in=["IC50", "Ki"],
                       standard_relation="=", pchembl_value__isnull=False)
        if lo is not None:
            q = q.filter(pchembl_value__gte=lo)
        if hi is not None:
            q = q.filter(pchembl_value__lt=hi)
        rows, seen = [], set()
        for r in q:
            cid = r.get("molecule_chembl_id"); smi = r.get("canonical_smiles")
            if not cid or not smi or cid in seen or cid in used:
                continue
            if Chem.MolFromSmiles(smi) is None:
                continue
            seen.add(cid)
            rows.append({"molecule_chembl_id": cid, "canonical_smiles": smi,
                         "standard_type": r.get("standard_type"),
                         "standard_value": float(r["standard_value"]),
                         "standard_units": r.get("standard_units"),
                         "pchembl_value": float(r["pchembl_value"]),
                         "stratum": name})
            if len(rows) >= 400:
                break
        rows.sort(key=lambda x: x["pchembl_value"])
        # 훈련과 같은 위치를 집지 않도록 반 칸 어긋나게 집는다
        take = evenly(rows[1:], k) if len(rows) > k + 1 else evenly(rows, k)
        picked += take
        log.append({"stratum": name, "candidates": len(rows), "taken": len(take)})
        print(f"  {name:7s} 후보 {len(rows):4d}(훈련 제외) → 선택 {len(take)}  "
              f"pChEMBL {take[0]['pchembl_value']:.2f}–{take[-1]['pchembl_value']:.2f}",
              file=sys.stderr)

    pv = [r["pchembl_value"] for r in picked]
    ids = {r["molecule_chembl_id"] for r in picked}
    checks = [
        ("훈련 세트와 겹치지 않음", not (ids & used)),
        ("중복 없음", len(ids) == len(picked)),
        ("역가 범위 >= 2.5 로그", (max(pv) - min(pv)) >= 2.5),
        ("목표 건수 확보", len(picked) == sum(want)),
    ]
    env = make_result(picked, "ChEMBL webresource_client",
                      f"target={TARGET}; 훈련 30건 제외; 층화 {want}",
                      checks,
                      notes=(f"held-out n={len(picked)}, pChEMBL {min(pv):.2f}–{max(pv):.2f}. "
                             "가중치 적합에 쓰이지 않는다."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n{a.out}  n={len(picked)}  pChEMBL {min(pv):.2f}–{max(pv):.2f}  "
          f"게이트 {'PASS' if env['verification']['passed'] else 'FAIL'}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
