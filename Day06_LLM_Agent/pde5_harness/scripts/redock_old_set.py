#!/usr/bin/env python3
"""옛 30건을 새 프로토콜로 다시 도킹해 '화합물 집합'과 '프로토콜'을 가른다.

§3.5 는 상관 붕괴의 원인을 특정하지 못한다고 적었다. 후보가 둘 남아 있었기 때문이다 —
화합물 집합이 다른 것과 프로토콜이 다른 것(탐색 깊이 16→64, 자세 규칙, 반복측정 집계).
공통 화합물이 3건뿐이라 비교할 수 없었는데, **옛 화합물을 새 프로토콜로 다시 돌리면**
그 비교가 가능해진다. 도킹 30건이면 몇 분이다. 못 한 게 아니라 안 한 것이었다.

  같은 화합물, 옛 프로토콜 rho ≈ -0.538
  같은 화합물, 새 프로토콜 rho = ?
    비슷하면 → 원인은 화합물 집합 (옛 30건이 특이했다)
    0 에 가까우면 → 원인은 프로토콜
"""
from __future__ import annotations
import argparse, json, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result                                    # noqa: E402
from analyze_controlled import spearman, perm_p, fisher_ci        # noqa: E402
from dock_controlled import one, WORK2                            # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=64)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=15)
    ap.add_argument("--out", default=str(SR / "old_set_new_protocol.json"))
    a = ap.parse_args()

    ds = json.loads((SR / "dataset30.json").read_text())["result"]
    reg = json.loads((SR / "regression.json").read_text()); reg = reg.get("result", reg)
    old_score = {r["chembl_id"]: r["top_pose_score"] for r in reg["rows"]}
    items = [(r["molecule_chembl_id"], r["canonical_smiles"], float(r["pchembl_value"]))
             for r in ds]
    WORK2.mkdir(parents=True, exist_ok=True)

    print(f"옛 {len(items)}건을 새 프로토콜(exh {a.exhaustiveness})로 재도킹…", file=sys.stderr)
    jobs = [(cid, smi, a.exhaustiveness, a.seed, 4) for cid, smi, _ in items]
    got, errs = {}, []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(one, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            if "error" in r: errs.append(r)
            else: got[r["chembl_id"]] = r
            if i % 10 == 0: print(f"  {i}/{len(jobs)}", file=sys.stderr)

    pv = {c: p for c, _, p in items}
    ids = [c for c, _, _ in items if c in got and c in old_score]
    y = [pv[c] for c in ids]
    s_old = [old_score[c] for c in ids]
    s_new = [got[c]["top_pose_score"] for c in ids]
    r_old, r_new = spearman(y, s_old), spearman(y, s_new)
    rmsd = [got[c].get("top_pose_mcs_rmsd") for c in ids
            if got[c].get("top_pose_mcs_rmsd") is not None]

    same = abs(r_new) >= 0.6 * abs(r_old)
    out = {"n": len(ids), "n_failed": len(errs),
           "exhaustiveness_old": 16, "exhaustiveness_new": a.exhaustiveness,
           "rho_old_protocol": round(r_old, 3), "rho_new_protocol": round(r_new, 3),
           "ci95_new": fisher_ci(r_new, len(ids)), "perm_p_new": perm_p(y, s_new),
           "spearman_between_protocols": round(spearman(s_old, s_new), 3),
           "mean_score_shift": round(sum(s_new) / len(s_new) - sum(s_old) / len(s_old), 3),
           "pose_validity_new": {
               "n": len(rmsd), "median_rmsd": round(sorted(rmsd)[len(rmsd) // 2], 3),
               "frac_under_2A": round(sum(1 for x in rmsd if x < 2) / len(rmsd), 3)} if rmsd else None,
           "rows": [{"chembl_id": c, "pchembl_value": pv[c],
                     "score_old_protocol": old_score[c],
                     "score_new_protocol": got[c]["top_pose_score"],
                     "top_pose_mcs_rmsd": got[c].get("top_pose_mcs_rmsd")} for c in ids],
           "verdict": ("원인은 화합물 집합이다: 같은 30건은 새 프로토콜에서도 상관을 유지한다"
                       if same else
                       "원인에 프로토콜이 포함된다: 같은 30건인데 새 프로토콜에서는 상관이 사라진다")}
    checks = [("전 화합물 재도킹 성공", len(ids) >= 0.9 * len(items)),
              ("두 프로토콜 점수 모두 확보", all(r["score_new_protocol"] is not None
                                                 for r in out["rows"])),
              ("동일 화합물 비교", len(ids) >= 25)]
    env = make_result(out, f"smina 재도킹 (exhaustiveness {a.exhaustiveness})",
                      f"옛 데이터셋 n={len(items)}", checks, notes=out["verdict"])
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n같은 {len(ids)}건")
    print(f"  옛 프로토콜 (exh 16)  rho {r_old:+.3f}")
    print(f"  새 프로토콜 (exh {a.exhaustiveness})  rho {r_new:+.3f}  "
          f"CI{out['ci95_new']}  p={out['perm_p_new']}")
    print(f"  두 프로토콜 점수 상관 {out['spearman_between_protocols']:+.3f}  "
          f"평균 점수 이동 {out['mean_score_shift']:+.3f} kcal/mol")
    if out["pose_validity_new"]:
        print(f"  새 프로토콜 자세 타당도: 중앙값 {out['pose_validity_new']['median_rmsd']} Å, "
              f"2 Å 미만 {out['pose_validity_new']['frac_under_2A']:.1%}")
    print(f"\n판정: {out['verdict']}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
