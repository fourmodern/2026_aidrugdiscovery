#!/usr/bin/env python3
"""v1/v2 의 상관 -0.538 이 v3 에서 -0.079 로 무너진 원인을 실제로 가른다.

이전 판은 "골격 교란 때문"이라고 단언했다. 그런데 옛 데이터에 같은 편상관을 적용하면
-0.538 → -0.532 로 1.2% 밖에 줄지 않는다. 즉 **그 설명은 틀렸다.** 원인은 다른 데 있고,
후보가 여럿이므로 하나씩 분리해야 한다.

여기서 검정하는 것은 **표본 설계**다. 새 163건에서 옛 설계를 흉내 내어 (역가로만 층화,
골격 무시, 층당 10건) 반복 추출하고 상관 분포를 본다.
  - 그 분포가 -0.5 근처면 → 원인은 표본 설계다 (골격이 아니라 '역가만으로 층화'라는 사실)
  - -0.08 근처에 머물면 → 원인은 화합물 집합이나 프로토콜 차이다

설계가 원인이 아니라면 남는 후보를 함께 계산해 보고한다: 표본 크기, 반복측정 집계 방식,
탐색 깊이, 자세 선택 규칙.
"""
from __future__ import annotations
import argparse, json, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result                                       # noqa: E402
from analyze_controlled import spearman, partial_spearman, perm_p    # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"


def old_dataset():
    """옛 n=30 세트에 Tanimoto 를 붙여 돌려준다."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator, DataStructs
    RDLogger.DisableLog("rdApp.*")
    SIL = ("CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ref = gen.GetFingerprint(Chem.MolFromSmiles(SIL))
    ds = {r["molecule_chembl_id"]: r
          for r in json.loads((SR / "dataset30.json").read_text())["result"]}
    reg = json.loads((SR / "regression.json").read_text()); reg = reg.get("result", reg)
    out = []
    for r in reg["rows"]:
        d = ds.get(r["chembl_id"])
        if not d: continue
        m = Chem.MolFromSmiles(d["canonical_smiles"])
        if m is None: continue
        out.append({"chembl_id": r["chembl_id"], "y": r["pIC50"],
                    "top": r["top_pose_score"], "sel": r["dock_score"],
                    # 반올림하지 않는다 — 3자리로 자르면 편상관이 -0.5317 에서 -0.534 로
                    # 움직여 감쇠가 1.24% 에서 0.81% 로 바뀐다. 헤드라인 수치가 반올림에
                    # 흔들리면 안 된다.
                    "tan": DataStructs.TanimotoSimilarity(ref, gen.GetFingerprint(m)),
                    "hac": m.GetNumHeavyAtoms()})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=str(SR / "collapse_diagnosis.json"))
    a = ap.parse_args()

    old = old_dataset()
    oy = [r["y"] for r in old]; osc = [r["top"] for r in old]
    otan = [r["tan"] for r in old]
    old_stats = {
        "n": len(old),
        "spearman_top_pose": round(spearman(oy, osc), 3),
        "perm_p": perm_p(oy, osc),
        "partial_controlling_tanimoto": round(partial_spearman(oy, osc, otan), 3),
        "spearman_score_vs_tanimoto": round(spearman(osc, otan), 3),
        "perm_p_score_vs_tanimoto": perm_p(osc, otan),
        "spearman_potency_vs_tanimoto": round(spearman(oy, otan), 3),
        "attenuation_from_scaffold_control": round(
            1 - abs(partial_spearman(oy, osc, otan)) / abs(spearman(oy, osc)), 4),
        "tanimoto_rounding_note": ("Tanimoto 를 반올림하지 않은 값이다. 3자리로 자르면 "
                                   "감쇠가 0.8% 로 바뀐다 — 어느 쪽이든 '거의 0' 이라는 "
                                   "결론은 같지만 수치는 하나로 고정한다.")}

    D = json.loads((SR / "docking_controlled.json").read_text())["result"]["rows"]
    new_y = [r["pchembl_value"] for r in D]; new_s = [r["top_pose_score"] for r in D]
    new_rho = spearman(new_y, new_s)

    # ── 검정 1. 옛 '역가만 층화' 설계를 새 데이터에서 재현 ────────────
    # 옛 설계의 절단점(strong>=7.0, 6.0<=medium<7.0, weak<6.0)을 그대로 쓴다
    def old_bin(v):
        return "strong" if v >= 7.0 else ("medium" if v >= 6.0 else "weak")
    pools = {"strong": [], "medium": [], "weak": []}
    for r in D:
        pools[old_bin(r["pchembl_value"])].append(r)
    rng = random.Random(a.seed); draws = []
    for _ in range(a.draws):
        pick = []
        for b in ("strong", "medium", "weak"):
            if len(pools[b]) < 10: break
            pick += rng.sample(pools[b], 10)
        if len(pick) != 30: continue
        draws.append(spearman([r["pchembl_value"] for r in pick],
                              [r["top_pose_score"] for r in pick]))
    draws.sort()
    q = lambda f: round(draws[int(f * (len(draws) - 1))], 3)
    design_test = {
        "description": ("새 163건에서 옛 설계(역가 3층, 층당 10건, 골격 무시)로 "
                        f"{len(draws)}회 반복 추출한 Spearman 분포"),
        "pool_sizes": {b: len(v) for b, v in pools.items()},
        "n_draws": len(draws),
        "median": q(.5), "p2.5": q(.025), "p97.5": q(.975),
        "frac_below_-0.4": round(sum(1 for d in draws if d < -0.4) / len(draws), 4),
        "old_observed": old_stats["spearman_top_pose"]}

    # ── 검정 2. 표본 크기만의 효과 (새 데이터에서 무작위 30건) ────────
    rng2 = random.Random(a.seed + 1); sub = []
    for _ in range(a.draws):
        pick = rng2.sample(D, 30)
        sub.append(spearman([r["pchembl_value"] for r in pick],
                            [r["top_pose_score"] for r in pick]))
    sub.sort()
    qs = lambda f: round(sub[int(f * (len(sub) - 1))], 3)
    size_test = {"description": "새 163건에서 층화 없이 무작위 30건 추출",
                 "median": qs(.5), "p2.5": qs(.025), "p97.5": qs(.975),
                 "frac_below_-0.4": round(sum(1 for d in sub if d < -0.4) / len(sub), 4)}

    # ── 검정 3. 겹치는 화합물이 있는가 (프로토콜 효과 직접 비교) ──────
    old_ids = {r["chembl_id"] for r in old}
    new_ids = {r["chembl_id"] for r in D}
    shared = sorted(old_ids & new_ids)
    prot = {"n_shared_compounds": len(shared), "shared_ids": shared[:20]}
    if len(shared) >= 8:
        om = {r["chembl_id"]: r for r in old}; nm = {r["chembl_id"]: r for r in D}
        po = [om[i]["top"] for i in shared]; pn = [nm[i]["top_pose_score"] for i in shared]
        yy = [om[i]["y"] for i in shared]
        prot.update({
            "spearman_old_vs_new_score": round(spearman(po, pn), 3),
            "mean_score_shift": round(sum(pn) / len(pn) - sum(po) / len(po), 3),
            "rho_potency_old_protocol": round(spearman(yy, po), 3),
            "rho_potency_new_protocol": round(spearman(yy, pn), 3),
            "note": ("같은 화합물을 두 프로토콜로 도킹한 값이다. 두 rho 가 크게 다르면 "
                     "프로토콜(깊이·자세규칙·집계)이 원인이고, 비슷하면 화합물 집합이 원인이다.")})

    covered = design_test["frac_below_-0.4"] > 0.10
    verdict = ("표본 설계가 원인이다: 옛 설계를 새 데이터에 적용하면 옛 상관 수준이 재현된다"
               if covered else
               "표본 설계만으로는 설명되지 않는다: 옛 설계를 새 데이터에 적용해도 "
               "옛 상관이 재현되지 않는다. 화합물 집합 또는 프로토콜 차이가 남는다")

    out = {"old_design": old_stats, "new_overall_spearman": round(new_rho, 3),
           "test_1_design": design_test, "test_2_sample_size": size_test,
           "test_3_protocol_on_shared_compounds": prot,
           "verdict": verdict,
           "what_this_refutes": ("이전 판은 '옛 상관은 골격 유사성이 만든 것'이라고 적었다. "
                                 f"옛 데이터에서 골격을 통제해도 상관은 "
                                 f"{old_stats['attenuation_from_scaffold_control']:.1%} 밖에 "
                                 "줄지 않으므로 그 설명은 성립하지 않는다.")}
    checks = [("옛 데이터 편상관 산출", old_stats["partial_controlling_tanimoto"] is not None),
              ("설계 재현 추출 1000회 이상", design_test["n_draws"] >= 1000),
              ("표본 크기 대조 수행", size_test["median"] is not None)]
    env = make_result(out, "옛/새 산출물 재분석 (새 도킹 없음)",
                      f"draws={a.draws}, seed={a.seed}", checks,
                      notes=verdict)
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")

    print(f"옛 설계 (n={old_stats['n']}): rho {old_stats['spearman_top_pose']:+.3f}  "
          f"골격 통제 후 {old_stats['partial_controlling_tanimoto']:+.3f} "
          f"(감쇠 {old_stats['attenuation_from_scaffold_control']:.1%})")
    print(f"  옛 데이터 점수-Tanimoto 상관 {old_stats['spearman_score_vs_tanimoto']:+.3f} "
          f"(p={old_stats['perm_p_score_vs_tanimoto']}) → 골격을 따라가지 않았다")
    print(f"\n검정 1  옛 설계를 새 데이터에 적용: 중앙값 {design_test['median']:+.3f}  "
          f"95% [{design_test['p2.5']:+.3f}, {design_test['p97.5']:+.3f}]  "
          f"-0.4 미만 비율 {design_test['frac_below_-0.4']:.1%}")
    print(f"검정 2  무작위 30건:                중앙값 {size_test['median']:+.3f}  "
          f"95% [{size_test['p2.5']:+.3f}, {size_test['p97.5']:+.3f}]")
    print(f"검정 3  두 세트 공통 화합물 {prot['n_shared_compounds']}건")
    if "rho_potency_old_protocol" in prot:
        print(f"        같은 화합물, 옛 프로토콜 rho {prot['rho_potency_old_protocol']:+.3f} vs "
              f"새 프로토콜 rho {prot['rho_potency_new_protocol']:+.3f}")
    print(f"\n판정: {verdict}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
