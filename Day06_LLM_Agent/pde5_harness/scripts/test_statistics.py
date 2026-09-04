#!/usr/bin/env python3
"""직접 짠 구현을 참조 구현과 대조한다 — 수식이나 파서가 틀리면 결론 전체가 무의미하다.

도킹 점수에는 동점이 많아 동점 처리가 특히 중요하다. 단순 정렬 순위를 쓰면 Spearman 값이
달라지고, 그 차이가 p 값의 유의성 판정을 뒤집을 수 있다. 따라서 동점을 일부러 만들어 검정한다.
"""
from __future__ import annotations
import json, math, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result                                   # noqa: E402
from analyze_controlled import (spearman, partial_spearman, partial_perm_p,  # noqa: E402
                                partial_spearman_multi, fisher_ci, auc)

ROOT = HERE.parent


def main() -> int:
    try:
        from scipy import stats
        import numpy as np
    except ImportError:
        print("scipy 없음 — 대조 불가"); return 1

    rng = random.Random(1)
    worst = {"spearman": 0.0, "partial": 0.0, "auc": 0.0}
    counts = {"spearman": 0, "partial": 0, "auc": 0}

    for _ in range(200):
        n = rng.randint(8, 60)
        x = [rng.gauss(0, 1) for _ in range(n)]
        y = [xi * .6 + rng.gauss(0, 1) for xi in x]
        for _ in range(n // 4):                       # 동점 주입
            i, j = rng.randrange(n), rng.randrange(n); y[i] = y[j]
        worst["spearman"] = max(worst["spearman"],
                                abs(spearman(x, y) - stats.spearmanr(x, y).statistic))
        counts["spearman"] += 1

    for _ in range(100):
        n = rng.randint(12, 60)
        z = [rng.gauss(0, 1) for _ in range(n)]
        x = [zi * .5 + rng.gauss(0, 1) for zi in z]
        y = [zi * .4 + xi * .3 + rng.gauss(0, 1) for zi, xi in zip(z, x)]
        rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
        def resid(a, b):
            s, i = np.polyfit(b, a, 1); return a - (s * b + i)
        ref = stats.pearsonr(resid(rx, rz), resid(ry, rz)).statistic
        worst["partial"] = max(worst["partial"], abs(partial_spearman(x, y, z) - ref))
        counts["partial"] += 1

    for _ in range(100):
        n = rng.randint(10, 50)
        lab = [rng.randint(0, 1) for _ in range(n)]
        if not (0 < sum(lab) < n): continue
        sc = [-(l * .8 + rng.gauss(0, 1)) for l in lab]
        pos = [-s for s, l in zip(sc, lab) if l]; neg = [-s for s, l in zip(sc, lab) if not l]
        ref = round(stats.mannwhitneyu(pos, neg, alternative="greater").statistic
                    / (len(pos) * len(neg)), 3)
        worst["auc"] = max(worst["auc"], abs(auc(sc, lab) - ref))
        counts["auc"] += 1

    # 다변량 편상관 — 여러 공변량 동시 통제
    for _ in range(150):
        nn = rng.randint(20, 120)
        z1 = [rng.gauss(0, 1) for _ in range(nn)]; z2 = [rng.gauss(0, 1) for _ in range(nn)]
        xc = [.5 * a1 + .3 * b1 + rng.gauss(0, 1) for a1, b1 in zip(z1, z2)]
        yc = [.4 * a1 - .2 * b1 + .25 * c1 + rng.gauss(0, 1) for a1, b1, c1 in zip(z1, z2, xc)]
        R = np.column_stack([np.ones(nn), stats.rankdata(z1), stats.rankdata(z2)])
        def _res(v):
            v = stats.rankdata(v); bb, *_ = np.linalg.lstsq(R, v, rcond=None)
            return v - R @ bb
        ref = stats.pearsonr(_res(xc), _res(yc)).statistic
        worst["partial_multi"] = max(worst.get("partial_multi", 0.0),
                                     abs(partial_spearman_multi(xc, yc, [z1, z2]) - ref))
        counts["partial_multi"] = counts.get("partial_multi", 0) + 1

    # SDF 파서 — 손으로 짠 다중 모델 파서를 RDKit 과 대조한다
    parser = {"files_checked": 0, "mismatches": 0, "skipped_in_flight": 0}
    try:
        import glob, time
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        sys.path.insert(0, str(HERE))
        from contact_concordance import sdf_models
        now = time.time()
        for f in sorted(glob.glob(str(ROOT / "sample_run" / "structures" /
                                      "work_controlled" / "*_poses.sdf"))):
            if now - Path(f).stat().st_mtime < 120:
                parser["skipped_in_flight"] += 1; continue
            try:
                mols = [m for m in Chem.SDMolSupplier(f, removeHs=True) if m]
            except OSError:
                continue
            if not mols:
                continue
            mine = sdf_models(Path(f))
            c = mols[0].GetConformer(); p0 = c.GetAtomPosition(0)
            ok = (len(mine) == len(mols)
                  and all(len(a) == m.GetNumHeavyAtoms() for a, m in zip(mine, mols))
                  and abs(mine[0][0][0] - p0.x) < 1e-3 and abs(mine[0][0][1] - p0.y) < 1e-3)
            parser["files_checked"] += 1
            if not ok:
                parser["mismatches"] += 1
            if parser["files_checked"] >= 25:
                break
    except Exception as exc:
        parser["error"] = f"{type(exc).__name__}: {exc}"

    # 편상관 순열 절차 — Freedman-Lane 대 단순 섞기가 결론을 바꾸는지 확인한다
    def _naive_p(x, y, zz, iters=4000, seed=42):
        obs = abs(partial_spearman(x, y, zz)); r = random.Random(seed)
        idx = list(range(len(y))); c = 0
        for _ in range(iters):
            r.shuffle(idx)
            if abs(partial_spearman(x, [y[i] for i in idx], zz)) >= obs - 1e-12:
                c += 1
        return (c + 1) / (iters + 1)
    rng3 = random.Random(0); perm_cmp = []
    for rho_zy in (0.1, 0.5, 0.8):
        nn = 120
        zc = [rng3.gauss(0, 1) for _ in range(nn)]
        yc = [rho_zy * v + rng3.gauss(0, 1) * (1 - rho_zy ** 2) ** .5 for v in zc]
        xc = [0.15 * v + rng3.gauss(0, 1) for v in yc]
        fl = partial_perm_p(xc, yc, zc, iters=4000)
        nv = _naive_p(xc, yc, zc)
        perm_cmp.append({"z_y_strength": rho_zy, "freedman_lane_p": round(fl, 4),
                         "naive_shuffle_p": round(nv, 4), "abs_diff": round(abs(fl - nv), 4)})

    z = 0.5 * math.log(1.5 / 0.5); se = 1 / math.sqrt(27)
    ci_ref = [round(math.tanh(z - 1.96 * se), 3), round(math.tanh(z + 1.96 * se), 3)]
    ci_mine = fisher_ci(0.5, 30)

    out = {"trials": counts, "max_abs_deviation": {k: float(f"{v:.3e}") for k, v in worst.items()},
           "reference": "scipy.stats (spearmanr / rankdata+pearsonr / mannwhitneyu)",
           "fisher_ci_check": {"input": "r=0.5, n=30", "mine": ci_mine, "analytic": ci_ref},
           "sdf_parser_vs_rdkit": parser,
           "partial_permutation_procedure": {
               "used": "Freedman-Lane (z 로 회귀한 잔차만 섞음)",
               "why": ("y 를 통째로 섞으면 y~z 관계까지 깨져 귀무가설이 달라진다. "
                       "본 연구는 올바른 절차를 쓰되, 단순 섞기와의 차이도 함께 확인했다."),
               "comparison": perm_cmp},
           "note": ("동점을 의도적으로 주입해 검정했다. 동점 보정을 빼면 도킹 점수처럼 "
                    "동점이 많은 데이터에서 값이 달라진다.")}
    checks = [
        ("Spearman scipy 일치 (동점 포함)", bool(worst["spearman"] < 1e-9)),
        ("편상관 순위잔차 방식 일치", bool(worst["partial"] < 1e-8)),
        ("다변량 편상관 최소제곱 잔차 일치", bool(worst.get("partial_multi", 1) < 1e-8)),
        ("AUC Mann-Whitney U 일치", bool(worst["auc"] < 1e-9)),
        ("Fisher CI 해석해 일치", bool(ci_mine == ci_ref)),
        ("SDF 파서 RDKit 일치", parser["files_checked"] > 0 and parser["mismatches"] == 0),
        ("편상관 순열 절차 두 방식 비교 수행", len(perm_cmp) == 3),
    ]
    env = make_result(out, "scipy · RDKit 대조 단위검정",
                      f"{sum(counts.values())} 회 무작위 시행 + SDF {parser['files_checked']}개 파일",
                      checks, notes=(f"최대 편차 {out['max_abs_deviation']}, "
                                     f"SDF 파서 불일치 {parser['mismatches']}건."))
    p = ROOT / "sample_run" / "statistics_validation.json"
    p.write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    for c in env["verification"]["checks"]:
        print(f"  {'OK ' if c['passed'] else 'NG '} {c['check']}")
    print(f"{p}  게이트 {'PASS' if env['verification']['passed'] else 'FAIL'}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
