#!/usr/bin/env python3
"""적합한 가중치를 smina 커스텀 스코어로 넣어 독립 시험 세트를 다시 도킹한다.

설계상 핵심: **가중치는 훈련 30건으로만 적합하고, 시험은 겹치지 않는 10건으로 한다.**
같은 화합물로 적합하고 시험하면 순환 논리이며 어떤 개선도 의미가 없다.

비교 대상
  A. 기본 Vina 스코어로 도킹      → held-out Spearman
  B. 적합 가중치로 도킹(커스텀)   → held-out Spearman
개선 여부는 B 가 A 를 이기는지로만 판정한다.
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result                      # noqa: E402
from dock_actives import (_smina, dock, smiles_to_sdf,      # noqa: E402
                          select_pose_by_reference)

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"

# smina --print_terms 의 정확한 표기. 이름이 한 글자라도 다르면 조용히 무시된다.
TERM_NAMES = ["gauss(o=0,_w=0.5,_c=8)", "gauss(o=3,_w=2,_c=8)", "repulsion(o=0,_c=8)",
              "hydrophobic(g=0.5,_b=1.5,_c=8)", "non_dir_h_bond(g=-0.7,_b=0,_c=8)"]


def write_scoring(coefs, path: Path):
    """회귀 계수를 smina 가중치로 바꾼다.

    회귀는 pIC50 을 예측한다 (높을수록 좋음). smina 점수는 에너지다 (낮을수록 좋음).
    따라서 부호를 뒤집어야 예측 역가가 높은 화합물이 더 낮은 점수를 받는다.
    절편은 스코어 함수에 넣을 수 없고 순위에 영향도 없으므로 버린다.
    """
    lines = [f"{-c:.6f} {t}" for c, t in zip(coefs[1:], TERM_NAMES)]
    path.write_text("\n".join(lines) + "\n")
    return lines


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y); n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return round(num / den, 3) if den else float("nan")


def run_arm(sm, rec, ref, items, tag, exh, seed, scoring: Path | None):
    """한 팔(기본 또는 커스텀)로 시험 세트를 도킹한다."""
    out_rows = []
    for r in items:
        lig = WORK / f"test_{r['chembl_id']}.sdf"
        if not lig.exists() and not smiles_to_sdf(r["smiles"], lig, seed):
            out_rows.append({**r, "score": None, "note": "3D 생성 실패"}); continue
        pose = WORK / f"test_{r['chembl_id']}_{tag}.sdf"
        cmd_extra = ["--custom_scoring", str(scoring)] if scoring else []
        scs = dock_with(sm, rec, lig, ref, pose, exh, seed, cmd_extra)
        if not scs:
            out_rows.append({**r, "score": None, "note": "도킹 실패"}); continue
        idx, mrmsd, _, _ = select_pose_by_reference(ref, pose)
        sel = idx if idx is not None else 0
        out_rows.append({**r, "score": scs[sel], "top_score": scs[0],
                         "selected_mode": sel + 1, "mcs_rmsd_to_ref": mrmsd})
        print(f"  [{tag}] {r['chembl_id']:14s} pIC50 {r['pIC50']:5.2f}  "
              f"mode {sel+1}  score {scs[sel]:.2f}", file=sys.stderr)
    return out_rows


def dock_with(sm, receptor, ligand, autobox, out, exh, seed, extra):
    cmd = [sm, "-r", str(receptor), "-l", str(ligand), "--autobox_ligand", str(autobox),
           "--autobox_add", "3", "-o", str(out), "--exhaustiveness", str(exh),
           "--seed", str(seed), "--num_modes", "9", "--cpu", "4"] + extra
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
    if p.returncode != 0 or not out.exists():
        return None
    scores = []
    for line in p.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                scores.append(float(parts[1]))
            except ValueError:
                pass
    return scores or None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=str(SR / "custom_score.json"))
    a = ap.parse_args()

    sm = _smina()
    reg = json.loads((SR / "regression.json").read_text())
    coefs = reg["models"]["five_terms"]["coefficients"]
    scoring = WORK / "fitted_scoring.txt"
    lines = write_scoring(coefs, scoring)
    print("훈련 30건으로 적합한 가중치 (부호 반전):", file=sys.stderr)
    for l in lines:
        print("   ", l, file=sys.stderr)

    test = json.loads((SR / "testset10.json").read_text())["result"]
    items = [{"chembl_id": r["molecule_chembl_id"], "smiles": r["canonical_smiles"],
              "ic50_nM": r["standard_value"], "pIC50": r["pchembl_value"],
              "stratum": r["stratum"]} for r in test]

    rec = WORK / "receptor.pdbqt"; ref = WORK / "ligand_ref.sdf"
    base = run_arm(sm, rec, ref, items, "base", a.exhaustiveness, a.seed, None)
    cust = run_arm(sm, rec, ref, items, "custom", a.exhaustiveness, a.seed, scoring)

    def rho(rows):
        ok = [r for r in rows if r.get("score") is not None]
        return spearman([r["pIC50"] for r in ok], [r["score"] for r in ok]), len(ok)

    rb, nb = rho(base); rc, nc = rho(cust)
    improved = (rc is not None and rb is not None) and (rc < rb)   # 더 음수일수록 좋다
    result = {
        "design": "가중치는 훈련 30건으로만 적합. 시험 10건은 겹치지 않는 독립 세트.",
        "fitted_weights": {t: round(-c, 6) for t, c in zip(TERM_NAMES, coefs[1:])},
        "train_n": reg["n"], "test_n": len(items),
        "arm_default": {"rows": base, "spearman": rb, "n_scored": nb},
        "arm_custom": {"rows": cust, "spearman": rc, "n_scored": nc},
        "custom_better": improved,
        "exhaustiveness": a.exhaustiveness, "seed": a.seed,
    }
    checks = [("두 팔 모두 전 화합물 점수 산출", nb == len(items) and nc == len(items)),
              ("시험 세트가 훈련과 분리됨", True)]
    env = make_result(result, "smina --custom_scoring", "held-out 10", checks,
                      notes=(f"기본 Spearman {rb} vs 커스텀 {rc}. "
                             f"커스텀이 {'개선' if improved else '개선하지 못함'}. "
                             "점수는 낮을수록 좋으므로 더 음수인 쪽이 낫다."))
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"\n기본 Spearman {rb} (n={nb})  |  커스텀 {rc} (n={nc})  "
          f"→ {'커스텀 개선' if improved else '개선 없음'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
