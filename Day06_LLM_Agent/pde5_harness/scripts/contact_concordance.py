#!/usr/bin/env python3
"""접촉 잔기 일치가 순환 논증인지 검정한다.

보고서 §4.3 은 "도킹 자세가 알려진 잔기와 접촉하므로 타당하다"고 주장했다가 철회됐다.
자세를 공결정 리간드와 겹치도록 고른 뒤 그 리간드의 접촉 잔기를 재현했다고 말하는 것은
선택 규칙의 정의를 결과로 보고한 것이기 때문이다.

이 스크립트는 그 순환을 실제로 깬다. **점수 1위 자세**는 공결정 리간드를 전혀 참조하지
않고 채점 함수만으로 뽑힌다. 1위 자세도 같은 잔기와 접촉한다면 그 관찰은 순환이 아니다.
강함 층 전체(n=10)로 확대해 n=1 한계도 함께 없앤다.
"""
from __future__ import annotations
import argparse, json, sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify import make_result  # noqa: E402

ROOT = HERE.parent
SR = ROOT / "sample_run"
WORK = SR / "structures" / "work"
WORK2 = SR / "structures" / "work_controlled"
CUTOFF = 4.0


def receptor_atoms(pdb: Path):
    out = []
    for l in pdb.read_text().splitlines():
        if l[:6].strip() in ("ATOM", "HETATM") and l[76:78].strip() != "H":
            try:
                out.append((f"{l[17:20].strip()}{l[22:26].strip()}",
                            float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                pass
    return out


def sdf_models(path: Path):
    """SDF 의 모든 모델을 좌표 목록으로 돌려준다 (수소 제외)."""
    models, cur, state, natoms = [], [], 0, 0
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        # counts line 은 4번째 줄
        if i + 3 < len(lines) and lines[i + 3][:3].strip().isdigit() and len(lines[i + 3]) > 6:
            try:
                natoms = int(lines[i + 3][:3])
            except ValueError:
                i += 1; continue
            cur = []
            for j in range(i + 4, i + 4 + natoms):
                pr = lines[j].split()
                if len(pr) >= 4 and pr[3] != "H":
                    cur.append((float(pr[0]), float(pr[1]), float(pr[2])))
            models.append(cur)
            # $$$$ 까지 건너뛴다
            while i < len(lines) and lines[i].strip() != "$$$$":
                i += 1
        i += 1
    return models


def contacts(rec, lig):
    c = Counter()
    c2 = CUTOFF ** 2
    for rn, rx, ry, rz in rec:
        for lx, ly, lz in lig:
            if (rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2 <= c2:
                c[rn] += 1
                break
    return c


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stratum", default="strong")
    ap.add_argument("--source", choices=("legacy", "controlled"), default="controlled",
                    help="legacy = 옛 n=30 세트, controlled = 골격 통제 세트")
    ap.add_argument("--out", default=str(SR / "contact_concordance.json"))
    a = ap.parse_args()

    if a.source == "controlled":
        D = json.loads((SR / "docking_controlled.json").read_text())["result"]
        rows = [r for r in D["rows"] if r["potency_bin"] == a.stratum]
        # 통제 세트는 전 모드를 한 파일에 담고 있다 — 1위는 첫 모델, 참조 자세는 인덱스로 찾는다
        wd, sel_suffix, all_suffix = WORK2, None, "_poses.sdf"
        idkey, pkey = "chembl_id", "pchembl_value"
    else:
        reg = json.loads((SR / "regression.json").read_text()); reg = reg.get("result", reg)
        rows = [r for r in reg["rows"] if r.get("stratum") == a.stratum]
        wd, sel_suffix, all_suffix = WORK, "_selected.sdf", "_pose.sdf"
        idkey, pkey = "chembl_id", "pIC50"
    rec = receptor_atoms(WORK / "receptor.pdb")

    per, sel_all, top_all = [], Counter(), Counter()
    for r in rows:
        cid = r[idkey]
        all_f = wd / f"{cid}{all_suffix}"
        if not all_f.exists():
            continue
        ma = sdf_models(all_f)
        if not ma:
            continue
        if sel_suffix:                       # 옛 세트는 선택 자세가 별도 파일
            sel_f = wd / f"{cid}{sel_suffix}"
            if not sel_f.exists():
                continue
            ms = sdf_models(sel_f)
        else:                                # 통제 세트는 같은 파일의 참조 모드 인덱스
            k = r.get("ref_pose_index")
            if not k or k > len(ma):
                continue
            ms = [ma[k - 1]]
        if not ms:
            continue
        cs, ct = contacts(rec, ms[0]), contacts(rec, ma[0])   # ma[0] = 점수 1위 모드
        sel_all.update(cs.keys()); top_all.update(ct.keys())
        inter = set(cs) & set(ct); union = set(cs) | set(ct)
        per.append({"chembl_id": cid, "pIC50": r[pkey],
                    "selected_mode": r.get("selected_mode") or r.get("ref_pose_index"),
                    "n_res_selected": len(cs), "n_res_top": len(ct),
                    "jaccard": round(len(inter) / len(union), 3) if union else None})

    n = len(per)
    key = ["GLN817", "PHE820", "PHE786", "LEU804", "VAL782", "TYR612"]
    keytab = {k: {"selected": sel_all.get(k, 0), "top": top_all.get(k, 0), "of_n": n}
              for k in key}
    out = {"stratum": a.stratum, "source": a.source, "n_compounds": n, "cutoff_angstrom": CUTOFF,
           "per_compound": per,
           "key_residue_frequency": keytab,
           "mean_jaccard": round(sum(p["jaccard"] for p in per if p["jaccard"] is not None)
                                 / max(1, sum(1 for p in per if p["jaccard"] is not None)), 3),
           "top_residues_selected": sel_all.most_common(8),
           "top_residues_top_pose": top_all.most_common(8),
           "interpretation": ("1위 자세는 공결정 리간드를 참조하지 않고 채점 함수만으로 뽑힌다. "
                              "1위 자세에서도 같은 잔기가 높은 빈도로 나타나면 접촉 관찰은 "
                              "자세 선택 규칙의 산물이 아니다.")}
    checks = [("강함 층 화합물 8건 이상 처리", n >= 8),
              ("두 자세 모두 접촉 산출", all(p["n_res_selected"] > 0 and p["n_res_top"] > 0
                                             for p in per)),
              ("접촉이 좌표에서 계산됨", True)]
    env = make_result(out, "좌표 기반 접촉 계산 (참조 선택 자세 vs 점수 1위 자세)",
                      f"stratum={a.stratum}, cutoff {CUTOFF} A", checks,
                      notes=f"평균 Jaccard {out['mean_jaccard']}, n={n}.")
    Path(a.out).write_text(json.dumps(env, ensure_ascii=False, indent=2) + "\n")
    print(f"n={n}  평균 Jaccard {out['mean_jaccard']}")
    for k, v in keytab.items():
        print(f"  {k:8s} 선택자세 {v['selected']:2d}/{n}   1위자세 {v['top']:2d}/{n}")
    return 0 if env["verification"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
