"""docking.py — 옵션 도킹 단계. smina 있으면 실행, 없으면 graceful 스킵(무-날조).

이 단계는 파이프라인의 **선택적** 단계다. smina 바이너리와 준비된 수용체(PDBQT)·
리간드가 모두 있어야 실제 도킹을 수행한다. 하나라도 없으면 도킹을 실행하지 않고
"도킹 미실행(구조·smina 필요)"을 표준 봉투로 반환한다. 도킹 스코어를 만들어내지 않는다.

주의: 도킹 스코어는 **가설**이며 결합 친화도의 실측이 아니다(CLAUDE.md 하드 규칙).

사용:
  python scripts/docking.py --receptor rec.pdbqt --ligand lig.pdbqt \
      --center X Y Z --size SX SY SZ
출력: 표준 결과 봉투(JSON).
"""
from __future__ import annotations
import sys, os, json, shutil, subprocess
from verify import make_result


def _arg(name, n=1):
    if name in sys.argv:
        i = sys.argv.index(name)
        vals = sys.argv[i + 1:i + 1 + n]
        return vals if n > 1 else (vals[0] if vals else None)
    return None if n == 1 else None


def find_smina():
    """PATH에서 smina(또는 vina) 탐색. 없으면 None."""
    for exe in ("smina", "smina.static", "vina"):
        p = shutil.which(exe)
        if p:
            return p, exe
    return None, None


def run_smina(smina, receptor, ligand, center, size, out_path):
    cmd = [smina, "-r", receptor, "-l", ligand,
           "--center_x", center[0], "--center_y", center[1], "--center_z", center[2],
           "--size_x", size[0], "--size_y", size[1], "--size_z", size[2],
           "--out", out_path, "--num_modes", "9"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return proc


def parse_scores(stdout: str):
    """smina/vina 표 출력에서 affinity(kcal/mol) 파싱. 실패 시 빈 리스트."""
    scores = []
    for line in stdout.splitlines():
        toks = line.split()
        if len(toks) >= 2 and toks[0].isdigit():
            try:
                scores.append(float(toks[1]))
            except ValueError:
                continue
    return scores


def main():
    receptor = _arg("--receptor")
    ligand = _arg("--ligand")
    center = _arg("--center", 3)
    size = _arg("--size", 3)
    out_path = _arg("--out") or "outputs/docking/docked.pdbqt"

    smina, exe = find_smina()

    # 전제 조건 점검 — 하나라도 없으면 실행하지 않고 정직하게 스킵.
    missing = []
    if smina is None:
        missing.append("smina/vina 바이너리 (PATH에 없음)")
    if not receptor or not os.path.exists(receptor):
        missing.append("수용체 파일(--receptor, PDBQT)")
    if not ligand or not os.path.exists(ligand):
        missing.append("리간드 파일(--ligand, PDBQT)")
    if not center or len(center) != 3:
        missing.append("도킹 박스 중심(--center X Y Z)")
    if not size or len(size) != 3:
        missing.append("도킹 박스 크기(--size SX SY SZ)")

    if missing:
        result = {
            "docking_executed": False,
            "reason": "전제 조건 미충족 — 도킹 미실행(구조·smina 필요)",
            "missing": missing,
            "smina_found": smina,
            "hint": ("Dockerfile로 smina 설치 + mol_utils.py로 구조 확보 + 수용체/리간드 "
                     "PDBQT 준비(prepare_receptor/prepare_ligand) 후 재실행."),
            "caveat": "도킹 미실행. 스코어를 만들어내지 않음(무-날조). 도킹 스코어는 실행 시에도 가설임.",
        }
        checks = [
            ("도킹 미실행을 정직하게 보고", result["docking_executed"] is False),
            ("미충족 전제 조건 명시", len(missing) > 0),
            ("스코어 날조 없음", "docking_executed" in result and not result["docking_executed"]),
        ]
        env = make_result(result, "docking (smina) — skipped", "preconditions unmet", checks,
                          notes="옵션 단계. 전제 미충족으로 graceful 스킵. 수치 날조 금지.")
        print(json.dumps(env, ensure_ascii=False, indent=2))
        return

    # 전제 충족 — 실제 도킹 시도.
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    try:
        proc = run_smina(smina, receptor, ligand, center, size, out_path)
        scores = parse_scores(proc.stdout)
        ok = proc.returncode == 0 and len(scores) > 0
        result = {
            "docking_executed": ok,
            "engine": exe,
            "engine_path": smina,
            "receptor": receptor,
            "ligand": ligand,
            "center": center,
            "size": size,
            "out_path": out_path if ok else None,
            "affinities_kcal_per_mol": scores,
            "best_affinity_kcal_per_mol": min(scores) if scores else None,
            "returncode": proc.returncode,
            "caveat": "도킹 스코어는 가설(결합 친화도 실측 아님). 실검증은 어세이 필요.",
        }
        checks = [
            ("도킹 프로세스 정상 종료", proc.returncode == 0),
            ("스코어 파싱됨", len(scores) > 0),
            ("스코어는 도구 실행 결과(날조 아님)", ok),
        ]
        env = make_result(result, f"{exe} docking", " ".join(sys.argv[1:]), checks,
                          notes=f"도킹 {'성공' if ok else '실패'}. 스코어는 가설로 취급.")
    except Exception as e:
        result = {
            "docking_executed": False,
            "reason": f"도킹 실행 중 예외: {type(e).__name__}",
            "detail": str(e)[:200],
            "caveat": "도킹 실패. 스코어 날조 금지.",
        }
        checks = [("도킹 실행", False), ("실패를 정직하게 보고", True)]
        env = make_result(result, f"{exe} docking (error)", " ".join(sys.argv[1:]), checks,
                          notes="도킹 실행 실패 — 수치 만들지 않음.")
    print(json.dumps(env, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
