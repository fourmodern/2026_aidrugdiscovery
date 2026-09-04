#!/usr/bin/env python3
"""PDE5 하네스 로컬 실행 헬퍼 (Docker 없이 .venv 모드).

수강생이 환경을 빠르게 점검하고, 원하면 파이프라인을 직접 순차 실행할 수 있게 한다.
무-날조 원칙은 각 scripts/*.py 안에서 강제된다(이 파일은 오케스트레이션만).

사용법:
    python run_harness.py check     # 환경/의존성/스크립트 점검 (rdkit, requests, 스크립트 존재)
    python run_harness.py run       # target → chembl → mol-properties → selectivity 순차 실행(표준 봉투 출력)
    python run_harness.py list      # 사용 가능한 스킬/스크립트 목록

권장 사용법은 Claude Code/Desktop 로 이 폴더를 열고 자연어로 지시하는 것:
    "SETUP.md, run.md(README), CLAUDE.md 읽고 harness 자율 실행해줘"
"""
from __future__ import annotations
import sys, os, json, subprocess, importlib.util

ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, "scripts")
PY = sys.executable


def _has(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


def cmd_check() -> int:
    print("== PDE5 harness 환경 점검 (.venv 로컬 모드) ==")
    ok = True
    # 1) Python
    print(f"[python] {sys.version.split()[0]}  ({PY})")
    # 2) 핵심 의존성
    for mod, need in [("rdkit", True), ("requests", True),
                      ("chembl_webresource_client", False)]:
        present = _has(mod)
        tag = "ok" if present else ("MISSING(필수)" if need else "missing(옵션)")
        print(f"[dep] {mod:28s} {tag}")
        if need and not present:
            ok = False
    # 3) 스크립트 존재
    for s in ["verify.py", "target_lookup.py", "chembl_actives.py",
              "mol_properties.py", "selectivity.py"]:
        p = os.path.join(SCRIPTS, s)
        exists = os.path.isfile(p)
        print(f"[script] {s:22s} {'ok' if exists else 'MISSING'}")
        ok = ok and exists
    # 4) rdkit 실계산 smoke test (무-날조: 실제 계산되는지)
    if _has("rdkit"):
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            m = Chem.MolFromSmiles("CCO")
            print(f"[rdkit] smoke test: ethanol MW={Descriptors.MolWt(m):.2f}  ok")
        except Exception as e:  # pragma: no cover
            print(f"[rdkit] smoke test FAILED: {e}")
            ok = False
    print("\n결과:", "PASS — 자율 실행 준비 완료" if ok else "FAIL — 위 MISSING 항목 설치 필요 (pip install -r requirements.txt)")
    return 0 if ok else 1


def _run(script: str, args=None, stdin_text=None) -> dict | None:
    args = args or []
    try:
        p = subprocess.run([PY, os.path.join(SCRIPTS, script), *args],
                           input=stdin_text, capture_output=True, text=True, timeout=120)
        sys.stderr.write(p.stderr)
        return json.loads(p.stdout)
    except Exception as e:
        print(f"[run_harness] {script} 실행 오류: {e}", file=sys.stderr)
        return None


def cmd_run() -> int:
    print("== PDE5 harness 순차 실행 (검증 게이트) ==\n")
    sys.path.insert(0, SCRIPTS)
    from verify import gate  # noqa

    # (a) target
    env = _run("target_lookup.py")
    if env is None or not gate(env, "target-lookup"):
        print("target-lookup 실패 → 중단"); return 1
    print(json.dumps(env["result"], ensure_ascii=False, indent=2))

    # (b) chembl actives (limit 10)
    env = _run("chembl_actives.py", ["10"])
    if env is None or not gate(env, "chembl-actives"):
        print("chembl-actives 실패 → 중단"); return 1
    actives_json = json.dumps(env)
    print(f"활성물질 {len(env['result'])}건 조회")

    # (c) mol-properties (stdin 파이프)
    env = _run("mol_properties.py", ["--stdin"], stdin_text=actives_json)
    if env is None or not gate(env, "mol-properties"):
        print("mol-properties 실패 → 중단"); return 1
    print(json.dumps(env["result"], ensure_ascii=False, indent=2))

    # (d) selectivity
    env = _run("selectivity.py")
    if env is not None:
        gate(env, "selectivity-check")
    print("\n순차 실행 완료. 보고서 작성은 Claude Code report-writer 스킬로 수행하세요.")
    return 0


def cmd_list() -> int:
    print("== 스킬 (.claude/skills/) ==")
    sk = os.path.join(ROOT, ".claude", "skills")
    if os.path.isdir(sk):
        for name in sorted(os.listdir(sk)):
            print("  -", name)
    print("\n== 스크립트 (scripts/) ==")
    for name in sorted(os.listdir(SCRIPTS)):
        if name.endswith(".py"):
            print("  -", name)
    return 0


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    return {"check": cmd_check, "run": cmd_run, "list": cmd_list}.get(cmd, cmd_check)()


if __name__ == "__main__":
    raise SystemExit(main())
