#!/usr/bin/env python3
"""항체 설계 하네스 로컬 실행 헬퍼 (Docker 없이 .venv 모드).

수강생이 환경을 빠르게 점검하고, 원하면 파이프라인을 직접 순차 실행할 수 있게 한다.
무-날조 원칙은 각 scripts/*.py 안에서 강제된다(이 파일은 오케스트레이션만).

사용법:
    python run_harness.py check              # 환경/의존성/스크립트 점검 + smoke test
    python run_harness.py list               # 사용 가능한 스킬/스크립트 목록
    python run_harness.py run [ACCESSION]    # 평가 파이프라인: antigen → antibody →
                                             # cdr → developability → humanness
                                             # (기본 P04626 = HER2, 결과는 outputs/ 에 저장)
    python run_harness.py design-check       # 설계 경로 A/B 의 GPU-free dry-run 점검

권장 사용법은 Claude Code 로 이 폴더를 열고 자연어로 지시하는 것:
    "CLAUDE.md 읽고 HER2 항체 하네스 자율 실행해줘"

로컬 .venv 안내 (Docker 불필요):
    python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt
    ./.venv/bin/python run_harness.py check

설계(경로 A/B)는 GPU 가 필요하다. `RUNPOD_가이드.md` 참조 — RunPod 제어는
`/home/hjpark/foundation_model_research/projects/_shared_infra/runpod_ctl.py` 를 재사용한다
(이 하네스는 RunPod 제어 코드를 자체 구현하지 않는다).
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, "scripts")
OUTPUTS = os.path.join(ROOT, "outputs")
PY = sys.executable

EVAL_STEPS = ["antigen_lookup.py", "antibody_search.py", "cdr_analysis.py",
              "developability.py", "humanness.py"]
DESIGN_STEPS = ["design_esmfold2.py", "design_rfantibody.py", "compare_designs.py"]
STEPS = EVAL_STEPS  # 하위호환
ALL_SCRIPTS = ["verify.py", "seq_utils.py"] + EVAL_STEPS + DESIGN_STEPS


def _has(mod: str) -> bool:
    try:
        return importlib.util.find_spec(mod) is not None
    except Exception:  # noqa: BLE001
        return False


# ---------------------------------------------------------------- check
def cmd_check() -> int:
    print("== 항체 설계 하네스 환경 점검 (.venv 로컬 모드) ==")
    ok = True

    print(f"[python] {sys.version.split()[0]}  ({PY})")
    if sys.version_info < (3, 8):
        print("[python] 경고: Python 3.8+ 권장")

    # 1) 의존성
    for mod, need, why in [
        ("requests", True, "UniProt/RCSB REST 조회"),
        ("Bio", True, "BioPython — ProtParam 물성 + PairwiseAligner germline identity"),
        ("abnumber", False, "정식 IMGT/Kabat 번호매김 (없으면 휴리스틱 폴백)"),
        ("anarci", False, "정식 IMGT/Kabat 번호매김 (없으면 휴리스틱 폴백)"),
    ]:
        present = _has(mod)
        tag = "ok" if present else ("MISSING(필수)" if need else "missing(옵션)")
        print(f"[dep] {mod:12s} {tag:16s} — {why}")
        if need and not present:
            ok = False

    # 2) 스크립트 존재
    for s in ALL_SCRIPTS:
        p = os.path.join(SCRIPTS, s)
        exists = os.path.isfile(p)
        print(f"[script] {s:22s} {'ok' if exists else 'MISSING'}")
        ok = ok and exists

    # 3) BioPython 실계산 smoke test (무-날조: 실제 계산되는지)
    if _has("Bio"):
        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            pa = ProteinAnalysis("EVQLVESGGGLVQPGGSLRLSCAAS")
            print(f"[biopython] ProtParam smoke: MW={pa.molecular_weight():.2f} "
                  f"pI={pa.isoelectric_point():.2f}  ok")
        except Exception as exc:  # noqa: BLE001
            print(f"[biopython] ProtParam smoke FAILED: {exc}")
            ok = False
        try:
            from Bio import Align
            from Bio.Align import substitution_matrices
            a = Align.PairwiseAligner()
            a.substitution_matrix = substitution_matrices.load("BLOSUM62")
            a.mode = "local"
            print(f"[biopython] Aligner smoke: BLOSUM62 local score="
                  f"{a.align('EVQLVESGGG', 'EVQLVESGGG')[0].score}  ok")
        except Exception as exc:  # noqa: BLE001
            print(f"[biopython] Aligner smoke FAILED: {exc}")
            ok = False

    # 4) CDR 휴리스틱 자기검증 (트라스투주맙 Kabat CDR 재현 — 문헌 공개값)
    sys.path.insert(0, SCRIPTS)
    try:
        from seq_utils import cdrs_heavy, cdrs_light, numbering_backend
        vh = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
              "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
        vl = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGT"
              "DFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
        h, lt = cdrs_heavy(vh), cdrs_light(vl)
        expect_h = {"CDR-H1": "DTYIH", "CDR-H2": "RIYPTNGYTRYADSVKG", "CDR-H3": "WGGDGFYAMDY"}
        expect_l = {"CDR-L1": "RASQDVNTAVA", "CDR-L2": "SASFLYS", "CDR-L3": "QQHYTTPPT"}
        good = all(h.get(k) == v for k, v in expect_h.items()) and \
               all(lt.get(k) == v for k, v in expect_l.items())
        print(f"[cdr] 휴리스틱 자기검증(트라스투주맙 Kabat CDR 6종 재현): "
              f"{'ok' if good else 'FAILED'}")
        print(f"[cdr] 번호매김 백엔드: {numbering_backend()}"
              + ("  (휴리스틱 근사 — IMGT 정확 번호 아님)"
                 if numbering_backend() == "heuristic" else ""))
        ok = ok and good
    except Exception as exc:  # noqa: BLE001
        print(f"[cdr] 자기검증 FAILED: {exc}")
        ok = False

    # 5) verify.py 계약 self-test
    try:
        from verify import gate, make_result
        env = make_result({"x": 1}, "self-test", "n/a", [("형식", True)])
        assert set(env) == {"result", "provenance", "verification"}
        assert gate(env, "self-test")
        bad = make_result(None, "self-test", "n/a", [("형식", False)])
        assert not gate(bad, "self-test")
        print("[verify] 표준 봉투 계약 self-test: ok")
    except Exception as exc:  # noqa: BLE001
        print(f"[verify] self-test FAILED: {exc}")
        ok = False

    # 6) 설계 경로 의존성 (옵션 — GPU 노드에서만 필요)
    print("\n-- 설계 경로 (GPU 필요, 로컬에서는 optional) --")
    for mod, why in [("torch", "GPU 점검 + ESMFold2"),
                     ("esm", "경로 A ESMFold2 inversion (pip install esm@git+...)")]:
        print(f"[design-dep] {mod:8s} {'ok' if _has(mod) else 'missing(옵션)'} — {why}")
    import shutil as _sh
    for exe in ("nvidia-smi", "rfdiffusion", "proteinmpnn", "rf2"):
        p = _sh.which(exe)
        print(f"[design-cli] {exe:12s} {'ok: ' + p if p else 'missing(옵션)'}")
    print("  → 설계는 GPU 노드에서 수행. 로컬에서는 "
          "`python run_harness.py design-check` 로 dry-run 만 검증.")

    # 7) outputs 디렉토리
    os.makedirs(OUTPUTS, exist_ok=True)
    print(f"\n[outputs] {OUTPUTS}  {'ok' if os.path.isdir(OUTPUTS) else 'MISSING'}")

    print("\n결과:", "PASS — 자율 실행 준비 완료" if ok
          else "FAIL — 위 MISSING 항목 설치 필요 (pip install -r requirements.txt)")
    return 0 if ok else 1


# ---------------------------------------------------------------- run
def _run(script: str, args=None, stdin_text=None):
    """스크립트를 실행해 표준 봉투(dict)를 반환. 실패 시 None."""
    args = args or []
    try:
        p = subprocess.run([PY, os.path.join(SCRIPTS, script), *args],
                           input=stdin_text, capture_output=True, text=True, timeout=300)
    except Exception as exc:  # noqa: BLE001
        print(f"[run_harness] {script} 실행 오류: {exc}", file=sys.stderr)
        return None
    if p.stderr.strip():
        sys.stderr.write(p.stderr)
    try:
        return json.loads(p.stdout)
    except Exception as exc:  # noqa: BLE001
        print(f"[run_harness] {script} 출력 JSON 파싱 실패: {exc}", file=sys.stderr)
        sys.stderr.write((p.stdout or "")[:500] + "\n")
        return None


def _save(name: str, env: dict) -> str:
    os.makedirs(OUTPUTS, exist_ok=True)
    path = os.path.join(OUTPUTS, name)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(env, fh, ensure_ascii=False, indent=2)
    return path


def cmd_run(acc: str = "P04626") -> int:
    print(f"== 항체 설계 하네스 순차 실행 (항원 {acc}, 검증 게이트) ==\n")
    sys.path.insert(0, SCRIPTS)
    from verify import gate  # noqa

    # (a) 항원 조사
    env = _run("antigen_lookup.py", [acc])
    if env is None or not gate(env, "antigen-lookup"):
        print("antigen-lookup 실패 → 중단 (무-날조: 항원 정보를 만들지 않음)")
        if env:
            _save("01_antigen.json", env)
        return 1
    _save("01_antigen.json", env)
    r = env["result"]
    print(f"항원: {r['gene']} / {r['protein_name']} / {r['organism']} / length={r['length']}")

    # (b) 알려진 항체 수집 (PDB)
    env = _run("antibody_search.py", [acc, "6"])
    if env is None:
        print("antibody-search 실행 실패 → 중단")
        return 1
    _save("02_antibodies.json", env)
    ab_json = json.dumps(env, ensure_ascii=False)
    n_entry = len(env.get("result") or [])
    n_chain = sum(1 for e in (env.get("result") or []) for c in e.get("chains", [])
                  if c.get("chain_type") in ("heavy", "light", "scfv"))
    print(f"항체 복합체 {n_entry} entry / 항체 사슬 {n_chain}건")
    if not gate(env, "antibody-search"):
        print("  → passed=false. 오프라인 캐시이거나 조회 실패입니다. "
              "결과를 '미확인'으로 플래그하고 보고서 한계에 명시하십시오.")
        if n_chain == 0:
            print("항체 사슬 0건 → 이후 단계 진행 불가. 중단.")
            return 1

    # (c) CDR 분석
    env = _run("cdr_analysis.py", ["--stdin"], stdin_text=ab_json)
    if env is None:
        print("cdr-analysis 실행 실패 → 중단")
        return 1
    _save("03_cdr.json", env)
    gate(env, "cdr-analysis")
    for d in (env.get("result") or [])[:4]:
        print(f"  {d['id']:12s} {d['chain_type']:6s} {d.get('cdrs')}")
    print(f"  method: {(env.get('result') or [{}])[0].get('method', 'n/a')}")

    # (d) Developability
    env = _run("developability.py", ["--stdin"], stdin_text=ab_json)
    if env is None:
        print("developability 실행 실패 → 중단")
        return 1
    _save("04_developability.json", env)
    gate(env, "developability")
    for d in (env.get("result") or [])[:4]:
        p = d["properties"]
        print(f"  {d['id']:12s} MW={p['molecular_weight_Da']:.1f} pI={p['isoelectric_point_pI']} "
              f"GRAVY={p['gravy_kyte_doolittle']} instab={p['instability_index']} "
              f"liab={d['summary']['n_liability_hits']}(CDR {d['summary']['n_in_cdr']})")

    # (e) Humanness
    env = _run("humanness.py", ["--stdin"], stdin_text=ab_json)
    if env is None:
        print("humanness 실행 실패 → 중단")
        return 1
    _save("05_humanness.json", env)
    gate(env, "humanness")
    for d in (env.get("result") or [])[:4]:
        b = d.get("nearest_germline")
        if b:
            print(f"  {d['id']:12s} nearest {b['germline_gene']} ({b['germline_accession']}) "
                  f"identity={b['germline_identity_percent']}%")

    print(f"\n순차 실행 완료. 산출물: {OUTPUTS}/0{{1..5}}_*.json")
    print("보고서 작성은 Claude Code 의 report-writer 스킬로 수행하세요 "
          "(수치는 위 JSON 봉투에 있는 값만 인용).")
    return 0


# ---------------------------------------------------------- design-check
def cmd_design_check() -> int:
    """설계 경로 A/B 의 GPU-free dry-run 점검 (강의실 노트북에서 사용)."""
    print("== 설계 경로 dry-run 점검 (GPU 불필요) ==\n")
    ok = True

    print("-- 경로 A: ESMFold2 inversion --")
    env = _run("design_esmfold2.py",
               ["--target-name", "pd-l1", "--binder-name", "minibinder", "--dry-run"])
    if env is None:
        print("  실행 실패"); ok = False
    else:
        r, v = env["result"], env["verification"]
        print(f"  passed={v['passed']}")
        for c in v["checks"]:
            print(f"    [{'ok' if c['passed'] else 'FAIL'}] {c['check']}")
        p = r.get("protocol", {})
        print(f"  공식 프로토콜: {p.get('note')} "
              f"(sha256 일치={p.get('matches_verified_sha256')})")
        print(f"  API 심볼: {r.get('api_symbol_check')}")
        print(f"  GPU: {r['environment']['gpu'].get('reason') or 'available'}")
        _save("design_a_dryrun.json", env)
        ok = ok and v["passed"]

    print("\n-- 경로 B: RFantibody --")
    print("  (--target / --framework HLT 파일이 필요하므로 여기서는 CLI·환경만 점검)")
    import shutil as _sh
    tools = {n: _sh.which(n) for n in ("rfdiffusion", "proteinmpnn", "rf2", "qvscorefile")}
    for n, p in tools.items():
        print(f"    {n:12s} {'ok: ' + p if p else 'missing — GPU 노드에서 설치 필요'}")
    print("    dry-run 예시:")
    print("      python scripts/design_rfantibody.py --target antigen.pdb \\")
    print("          --framework framework_HLT.pdb --out outputs/design_rfab \\")
    print('          --loops "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11" \\')
    print('          --hotspots "T570,T593" -n 10 --dry-run')

    print("\n결과:", "PASS — 설계 경로 설정 검증 완료 (실제 설계는 GPU 노드 필요)" if ok
          else "FAIL — 위 항목 확인")
    print("GPU 실행 절차: RUNPOD_가이드.md")
    return 0 if ok else 1


# ---------------------------------------------------------------- list
def cmd_list() -> int:
    print("== 스킬 (.claude/skills/) ==")
    sk = os.path.join(ROOT, ".claude", "skills")
    if os.path.isdir(sk):
        for name in sorted(os.listdir(sk)):
            desc = ""
            p = os.path.join(sk, name, "SKILL.md")
            if os.path.isfile(p):
                with open(p, "r", encoding="utf-8") as fh:
                    for line in fh:
                        if line.startswith("description:"):
                            desc = line.split(":", 1)[1].strip()
                            break
            print(f"  - {name:18s} {desc}")

    print("\n== 스크립트 (scripts/) ==")
    for name in sorted(os.listdir(SCRIPTS)):
        if not name.endswith(".py"):
            continue
        if name in EVAL_STEPS:
            mark = " (평가 파이프라인)"
        elif name == "compare_designs.py":
            mark = " (설계 비교 — GPU 불필요)"
        elif name in DESIGN_STEPS:
            mark = " (설계 — GPU 필요)"
        else:
            mark = " (공용 유틸)"
        print(f"  - {name}{mark}")

    print("\n== 파이프라인 ==")
    print("  [설계]  경로 A design_esmfold2  ─┐")
    print("          경로 B design_rfantibody ┴→ compare_designs")
    print("  [평가]  antigen_lookup → antibody_search → cdr_analysis →")
    print("          developability → humanness → report")
    print("  설계 산출 서열도 같은 평가 3종(cdr/developability/humanness)으로 채점된다.")
    return 0


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    if cmd == "run":
        acc = sys.argv[2] if len(sys.argv) > 2 else "P04626"
        return cmd_run(acc)
    return {"check": cmd_check, "list": cmd_list,
            "design-check": cmd_design_check}.get(cmd, cmd_check)()


if __name__ == "__main__":
    raise SystemExit(main())
