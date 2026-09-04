"""경로 B (고전 비교군) — RFantibody 3단계 파이프라인 오케스트레이션.

RFdiffusion(항체 튜닝) → ProteinMPNN(서열) → RoseTTAFold2(구조 검증·필터)
+ 선택적 ESMFold2 refolding self-consistency (고전 바인더 설계 워크플로 그대로).

무-날조 정책:
- CLI 인자를 **지어내지 않았다.** 아래 인자는 공식 README 에서 직접 읽은 것이다.
- 도구가 없거나 GPU 가 없으면 **실행하지 않고** verification.passed=false 로 보고한다.
  설계 서열·pAE·RMSD 를 만들어내지 않는다.
- 필터 임계값은 저자가 문서에 명시한 값을 그대로 쓴다(발명 금지).

검증된 출처 (2026-09-01 확인):
- repo    : https://github.com/RosettaCommons/RFantibody  (license: MIT, default_branch: main)
- README  : https://raw.githubusercontent.com/RosettaCommons/RFantibody/main/README.md
            sha256(당시): 1a1bd53922d063af09148e01d7486bb6829f233bc8559e6112f681f6b970ed5a
- 프리프린트: https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1
- 가중치  : `bash include/download_weights.sh`
- 설치    : uv (`uv sync`) / Docker / Apptainer. NVIDIA GPU + CUDA 11.8+, Ubuntu 22.04 권장

공식 README 에서 확인한 CLI (추측 아님):
  rfdiffusion  -t/--target <antigen.pdb>  -f/--framework <HLT.pdb>
               -o/--output <prefix> | -q/--output-quiver <file.qv>
               -n/--num-designs N
               -l/--design-loops "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11"
               -h/--hotspots "T305,T456"      # 체인 접두어 + 잔기번호
               --deterministic
  proteinmpnn  -i/--input-dir | -q/--input-quiver
               -o/--output-dir | --output-quiver
               -l/--loops "H1,H2,H3,L1,L2,L3"  (default)
               -n/--seqs-per-struct N   -t/--temperature 0.1   --deterministic
  rf2          -p/--input-pdb | -i/--input-dir | -q/--input-quiver
               -o/--output-dir | --output-quiver
               -r/--num-recycles 10   -s/--seed   --hotspot-show-prop 0.1
  qvscorefile <file.qv>   # 점수를 TSV 로 추출  (그 외 qvls/qvextract/... 유틸)

HLT 포맷 (공식 README): PDB 변형. 체인 ID Heavy='H', Light='L', Target='T',
체인 순서 H→L→T. CDR 은 PDB remark 로 1-indexed 절대 잔기번호 주석:
    REMARK PDBinfo-LABEL:   32 H1

저자가 명시한 한계 (그대로 인용해야 함):
- "The lack of an effective filter is the main limitation of the RFantibody pipeline
   at the moment."
- 일반적으로 hit 를 찾으려면 **10k 규모** 설계 캠페인이 필요.
- 권장 최소 필터: RF2 pAE < 10, RMSD(design vs RF2 예측) < 2 Å, (선택) Rosetta ddG < -20.
- hotspot 선택에 민감 → 대규모 캠페인 전에 pilot run 권장.

사용:
    # GPU 없이 검증 — 명령어·입력만 확인
    python scripts/design_rfantibody.py --target antigen.pdb --framework hu-4D5-8_Fv.pdb \
        --loops "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11" --hotspots "T570,T593" \
        -n 10 --out outputs/design_rfab --dry-run

    # RunPod GPU 에서 실제 실행
    python scripts/design_rfantibody.py --target antigen.pdb --framework fw.pdb \
        --hotspots "T570,T593" -n 100 --out outputs/design_rfab

출력: 표준 결과 봉투(JSON).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time

from verify import emit, make_result

REPO_URL = "https://github.com/RosettaCommons/RFantibody"
README_SHA256_VERIFIED = "1a1bd53922d063af09148e01d7486bb6829f233bc8559e6112f681f6b970ed5a"
PREPRINT = "https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1"

# 저자가 문서에 명시한 권장 최소 필터 (발명 아님)
FILTER_PAE_MAX = 10.0
FILTER_RMSD_MAX_A = 2.0
FILTER_DDG_MAX_OPTIONAL = -20.0

AUTHOR_LIMITATIONS = [
    "저자 문서 인용: \"The lack of an effective filter is the main limitation of the "
    "RFantibody pipeline at the moment.\" — RF2 필터가 binder/non-binder 를 신뢰성 있게 "
    "구분한다는 근거가 아직 부족하다.",
    "저자 문서 인용: 일부 타깃에서는 95 designs 로 VHH binder 를 찾았으나, 일반적으로는 "
    "\"design campaigns in the 10k range\" 가 필요할 것으로 예상된다.",
    "hotspot 선택에 vanilla RFdiffusion 보다 민감 → 대규모 캠페인 전 pilot run 권장.",
    "하전·극성 부위, 글라이칸 인접부, 비구조화 루프 표적은 여전히 어렵다.",
    "런타임이 O(N²) 로 증가 → 큰 항원은 절단(truncation)이 필요하다.",
]

LOOPS_RE = re.compile(r"^(?:[HL][123]:\d+(?:-\d+)?)(?:,[HL][123]:\d+(?:-\d+)?)*$")
HOTSPOT_RE = re.compile(r"^[A-Za-z]\d+$")


# --------------------------------------------------------------- 환경/입력 점검
def tool_status() -> dict:
    """RFantibody CLI 실행 파일 존재 여부 (실행하지 않고 which 만)."""
    return {name: shutil.which(name) for name in
            ("rfdiffusion", "proteinmpnn", "rf2", "qvscorefile", "qvls")}


def gpu_status() -> dict:
    """nvidia-smi 로 GPU 확인 (torch 의존 없음)."""
    exe = shutil.which("nvidia-smi")
    if not exe:
        return {"nvidia_smi": False, "gpus": [], "reason": "nvidia-smi 없음"}
    try:
        p = subprocess.run([exe, "--query-gpu=name,memory.total,driver_version",
                            "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=30)
    except Exception as exc:  # noqa: BLE001
        return {"nvidia_smi": True, "gpus": [], "reason": f"{type(exc).__name__}: {exc}"}
    if p.returncode != 0:
        return {"nvidia_smi": True, "gpus": [], "reason": p.stderr.strip()[:200]}
    gpus = [ln.strip() for ln in p.stdout.splitlines() if ln.strip()]
    return {"nvidia_smi": True, "gpus": gpus, "reason": None if gpus else "GPU 미검출"}


def inspect_hlt(path: str) -> dict:
    """HLT 프레임워크 파일의 형식 점검 (공식 README 규격 기준)."""
    out = {"path": path, "exists": os.path.isfile(path), "chains": [],
           "cdr_labels": [], "chain_order_ok": None, "issues": []}
    if not out["exists"]:
        out["issues"].append("파일 없음")
        return out
    order, seen = [], set()
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK PDBinfo-LABEL:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        out["cdr_labels"].append({"resi": parts[2], "loop": parts[3]})
                elif line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    ch = line[21]
                    if ch not in seen:
                        seen.add(ch)
                        order.append(ch)
    except Exception as exc:  # noqa: BLE001
        out["issues"].append(f"읽기 실패: {type(exc).__name__}: {exc}")
        return out
    out["chains"] = order
    expected = [c for c in ("H", "L", "T") if c in order]
    out["chain_order_ok"] = ([c for c in order if c in ("H", "L", "T")] == expected
                             and bool(expected))
    if not out["chain_order_ok"]:
        out["issues"].append(f"체인 순서가 H→L→T 가 아님 (관측: {order}). "
                             f"공식 HLT 규격 위반 가능")
    unknown = [c for c in order if c not in ("H", "L", "T")]
    if unknown:
        out["issues"].append(f"HLT 규격 외 체인 ID: {unknown}")
    if not out["cdr_labels"]:
        out["issues"].append("REMARK PDBinfo-LABEL CDR 주석 없음 — "
                             "RFdiffusion 이 CDR 을 인식하지 못할 수 있음")
    return out


# ------------------------------------------------------------ self-consistency
def ca_rmsd(pdb_a: str, pdb_b: str, chain: str | None = None):
    """두 구조의 Cα RMSD (Kabsch superposition, BioPython).

    설계 구조 vs 예측(refold) 구조의 self-consistency 지표.
    잔기 수가 다르면 앞에서부터 공통 길이만 사용하고 그 사실을 함께 반환한다.
    """
    try:
        from Bio.PDB import PDBParser, Superimposer
    except Exception as exc:  # noqa: BLE001
        return None, f"BioPython 없음: {exc}"
    try:
        p = PDBParser(QUIET=True)
        sa = p.get_structure("a", pdb_a)
        sb = p.get_structure("b", pdb_b)

        def cas(st):
            out = []
            for model in st:
                for ch in model:
                    if chain and ch.id != chain:
                        continue
                    for res in ch:
                        if "CA" in res:
                            out.append(res["CA"])
                break  # 첫 모델만
            return out

        A, B = cas(sa), cas(sb)
        if not A or not B:
            return None, "Cα 원자를 찾지 못함"
        n = min(len(A), len(B))
        sup = Superimposer()
        sup.set_atoms(A[:n], B[:n])
        note = None if len(A) == len(B) else f"잔기 수 불일치({len(A)} vs {len(B)}) → 앞 {n}개만 사용"
        return round(float(sup.rms), 3), note
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


# ------------------------------------------------------------------ 실행 헬퍼
def build_commands(a) -> list[dict]:
    """3단계 명령을 공식 README 인자 그대로 구성 (실행은 하지 않음)."""
    out = a.out
    qv1 = os.path.join(out, "1_rfdiffusion.qv")
    qv2 = os.path.join(out, "2_proteinmpnn.qv")
    qv3 = os.path.join(out, "3_rf2.qv")

    c1 = ["rfdiffusion", "-t", a.target, "-f", a.framework,
          "--output-quiver", qv1, "-n", str(a.num_designs)]
    if a.loops:
        c1 += ["--design-loops", a.loops]
    if a.hotspots:
        c1 += ["--hotspots", a.hotspots]
    if a.deterministic:
        c1 += ["--deterministic"]

    c2 = ["proteinmpnn", "--input-quiver", qv1, "--output-quiver", qv2,
          "--seqs-per-struct", str(a.seqs_per_struct), "--temperature", str(a.temperature)]
    if a.mpnn_loops:
        c2 += ["--loops", a.mpnn_loops]
    if a.deterministic:
        c2 += ["--deterministic"]

    c3 = ["rf2", "--input-quiver", qv2, "--output-quiver", qv3,
          "--num-recycles", str(a.num_recycles), "--seed", str(a.seed)]

    c4 = ["qvscorefile", qv3]

    return [
        {"step": 1, "name": "RFdiffusion (backbone)", "cmd": c1, "output": qv1},
        {"step": 2, "name": "ProteinMPNN (sequence)", "cmd": c2, "output": qv2},
        {"step": 3, "name": "RoseTTAFold2 (structure + confidence)", "cmd": c3, "output": qv3},
        {"step": 4, "name": "qvscorefile (scores → TSV)", "cmd": c4,
         "output": qv3.replace(".qv", ".sc")},
    ]


def run_cmd(cmd: list[str], timeout: int) -> dict:
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        return {"ok": False, "returncode": None, "error": f"실행 파일 없음: {cmd[0]}",
                "elapsed_s": 0.0}
    except subprocess.TimeoutExpired:
        return {"ok": False, "returncode": None, "error": f"timeout {timeout}s",
                "elapsed_s": round(time.time() - t0, 1)}
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "returncode": None,
                "error": f"{type(exc).__name__}: {exc}", "elapsed_s": 0.0}
    return {"ok": p.returncode == 0, "returncode": p.returncode,
            "stdout_tail": (p.stdout or "")[-1500:], "stderr_tail": (p.stderr or "")[-1500:],
            "elapsed_s": round(time.time() - t0, 1)}


def read_scorefile(path: str):
    """qvscorefile 이 만든 TSV 를 읽는다. 컬럼명은 파일에 있는 그대로 유지."""
    if not os.path.isfile(path):
        return None, f"score 파일 없음: {path}"
    try:
        with open(path, "r", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows, None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def apply_filters(rows: list) -> dict:
    """저자 권장 최소 필터 적용. 컬럼이 없으면 '판정 불가'로 남긴다(추측 금지)."""
    if not rows:
        return {"applied": False, "reason": "score row 0건", "passed": [], "n_passed": 0}
    cols = set(rows[0].keys())
    pae_col = next((c for c in cols if c.lower() in
                    ("pae", "pae_interaction", "rf2_pae", "mean_pae")), None)
    rmsd_col = next((c for c in cols if "rmsd" in c.lower()), None)
    if pae_col is None and rmsd_col is None:
        return {"applied": False,
                "reason": f"pAE/RMSD 컬럼을 찾지 못함. 실제 컬럼: {sorted(cols)}",
                "available_columns": sorted(cols), "passed": [], "n_passed": 0}
    passed, undecidable = [], 0
    for r in rows:
        ok, decided = True, False
        if pae_col:
            try:
                ok = ok and float(r[pae_col]) < FILTER_PAE_MAX
                decided = True
            except (TypeError, ValueError):
                pass
        if rmsd_col:
            try:
                ok = ok and float(r[rmsd_col]) < FILTER_RMSD_MAX_A
                decided = True
            except (TypeError, ValueError):
                pass
        if not decided:
            undecidable += 1
            continue
        if ok:
            passed.append(r)
    return {"applied": True, "pae_column": pae_col, "rmsd_column": rmsd_col,
            "thresholds": {"pAE_max": FILTER_PAE_MAX, "rmsd_max_A": FILTER_RMSD_MAX_A,
                           "optional_rosetta_ddG_max": FILTER_DDG_MAX_OPTIONAL},
            "n_input": len(rows), "n_passed": len(passed),
            "n_undecidable": undecidable, "passed": passed[:200]}


# ------------------------------------------------------------------------ main
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="design_rfantibody.py",
        description="RFantibody 3단계 파이프라인 (경로 B, 고전 비교군)")
    p.add_argument("--target", required=True, help="항원 PDB (rfdiffusion -t)")
    p.add_argument("--framework", required=True,
                   help="HLT 포맷 항체 프레임워크 PDB (rfdiffusion -f)")
    p.add_argument("--out", required=True, help="출력 디렉토리")
    p.add_argument("-n", "--num-designs", type=int, default=10)
    p.add_argument("--loops", default=None,
                   help='RFdiffusion CDR 길이. 예: "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11"')
    p.add_argument("--hotspots", default=None, help='예: "T305,T456" (체인+잔기번호)')
    p.add_argument("--mpnn-loops", default=None,
                   help='ProteinMPNN 설계 루프 (기본 "H1,H2,H3,L1,L2,L3")')
    p.add_argument("--seqs-per-struct", type=int, default=1)
    p.add_argument("--temperature", type=float, default=0.1)
    p.add_argument("--num-recycles", type=int, default=10)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--deterministic", action="store_true")
    p.add_argument("--step-timeout", type=int, default=21600, help="단계별 타임아웃(초)")
    p.add_argument("--self-consistency", action="store_true",
                   help="RF2 예측 vs 설계 백본 Cα RMSD 를 추가 계산 (BioPython)")
    p.add_argument("--dry-run", action="store_true",
                   help="GPU/도구 없이 명령어·입력 형식만 검증")
    return p


def main() -> int:
    a = build_parser().parse_args()

    # 1) 입력 형식 검증 (GPU 불필요 — 강의실 노트북에서도 돌아감)
    input_checks = {
        "target_exists": os.path.isfile(a.target),
        "framework_hlt": inspect_hlt(a.framework),
        "loops_syntax_ok": (a.loops is None) or bool(LOOPS_RE.match(a.loops)),
        "hotspots_syntax_ok": (a.hotspots is None) or all(
            HOTSPOT_RE.match(h.strip()) for h in a.hotspots.split(",") if h.strip()),
        "num_designs_positive": a.num_designs > 0,
    }
    if a.loops and not input_checks["loops_syntax_ok"]:
        input_checks["loops_hint"] = ('형식: "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11" '
                                      '(단일값=고정 길이, 범위=균등 샘플링)')
    if a.hotspots and not input_checks["hotspots_syntax_ok"]:
        input_checks["hotspots_hint"] = '형식: "T305,T456" (체인 문자 + 잔기번호)'

    cmds = build_commands(a)
    env_info = {"tools": tool_status(), "gpu": gpu_status()}
    refs = {"repo": REPO_URL, "preprint": PREPRINT, "license": "MIT",
            "readme_sha256_verified": README_SHA256_VERIFIED,
            "weights": "bash include/download_weights.sh",
            "requirements": "NVIDIA GPU, CUDA 11.8+, Ubuntu 22.04 권장"}

    base_ok = (input_checks["target_exists"]
               and input_checks["framework_hlt"]["exists"]
               and input_checks["loops_syntax_ok"]
               and input_checks["hotspots_syntax_ok"]
               and input_checks["num_designs_positive"])

    # 2) dry-run
    if a.dry_run:
        checks = [
            ("항원 PDB 존재", input_checks["target_exists"]),
            ("프레임워크 HLT 파일 존재", input_checks["framework_hlt"]["exists"]),
            ("HLT 체인 순서 H→L→T", bool(input_checks["framework_hlt"].get("chain_order_ok"))),
            ("--loops 문법 유효", input_checks["loops_syntax_ok"]),
            ("--hotspots 문법 유효", input_checks["hotspots_syntax_ok"]),
        ]
        emit(make_result(
            {"mode": "dry-run", "commands": [{**c, "cmd_str": " ".join(c["cmd"])}
                                             for c in cmds],
             "input_checks": input_checks, "environment": env_info,
             "filters": {"pAE_max": FILTER_PAE_MAX, "rmsd_max_A": FILTER_RMSD_MAX_A,
                         "optional_rosetta_ddG_max": FILTER_DDG_MAX_OPTIONAL,
                         "source": "RFantibody README §Filtering Strategies"},
             "author_limitations": AUTHOR_LIMITATIONS, "references": refs,
             "gpu_required": True,
             "next_step": "RFantibody 설치된 GPU 노드에서 --dry-run 없이 재실행 "
                          "(RUNPOD_가이드.md)."},
            "design_rfantibody dry-run (no GPU, no tool execution)",
            f"target={a.target} framework={a.framework} n={a.num_designs}",
            checks,
            notes=("DRY-RUN — 아무 도구도 실행하지 않았습니다. 명령어 구성과 입력 형식만 검증. "
                   "설계 서열·pAE·RMSD 를 생성하지 않습니다(무-날조).")))
        return 0 if all(ok for _, ok in checks) else 1

    # 3) 실제 실행 — 도구 + GPU 필수
    missing = [k for k, v in env_info["tools"].items()
               if v is None and k in ("rfdiffusion", "proteinmpnn", "rf2")]
    if missing or not env_info["gpu"]["gpus"]:
        emit(make_result(
            {"mode": "aborted", "commands": [{**c, "cmd_str": " ".join(c["cmd"])}
                                             for c in cmds],
             "input_checks": input_checks, "environment": env_info,
             "author_limitations": AUTHOR_LIMITATIONS, "references": refs},
            "design_rfantibody (PRECONDITION FAILED)",
            f"target={a.target} framework={a.framework}",
            [("RFantibody CLI 설치 (rfdiffusion/proteinmpnn/rf2)", not missing),
             ("NVIDIA GPU 검출", bool(env_info["gpu"]["gpus"]))],
            notes=(f"실행 전제 미충족 — 누락 CLI: {missing or '없음'}, "
                   f"GPU: {env_info['gpu'].get('reason') or env_info['gpu']['gpus']}. "
                   f"설계 결과를 생성하지 않습니다. RunPod GPU 노드에서 "
                   f"`{REPO_URL}` 설치 후 실행하십시오 (RUNPOD_가이드.md).")))
        return 1
    if not base_ok:
        emit(make_result({"mode": "aborted", "input_checks": input_checks},
                         "design_rfantibody (INPUT INVALID)", a.framework,
                         [("입력 파일/문법 유효", False)],
                         notes=f"입력 검증 실패: {input_checks}"))
        return 1

    os.makedirs(a.out, exist_ok=True)
    step_logs, aborted = [], False
    for c in cmds:
        res = run_cmd(c["cmd"], timeout=a.step_timeout)
        step_logs.append({"step": c["step"], "name": c["name"],
                          "cmd_str": " ".join(c["cmd"]), "output": c["output"], **res})
        if not res["ok"]:
            aborted = True
            break

    score_path = cmds[3]["output"]
    rows, srr = read_scorefile(score_path)
    filt = apply_filters(rows or [])

    # 4) 선택적 self-consistency (RF2 예측 vs 설계 백본)
    sc = None
    if a.self_consistency and not aborted:
        sc = {"method": "Cα RMSD (BioPython Superimposer, Kabsch)",
              "note": "RF2 예측 구조와 설계 백본을 비교. "
                      "Quiver 에서 PDB 추출(qvextract) 후 쌍을 지어 계산해야 하며, "
                      "추출 파일이 없으면 수행하지 않는다.",
              "pairs": [], "skipped_reason": None}
        d1 = os.path.join(a.out, "extracted_rfdiffusion")
        d3 = os.path.join(a.out, "extracted_rf2")
        if os.path.isdir(d1) and os.path.isdir(d3):
            for fn in sorted(os.listdir(d3)):
                if not fn.endswith(".pdb"):
                    continue
                p3, p1 = os.path.join(d3, fn), os.path.join(d1, fn)
                if not os.path.isfile(p1):
                    continue
                rms, note = ca_rmsd(p1, p3)
                sc["pairs"].append({"design": p1, "prediction": p3,
                                    "ca_rmsd_A": rms, "note": note,
                                    "passes_rmsd_filter": (rms is not None
                                                           and rms < FILTER_RMSD_MAX_A)})
        else:
            sc["skipped_reason"] = (
                f"추출 디렉토리 없음 ({d1}, {d3}). "
                f"`qvextract` 로 각 Quiver 를 PDB 로 풀어 놓은 뒤 재실행하십시오. "
                f"RMSD 를 추정값으로 만들어내지 않습니다(무-날조).")

    checks = [
        ("RFantibody CLI + GPU 확인", True),
        ("4 단계 전부 성공", not aborted and all(s["ok"] for s in step_logs)),
        ("score 파일 파싱", rows is not None and len(rows) > 0),
        ("저자 권장 필터 적용됨", bool(filt.get("applied"))),
    ]
    notes = (f"RFantibody 3단계 + qvscorefile 실행. "
             f"score row {len(rows or [])}건, 필터 통과 {filt.get('n_passed', 0)}건 "
             f"(pAE<{FILTER_PAE_MAX}, RMSD<{FILTER_RMSD_MAX_A}Å — 저자 문서 권장값). "
             f"필터 신뢰도에 대한 저자 자체 경고를 반드시 함께 보고할 것. "
             + (f"score 파일 문제: {srr}" if srr else ""))
    emit(make_result(
        {"mode": "design", "steps": step_logs, "score_file": score_path,
         "score_rows": (rows or [])[:200], "filter_result": filt,
         "self_consistency": sc, "input_checks": input_checks,
         "environment": env_info, "author_limitations": AUTHOR_LIMITATIONS,
         "references": refs},
        "RFantibody (RFdiffusion → ProteinMPNN → RoseTTAFold2)",
        f"target={a.target} framework={a.framework} n={a.num_designs} "
        f"loops={a.loops} hotspots={a.hotspots}",
        checks, notes=notes))
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result(None, "design_rfantibody (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
