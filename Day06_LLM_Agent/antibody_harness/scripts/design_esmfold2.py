"""경로 A (주축) — ESMFold2 inversion 으로 de novo 바인더/항체 설계.

무-날조 정책:
- 이 스크립트는 API 를 **재구현하지 않는다.** Biohub 공식 저장소의 프로토콜 스크립트
  `cookbook/tutorials/binder_design.py` 를 그대로 내려받아 import 한다.
  (직접 베낀 코드가 없으므로 함수명·인자를 지어낼 여지가 없다.)
- 설계 서열·점수는 모델이 반환한 값만 기록한다. GPU 가 없으면 **실행하지 않고**
  `verification.passed=false` + "RunPod GPU 필요" 로 정직 보고한다.

검증된 출처 (2026-09-01 확인):
- repo      : https://github.com/Biohub/esm  (MIT, LICENSE.md)
- 프로토콜  : https://github.com/Biohub/esm/blob/main/cookbook/tutorials/binder_design.py
              raw  : PROTOCOL_URL 상수
              커밋 : 827ec128e4cdaf80f7d6f95fb367a08980b34918 (해당 경로 최신)
              sha256(당시): 28d672a3b1ff722e6d3d50f7538b806bbaf04291c70d83454fa6c912869cd3d3
- 가중치    : https://huggingface.co/biohub/ESMFold2 (license:mit), biohub/ESMC-6B
- 설치      : binder_design.py 헤더가 명시한 의존성
              esm@git+https://github.com/Biohub/esm.git@main , abnumber , modal
- 프리프린트: https://www.biorxiv.org/content/10.64898/2026.06.03.729735

공식 프로토콜에서 확인한 API (추측 아님, 파일에서 직접 읽음):
    class ESMFold2Design:
        lm_name = "biohub/ESMC-6B"
        inversion_model_names = ["ESMFold2-Experimental-Fast",
                                 "ESMFold2-Experimental-Fast-Cutoff2025"]
        def load(self, use_scaling_critics: bool) -> None
        def design(self, target_name, binder_name, target_sequence=None,
                   binder_sequence=None, is_antibody=None, seed=0, batch_size=1,
                   target_hotspot_ids=None, epitope_contact_distance=12.0)
            -> tuple[list[str], Trajectory, list[dict]]
    TARGET_SEQUENCES      : preset 항원 5종 (cd45, ctla4, egfr, pd-l1, pdgfr)
    BINDER_PROMPT_FACTORIES: minibinder / trastuzumab_framework_vhvl /
                             atezolizumab_framework_vhvl / ocankitug_framework_vhvl
    critic 결과 dict 키   : designed_sequence, complex, final_loss, iptm, critic_name,
                            is_scaling_critic, batch_idx, is_antibody, logits,
                            distogram_iptm_proxy, cdr_distogram_iptm_proxy

사용:
    # (1) GPU 없이 검증 가능 — 설정 확인만, 모델 로드 안 함
    python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder --dry-run

    # (2) 커스텀 항원(예: HER2 ECD 일부) + 항체 프레임워크
    python scripts/design_esmfold2.py --target-name her2_ecd4 \
        --target-sequence-file outputs/her2_epitope.fasta \
        --binder-name trastuzumab_framework_vhvl --hotspots T570,T593 --dry-run

    # (3) RunPod A100/H100 에서 실제 설계
    python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder \
        --seed 0 --batch-size 1 --num-seeds 4 --out outputs/design_esmfold2

출력: 표준 결과 봉투(JSON).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

from verify import emit, make_result

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VENDOR = os.path.join(ROOT, "vendor")
PROTOCOL_URL = ("https://raw.githubusercontent.com/Biohub/esm/main/"
                "cookbook/tutorials/binder_design.py")
PROTOCOL_PATH = os.path.join(VENDOR, "binder_design.py")
PROTOCOL_SHA256_VERIFIED = "28d672a3b1ff722e6d3d50f7538b806bbaf04291c70d83454fa6c912869cd3d3"
PROTOCOL_COMMIT_VERIFIED = "827ec128e4cdaf80f7d6f95fb367a08980b34918"

REPO_URL = "https://github.com/Biohub/esm"
HF_MODELS = ["biohub/ESMFold2", "biohub/ESMC-6B"]
PIP_SPEC = "esm@git+https://github.com/Biohub/esm.git@main"

# 공식 프로토콜에 하드코딩된 preset 이름 (참조용 — 실제 값은 vendor 파일에서 읽는다)
KNOWN_TARGET_PRESETS = ["cd45", "ctla4", "egfr", "pd-l1", "pdgfr"]
KNOWN_BINDER_PRESETS = ["minibinder", "trastuzumab_framework_vhvl",
                        "atezolizumab_framework_vhvl", "ocankitug_framework_vhvl"]


# ------------------------------------------------------------------ 환경 점검
def gpu_status() -> dict:
    """GPU 가용성 점검. torch 미설치도 정직하게 보고."""
    try:
        import torch
    except Exception as exc:  # noqa: BLE001
        return {"torch": False, "cuda": False, "reason": f"torch 미설치: {exc}",
                "devices": []}
    try:
        avail = bool(torch.cuda.is_available())
        devs = []
        if avail:
            for i in range(torch.cuda.device_count()):
                p = torch.cuda.get_device_properties(i)
                devs.append({"index": i, "name": p.name,
                             "total_memory_GiB": round(p.total_memory / 1024**3, 1)})
        return {"torch": True, "torch_version": torch.__version__, "cuda": avail,
                "reason": None if avail else "CUDA 사용 불가", "devices": devs}
    except Exception as exc:  # noqa: BLE001
        return {"torch": True, "cuda": False, "reason": f"{type(exc).__name__}: {exc}",
                "devices": []}


def esm_status() -> dict:
    """esm 패키지와 experimental 모듈 존재 점검 (import 만, 가중치 다운로드 없음)."""
    out = {"esm_installed": False, "experimental_module": False, "abnumber": False,
           "reason": None}
    try:
        import esm  # noqa: F401
        out["esm_installed"] = True
        out["esm_version"] = getattr(esm, "__version__", "unknown")
    except Exception as exc:  # noqa: BLE001
        out["reason"] = f"esm 미설치: {type(exc).__name__}: {exc}"
        return out
    try:
        from esm.models.esmfold2.experimental import EsmFold2ExperimentalModel  # noqa: F401
        out["experimental_module"] = True
    except Exception as exc:  # noqa: BLE001
        out["reason"] = (f"esm.models.esmfold2.experimental import 실패 "
                         f"({type(exc).__name__}: {exc}) — PyPI 판이 아니라 "
                         f"`pip install {PIP_SPEC}` 로 설치해야 할 수 있음")
    try:
        import abnumber  # noqa: F401
        out["abnumber"] = True
    except Exception:  # noqa: BLE001
        pass
    return out


# --------------------------------------------------- 공식 프로토콜 스크립트 확보
def fetch_protocol(force: bool = False) -> tuple[str | None, str | None, str]:
    """공식 binder_design.py 를 vendor/ 로 내려받는다. (path, sha256, note)."""
    if os.path.isfile(PROTOCOL_PATH) and not force:
        with open(PROTOCOL_PATH, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        return PROTOCOL_PATH, digest, "기존 vendor 사본 사용"
    try:
        import requests
    except Exception as exc:  # noqa: BLE001
        return None, None, f"requests 미설치: {exc}"
    try:
        r = requests.get(PROTOCOL_URL, timeout=60)
    except Exception as exc:  # noqa: BLE001
        return None, None, f"다운로드 네트워크 오류: {type(exc).__name__}: {exc}"
    if r.status_code != 200:
        return None, None, f"다운로드 HTTP {r.status_code}"
    os.makedirs(VENDOR, exist_ok=True)
    with open(PROTOCOL_PATH, "wb") as fh:
        fh.write(r.content)
    digest = hashlib.sha256(r.content).hexdigest()
    note = "새로 다운로드"
    if digest != PROTOCOL_SHA256_VERIFIED:
        note += (f" — 주의: sha256 이 검증 시점({PROTOCOL_SHA256_VERIFIED[:12]}…)과 다름. "
                 f"upstream main 이 갱신된 것이며, API 가 바뀌었을 수 있으니 "
                 f"vendor/binder_design.py 를 직접 확인하십시오.")
    return PROTOCOL_PATH, digest, note


def load_protocol(path: str):
    """vendor 사본을 모듈로 로드."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("esm_binder_design", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["esm_binder_design"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def read_sequence_file(path: str) -> str | None:
    from seq_utils import parse_fasta, clean_seq

    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as fh:
        txt = fh.read()
    if ">" in txt:
        recs = parse_fasta(txt)
        return recs[0]["sequence"] if recs else None
    return clean_seq(txt) or None


# ------------------------------------------------------------------------ main
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="design_esmfold2.py",
        description="ESMFold2 inversion 기반 de novo 바인더/항체 설계 (경로 A)")
    p.add_argument("--target-name", required=True,
                   help=f"항원 이름. preset: {', '.join(KNOWN_TARGET_PRESETS)}. "
                        f"그 외 이름이면 --target-sequence(-file) 필수")
    p.add_argument("--target-sequence", default=None, help="항원(에피토프) 서열 문자열")
    p.add_argument("--target-sequence-file", default=None,
                   help="항원 서열 FASTA/텍스트 파일")
    p.add_argument("--binder-name", default="minibinder",
                   help=f"바인더 프롬프트. preset: {', '.join(KNOWN_BINDER_PRESETS)}")
    p.add_argument("--binder-sequence", default=None,
                   help="커스텀 바인더 프롬프트 서열 ('#' = 가변 위치)")
    p.add_argument("--is-antibody", choices=["true", "false"], default=None,
                   help="항체 포맷 여부. preset 사용 시 생략(자동)")
    p.add_argument("--hotspots", default=None,
                   help="에피토프 hotspot, 쉼표 구분. 공식 API 는 target 서열 기준 "
                        "1-indexed 잔기 ID 리스트 (예: 'L150' 또는 'T305,T456')")
    p.add_argument("--epitope-contact-distance", type=float, default=12.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--num-seeds", type=int, default=1,
                   help="서로 다른 seed 로 몇 번 설계할지 (seed, seed+1, ...)")
    p.add_argument("--batch-size", type=int, default=1)
    p.add_argument("--no-scaling-critics", action="store_true",
                   help="scaling critic 비활성 (VRAM/시간 절약; 로컬 uv 실행 시 필요)")
    p.add_argument("--out", default=None, help="설계 산출물 디렉토리 (PDB/mmCIF 저장)")
    p.add_argument("--fetch-protocol", action="store_true",
                   help="공식 binder_design.py 를 강제로 다시 내려받음")
    p.add_argument("--dry-run", action="store_true",
                   help="GPU/모델 없이 설정·환경만 검증 (강의실 노트북에서 사용)")
    return p


def main() -> int:
    a = build_parser().parse_args()

    # 1) 항원 서열 결정
    tseq = a.target_sequence
    if tseq is None and a.target_sequence_file:
        tseq = read_sequence_file(a.target_sequence_file)
        if tseq is None:
            emit(make_result(None, "input validation", a.target_sequence_file,
                             [("항원 서열 파일 읽기", False)],
                             notes=f"'{a.target_sequence_file}' 를 읽지 못했습니다."))
            return 1
    is_preset_target = a.target_name in KNOWN_TARGET_PRESETS
    if is_preset_target and tseq is not None:
        emit(make_result(None, "input validation", a.target_name,
                         [("preset 항원에 서열 중복 지정 안 함", False)],
                         notes=f"'{a.target_name}' 은 공식 프로토콜의 preset 항원입니다. "
                               f"--target-sequence 를 생략하십시오 "
                               f"(design_binder() 가 ValueError 를 던짐)."))
        return 1
    if not is_preset_target and tseq is None:
        emit(make_result(None, "input validation", a.target_name,
                         [("커스텀 항원에 서열 제공", False)],
                         notes=f"'{a.target_name}' 은 preset 이 아니므로 "
                               f"--target-sequence 또는 --target-sequence-file 이 필요합니다. "
                               f"preset: {KNOWN_TARGET_PRESETS}"))
        return 1

    hotspots = [h.strip() for h in (a.hotspots or "").split(",") if h.strip()] or None
    is_ab = None if a.is_antibody is None else (a.is_antibody == "true")

    plan = {
        "track": "A — ESMFold2 inversion",
        "target_name": a.target_name,
        "target_is_preset": is_preset_target,
        "target_sequence_length": len(tseq) if tseq else None,
        "binder_name": a.binder_name,
        "binder_is_preset": a.binder_name in KNOWN_BINDER_PRESETS,
        "is_antibody": is_ab,
        "target_hotspot_ids": hotspots,
        "epitope_contact_distance": a.epitope_contact_distance,
        "seeds": list(range(a.seed, a.seed + a.num_seeds)),
        "batch_size": a.batch_size,
        "use_scaling_critics": not a.no_scaling_critics,
        "output_dir": a.out,
    }
    env_info = {"gpu": gpu_status(), "esm": esm_status()}
    refs = {
        "repo": REPO_URL, "protocol_url": PROTOCOL_URL,
        "protocol_commit_verified": PROTOCOL_COMMIT_VERIFIED,
        "hf_models": HF_MODELS, "pip_spec": PIP_SPEC, "license": "MIT",
    }

    # 2) 공식 프로토콜 확보
    ppath, pdigest, pnote = fetch_protocol(force=a.fetch_protocol)
    protocol = {"path": ppath, "sha256": pdigest, "note": pnote,
                "matches_verified_sha256": pdigest == PROTOCOL_SHA256_VERIFIED}

    # 3) dry-run: GPU 없이 여기까지만
    if a.dry_run:
        api_ok, api_note = False, "vendor 프로토콜 미확보"
        if ppath:
            try:
                src = open(ppath, "r", encoding="utf-8").read()
                needed = ["class ESMFold2Design", "def design_binder(",
                          "TARGET_SEQUENCES", "BINDER_PROMPT_FACTORIES"]
                missing = [n for n in needed if n not in src]
                api_ok = not missing
                api_note = ("공식 프로토콜에서 기대 심볼 전부 확인" if api_ok
                            else f"기대 심볼 누락: {missing} — upstream API 변경 가능")
            except Exception as exc:  # noqa: BLE001
                api_note = f"프로토콜 파일 읽기 실패: {exc}"
        checks = [
            ("입력 설정 유효 (항원/바인더 조합)", True),
            ("공식 프로토콜 스크립트 확보", ppath is not None),
            ("프로토콜 내 기대 API 심볼 존재", api_ok),
        ]
        emit(make_result(
            {"mode": "dry-run", "plan": plan, "environment": env_info,
             "protocol": protocol, "api_symbol_check": api_note, "references": refs,
             "gpu_required": True,
             "next_step": "GPU 노드(A100 40GB+)에서 --dry-run 없이 재실행. "
                          "RUNPOD_가이드.md 참조."},
            "design_esmfold2 dry-run (no GPU, no model load)",
            f"target={a.target_name} binder={a.binder_name} seeds={plan['seeds']}",
            checks,
            notes=("DRY-RUN — 설계를 수행하지 않았습니다. 설정·환경·공식 API 심볼만 검증. "
                   "설계 서열/점수는 생성하지 않습니다(무-날조). "
                   f"GPU: {env_info['gpu'].get('reason') or 'available'}")))
        return 0 if all(ok for _, ok in checks) else 1

    # 4) 실제 실행 — GPU 필수
    g = env_info["gpu"]
    if not g.get("cuda"):
        emit(make_result(
            {"mode": "aborted", "plan": plan, "environment": env_info,
             "protocol": protocol, "references": refs},
            "design_esmfold2 (GPU UNAVAILABLE)",
            f"target={a.target_name} binder={a.binder_name}",
            [("CUDA GPU 사용 가능", False)],
            notes=(f"CUDA GPU 없음 ({g.get('reason')}). ESMFold2 inversion 은 "
                   f"A100 40GB 이상을 요구합니다(공식 프로토콜 주석: batch_size=1, "
                   f"REUSE_ESMC=True 기준 27GB, False 기준 51GB VRAM). "
                   f"설계 결과를 생성하지 않습니다 — RunPod 에서 실행하십시오 "
                   f"(RUNPOD_가이드.md).")))
        return 1
    if ppath is None:
        emit(make_result({"mode": "aborted", "protocol": protocol},
                         "design_esmfold2 (PROTOCOL FETCH FAILED)", PROTOCOL_URL,
                         [("공식 프로토콜 스크립트 확보", False)],
                         notes=f"공식 binder_design.py 를 가져오지 못했습니다: {pnote}. "
                               f"API 를 추측해 재구현하지 않습니다(무-날조)."))
        return 1

    try:
        bd = load_protocol(ppath)
    except Exception as exc:  # noqa: BLE001
        emit(make_result({"mode": "aborted", "protocol": protocol, "environment": env_info},
                         "design_esmfold2 (PROTOCOL IMPORT FAILED)", ppath,
                         [("공식 프로토콜 import", False)],
                         notes=f"import 실패: {type(exc).__name__}: {exc}. "
                               f"`pip install {PIP_SPEC}` 및 abnumber 설치를 확인하십시오."))
        return 1

    app = bd.ESMFold2Design()
    t0 = time.time()
    app.load(use_scaling_critics=not a.no_scaling_critics)
    load_s = round(time.time() - t0, 1)

    if a.out:
        os.makedirs(a.out, exist_ok=True)

    designs, errors = [], []
    for seed in plan["seeds"]:
        try:
            t1 = time.time()
            seqs, trajectory, results = app.design(
                target_name=a.target_name,
                target_sequence=tseq,
                binder_name=a.binder_name,
                binder_sequence=a.binder_sequence,
                is_antibody=is_ab,
                seed=seed,
                batch_size=a.batch_size,
                target_hotspot_ids=hotspots,
                epitope_contact_distance=a.epitope_contact_distance,
            )
            elapsed = round(time.time() - t1, 1)
            for r in results:
                rec = {
                    "seed": seed,
                    "designed_sequence": r.get("designed_sequence"),
                    "binder_sequence": (r.get("designed_sequence") or "").split("|")[-1],
                    "critic_name": r.get("critic_name"),
                    "is_scaling_critic": r.get("is_scaling_critic"),
                    "is_antibody": r.get("is_antibody"),
                    "batch_idx": r.get("batch_idx"),
                    "final_loss": r.get("final_loss"),
                    "iptm": r.get("iptm"),
                    "distogram_iptm_proxy": r.get("distogram_iptm_proxy"),
                    "cdr_distogram_iptm_proxy": r.get("cdr_distogram_iptm_proxy"),
                    "elapsed_s": elapsed,
                }
                cx = r.get("complex")
                if cx is not None and a.out:
                    fn = (f"design_seed{seed}_b{r.get('batch_idx')}_"
                          f"{r.get('critic_name')}.cif")
                    fp = os.path.join(a.out, fn)
                    try:
                        with open(fp, "w", encoding="utf-8") as fh:
                            fh.write(cx.to_mmcif())
                        rec["structure_file"] = fp
                    except Exception as exc:  # noqa: BLE001
                        rec["structure_file"] = None
                        rec["structure_error"] = f"{type(exc).__name__}: {exc}"
                designs.append(rec)
            if a.out:
                tp = os.path.join(a.out, f"trajectory_seed{seed}.json")
                try:
                    with open(tp, "w", encoding="utf-8") as fh:
                        json.dump({str(k): {kk: (vv if isinstance(vv, (int, float, str))
                                                 else str(vv))
                                            for kk, vv in v.items()}
                                   for k, v in trajectory.items()}, fh)
                except Exception:  # noqa: BLE001
                    pass
        except Exception as exc:  # noqa: BLE001
            errors.append(f"seed={seed}: {type(exc).__name__}: {exc}")

    # ranked (iptm 우선, 없으면 distogram proxy, 그다음 loss 낮은 순)
    def rank_key(d):
        return (-(d.get("iptm") if d.get("iptm") is not None else -1),
                -(d.get("distogram_iptm_proxy") or -1),
                d.get("final_loss") if d.get("final_loss") is not None else 1e9)

    ranked = sorted([d for d in designs if not d.get("is_scaling_critic")], key=rank_key)

    checks = [
        ("CUDA GPU 사용", True),
        ("공식 프로토콜 로드", True),
        ("설계 결과 ≥ 1", len(designs) >= 1),
        ("모든 seed 정상 완료 (오류 0건)", len(errors) == 0),
        ("설계 서열이 비어있지 않음",
         all(d.get("binder_sequence") for d in designs) and len(designs) >= 1),
    ]
    notes = (f"ESMFold2 inversion 설계 {len(designs)}건 "
             f"(seed {plan['seeds']}, batch_size={a.batch_size}, "
             f"모델 로드 {load_s}s). ranked 는 iptm → distogram proxy → final_loss 순. "
             f"점수는 모두 공식 프로토콜이 반환한 값 (무-날조). "
             + (f"오류: {errors}" if errors else ""))
    emit(make_result(
        {"mode": "design", "plan": plan, "environment": env_info, "protocol": protocol,
         "references": refs, "designs": designs, "ranked": ranked,
         "model_load_seconds": load_s},
        "ESMFold2 inversion (Biohub/esm cookbook/tutorials/binder_design.py)",
        f"target={a.target_name} binder={a.binder_name} seeds={plan['seeds']}",
        checks, notes=notes))
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result(None, "design_esmfold2 (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
