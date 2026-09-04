"""경로 A(ESMFold2 inversion) vs 경로 B(RFantibody) 산출물을 동일 지표로 비교.

무-날조 정책:
- **두 경로의 점수는 서로 다른 척도다.** ESMFold2 의 ipTM 과 RF2 의 pAE 를 하나의
  "종합 점수"로 합치지 않는다. 그런 통합 점수는 근거 없는 발명이다.
- 대신 (a) 각 경로가 실제로 낸 신뢰도 지표를 **원본 이름 그대로** 병기하고,
  (b) **양쪽에 동일하게 적용 가능한 서열 기반 지표**(developability, humanness,
  CDR 길이)만 공통 비교축으로 삼는다.
- 한쪽 경로 결과가 없으면 그 칸은 null 로 남긴다. 채워 넣지 않는다.

공통 비교축 (서열만 있으면 GPU 없이 계산 가능 — 실제로 검증됨):
    chain_type, CDR-H1/H2/H3·L1/L2/L3 및 길이   (cdr_analysis)
    MW / pI / GRAVY / instability index          (developability, BioPython ProtParam)
    liability 히트 수 (전체 / CDR 내부)          (developability, 규칙 기반)
    germline identity (%) + nearest germline     (humanness, UniProt + BLOSUM62)

경로별 고유 지표 (원본 이름 유지, 교차 비교 금지):
    A: iptm, distogram_iptm_proxy, cdr_distogram_iptm_proxy, final_loss, critic_name
    B: RF2 score 파일의 컬럼 그대로(pAE 등), Cα RMSD self-consistency

사용:
    python scripts/compare_designs.py \
        --track-a outputs/design_esmfold2.json \
        --track-b outputs/design_rfantibody.json \
        --out outputs/06_comparison.json

    # 한쪽만 있어도 동작 (없는 쪽은 null)
    python scripts/compare_designs.py --track-a outputs/design_esmfold2.json

출력: 표준 결과 봉투(JSON).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

from seq_utils import clean_seq
from verify import emit, make_result

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def load_envelope(path: str | None):
    if not path:
        return None, "미지정"
    if not os.path.isfile(path):
        return None, f"파일 없음: {path}"
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh), None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def extract_track_a(env) -> list:
    """design_esmfold2.py 봉투 → 공통 레코드."""
    if not env:
        return []
    r = env.get("result") or {}
    if not isinstance(r, dict) or r.get("mode") != "design":
        return []
    out = []
    for d in r.get("designs", []):
        seq = clean_seq(d.get("binder_sequence") or "")
        if not seq:
            continue
        out.append({
            "track": "A_esmfold2_inversion",
            "id": f"A_seed{d.get('seed')}_b{d.get('batch_idx')}_{d.get('critic_name')}",
            "sequence": seq,
            "is_antibody": d.get("is_antibody"),
            "structure_file": d.get("structure_file"),
            "track_specific_metrics": {
                "iptm": d.get("iptm"),
                "distogram_iptm_proxy": d.get("distogram_iptm_proxy"),
                "cdr_distogram_iptm_proxy": d.get("cdr_distogram_iptm_proxy"),
                "final_loss": d.get("final_loss"),
                "critic_name": d.get("critic_name"),
                "is_scaling_critic": d.get("is_scaling_critic"),
            },
        })
    return out


def extract_track_b(env) -> list:
    """design_rfantibody.py 봉투 → 공통 레코드.

    RF2 score 파일에 서열 컬럼이 있으면 서열 기반 공통 지표도 계산 가능하다.
    없으면 서열 없이 track-specific 지표만 남긴다(서열을 만들어내지 않는다).
    """
    if not env:
        return []
    r = env.get("result") or {}
    if not isinstance(r, dict) or r.get("mode") != "design":
        return []
    rows = r.get("score_rows") or []
    sc_by_pred = {}
    for p in ((r.get("self_consistency") or {}).get("pairs") or []):
        sc_by_pred[os.path.basename(p.get("prediction") or "")] = p.get("ca_rmsd_A")

    out = []
    for i, row in enumerate(rows):
        seq_col = next((c for c in row
                        if c.lower() in ("sequence", "seq", "designed_sequence")), None)
        seq = clean_seq(row.get(seq_col, "")) if seq_col else ""
        tag = row.get("tag") or row.get("description") or f"row{i}"
        out.append({
            "track": "B_rfantibody",
            "id": f"B_{tag}",
            "sequence": seq or None,
            "sequence_available": bool(seq),
            "is_antibody": True,
            "structure_file": None,
            "track_specific_metrics": {
                **{k: v for k, v in row.items() if k != seq_col},
                "ca_rmsd_self_consistency_A": sc_by_pred.get(f"{tag}.pdb"),
            },
        })
    return out


def run_evaluator(script: str, records: list) -> tuple[dict | None, str | None]:
    """공통 평가 스크립트를 --stdin 으로 호출 (기존 검증된 코드 재사용)."""
    payload = json.dumps({
        "result": [{"id": r["id"], "sequence": r["sequence"]}
                   for r in records if r.get("sequence")],
        "provenance": {"source": "compare_designs", "query": "", "timestamp": ""},
        "verification": {"passed": True, "checks": [], "notes": ""},
    }, ensure_ascii=False)
    try:
        p = subprocess.run([PY, os.path.join(SCRIPTS, script), "--stdin"],
                           input=payload, capture_output=True, text=True, timeout=900)
    except Exception as exc:  # noqa: BLE001
        return None, f"{script} 실행 실패: {type(exc).__name__}: {exc}"
    try:
        return json.loads(p.stdout), None
    except Exception as exc:  # noqa: BLE001
        return None, f"{script} 출력 파싱 실패: {exc}; stderr={p.stderr[:300]}"


def main() -> int:
    ap = argparse.ArgumentParser(
        prog="compare_designs.py",
        description="경로 A/B 설계 산출물의 동일-지표 비교")
    ap.add_argument("--track-a", default=None, help="design_esmfold2.py 출력 JSON")
    ap.add_argument("--track-b", default=None, help="design_rfantibody.py 출력 JSON")
    ap.add_argument("--out", default=None, help="비교 결과 JSON 저장 경로")
    a = ap.parse_args()

    if not a.track_a and not a.track_b:
        emit(make_result(None, "input validation", "no input",
                         [("최소 한 경로의 결과 제공", False)],
                         notes="--track-a 또는 --track-b 중 하나 이상이 필요합니다."))
        return 1

    env_a, err_a = load_envelope(a.track_a)
    env_b, err_b = load_envelope(a.track_b)
    recs = extract_track_a(env_a) + extract_track_b(env_b)
    with_seq = [r for r in recs if r.get("sequence")]

    # 공통 지표 계산 (기존 검증된 평가 스크립트 재사용)
    common, eval_errors = {}, []
    if with_seq:
        for key, script in (("cdr", "cdr_analysis.py"),
                            ("developability", "developability.py"),
                            ("humanness", "humanness.py")):
            env, err = run_evaluator(script, with_seq)
            if err:
                eval_errors.append(err)
                continue
            if not env.get("verification", {}).get("passed"):
                eval_errors.append(
                    f"{script}: verification.passed=false — "
                    f"{env.get('verification', {}).get('notes', '')[:200]}")
            common[key] = env

    # id 기준으로 공통 지표를 각 레코드에 병합
    def index_by_id(env, id_key="id"):
        idx = {}
        for r in ((env or {}).get("result") or []):
            base = str(r.get(id_key, "")).split(":")[0]
            idx.setdefault(base, []).append(r)
        return idx

    cdr_idx = index_by_id(common.get("cdr"))
    dev_idx = index_by_id(common.get("developability"))
    hum_idx = index_by_id(common.get("humanness"))

    rows = []
    for r in recs:
        rid = r["id"]
        cdr = (cdr_idx.get(rid) or [None])[0]
        dev = (dev_idx.get(rid) or [None])[0]
        hum = (hum_idx.get(rid) or [None])[0]
        nearest = (hum or {}).get("nearest_germline") or {}
        props = (dev or {}).get("properties") or {}
        summ = (dev or {}).get("summary") or {}
        rows.append({
            "id": rid,
            "track": r["track"],
            "sequence": r.get("sequence"),
            "sequence_length": len(r["sequence"]) if r.get("sequence") else None,
            "structure_file": r.get("structure_file"),
            # --- 공통 비교축 (동일 방법으로 계산됨) ---
            "common_metrics": {
                "chain_type": (cdr or {}).get("chain_type"),
                "cdrs": (cdr or {}).get("cdrs"),
                "cdr_lengths": (cdr or {}).get("cdr_lengths"),
                "cdr_method": (cdr or {}).get("method"),
                "molecular_weight_Da": props.get("molecular_weight_Da"),
                "isoelectric_point_pI": props.get("isoelectric_point_pI"),
                "gravy_kyte_doolittle": props.get("gravy_kyte_doolittle"),
                "instability_index": props.get("instability_index"),
                "n_liability_hits": summ.get("n_liability_hits"),
                "n_liabilities_in_cdr": summ.get("n_in_cdr"),
                "n_hydrophobic_patches": summ.get("n_hydrophobic_patches"),
                "germline_identity_percent": nearest.get("germline_identity_percent"),
                "nearest_germline": nearest.get("germline_gene"),
            },
            # --- 경로 고유 지표 (원본 이름 유지, 교차 비교 금지) ---
            "track_specific_metrics": r["track_specific_metrics"],
        })

    n_a = sum(1 for r in rows if r["track"].startswith("A"))
    n_b = sum(1 for r in rows if r["track"].startswith("B"))

    comparison = {
        "track_a": {"source_file": a.track_a, "load_error": err_a, "n_designs": n_a,
                    "confidence_metrics": ["iptm", "distogram_iptm_proxy",
                                           "cdr_distogram_iptm_proxy", "final_loss"]},
        "track_b": {"source_file": a.track_b, "load_error": err_b, "n_designs": n_b,
                    "confidence_metrics": ["RF2 score 파일 컬럼 (pAE 등)",
                                           "ca_rmsd_self_consistency_A"]},
        "common_axis": ["chain_type", "cdr_lengths", "molecular_weight_Da",
                        "isoelectric_point_pI", "gravy_kyte_doolittle",
                        "instability_index", "n_liability_hits",
                        "n_liabilities_in_cdr", "germline_identity_percent"],
        "designs": rows,
        "interpretation_rules": [
            "경로 A 의 ipTM 과 경로 B 의 pAE 는 서로 다른 모델·다른 척도다. "
            "하나의 종합 점수로 합치지 말 것 — 그런 통합 지표는 근거 없는 발명이다.",
            "공통 비교축(서열 기반)만 두 경로 간 직접 비교가 정당하다.",
            "germline identity 와 liability 플래그는 '개발 가능성 위험 신호'이지 "
            "결합력·면역원성 예측이 아니다.",
            "어느 경로도 wet-lab 검증을 대체하지 않는다. 실측(발현·SPR/BLI·DSF·SEC)만이 "
            "hit 여부를 결정한다.",
            "RFantibody 저자는 신뢰할 만한 필터의 부재를 최대 한계로 명시했고, "
            "일반적으로 10k 규모 캠페인을 예상한다 — 소수 설계의 순위는 약한 증거다.",
        ],
        "evaluation_errors": eval_errors,
    }

    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8") as fh:
            json.dump(comparison, fh, ensure_ascii=False, indent=2)

    checks = [
        ("최소 한 경로의 설계 결과 로드", (n_a + n_b) >= 1),
        ("서열이 있는 설계 ≥ 1 (공통 지표 계산 가능)", len(with_seq) >= 1),
        ("공통 평가 3종 모두 성공 (cdr/developability/humanness)",
         len(common) == 3 and not eval_errors),
    ]
    notes = (f"경로 A {n_a}건, 경로 B {n_b}건 비교. "
             f"공통 지표는 서열 기반 3종(cdr_analysis / developability / humanness)만 사용. "
             f"경로별 신뢰도 지표는 원본 이름 그대로 병기하며 통합 점수를 만들지 않는다. "
             + (f"평가 오류: {eval_errors}" if eval_errors else "")
             + (f" track-a 로드: {err_a}." if err_a else "")
             + (f" track-b 로드: {err_b}." if err_b else ""))
    emit(make_result(comparison, "compare_designs (track A vs B, common sequence axis)",
                     f"track_a={a.track_a} track_b={a.track_b}", checks, notes=notes))
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        emit(make_result(None, "compare_designs (CRASH)", " ".join(sys.argv[1:]),
                         [("스크립트 정상 종료", False)],
                         notes=f"예기치 못한 오류: {type(exc).__name__}: {exc}"))
        raise SystemExit(1)
