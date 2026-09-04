---
name: design-esmfold2
description: 경로 A (주축) — ESMFold2 inversion 으로 타깃 서열에서 de novo 미니바인더/항체 scFv 를 설계한다. GPU 필요.
---
# design-esmfold2

## 언제 사용하는가
항원(에피토프) 서열이 정해진 뒤, **새 바인더를 설계할 때**. 이 하네스의 주축 설계 경로다.

## 핵심 아이디어
구조 예측 모델 ESMFold2 를 **역방향으로 미분(inversion)** 해서, 타깃에 잘 결합하는 바인더 서열을
gradient 로 최적화한다. 별도의 생성 모델 없이 folding model 하나로 설계한다.

## 출처 (모두 확인됨 — 추측 아님)
- 저장소: <https://github.com/Biohub/esm> (MIT)
- 프로토콜: `cookbook/tutorials/binder_design.py` (타깃 서열 → ranked binder 엔드투엔드)
- 가중치: `biohub/ESMFold2`, `biohub/ESMC-6B` (HuggingFace, MIT)
- 설치: `pip install 'esm@git+https://github.com/Biohub/esm.git@main' abnumber modal`
- 프리프린트: <https://www.biorxiv.org/content/10.64898/2026.06.03.729735>

**wet-lab 검증 (보고서에 인용 가능)**: 5개 타깃(EGFR·PDGFRβ·PD-L1·CTLA-4·CD45)에서
항체 포맷 hit rate **15–29%**, 미니바인더 **36–88%**, nM 친화도.
FoldBench 항체-항원 DockQ pass-rate 에서 **AF3 상회**.

## 어떻게 (호출 명령)
```bash
# GPU 없이 설정·환경·API 심볼만 검증 (강의실 노트북)
python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder --dry-run

# 커스텀 항원 (HER2 등) + 항체 프레임워크
python scripts/design_esmfold2.py --target-name her2_ecd4 \
    --target-sequence-file inputs/her2_epitope.fasta \
    --binder-name trastuzumab_framework_vhvl --hotspots T570,T593 --dry-run

# GPU 노드에서 실제 설계
python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder \
    --seed 0 --num-seeds 4 --out outputs/design_esmfold2 > outputs/design_a.json
```

## preset (공식 프로토콜에 하드코딩)
- **항원**: `cd45`, `ctla4`, `egfr`, `pd-l1`, `pdgfr`.
  preset 이면 `--target-sequence` 를 **주면 안 된다** (프로토콜이 ValueError 를 던진다).
- **바인더**: `minibinder`, `trastuzumab_framework_vhvl`,
  `atezolizumab_framework_vhvl`, `ocankitug_framework_vhvl`.
- preset 밖 항원은 `--target-sequence(-file)` 로 에피토프 서열을 직접 준다.

## 반환 검증
- **dry-run** checks: 입력 조합 유효 / 공식 프로토콜 확보 / 기대 API 심볼 존재.
  `protocol.matches_verified_sha256` 가 `false` 면 upstream 이 갱신된 것이니
  `vendor/binder_design.py` 를 직접 열어 API 변경 여부를 확인하라.
- **실행** checks: CUDA GPU / 프로토콜 로드 / 설계 ≥1 / 전 seed 완료 / 서열 비어있지 않음.
- `ranked` 는 `iptm` → `distogram_iptm_proxy` → `final_loss` 순 정렬. 모두 모델 반환값이다.

## 무-날조 (중요)
- 이 스크립트는 공식 `binder_design.py` 를 `vendor/` 로 내려받아 **그대로 import** 한다.
  API 를 재구현하거나 함수명을 추측하지 않는다.
- **GPU 없으면 설계하지 않는다.** `--dry-run` 없이 CUDA 미보유로 실행하면
  `mode: aborted` + `passed: false` 를 내고 설계 서열/점수를 생성하지 않는다.
  이 상태를 보고할 때는 **"RunPod GPU 에서 검증 필요"** 라고 명시하라.
- ipTM 은 **구조 신뢰도**이지 결합상수(KD)가 아니다. 설계 결과는 **가설**이며
  발현·SPR/BLI 실측만이 hit 여부를 결정한다.

## GPU 요구
공식 프로토콜 주석 기준: `batch_size=1` 에서 `REUSE_ESMC=True` 27GB / `False` 51GB VRAM.
**A100 40GB 이상 권장(80GB 이면 batch_size 6까지).** 실행 절차는 `RUNPOD_가이드.md`.
