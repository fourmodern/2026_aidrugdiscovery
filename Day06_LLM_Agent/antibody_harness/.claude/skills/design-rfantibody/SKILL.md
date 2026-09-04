---
name: design-rfantibody
description: 경로 B (고전 비교군) — RFdiffusion → ProteinMPNN → RoseTTAFold2 3단계로 항체를 설계하고 저자 권장 필터를 적용한다. GPU 필요.
---
# design-rfantibody

## 언제 사용하는가
경로 A(ESMFold2 inversion)와 **비교할 고전 baseline** 이 필요할 때.
"2024년 방식 vs 2026년 방식" 대비를 보여주는 것이 이 경로의 교육적 목적이다.

## 3단계
1. **RFdiffusion** (항체 튜닝 버전) — 항원 hotspot 과 CDR 길이를 지정해 백본 생성
2. **ProteinMPNN** — 그 백본에 맞는 서열 설계
3. **RoseTTAFold2** — 구조 재예측 + pLDDT/PAE 로 필터

## 출처 (모두 확인됨)
- 저장소: <https://github.com/RosettaCommons/RFantibody> (MIT, default branch `main`)
- 가중치: `bash include/download_weights.sh`
- 설치: `uv sync` / Docker / Apptainer. **NVIDIA GPU + CUDA 11.8+**, Ubuntu 22.04 권장
- 프리프린트: <https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1>

## 어떻게 (호출 명령)
```bash
# GPU 없이 명령어·입력 형식만 검증
python scripts/design_rfantibody.py \
    --target antigen.pdb --framework framework_HLT.pdb --out outputs/design_rfab \
    --loops "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11" --hotspots "T570,T593" \
    -n 10 --dry-run

# GPU 노드에서 실제 실행
python scripts/design_rfantibody.py --target antigen.pdb --framework framework_HLT.pdb \
    --hotspots "T570,T593" -n 100 --out outputs/design_rfab \
    --self-consistency > outputs/design_b.json
```

스크립트가 구성해 실행하는 명령 (공식 README 인자 그대로):
```
rfdiffusion  -t <antigen.pdb> -f <HLT.pdb> --output-quiver 1_rfdiffusion.qv -n N
             --design-loops "H1:7,..." --hotspots "T305,T456" [--deterministic]
proteinmpnn  --input-quiver 1_rfdiffusion.qv --output-quiver 2_proteinmpnn.qv
             --seqs-per-struct N --temperature 0.1 [--loops "H1,H2,H3,L1,L2,L3"]
rf2          --input-quiver 2_proteinmpnn.qv --output-quiver 3_rf2.qv
             --num-recycles 10 --seed S
qvscorefile  3_rf2.qv
```

## HLT 입력 포맷 (공식 규격)
PDB 변형. 체인 ID **Heavy=`H`, Light=`L`, Target=`T`**, 체인 순서 **H→L→T**.
CDR 은 PDB remark 로 1-indexed 절대 잔기번호 주석:
```
REMARK PDBinfo-LABEL:   32 H1
REMARK PDBinfo-LABEL:   52 H2
```
`design_rfantibody.py` 가 이 형식(체인 목록·순서·CDR 라벨)을 실제로 검사하고
위반 시 `input_checks.framework_hlt.issues` 에 기록한다.

## 필터 (저자 문서 값 — 발명 금지)
`RF2 pAE < 10`, `RMSD(design vs RF2 예측) < 2 Å`, (선택) `Rosetta ddG < -20`.
score 파일에 해당 컬럼이 없으면 `applied: false` + 실제 컬럼 목록을 보고한다 (추측하지 않는다).

## self-consistency
`--self-consistency` 로 설계 백본 vs RF2 예측의 **Cα RMSD** 를 계산한다
(BioPython `Superimposer`, Kabsch). `qvextract` 로 Quiver 를 PDB 로 풀어
`<out>/extracted_rfdiffusion/` 과 `<out>/extracted_rf2/` 에 둔 뒤 실행해야 한다.
디렉토리가 없으면 RMSD 를 추정하지 않고 `skipped_reason` 을 남긴다.

## 저자가 명시한 한계 — 반드시 보고서에 인용
> "The lack of an effective filter is the main limitation of the RFantibody pipeline at the moment."

- 일부 타깃은 95 designs 로 VHH binder 를 찾았으나, 일반적으로 **10k 규모** 캠페인 예상.
- hotspot 선택에 vanilla RFdiffusion 보다 민감 → 대규모 전 **pilot run** 권장.
- 하전·극성 부위, 글라이칸 인접부, 비구조화 루프는 여전히 어렵다.
- 런타임 O(N²) → 큰 항원은 절단(truncation) 필요.

**소수 설계의 순위는 약한 증거다.** 이 문장을 결론에 반드시 포함하라.

## 무-날조
CLI 3종이 없거나 GPU 미검출이면 `mode: aborted` + `passed: false` 로 보고하고
설계 서열·pAE·RMSD 를 생성하지 않는다. 실행 절차는 `RUNPOD_가이드.md`.
