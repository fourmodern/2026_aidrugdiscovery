# RunPod 실행 가이드 — 항체 설계 하네스 (GPU 단계)

> ## ⚠️ **비용 경고 — 작업이 끝나면 반드시 pod 을 종료하라**
> ```bash
> python <RUNPOD_CTL> rm <POD_ID>
> ```
> **`rm` 을 하지 않으면 pod 이 켜져 있는 동안 계속 과금된다.** A100 80GB 급은 시간당 과금이
> 크므로, 설계가 끝나면 `get` 으로 결과를 회수한 즉시 `rm` 한다.
> `list` 로 살아 있는 pod 이 없는지 항상 마지막에 확인할 것.

---

## 0. 원칙 — RunPod 제어 코드를 새로 만들지 않는다

이 하네스는 RunPod API 를 직접 호출하지 않는다. 이미 검증된 헬퍼를 **그대로 재사용**한다:

```
/home/hjpark/foundation_model_research/projects/_shared_infra/runpod_ctl.py
```

이 도구의 방침(파일 docstring 인용): **"ad-hoc curl 금지 — pod 배포·상태·대기·ssh·전송·종료를
이 도구로만 한다."** REST 엔드포인트는 `https://rest.runpod.io/v1`, SSH 키는 `~/.ssh/id_ed25519`
(배포 시 `.pub` 가 `PUBLIC_KEY` 환경변수로 주입된다).

편의상 아래에서는 이렇게 줄여 쓴다:

```bash
RUNPOD_CTL=/home/hjpark/foundation_model_research/projects/_shared_infra/runpod_ctl.py
HARNESS=/home/hjpark/lecture_drug1/0905_agent/antibody_harness
```

---

## 1. 환경변수 (이름만 — 값은 절대 문서·로그·코드에 남기지 말 것)

| 변수 | 용도 |
|------|------|
| `RUNPOD_API_KEY` | RunPod REST 인증. `runpod_ctl.py` 가 환경변수 또는 `.env` 에서 읽는다 |
| `HF_TOKEN` | HuggingFace 가중치 다운로드 (`biohub/ESMC-6B`, `biohub/ESMFold2`) |

```bash
export RUNPOD_API_KEY=...      # 값은 여기 적지 않는다
export HF_TOKEN=...            # 값은 여기 적지 않는다
```

`/home/hjpark/foundation_model_research/.env` 에 정의되어 있다.

> 🔒 **이 하네스는 공개 GitHub 리포에 푸시된다.** `.env` 파일이나 키 값을 복사·기록·출력하지 말 것.
> `.gitignore` 에 `.env` 가 포함되어 있다. 문서에는 **변수 이름과 사용법만** 적는다.

---

## 2. 권장 GPU

| 경로 | 권장 GPU | 근거 |
|------|----------|------|
| **A. ESMFold2 inversion** | **A100 80GB** (최소 A100 40GB) | 공식 프로토콜 `binder_design.py` 주석: `target_name="cd45", binder_name="trastuzumab_framework_vhvl", batch_size=1` 기준 **`REUSE_ESMC=True` 에서 27GB, False 에서 51GB VRAM**. 80GB 면 batch_size 를 최대 6까지 올릴 수 있다고 같은 주석에 명시 |
| **B. RFantibody** | CUDA 11.8+ NVIDIA GPU | 공식 README 요구사항. VRAM 요구치는 문서에 수치로 명시되어 있지 않음 — **미확인** |

경로 A 는 ESMC-6B(60억 파라미터) 를 포함하므로 24GB 급(RTX 3090/4090)에서는
`REUSE_ESMC=True` + `--no-scaling-critics` 로도 빠듯하다. **A100 이상을 권장한다.**

> `--gpu` 인자에는 RunPod 의 GPU 타입 ID 문자열을 쉼표로 나열한다
> (`runpod_ctl.py` 기본값은 `"NVIDIA RTX A4000"`). 정확한 ID 문자열은 RunPod 콘솔에서 확인할 것 —
> 이 문서에 A100 의 ID 문자열을 추측해 적지 않는다.

---

## 3. 표준 흐름

### (1) 배포

```bash
python $RUNPOD_CTL deploy --name ab-design \
    --gpu "<RunPod A100 GPU type ID>" \
    --disk 200 \
    --image runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04
# → 마지막 줄에 POD_ID 가 출력된다
```

`runpod_ctl.py deploy` 인자 (파일에서 직접 확인):

| 인자 | 기본값 | 설명 |
|------|--------|------|
| `--name` | (필수) | pod 이름 |
| `--gpu` | `NVIDIA RTX A4000` | GPU 타입 ID, 쉼표 구분 다중 지정 가능 |
| `--cpu` | None | CPU flavor (지정 시 CPU pod) |
| `--vcpu` | 4 | CPU pod 의 vCPU 수 |
| `--dc` | None | 데이터센터 ID |
| `--volume` | None | network volume ID |
| `--mount` | `/workspace` | 볼륨 마운트 경로 |
| `--image` | `runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04` | 컨테이너 이미지 |
| `--disk` | 60 | 컨테이너 디스크 GB |

배포는 `cloudType: SECURE`, `ports: ["22/tcp"]`, `supportPublicIp: true` 로 고정되어 있다.

### (2) 대기 — 인내심 설계를 이해할 것

```bash
python $RUNPOD_CTL wait <POD_ID> --timeout 1500
```

`wait` 는 **최대 25분(기본 1500초)** 동안 30초 간격으로 폴링한다. 중요한 점:

- **`runtime` 이 `null` 이어도 계속 SSH 를 시도한다.** RunPod API 의 `runtime` 필드는
  실제 SSH 준비보다 늦게 채워지는 경우가 있다.
- 판정 기준은 API 상태가 아니라 **`ssh ... echo READY` 가 실제로 성공하는지**다.
- 10GB 급 이미지 pull 이 느릴 수 있어서 기본 타임아웃이 길게 잡혀 있다.
  진행 중에는 `[{초}s] {ip}:{port} runtime={up|null} — 대기중` 로그가 찍힌다.

즉, 몇 분간 `runtime=null` 이 계속 나와도 **정상이며 기다려야 한다.**

### (3) 하네스 전송

```bash
tar czf /tmp/ab_harness.tgz -C $(dirname $HARNESS) antibody_harness
python $RUNPOD_CTL put <POD_ID> /tmp/ab_harness.tgz /workspace/ab_harness.tgz
python $RUNPOD_CTL ssh <POD_ID> "cd /workspace && tar xzf ab_harness.tgz && ls antibody_harness"
```

`put`/`get` 은 내부적으로 `scp` 를 쓰며 기본 타임아웃 1800초다 (`--timeout` 으로 조정).

### (4) 가중치 다운로드 + 설계 실행

**경로 A — ESMFold2 inversion**

```bash
python $RUNPOD_CTL ssh <POD_ID> "
  export HF_TOKEN='$HF_TOKEN' &&
  cd /workspace/antibody_harness &&
  pip install -q 'esm@git+https://github.com/Biohub/esm.git@main' abnumber requests biopython &&
  python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder --dry-run
" --timeout 3600
```

dry-run 이 통과하면 실제 설계:

```bash
python $RUNPOD_CTL ssh <POD_ID> "
  export HF_TOKEN='$HF_TOKEN' &&
  cd /workspace/antibody_harness &&
  python scripts/design_esmfold2.py --target-name pd-l1 --binder-name minibinder \
      --seed 0 --num-seeds 4 --batch-size 1 \
      --out outputs/design_esmfold2 > outputs/design_a.json
" --timeout 21600
```

> `abnumber` 는 conda 패키지에 의존한다. 공식 `binder_design.py` 의 `main()` 에는
> *"'abnumber' will fail if running this script with uv run. It requires conda packages."*
> 라는 assert 주석이 있고, 로컬 실행 시 `use_scaling_critics=False` 를 요구한다.
> pip 설치가 실패하면 `--no-scaling-critics` 로 실행하라.

**경로 B — RFantibody**

```bash
python $RUNPOD_CTL ssh <POD_ID> "
  cd /workspace &&
  git clone https://github.com/RosettaCommons/RFantibody.git &&
  cd RFantibody &&
  bash include/download_weights.sh &&
  curl -LsSf https://astral.sh/uv/install.sh | sh &&
  \$HOME/.local/bin/uv sync
" --timeout 7200

python $RUNPOD_CTL ssh <POD_ID> "
  cd /workspace/RFantibody && source .venv/bin/activate &&
  cd /workspace/antibody_harness &&
  python scripts/design_rfantibody.py \
      --target inputs/antigen.pdb --framework inputs/framework_HLT.pdb \
      --loops 'H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11' --hotspots 'T570,T593' \
      -n 100 --out outputs/design_rfab > outputs/design_b.json
" --timeout 21600
```

**비교**

```bash
python $RUNPOD_CTL ssh <POD_ID> "
  cd /workspace/antibody_harness &&
  python scripts/compare_designs.py --track-a outputs/design_a.json \
      --track-b outputs/design_b.json --out outputs/06_comparison.json
" --timeout 1800
```

### (5) 결과 회수

```bash
python $RUNPOD_CTL ssh <POD_ID> "cd /workspace/antibody_harness && tar czf /workspace/results.tgz outputs"
python $RUNPOD_CTL get <POD_ID> /workspace/results.tgz /tmp/results.tgz
tar xzf /tmp/results.tgz -C $HARNESS
```

### (6) 🔴 종료 — 잊지 말 것

```bash
python $RUNPOD_CTL rm <POD_ID>
python $RUNPOD_CTL list          # 살아 있는 pod 이 없는지 반드시 확인
```

---

## 4. 스토리지 주의

- **가중치가 크다.** ESMC-6B(60억 파라미터)에 더해 프로토콜이 여러 critic 체크포인트
  (`ESMFold2-Experimental-Fast`, `-Cutoff2025`, `ESMFold2-Experimental` 등 4종 + 선택적
  scaling critic 15종)를 내려받는다. `--disk` 를 **넉넉히**(200GB 이상) 잡아라.
  컨테이너 디스크가 차면 다운로드가 조용히 실패한다.
- **반복 사용하면 network volume 을 써라.** 매번 다시 받으면 시간과 대역폭이 낭비된다.

  ```bash
  python $RUNPOD_CTL vol-list
  python $RUNPOD_CTL vol-create --name ab-weights --size 300 --dc <DATACENTER_ID>
  python $RUNPOD_CTL deploy --name ab-design --gpu "<GPU_ID>" \
      --volume <VOLUME_ID> --mount /workspace --disk 100
  ```

  볼륨을 붙였으면 `HF_HOME=/workspace/hf_cache` 를 설정해 가중치 캐시를 볼륨에 남긴다:

  ```bash
  python $RUNPOD_CTL ssh <POD_ID> "export HF_HOME=/workspace/hf_cache && ..."
  ```

  > `vol-delete` 는 되돌릴 수 없다. 안전장치로 `--confirm <VOLUME_ID>` 를 정확히 재입력해야 실행된다.
- **network volume 은 데이터센터에 묶인다.** 볼륨과 pod 의 `--dc` 가 같아야 한다.

---

## 5. 예상 실행 시간 — 개념 (수치는 미측정)

> **정직 고지**: 아래 표는 **개념적 안내**이며, 이 하네스로 **실측한 값이 아니다.**
> 공식 문서도 벽시계 시간을 명시하지 않는다. 실제 시간은 GPU 종류·배치 크기·타깃 길이·
> 네트워크 속도에 따라 크게 달라진다. **첫 실행 때 직접 측정해 이 표를 채워라.**

| 단계 | 비용 요인 | 측정값 |
|------|-----------|--------|
| pod 배포 + 이미지 pull | 이미지 크기(10GB 급), 데이터센터 | 미측정 (`wait` 기본 25분 여유) |
| HF 가중치 다운로드 | ESMC-6B + critic 다수, 볼륨 캐시 유무 | 미측정 |
| 경로 A 모델 로드 | critic 개수, `torch.compile` | 봉투의 `model_load_seconds` 에 기록됨 |
| 경로 A 설계 1 seed | 공식 프로토콜 `STEPS = 150` 최적화 스텝 + critic 폴딩 | 봉투의 `elapsed_s` 에 기록됨 |
| 경로 B RFdiffusion | `-n` 설계 수, 항원 길이(런타임 O(N²)) | 미측정 |
| 경로 B ProteinMPNN | 백본 수 × `--seqs-per-struct` | 미측정 |
| 경로 B RF2 | 설계 수 × `--num-recycles` | 미측정 |

`design_esmfold2.py` 는 `model_load_seconds` 와 설계별 `elapsed_s` 를 봉투에 기록한다.
경로 B 는 `steps[].elapsed_s` 에 단계별 실측 시간이 남는다. **첫 실행 후 그 값을 인용하라.**

**비용 개념**: RunPod 는 pod 이 켜져 있는 시간에 비례해 과금한다(`runpod_ctl.py list` 가
`costPerHr` 를 출력한다). 따라서 비용 = `costPerHr` × (배포~`rm` 까지의 시간). 가중치 다운로드
시간도 과금 대상이므로, 반복 실행 시 network volume 캐시가 실질적인 비용 절감 수단이다.

---

## 6. 강의 데모 팁

- 강의실에서는 **`python run_harness.py design-check` 로 dry-run 만** 보여주는 것으로 충분하다
  (GPU 불필요, 즉시 실행, 공식 프로토콜 sha256 검증까지 보여준다).
- 실제 설계는 **사전에 RunPod 에서 돌려 결과 JSON 을 받아두고**, 강의 중에는
  `compare_designs.py` 로 비교·해석하는 데 시간을 쓰는 편이 안전하다.
- 데모 중 pod 을 띄웠다면 **강의 종료 직후 `rm`** 을 체크리스트에 넣어라.

---

## 7. 문제 해결

| 증상 | 확인 |
|------|------|
| `wait` 가 계속 `runtime=null` | 정상일 수 있음 — 이미지 pull 중. 25분까지 기다린다 |
| `wait` TIMEOUT | `status <id>` 로 `desiredStatus`·IP·포트 확인, `list` 로 pod 존재 확인 |
| `ssh` "엔드포인트 없음" | pod 이 아직 준비 안 됨 → `wait` 먼저 |
| 가중치 다운로드 중단 | 디스크 부족 의심 → `ssh <id> "df -h"`, `--disk` 늘려 재배포 |
| `esm` import 실패 | PyPI 판이 아니라 `pip install 'esm@git+https://github.com/Biohub/esm.git@main'` 필요 |
| `abnumber` 설치 실패 | `--no-scaling-critics` 로 실행 (공식 스크립트도 로컬 실행 시 동일 제약) |
| CUDA OOM (경로 A) | `--batch-size 1`, `--no-scaling-critics`, 더 큰 VRAM GPU |
| RFantibody CLI 없음 | `uv sync` 후 `.venv` 활성화했는지 확인 (`rfdiffusion --help`) |
