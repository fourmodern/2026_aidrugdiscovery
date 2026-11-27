# Claude Science — 2026년 11월 확인 사실

> 11/28 실습(수강생 전원 Claude Science)을 설계하기 위해 **공식 문서만** 보고 정리했다.
> 확인일 **2026-11-27**. 문서에 없는 것은 "미확인"으로 적었다. 추정·날조 금지.
> 출처는 모두 `https://claude.com/docs/claude-science/<page>.md` 또는 명시한 URL.

---

## 1. 접근 조건 · 설치

| 항목 | 사실 | 출처 |
|---|---|---|
| 정체 | 데스크톱 앱(웹사이트 아님). Claude 데스크톱 앱과 **별도 설치**. 공개 URL 없음 | overview, get-started |
| 베타 | **macOS · Windows · Linux 베타** | overview |
| 플랜 | Pro · Max · Ultra · Team · Enterprise 포함, **별도 구매 없음**. **Free 플랜 불가** | overview §Plans and usage |
| Team/Enterprise | Owner가 조직 단위로 **먼저 활성화**해야 함. `Claude Team plan for scientists`는 기본 on | overview, enable-claude-science |
| OS 요건 | macOS 13+ (Apple silicon/Intel), **Windows 11 (x64)**, Linux x64 glibc | overview §Requirements |
| 디스크 | 런타임·기본 환경에 **약 5 GB** | overview |
| 지역 | claude.ai 제공 국가 (한국 포함) | overview |
| 설치 파일 | https://claude.com/product/claude-science | get-started |
| 로그인 | Claude 계정으로 웹 로그인. **API 키 불필요** | get-started |
| 데이터 위치 | `%USERPROFILE%\.claude-science` (Win) / `~/.claude-science` (mac·Linux) | get-started |
| 최신 버전 | **0.1.75 (2026-11-24)** | changelog |

### OS별 주의 (실습 당일 사고 예방)
- **Windows**: 관리자 권한 없이 사용자 계정에 설치. 첫 실행에 **앱 윈도 엔진 약 150 MB를 `downloads.claude.ai`에서 내려받음** → 사내망에서 막히면 "couldn't set up the app window engine" 오류. 파이썬·R 셀은 추가 권한 없이 동작하고, git·cmd 스크립트만 1회 Windows 권한이 필요. 프록시·사내 루트 인증서는 Windows 설정을 그대로 따름.
- **Windows 제한**: 폴더 권한은 **로컬 NTFS/ReFS만** — 네트워크 드라이브·USB(FAT/exFAT) 폴더는 지정 불가(복사해서 쓰기). **Settings > Compute의 Model endpoints(NVIDIA BioNeMo NIM) 없음**.
- **macOS**: Intel/Apple Silicon 빌드를 **골라서** 받아야 함(잘못 받으면 "지원되지 않는 앱"). 브라우저 탭으로 열리며, 업데이트 후 자체 창으로 바뀔 수 있음.
- **Linux**: `bubblewrap 0.8.0+`와 `socat`을 **시스템 폴더에** sudo로 먼저 설치. Ubuntu 22.04 저장소의 bubblewrap은 **버전 미달**, 24.04는 충족.
- 네트워크 허용 도메인은 `network-requirements.md` 참조(사내망이면 IT에 미리 요청).

---

## 2. 작업 방식 · 재현성

- **프로젝트**가 세션과 산출물(artifact)을 묶는 단위. 결과는 **버전이 붙은 artifact**로 저장되고 provenance 기록이 남는다. (core-concepts, artifacts)
- Claude가 **Python · R · shell**을 샌드박스 안에서 직접 작성·실행. 커널은 세션 중 변수를 유지하고 **약 30분 유휴 시 종료**. (tools-and-environments)
- 기본 환경 2개(읽기 전용): Python(numpy·pandas·scipy·matplotlib·seaborn·pillow), R(tidyverse·ggplot2·jsonlite). 추가 패키지는 **task environment**를 만들어 영구 설치. (tools-and-environments)
- 패키지 출처: **conda-forge·bioconda·defaults·pytorch(micromamba), PyPI(pip), CRAN·Bioconductor**. 패키지가 없으면 소스를 받아 샌드박스에서 빌드. `apt`·`sudo`·root 없음. (tools-and-environments)
- pip는 **휠(prebuilt)만** 설치. git URL·로컬 폴더 설치 불가 → 없으면 conda-forge를 쓰거나 셀 안에서 임시 설치.
- **샌드박스**: 네트워크는 기본 차단(deny-by-default). 패키지 매니저와 **기본 커넥터가 쓰는 과학 DB**, 그리고 사용자가 승인한 호스트만 허용. 파일은 승인한 폴더를 **읽기 전용**(mac·Linux)으로만 접근. (core-concepts §Sandbox)
- **reviewer(검토 에이전트)**: 응답·승인된 계획·artifact·실행 기록을 다시 읽고 **주장이 실제 실행과 맞는지** 검사. 계산을 다시 돌리지는 않는다. 자동 검토는 **Max·Ultra·Team·Enterprise는 기본 on, Pro는 기본 off**(세션별로 켜기). `Request review`로 수동 실행. (the-reviewer)
- **auto mode**: 켜면 매 단계 승인 카드 없이 코드 실행·커넥터 호출을 진행. (auto-mode)
- **GPU**: Linux + NVIDIA에서 `Settings > Compute`의 Host GPU 스위치로만 사용. **Windows 노트북에서 로컬 GPU 경로는 문서에 없음.** 원격은 SSH(워크스테이션·Slurm) 또는 **본인 Modal 계정**(Modal이 직접 과금). (tools-and-environments, remote-compute-clusters, compute-providers)

---

## 3. 커넥터 · 스킬 (실습에 쓰는 것만)

**기본 제공(Featured, 기본 on)** — connectors-and-skills
- 문헌·지원: **Literature Graph (OpenAlex·arXiv)**, Research Resources
- 화학: **PubChem · ChEBI · Rhea · BindingDB**, 추가로 **ZINC**(구매 가능 화학공간), **Ketcher**(2D 분자 스케치)
- 구조: **PDB · AlphaFold · EMDB · Complex Portal · IntAct**
- 단백질 주석: InterPro·Pfam·Human Protein Atlas·STRING / 유전자: MyGene·UniProt·GO·Reactome
- 규제: FDA drug data · openFDA

**디렉터리 커넥터(별도 추가 필요)**: **PubMed · ClinicalTrials · ChEMBL · bioRxiv**
→ `Settings > Connectors > Add connector > Browse Connectors Directory`. **Team/Enterprise는 조직에 추가돼야 목록에 나온다** → 실습 전에 반드시 준비.

**기본 과학 스킬**: literature review, indication dossier, 그리고 모델별 스킬 — **AlphaFold2 · Boltz-2 · Chai-1 · ESMFold2 · OpenFold3 · ProteinMPNN(LigandMPNN·SolubleMPNN) · DiffDock · ESM-2 · Evo 2 · Borzoi · scGPT · scvi-tools**
- 정렬(MSA)은 공개 **ColabFold 서버(api.colabfold.com)** 를 쓸 수 있고, 이때 서열은 Anthropic을 거치지 않고 해당 서버로 **직접** 전송된다 → 사내 서열은 주의.
- `/`를 입력해 스킬을 직접 지정 가능. 커스텀 스킬은 **Chat with Claude · Write from scratch · Upload · Import from GitHub**으로 추가(GitHub 가져오기는 자동 업데이트 안 됨).
- **OpenAlex는 요청마다 키가 필요** → `Settings > Credentials`에 무료 OpenAlex API 키를 넣어야 OpenAlex 기반 문헌 검색이 동작. (literature-access)
- 문헌 전문: open access는 키 없이(Unpaywall·Semantic Scholar·PMC), 유료는 출판사 키나 기관 EZproxy를 넣어야 접근. NCBI 키를 넣으면 PubMed 속도 3→10 req/s. (literature-access)

### 도킹 관련 — 실습 설계의 핵심 판단
- **DiffDock 스킬은 기본 제공**이지만 딥러닝 도킹이라 실질적으로 GPU가 필요하다. **수강생 Windows 노트북에는 GPU 경로가 없다**(위 §2).
- 따라서 **고전 도킹(smina)** 을 샌드박스에 conda로 설치해 CPU로 돌리는 경로를 쓴다. conda-forge에 **win-64 · osx-arm64 · linux-64 빌드 모두 존재(smina 2020.12.10)** — 2026-11-27 anaconda.org API로 확인.
  - `conda-forge/vina`에는 **win-64 빌드가 없다**(linux·osx만) → Windows 실습에는 vina 대신 **smina**.
  - smina는 수용체 PDB와 리간드 SDF를 직접 읽어 pdbqt 변환이 불필요해 수업에 유리하다.
  - `conda-forge/rdkit`은 win-64 2026.09.3까지 있음.
- 이 선택은 우리 PDE5 하네스(smina, exhaustiveness 64, PDB 1UDT)와 동일해 이론(3·4교시)→실습 연결이 유지된다.

---

## 4. 사용량 · 비용 · 데이터

- Claude Science 사용량은 **claude.ai·Claude Code·Cowork와 같은 플랜 한도를 공유**. 앱의 `Settings > Usage`에서 확인. (overview)
- 한도는 **5시간 단위 + 주간** 두 가지. 한도에 닿으면 세션이 대기하거나(리셋 8시간 이내) 멈춘다. 대화·artifact는 보존. (overview)
- 추가 비용이 나는 경우: ① usage credits/extra usage를 켠 경우 ② 사용량 기반 Enterprise ③ **본인 계정으로 연결한 외부 서비스(Modal 등 — 직접 과금)**.
- 일부 모델(Claude Fable 계열)은 Pro·일반 Team 좌석에서 **첫 메시지부터 usage credits로만** 동작 → 0.1.73+는 세션 첫 메시지 전에 물어본다.
- 데이터: 파일은 내 컴퓨터에 남고 코드는 샌드박스에서 실행. 무엇이 Anthropic으로 가고 얼마나 보관되는지, 학습 미사용·소유권 약정은 `how-claude-science-works-with-your-data.md` 참조. **회사 데이터 업로드는 사내 규정 확인 후** 진행.
- **안전장치**: 생명과학 관련 내용은 광범위한 safeguards에 걸릴 수 있고, 걸리면 세션이 멈추거나 **자동으로 다른 모델로 재시도**(기본 on). Team/Enterprise는 **Life Sciences Verification Program(베타)** 으로 완화 가능. (safeguards)

---

## 5. 2026년 9~11월 변경점 (changelog)

| 버전 | 날짜 | 실습에 영향 있는 변경 |
|---|---|---|
| 0.1.75 | 11-24 | mac·Linux에서 Claude가 다시 git 커밋 가능(자체 사본에 커밋 후 pull) |
| 0.1.73 | 11-18 | Read & write 폴더를 직접 고치지 않고 **앱 자체 폴더에서 작업 후 반영**. 폴더 승인 카드가 **Read-only에서 시작**. literature review 스킬의 철회 논문 탐지 강화. Pro에서 credits 모델 사전 확인. ESMFold2는 **NVIDIA 드라이버 580+** 필요 |
| 0.1.72 | 11-05 | LaTeX 미리보기 개선. mac에서 자체 창으로 전환 가능. Linux 디스크 약 5 GB 추가 |
| 0.1.71 | 10-27 | 한도 도달 시 credits 소비 전 확인. Linux 샌드박스는 **시스템 폴더의 bwrap·socat만** 사용 |

Windows 지원 충돌에 대한 판정: **공식 `overview.md`·`get-started.md`가 Windows 앱(Windows 11 x64, 자체 창, 서명된 설치 파일)을 명시**하므로 Windows 네이티브 앱은 **있다**. `run-on-windows-wsl` 경로는 "Windows에서 **Linux 판**을 쓰고 싶을 때"의 대안이고, "네이티브 빌드가 없다"는 3자 저장소 설명은 낡았다.

---

## 6. 11/28 실습에 그대로 쓰는 결론

1. 수강생 준비: **Pro 이상 계정 + 앱 설치 + 로그인**까지 **사전 과제**로 끝내 둔다(첫 실행 5 GB·150 MB 다운로드 때문에 수업 중 설치는 위험).
2. 수업 시작 전 확인 3개: ① `Settings > Connectors`에 **ChEMBL·PubMed 디렉터리 커넥터 추가** ② `Settings > Credentials`에 **OpenAlex 키**(무료) ③ `Settings > General`에 **contact email**(Unpaywall용).
3. Pro 플랜은 **자동 검토가 꺼져 있다** → 실습 1단계에서 **auto-review를 켜게** 안내(무-날조 교육의 핵심 장치).
4. 도킹은 **smina(CPU) + PDB 1UDT**. DiffDock·Boltz-2는 "GPU 있으면"으로 선택 과제.
5. 사내 데이터는 올리지 않는다 — 전 과정을 **공개 DB(PDB·ChEMBL·PubChem·BindingDB·PubMed)** 로만 수행.
6. 한도 관리: 5시간 한도를 공유하므로 수강생에게 **불필요한 재실행을 줄이고 결과를 artifact로 저장**하게 안내.
