---
name: cdr-analysis
description: 항체 가변영역에서 CDR-H1/H2/H3, L1/L2/L3 을 추출한다. anarci/abnumber 가 있으면 정식 번호매김, 없으면 문서화된 휴리스틱 근사.
---
# cdr-analysis

## 언제 사용하는가
antibody-search 또는 설계 경로가 항체 사슬 서열을 낸 뒤. **알려진 항체와 설계 항체 모두**에 적용한다.

## 어떻게 (호출 명령)
```bash
python scripts/antibody_search.py | python scripts/cdr_analysis.py --stdin
python scripts/cdr_analysis.py "EVQLVESGG..."          # 서열 직접
python scripts/cdr_analysis.py chains.fasta            # FASTA
python scripts/cdr_analysis.py --stdin --scheme kabat  # anarci 설치 시 scheme 지정
```

## 방법 (반드시 구분해서 보고할 것)
1. **정식 번호매김**: `abnumber` 또는 `anarci` 가 import 되면 IMGT/Kabat 번호로 추출.
   `method` 필드가 `abnumber/ANARCI (IMGT numbering)` 등으로 표시된다.
2. **휴리스틱 폴백**: 둘 다 없으면 보존 모티프 정규식으로 Kabat 근사 추출.
   `method` 가 `heuristic (approximate) — Kabat-like, conserved-motif regex` 로 표시된다.
   - CDR-H1: FR1 보존 Cys +9 ~ FR2 Trp(`W[VIAFGLM][RKQGNHS]Q`) 직전
   - CDR-H2: FR2 Trp +14 ~ FR3 시작 모티프(`[RK][FVLIATM][TSVA][IFMLV][ST]`) 직전
   - CDR-H3: FR3 끝 보존 Cys(`[YFHVA][YFHCVL]C`) +3 ~ FR4(`WG[..]G[TSA]`) 직전
   - CDR-L1: FR1 보존 Cys +1 ~ FR2 Trp 직전
   - CDR-L2: FR2 Trp +15 부터 7 잔기 (Kabat L50-56 고정 길이 **가정**)
   - CDR-L3: FR3 끝 보존 Cys +1 ~ FR4(`FG[..]G[TS]`) 직전

## 반환 검증
- checks: 가변영역 ≥1 / CDR 3종 완전 추출 ≥1 / CDR 길이 1-40 범위 / 스킵 0건.
- `warnings` 를 반드시 읽을 것 — CDR-L2 고정 길이 가정, Cys 폴백 사용, scFv 분할 등이 기록된다.
- 경계를 못 찾은 CDR 은 `null`. 채워 넣지 말 것.

## 무-날조 (중요)
`method` 에 `heuristic (approximate)` 가 있으면 보고서에도 **근사임을 그대로 명시**해야 한다.
**"IMGT 번호매김으로 추출했다"고 쓰면 안 된다.** 정확한 번호가 필요하면
`conda install -c bioconda abnumber` (또는 anarci) 후 재실행하고, 그때 method 가 바뀐 것을 확인한다.

## 참고 (개발 시 확인한 사실)
이 휴리스틱은 트라스투주맙(PDB 1N8Z)과 퍼투주맙(PDB 1S78) 사슬에서 Kabat CDR 6종을
문헌 공개값과 일치하게 재현했다. 그러나 모든 항체(비정형 프레임워크·삽입·비인간)에 대해
보장되지 않는다.
