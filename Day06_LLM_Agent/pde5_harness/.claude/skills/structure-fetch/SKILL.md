---
name: structure-fetch
description: 표적의 실험 결정 구조를 RCSB 에서 가져온다. 공결정 리간드가 있는 구조를 고르는 것이 중요하다.
---
# structure-fetch

## 언제 사용하는가
도킹 전. 수용체 구조와 **결합 부위를 정의할 공결정 리간드**가 함께 필요하다.

## 어떻게
```bash
# UniProt accession 으로 구조 검색
curl -s "https://search.rcsb.org/rcsbsearch/v2/query?json=..." # accession 매칭
curl -sL "https://files.rcsb.org/download/<PDBID>.pdb" -o sample_run/structures/<PDBID>.pdb
```

## 무엇을 보고 고르는가
- **공결정 리간드 유무** — 있으면 결합 부위를 좌표로 정의할 수 있고 재도킹 대조가 가능하다
- 해상도 — 낮을수록 좋다
- 결손 잔기 — 결합 부위 근처에 있으면 피한다
- 촉매 금속 — PDE 계열은 Zn/Mg 를 남겨야 한다

## 반환 검증
- 수용체 원자 수가 기대 범위인가
- 공결정 리간드가 실제로 존재하는가 (HETATM 잔기명 확인)
- 물은 제거했고 금속은 남겼는가
