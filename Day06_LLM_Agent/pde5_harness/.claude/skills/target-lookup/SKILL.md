---
name: target-lookup
description: PDE5A 표적의 실재·기능을 UniProt(O76074)로 확인한다. 활성물질 조회 전 표적 검증 단계.
---
# target-lookup

## 언제 사용하는가
PLAN 직후, 첫 실행 단계(EXECUTE (a)). 후속 단계(ChEMBL 조회 등)로 진행하기 전에
표적 PDE5A가 실재하고 기능이 확인되는지 검증한다.

## 어떻게 (호출 명령)
```bash
python scripts/target_lookup.py            # 기본 accession O76074
python scripts/target_lookup.py O76074     # 명시적
```
- UniProt REST(`rest.uniprot.org`)에서 O76074를 조회 → gene(PDE5A)/protein/function/length 추출.
- 네트워크 불가 시 오프라인 데모 폴백(공개 사실, "실데이터 아님" 표기)으로 진행.

## 반환 검증
표준 봉투 `{result, provenance, verification}` 를 반환한다. 다음을 확인:
- `verification.passed == true` — 아니면 다음 단계 진행 금지(최대 2회 재시도 후 "미확인" 플래그).
- checks: `accession == O76074`, `phosphodiesterase 5 / PDE5A 언급`, `기능 서술 존재`.
- `provenance.source` 가 UniProt(또는 OFFLINE DEMO 명시)인지 확인.
- 표적 사실: PDE5A = cGMP 분해효소, 억제 시 cGMP↑ → NO-cGMP 평활근 이완·혈관확장.
