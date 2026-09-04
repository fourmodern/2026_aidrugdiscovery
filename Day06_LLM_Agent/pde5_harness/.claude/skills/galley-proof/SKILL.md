---
name: galley-proof
description: 출판 직전 교정. 8개 카테고리 전부 CLEAR 여야 통과. 문장이 아니라 정합성을 본다.
---
# galley-proof

## 언제 사용하는가
`devils-advocate` 통과 후, docx/pptx 최종본을 내보내기 직전.

## 8개 카테고리
| key | 확인 |
|-----|------|
| `text_integrity` | 오타·띄어쓰기·용어 일관성 |
| `citation_reference` | 본문 인용과 참고문헌 1:1, 식별자 실재 |
| `figure_table` | 번호·캡션·본문 참조 순서 일치, 해상도 |
| `equation_math` | 수식 번호와 상호참조 |
| `metadata` | 제목·저자·일시·버전·라이선스 |
| `supplementary` | 부속 자료 경로 유효 |
| `license_copyright` | 데이터·이미지 출처와 라이선스 |
| `numeric_traceability` | **본문의 모든 수치가 도구 산출 파일에 존재** |

## 외부 평가 의무
`Agent(subagent_type="oh-my-claudecode:verifier")` 로 검증한다.
특히 `numeric_traceability` 는 본문 수치를 `run_stdout.json` 과 기계적으로 대조한다.

## 반환 형식
```json
{"per_category": {"text_integrity": "CLEAR", "...": "UNCLEAR"},
 "findings": ["..."], "pass": false}
```
