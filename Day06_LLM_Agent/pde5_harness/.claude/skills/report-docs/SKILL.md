---
name: report-docs
description: 보고서 마크다운을 Word(.docx)와 발표용 PowerPoint(.pptx)로 내보낸다. 그림과 표를 실제로 삽입한다.
---
# report-docs

## 언제 사용하는가
`galley-proof` 통과 후 마지막 단계.

## 어떻게 (호출 명령)
```bash
python scripts/make_docs.py --md outputs/report_pde5.md \
    --figures outputs/figures --out outputs/docs
```

## 산출물
| 파일 | 내용 |
|------|------|
| `report_pde5.docx` | IMRAD 본문 + 그림 5장 + 표 3종 실제 삽입 |
| `report_pde5.pptx` | 발표용 12장 내외, 표지에 graphical abstract |

## 하드 규칙
- 그림은 **링크가 아니라 실제 이미지로 삽입**한다 (`add_picture`).
- 표는 텍스트 나열이 아니라 실제 표 객체로 넣는다 (`add_table`).
- 변환 과정에서 수치를 재계산하거나 반올림해 바꾸지 않는다.

## 반환 검증
- docx 안에 이미지 5개, 표 3개 존재.
- pptx 슬라이드 수 ≥ 10, 표지에 그림 1장.
