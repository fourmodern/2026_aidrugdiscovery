# 실습용 과학 스킬 5종

11/28 실습(Claude Science)의 ⑥ 보고서·⑦ 리뷰 단계에서 쓰는 스킬 모음이다.

| 스킬 | 쓰는 단계 | 하는 일 |
|------|----------|---------|
| `scientific-writing` | ⑥ 보고서 | IMRAD 구조, 문단 서술(불릿 금지), 인용 형식 |
| `scientific-visualization` | ⑥ 보고서 | 출판용 그림 — 다중 패널, 색각 이상 안전 팔레트, 저널 규격 |
| `citation-management` | ⑥ 보고서 | 서지 검증, BibTeX 생성, DOI→인용 변환 |
| `peer-review` | ⑦ 리뷰 | 체크리스트 기반 심사, 방법론·보고 기준(CONSORT/STROBE) |
| `scientific-critical-thinking` | ⑦ 리뷰 | 근거 품질 평가, 편향·교란 식별 |

## 왜 여기에 두었나

원본은 [`K-Dense-AI/scientific-agent-skills`](https://github.com/K-Dense-AI/scientific-agent-skills)
(MIT, skill-author: K-Dense Inc.)이다. 그 저장소는 **약 471 MB** 라서 Claude Science 의
GitHub 가져오기 제한(저장소 다운로드 100 MB · 스킬 폴더 합계 100 MB)을 넘는다.
그래서 실습에 필요한 5종(합계 약 1.1 MB)만 추려 두었다.

## 수강생이 쓰는 법 — 두 가지

**A. 이 저장소에서 가져오기 (빠름)**

Claude Science → `Settings > Skills > Add skill > Import from GitHub` →
`fourmodern/2026_aidrugdiscovery` 를 지정한다. 가져온 뒤 목록에서 위 5종을 켜 둔다.
(같은 저장소의 `pde5_harness/.claude/skills/` 19종도 함께 보인다 — 강사 데모 하네스의 스킬이다.)

**B. 직접 만들기 (제한 없음 · 4교시 Agent Skills 실습)**

`Settings > Skills > Add skill > Chat with Claude` 로 "PDE5 가상 스크리닝 보고서" 스킬을
직접 만든다. 핸드북 ⑥-0 장의 프롬프트를 그대로 쓰면 된다. 가져오기가 막힌 환경(조직이
커스텀 스킬을 제한한 경우 등)에서도 이 경로가 열려 있다면 쓸 수 있다.

> 스킬은 Claude 가 필요할 때 스스로 불러오는 지침 파일이다. 대화창에서 `/` 를 누르면
> 직접 지정할 수도 있다. 4교시에서 다룬 Agent Skills 가 이것이다.

## 라이선스

각 스킬의 `SKILL.md` frontmatter 에 적힌 대로 **MIT**, skill-author **K-Dense Inc.**
원본 저장소의 라이선스를 따른다. 교육 목적의 부분 재배포이며, 내용은 수정하지 않았다.
