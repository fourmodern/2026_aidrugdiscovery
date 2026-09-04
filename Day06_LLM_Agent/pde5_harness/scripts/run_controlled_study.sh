#!/usr/bin/env bash
# 골격 통제 연구 전체 파이프라인 — 도킹 이후 단계를 순서대로 돌린다.
# 각 단계는 게이트가 FAIL 이면 비영(非零)으로 끝나지만, 파이프라인은 계속 진행하고
# 마지막에 어느 단계가 실패했는지 요약한다. 실패를 조용히 넘기지 않는다.
set -u
cd "$(dirname "$0")/.."
PY=/home/hjpark/lecture_drug1/.venv_gpu/bin/python
LOG=/tmp/claude-1001/t004
mkdir -p "$LOG"
declare -a FAILED=()

step () {                       # step <이름> <로그파일> <명령...>
  local name="$1" logf="$2"; shift 2
  echo "── $name ──"
  if "$@" > "$LOG/$logf" 2>&1; then
    echo "   OK   ($(wc -l < "$LOG/$logf") 줄)"
  else
    echo "   FAIL (exit $?) — $LOG/$logf"
    FAILED+=("$name")
  fi
  tail -4 "$LOG/$logf" | sed 's/^/     /'
}

# 도킹이 끝날 때까지 기다린다
while pgrep -f "dock_controlled.py" > /dev/null; do sleep 20; done
echo "도킹 완료. 후속 단계 시작 $(date -u +%H:%M:%S)"

step "탐색 깊이 스윕"   sweep.log      $PY scripts/exhaustiveness_sweep.py
step "통제 분석"        analyze.log    $PY scripts/analyze_controlled.py
step "항 분석"          terms.log      $PY scripts/terms_controlled.py
step "그림 생성"        figs.log       $PY scripts/figures_controlled.py
step "보고서"           report.log     $PY scripts/report_controlled.py

echo
if [ ${#FAILED[@]} -eq 0 ]; then
  echo "전 단계 통과"
else
  echo "실패 단계: ${FAILED[*]}"
fi
