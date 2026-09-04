#!/usr/bin/env bash
# 3단계 — 커스텀 채점 독립 검증을 돌리고, 그림·보고서·문서를 최종 재생성한다.
set -u
cd "$(dirname "$0")/.."
PY=/home/hjpark/lecture_drug1/.venv_gpu/bin/python
LOG=/tmp/claude-1001/t004
declare -a FAILED=()
step () {
  local name="$1" logf="$2"; shift 2
  echo "── $name ──"
  if "$@" > "$LOG/$logf" 2>&1; then echo "   OK"; else
    echo "   FAIL (exit $?) — $LOG/$logf"; FAILED+=("$name"); fi
  tail -4 "$LOG/$logf" | sed 's/^/     /'
}
while pgrep -f "run_stage2.sh" > /dev/null; do sleep 20; done
echo "2단계 완료. 3단계 시작 $(date -u +%H:%M:%S)"

step "커스텀 채점 검증"  s3_custom.log  $PY scripts/custom_scoring_controlled.py
step "그림 최종"         s3_figs.log    $PY scripts/figures_controlled.py
step "보고서 최종"       s3_report.log  $PY scripts/report_controlled.py
step "문서 최종"         s3_docs.log    $PY scripts/make_docs.py \
                             --md sample_run/report/report_controlled.md \
                             --out sample_run/report/docs --stem report_controlled

echo
if [ ${#FAILED[@]} -eq 0 ]; then echo "3단계 전부 통과"; else echo "실패: ${FAILED[*]}"; fi
