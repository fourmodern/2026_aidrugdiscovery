#!/usr/bin/env bash
# 2단계 — 결합양상·접촉 분석·통계검증을 만든 뒤 그림과 보고서를 다시 만들고 문서로 굽는다.
# 1단계(run_controlled_study.sh)가 끝난 뒤 실행된다.
set -u
cd "$(dirname "$0")/.."
PY=/home/hjpark/lecture_drug1/.venv_gpu/bin/python
LOG=/tmp/claude-1001/t004
FIGD=sample_run/report/figures_controlled
declare -a FAILED=()

step () {
  local name="$1" logf="$2"; shift 2
  echo "── $name ──"
  if "$@" > "$LOG/$logf" 2>&1; then echo "   OK"; else
    echo "   FAIL (exit $?) — $LOG/$logf"; FAILED+=("$name"); fi
  tail -3 "$LOG/$logf" | sed 's/^/     /'
}

while pgrep -f "run_controlled_study.sh" > /dev/null; do sleep 20; done
echo "1단계 완료. 2단계 시작 $(date -u +%H:%M:%S)"

step "통계 함수 검증"  s2_stats.log   $PY scripts/test_statistics.py
step "결합양상 렌더"    s2_bind.log    $PY scripts/render_binding_mode.py \
                                          --source controlled --pose top --out "$FIGD"
step "접촉 일치 검정"   s2_conc.log    $PY scripts/contact_concordance.py --source controlled
step "그림 재생성"      s2_figs.log    $PY scripts/figures_controlled.py
step "보고서 재생성"    s2_report.log  $PY scripts/report_controlled.py
step "문서 변환"        s2_docs.log    $PY scripts/make_docs.py \
                                          --md sample_run/report/report_controlled.md \
                                          --out sample_run/report/docs \
                                          --stem report_controlled

echo
if [ ${#FAILED[@]} -eq 0 ]; then echo "2단계 전부 통과"; else echo "실패: ${FAILED[*]}"; fi
