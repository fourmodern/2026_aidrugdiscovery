#!/usr/bin/env bash
# 4단계 — 자체 감사 → 강의 덱 재빌드 → 레포 동기화 준비.
# 감사가 실패하면 덱·동기화를 진행하지 않는다. 틀린 결과를 배포하지 않기 위해서다.
set -u
cd "$(dirname "$0")/.."
PY=/home/hjpark/lecture_drug1/.venv_gpu/bin/python
APY=/share/anaconda3/bin/python
LOG=/tmp/claude-1001/t004

while pgrep -f "run_stage3.sh" > /dev/null; do sleep 20; done
echo "3단계 완료. 4단계 시작 $(date -u +%H:%M:%S)"

echo "── 자체 감사 ──"
if $PY scripts/audit_release.py > "$LOG/s4_audit.log" 2>&1; then
  echo "   감사 통과"
  cat "$LOG/s4_audit.log"
else
  echo "   감사 실패 — 덱·동기화 중단"
  cat "$LOG/s4_audit.log"
  exit 1
fi

echo "── 강의 덱 재빌드 ──"
cd /home/hjpark/lecture_drug1/lecture_materials
if $APY build_0912_harness.py > "$LOG/s4_deck.log" 2>&1; then
  tail -1 "$LOG/s4_deck.log"
else
  echo "   덱 빌드 실패"; tail -5 "$LOG/s4_deck.log"; exit 1
fi

echo "── 레포 동기화 ──"
bash /tmp/claude-1001/t004/sync_repos.sh
echo "4단계 완료 — 커밋은 수동 확인 후"
