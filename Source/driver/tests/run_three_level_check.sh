#!/bin/bash
# D-053 three-level hierarchy check: bind 3 levels (registry + TimeLoop::bind_level), then a short advance with EXACT_SUMS on.
# usage: run_three_level_check.sh <build-dir (absolute)> [work-dir] [steps=10]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/tl3}; STEPS=${3:-10}; rc=0; rm -rf "$W"; mkdir -p "$W"
CASE="$HERE/cases/ns2d_16_l0.fds"
[ -f "$CASE" ] || CASE="$HERE/cases/blob_16_l0.fds"
cp "$CASE" "$W/"; c=$(basename "$CASE" .fds)
# 1) bind-only: hierarchy + bind levels 1 and 2, no advance
(cd "$W" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" "$c.fds" --two-level-run --levels 3 --bind-only --patch 4 11 4 11 --patch2 6 9 6 9 --blocking 4 --maxsize 16 --outdir . --chid bind3 > bind3.log 2>&1) || { echo "RUN FAILED bind-only"; rc=1; }
grep -q "MULTI-LEVEL BIND-ONLY PASS" "$W/bind3.log" && echo "PASS bind-only three levels" || { echo "FAIL bind-only"; grep -E "MULTI-LEVEL|expected finest|Abort|error" "$W/bind3.log" | head -20; rc=1; }
# 2) short advance with exact zone sums (composite over 3 levels)
(cd "$W" && timeout 1800 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" "$c.fds" --two-level-run --levels 3 --steps "$STEPS" --exact-zone-sums --patch 4 11 4 11 --patch2 6 9 6 9 --blocking 4 --maxsize 16 --outdir . --chid run3 --log-every "$STEPS" > run3.log 2>&1) || { echo "RUN FAILED advance"; rc=1; }
if grep -q "MULTI-LEVEL RUN PASS" "$W/run3.log"; then
  echo "PASS three-level advance ($STEPS steps)"
  grep "MULTI-LEVEL RESULT" "$W/run3.log" | head -8
else
  echo "FAIL three-level advance (composite or step loop)"
  grep -E "MULTI-LEVEL|composite pressure|Abort|NotBuilt|error|STOP" "$W/run3.log" | head -30
  rc=1
fi
exit $rc
