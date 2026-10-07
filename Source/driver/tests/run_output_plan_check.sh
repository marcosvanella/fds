#!/bin/bash
# OutputPlan + D6 check (notes/output-plan.md §4–5): driver_unit_tests test_output_plan.
# TwoLevelRun FDSTL_OUTPLAN diagnostic hook deferred (Role 3 owns TwoLevelRun).
# usage: run_output_plan_check.sh <build-dir>
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1
export FDSTL_SRC="$(cd "$HERE/.." && pwd)"
rc=0
for np in 1 4; do
  out=$(mpirun --bind-to none --oversubscribe -np $np "$B/driver_unit_tests" 2>&1) || true
  echo "$out" | grep -E "OutputPlan|RESULT D6|ALL PASS|SOME FAILED" | sed "s/^/[outplan np=$np] /"
  echo "$out" | grep -q "RESULT D6 exact composite mass" || { echo "[outplan np=$np] FAIL no D6 line"; rc=1; continue; }
  echo "$out" | grep -q "OutputPlan groups" || { echo "[outplan np=$np] FAIL no OutputPlan groups"; rc=1; continue; }
  echo "$out" | grep -q "^PASS OutputPlan layout/fill + D6" || { echo "[outplan np=$np] FAIL OutputPlan test"; rc=1; continue; }
  python3 -c "
import re,sys
t=sys.stdin.read()
m=re.search(r'D6 exact composite mass ([0-9.eE+-]+) \(FP composite ([0-9.eE+-]+), covered-coarse extra ([0-9.eE+-]+), naive with double-count ([0-9.eE+-]+)\)',t)
assert m, 'no D6 parse'
ex,fp,exx,naive=map(float,m.groups())
print(f'  D6 exact {ex} FP {fp} covered_extra {exx} naive {naive}')
assert abs(ex-fp)<=1e-9*max(1.0,abs(fp)), 'exact vs FP'
assert exx>0 and ex<naive, 'exact must beat double-count'
print('PASS D6 uncovered')
" <<<"$out" || rc=1
done
[ $rc = 0 ] && echo "PASS output_plan check" || echo "FAIL output_plan check"
exit $rc
