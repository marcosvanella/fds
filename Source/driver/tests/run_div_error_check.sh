#!/bin/bash
# D-058 diagnostic FluxStages::max_divergence_error (max |div u - D| over the valid cells of level 0), printed by FDSTL_DIVERR=1 at the end of a run.
# After a full step with the pressure solve D equals div u to the solver tolerance: the value must be below 1e-8 x max(1, max|D|), and max|D| must be non-zero (the cases carry
# species and so a real divergence). Results are not changed by the flag (final fields bitwise equal to a plain run).
# usage: run_div_error_check.sh <build-dir> [steps] [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; STEPS=${2:-10}; WORK=${3:-$BLD/run/diverr}; rm -rf "$WORK"; mkdir -p "$WORK"; rc=0
for T in "dec1 1" "dec4_np4 4"; do
  set -- $T; C=$1; NP=$2
  for v in plain diverr; do
    D="$WORK/${C}_$v"; mkdir -p "$D"; cd "$D"; cp "$HERE/cases/$C.fds" .
    if [ $v = diverr ]; then export FDSTL_DIVERR=1; else unset FDSTL_DIVERR; fi
    timeout 600 mpirun --bind-to none --oversubscribe -np $NP "$BLD/fds_amr" $C.fds --run --outdir . --chid $C --steps $STEPS --exact-zone-sums > run.out 2> run.err || { echo "FAIL $C $v: run failed"; rc=1; }
  done
  unset FDSTL_DIVERR
  line=$(grep -a "^DIVERR" "$WORK/${C}_diverr/run.out" | head -1)
  e=$(echo "$line" | sed -n 's/.*cells = \([^ ]*\) (1\/s).*/\1/p'); d=$(echo "$line" | sed -n 's/.*max |D| = \(.*\)$/\1/p')
  ok=$(python3 -c "e=float('${e:-1}'); d=float('${d:-0}'); print(1 if (d>0 and e < 1e-8*max(1.0,d)) else 0)")
  same=1; for f in "$WORK/${C}_plain/${C}"_final_*.bin; do cmp -s "$f" "$WORK/${C}_diverr/$(basename $f)" || same=0; done
  if [ "$ok" = 1 ] && [ $same = 1 ]; then echo "PASS $C np=$NP: max|div u - D| = $e (max|D| = $d), results bitwise equal to the plain run"
  else echo "FAIL $C np=$NP: '$line' same=$same"; rc=1; fi
done
exit $rc
