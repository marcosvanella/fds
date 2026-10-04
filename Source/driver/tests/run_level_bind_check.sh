#!/bin/bash
# S12: TimeLoop::bind_level. (1) `--level-bind-check`: a ratio-1 copy of level 0 is bound as level 1 and the first stage of a step is run on both levels position by position;
# every field must agree bitwise (documented exceptions in tests/level_bind_check.cpp: FVY 1e-16, D 4e-15, seam-side FVX/FVZ 1e-9, and the known DS gap of the predictor DIVERGENCE_PART_1 for non-uniform species), 1 rank and 4 ranks (box-box faces of the fine level exchanged by the per-level BcStep).
# (2) optional single-level regression: with a reference executable (built from the tree before bind_level) the final fields and the step log of a normal run must be bitwise equal.
# usage: run_level_bind_check.sh <build-dir> [reference fds_amr] [steps] [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; REF=${2:-}; STEPS=${3:-10}; WORK=${4:-$BLD/run/levelbind}; rm -rf "$WORK"; mkdir -p "$WORK"; rc=0
export OMP_NUM_THREADS=1
for T in "ns2d_16_l0 1" "dec1 1" "dec4_np4 4"; do
  set -- $T; C=$1; NP=$2; D="$WORK/$C"; mkdir -p "$D"; cd "$D"; cp "$HERE/cases/$C.fds" .
  timeout 600 mpirun --bind-to none --oversubscribe -np $NP "$BLD/fds_amr" $C.fds --level-bind-check > lb.out 2> lb.er
  if grep -q "LEVEL-BIND-CHECK PASS" lb.out; then echo "PASS $C np=$NP level-bind-check ($(grep -c '  ok ' lb.out) fields equal to level 0 [bitwise except the documented rounding-level ones], $(grep -c '  GAP ' lb.out) known gap)"
  else echo "FAIL $C np=$NP level-bind-check"; grep -a "DIFF\|Abort\|ERROR" lb.out lb.er | head -5; rc=1; fi
done
if [ -n "$REF" ] && [ -x "$REF" ]; then
  for T in "ns2d_16_l0 1" "dec1 1" "dec4_np4 4"; do
    set -- $T; C=$1; NP=$2
    for V in new ref; do
      X=$BLD/fds_amr; [ $V = ref ] && X=$REF
      D="$WORK/${C}_$V"; mkdir -p "$D"; cd "$D"; cp "$HERE/cases/$C.fds" .
      timeout 600 mpirun --bind-to none --oversubscribe -np $NP "$X" $C.fds --run --outdir . --chid $C --steps $STEPS --exact-zone-sums > run.out 2> run.er
    done
    ok=1; n=0
    for f in "$WORK/${C}_new/$C"_final_*.bin; do n=$((n+1)); cmp -s "$f" "$WORK/${C}_ref/$(basename $f)" || ok=0; done
    cmp -s "$WORK/${C}_new/${C}_driver_steps.csv" "$WORK/${C}_ref/${C}_driver_steps.csv" || ok=0
    if [ $ok = 1 ] && [ $n -gt 0 ]; then echo "PASS $C np=$NP single-level run bitwise equal to the reference executable ($n final-field files + step log, $STEPS steps)"
    else echo "FAIL $C np=$NP single-level run differs from the reference executable (files: $n)"; rc=1; fi
  done
else echo "NOTE no reference executable given: single-level regression not run"; fi
exit $rc
