#!/bin/bash
# Wall-state seam (ADR-001 W1/W2, notes/wall-seam-design.md), host-side test:
#   FDSTL_WSEAM=1  before the density update and before WALL_BC of every stage UVW_SAVE/U_GHOST/V_GHOST/W_GHOST are staged over the box's WLIST_EXT with a checksum assertion, and
#                  after the stage the host arrays are asserted unchanged: the run must finish and be bitwise equal to the plain run (final fields and step log).
#   FDSTL_WSEAM=2  negative control: a host value is perturbed after the upload; the run must STOP with the checksum-assertion message.
# usage: run_wall_seam_check.sh <build-dir> [steps] [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; STEPS=${2:-10}; WORK=${3:-$BLD/run/wallseam}; rm -rf "$WORK"; mkdir -p "$WORK"; rc=0
run() {  # name case np mode
  local D="$WORK/$1"; mkdir -p "$D"; cd "$D"; cp "$HERE/cases/$2.fds" .
  FDSTL_WSEAM=$4 timeout 900 mpirun --bind-to none --oversubscribe -np $3 "$BLD/fds_amr" $2.fds --run --outdir . --chid $2 --steps $STEPS --exact-zone-sums > run.out 2> run.err
  return $?
}
same() { local ok=1; for f in "$WORK/$1/$3"_final_*.bin; do cmp -s "$f" "$WORK/$2/$(basename $f)" || ok=0; done; cmp -s "$WORK/$1/$3_driver_steps.csv" "$WORK/$2/$3_driver_steps.csv" || ok=0; echo $ok; }
for T in "dec1 1" "dec4_np4 4"; do
  set -- $T; C=$1; NP=$2
  run ${C}_plain $C $NP 0 || { echo "FAIL $C np=$NP plain run"; rc=1; continue; }
  if run ${C}_ws1 $C $NP 1 && [ "$(same ${C}_plain ${C}_ws1 $C)" = 1 ]; then
    echo "PASS $C np=$NP FDSTL_WSEAM=1: stage + checksum assertion, final fields and step log bitwise equal to the plain run ($(grep -a 'wall seam check' $WORK/${C}_ws1/run.out | head -1 | sed 's/.*(W1/W1/'))"
  else echo "FAIL $C np=$NP FDSTL_WSEAM=1"; rc=1; fi
  if run ${C}_ws2 $C $NP 2; then echo "FAIL $C np=$NP negative control: the run did not stop"; rc=1
  elif grep -aq "checksum assertion of ADR-001 W2" $WORK/${C}_ws2/run.out $WORK/${C}_ws2/run.err; then echo "PASS $C np=$NP negative control (host value perturbed after the upload): the run stops with the checksum-assertion message"
  else echo "FAIL $C np=$NP negative control: stopped without the expected message"; rc=1; fi
done
exit $rc
