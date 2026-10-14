#!/bin/bash
# run_open_vent_check.sh <build-dir> [<run-dir>]: two-level regression for an OPEN domain face (Dirichlet pressure) with a burner kept on level 0
# (simple_test + AMR BLOCKING_FACTOR=2, tag-vort auto, regrid every 5). Two faults are covered: the fine H ghost on a Dirichlet face was the even
# mirror (Neumann) while the composite solve folds the level-0 BXF as Dirichlet data (L1 div 0.01 -> 3.5 between regrids, st_tracer), and the
# post-regrid projection was all-Neumann (an open case's net inflow became a removed mean: L0 = L1 div 9e-9 .. 1.7e-3 on changed regrids).
# Cases: open_tracer_two_level (st_tracer: MASS_FLUX tracer burner, 5 steps) and open_burner_two_level (20 steps). Pass: every step max|div u - D| and L1 <= 1e-9.
set -u
BLD=$1; RUN=${2:-$BLD/run/openvent}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
rc=0
for spec in "open_tracer_two_level 5" "open_burner_two_level 20"; do   # tracer: the H-ghost fault (L1 0.012 at step 1); burner: the projection fault (regrid steps 5, 20)
  set -- $spec; C=$1; N=$2
  cp "$HERE/cases/$C.fds" .
  OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" $C.fds --two-level-run --steps $N --regrid-every 5 --plot-every 0 --tag-vort auto --outdir . --chid $C > $C.log 2>&1
  n=$(grep -c "^TWO-LEVEL step " $C.log)
  bad=$(awk '/^TWO-LEVEL step /{for(i=1;i<=NF;i++) if($i=="D|"||$i=="L1"){ if ($(i+1)+0 > 1e-9) {print; break} }}' $C.log | wc -l)
  if [ "$n" -eq "$N" ] && [ "$bad" -eq 0 ]; then echo "PASS open_vent $C: $N steps, max|div u - D| and L1 <= 1e-9"; else echo "FAIL open_vent $C: steps $n of $N, steps with div > 1e-9: $bad (see $RUN/$C.log)"; rc=1; fi
done
exit $rc
