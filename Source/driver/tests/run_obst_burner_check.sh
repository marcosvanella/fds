#!/bin/bash
# run_obst_burner_check.sh <build-dir> [<run-dir>]: regression for an L0-only OBST with a MASS_FLUX burner on its top face, two levels
# (fine patch away from the OBST; work-species simple_test_s1). Level 0 keeps FDS's own OBST treatment in the composite solve (no L0 mask,
# no zeroed L0 solid-face FV): the burner inflow stays the wall normal velocity. Before the fix L0 div was 0.075 at step 1 and 4679 at step 6.
# Pass: every step max|div u - D| <= 1e-10 and no L0 clipped cells.
set -u
BLD=$1; RUN=${2:-$BLD/run/obstburner}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/obst_burner_two_level.fds" .
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" obst_burner_two_level.fds --two-level-run --steps 10 --patch 8 23 0 7 --outdir . --chid ob > run.log 2>&1
n=$(grep -c "^TWO-LEVEL step " run.log)
bad=$(awk '/^TWO-LEVEL step /{for(i=1;i<=NF;i++) if($i=="D|"){ if ($(i+1)+0 > 1e-10) print }}' run.log | wc -l)
clip=$(grep -cE "^CLIP_CELLS N=[0-9]+ L=0 count=[1-9]" run.log)
if [ "$n" -eq 10 ] && [ "$bad" -eq 0 ] && [ "$clip" -eq 0 ]; then echo "PASS obst_burner: 10 steps, max|div u - D| <= 1e-10, no L0 clipping"; exit 0; fi
echo "FAIL obst_burner: steps $n, steps with div > 1e-10: $bad, L0 clip lines: $clip (see $RUN/run.log)"; exit 1
