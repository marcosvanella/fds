#!/bin/bash
# run_obst_offgrid_check.sh <build-dir> [levels=3] [<run-dir>]: repro of the N-level fine OBST wall gap. The OBST and its burner (XB 0.03,0.43,1.02,1.38,0,0.18) are off the level-0 grid,
# so the finest level snaps them differently from level 0 (every level follows FDS's NINT at its own dx), and the finest level has solid faces with no level-0 parent wall.
# Before the per-level fine walls (5d56a9010f, c1d59c514a): no crash, but max|div u - D| on the finest level was 0.065 (3 levels) / 0.41 (2 levels) from step 1, levels below 1e-15.
# Passes when every step stays <= 1e-9 (since 5d56a9010f: 6.6e-15 / 2.4e-15; notes/fine-obst-walls.md R12, M1).
set -u
BLD=$1; LEV=${2:-3}; RUN=${3:-$BLD/run/obstoffgrid}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/obst_offgrid_levels.fds" .
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" obst_offgrid_levels.fds --two-level-run --levels $LEV --steps 6 --tag-obst --nofine -1 --tag-temp 400 --regrid-every 3 --blocking 4 --maxsize 8 --outdir . --chid og --log-every 1 > run.log 2>&1
n=$(grep -cE "^(TWO|MULTI)-LEVEL step " run.log)
worst=$(awk '/^(TWO|MULTI)-LEVEL step /{for(i=1;i<=NF;i++) if($i=="D|" && $(i+1)+0 > w) w=$(i+1)+0} END{printf "%.3g", w}' run.log)
if [ "$n" -eq 6 ] && awk -v w="$worst" 'BEGIN{exit !(w <= 1e-9)}'; then echo "PASS obst_offgrid ($LEV levels): 6 steps, worst max|div u - D| $worst"; exit 0; fi
echo "FAIL obst_offgrid ($LEV levels): $n of 6 steps, worst max|div u - D| $worst (see $RUN/run.log)"; exit 1
