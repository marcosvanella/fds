#!/bin/bash
# run_w12_smoke.sh <fds_amr> <c10 cut-down C064 input (4 meshes at 1 m, &AMR line, vegetation_model.txt beside it)> [work-dir] [ranks]
# P0 W12 smoke (G-C): 2 levels, bed particles, regrid every 5 steps, to T_END of the input. PASS when the run ends at T_END with the TWO-LEVEL gate PASS
# (fine div, mass and species balance, species check 7) and every W8 regrid check keeps the particle count and record hash (FDSTL_W8_CHECK SAME).
# W6, W9, W10, W11 have no gated tests yet (PROGRESS.md proposals); restart-and-continue is not implemented in the driver.
set -u
BIN=$1; IN=$2; W=${3:-$PWD/w12_smoke}; NP=${4:-1}
rm -rf "$W"; mkdir -p "$W"; cp "$IN" "$W/c10.fds"; cp "$(dirname "$IN")/vegetation_model.txt" "$W/"
cd "$W" && FDSTL_W8_CHECK=1 OMP_NUM_THREADS=1 timeout 3600 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" c10.fds --two-level-run --bed -10 30 -20 20 0 20 \
    --regrid-every 5 --tag-temp 330 > log 2>&1; rc=$?
rm -f *.s3d* *.sf *.bf *.prt5* *.restart
why=""
[ $rc = 0 ] || why="rc $rc"
grep -aq "RESULT gate PASS" log || why="$why; gate: $(grep -a 'RESULT gate' log | cut -c1-200)"
grep -a "RESULT species check 7" log | grep -q PASS || why="$why; species check 7"
nw=$(grep -ac "W8 CHECK" log); nb=$(grep -a "W8 CHECK" log | grep -vc SAME)
[ "$nw" -gt 0 ] && [ "$nb" = 0 ] || why="$why; W8 checks $nw, not SAME $nb"
t=$(grep -a "RESULT composite HRR final" log | grep -o "t [0-9.]* (T_END [0-9.]*)")
if [ -z "$why" ]; then echo "PASS w12_smoke np$NP ($t, W8 checks $nw SAME)"; else echo "FAIL w12_smoke np$NP:$why ($t)"; exit 1; fi
