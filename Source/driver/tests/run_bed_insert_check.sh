#!/bin/bash
# run_bed_insert_check.sh <build-dir> [<run-dir>]: cases/bed_insert_c064.fds (C064 cut-down c10, &INIT grass bed XB 0,25,-15,15,0,0.21, T_END 0.3 s) as a two-level
# run with --bed on the same XB. The &INIT bed must be inserted on the bed level L_ins = 1 (0.5 m: 50 x 60 x 1 cells, 1 particle per cell = 3000) and the level-0
# particles under level 1 removed (leaf rule); before the fix --bed skipped the insert and the 750 level-0 particles stayed on covered cells (mass gate 1.40, rc 6).
# Pass: exit 0 (TWO-LEVEL gate PASS) and "particles at t=0: 3000".
BLD=$1; RUN=${2:-$BLD/run/bed_insert}; HERE=$(cd "$(dirname "$0")" && pwd)
rm -rf "$RUN"; mkdir -p "$RUN"; cd "$RUN" || exit 1
cp "$HERE/cases/bed_insert_c064.fds" c.fds; cp "$HERE/../../../Validation/CSIRO_Grassland_Fires/FDS_Input_Files/vegetation_model.txt" . || exit 1
OMP_NUM_THREADS=1 mpirun --oversubscribe -np 4 "$BLD/fds_amr" c.fds --two-level-run --bed 0 25 -15 15 0 0.21 > run.log 2>&1; rc=$?
n=$(grep -a "particles at t=0:" run.log | head -1 | sed 's/.*t=0: \([0-9]*\),.*/\1/')
g=$(grep -a "RESULT gate" run.log | grep -o "gate [A-Z]*: .*mass balance| / limit [0-9.e+-]*" | sed 's/: .*limit / mass\/limit /')
echo "rc=$rc particles at t=0: ${n:-none}; $g"
[ $rc -eq 0 ] && [ "${n:-0}" = 3000 ] && echo "PASS bed_insert" && exit 0
echo "FAIL bed_insert"; exit 1
