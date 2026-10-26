#!/bin/bash
# run_drag_clip_check.sh <fds_amr> <c10 input as run_w12_smoke.sh> [work-dir] [ranks] [t_end=0.5]
# Role 3 design test 5 (PV-4 drag stencil clip): 2 levels, bed particles, FDSTL_DRAG_CHECK=1. PASS when every particle's deposited momentum equals its drag to 1e-15
# relative and none is left without a kept face; reports the weights dropped on fine boxes (seam faces of the box below, ghost faces) and on level 0 (coarse/fine faces).
set -u
BIN=$1; IN=$2; W=${3:-$PWD/drag_clip}; NP=${4:-1}; TE=${5:-0.5}
rm -rf "$W"; mkdir -p "$W"; sed "s/T_END=[0-9.]*/T_END=$TE/" "$IN" > "$W/c10.fds"; cp "$(dirname "$IN")/vegetation_model.txt" "$W/"
cd "$W" && FDSTL_DRAG_CHECK=1 OMP_NUM_THREADS=1 timeout 1800 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" c10.fds --two-level-run --bed -10 30 -20 20 0 20 > log 2>&1; rc=$?
r=$(grep -a "RESULT drag check 5" log | sed 's/.*drag check 5: //')
if [ $rc = 0 ] && [[ "$r" == *PASS ]]; then rm -rf "$W"; echo "PASS drag_clip np$NP: $r"; else echo "FAIL drag_clip np$NP (rc $rc): $r"; exit 1; fi
