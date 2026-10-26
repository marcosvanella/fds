#!/bin/bash
# run_checkpoint_m1_check.sh <fds_amr> <c10 input as run_w12_smoke.sh> [work-dir] [ranks] [step=30]
# Restart M1 smoke: c10 2 levels, bed particles, regrid every 5, --checkpoint-at STEP; the same run without the flag. PASS when the checkpoint has its Header, a
# field file per level and registry field, part_r<rank> and walls_r<rank> per rank and the level-0 stock restart files, and writing it leaves the run unchanged
# (two_level_multi_level.csv identical).
set -u
BIN=$1; IN=$2; W=${3:-$PWD/chk_m1}; NP=${4:-1}; S=${5:-30}
run() { rm -rf "$W/$1"; mkdir -p "$W/$1"; cp "$IN" "$W/$1/c10.fds"; cp "$(dirname "$IN")/vegetation_model.txt" "$W/$1/"
  (cd "$W/$1" && OMP_NUM_THREADS=1 timeout 3600 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" c10.fds --two-level-run --bed -10 30 -20 20 0 20 \
    --regrid-every 5 --tag-temp 330 $2 > log 2>&1; echo $? > rc; rm -f *.s3d* *.sf *.bf *.prt5*); }
run on "--checkpoint-at $S"; run off ""
C=$(ls -d "$W"/on/chk_*_s$S 2>/dev/null | head -1); why=""
[ "$(cat "$W/on/rc")$(cat "$W/off/rc")" = 00 ] || why="rc"
[ -n "$C" ] && [ -s "$C/Header" ] || why="$why; no Header"
for r in $(seq 0 $((NP-1))); do [ -s "$C/part_r$r" ] && [ -s "$C/walls_r$r" ] || why="$why; rank $r files"; done
nf=$(ls "$C" | grep -c '_H$'); nr=$(ls "$W/on" | grep -c '\.restart$')
[ "$nf" -gt 0 ] && [ "$nr" -gt 0 ] || why="$why; fields $nf restart $nr"
cmp -s "$W/on/two_level_multi_level.csv" "$W/off/two_level_multi_level.csv" || why="$why; run changed by the checkpoint"
sz=$(du -sk "$C" | cut -f1); szr=$(du -ck "$W"/on/*.restart | tail -1 | cut -f1)
if [ -z "$why" ]; then echo "PASS checkpoint_m1 np$NP: step $S, $nf field files, $nr restart files, checkpoint ${sz} kB + restart ${szr} kB, run unchanged"; rm -rf "$W"
else echo "FAIL checkpoint_m1 np$NP:$why"; exit 1; fi
