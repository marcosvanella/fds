#!/bin/bash
# D-056 (option B) guard: a kernel wrapper called with a fine-level mesh number (above NMESHES) must stop the run with a clear message and a non-zero exit code
# (nothing may read MESHES(NM) for such a number). Uses dec1.fds and the mode `--fine-guard-test` of fds_amr.
# With a build that has patch 0007 and -DFDS_AMR_FINE_B_DRAFT=ON (scratch tree, see patches/0007-mesh-fine-level-boxes.md) the draft checks run too (arg 3 = 1):
#  `--fine-b-selftest` (pointer routing of POINT_TO_BOX to FINE_LEVEL, level 0 unchanged) and `--fine-b-abort 1|2|3` (POINT_TO_MESH on a fine number aborts with the message,
#  POINT_TO_BOX on an unknown number aborts, POINT_TO_BOX on a fine box works).
# usage: run_fine_guard_check.sh <build-dir> [work-dir] [draft=0|1]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; WORK=${2:-$BLD/run/fineguard}; DRAFT=${3:-0}; rm -rf "$WORK"; mkdir -p "$WORK"; cp "$HERE/cases/dec1.fds" "$WORK/"; rc=0
cd "$WORK"
run() { timeout 300 mpirun --bind-to none --oversubscribe -np 1 "$BLD/fds_amr" dec1.fds "$@" > "$OUT.out" 2> "$OUT.err"; }
OUT=guard; run --fine-guard-test; r=$?
if [ $r -ne 0 ] && grep -q "mesh number 1000000 is not a level-0 FDS mesh" guard.err && ! grep -q "FAIL: returned" guard.out; then echo "PASS wrapper guard (exit $r, message printed)"; else echo "FAIL wrapper guard (exit $r)"; rc=1; fi
if [ "$DRAFT" = 1 ]; then
  OUT=bsel; run --fine-b-selftest; r=$?; if [ $r -eq 0 ] && grep -q "FINE-B SELFTEST PASS" bsel.out; then echo "PASS fine-b selftest"; else echo "FAIL fine-b selftest (exit $r)"; cat bsel.err | head; rc=1; fi
  for k in 1 2; do OUT=babort$k; run --fine-b-abort $k; r=$?
    if [ $r -ne 0 ] && grep -q "is not a level-0 FDS mesh" babort$k.err && ! grep -q "returned" babort$k.out; then echo "PASS abort case $k (exit $r)"; else echo "FAIL abort case $k (exit $r)"; rc=1; fi; done
  OUT=babort3; run --fine-b-abort 3; r=$?; if [ $r -eq 0 ] && grep -q "returned" babort3.out; then echo "PASS fine box reachable through POINT_TO_BOX"; else echo "FAIL case 3 (exit $r)"; rc=1; fi
fi
exit $rc
