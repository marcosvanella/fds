#!/bin/bash
# S9.3 (A)(i): the FDSTL_STAGES dump of CELL%SOLID (scr_static_SOLID_b<box>.bin). dec2_obst.fds = dec2 (4 boxes of 16x1x16) with one obstruction -0.25..0.25 x -0.25..0.25
# straddling all four boxes: every box holds 4x4 solid cells at its inner corner plus (FDS index 0 or IBP1) a ghost layer from the neighbour box -> 25 values of 1 per box.
# usage: run_solid_dump_check.sh <build-dir> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; WORK=${2:-$BLD/run/solid}; rm -rf "$WORK"; mkdir -p "$WORK"; cp "$HERE/cases/dec2_obst.fds" "$WORK/"; rc=0
(cd "$WORK" && FDSTL_STAGES=1 timeout 300 mpirun --bind-to none --oversubscribe -np 1 "$BLD/fds_amr" dec2_obst.fds --run --outdir . --chid dec2_obst --quiet > stdout.txt 2> stderr.txt) || { echo "FAIL run"; exit 1; }
python3 - "$WORK" <<'P' || rc=1
import numpy as np, sys
w = sys.argv[1]; bad = 0
exp = {0: ((13, 13), (17, 17)), 1: ((0, 13), (4, 17)), 2: ((13, 0), (17, 4)), 3: ((0, 0), (4, 4))}
for b in range(4):
    r = open('%s/scr_static_SOLID_b%d.bin' % (w, b), 'rb').read()
    h = np.frombuffer(r[:64], dtype=np.int32); d = np.frombuffer(r[64:], dtype=np.float64)
    ok = list(h[:9]) == [3, 0, 0, 0, 1, 17, 2, 17, 1] and d.size == 18 * 3 * 18
    d = d.reshape(18, 3, 18, order='F')[:, 1, :]; nz = np.argwhere(d > 0)
    ok = ok and set(np.unique(d)) <= {0.0, 1.0} and int(d.sum()) == 25 and tuple(nz.min(0)) == exp[b][0] and tuple(nz.max(0)) == exp[b][1]
    print('box', b, 'header/shape/values', 'ok' if ok else 'FAIL'); bad += not ok
print('PASS' if not bad else 'FAIL'); sys.exit(bad != 0)
P
exit $rc
