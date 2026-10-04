#!/bin/bash
# Patch 0011 check (needs main.f90 patch 0011 applied and a USE_AMREX=ON build): fds_amr --run with FDSTL_MESH_DUMPS=1 writes the per-mesh FDS dumps
# (slice files, boundary files) through FDS_SETUP(MODE=3). Compared with the USE_AMREX=OFF baseline ($BASELINE) of shunn3_32 (1 rank) and shunn3_4mesh_32 (4 ranks):
# same list of .sf/.sf.bnd files, same sizes, and the data equal to the baseline within single-precision round-off (bitwise count reported).
# usage: run_mesh_dumps_check.sh <build-dir (absolute)> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/mdchk}; rc=0; mkdir -p "$W"
for spec in "shunn3_32 1" "shunn3_4mesh_32 4"; do
  set -- $spec; c=$1; np=$2; d=$W/$c; rm -rf "$d"; mkdir -p "$d"; cp "$BASELINE/$c/$c.fds" "$d/"
  (cd "$d" && FDSTL_MESH_DUMPS=1 timeout 1500 mpirun --bind-to none --oversubscribe -np $np "$B/fds_amr" "$c.fds" --run --outdir . --chid "$c" --quiet > stdout.txt 2> stderr.txt) || { echo "RUN FAILED $c"; rc=1; continue; }
  python3 - "$d" "$BASELINE/$c" "$c" <<'PY' || rc=1
import sys, os, glob, numpy as np
d, ref, c = sys.argv[1:4]
names = sorted(os.path.basename(f) for f in glob.glob(ref + '/*.sf') + glob.glob(ref + '/*.sf.bnd'))
miss = [n for n in names if not os.path.exists(d + '/' + n)]
size = [n for n in names if n not in miss and os.path.getsize(d + '/' + n) != os.path.getsize(ref + '/' + n)]
same = 0; worst = 0.0
for n in names:
    if n in miss or n in size: continue
    if open(d + '/' + n, 'rb').read() == open(ref + '/' + n, 'rb').read(): same += 1; continue
    if n.endswith('.bnd'): worst = max(worst, 1.0); continue
    a = np.fromfile(d + '/' + n, dtype=np.float32); b = np.fromfile(ref + '/' + n, dtype=np.float32)
    s = max(np.abs(b).max(), 1e-30); worst = max(worst, float(np.abs(a - b).max() / s))
ok = len(names) > 0 and not miss and not size and worst < 1e-5
print(('PASS ' if ok else 'FAIL ') + c + ': %d slice/boundary files in the baseline, %d missing, %d with other size, %d bitwise equal, worst relative difference %.2e' % (len(names), len(miss), len(size), same, worst))
sys.exit(0 if ok else 1)
PY
done
exit $rc
