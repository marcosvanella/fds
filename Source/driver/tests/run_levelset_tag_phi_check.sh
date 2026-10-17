#!/bin/bash
# run_levelset_tag_phi_check.sh <build-dir> [<run-dir>]: 2-level --tag-phi run of the CSIRO C064 level-set crude case (60 s): the fine boxes must sit over the
# columns where PHI crosses the threshold and reach the lowest allowed fine layer (k=2) and the top (k=7) above the no-fine layer of the stateful ground VENTs (D-010, --nofine 1: level-0 layer k=0 stays coarse,
# fine k >= 2); a box list that floats higher (the tag column mis-indexed) or misses the front column fails.
# Open: with --nofine 0 (fine k=0 over the stateful vents) the hidden ground vent walls are not supported yet (crash / gate failure); that run is not asserted here.
set -u
BLD=$1; RUN=${2:-$BLD/run/levelset_tag_phi}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
sed -i 's/^&HEAD CHID=[^,]*,/\&HEAD CHID='"'c'"',/' c.fds
sed -i '1a &AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=12 /' c.fds
OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --nofine 1 --plotfile plt --plot-every 20 --outdir . --chid c > run.log 2>&1
grep -aq "gate PASS" run.log || { echo "FAIL levelset_tag_phi (gate, see $RUN/run.log)"; exit 1; }
python3 - "$RUN" <<'P'
import sys, re, glob, numpy as np
run = sys.argv[1]; bad = 0; n = 0
for d in sorted(glob.glob(run + '/plt/plt_*_s*')):
    boxes = [tuple(map(int, m.groups())) for m in (re.match(r'\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)', l) for l in open(d + '/Level_1/Cell_H')) if m]
    if not boxes: print('no fine box in', d); bad += 1; continue
    if min(b[2] for b in boxes) != 2 or max(b[5] for b in boxes) != 7: print('box k range', min(b[2] for b in boxes), max(b[5] for b in boxes), 'in', d, '(expected 2..7)'); bad += 1
    # front column: level-0 PHI (component 6) crosses 0 along the mid-y row; some fine box must cover its x index
    raw = open(d + '/Level_0/Cell_D_00000', 'rb').read(); h = raw[:raw.index(b'\n')].decode()
    g = re.search(r'\(\((\d+),(\d+),(\d+)\) \((\d+),(\d+),(\d+)\)', h); dims = [int(g.group(k + 3)) - int(g.group(k)) + 1 for k in (1, 2, 3)]
    a = np.frombuffer(raw[len(h) + 1:len(h) + 1 + 8 * int(np.prod(dims)) * 7], '<f8').reshape(7, dims[2], dims[1], dims[0])
    row = a[6, 0, dims[1] // 2, :]
    cross = [i for i in range(len(row) - 1) if row[i] > 0 >= row[i + 1]]
    for i in cross:
        if not any(b[2] == 2 and b[0] // 2 <= i <= b[3] // 2 for b in boxes): print('front column', i, 'not covered in', d); bad += 1
    n += 1
print(('FAIL' if bad or n == 0 else 'PASS') + ' levelset_tag_phi: %d plotfiles checked, %d violation(s)' % (n, bad))
sys.exit(1 if bad or n == 0 else 0)
P
