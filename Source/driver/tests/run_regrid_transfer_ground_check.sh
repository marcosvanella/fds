#!/bin/bash
# run_regrid_transfer_ground_check.sh <build-dir> [<run-dir>]: boxes created, moved and removed over the burning ground with --nofine 0 (fine k=0 over the stateful vents), CSIRO C064 level-set crude case,
# --regrid-every 2 --tag-phi (60 s). Asserts: boxes are created (box count grows at a regrid) and removed (it shrinks at a later one) while the fire burns; the fine div u - D stays <= 1e-9 in every step,
# also in the step after a regrid that created or re-made boxes (the vent faces of the new layout, the velocity re-apply); the HRR of the step after such a regrid keeps within 35% of the step before
# (a box re-created at the ground over the fire starts without the fine QR, REGRID-TRANSFER-MAP.md gap b, and must not lose the fire). Mass balance of --nofine 0 is not asserted here.
set -u
BLD=$1; RUN=${2:-$BLD/run/regrid_transfer_ground}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
sed -i 's/^&HEAD CHID=[^,]*,/\&HEAD CHID='"'c'"',/' c.fds
sed -i '1a &AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=12 /' c.fds
OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --nofine 0 --outdir . --chid c > run.log 2>&1
python3 - run.log <<'P'
import sys, re
boxes = {}; chg = {}; div = {}; hrr = {}
for l in open(sys.argv[1], errors='replace'):
    m = re.match(r'TWO-LEVEL mid-regrid after step (\d+) .* changed (\d+) finest \d+ boxes L1 (\d+)', l)
    if m: s = int(m.group(1)); chg[s] = int(m.group(2)); boxes[s] = int(m.group(3)); continue
    m = re.match(r'TWO-LEVEL step (\d+) .* L1 (\S+)', l)
    if m: div[int(m.group(1))] = float(m.group(2)); continue
    m = re.match(r'TWO-LEVEL HRR at step (\d+) composite (\S+) kW', l)
    if m: hrr[int(m.group(1))] = float(m.group(2))
bad = 0
if not div: print('no steps'); bad += 1
worst = max(div.values(), default=0.0)
if worst > 1e-9: print('fine div', worst); bad += 1
seq = [boxes[s] for s in sorted(boxes)]
grow = any(b > a for a, b in zip(seq, seq[1:])); shrink = any(b < a for a, b in zip(seq, seq[1:]))
if not grow: print('no box created at a regrid'); bad += 1
if not shrink: print('no box removed at a regrid'); bad += 1
for s in sorted(chg):
    if chg[s] and s - 1 in hrr and s in hrr and hrr[s - 1] > 100.0 and hrr[s] < 0.65 * hrr[s - 1]:
        print('HRR drop after regrid at step', s, hrr[s - 1], '->', hrr[s]); bad += 1
print(('FAIL' if bad else 'PASS') + ' regrid_transfer_ground: %d steps, worst fine div %.2e, %d regrids (%s, %s), %d violation(s)' % (len(div), worst, len(chg), 'created' if grow else 'no create', 'removed' if shrink else 'no remove', bad))
sys.exit(1 if bad else 0)
P
