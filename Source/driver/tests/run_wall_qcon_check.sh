#!/bin/bash
# run_wall_qcon_check.sh <fds_amr binary> [<run-dir>]: gas-side convective wall flux of the hidden fine walls (FDS_FINE_WALLS_QCON, read by the fine DIVERGENCE_PART_1).
# cases/t0_state.fds (hot block on the floor). FDSTL_T0RAD / FDSTL_QCON print per level, over the SOLID walls whose gas cell is leaf: solid side sum(Q_CON_F*AREA),
# gas side sum(AREA_ADJUST*Q_CON_F*RDN*V) (what DIVERGENCE_PART_1 adds), and sum|Q_CON_F*AREA| of the non-leaf walls.
# 1. Leaf-only energy at t=0: level 1 over the whole domain vs 1 level at the level-1 dx: the level-1 leaf sum equals the 1-level sum (1e-12), the level-0 leaf sum is 0.
# 2. Closure, level 1 over half the domain, steps 1-6: on every level gas side = solid side (1e-12 relative), and the non-leaf hidden walls carry no flux (a face counts once).
# 3. Thermally thick floor (cases/qcon_thick.fds, level 1 over half the domain, 8 steps, WALL_INCREMENT 2): the covered level-0 solids take the area-weighted leaf fine
#    Q_CON_F (L0QCON: fine energy committed vs energy the solids took, cumulative); at every solid update they agree to 1e-12 relative, and the solids took heat.
set -u
BIN=$1; RUN=${2:-$PWD/wall_qcon_run}; SRC=$(cd "$(dirname "$0")" && pwd)/cases/t0_state.fds
mkdir -p "$RUN" && cd "$RUN" || exit 1
sed 's/IJK=8,8,4/IJK=16,16,8/' "$SRC" > r1.fds; sed 's/T_END=0.1/T_END=2./' "$SRC" > r2.fds
run() { local L=$1 log=$2; shift 2; env FDSTL_T0RAD=1 FDSTL_QCON=1 OMP_NUM_THREADS=1 timeout 600 "$BIN" r$L.fds --two-level-run --levels $L --maxsize 16 --nofine -1 --outdir . --chid r$L "$@" > $log 2>&1 \
    || { echo "FAIL wall_qcon: run rc (see $RUN/$log)"; exit 1; }; }
run 1 r1.log --patch 0 7 0 3 --steps 1
run 2 r2.log --patch 0 7 0 3 --steps 1
run 2 r2p.log --patch 0 3 0 3 --steps 6
cp "$(dirname "$SRC")/qcon_thick.fds" r3.fds
env FDSTL_QCON=1 OMP_NUM_THREADS=1 timeout 600 "$BIN" r3.fds --two-level-run --levels 2 --maxsize 16 --nofine -1 --outdir . --chid r3 --patch 0 3 0 3 --steps 8 > r3.log 2>&1 \
    || { echo "FAIL wall_qcon: thick run rc (see $RUN/r3.log)"; exit 1; }
rm -f ./*.s3d ./*.smv ./*.sf ./*.bf ./*.prt5 ./*.csv
grep -a "^T0QCON" r1.log r2.log; grep -a "^QCON" r2p.log; grep -a "^L0QCON" r3.log
python3 - <<'P'
import sys
def rd(f, tag):
    return [l.split() for l in open(f, errors='replace') if l.startswith(tag)]
fail = 0
one = {r[1]: list(map(float, (r[3], r[5], r[7]))) for r in rd('r1.log', 'T0QCON ')}
two = {r[1]: list(map(float, (r[3], r[5], r[7]))) for r in rd('r2.log', 'T0QCON ')}
e = abs(two['L1'][0] - one['L0'][0]) / abs(one['L0'][0])
print('t=0 leaf Q_CON: 1L %.17g W, 2L level 1 %.17g W, rel %.2e; 2L level-0 leaf %.3g W' % (one['L0'][0], two['L1'][0], e, two['L0'][0]))
if not (one['L0'][0] != 0 and e <= 1e-12 and two['L0'][0] == 0): print('  leaf-only check failed'); fail = 1
rows = rd('r2p.log', 'QCON step')
if len(rows) < 12: print('  closure: missing QCON lines'); fail = 1
w = 0.0
for r in rows:
    s, g, n = float(r[5]), float(r[7]), float(r[9])
    w = max(w, abs(g - s) / max(abs(s), 1e-300))
    if r[3] != 'L0' and n != 0: print('  non-leaf hidden walls carry flux:', ' '.join(r)); fail = 1
print('closure worst |gas - solid|/|solid| %.2e over %d level-steps' % (w, len(rows)))
if w > 1e-12: fail = 1
th = [(int(r[2]), float(r[4]), float(r[6])) for r in rd('r3.log', 'L0QCON ')]
upd = [(n, f, s) for (n, f, s) in th if n % 2 == 0]
wt = max([abs(f - s) / max(abs(f), 1e-300) for (n, f, s) in upd] or [1.0])
print('thick floor: %d solid updates, fine %.17g J vs solid %.17g J at the last, worst rel %.2e' % (len(upd), upd[-1][1] if upd else 0, upd[-1][2] if upd else 0, wt))
if len(upd) < 4 or wt > 1e-12 or upd[-1][2] <= 0: print('  thick closure failed'); fail = 1
print('PASS wall_qcon' if fail == 0 else 'FAIL wall_qcon'); sys.exit(fail)
P
