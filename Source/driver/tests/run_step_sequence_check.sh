#!/bin/bash
# Step sequence (T, DT per step) of the driver against plain FDS on a multi-mesh input: shunn3_4mesh_32 (2x2 meshes), driver at 1 and 4 ranks, first 40 steps, against the
# USE_AMREX=OFF baseline $BASELINE/shunn3_4mesh_32/*_steps.csv (FDS writes the diagnostic steps only, so the steps are matched by number).
# Checks: (1) driver np 1 and np 4 give the same T and DT sequence; (2) DT of every compared step is within 2e-3 of the baseline (3 printed digits);
# (3) T within 2e-3 of the baseline.
# Known, not a driver defect (notes/multi-mesh-step-sequence.md): plain FDS on the 2x2 arrangement gives a step-1 DT 1.0e-3 larger than FDS on the same domain as 1 mesh or as
# 2 meshes (0.6195e-3 vs 0.6188e-3); the driver follows the 1-mesh value (0.6189e-3), so T is offset by -8.5e-4 (relative) after step 1.
# usage: run_step_sequence_check.sh <build-dir (absolute)> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/seqchk}; c=shunn3_4mesh_32; rc=0; rm -rf "$W"; mkdir -p "$W/np1" "$W/np4"
for np in 1 4; do
  cp "$BASELINE/$c/$c.fds" "$W/np$np/"
  (cd "$W/np$np" && timeout 1500 mpirun --bind-to none --oversubscribe -np $np "$B/fds_amr" $c.fds --run --steps 40 --outdir . --chid $c --quiet > stdout.txt 2> stderr.txt) || { echo "RUN FAILED np=$np"; rc=1; }
done
python3 - "$W" "$BASELINE/$c/${c}_steps.csv" $c <<'PY' || rc=1
import sys, numpy as np
W, ref, c = sys.argv[1:4]
rd = lambda f: np.genfromtxt(f, delimiter=',', skip_header=2, usecols=(0, 2, 3))
a1, a4, b = rd(f'{W}/np1/{c}_steps.csv'), rd(f'{W}/np4/{c}_steps.csv'), rd(ref)
ok = True
def chk(name, cond, msg=''):
    global ok; print(('PASS ' if cond else 'FAIL ') + name + (' ' + msg if msg else '')); ok = ok and cond
chk('np1 and np4 give the same step sequence', a1.shape == a4.shape and np.allclose(a1, a4, rtol=1e-12, atol=0) and len(a1) > 5)
common = [int(s) for s in a1[:, 0] if s in set(b[:, 0])]
ia = {int(s): i for i, s in enumerate(a1[:, 0])}; ib = {int(s): i for i, s in enumerate(b[:, 0])}
edt = max(abs(a1[ia[s], 1] - b[ib[s], 1]) / b[ib[s], 1] for s in common); et = max(abs(a1[ia[s], 2] - b[ib[s], 2]) / b[ib[s], 2] for s in common)
chk('DT per step within 2e-3 of plain FDS', edt < 2e-3, '(%d steps compared, max %.2e)' % (len(common), edt))
chk('T per step within 2e-3 of plain FDS', et < 2e-3, '(max %.2e; step 1: driver %.7g, FDS %.7g)' % (et, a1[0, 2], b[0, 2]))
sys.exit(0 if ok else 1)
PY
exit $rc
