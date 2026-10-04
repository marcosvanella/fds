#!/bin/bash
# Step sequence (T, DT per step) of the driver against SINGLE-MESH plain FDS (ground truth: a multi-mesh run must equal the single-mesh run).
# Reference: $BASELINE/shunn3_4mesh_32__1mesh/*_steps.csv (the same domain as one 32x1x32 mesh; FDS writes the diagnostic steps only, so steps are matched by number).
# Driver runs, first 40 steps: the 4-mesh input shunn3_4mesh_32 at 1 and 4 ranks and the 1-mesh input (the reference .fds) at 1 rank.
# Checks: (1) the three driver runs give the same T/DT sequence; (2) DT per compared step within 6e-3 of the single-mesh FDS (3 printed digits, i.e. rounding of the .csv); (3) T within 5e-4.
# Known reference quirk (notes/multi-mesh-step-sequence.md): plain FDS on the 2x2 four-mesh input has a step-1 DT larger than single-mesh FDS (T after step 1 6.195e-4 s
# against 6.188e-4 s). That multi-mesh FDS run is NOT used as reference here. Candidate cause (unconfirmed): interpolated-boundary UVW_SAVE data in the first trial CFL pass.
# usage: run_step_sequence_check.sh <build-dir (absolute)> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/seqchk}; c=shunn3_4mesh_32; R=$BASELINE/${c}__1mesh; rc=0; rm -rf "$W"; mkdir -p "$W/m4np1" "$W/m4np4" "$W/m1np1"
run() { # dir input np
  cp "$2" "$1/"; (cd "$1" && timeout 1500 mpirun --bind-to none --oversubscribe -np $3 "$B/fds_amr" $c.fds --run --steps 40 --outdir . --chid $c --quiet > stdout.txt 2> stderr.txt) || { echo "RUN FAILED $1"; rc=1; }
}
run "$W/m4np1" "$BASELINE/$c/$c.fds" 1; run "$W/m4np4" "$BASELINE/$c/$c.fds" 4; run "$W/m1np1" "$R/$c.fds" 1
python3 - "$W" "$R/${c}_steps.csv" $c <<'PY' || rc=1
import sys, numpy as np
W, ref, c = sys.argv[1:4]
rd = lambda f: np.genfromtxt(f, delimiter=',', skip_header=2, usecols=(0, 2, 3))
a = {k: rd(f'{W}/{k}/{c}_steps.csv') for k in ('m4np1', 'm4np4', 'm1np1')}; b = rd(ref)
ok = True
def chk(name, cond, msg=''):
    global ok; print(('PASS ' if cond else 'FAIL ') + name + (' ' + msg if msg else '')); ok = ok and cond
chk('4-mesh np1, 4-mesh np4 and 1-mesh driver runs give the same step sequence', all(a[k].shape == a['m1np1'].shape and np.allclose(a[k], a['m1np1'], rtol=1e-9, atol=0) for k in a) and len(a['m1np1']) > 5)
ib = {int(s): i for i, s in enumerate(b[:, 0])}
for k in a:
    ia = {int(s): i for i, s in enumerate(a[k][:, 0])}; common = [s for s in ia if s in ib]
    edt = max(abs(a[k][ia[s], 1] - b[ib[s], 1]) / b[ib[s], 1] for s in common); et = max(abs(a[k][ia[s], 2] - b[ib[s], 2]) / b[ib[s], 2] for s in common)
    chk(k + ' against single-mesh FDS', edt < 6e-3 and et < 5e-4, '(%d steps, max DT diff %.2e, max T diff %.2e; step 1 T %.7g, FDS %.7g)' % (len(common), edt, et, a[k][0, 2], b[0, 2]))
sys.exit(0 if ok else 1)
PY
exit $rc
