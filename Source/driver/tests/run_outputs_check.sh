#!/bin/bash
# S7 writer check: the FDS output files written by the driver (CHID_devc.csv, _hrr.csv, _mass.csv, _steps.csv, _cpu.csv, .out; main.f90 patch 0006) on the M2a cases.
# usage: [REF_FDS=<same-compiler USE_AMREX=OFF fds> [REF_MPI_OPTS=<opts>]] [MPI_OPTS=<opts>] run_outputs_check.sh <driver-build> <work>
# Same-compiler reference (optional): the gfortran baseline is only a valid reference for a gfortran driver build; for any other compiler (oneAPI, ...) the 1e-8 _hrr.csv gate and the
# bitwise _mass.csv / pressure-iteration gates compare against round-off that is the compiler's, not the driver's. With REF_FDS set, the shunn3_32 and csmag_32 references are
# produced here by running that USE_AMREX=OFF fds (same compiler and flags as the driver build, 1 rank) on the same inputs in <work>/ref/, and every gate below compares against them.
# Without REF_FDS nothing changes (the $BASELINE / $CSMAG_DIR references are used). MPI_OPTS (default '--bind-to none --oversubscribe', Open MPI) is the mpirun option string for every
# run, set it to '--bind-to none' for Intel MPI (Hydra has no --oversubscribe). REF_MPI_OPTS defaults to MPI_OPTS.
# Checks (no tolerance is loosened for the physics gates; the HRR tolerance below is only for the writer check):
#   shunn3_32 (1 rank): _mass.csv bitwise equal to the baseline; _hrr.csv time column and all columns equal to the baseline within 1e-8 of the column scale (the
#       full-step solution differs from FDS at round-off, T2); .out pressure-iteration lines equal to the baseline; _steps.csv has the same number of rows as the baseline.
#   shunn3_4mesh_32 (1 and 4 ranks): _mass.csv and _hrr.csv of the two rank counts equal (mass bitwise, hrr within 1e-8 of the scale); total mass constant to 1e-12.
#   csmag_32 (1 rank, scratch input, no official baseline): _devc.csv has the KE device rows, t=0 value equal to the scratch FDS reference's t=0 value.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=$2; rc=0
MPI_OPTS=${MPI_OPTS---bind-to none --oversubscribe}; REF_MPI_OPTS=${REF_MPI_OPTS:-$MPI_OPTS}; REF_FDS=${REF_FDS:-}
CSDIR=${CSMAG_DIR:-/workspace/fds-amr/scratch/role1-s3-work/ref_runs/csmag_32}
run() { # dir src np
  rm -rf "$1"; mkdir -p "$1"; local c; c=$(basename "$2" .fds); cp "$2" "$1/"
  for f in "$(dirname "$2")"/*uvw*.csv; do [ -f "$f" ] && cp "$f" "$1/"; done
  (cd "$1" && mpirun $MPI_OPTS -np "$3" "$B/fds_amr" "$c.fds" --run --outdir . --chid "$c" --quiet > stdout.txt 2> stderr.txt) || { echo "RUN FAILED $1"; rc=1; }
}
BLREF=$BASELINE; CSREF=$CSDIR
if [ -n "$REF_FDS" ]; then   # same-compiler reference runs (USE_AMREX=OFF fds), OMP_NUM_THREADS from the caller (1 in the oneAPI set-up)
  refrun() { # case dir-with-inputs
    local d="$W/ref/$1"; rm -rf "$d"; mkdir -p "$d"; cp "$2/$1.fds" "$d/"; for f in "$2"/*uvw*.csv; do [ -f "$f" ] && cp "$f" "$d/"; done
    (cd "$d" && mpirun $REF_MPI_OPTS -np 1 "$REF_FDS" "$1.fds" > stdout.txt 2> stderr.txt) || { echo "REFERENCE RUN FAILED $d"; rc=1; }
  }
  refrun shunn3_32 "$BASELINE/shunn3_32"; refrun csmag_32 "$CSDIR"
  BLREF=$W/ref; CSREF=$W/ref/csmag_32
  echo "INFO reference: same-compiler FDS $REF_FDS (runs in $W/ref); the $BASELINE / $CSDIR references are not used for the gates"
fi
run "$W/sh1" "$BASELINE/shunn3_32/shunn3_32.fds" 1
run "$W/m1" "$BASELINE/shunn3_4mesh_32/shunn3_4mesh_32.fds" 1
run "$W/m4" "$BASELINE/shunn3_4mesh_32/shunn3_4mesh_32.fds" 4
run "$W/cs" "$CSDIR/csmag_32.fds" 1
python3 - "$W" "$BLREF" "$CSREF" <<'PY' || rc=1
import sys, re, numpy as np
W, BL, CS = sys.argv[1:4]
ok = True
def rd(f): return np.genfromtxt(f, delimiter=',', skip_header=2)
def chk(name, cond, msg=''):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name + (' ' + msg if msg else ''))
    ok = ok and cond
def scaled(a, b):
    return max((np.abs(a[:, j] - b[:, j]).max() / max(np.abs(b[:, j]).max(), 1e-300) if np.abs(b[:, j]).max() > 0 else np.abs(a[:, j]).max()) for j in range(a.shape[1]))
c = 'shunn3_32'
chk(c + ' mass.csv bitwise', open(f'{W}/sh1/{c}_mass.csv').read() == open(f'{BL}/{c}/{c}_mass.csv').read())
a, b = rd(f'{W}/sh1/{c}_hrr.csv'), rd(f'{BL}/{c}/{c}_hrr.csv')
chk(c + ' hrr.csv shape+time', a.shape == b.shape and np.array_equal(a[:, 0], b[:, 0]))
e = scaled(a, b); chk(c + ' hrr.csv columns within 1e-8 of scale', e < 1e-8, '(max %.2e)' % e)
pi = lambda f: re.findall(r'Pressure Iterations: (\d+)', open(f).read())
chk(c + ' .out pressure iterations equal baseline', pi(f'{W}/sh1/{c}.out') == pi(f'{BL}/{c}/{c}.out') and len(pi(f'{W}/sh1/{c}.out')) > 0)
n = sum(1 for _ in open(f'{W}/sh1/{c}_steps.csv')) - 2
nb = sum(1 for _ in open(f'{BL}/{c}/{c}_steps.csv')) - 2
chk(c + ' _steps.csv rows = baseline rows (FDS writes the diagnostic steps only)', n == nb, '(%d rows, baseline %d)' % (n, nb))
chk(c + ' _cpu.csv present', len(open(f'{W}/sh1/{c}_cpu.csv').read().splitlines()) >= 2)
c = 'shunn3_4mesh_32'
a, b = rd(f'{W}/m1/{c}_mass.csv'), rd(f'{W}/m4/{c}_mass.csv')
chk(c + ' mass.csv np1 == np4 (bitwise)', open(f'{W}/m1/{c}_mass.csv').read() == open(f'{W}/m4/{c}_mass.csv').read())
chk(c + ' total mass constant to 1e-12', np.abs(a[:, 1] - a[0, 1]).max() < 1e-12, '(max %.2e)' % np.abs(a[:, 1] - a[0, 1]).max())
a, b = rd(f'{W}/m1/{c}_hrr.csv'), rd(f'{W}/m4/{c}_hrr.csv')
e = scaled(a, b); chk(c + ' hrr.csv np1 vs np4 within 1e-8 of scale', e < 1e-8, '(max %.2e)' % e)
c = 'csmag_32'
d = open(f'{W}/cs/{c}_devc.csv').read().splitlines()
ref = open(f'{CS}/{c}_devc.csv').read().splitlines()
chk(c + ' devc.csv header (KE) and rows', d[1].split(',')[1] == 'KE' and len(d) == len(ref), '(%d rows, ref %d)' % (len(d), len(ref)))
chk(c + ' devc.csv t=0 equals the FDS reference t=0', d[2] == ref[2])
v = float(d[-1].split(',')[1]); r = float(ref[-1].split(',')[1])
print('INFO csmag_32 KE at T_END driver %.7E scratch reference %.7E (rel %.1e; the reference is not an official baseline and its pressure solve differs, see README S6b)' % (v, r, abs(v - r) / r))
print('OUTPUTS CHECK PASS' if ok else 'OUTPUTS CHECK FAIL')
sys.exit(0 if ok else 1)
PY
exit $rc
