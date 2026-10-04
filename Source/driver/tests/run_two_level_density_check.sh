#!/bin/bash
# Variable-density two-level run with discriminating controls (D-050, vv/test-plan.md 5.12.6): a heavy tracer blob (density contrast 1.6, MW 44 in MW 28) is advected across the
# coarse-fine interfaces of a fine patch (coarse cells 4..11). `fds_amr <case> --two-level-run`:
#   ON  (interface flux overwrite): composite mass and rho*Z of both species conserved to round-off                                 -> relative change <= TOL_ON (default 1e-12)
#   OFF (--no-overwrite, negative control): the conservation drift returns                                                             -> relative change of mass or a rho*Z >= MIN_OFF (default 1e-5)
#   REF (patch = the whole domain, i.e. a uniform fine run): max|div u - D| of the last step. With the FFT-type Poisson solve and no pressure iteration the baroclinic term is lagged, so
#   div u - D is not round-off for a variable-density flow (single-level FDS has the same); the check is that the two-level value is within a factor DIV_FACTOR (default 3) of the uniform fine run.
# Runs on 1 rank (blob_16_l0.fds) and on 4 ranks (blob_16_4m.fds, one level-0 mesh per rank).
# usage: run_two_level_density_check.sh <build-dir> [work-dir] [steps1=100] [steps4=40]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; WORK=${2:-$BLD/run/density}; S1=${3:-100}; S4=${4:-40}
TOL_ON=${TOL_ON:-1e-12}; MIN_OFF=${MIN_OFF:-1e-5}; DIV_FACTOR=${DIV_FACTOR:-3}
rm -rf "$WORK"; mkdir -p "$WORK"; cp "$HERE/cases/blob_16_l0.fds" "$HERE/cases/blob_16_4m.fds" "$WORK/"; cd "$WORK"; rc=0
run() { # tag np case steps patch args...
  local tag=$1 np=$2 c=$3 st=$4 pt=$5; shift 5; mkdir -p out_$tag
  timeout 3000 mpirun --bind-to none --oversubscribe -np $np "$BLD/fds_amr" $c --two-level-run --steps $st --patch $pt --outdir out_$tag --chid $tag "$@" > $tag.log 2>&1
}
maxabs() { grep -a "TWO-LEVEL RESULT" $1 | grep -a "relative change" | sed 's/.*relative change \([-0-9.e+]*\).*/\1/' | python3 -c 'import sys; print(max(abs(float(x)) for x in sys.stdin))'; }
lastdiv() { grep -a "TWO-LEVEL step" $1 | tail -1 | sed 's/.*max|div u - D| \([-0-9.e+]*\) .*/\1/'; }
check() { # name np case steps
  local n=$1 np=$2 c=$3 st=$4
  run ${n}_on $np $c $st "4 11 4 11"; run ${n}_off $np $c $st "4 11 4 11" --no-overwrite; run ${n}_ref $np $c $st "0 15 0 15"
  local on off dv dref
  on=$(maxabs ${n}_on.log); off=$(maxabs ${n}_off.log); dv=$(lastdiv ${n}_on.log); dref=$(lastdiv ${n}_ref.log)
  echo "  $n: steps $st, ON drift $on, OFF drift $off, last-step max|div u - D| two-level $dv, uniform fine $dref"
  python3 - "$on" "$off" "$dv" "$dref" "$TOL_ON" "$MIN_OFF" "$DIV_FACTOR" "$n" <<'PY' || rc=1
import sys
on,off,dv,dref,tol,mn,fac=map(float,sys.argv[1:8]); n=sys.argv[8]
ok=True
if not on<=tol: print("FAIL %s: overwrite ON drift %g > %g"%(n,on,tol)); ok=False
if not off>=mn: print("FAIL %s: overwrite OFF control drift %g < %g (the control does not discriminate)"%(n,off,mn)); ok=False
if not dv<=fac*dref: print("FAIL %s: max|div u - D| %g > %g x uniform fine %g"%(n,dv,fac,dref)); ok=False
if ok: print("PASS %s: ON %g <= %g, OFF %g >= %g, div %g within %g x %g"%(n,on,tol,off,mn,dv,fac,dref))
sys.exit(0 if ok else 1)
PY
}
check np1 1 blob_16_l0.fds $S1
check np4 4 blob_16_4m.fds $S4
exit $rc
