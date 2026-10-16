#!/bin/bash
# run_burner_inject_check.sh <build-dir> [<run-dir>]: fuel injected by the prescribed MASS_FLUX burner must not depend on where the fine levels sit.
# KNOWN ISSUE - EXPECTED TO FAIL on 215c55bbd8 (near case): see scratch/role3-walls/BURNER-NEAR-EVIDENCE.md.
# Case cases/burner_obst_3level.fds (= scratch/m1/b3 simple_test: 3.6x2.4x2.4 m, 36x24x24, MASS_FLUX(1)=0.02 on a 0.4x0.4 burner over an OBST, OPEN vent at x=3.6).
# Runs 30 steps (t = 3.0919 s) of  (near) --patch 0 7 0 7 --patch2 0 7 0 3   and  (far) --patch2 0 7 4 7  [and (two) --levels 2 --patch 0 7 0 7 with TWO=1], 1 OMP thread each, in parallel.
# Injected mass = "composite rho*Z2 0 -> X" of the RESULT line (propane mass gained by the composite, equals the counted burner inflow when the gate passes).
# Passes when |near/far - 1| < 2% and near and far each within 2% of the real Fortran FDS value REF_Z2 below.
# REF_Z2 = 7.5858908e-3 kg: real FDS (work-species/build_fds/fds, single mesh, same simple_test.fds with T_END=4.1, no AMR, 40 steps with dt = 0.1030633 .. ), the DEVC
#   &DEVC ID='m_prop', XB=0,3.6,0,2.4,0,2.4, QUANTITY='MASS FRACTION', SPEC_ID='PROPANE', STATISTICS='MASS INTEGRAL'  read at t = 3.0918990 s (= AMR step 30, same dt sequence).
#   (at t = 4.1 s, step 40, the same device reads 1.0790492e-2; nominal 3.2e-3 kg/s * t = 0.0132 kg is higher because MASS_FLUX ramps up with the default TAU_MF.)
# Measured at 215c55bbd8: 30 steps near 5.553e-3, far 7.611e-3, ref 7.586e-3 (~2 min wall, 2 runs in parallel); 40 steps near 5.456e-3, far 1.0873e-2, two-level 1.0862e-2 (ref 1.0790e-2). Cause: scratch/role3-walls/BURNER-NEAR-EVIDENCE.md.
set -u
BLD=$1; RUN=${2:-$BLD/run/burnerinject}
STEPS=${STEPS:-30}; REF_Z2=${REF_Z2:-7.5858908e-3}; TOL=${TOL:-0.02}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
run() {  # <tag> <args...>
  local t=$1; shift; mkdir -p $t; cp "$HERE/cases/burner_obst_3level.fds" $t/simple_test.fds
  (cd $t && OMP_NUM_THREADS=1 "$BLD/fds_amr" simple_test.fds --two-level-run "$@" --steps $STEPS --blocking 4 --maxsize 8 --outdir . --chid c > run.log 2>&1)
}
run near --levels 3 --patch 0 7 0 7 --patch2 0 7 0 3 &
run far  --levels 3 --patch 0 7 0 7 --patch2 0 7 4 7 &
[ "${TWO:-0}" = 1 ] && { run two --levels 2 --patch 0 7 0 7 & }
wait
z() { grep -a -m1 "composite rho\*Z2" $1/run.log | sed -E 's/.*rho\*Z2 0 -> ([-0-9.e+]+) .*/\1/'; }
N=$(z near); F=$(z far)
if [ -z "$N" ] || [ -z "$F" ]; then echo "FAIL burner_inject (no RESULT line, see $RUN/*/run.log)"; exit 1; fi
rc=0
msg=$(awk -v n=$N -v f=$F -v r=$REF_Z2 -v t=$TOL 'BEGIN{ a=(n/f-1); b=(f/r-1); c=(n/r-1); if(a<0)a=-a; if(b<0)b=-b; if(c<0)c=-c;
  printf "near %.5e far %.5e ref %.5e: |near/far-1| %.3f |far/ref-1| %.3f |near/ref-1| %.3f (tol %.2f)", n,f,r,a,b,c,t; exit !(a<t && b<t && c<t) }') || rc=1
if [ "${TWO:-0}" = 1 ]; then T2=$(z two); msg="$msg; two-level $T2"; fi
if [ $rc = 0 ]; then echo "PASS burner_inject: $msg"; else echo "FAIL burner_inject (KNOWN ISSUE, BURNER-NEAR-EVIDENCE.md): $msg"; fi
exit $rc
