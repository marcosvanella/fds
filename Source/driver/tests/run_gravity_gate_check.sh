#!/bin/bash
# run_gravity_gate_check.sh <build-dir> [<run-dir>]: quiescent tall open domain under gravity (default STRATIFICATION): the composite air balance must close with fine levels.
# Case cases/gravity_tall_2level.fds (1000x1000x400 m, 40x40x16, all sides OPEN, no fire, no species), 20 steps, static patches, 2 OMP threads each.
# Fine children are injected from the parent (CellTransfer), so a fine box needs the parent's background profile (P_0, TMP_0, RHO_0, PBAR piecewise constant per parent cell);
# with its own profile at the fine cell centres the injected rho is off the fine hydrostatic state and the balance drifts with the fine cell count (scratch/role3-walls: fine_strat_background.patch).
# Before the patch: 2-level residual/limit 4.1e-3 (GVEC=0 or STRATIFICATION=F: 1e-10), 3-level 1.24 (gate FAIL). After: 2-level 1.5e-5, 3-level 0.076.
# Passes when 2-level residual/limit < 2e-4 and 3-level < 0.3.
set -u
BLD=$1; RUN=${2:-$BLD/run/gravity}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
run() {  # <tag> <args...>
  local t=$1; shift; mkdir -p $t; cp "$HERE/cases/gravity_tall_2level.fds" $t/c.fds
  (cd $t && OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run "$@" --blocking 4 --maxsize 8 --outdir . --chid c > run.log 2>&1)
}
run two   --levels 2 --init-patch --patch 2 5 2 5
run three --levels 3 --init-patch --patch 2 5 2 5 --patch2 3 4 3 4
r() { grep -a -m1 "composite rho\*Z1" $1/run.log | sed -E 's/.*\(\/ limit ([-0-9.e+]+)\).*/\1/'; }
A=$(r two); B=$(r three)
if [ -z "$A" ] || [ -z "$B" ]; then echo "FAIL gravity_gate (no RESULT line, see $RUN/*/run.log)"; exit 1; fi
awk -v a=$A -v b=$B 'BEGIN{ if(a<0)a=-a; if(b<0)b=-b; ok=(a<2e-4 && b<0.3);
  printf "%s gravity_gate: 2-level residual/limit %.3e (< 2e-4), 3-level %.3e (< 0.3)\n", ok?"PASS":"FAIL", a, b; exit !ok }'
