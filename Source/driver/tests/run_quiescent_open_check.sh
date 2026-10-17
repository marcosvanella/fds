#!/bin/bash
# run_quiescent_open_check.sh <build-dir> [<run-dir>]: quiescent uniform air (GVEC=0, radiation on, all sides OPEN, ground wall) must stay quiescent with fine boxes
# that touch the ground, the x/y open faces and (patch 0 15) a lateral open face; 2 levels and 3 levels, static patches.
# Cause of the earlier failure (net inflow ~190 kg per step in the 1000 m case): a hidden fine wall started its radiation intensities ILW at 0 (ALLOCATE_STORAGE),
# so the wall reflected nothing, the gas next to a ground wall saw a 0.6% intensity deficit, absorbed less than it emitted, and cooled (D < 0, mass drawn in through the open faces).
# Level 0 reaches the ambient equilibrium in its initialization sweeps. Fix: the hidden-wall ILW starts at SIGMA*TMPA^4/PI.
# Passes when |net boundary inflow| / initial mass < 1e-12 for both runs (before the fix: 2-level ~1e-6).
set -u
BLD=$1; RUN=${2:-$BLD/run/quiescent}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
run() {  # <tag> <args...>
  local t=$1; shift; mkdir -p $t; cp "$HERE/cases/quiescent_open_2level.fds" $t/c.fds
  (cd $t && OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run "$@" --blocking 4 --maxsize 8 --outdir . --chid c > run.log 2>&1)
}
run two   --levels 2 --init-patch --patch 0 3 0 3
run three --levels 3 --init-patch --patch 0 3 0 3 --patch2 1 2 1 2
r() { grep -a -m1 "LEVEL RESULT steps" $1/run.log | sed -E 's/.*composite mass ([0-9.e+]+) -> [0-9.e+]+ boundary inflow ([-0-9.e+]+) .*/\1 \2/'; }
A=$(r two); B=$(r three)
if [ -z "$A" ] || [ -z "$B" ]; then echo "FAIL quiescent_open (no RESULT line, see $RUN/*/run.log)"; exit 1; fi
echo "$A $B" | awk '{ a=$2/$1; b=$4/$3; if(a<0)a=-a; if(b<0)b=-b; ok=(a<1e-12 && b<1e-12);
  printf "%s quiescent_open: |net inflow|/mass 2-level %.3e, 3-level %.3e (< 1e-12)\n", ok?"PASS":"FAIL", a, b; exit !ok }'
