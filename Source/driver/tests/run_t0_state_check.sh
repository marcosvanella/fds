#!/bin/bash
# run_t0_state_check.sh <fds_amr binary> [<run-dir>]: t=0 divergence D of a fine level (TwoLevelRun t0_stock_rho_levels / t0_divergence_levels).
# cases/t0_state.fds (hot &INIT block on the floor, open top) as 2 levels with level 1 over the whole domain, and as 1 level at the level-1 cell size (IJK doubled),
# with radiation off and on. D needs the stock start state of the fine level: RHO from its own P_0 profile, U = U_WIND plus INITIAL_NOISE, US = U before the noise,
# MU from COMPUTE_VISCOSITY before the time-step ghost rules, then the initial radiation and DIVERGENCE_PART_1(T_BEGIN) in corrector form.
# Pass: sum, norm1, min, max of the level-1 D (FDSTL_T0RAD) within 1e-12 of max|D| of the 1-level ones (bitwise apart from summation order).
# The hot block lies on the floor: the floor cells under it need the gas-side wall flux of the hidden fine walls (FDS_FINE_WALLS_QCON; without it -0.0122 there).
set -u
BIN=$1; RUN=${2:-$PWD/t0_state_run}; SRC=$(cd "$(dirname "$0")" && pwd)/cases/t0_state.fds
mkdir -p "$RUN" && cd "$RUN" || exit 1
fail=0
for rad in off on; do
    sed 's/IJK=8,8,4/IJK=16,16,8/' "$SRC" > r1.fds; cp "$SRC" r2.fds
    [ $rad = off ] && for f in r1.fds r2.fds; do sed -i 's#^&TAIL /#\&RADI RADIATION=.FALSE. /\n\&TAIL /#' $f; done
    for L in 1 2; do
        FDSTL_T0RAD=1 OMP_NUM_THREADS=1 timeout 600 "$BIN" r$L.fds --two-level-run --levels $L --patch 0 7 0 3 --maxsize 16 --nofine -1 --steps 1 --outdir . --chid r$L > r${L}_$rad.log 2>&1 \
            || { echo "FAIL t0_state: radiation $rad, $L-level run rc (see $RUN/r${L}_$rad.log)"; exit 1; }
    done
    rm -f ./*.s3d ./*.smv ./*.sf ./*.bf ./*.prt5 ./*.csv
    a=$(grep -a "^T0RAD L0 D " r1_$rad.log); b=$(grep -a "^T0RAD L1 D " r2_$rad.log)
    echo "radiation $rad 1L: ${a#T0RAD L0 D }"; echo "radiation $rad 2L: ${b#T0RAD L1 D }"
    python3 - "$a" "$b" <<'P' || fail=1
import sys
a, b = sys.argv[1].split(), sys.argv[2].split()
if len(a) < 11 or len(b) < 11: print("  missing T0RAD line"); sys.exit(1)
x = [float(a[i]) for i in (4, 6, 8, 10)]; y = [float(b[i]) for i in (4, 6, 8, 10)]
s = max(abs(x[2]), abs(x[3])); e = max(abs(p - q) for p, q in zip(x, y)) / s
print("  max diff / max|D| %.3e" % e); sys.exit(0 if e <= 1e-12 else 1)
P
done
[ $fail -eq 0 ] && echo "PASS t0_state" || echo "FAIL t0_state (see $RUN)"
exit $fail
