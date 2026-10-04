#!/bin/bash
# S6 decomposition independence (USE_AMREX=ON build, threads=1): the 32x1x32 periodic shunn3 problem as 1 box (cases/dec1.fds), 4 boxes of 16x16 (dec2) and 16 boxes of 8x8 (dec4),
# at 1 and 4 ranks. Kernel-facing rules (M2a): passive scalars are carried by Fields.cpp; only uniform Cartesian metrics are used.
# usage: run_decomp_check.sh <build-dir> [work-dir]
#  (1) first step, pass 1 (FDSTL_STAGE=1): the reduction-free stages RHOS ZZS (density), FVX FVZ MU (VELOCITY_FLUX), DS MU KRES TMP RSUM (DIVERGENCE_PART_1) must be
#      BITWISE equal to the 1-box run (DS at most 1 cell, <= 1e-15, for the 16-box split: see README "Decomposition")
#      with the exact zone sums (--exact-zone-sums, here env FDSTL_EXACT_ZONES=1, D-053: FDS order is the default); the default mode is checked on dec4/1 rank and dec2/4 ranks:
#      reduction-free stages bitwise, D/DS/DDDT within 1e-12 absolute of the 1-box run (zone sums in FDS order differ by rounding with the split)
#  (2) the whole run (T = 1): final U W H HS D DS RHO TMP ZZ within 2e-12 of max|field| of the 1-box run (the pressure solve sums in a layout dependent order)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; WORK=${2:-$BLD/run/decomp}; mkdir -p "$WORK"; rc=0
run() { # label case np [env...]
  local d="$WORK/$1"; rm -rf "$d"; mkdir -p "$d"; cp "$HERE/cases/$2.fds" "$d/"
  (cd "$d" && env "${@:4}" timeout 900 mpirun --bind-to none --oversubscribe -np $3 "$BLD/fds_amr" "$2.fds" --run --outdir . --chid $2 --quiet > stdout.txt 2> stderr.txt) || { echo "FAIL run $1"; rc=1; }
}
for c in "dec1 dec1 1" "dec2 dec2 1" "dec4 dec4 1" "dec2 dec2_np4 4" "dec4 dec4_np4 4"; do set -- $c; run "s_$1_np$3" $2 $3 FDSTL_STAGE=1 FDSTL_EXACT_ZONES=1 ; run "f_$1_np$3" $2 $3; done
# D-053 default mode (zone sums in FDS order): the same first-step stages, not exact for D/DS/DDDT (the zone sums change with the split by rounding)
run "sd_dec4_np1" dec4 1 FDSTL_STAGE=1 ; run "sd_dec2_np4" dec2_np4 4 FDSTL_STAGE=1
grep -q 'EXTERNAL_GHOSTS_FILLED=1' "$WORK/s_dec2_np1/stdout.txt" && echo 'EXTERNAL_GHOSTS_FILLED path exercised (flag set by the driver)' || { echo 'FAIL EXTERNAL_GHOSTS_FILLED not set'; rc=1; }
python3 - "$WORK" <<'P' || rc=1
import numpy as np, sys, glob, os
w = sys.argv[1]; bad = 0
ref_s = w + '/s_dec1_np1'; ref_f = w + '/f_dec1_np1'
exact = ['p1_dens_RHOS', 'p1_dens_ZZS', 'p1_vflux_FVX', 'p1_vflux_FVZ', 'p1_vflux_MU', 'p1_div1_MU', 'p1_div1_KRES', 'p1_div1_TMP', 'p1_div1_RSUM']
near = ['p1_div1_DS', 'p1_div_DS', 'p1_div_D', 'p1_div_DDDT']
for run in ['dec2_np1', 'dec4_np1', 'dec2_np4', 'dec4_np4']:
    d = w + '/s_' + run
    nb = [t for t in exact if not np.array_equal(np.fromfile('%s/stage_%s.bin' % (ref_s, t)), np.fromfile('%s/stage_%s.bin' % (d, t)))]
    nn = []
    for t in near:
        a = np.fromfile('%s/stage_%s.bin' % (ref_s, t)); b = np.fromfile('%s/stage_%s.bin' % (d, t))
        nd = int((a != b).sum())
        # tolerance: 2e-15 absolute, or 1 ulp of the larger value where that is bigger (D/DDDT reach |x| > 16 where 1 ulp = 3.6e-15)
        tol = np.maximum(2e-15, np.spacing(np.maximum(np.abs(a), np.abs(b))))
        if nd > 1 or (np.abs(a - b) > tol).any(): nn.append((t, nd, float(np.abs(a - b).max())))
    print('%-9s reduction-free stages: %d of %d bitwise%s; DS/D/DDDT: %s' % (run, len(exact) - len(nb), len(exact), (' DIFFER ' + str(nb)) if nb else '', 'within 1 cell / max(2e-15, 1 ulp)' if not nn else 'FAIL ' + str(nn)))
    bad += len(nb) + len(nn)
    worst = 0.0
    for n in ['U', 'W', 'H', 'HS', 'D', 'DS', 'RHO', 'TMP', 'ZZ2']:
        a = np.fromfile(glob.glob('%s/*_final_%s.bin' % (ref_f, n))[0]); b = np.fromfile(glob.glob('%s/f_%s/*_final_%s.bin' % (w, run, n))[0])
        worst = max(worst, np.abs(a - b).max() / max(np.abs(a).max(), 1e-300))
    print('          final fields vs 1 box: max relative difference %.2e %s' % (worst, 'ok' if worst < 2e-12 else 'FAIL'))
    if worst >= 2e-12: bad += 1
for run in ['dec4_np1', 'dec2_np4']:
    d = w + '/sd_' + run
    nb = [t for t in exact if not np.array_equal(np.fromfile('%s/stage_%s.bin' % (ref_s, t)), np.fromfile('%s/stage_%s.bin' % (d, t)))]
    mx = max(float(np.abs(np.fromfile('%s/stage_%s.bin' % (ref_s, t)) - np.fromfile('%s/stage_%s.bin' % (d, t))).max()) for t in near)
    ok = not nb and mx < 1e-12
    print('%-9s default (FDS-order zone sums): reduction-free stages %d of %d bitwise; D/DS/DDDT max abs difference %.2e %s' % (run, len(exact) - len(nb), len(exact), mx, 'ok' if ok else 'FAIL'))
    bad += 0 if ok else 1
sys.exit(1 if bad else 0)
P
[ $rc = 0 ] && echo "DECOMPOSITION CHECK PASS" || echo "DECOMPOSITION CHECK FAIL"
exit $rc
