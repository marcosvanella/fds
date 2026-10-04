#!/bin/bash
# Driver-level end-to-end checks of the multi-level transport (R2b): Role 3's FluxStageRunner between the stage pieces of Role 1's TimeLoop, a level 1 bound by TimeLoop::bind_level
# (patches 0007-0009, fine-level FDS kernels) over the FDS cases of tests/cases. Run with the gfortran and the oneAPI (ifx) builds, same numbers.
# usage: run_e2e_driver.sh <driver-build> [work-dir]        Exit 0 all pass, 1 a check failed, 77 skipped (no driver / driver without the --rt-e2e hook / no bind_level).
# The driver needs the hook of docs/upstream-patches/0006-r2b-driver-e2e-hook.patch (Source/driver/main.cpp and CMakeLists.txt include DriverSources.cmake): `fds_amr <case> --rt-e2e <test>`.
# The velocity is prescribed (constant) and the pressure solve is not run: the FDS stages of one step are viscosity, ADV read-out/overrides, density, exchange 1/4, boundary, velocity flux,
# WALL_BC, DIF read-out/overrides, DIVERGENCE_PART_1, run level by level (DriverModes.cpp). Two species of equal molecular weight; the second is a passive tracer.
#   E1 full-domain level 1 equals the uniform-fine run (RHO and both ZZ, field by field), 1 and 4 ranks
#   E2 uniform state stays uniform across the interface (patch of level 1, a few steps), 1 and 4 ranks
#   E3 static two-level layout of ns2d_16_int_1to2_refinement (patch = coarse cells 4..11): composite mass and species conserved to round-off with the overwrite ON; OFF drifts (negative control)
#   E4 FR-016 ghost check through a real TimeLoop stage (needs the reference dumps: RT_FR016_DUMP=<prefix>; skipped with a note otherwise)
#   E5 D-059: shared ghost cell count printed; corner blobs (blobs astride the patch corners) with the overwrite ON
# Tolerances (set from the real runs, see notes/flux-stage-wiring.md):
#   E1 RHO difference <= 1e-12; ZZ difference <= E1_ZZ_TOL (default 1e-4, observed 4e-5): the unbound y domain edge of a 2D fine box (Role 1, see E3b below and E1c); NOT round-off, tighten when Role 1 closes it.
#   E3a advection only (species diffusivity 1e-12): overwrite ON conserves the composite species to round-off (< E3A_TOL, default 1e-12, observed 2e-14), OFF drifts (> 1e-6 and > 1e3 x ON)
#   E3b default diffusivity: overwrite ON drift < MASS_TOL (default 5e-6, observed 4e-7 after 6 steps, NOT round-off). CAUSE (found, not ours): the 2D fine box has no wall cells at its y domain
#      edge and its y ghost layer is never filled, so DIVERGENCE_PART_1 adds a spurious y diffusive flux (Role 1, level binding gap "fine-box domain-edge wall cells"); total mass exact (< 1e-12);
#      OFF drifts > 1e-5 and > 30 x ON (observed 4e-5).
#   E1c / E3c the same on a triply periodic 3D case where that edge does not exist: E1c full-domain level 1 equals the uniform fine run to E1C_TOL (default 1e-12, observed 3e-16); E3c interior fine
#      box (y range 4..11), default diffusivity: ON drift < E3A_TOL (observed 1e-14), OFF drifts (negative control, observed 1.6e-5).
set -u
BLD=${1:-}; WORK=${2:-/tmp/rt_e2e_driver}
HERE=$(cd "$(dirname "$0")" && pwd); CASES=$HERE/cases
E1_ZZ_TOL=${E1_ZZ_TOL:-1e-4}; MASS_TOL=${MASS_TOL:-5e-6}; E3A_TOL=${E3A_TOL:-1e-12}; E1C_TOL=${E1C_TOL:-1e-12}; STEPS=${E2E_STEPS:-6}
if [ -z "$BLD" ] || [ ! -x "$BLD/fds_amr" ]; then echo "SKIP all: driver executable fds_amr not found (argument 1)"; exit 77; fi
mkdir -p "$WORK" && cd "$WORK" || exit 1
cp "$CASES"/rt_blob2d_*.fds "$CASES"/rt_blob3d_*.fds . 2>/dev/null
for c in rt_blob2d_16 rt_blob2d_16_np4; do sed "s/VISCOSITY=0.002,/VISCOSITY=0.002,DIFFUSIVITY=1E-12,/" $c.fds > ${c}_nodif.fds; done   # advection-only variants (species diffusivity switched off)
export OMP_NUM_THREADS=1
run() {  # run <np> <case> <args...>  -> stdout of the run, stored in $LAST
  local np=$1 case=$2; shift 2
  LAST=$(timeout ${E2E_TIMEOUT:-900} mpirun --bind-to none --oversubscribe -np "$np" "$BLD/fds_amr" "$case.fds" --rt-e2e "$@" 2>&1)
  RC=$?
}
run 1 rt_blob2d_16 list
if ! echo "$LAST" | grep -q "RTE2E modes"; then
  echo "SKIP all: this driver has no --rt-e2e mode (apply docs/upstream-patches/0006-r2b-driver-e2e-hook.patch: Source/driver/main.cpp + CMakeLists.txt include DriverSources.cmake)"; exit 77
fi
fails=0
pass() { echo "PASS $1"; }
fail() { echo "FAIL $1"; fails=$((fails + 1)); }
val() { echo "$LAST" | sed -n "s/.*$1 \([-+0-9.eE]*\).*/\1/p" | head -1; }   # first number after the keyword
lt() { python3 -c "import sys; sys.exit(0 if float('$1') < float('$2') else 1)"; }
BLOB="--blob 0.25 0.5 0.375 0.625"          # tracer blob aligned with the coarse cells of the 16 x 16 case
DT=0.005                                      # CFL 0.24 on the 32 x 32 grid for (u, w) = (1, 0.5)
for np in 1 4; do
  S=""; MS=100000; MS1=100000; [ $np = 4 ] && { S="_np4"; MS=8; MS1=16; }   # one fine box on 1 rank; 4 ranks: one fine box per rank (several boxes on one rank need docs/upstream-patches/0007-r2b-fill-om-bounds.patch)
  C16=rt_blob2d_16$S; C32=rt_blob2d_32$S
  # ---- E1
  rm -f E1A_np$np.* E1B_np$np.*
  run $np $C32 transport --none --steps $STEPS --dt $DT --smooth 0.3 --label E1fine --dump E1A_np$np; a=$RC; echo "$LAST" | grep -E "RTE2E E1fine (levels|WORST|RHO range)" | sed 's/^/  /'
  run $np $C16 transport --full --maxsize $MS1 --steps $STEPS --dt $DT --smooth 0.3 --label E1full --dump E1B_np$np; b=$RC; echo "$LAST" | grep -E "RTE2E E1full (levels|WORST|RHO range|runner)" | sed 's/^/  /'
  if [ $a = 0 ] && [ $b = 0 ]; then
    out=$(python3 "$HERE/compare_level_dumps.py" E1A_np$np 0 E1B_np$np 1 "$E1_ZZ_TOL"); rc=$?
    echo "$out" | sed 's/^/  /'
    rho_bad=$(python3 "$HERE/compare_level_dumps.py" E1A_np$np 0 E1B_np$np 1 1e-12 | sed -n 's/^RTE2E-CMP RHO:.*cells above tol \([0-9]*\).*/\1/p')
    { [ $rc = 0 ] && [ "${rho_bad:-1}" = 0 ]; } && pass "E1 np=$np full-domain level 1 vs uniform fine run after $STEPS steps: RHO equal to 1e-12, ZZ within $E1_ZZ_TOL (known binding gap, see header)" || fail "E1 np=$np full-domain level 1 vs uniform fine run"
  else fail "E1 np=$np runs did not complete (rc $a, $b)"; echo "$LAST" | tail -5; fi
  # ---- E2
  run $np $C16 transport --patch 4 11 4 11 --maxsize $MS --no-blob --steps $STEPS --dt $DT --label E2
  echo "$LAST" | grep -E "RTE2E E2 (levels|CHECK|RHO range|runner)" | sed 's/^/  /'
  echo "$LAST" | grep -q "E2 CHECK-OK" && [ $RC = 0 ] && pass "E2 np=$np uniform state stays uniform across the interface" || fail "E2 np=$np uniform state"
  # ---- E3a advection only, E3b with diffusion; each with the negative control (overwrite OFF)
  for kind in a b; do
    CASEX=$C16; [ $kind = a ] && CASEX=${C16}_nodif
    for ov in 1 0; do
      run $np $CASEX transport --patch 4 11 4 11 --maxsize $MS --smooth 0.3 --steps $STEPS --dt $DT --overwrite $ov --label E3${kind}ov$ov
      echo "$LAST" | grep -E "RTE2E E3${kind}ov$ov (levels|WORST|composite|runner)" | sed 's/^/  /'
      d=$(val "E3${kind}ov$ov WORST-MASS-DRIFT")
      if [ $kind = a ]; then tol=$E3A_TOL; offmin=1e-6; ratio=1000; else tol=$MASS_TOL; offmin=1e-5; ratio=30; fi
      if [ $ov = 1 ]; then
        don=$d; { [ $RC = 0 ] && [ -n "$d" ] && lt "$d" "$tol"; } && pass "E3$kind np=$np composite mass and species, overwrite ON ($d < $tol)" || fail "E3$kind np=$np overwrite ON drift '$d' (limit $tol)"
      else
        { [ $RC = 0 ] && [ -n "$d" ] && [ -n "${don:-}" ] && lt "$offmin" "$d" && lt "$(python3 -c "print($ratio*max(float('$don'),1e-16))")" "$d"; } && pass "E3$kind np=$np negative control, overwrite OFF drifts ($d > $offmin and > $ratio x ON drift $don)" || fail "E3$kind np=$np overwrite OFF drift '$d' (ON drift '${don:-}')"
      fi
    done
  done
  # ---- E1c/E3c: the same two checks on a triply periodic 3D case (1 rank). The fine box of the 2D cases (one cell in y) has domain-edge faces in y that the level binding does not bind
  # (fine box without wall cells: its y ghost layer is never filled and DIVERGENCE_PART_1 differences against it), which is the whole source of the E3b residual and of the E1 ZZ gap;
  # with periodic y (E1c: full level 1) or a fine box that touches no domain edge (E3c: y range 4..11) neither exists and the checks are round-off (notes/flux-stage-wiring.md).
  if [ $np = 1 ]; then
    rm -f E1A3.* E1B3.*
    run 1 rt_blob3d_32 transport --none --steps $STEPS --dt $DT --smooth 0.3 --label E1cfine --dump E1A3; a=$RC
    run 1 rt_blob3d_16 transport --full --maxsize 100000 --steps $STEPS --dt $DT --smooth 0.3 --label E1cfull --dump E1B3; b=$RC
    if [ $a = 0 ] && [ $b = 0 ]; then
      out=$(python3 "$HERE/compare_level_dumps.py" E1A3 0 E1B3 1 "$E1C_TOL"); rc=$?; echo "$out" | sed 's/^/  /'
      [ $rc = 0 ] && pass "E1c 3D periodic: full-domain level 1 vs uniform fine run, RHO and ZZ within $E1C_TOL (the 2D E1 gap is the unbound y edge)" || fail "E1c 3D full-domain level 1 vs uniform fine run"
    else fail "E1c runs did not complete (rc $a, $b)"; fi
    for ov in 1 0; do
      run 1 rt_blob3d_16 transport --patch 4 11 4 11 --py 4 11 --maxsize 100000 --smooth 0.3 --steps $STEPS --dt $DT --velocity 1 0 --overwrite $ov --label E3cov$ov
      echo "$LAST" | grep -E "RTE2E E3cov$ov (levels|WORST|composite|runner)" | sed 's/^/  /'
      d=$(val "E3cov$ov WORST-MASS-DRIFT")
      if [ $ov = 1 ]; then
        donc=$d; { [ $RC = 0 ] && [ -n "$d" ] && lt "$d" "$E3A_TOL"; } && pass "E3c 3D interior fine box, default diffusion, overwrite ON: composite drift $d < $E3A_TOL (round-off)" || fail "E3c overwrite ON drift '$d' (limit $E3A_TOL)"
      else
        { [ $RC = 0 ] && [ -n "$d" ] && [ -n "${donc:-}" ] && lt 1e-6 "$d" && lt "$(python3 -c "print(1000*max(float('$donc'),1e-16))")" "$d"; } && pass "E3c negative control, overwrite OFF drifts ($d > 1e-6 and > 1000 x ON drift $donc)" || fail "E3c overwrite OFF drift '$d' (ON drift '${donc:-}')"
      fi
    done
  fi
  # ---- E5 corner blobs
  run $np $C16 transport --patch 4 11 4 11 --maxsize $MS --smooth 0 $BLOB --blob 0.6 0.9 0.6 0.9 --blob 0.1 0.3 0.1 0.3 --blob 0.6 0.9 0.1 0.3 --blob 0.1 0.3 0.6 0.9 --steps $STEPS --dt $DT --label E5
  echo "$LAST" | grep -E "RTE2E E5 (levels|D-059|WORST|tracer range|runner)" | sed 's/^/  /'
  [ $RC = 0 ] && echo "$LAST" | grep -q "E5 D-059 shared ghost cells" && pass "E5 np=$np corner blobs ran; D-059 count printed" || fail "E5 np=$np corner blobs"
done
# ---- E4 (1 rank, ns2d_16_l0 geometry)
if [ -n "${RT_FR016_DUMP:-}" ] && [ -e "${RT_FR016_DUMP}" ]; then
  cp "$HERE/../../driver/tests/cases/ns2d_16_l0.fds" .
  run 1 ns2d_16_l0 ghost --dumps "$RT_FR016_DUMP" --step ${RT_FR016_STEP:-3} --label E4
  echo "$LAST" | grep -E "RTE2E ghost|RTE2E E4" | sed 's/^/  /'
  [ $RC = 0 ] && echo "$LAST" | grep -q "RTE2E ghost .*covered cells written" && pass "E4 ghost cells written by a real stage_exchange agree with the FDS dumps (non-conflict cells)" || fail "E4 FR-016 ghost check through a stage"
else echo "SKIP E4 FR-016 ghost check through a stage: set RT_FR016_DUMP=<prefix of the reference dumps of ns2d_16_int_1to2_refinement> (instrumented FDS, FDSREF_FILE=<prefix> mpirun -np 13)"; fi
# ---- Role 1's two-level prescribed-velocity test with the FluxStageRunner override lists
if [ -x "$BLD/driver_unit_tests" ]; then
  for np in 1 4; do
    out=$(mpirun --bind-to none --oversubscribe -np $np "$BLD/driver_unit_tests" 2>&1 | grep -E "TWOLEVEL|^ALL PASS|^SOME FAILED")
    echo "$out" | sed "s/^/  [np=$np] /"
    echo "$out" | grep -q "ALL PASS" && pass "driver_unit_tests np=$np (two-level prescribed-velocity transport with the FluxStageRunner lists)" || fail "driver_unit_tests np=$np"
  done
fi
[ $fails = 0 ] && echo "E2E DRIVER PASS" || echo "E2E DRIVER FAIL: $fails check(s)"
[ $fails = 0 ]
