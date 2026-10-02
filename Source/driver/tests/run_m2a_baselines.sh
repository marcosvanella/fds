#!/bin/bash
# M2a gate items against the official baselines (S8): csmag_32 plain and periodic variant (+ their __sf17 copies), shunn3_4mesh_32 against the GLMAT 4-mesh baseline.
# usage: run_m2a_baselines.sh <driver-build> <work>
# The csmag_32 baseline directories keep no restart file, so the reference binary (vv-runs/refbin) is rerun on the same input in <work>/ref; the rerun's _devc.csv is
# checked byte-identical to the baseline's before its restart file is used. The driver refuses SOLVER='GLMAT' (its pressure solve is the level FFT), so the
# GLMAT baseline is compared with a driver run of the same input without that &PRES line (the default-solver input, directory shunn3_4mesh_32).
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=$2; BL=$BASELINE; REF=${REFBIN:-$(dirname "$BL")/../refbin/gnu_ompi_firex-36975d7/fds}; REF=$(readlink -f "$REF"); rc=0
[ -x "$REF" ] || { echo "reference binary not found: $REF"; exit 2; }
mkdir -p "$W"
drv() { # outdir srcdir np
  local out=$1 src=$2 np=$3 f; f=$(basename "$(ls "$src"/*.fds | head -1)")
  rm -rf "$out"; mkdir -p "$out"; cp "$src/$f" "$out/"
  for c in "$src"/*.csv; do case $c in *_hrr.csv|*_mass.csv|*_mms.csv|*_cpu.csv|*_steps.csv|*_devc.csv|*_uvw_t*_m*.csv) ;; *) cp "$c" "$out/";; esac; done
  (cd "$out" && mpirun --bind-to none --oversubscribe -np "$np" "$B/fds_amr" "$f" --run --outdir . --chid "${f%.fds}" --quiet > stdout.txt 2> stderr.txt) || { echo "DRIVER RUN FAILED $out"; rc=1; }
}
refrun() { # name
  local d=$W/ref/$1; rm -rf "$d"; mkdir -p "$d"; cp "$BL/$1/csmag_32.fds" "$BL/$1/cbc32_uvw.csv" "$d/"
  (cd "$d" && mpirun --bind-to none -np 1 "$REF" csmag_32.fds > stdout.txt 2> stderr.txt) || { echo "REFERENCE RUN FAILED $1"; rc=1; }
  cmp -s "$d/csmag_32_devc.csv" "$BL/$1/csmag_32_devc.csv" || { echo "reference rerun differs from the baseline devc: $1"; rc=1; }
}
for c in csmag_32 csmag_32__fishpak_bc000 csmag_32__sf17 csmag_32__sf17_fishpak_bc000; do refrun $c; drv "$W/run/$c" "$BL/$c" 1; done
echo "== csmag_32 periodic variant (FISHPAK_BC=0,0,0): same discretisation as the driver; KE rows equal to the printed 8 digits"
python3 "$HERE/m2a_compare.py" csmag "$W/run/csmag_32__fishpak_bc000" "$BL/csmag_32__fishpak_bc000" csmag_32 --restart "$W/ref/csmag_32__fishpak_bc000/csmag_32_1.restart" | tail -n 16 || rc=1
echo "== csmag_32 periodic variant, SIG_FIGS=17 copy: KE rows to 1e-12"
python3 "$HERE/m2a_compare.py" csmag "$W/run/csmag_32__sf17_fishpak_bc000" "$BL/csmag_32__sf17_fishpak_bc000" csmag_32 --tol-ke 1e-12 --restart "$W/ref/csmag_32__sf17_fishpak_bc000/csmag_32_1.restart" | grep -E "KE relative|COMPARE" || rc=1
echo "== csmag_32 plain baseline (coupled-Dirichlet iteration in FDS): T2 form against the baseline's own solver spread"
python3 "$HERE/m2a_compare.py" csmag-t2 "$W/run/csmag_32" "$BL/csmag_32" csmag_32 --alt "$BL/csmag_32__fishpak_bc000" --restart "$W/ref/csmag_32/csmag_32_1.restart" --alt-restart "$W/ref/csmag_32__fishpak_bc000/csmag_32_1.restart" || rc=1
echo "== csmag_32 plain baseline, SIG_FIGS=17 copy"
python3 "$HERE/m2a_compare.py" csmag-t2 "$W/run/csmag_32__sf17" "$BL/csmag_32__sf17" csmag_32 --alt "$BL/csmag_32__sf17_fishpak_bc000" --restart "$W/ref/csmag_32__sf17/csmag_32_1.restart" --alt-restart "$W/ref/csmag_32__sf17_fishpak_bc000/csmag_32_1.restart" | grep -E "KE series|COMPARE" || rc=1
for np in 1 4; do
  drv "$W/run/glmat_np$np" "$BL/shunn3_4mesh_32" $np
  echo "== shunn3_4mesh_32 vs the GLMAT 4-mesh baseline, $np rank(s)"
  python3 "$HERE/m2a_compare.py" glmat "$W/run/glmat_np$np" "$BL/shunn3_4mesh_32__glmat" shunn3_4mesh_32 || rc=1
done
[ $rc = 0 ] && echo "M2A BASELINES PASS" || echo "M2A BASELINES FAIL"
exit $rc
