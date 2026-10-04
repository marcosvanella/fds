#!/bin/bash
# M2a gate items against the official baselines (S8): csmag_32 plain and periodic variant (+ their __sf17 copies), shunn3_4mesh_32 against the GLMAT 4-mesh baseline.
# usage: run_m2a_baselines.sh <driver-build> <work>
# The csmag_32 baseline directories keep no restart file, so the reference binary (vv-runs/refbin) is rerun on the same input in <work>/ref; the rerun's _devc.csv is
# checked byte-identical to the baseline's before its restart file is used. The driver refuses SOLVER='GLMAT' (its pressure solve is the level FFT), so the
# GLMAT baseline is compared with a driver run of the same input without that &PRES line (the default-solver input, directory shunn3_4mesh_32).
# Reference-rerun check: the stored baselines were produced by gfortran. REF_COMPILER names the compiler of the reference binary (default gfortran, the Open MPI reference). With
# REF_COMPILER=gfortran the rerun's csmag_32_devc.csv must be byte-identical to the baseline's (hard failure). With any other value (e.g. REF_COMPILER=ifx for the oneAPI reference,
# REFBIN=<that fds>) the same comparison is still made but only reported (INFO line with the largest relative difference over the numeric columns, relative to the column maximum) and
# never fails: round-off differs between compilers. Only the compiler of REFBIN is declared; the driver under test is compared with that same reference's restart file either way.
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
REF_COMPILER=${REF_COMPILER:-gfortran}
devc_diff() { # rerun baseline: largest |a-b| / (column max |value|) over the numeric columns of two devc CSV files (2 header rows), or why they cannot be compared
  python3 - "$1" "$2" <<'PYEOF'
import sys
def rd(f):
    L = open(f).read().splitlines()[2:]
    return [[float(x) for x in l.split(',')] for l in L if l.strip()]
try:
    a, b = rd(sys.argv[1]), rd(sys.argv[2])
    if len(a) != len(b): print("row count %d vs %d" % (len(a), len(b))); sys.exit()
    worst = 0.0
    for c in range(len(a[0])):
        sc = max(abs(r[c]) for r in b) or 1.0
        worst = max(worst, max(abs(x[c] - y[c]) for x, y in zip(a, b)) / sc)
    print("max relative difference %.3e over %d rows x %d columns" % (worst, len(a), len(a[0])))
except Exception as e:
    print("cannot compare (%s)" % e)
PYEOF
}
refrun() { # name
  local d=$W/ref/$1; rm -rf "$d"; mkdir -p "$d"; cp "$BL/$1/csmag_32.fds" "$BL/$1/cbc32_uvw.csv" "$d/"
  (cd "$d" && mpirun --bind-to none -np 1 "$REF" csmag_32.fds > stdout.txt 2> stderr.txt) || { echo "REFERENCE RUN FAILED $1"; rc=1; }
  if cmp -s "$d/csmag_32_devc.csv" "$BL/$1/csmag_32_devc.csv"; then echo "reference rerun devc byte-identical to the baseline: $1"
  elif [ "$REF_COMPILER" = gfortran ]; then echo "reference rerun differs from the baseline devc: $1"; rc=1
  else echo "INFO reference rerun devc ($REF_COMPILER) differs from the gfortran baseline devc: $1; $(devc_diff "$d/csmag_32_devc.csv" "$BL/$1/csmag_32_devc.csv") (informational, not a gate)"; fi
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
