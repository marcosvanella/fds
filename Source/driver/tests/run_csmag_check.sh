#!/bin/bash
# csmag_32 (3-D, all six faces periodic) against a reference-dump run of the same input with &PRES FISHPAK_BC=0,0,0 (cases/csmag_32_fishpak.fds).
# Why that reference: with the plain csmag_32.fds FDS treats the periodic vents as mesh-to-self interpolated boundaries and iterates the Dirichlet-coupled pressure
# problem to the velocity tolerance only (README, "csmag_32"), which the level-wide periodic FFT solve does not do; FISHPAK_BC=0,0,0 makes FDS itself solve the periodic problem.
# usage: run_csmag_check.sh <build-dir> <ref-dir> [work-dir]   <ref-dir> holds ref.dump (of cases/csmag_32_fishpak.fds, steps 1 and 2) and cbc32_uvw.csv
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; REF=$2; WORK=${3:-$BLD/run/csmag}; mkdir -p "$WORK"; rc=0
for s in 1 2; do
  d="$WORK/s$s"; rm -rf "$d"; mkdir -p "$d"; cp "$HERE/cases/csmag_32_fishpak.fds" "$d/csmag_32.fds"; cp "$REF"/*.csv "$d/"
  (cd "$d" && FDSTL_STAGE=$s timeout 900 mpirun --bind-to none -np 1 "$BLD/fds_amr" csmag_32.fds --run --outdir . --chid csmag_32 --quiet > stdout.txt 2> stderr.txt) || { echo "FAIL run s$s"; rc=1; }
done
python3 "$HERE/csmag_compare.py" "$WORK" "$REF/ref.dump" || rc=1
exit $rc
