#!/bin/bash
# IR-006 check: a USE_AMREX=OFF build of a (patched) source tree gives bitwise-identical FDS output to the FireX baseline.
# usage: check_off_bitwise.sh <source-tree> <build-dir> <case> [np]     (case: shunn3_32 (np 1), shunn3_4mesh_32 (np 4), ...)
# <source-tree> must be a copy OUTSIDE the reference tree with the patches applied; <build-dir> outside the source tree.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
SRC=$1; BLD=$2; CASE=$3; NP=${4:-1}
[ -x "$BLD/fds" ] || { cmake -S "$SRC" -B "$BLD" $FDS_CMAKE_COMMON > "$BLD.cfg.log" 2>&1 && cmake --build "$BLD" -j4 > "$BLD.build.log" 2>&1; } || { echo "BUILD FAILED"; exit 2; }
RUN=$BLD/run/$CASE; rm -rf "$RUN"; mkdir -p "$RUN"; FDSF=$(basename "$(ls "$BASELINE/$CASE"/*.fds | head -1)")   # derived directories (e.g. csmag_32__fishpak_bc000) keep the CHID of the source input
cp "$BASELINE/$CASE/$FDSF" "$RUN/"
# input csv files of the case (e.g. csmag_32 UVW_FILE data) are copied; the FDS output csv files are the ones compared below
for f in "$BASELINE/$CASE"/*.csv; do case "$f" in *_hrr.csv|*_mass.csv|*_mms.csv|*_devc.csv|*_cpu.csv|*_steps.csv|*_uvw_t*_m*.csv) ;; *) cp "$f" "$RUN/";; esac; done 2>/dev/null
(cd "$RUN" && timeout 600 mpirun --bind-to none -np "$NP" "$BLD/fds" "$FDSF" > stdout.txt 2> stderr.txt) || { echo "RUN FAILED"; exit 2; }
rc=0; n=0
for f in $(cd "$BASELINE/$CASE" && ls *.csv *.restart *.sf 2>/dev/null | grep -v -e _cpu.csv -e _steps.csv); do
  n=$((n+1)); cmp -s "$BASELINE/$CASE/$f" "$RUN/$f" || { echo "DIFFER $f"; rc=1; }
done
[ $rc = 0 ] && echo "PASS $CASE np=$NP: $n output files bitwise identical to baseline" || echo "FAIL $CASE"
exit $rc
