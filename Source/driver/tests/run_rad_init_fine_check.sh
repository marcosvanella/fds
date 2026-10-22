#!/bin/bash
# run_rad_init_fine_check.sh <fds_amr binary> [<run-dir>]: initial radiation of a fine level (main.f90 ~599 INITIAL_RADIATION_ITERATIONS, TimeLoop rad_init_levels).
# cases/rad_init_fine.fds (hot &INIT block aligned to level-0 cells, black walls, open top) as 2 levels with level 1 over the whole domain, and as 1 level at the level-1 cell
# size (IJK doubled). FDSTL_T0RAD prints sum, norm1, min, max of QR and D of the finest level at t=0, after the initial projection.
# Pass: the 2-level QR line equals the 1-level one digit for digit, and the negative control (FDSTL_NO_RAD_INIT, no initial iterations on level 1) differs.
# Black walls: a fine wall takes its radiative state from the level-0 wall (fds_fine_walls_rad_xchg), whose reflected part comes from the coarse solve; with reflecting
# walls the 2-level QR is ~5e-6 off (relative), see NOTE_fine_init_rad.md. D is printed for information: the level-1 D at t=0 is the level-0 D prolonged (not radiation).
set -u
BIN=$1; RUN=${2:-$PWD/rad_init_fine_run}; SRC=$(cd "$(dirname "$0")" && pwd)/cases/rad_init_fine.fds
mkdir -p "$RUN" && cd "$RUN" || exit 1
sed 's/IJK=8,8,4/IJK=16,16,8/' "$SRC" > r1.fds; cp "$SRC" r2.fds
run() { env "$@" FDSTL_T0RAD=1 OMP_NUM_THREADS=1 timeout 600 "$BIN" r$L.fds --two-level-run --levels $L --patch 0 7 0 3 --maxsize 16 --nofine -1 --steps 1 --outdir . --chid r$L 2>&1; }
L=1; run X=1 > r1.log || { echo "FAIL rad_init_fine: 1-level run rc (see $RUN/r1.log)"; exit 1; }
L=2; run X=1 > r2.log || { echo "FAIL rad_init_fine: 2-level run rc (see $RUN/r2.log)"; exit 1; }
L=2; run FDSTL_NO_RAD_INIT=1 > r2n.log || { echo "FAIL rad_init_fine: control run rc (see $RUN/r2n.log)"; exit 1; }
rm -f ./*.s3d ./*.smv ./*.sf ./*.bf ./*.prt5
q() { grep -a "^T0RAD L[0-9] $2 " "$1" | cut -d' ' -f3-; }
fail=0
echo "1L QR: $(q r1.log QR)"; echo "2L QR: $(q r2.log QR)"; echo "2L QR without initial iterations (control): $(q r2n.log QR)"
echo "1L D:  $(q r1.log D)"; echo "2L D:  $(q r2.log D) (info)"
[ -n "$(q r1.log QR)" ] && [ "$(q r1.log QR)" = "$(q r2.log QR)" ] || { echo "  QR at t=0: 2-level differs from 1-level"; fail=1; }
[ "$(q r2n.log QR)" != "$(q r1.log QR)" ] || { echo "  control: QR equal without the initial iterations (test does not catch)"; fail=1; }
[ $fail -eq 0 ] && echo "PASS rad_init_fine" || echo "FAIL rad_init_fine (see $RUN)"
exit $fail
