#!/bin/bash
# run_rad_init_csiro_check.sh <fds_amr binary> [<run-dir>]: CSIRO-like t=0 radiation (cases/rad_init_csiro.fds: TMPA 32 C with the &WIND L profile, reflecting DIRT ground,
# open sides and top, static grass bed in 0..0.25 m) as 2 levels with level 1 over the whole domain and the bed on level 1 (--bed), and as 1 level at the level-1 cell size.
# Stock (main.f90) runs the initial radiation iterations (:597-622) before it inserts the particles (:634), so the bed takes no part in the t=0 QR.
# Pass: the 2-level QR line (FDSTL_T0RAD, finest level at t=0) equals the 1-level one digit for digit, and the negative control (FDSTL_BED_BEFORE_RAD_INIT: bed inserted
# before rad_init_levels) differs. D is printed for information (level-1 t=0 D: Integration).
set -u
BIN=$1; RUN=${2:-$PWD/rad_init_csiro_run}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" vegetation_model.txt
sed 's/IJK=8,8,4/IJK=16,16,8/' "$HERE/cases/rad_init_csiro.fds" > r1.fds; cp "$HERE/cases/rad_init_csiro.fds" r2.fds
run() { env "$@" FDSTL_T0RAD=1 OMP_NUM_THREADS=1 timeout 600 "$BIN" r$L.fds --two-level-run --levels $L --patch 0 7 0 3 --maxsize 16 --nofine -1 --bed 0.5 3.5 0.5 3.5 0 0.25 --steps 1 --outdir . --chid r$L 2>&1; }
L=1; run X=1 > r1.log || { echo "FAIL rad_init_csiro: 1-level run rc (see $RUN/r1.log)"; exit 1; }
L=2; run X=1 > r2.log || { echo "FAIL rad_init_csiro: 2-level run rc (see $RUN/r2.log)"; exit 1; }
L=2; run FDSTL_BED_BEFORE_RAD_INIT=1 > r2n.log || { echo "FAIL rad_init_csiro: control run rc (see $RUN/r2n.log)"; exit 1; }
rm -f ./*.s3d ./*.smv ./*.sf ./*.bf ./*.prt5
q() { grep -a "^T0RAD L[0-9] $2 " "$1" | cut -d' ' -f3-; }
fail=0
echo "1L QR: $(q r1.log QR)"; echo "2L QR: $(q r2.log QR)"; echo "2L QR with the bed inserted before the initial radiation (control): $(q r2n.log QR)"
echo "1L D:  $(q r1.log D)"; echo "2L D:  $(q r2.log D) (info)"
grep -aq "particles at t=0: [1-9]" r2.log || { echo "  no bed particles on level 1"; fail=1; }
[ -n "$(q r1.log QR)" ] && [ "$(q r1.log QR)" = "$(q r2.log QR)" ] || { echo "  QR at t=0: 2-level differs from 1-level"; fail=1; }
[ "$(q r2n.log QR)" != "$(q r1.log QR)" ] || { echo "  control: QR equal with the bed inserted first (test does not catch)"; fail=1; }
[ $fail -eq 0 ] && echo "PASS rad_init_csiro" || echo "FAIL rad_init_csiro (see $RUN)"
exit $fail
