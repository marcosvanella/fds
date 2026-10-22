#!/bin/bash
# run_rad_init_csiro_check.sh <fds_amr binary> [<run-dir>]: CSIRO-like t=0 radiation (cases/rad_init_csiro.fds: TMPA 32 C with the &WIND L profile, open sides and top),
# 2 levels with level 1 over the whole domain vs 1 level at the level-1 cell size. FDSTL_T0RAD prints QR and D of the finest level at t=0.
# Variant bed: reflecting DIRT ground and a static grass bed in 0..0.25 m on level 1 (--bed). Stock runs the initial radiation iterations (main.f90:597-622) before it
#   inserts the particles (:634). Control FDSTL_BED_BEFORE_RAD_INIT: the bed inserted before rad_init_levels.
# Variant fixtmp: fixed-temperature ground (TMP_FRONT=32), no bed, NOISE off. Its TMP_F is TMP_0 at the face's own first-cell height (wall.f90 SPECIFIED_TEMPERATURE), so the
#   fine wall must not take the level-0 parent's (fds_fine_level.f90 FRAD_WALLS). Control FDSTL_FINE_TMPF_PULL: the parent's TMP_F, EMISSIVITY, Q_RAD_OUT.
# Pass per variant: the 2-level QR line equals the 1-level one digit for digit, and the control differs. D is printed for information (level-1 t=0 D: Integration).
set -u
BIN=$1; RUN=${2:-$PWD/rad_init_csiro_run}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" vegetation_model.txt
fail=0
q() { grep -a "^T0RAD L[0-9] $2 " "$1" | cut -d' ' -f3-; }
variant() {   # name control-env bed-args sed-expr
  local v=$1 ctl=$2 bed=$3 L
  sed "$4" "$HERE/cases/rad_init_csiro.fds" > ${v}2.fds; sed 's/IJK=8,8,4/IJK=16,16,8/' ${v}2.fds > ${v}1.fds
  run() { env "$@" FDSTL_T0RAD=1 OMP_NUM_THREADS=1 timeout 600 "$BIN" $v$L.fds --two-level-run --levels $L --patch 0 7 0 3 --maxsize 16 --nofine -1 $bed --steps 1 --outdir . --chid $v$L 2>&1; }
  L=1; run X=1 > ${v}1.log || { echo "FAIL rad_init_csiro $v: 1-level run rc (see $RUN/${v}1.log)"; fail=1; return; }
  L=2; run X=1 > ${v}2.log || { echo "FAIL rad_init_csiro $v: 2-level run rc (see $RUN/${v}2.log)"; fail=1; return; }
  L=2; run $ctl=1 > ${v}2n.log || { echo "FAIL rad_init_csiro $v: control run rc (see $RUN/${v}2n.log)"; fail=1; return; }
  echo "[$v] 1L QR: $(q ${v}1.log QR)"; echo "[$v] 2L QR: $(q ${v}2.log QR)"; echo "[$v] control $ctl: $(q ${v}2n.log QR)"
  echo "[$v] 1L D:  $(q ${v}1.log D)"; echo "[$v] 2L D:  $(q ${v}2.log D) (info)"
  [ -n "$(q ${v}1.log QR)" ] && [ "$(q ${v}1.log QR)" = "$(q ${v}2.log QR)" ] || { echo "  [$v] QR at t=0: 2-level differs from 1-level"; fail=1; }
  [ "$(q ${v}2n.log QR)" != "$(q ${v}1.log QR)" ] || { echo "  [$v] control: QR equal (test does not catch)"; fail=1; }
}
variant bed FDSTL_BED_BEFORE_RAD_INIT "--bed 0.5 3.5 0.5 3.5 0 0.25" "s/^//"
grep -aq "particles at t=0: [1-9]" bed2.log || { echo "  [bed] no bed particles on level 1"; fail=1; }
variant fixtmp FDSTL_FINE_TMPF_PULL "" "/&INIT PART/d;s/MATL_ID='DIRT', THICKNESS=0.2 \//TMP_FRONT=32. \//;s/TMPA=32\./TMPA=32., NOISE=.FALSE./"
grep -aq "TMP_FRONT=32" fixtmp2.fds || { echo "  [fixtmp] variant edit failed"; fail=1; }
rm -f ./*.s3d ./*.smv ./*.sf ./*.bf ./*.prt5
[ $fail -eq 0 ] && echo "PASS rad_init_csiro" || echo "FAIL rad_init_csiro (see $RUN)"
exit $fail
