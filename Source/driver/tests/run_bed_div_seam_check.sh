#!/bin/bash
# run_bed_div_seam_check.sh <fds_amr binary> [<run-dir>]: grass bed (&INIT PART, x 0..25 m) on a partial-coverage level 1 (--bed -10 10 -20 20 0 10, 4 ranks): the bed crosses
# the coarse-fine face at x = 10 m and the faces between the level-1 boxes. Corrector particle drag on a face shared by two fine boxes must carry the drag of both boxes
# (stock PARTICLE_MOMENTUM_TRANSFER at an interpolated mesh boundary); with only the box's own share the copies of the face disagree and fine div u - D ~ 1e-2.
# Pass: the run's own gate PASS (worst fine div u - D and predictor div US - DS <= 1e-9) and bed particles present.
set -u
BIN=$1; RUN=${2:-$PWD/bed_div_seam_run}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" vegetation_model.txt
cp "$HERE/cases/bed_div_seam.fds" s.fds
OMP_NUM_THREADS=1 timeout 600 mpirun -np 4 "$BIN" s.fds --two-level-run --bed -10 10 -20 20 0 10 --log-every 1 > log 2>&1; rc=$?
fail=0
grep -a "worst fine-level\|RESULT gate" log | sed 's/, worst |mass.*//'
grep -aq "particles at t=0: [1-9]" log || { echo "  no bed particles"; fail=1; }
grep -aq "RESULT gate PASS" log && [ $rc -eq 0 ] || { echo "  gate FAIL or rc $rc"; fail=1; }
rm -f ./*.s3d ./*.s3d.sz ./*.smv ./*.sf ./*.bf ./*.prt5 ./*.prt5.bnd
[ $fail -eq 0 ] && echo "PASS bed_div_seam" || echo "FAIL bed_div_seam (see $RUN)"
exit $fail
