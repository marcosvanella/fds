#!/bin/bash
# run_regrid_burner_check.sh <bin> <abs-path-case.fds> [steps=70]: two-level burner-on-fine run with --tag-obst --tag-temp 400 and regrids every 5 steps
# (no --tag-vort, no --patch: the default TRACER criterion must not apply). Regression of the "flux override rejected" abort: the burner face lists of a fine level
# must follow the box layout of every regrid. Passes when the run reaches the last step without the flux-override abort (the end-of-run mass-balance gate is not part of this check). Run in a scratch dir.
set -u
steps=${3:-70}
rm -rf rb_run; mkdir rb_run; cd rb_run
OMP_NUM_THREADS=1 $1 $2 --two-level-run --levels 2 --steps $steps --tag-obst --nofine -1 --tag-temp 400 --regrid-every 5 --blocking 4 --maxsize 8 --outdir . --chid st --log-every 1 > run.log 2>&1
last=$(grep -c "^TWO-LEVEL step" run.log)
if grep -q "flux override rejected" run.log || [ "$last" -lt "$steps" ]; then echo "FAIL regrid burner: $last of $steps steps"; grep -m1 "Abort" run.log; exit 1; fi
echo "PASS regrid burner: $last steps, $(grep -c 'mid-regrid' run.log) regrids"
