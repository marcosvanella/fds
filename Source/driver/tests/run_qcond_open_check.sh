#!/bin/bash
# Q_COND OPEN-face check (FDSTL_QOPEN per-step QOPEN line, kW): leaf Q_COND summed over the OPEN boundary faces only.
# (a) 2 levels (level 1 at 0.5 m spans y, so it reaches the YMIN/YMAX OPEN faces) against 1 level at like dx (0.5 m): PASS within 2% integrated over the window
#     or |diff| <= 1e-4 x the 1-level integrated HRR over the same window (about 0.15 kJ for c10 over 30 steps); (b) control: the same pair with YMIN/YMAX closed (default walls), PASS likewise.
# The floor: open-face Q_COND without the defect is about 0.1 kJ over the window, so 2% of it is below the noise between resolutions and a relative rule alone
# is noise; the floor ties the absolute tolerance to the energy scale of the case (1e-4 of the energy released).
# usage: run_qcond_open_check.sh <build-dir (absolute)> <c10 cut-down C064 input at 1 m, IJK=20,20,20, with &AMR line, vegetation_model.txt beside it (absolute)> [work-dir] [steps]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${3:-$PWD/qcopen}; N=${4:-30}; rc=0; rm -rf "$W"; mkdir -p "$W"
cp "$2" "$W/c1.fds"; cp "$(dirname "$2")/vegetation_model.txt" "$W/"; sed 's/IJK=20,20,20/IJK=40,40,40/' "$2" > "$W/c05.fds"; grep -v "MB='Y" "$W/c1.fds" > "$W/y1.fds"; grep -v "MB='Y" "$W/c05.fds" > "$W/y05.fds"
run() {  # tag input levels
  (cd "$W" && FDSTL_QOPEN=1 timeout 1800 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" "$2.fds" --two-level-run --levels "$3" --blocking 4 \
     --outdir . --steps "$N" --chid "$1" > "$1.log" 2>&1; mv c_cat_hrr.csv "$1_hrr.csv")   # input CHID c
}
run l2 c1 2; run l1 c05 1; run y2 y1 2; run y1 y05 1
python3 "$HERE/qcond_open_cmp.py" "$W/l2.log" "$W/l1.log" 0.02 1e-4 "$W/l1_hrr.csv" "2L leaf OPEN vs 1L 0.5 m" || rc=1
python3 "$HERE/qcond_open_cmp.py" "$W/y2.log" "$W/y1.log" 0.02 1e-4 "$W/y1_hrr.csv" "control, YMIN/YMAX closed: 2L vs 1L 0.5 m" || rc=1
exit $rc
