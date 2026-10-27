#!/bin/bash
# run_drag_total_check.sh <fds_amr> [work-dir] [ranks=1|2]
# W10 test (PV-4 drag clip and renormalization at level and box edges). Case cases/drag_total.fds: grass bed x 0..25 m, uniform U0 = 2 m/s, free-slip ground, no fire,
# one step. Test run: level 1 over x -10..10 m in 8-cell boxes (the bed crosses box seams and the coarse/fine face at x = 10 m). Reference: one level-1 box over the
# whole bed (no seam, no coarse/fine face at the bed). The gas at the grass is uniform in the first step, so the drag the particles deposit on the gas (FDSTL_DRAG_CHECK,
# summed over all MOVE_IN_GAS calls) must match: |Fx/Fx_ref - 1| <= 1e-3 and |Fy|, |Fz| <= 1e-3 |Fx_ref| (the one-step flow differs at the 1e-4 level), plus test 5 PASS
# in both runs. A face weight dropped without renormalization (or a particle deposited twice) moves Fx by about 1e-2 or more.
set -u
BIN=$1; W=${2:-$PWD/drag_total}; NP=${3:-1}; HERE=$(cd "$(dirname "$0")" && pwd)
run() {   # run <dir> <bed x-high> <maxsize>
  mkdir -p "$W/$1" && sed "s/MPI_PROCESS=1/MPI_PROCESS=$((NP > 1 ? 1 : 0))/" "$HERE/cases/drag_total.fds" > "$W/$1/s.fds" && cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" "$W/$1/vegetation_model.txt"
  (cd "$W/$1" && FDSTL_DRAG_CHECK=1 OMP_NUM_THREADS=1 timeout 600 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" s.fds --two-level-run --bed -10 "$2" -20 20 0 10 --maxsize "$3" --steps 1 > log 2>&1)
}
rm -rf "$W"; run test 10 8; rt=$?; run ref 30 64; rr=$?
python3 - "$W/test/log" "$W/ref/log" "$rt" "$rr" "$NP" <<'P' || { echo "  logs in $W"; exit 1; }
import sys, re
def get(f):
    s = open(f, errors='replace').read()
    t = re.findall(r"RESULT drag total \(deposited force, summed over calls\): (\S+) (\S+) (\S+)", s)
    c = re.findall(r"RESULT drag check 5: (\S+) particle drags.*?(PASS|FAIL)\n", s)
    return ([float(x) for x in t[-1]] if t else None), (c[-1] if c else None)
(ft, ct), (fr, cr) = get(sys.argv[1]), get(sys.argv[2])
np = sys.argv[5]
if sys.argv[3] != '0' or sys.argv[4] != '0' or not ft or not fr or not ct or not cr:
    print(f"FAIL drag_total np{np}: run failed or no RESULT lines (rc test {sys.argv[3]}, ref {sys.argv[4]})"); sys.exit(1)
ex = abs(ft[0] / fr[0] - 1); ey = abs(ft[1]) / abs(fr[0]); ez = abs(ft[2]) / abs(fr[0])
ok = ex <= 1e-3 and ey <= 1e-3 and ez <= 1e-3 and ct[1] == 'PASS' and cr[1] == 'PASS' and ct[0] == cr[0]
msg = f"drag_total np{np}: Fx {ft[0]:.10g} vs ref {fr[0]:.10g} (rel {ex:.2e}), |Fy|/|Fx| {ey:.1e}, |Fz|/|Fx| {ez:.1e} (each <= 1e-3); drags {ct[0]}/{cr[0]}, test 5 {ct[1]}/{cr[1]}"
print(("PASS " if ok else "FAIL ") + msg); sys.exit(0 if ok else 1)
P
rc=$?; [ $rc = 0 ] && rm -rf "$W"; exit $rc
