#!/bin/bash
# D-053: pressure-zone sums USUM/DSUM/PSUM. Default = FDS summation order; --exact-zone-sums = the D-028 fixed-point sums (GPU option, off by default). Separate test of the option.
# usage: run_zone_sum_check.sh <driver-build> [work-dir]
# The first step is run (--steps 1, FDSTL_ZONES=1 prints the zone sums of the predictor and corrector pass) on dec1 (1 box, 1 rank), dec4 (16 boxes, 1 rank), dec2_np4 and dec4_np4 (4 ranks):
#  (1) exact mode: the printed predictor-pass sums (first line: identical input state, reduction-free stages are bitwise) of every layout are bitwise equal to the 1-box ones;
#      later passes start from a pressure solution whose summation order depends on the layout, so only the first line is compared;
#  (2) default mode, 1 box / 1 rank: the applied predictor sums are bitwise the "FDS order" sums that the exact mode prints next to its own (same accumulation, only the write-back differs);
#  (3) the two modes differ by rounding only: |exact - FDS order| <= 1e-9 * max(|exact|, 1) for every sum;
#  The exact mode reproduces the pre-D-053 driver bitwise (checked once at S9: final fields and steps files byte-identical for dec1, dec2, dec4/4 ranks).
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; WORK=${2:-$BLD/run/zonesum}; mkdir -p "$WORK"; rc=0
run() { # label case np [args]
  local d="$WORK/$1"; rm -rf "$d"; mkdir -p "$d"; cp "$HERE/cases/$2.fds" "$d/"
  (cd "$d" && FDSTL_ZONES=1 timeout 600 mpirun --bind-to none --oversubscribe -np $3 "$BLD/fds_amr" "$2.fds" --run --steps 1 --outdir . --chid $2 --quiet "${@:4}" > stdout.txt 2> stderr.txt) || { echo "FAIL run $1"; rc=1; }
}
for c in "dec1 dec1 1" "dec4 dec4 1" "dec2_np4 dec2_np4 4" "dec4_np4 dec4_np4 4"; do set -- $c; run "x_$1" $2 $3 --exact-zone-sums; done
run "d_dec1" dec1 1
run "d_dec4_np4" dec4_np4 4
python3 - "$WORK" <<'P' || rc=1
import re, sys
w = sys.argv[1]; bad = 0
def zl(path):   # lines "[zone] icyc N P|C DSUM a PSUM b USUM c (FDS order: a b c)" or "... FDS order DSUM a PSUM b USUM c"
    ex, fo = [], []
    for l in open(path + '/stdout.txt'):
        if '[zone]' not in l: continue
        m = re.search(r'icyc (\d+) (P|C) DSUM (\S+) PSUM (\S+) USUM (\S+) \(FDS order: (\S+) (\S+) (\S+)\)', l)
        if m: ex.append(m.group(1, 2, 3, 4, 5)); fo.append(m.group(1, 2, 6, 7, 8)); continue
        m = re.search(r'icyc (\d+) (P|C) FDS order DSUM (\S+) PSUM (\S+) USUM (\S+)', l)
        if m: fo.append(m.group(1, 2, 3, 4, 5))
    return ex, fo
ref, ref_fo = zl(w + '/x_dec1')
if not ref: print('FAIL: no zone lines in the exact-mode 1-box run'); sys.exit(1)
for r in ['dec4', 'dec2_np4', 'dec4_np4']:
    ex, _ = zl(w + '/x_' + r)
    same = ex[:1] == ref[:1]
    print('exact mode %-9s predictor-pass zone sums bitwise equal to the 1-box run: %s (%d lines)' % (r, 'yes' if same else 'NO', len(ex)))
    bad += 0 if same else 1
_, d1 = zl(w + '/d_dec1')
same = d1[:1] == ref_fo[:1]
print('default mode dec1: applied predictor sums == FDS-order sums printed by the exact mode: %s' % ('yes' if same else 'NO'))
bad += 0 if same else 1
worst = 0.0
for a, b in zip(ref, ref_fo):
    for q in range(2, 5):
        e, f = float(a[q]), float(b[q]); worst = max(worst, abs(e - f) / max(abs(e), 1.0))
print('exact vs FDS-order sums (1 box): max |diff| / max(|exact|, 1) = %.2e %s' % (worst, 'ok' if worst <= 1e-9 else 'FAIL'))
bad += 0 if worst <= 1e-9 else 1
_, d4 = zl(w + '/d_dec4_np4')
print('information: default mode dec4_np4 predictor sums %s the 1-box FDS-order sums (rounding differences between layouts are the reason for the exact option)' % ('equal' if d4[:1] == d1[:1] else 'differ from'))
sys.exit(1 if bad else 0)
P
[ $rc = 0 ] && echo "ZONE SUM CHECK PASS" || echo "ZONE SUM CHECK FAIL"
exit $rc
