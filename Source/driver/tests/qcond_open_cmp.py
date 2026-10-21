#!/usr/bin/env python3
# Q_COND OPEN-face comparison (run_qcond_open_check.sh): QOPEN lines (t, kW) of two logs, integrated (trapezoid) over the common time window;
# PASS if the relative difference is <= TOL or the absolute one <= FLOOR kJ (open-face Q_COND without the defect is ~0.1 kJ, below the 2%
# resolution of the run-to-run noise); also prints the max abs difference per output (b interpolated to the times of a).
import sys, numpy as np
def rd(f): return np.array([[float(x) for x in l.split()[1:3]] for l in open(f, errors='replace') if l.startswith('QOPEN ')])
a, b, tol, fl, tag = rd(sys.argv[1]), rd(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
t1 = min(a[-1, 0], b[-1, 0]); a, b = a[a[:, 0] <= t1], b[b[:, 0] <= t1]
ia, ib = np.trapezoid(a[:, 1], a[:, 0]), np.trapezoid(b[:, 1], b[:, 0])
ok = abs(ia - ib) <= max(tol * abs(ib), fl); rel = abs(ia - ib) / max(abs(ib), 1e-30); mx = np.max(np.abs(a[:, 1] - np.interp(a[:, 0], b[:, 0], b[:, 1])))
print(f"{'PASS' if ok else 'FAIL'} {tag}: t 0..{t1:.4g} s, int {ia:.6g} vs {ib:.6g} kJ, rel {rel:.3g} (tol {tol}, floor {fl} kJ), max abs per output {mx:.4g} kW")
sys.exit(0 if ok else 1)
