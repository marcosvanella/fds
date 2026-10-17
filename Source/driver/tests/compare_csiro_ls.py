#!/usr/bin/env python3
"""Compare CSIRO LS runs: usage cmp_csiro.py CHID x0,x1,... T1 [run dirs: stock l0 st rg]
Prints 40 s window means and totals of HRR and MLR_FUEL VAPOR, and the level-set front
position (phi=0 crossing along y=0) from the phi<x> point devices."""
import sys, numpy as np
chid, xs, tend = sys.argv[1], np.array([float(v) for v in sys.argv[2].split(',')]), float(sys.argv[3])
runs = sys.argv[4:] or ['stock', 'l0', 'st', 'rg']
def rd(p):
    f = open(p); f.readline(); h = f.readline().strip().split(','); return h, np.loadtxt(f, delimiter=',')
H = {r: rd(f'{r}/{chid}_hrr.csv') for r in runs}; V = {r: rd(f'{r}/{chid}_devc.csv') for r in runs}
h = H[runs[0]][0]
for nm, i in (('HRR kW', h.index('HRR')), ('MLR_FUEL VAPOR kg/s', [k for k, x in enumerate(h) if 'FUEL' in x][0])):
    print(nm)
    w = [(a, min(a + 40, tend)) for a in np.arange(0, tend, 40)]
    for lo, hi in w:
        row = [H[r][1][(H[r][1][:, 0] >= lo) & (H[r][1][:, 0] <= hi), i].mean() for r in runs]
        print(f'  {lo:.0f}-{hi:.0f}s ' + ' '.join(f'{r}={v:.4g}' for r, v in zip(runs, row)) + '  ratio to ' + runs[1] + ': ' + ' '.join(f'{v/row[1]:.3f}' for v in row))
    print('  total ' + ' '.join(f'{r}={np.trapezoid(H[r][1][:, i], H[r][1][:, 0]):.4g}' for r in runs))
def pos(p):
    for i in range(len(p) - 1):
        if p[i] > 0 >= p[i + 1]: return xs[i] + (xs[i + 1] - xs[i]) * p[i] / (p[i] - p[i + 1])
    return np.nan
print('front x (phi=0)')
for t in np.arange(tend / 4, tend + 1, tend / 4):
    print(f'  t={t:.0f} ' + ' '.join(f'{r}={pos(V[r][1][np.argmin(abs(V[r][1][:,0]-t)), 3:3+len(xs)]):.2f}' for r in runs))
