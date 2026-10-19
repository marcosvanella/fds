#!/usr/bin/env python3
"""Radiation validation tables: fds_amr (FDS_RAD_LEAF=1 RADLEAF lines in run.log, leaf-only sums over all boxes) against plain FDS devc.csv
(devices QR_int, W_net, W_inc, W_q, U_reg, Q_reg from Source/driver/tests/radiation/*.fds-like inputs: RADIATION LOSS / RADIATIVE HEAT FLUX /
INCIDENT HEAT FLUX surface and volume integrals) on the coarse and the fine uniform grid.
The leaf rule (finest level covering each location) is applied in the driver diagnostic (RAD_COVERED skips cells and wall faces under the next level; the
solid and OPEN hidden walls of fine boxes are summed at their level), so the AMR numbers here are plain sums. --mode mean is a trapezoid time mean.
usage: val_compare.py --mode last|mean [--t0 T] --fds c=DIR f1=DIR [f2=DIR] --amr 1=DIR 2=DIR 3=DIR [--hrr]"""
import argparse, csv, re
import numpy as np

def fds_series(d, chid='in'):
    r = list(csv.reader(open(f'{d}/{chid}_devc.csv'))); h = [x.strip() for x in r[1]]
    a = np.array([[float(x) for x in row] for row in r[2:] if row])
    out = {c: (a[:, 0], a[:, i]) for i, c in enumerate(h) if i}
    try:
        r = list(csv.reader(open(f'{d}/{chid}_hrr.csv'))); h = [x.strip() for x in r[1]]
        a = np.array([[float(x) for x in row] for row in r[2:] if row])
        for c in ('HRR', 'Q_RADI'): out[c] = (a[:, 0], a[:, h.index(c)])
    except (FileNotFoundError, ValueError, IndexError): pass
    return out

def amr_series(d, chid='in'):
    """groups of RADLEAF lines (one radiation call of every box) -> time, sums"""
    g, cur, prev = [], None, None
    for l in open(f'{d}/run.log', errors='replace'):
        if not l.startswith('RADLEAF'): continue
        f = l.split(); t, nm, lev = float(f[1]), int(f[2]), int(f[3]); v = [float(x) for x in f[4:9]] + [float(f[10]), int(f[9])]; v[1] = float(f[13]) + (float(f[15]) if len(f) > 15 else 0.)   # W_net = SNA (net flux of every leaf solid face incl. OBST) + SNO (open boundary), as the FDS RADIATIVE HEAT FLUX surface integral
        if cur is None or (lev, nm) <= prev: cur = [t, np.zeros(6), 0, []]; g.append(cur)
        prev = (lev, nm); cur[1] += np.array(v[:5] + [v[5]]); cur[3].append((nm, lev, v[0], v[4])); cur[2] += 1
    t = np.array([x[0] for x in g]); S = np.array([x[1] for x in g])
    out = {c: (t, S[:, i]) for i, c in enumerate(['QR_int', 'W_net', 'W_inc', 'W_q', 'U_reg', 'Q_reg'])}
    try:
        r = list(csv.reader(open(f'{d}/{chid}_hrr.csv'))); h = [x.strip() for x in r[1]]
        a = np.array([[float(x) for x in row] for row in r[2:] if row])
        for c in ('HRR', 'Q_RADI'): out[c] = (a[:, 0], a[:, h.index(c)])
    except (FileNotFoundError, ValueError, IndexError): pass
    return out, g

def value(s, mode, t0):
    t, v = s
    if mode == 'last': return v[-1]
    k = t >= t0; t, v = t[k], v[k]
    if len(t) < 2: return v[-1]
    return np.trapezoid(v, t) / (t[-1] - t[0])

def main():
    p = argparse.ArgumentParser(); p.add_argument('--mode', default='last'); p.add_argument('--t0', type=float, default=1.0)
    p.add_argument('--fds', nargs='+', required=True); p.add_argument('--amr', nargs='+', required=True)
    o = p.parse_args()
    F = {k: fds_series(v) for k, v in (x.split('=') for x in o.fds)}
    A = {k: amr_series(v)[0] for k, v in (x.split('=') for x in o.amr)}
    cols = [c for c in ('HRR', 'Q_RADI', 'QR_int', 'W_net', 'W_inc', 'W_q', 'U_reg', 'Q_reg') if all(c in s for s in list(F.values()) + list(A.values()))]
    names = [f'FDS {k}' for k in F] + [f'AMR {k}-lev' for k in A]
    print(f'| quantity ({o.mode}) | ' + ' | '.join(names) + ' |'); print('|---' * (len(names) + 1) + '|')
    V = {}
    for c in cols:
        row = [value(s[c], o.mode, o.t0) for s in list(F.values()) + list(A.values())]; V[c] = row
        print(f'| {c} | ' + ' | '.join(f'{x:.5g}' for x in row) + ' |')
    nf = len(F); keys = list(F)
    for ref in [x for x in ('f2', 'f1', 'c') if x in F]:
        i = keys.index(ref)
        print(f'\nrelative error vs FDS {ref}: |x - ref| / |ref|'); print('| quantity | ' + ' | '.join(n for j, n in enumerate(names) if j != i) + ' |'); print('|---' * len(names) + '|')
        for c in cols:
            r = V[c][i]; print(f'| {c} | ' + ' | '.join(f'{abs(x - r) / (abs(r) or 1):.3f}' for j, x in enumerate(V[c]) if j != i) + ' |')

if __name__ == '__main__': main()
