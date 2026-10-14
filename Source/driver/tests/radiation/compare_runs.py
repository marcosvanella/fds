#!/usr/bin/env python3
"""Compare FDS / fds_amr runs of one case: <CHID>_hrr.csv and <CHID>_devc.csv per run directory.
Means over t >= T0 and the max relative difference (aligned by linear interpolation onto the reference times,
divided by max|reference|) against reference runs. Optional plot of HRR(t) and Q_RADI(t).
usage: compare_runs.py --chid simple_test --ref fds_f1 --ref2 fds_c --t0 1 --png out.png name=dir [name=dir ...]"""
import argparse, csv, sys
import numpy as np

def read(path):
    r = list(csv.reader(open(path)))
    h = [x.strip() for x in r[1]]
    a = np.array([[float(x) for x in row] for row in r[2:] if row])
    return h, a

def load(d, chid):
    out = {}
    for k in ('hrr', 'devc'):
        try:
            h, a = read(f'{d}/{chid}_{k}.csv')
            for i, c in enumerate(h[1:], 1): out[c] = (a[:, 0], a[:, i])
        except (FileNotFoundError, IndexError):
            pass
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--chid', default='simple_test'); p.add_argument('--ref', required=True); p.add_argument('--ref2')
    p.add_argument('--t0', type=float, default=1.0); p.add_argument('--png'); p.add_argument('runs', nargs='+')
    o = p.parse_args()
    runs = dict(x.split('=', 1) for x in o.runs)
    D = {n: load(d, o.chid) for n, d in runs.items()}
    cols = [c for c in D[o.ref] if all(c in D[n] for n in D)]
    names = list(runs)
    print('| quantity | ' + ' | '.join(names) + ' |'); print('|---' * (len(names) + 1) + '|')
    for c in cols:
        row = []
        for n in names:
            t, v = D[n][c]; m = v[t >= o.t0]
            row.append(f'{m.mean():.4g}' if len(m) else 'n/a')
        print(f'| {c} mean t>={o.t0:g} | ' + ' | '.join(row) + ' |')
    for ref in [o.ref] + ([o.ref2] if o.ref2 else []):
        print(f'\nmax |run - {ref}| / max|{ref}| over the common time range (t >= {o.t0:g})')
        print('| quantity | ' + ' | '.join(n for n in names if n != ref) + ' |'); print('|---' * len(names) + '|')
        for c in cols:
            tr, vr = D[ref][c]; sc = np.max(np.abs(vr[tr >= o.t0])) or 1.0
            row = []
            for n in names:
                if n == ref: continue
                t, v = D[n][c]; k = (tr >= max(o.t0, t[0])) & (tr <= t[-1])
                row.append(f'{np.max(np.abs(np.interp(tr[k], t, v) - vr[k])) / sc:.3f}' if k.any() else 'n/a')
            print(f'| {c} | ' + ' | '.join(row) + ' |')
    if o.png:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        for n in names:
            for a, c in zip(ax, ('HRR', 'Q_RADI')):
                if c in D[n]: a.plot(*D[n][c], label=n, lw=1)
        ax[0].set_ylabel('HRR (kW)'); ax[1].set_ylabel('Q_RADI (kW)'); ax[1].set_xlabel('t (s)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3); ax[1].grid(alpha=.3)
        fig.tight_layout(); fig.savefig(o.png, dpi=120); print('\nplot', o.png)

if __name__ == '__main__':
    sys.exit(main())
