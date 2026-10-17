#!/usr/bin/env python3
"""Compare stock FDS and fds_amr (unrefined, and refined if present) output of one level-set case: final HRR and the last row of every devc column
(first column t excluded), relative to stock. usage: wui_levelset_compare.py CASEDIR CHID"""
import csv, os, sys
d, c = sys.argv[1], sys.argv[2]; tol = float(os.environ.get('TOL', '1e-6'))
def last(path):
    if not os.path.exists(path): return None
    r = list(csv.reader(open(path)))
    if len(r) < 3: return None
    return dict(zip([x.strip() for x in r[1]], [float(x) for x in r[-1]]))
def load(k):
    w = f'{d}/{k}'
    if not os.path.isdir(w): return None
    rc = open(f'{w}/rc').read().strip() if os.path.exists(f'{w}/rc') else '?'
    return rc, last(f'{w}/{c}_hrr.csv'), last(f'{w}/{c}_devc.csv')
S = load('stock'); ok = True
print(f'== {c}  stock rc={S[0]}')
for k in ('amr', 'ref'):
    X = load(k)
    if X is None: continue
    if X[2] is None or S[2] is None: print(f'  {k}: rc={X[0]} no devc output'); ok = ok and k == 'ref'; continue
    worst = 0.; line = []
    for n, v in S[2].items():
        if n == 'Time' or n not in X[2]: continue
        e = abs(X[2][n] - v) / max(abs(v), 1e-12) if abs(v) > 1e-12 else abs(X[2][n] - v)
        worst = max(worst, e)
        if e > 0.01: line.append(f'{n}: {v:.5g} -> {X[2][n]:.5g}')
    h = ''
    if S[1] and X[1]: h = f" HRR stock {S[1].get('HRR', float('nan')):.5g} {k} {X[1].get('HRR', float('nan')):.5g}"
    print(f'  {k}: rc={X[0]} t_end stock {S[2]["Time"]:.4g} {k} {X[2]["Time"]:.4g} devices max rel diff {worst:.3g}{h}' + ('  >1%: ' + '; '.join(line[:6]) if line else ''))
    if k == 'amr' and worst > tol: ok = False
sys.exit(0 if ok else 1)
