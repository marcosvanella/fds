#!/usr/bin/env python3
"""Leaf energy balance from RADLEAF lines (FDS_RAD_LEAF=1) of the last radiation time: leaf gas loss, leaf wall absorbed (solid SNA), open boundary net (SNO), residual.
usage: bal2.py DIR...   (DIR/run.log)   kW"""
import sys,collections
for d in sys.argv[1:]:
    L=[l.split() for l in open(d+'/run.log',errors='replace') if l.startswith('RADLEAF')]
    if not L: print(d,'no RADLEAF'); continue
    last=max(float(f[1]) for f in L); L=[f for f in L if float(f[1])==last]
    q=collections.defaultdict(float); w=collections.defaultdict(float); o=0.
    for f in L: q[int(f[3])]+=float(f[4]); w[int(f[3])]+=float(f[13]); o+=float(f[15]) if len(f)>15 else 0
    Q=-sum(q.values()); W=sum(w.values())
    print(f"{d:18s} t={last:.3e} gas_loss={Q:9.4f} walls={W:9.4f} open={o:8.4f} resid(loss-walls-open)={Q-W-o:8.4f}  per-level loss {[round(-q[k],4) for k in sorted(q)]} walls {[round(w[k],4) for k in sorted(w)]}")
