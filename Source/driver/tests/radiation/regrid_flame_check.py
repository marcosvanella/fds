#!/usr/bin/env python3
"""Radiation continuity across regrids: reads run.log of run_regrid_flame.sh (RADLEAF lines, FDS_RAD_LEAF=1 FDS_RAD_REG=whole domain; 'mid-regrid ... changed 1' lines).
Per regrid: leaf QR*V (SQ) and UII*V (PU) of the last radiation group before and the first after (relative change, against the median step-to-step change of the run),
and the fine boxes' own SQ/PU in the first and the second group after the regrid (a box that radiates at once has first ~ second; one started at ambient would not)."""
import sys
g, cur, prev, ng, marks = [], None, None, 0, []
for l in open(sys.argv[1], errors='replace'):
    if ('mid-regrid' in l and ' changed 1 ' in l) or ('regrid_dynamic' in l and 'changed 1' in l): marks.append(len(g))
    if not l.startswith('RADLEAF'): continue
    f = l.split(); t, nm, lev = float(f[1]), int(f[2]), int(f[3]); sq, pu = float(f[4]), float(f[8])
    if cur is None or (lev, nm) <= prev: cur = {'t': t, 'SQ': 0.0, 'PU': 0.0, 'FSQ': 0.0, 'FPU': 0.0, 'nb': 0}; g.append(cur)
    prev = (lev, nm); cur['SQ'] += sq; cur['PU'] += pu; cur['nb'] += 1
    if lev > 0: cur['FSQ'] += sq; cur['FPU'] += pu
rel = lambda a, b: abs(a - b) / (abs(b) or 1.0)
step = sorted(rel(g[i]['SQ'], g[i-1]['SQ']) for i in range(1, len(g)))
stepu = sorted(rel(g[i]['PU'], g[i-1]['PU']) for i in range(1, len(g)))
print(f'{len(g)} radiation groups, {len(marks)} regrids; median group-to-group change SQ {step[len(step)//2]:.2e}, PU {stepu[len(stepu)//2]:.2e}')
print('| regrid | t | boxes before->after | SQ before | SQ after | rel | PU rel | fine SQ 1st / 2nd after | fine PU 1st / 2nd after |'); print('|---' * 9 + '|')
for k, m in enumerate(marks):
    if m < 1 or m + 1 >= len(g): continue
    b, a1, a2 = g[m-1], g[m], g[m+1]
    print(f"| {k} | {a1['t']:.3f} | {b['nb']}->{a1['nb']} | {b['SQ']:.4g} | {a1['SQ']:.4g} | {rel(a1['SQ'], b['SQ']):.2e} | {rel(a1['PU'], b['PU']):.2e} | {a1['FSQ']:.4g} / {a2['FSQ']:.4g} | {a1['FPU']:.4g} / {a2['FPU']:.4g} |")
import statistics as st
R = [(rel(g[m]['SQ'], g[m-1]['SQ']), rel(g[m]['PU'], g[m-1]['PU']), g[m+1]['FPU'] / g[m]['FPU'] if g[m]['FPU'] > 0 and g[m]['nb'] == g[m+1]['nb'] else None) for m in marks if 1 <= m < len(g) - 1 and g[m]['FPU'] > 0]
if R:
    ramp = [x[2] for x in R if x[2]]
    print(f"summary over {len(R)} regrids with fine boxes: median change across a regrid SQ {st.median(x[0] for x in R):.3f}, PU {st.median(x[1] for x in R):.3f} (median over all steps SQ {step[len(step)//2]:.3f}, PU {stepu[len(stepu)//2]:.3f}); "
          f"median fine PU 2nd/1st update after a regrid (same boxes) {st.median(ramp):.3f} (ambient start would give a large ramp)" if ramp else '')
