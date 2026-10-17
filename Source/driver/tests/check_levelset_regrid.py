#!/usr/bin/env python3
"""check_levelset_regrid.py <run-dir>: invariants of the regrid transfer on the plotfiles of a --tag-phi level-set run (plot-every = regrid-every, so every plot
is written right after a regrid). Per plot: PHI is constant in k on both levels; the fine PHI is the injected level-0 PHI of its parent cell (fill_phi at the regrid);
Over the series: fine cells were created and the layout moved. Informational: the largest relative gap between the level-0 RHO under a fine box and the average of its
8 children at the end of the step (not asserted: covered level-0 scalars are only averaged down at regrids, the gap is up to a few percent next to the box edges)."""
import sys, re, glob, numpy as np

def fabs(level_dir, ncomp):
    out = []
    lines = open(level_dir + '/Cell_H').read().split('\n')
    boxes = [tuple(map(int, m.groups())) for m in (re.match(r'\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)', l) for l in lines) if m]
    offs = [(l.split()[1], int(l.split()[2])) for l in lines if l.startswith('FabOnDisk')]
    for b, (fn, off) in zip(boxes, offs):
        raw = open(level_dir + '/' + fn, 'rb')
        raw.seek(off); h = raw.readline()
        dims = [b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1]
        n = int(np.prod(dims)) * ncomp
        a = np.frombuffer(raw.read(8 * n), '<f8').reshape(ncomp, dims[2], dims[1], dims[0])
        out.append((b, a))
    return out

def names(plot):
    L = open(plot + '/Header').read().split('\n'); n = int(L[1]); return L[2:2 + n]

def main(run):
    plots = sorted(glob.glob(run + '/plt/plt_*_s*'), key=lambda d: int(re.search(r'_s(\d+)$', d).group(1)))
    bad = 0; gap = 0.0; checked = 0; created = moved = removed = 0; prev = None
    for d in plots:
        nm = names(d); ip = nm.index('PHI'); ir = nm.index('RHO'); it = nm.index('TMP'); nc = len(nm)
        f0 = fabs(d + '/Level_0', nc); f1 = fabs(d + '/Level_1', nc)
        if not f1: continue
        b0, a0 = f0[0]
        for lev, fl in ((0, f0), (1, f1)):
            for b, a in fl:
                if np.ptp(a[ip], axis=0).max() > 1e-12 * max(1.0, np.abs(a[ip]).max()): print('PHI not constant in k, level', lev, b, 'in', d); bad += 1
        cov = {}
        for b, a in f1:
            for k in range(a.shape[1]):
                for j in range(a.shape[2]):
                    for i in range(a.shape[3]):
                        pi, pj, pk = (b[0] + i) // 2 - b0[0], (b[1] + j) // 2 - b0[1], (b[2] + k) // 2 - b0[2]
                        if abs(a[ip, k, j, i] - a0[ip, pk, pj, pi]) > 1e-12 * max(1.0, abs(a0[ip, pk, pj, pi])): bad += 1; print('fine PHI differs from parent at', (b[0] + i, b[1] + j, b[2] + k), 'in', d); break
                        cov.setdefault((pi, pj, pk), []).append((a[ir, k, j, i], a[it, k, j, i]))
        for (pi, pj, pk), v in cov.items():
            if len(v) == 8: gap = max(gap, abs(a0[ir, pk, pj, pi] - np.mean([x[0] for x in v])) / a0[ir, pk, pj, pi])
        cells = set(cov)
        if prev is not None:
            created += len(cells - prev); removed += len(prev - cells); moved += int(cells != prev)
        prev = cells; checked += 1
    print('layout events over %d plots: %d cells created, %d removed, %d plot-to-plot layout changes' % (checked, created, removed, moved))
    print('largest level-0 RHO vs children-average gap at the end of a step: %.2e (informational)' % gap)
    if checked < 5 or created == 0 or moved < 3: print('too few layout events to test'); bad += 1
    print(('FAIL' if bad else 'PASS') + ' levelset_regrid plots: %d checked, %d violation(s)' % (checked, bad))
    return 1 if bad else 0

sys.exit(main(sys.argv[1]))
