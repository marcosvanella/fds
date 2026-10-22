#!/usr/bin/env python3
"""lsfine_compare.py REF_DIR RUN_DIR [--mode bitid|crossmesh] [--times 10,20,30] [--chid c]
Tolerance AGREED (Role 1, 2026-10-19); see Source/driver/notes/vv-lsfine-multirank-matrix.md. Exit 1 on any FAIL.
--mode bitid (same mesh count on more ranks, vs the 1-rank run of that mesh count; also the 1/1 determinism repeat):
   Owner ruling 2026-10-20: c_hrr.csv Q_* columns may differ by <= HRR_Q_RANK_TOL (1e-11) relative (round-off); devc, front x, gate strict.
   NEEDS-RULING (exit 3 if nothing else fails): leaf PHI differences, and c_multi_level.csv differences confined to mass / mass_rel_change.
   otherwise strict bit-identity: every <chid>_*.csv byte-identical (timing files _cpu, _driver_perf, _steps excluded), front x identical,
   final-plotfile PHI on leaf cells identical, gate PASS.
--mode crossmesh (any mesh count vs 1 mesh / 1 rank):
   front x identical (DEVC 'x', MAXLOC X of burning rate, at --times), gate verdict identical to the reference,
   mean HRR for t >= HRR_MEAN_T0 (default 10 s) within HRR_MEAN_TOL (default 0.02) relative. CSV and PHI differences are printed as INFO only
   (PHI also flags non-finite / garbage values > 1e3).
Leaf cells only: the csv integrals are the driver's leaf sums; PHI is compared with each finest-grid cell holding its leaf cell value.
"""
import sys, os, glob, re, argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('ref'); ap.add_argument('run')
ap.add_argument('--times', default='10,20,30'); ap.add_argument('--chid', default='c')
ap.add_argument('--mode', default='bitid', choices=['bitid', 'crossmesh'])
a = ap.parse_args()
TOL_CSV = float(os.environ.get('TOL_CSV', '1e-8')); TOL_FRONT = float(os.environ.get('TOL_FRONT', '0'))
BITID = a.mode == 'bitid' or os.environ.get('BITID', '0') == '1'
CROSS = a.mode == 'crossmesh'
HRR_Q_RANK_TOL = float(os.environ.get('HRR_Q_RANK_TOL', '2e-11'))   # ruling (a): Q_* sums of _hrr.csv 2e-11 relative across rank counts   # Owner ruling 2026-10-20: c_hrr.csv Q_* columns only, same-mesh rank comparisons only
HRR_Q_COLS_EXTRA = tuple(c for c in os.environ.get('HRR_Q_RANK_EXTRA', 'HRR').split(',') if c)   # HRR column counted with the Q_* columns (interpretation, to be confirmed)
RANK_SUMMARY_TOL_ALL = os.environ.get('RANK_SUMMARY_TOL_ALL', '1') not in ('0', 'off', 'OFF', 'no', 'false')   # Ruling 2026-10-20, extended: every c_hrr/c_multi_level column + leaf PHI <= RANK_SUMMARY_TOL (same-mesh rank comparisons)
RANK_SUMMARY_TOL = float(os.environ.get('RANK_SUMMARY_TOL', '1e-11'))
RANK_ABS_FLOOR = float(os.environ.get('RANK_ABS_FLOOR', '1e-12'))   # Ruling 2026-10-20, round-off floor: multi_level mass_rel_change, div*: |d| <= max(tol*|ref|, floor) per entry
def coltol(base, name): return HRR_Q_RANK_TOL if base.endswith('_hrr.csv') and name.startswith('Q_') else RANK_SUMMARY_TOL
def floor_col(base, name): return base.endswith('_multi_level.csv') and (name == 'mass_rel_change' or name.startswith('div'))
RULING_MASS_COLS = ('mass', 'mass_rel_change')   # c_multi_level.csv composite-mass round-off: NEEDS RULING (not covered by the 1e-11 rule)
ruling = 0
def needs_ruling(name, msg):
    global ruling
    print(f"  NEEDS-RULING {name}: {msg}"); ruling += 1
HRR_MEAN_TOL = float(os.environ.get('HRR_MEAN_TOL', '0.02')); HRR_MEAN_T0 = float(os.environ.get('HRR_MEAN_T0', '10'))
TOL_PHI = 0.0
SKIP = ('_cpu.csv', '_driver_perf.csv', '_steps.csv')   # wall/CPU timing only (c_driver_steps.csv is kept: solver iterations)
fails = 0

def report(ok, name, msg):
    global fails
    print(f"  {'PASS' if ok else 'FAIL'} {name}: {msg}")
    if not ok: fails += 1

def load_csv(p):
    lines = open(p).read().splitlines()
    first_num = next(i for i, l in enumerate(lines) if l.strip() and re.match(r'^\s*[-+0-9.]', l))
    names = [c.strip() for c in lines[first_num - 1].split(',')]
    data = np.genfromtxt(lines[first_num:], delimiter=',')
    return names, np.atleast_2d(data), lines[first_num:]

# ---- csv files
for rp in sorted(glob.glob(os.path.join(a.ref, f'{a.chid}_*.csv'))):
    base = os.path.basename(rp)
    if base.endswith(SKIP) and not base.endswith('_driver_steps.csv'): continue
    pp = os.path.join(a.run, base)
    if not os.path.exists(pp): report(False, f'csv:{base}', 'missing in run'); continue
    if open(rp).read() == open(pp).read(): report(True, f'csv:{base}', 'BITID'); continue
    try:
        nr, dr, _ = load_csv(rp); npp, dp, _ = load_csv(pp)
    except Exception as e:
        report(False, f'csv:{base}', f'unreadable ({e})'); continue
    if dr.shape != dp.shape or nr != npp:
        if CROSS: print(f"  INFO csv:{base}: shape differs {dr.shape} vs {dp.shape} (not gated across mesh counts)")
        else:
            n = min(len(dr), len(dp)); dif = np.argwhere(~((dp[:n] == dr[:n]) | (np.isnan(dp[:n]) & np.isnan(dr[:n])))) if dr.shape[1] == dp.shape[1] else []
            first = f'; first diff row {int(dif[0][0])} (col0={dr[int(dif[0][0]),0]:g}) column {nr[int(dif[0][1])]}' if len(dif) else ''
            report(False, f'csv:{base}', f'shape/columns differ {dr.shape} vs {dp.shape}{first}')
        continue
    scale = np.maximum(np.nanmax(np.abs(dr), axis=0), 1e-300)
    rel = np.nanmax(np.abs(dp - dr), axis=0) / scale
    j = int(np.nanargmax(rel)); worst = float(rel[j])
    if CROSS: print(f"  INFO csv:{base}: max rel diff {worst:.3e} (column {nr[j]}) (not gated across mesh counts)")
    else:
        dif = np.argwhere(~((dp == dr) | (np.isnan(dp) & np.isnan(dr))))
        cols = sorted(set(int(c) for _, c in dif))
        if RANK_SUMMARY_TOL_ALL and cols and (base.endswith('_hrr.csv') or base.endswith('_multi_level.csv')):
            r0, c0 = int(dif[0][0]), int(dif[0][1]); bad = []
            for c in cols:
                d = np.abs(dp[:, c] - dr[:, c])
                if floor_col(base, nr[c]):
                    lim = np.maximum(RANK_SUMMARY_TOL * np.abs(dr[:, c]), RANK_ABS_FLOOR); ok = bool(np.all((d <= lim) | (np.isnan(d) & np.isnan(dr[:, c]))))
                    print(f"  INFO csv:{base}: {nr[c]} max|diff| {np.nanmax(d):.3e} vs max(1e-11*|ref|, {RANK_ABS_FLOOR:g}) -> {'ok' if ok else 'EXCEEDS'} (round-off floor)")
                else:
                    ok = float(rel[c]) <= coltol(base, nr[c])
                if not ok: bad.append(nr[c])
            rc_ = [c for c in cols if not floor_col(base, nr[c])]
            wc = max(rc_, key=lambda c: float(rel[c]) / coltol(base, nr[c])) if rc_ else cols[0]; wq = float(rel[wc]) if rc_ else 0.0
            report(not bad, f'csv:{base}', f"summary columns ({', '.join(nr[c] for c in cols)}), max rel (to column max) {wq:.3e} in {nr[wc]} vs {coltol(base, nr[wc]):g} (RANK_SUMMARY_TOL; _hrr.csv Q_* HRR_Q_RANK_TOL; floor cols per entry){' FAILING: '+', '.join(bad) if bad else ''}; first diff row {r0} (col0={dr[r0,0]:g}) {nr[c0]}")
            continue
        if base.endswith('_hrr.csv') and cols and all(nr[c].startswith('Q_') or nr[c] in HRR_Q_COLS_EXTRA for c in cols):
            wq = max(float(rel[c]) for c in cols)
            r0, c0 = int(dif[0][0]), int(dif[0][1])
            report(wq <= HRR_Q_RANK_TOL, f'csv:{base}', f"Q columns only ({', '.join(nr[c] for c in cols)}), max rel {wq:.3e} vs HRR_Q_RANK_TOL {HRR_Q_RANK_TOL:g} (Owner ruling 2026-10-20); first diff row {r0} (t={dr[r0,0]:g}) {nr[c0]}")
            continue
        if base.endswith('_multi_level.csv') and cols and all(nr[c] in RULING_MASS_COLS for c in cols):
            r0, c0 = int(dif[0][0]), int(dif[0][1])
            needs_ruling(f'csv:{base}', f"only {', '.join(nr[c] for c in cols)} differ: max rel (to column max) {max(float(rel[c]) for c in cols):.3e}; first diff step {dr[r0,0]:g} {nr[c0]}")
            continue
        r0, c0 = (int(dif[0][0]), int(dif[0][1])) if len(dif) else (-1, -1)
        first = f'first diff row {r0} (t/step col0={dr[r0,0]:g}) column {nr[c0]}' if r0 >= 0 else 'text differs only'
        report(False, f'csv:{base}', f'not bit-identical: max rel diff {worst:.3e} (column {nr[j]}); {first}')

# ---- front x
try:
    nr, dr, _ = load_csv(os.path.join(a.ref, f'{a.chid}_devc.csv')); npp, dp, _ = load_csv(os.path.join(a.run, f'{a.chid}_devc.csv'))
    jx = nr.index('x')
    for t in [float(x) for x in a.times.split(',')]:
        ir = np.where(dr[:, 0] <= t + 1e-9)[0]; ip = np.where(dp[:, 0] <= t + 1e-9)[0]
        if len(ir) == 0 or len(ip) == 0 or dr[ir[-1], 0] < t - 1.0:
            print(f"  INFO front_x t={t:g}: not available (run shorter than t)"); continue
        xr, xp = dr[ir[-1], jx], dp[ip[-1], jx]
        rel = abs(xp - xr) / 120.0
        report(xr == xp, f'front_x t={t:g}', f'ref {xr:g} m, run {xp:g} m (identical required)')
except Exception as e:
    report(False, 'front_x', f'error {e}')

# ---- PHI on leaf cells
def read_plt(d):
    H = open(os.path.join(d, 'Header')).read().split('\n')
    nv = int(H[1]); names = H[2:2 + nv]; i = 2 + nv
    dim = int(H[i]); i += 1; i += 1; finest = int(H[i]); i += 1
    i += 2  # prob lo/hi
    ratios = [int(x) for x in H[i].split()] if finest > 0 else []; i += 1
    dom = [list(map(int, re.findall(r'-?\d+', s)))[:6] for s in re.findall(r'\(\([^)]*\) \([^)]*\) \([^)]*\)\)', H[i])]
    comp = names.index('PHI'); levels = []
    for l in range(finest + 1):
        ch = open(os.path.join(d, f'Level_{l}', 'Cell_H')).read()
        boxes = [list(map(int, re.findall(r'-?\d+', b)))[:6] for b in re.findall(r'\(\(\s*-?\d+,-?\d+,-?\d+\) \(-?\d+,-?\d+,-?\d+\) \(\d+,\d+,\d+\)\)', ch)]
        fods = re.findall(r'FabOnDisk: (\S+) (\d+)', ch)
        fabs = []
        for b, (fn, off) in zip(boxes, fods):
            with open(os.path.join(d, f'Level_{l}', fn), 'rb') as f:
                f.seek(int(off)); hdr = f.readline().decode()
                lo = b[:3]; hi = b[3:]; n = [hi[k] - lo[k] + 1 for k in range(3)]; nc = int(hdr.split()[-1])
                arr = np.fromfile(f, dtype='<f8', count=n[0] * n[1] * n[2] * nc).reshape(nc, n[2], n[1], n[0])
                fabs.append((lo, hi, arr[comp].transpose(2, 1, 0)))
        levels.append(fabs)
    return dom, ratios, levels

def leaf_field(d):
    dom, ratios, levels = read_plt(d)
    R = int(np.prod(ratios)) if ratios else 1
    n0 = [dom[0][3 + k] - dom[0][k] + 1 for k in range(3)]
    F = np.full([n0[k] * R for k in range(3)], np.nan)
    rr = 1
    for l, fabs in enumerate(levels):   # coarse to fine: finer levels overwrite => each point holds its leaf cell value
        s = R // rr
        for lo, hi, a3 in fabs:
            F[lo[0]*s:(hi[0]+1)*s, lo[1]*s:(hi[1]+1)*s, lo[2]*s:(hi[2]+1)*s] = np.repeat(np.repeat(np.repeat(a3, s, 0), s, 1), s, 2)
        if l < len(ratios): rr *= ratios[l]
    return F

pr = sorted(glob.glob(os.path.join(a.ref, '**', 'plt_*_final'), recursive=True))
pp = sorted(glob.glob(os.path.join(a.run, '**', 'plt_*_final'), recursive=True))
if pr and pp:
    try:
        Fr, Fp = leaf_field(pr[0]), leaf_field(pp[0])
        if Fr.shape != Fp.shape: report(CROSS, 'phi', f'leaf grid shape differs {Fr.shape} vs {Fp.shape}')
        else:
            dmax = float(np.nanmax(np.abs(Fp - Fr))); sc = float(np.nanmax(np.abs(Fr)))
            if CROSS: print(f"  INFO phi: leaf max|dPHI| {dmax:.3e}, rel {dmax/sc:.3e} (not gated across mesh counts){' NON-FINITE/GARBAGE (>1e3)' if not np.isfinite(dmax) or np.nanmax(np.abs(Fp)) > 1e3 else ''}")
            elif dmax == 0: report(True, 'phi', 'leaf PHI identical')
            elif RANK_SUMMARY_TOL_ALL: report(dmax / sc <= RANK_SUMMARY_TOL, 'phi', f'leaf max|dPHI| {dmax:.3e}, rel {dmax/sc:.3e} vs RANK_SUMMARY_TOL {RANK_SUMMARY_TOL:g} (Ruling 2026-10-20, extended)')
            else: needs_ruling('phi', f'leaf max|dPHI| {dmax:.3e}, rel {dmax/sc:.3e} (not covered by the 1e-11 rule; strict bit-identity otherwise)')
    except Exception as e:
        report(False, 'phi', f'error {e}')
else:
    print('  INFO phi: no final plotfile (1-level --run writes none); PHI covered by the devc/front-x and csv checks')

# ---- gates: run must complete; crossmesh: same verdict as the reference; bitid: verdict PASS (csv bit-identity covers the rest)
def verdict(d):
    log = open(os.path.join(d, 'run.log'), errors='replace').read()
    if 'LEVEL RESULT' in log:
        g = re.search(r'LEVEL RESULT gate (PASS|FAIL)[^\n]*', log)
        return ('PASS' if ('LEVEL RUN PASS' in log and g and g.group(1) == 'PASS') else 'FAIL'), (g.group(0)[:160] if g else 'no RESULT gate line')
    return ('PASS' if 'completed successfully' in log else 'FAIL'), ('FDS completed' if 'completed successfully' in log else 'run did not complete')
vr, mr = verdict(a.ref); vp, mp = verdict(a.run)
if CROSS: report(vr == vp, 'gate verdict', f'ref {vr}, run {vp}: {mp}')
else: report(vp == 'PASS', 'gate', mp)
if CROSS:
    try:
        hr = np.genfromtxt(os.path.join(a.ref, f'{a.chid}_hrr.csv'), delimiter=',', skip_header=2)
        hp = np.genfromtxt(os.path.join(a.run, f'{a.chid}_hrr.csv'), delimiter=',', skip_header=2)
        mr_ = hr[hr[:, 0] >= HRR_MEAN_T0, 1].mean(); mp_ = hp[hp[:, 0] >= HRR_MEAN_T0, 1].mean()
        rel = abs(mp_ - mr_) / abs(mr_)
        report(rel <= HRR_MEAN_TOL, 'mean HRR', f't>={HRR_MEAN_T0:g} s ref {mr_:.5g} kW, run {mp_:.5g} kW, rel {rel:.3e}, tol {HRR_MEAN_TOL:g}')
    except Exception as e:
        report(False, 'mean HRR', f'error {e}')
sys.exit(1 if fails else (3 if ruling else 0))   # 3 = only NEEDS-RULING items
