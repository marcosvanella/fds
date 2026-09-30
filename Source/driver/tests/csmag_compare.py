#!/usr/bin/env python3
"""csmag_32 (3-D, six periodic faces) step check against a reference-dump run of the SAME input with &PRES FISHPAK_BC=0,0,0 (true periodic Crayfishpak FFT solve).
usage: csmag_compare.py <driver-run-dir> <reference-dump> [tol_rel]
The driver run must have been made with FDSTL_STAGE=1 and FDSTL_STAGE=2 output in <dir>/s1 and <dir>/s2 (run_csmag_check.sh does that).
Gate: FVX, H, HS, DS of step 1 (pass 1, pass 2, corrector) and step 2 agree with the reference to tol_rel*max|ref| (default 1e-12, absolute floor 1e-12 for the ~1e-15 divergence)."""
import sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, 'refdump'))
import rd

run, dump = sys.argv[1], sys.argv[2]
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-12
recs = rd.read(dump, True)

def interior(a):
    d = a[3]; d = d[..., 0] if d.ndim == 4 else d
    o = [1 - lb for lb in a[0][:3]]
    return d[o[0]:o[0] + 32, o[1]:o[1] + 32, o[2]:o[2] + 32]

def get(name, ic, first, lst, f):
    for r in recs:
        if r['name'] == name and r['icyc'] == ic and r['first'] == first:
            return r[lst][f]

def ld(d, s):
    return np.fromfile('%s/stage_%s.bin' % (d, s)).reshape((32, 32, 32), order='F')

rows = [('step1 pass1 FVX', 'VFLUX_P', 1, 1, 'aft', 'FVX', 's1', 'p1_vflux_FVX'),
        ('step1 pass1 H', 'VPRED', 1, 1, 'bef', 'H', 's1', 'p1_press_H'),
        ('step1 pass2 H', 'VPRED', 1, 0, 'bef', 'H', 's1', 'p2_press_H'),
        ('step1 pass2 DS', 'VPRED', 1, 0, 'bef', 'DS', 's1', 'p2_div_DS'),
        ('step1 corrector HS', 'VCORR', 1, 0, 'bef', 'HS', 's1', 'c_press_HS'),
        ('step1 corrector U', 'VCORR', 1, 0, 'aft', 'U', 's1', 'c_end_U'),
        ('step1 corrector W', 'VCORR', 1, 0, 'aft', 'W', 's1', 'c_end_W'),
        ('step2 FVX', 'VFLUX_P', 2, 1, 'aft', 'FVX', 's2', 'p1_vflux_FVX'),
        ('step2 DS', 'VPRED', 2, 1, 'bef', 'DS', 's2', 'p1_div_DS'),
        ('step2 H', 'VPRED', 2, 1, 'bef', 'H', 's2', 'p1_press_H'),
        ('step2 corrector HS', 'VCORR', 2, 1, 'bef', 'HS', 's2', 'c_press_HS'),
        ('step2 corrector U', 'VCORR', 2, 1, 'aft', 'U', 's2', 'c_end_U'),
        ('step2 corrector W', 'VCORR', 2, 1, 'aft', 'W', 's2', 'c_end_W')]
bad = 0
for tag, rec, ic, fi, lst, f, d, s in rows:
    a = get(rec, ic, fi, lst, f)
    if a is None:
        print('%-22s NOT IN DUMP' % tag); bad += 1; continue
    ref = interior(a); mine = ld(os.path.join(run, d), s)
    e = float(np.abs(mine - ref).max()); lim = max(tol * float(np.abs(ref).max()), 1e-12)
    ok = e <= lim; bad += (not ok)
    print('%-22s max|d| %.3e  (ref max %.3e, limit %.1e) %s' % (tag, e, np.abs(ref).max(), lim, 'OK' if ok else 'FAIL'))
print('CSMAG COMPARE %s' % ('PASS' if bad == 0 else 'FAIL'))
sys.exit(1 if bad else 0)
