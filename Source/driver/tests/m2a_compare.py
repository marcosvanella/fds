#!/usr/bin/env python3
"""T2 comparison of a driver run with an official baseline directory (M2a, S8). Kernel-facing rules: (a) passive scalars not involved; (b) uniform Cartesian metrics only.

usage: m2a_compare.py csmag <driver_outdir> <baseline_dir> <chid> [--restart <baseline-rerun>_1.restart] [--tol-ke REL] [--tol-field REL]
       m2a_compare.py glmat <driver_outdir> <baseline_dir_4mesh> <chid>

csmag (no analytic solution): (1) T and DT of every step that <chid>_steps.csv lists equal to its printed digits (DT 3 digits, T as many digits as the file prints, the existing gate; FDS lists the diagnostic steps only);
  (2) the KE device series (SPATIAL_STATISTIC='MEAN') of <chid>_devc.csv: relative difference per row; (3) final fields against a restart file when given (the baseline
  directory keeps no restart file, so the file comes from a rerun of the reference binary on the same input, whose _devc.csv is checked byte-identical to the baseline by the caller).
  The gates are --tol-ke (default 1e-8) and --tol-field (default 1e-8, relative to the field maximum, absolute floor 1e-12); they are the same for the plain and the periodic variant, the verdict
  line says which variant it is. Exit code 0 only if every gate holds.
csmag-t2 <driver_outdir> <baseline_plain> <chid> --alt <baseline_periodic_variant> --restart <plain rerun restart> --alt-restart <periodic rerun restart>:
  T2 form e_driver <= 1.05 e_base (+ absolute floor 1e-12) for a case without analytic solution or Verification Guide metric: e is the deviation from the plain baseline
  of (a) the KE device series (linear interpolation in time onto the baseline rows, relative to KE) and (b) each final field (relative to the field maximum); e_base is the same
  deviation of the baseline's own periodic variant (`&PRES FISHPAK_BC=0,0,0`, the same FDS binary and input apart from that line), i.e. the spread FDS itself shows between two
  pressure-solver settings of this input. The driver solves the periodic problem, so e_driver equals e_base up to round-off by construction: this is a "no worse than the
  solver spread" statement, and no fixed numerical bar for csmag_32 exists in the requirements.
glmat (shunn3_4mesh_32, 4 meshes 16x1x16 in a 2x2 layout of x and z): (1) T/DT as above; (2) the MMS error of mesh 1 (the only mesh <chid>_mms.csv of the baseline holds)
  for the driver (quadrant of its level-wide MMS file) must satisfy e <= 1.05 e_base (T2 rule); (3) final fields against the four restart files, assembled.
"""
import os, struct, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.dont_write_bytecode = True
import compare_run as cr

def steps_check(out, base, chid):
    drows = {int(l.split(',')[0]): l.strip().split(',') for l in list(open('%s/%s_driver_steps.csv' % (out, chid)))[1:]}
    brows = {}
    for l in list(open('%s/%s_steps.csv' % (base, chid)))[2:]:
        p = [x.strip() for x in l.split(',')]; brows[int(p[0])] = p
    bad = []
    for k, p in sorted(brows.items()):
        if k not in drows: bad.append(('missing', k)); continue
        dt = float(drows[k][2]); t = float(drows[k][1])
        m, e = ('%.2E' % dt).split('E'); fm = '0.%s%sE%+03d' % (m[0], m[2:4], int(e) + 1)
        if fm != p[2]: bad.append(('DT', k, dt, p[2]))
        nd = len(p[3].split('.')[1])   # FDS prints T with F10.<digits>, digits = 7 in DNS mode and at most 5 in LES mode, fewer for larger T (dump.f90 WRITE_DIAGNOSTICS)
        if ('%.' + str(nd) + 'f') % t != p[3]: bad.append(('T', k, t, p[3]))
    n = len(brows)
    print('T/DT: %d baseline rows, driver rows %d, mismatches %d %s' % (n, len(drows), len(bad), bad[:3]))
    return len(bad) == 0 and all(k in drows for k in brows)

def devc(fn):
    rows = [l.strip().split(',') for l in open(fn)][2:]
    return np.array([[float(x) for x in r] for r in rows if r and r[0]])

def csmag_t2(out, base, chid, alt, rs, rs_alt):
    n = tuple(int(x) for x in open('%s/%s.fds' % (base, chid)).read().split('IJK=')[1].split('/')[0].split(',')[:3])
    b = devc('%s/%s_devc.csv' % (base, chid))
    def ke_dev(fn):
        d = devc(fn)
        ki = np.interp(b[:, 0], d[:, 0], d[:, 1])
        return (np.abs(ki - b[:, 1]) / np.abs(b[:, 1])).max()
    e_ke_drv = ke_dev('%s/%s_devc.csv' % (out, chid)); e_ke_base = ke_dev('%s/%s_devc.csv' % (alt, chid))
    print('KE series (interpolated onto the %d baseline rows), max relative deviation from the plain baseline: driver %.4e, baseline periodic variant %.4e, ratio %.6f' % (len(b), e_ke_drv, e_ke_base, e_ke_drv / e_ke_base))
    ok = e_ke_drv <= 1.05 * e_ke_base + 1e-12
    print('  %s' % ('T2 pass (e <= 1.05 e_base)' if e_ke_drv <= 1.05 * e_ke_base + 1e-12 else 'T2 FAIL'))
    ref = cr.restart_fields(rs, n); ra = cr.restart_fields(rs_alt, n)
    print('final fields, max deviation from the plain baseline rerun restart, relative to the field maximum (absolute floor 1e-12):')
    for nm in ['U', 'V', 'W', 'H', 'HS', 'US', 'VS', 'WS', 'D', 'DS', 'RHO', 'TMP']:
        sc = max(np.abs(ref[nm]).max(), 1e-300)
        ed = np.abs(cr.load_bin(out, chid, nm, n) - ref[nm]).max(); eb = np.abs(ra[nm] - ref[nm]).max()
        good = ed <= 1.05 * eb + 1e-12; ok = ok and good
        print('  %-4s driver %.3e (rel %.2e)  baseline periodic variant %.3e (rel %.2e)  ratio %.5f  %s' % (nm, ed, ed / sc, eb, eb / sc, ed / max(eb, 1e-300), 'T2 pass' if good else 'T2 FAIL'))
    print('M2A COMPARE csmag-t2 %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1

def main():
    mode, out, base, chid = sys.argv[1:5]
    args = sys.argv[5:]
    opt = lambda k, d: args[args.index(k) + 1] if k in args else d
    if mode == 'csmag-t2':
        sys.exit(csmag_t2(out, base, chid, opt('--alt', None), opt('--restart', None), opt('--alt-restart', None)))
    ok = steps_check(out, base, chid) if mode != 'glmat' else True
    if mode == 'csmag':
        tol_ke = float(opt('--tol-ke', 1e-8)); tol_f = float(opt('--tol-field', 1e-8))
        b = devc('%s/%s_devc.csv' % (base, chid)); d = devc('%s/%s_devc.csv' % (out, chid))
        same = b.shape == d.shape and np.array_equal(b[:, 0], d[:, 0])
        print('KE series: baseline rows %d driver rows %d, times equal: %s' % (len(b), len(d), same))
        if same:
            rel = np.abs(d[:, 1] - b[:, 1]) / np.abs(b[:, 1])
            print('  KE relative difference per row: max %.3e, final row %.3e (baseline %.7E driver %.7E); gate %.1e' % (rel.max(), rel[-1], b[-1, 1], d[-1, 1], tol_ke))
            ok = ok and rel.max() <= tol_ke
        else:
            ok = False
        rs = opt('--restart', None)
        if rs:
            n = tuple(int(x) for x in open('%s/%s.fds' % (base, chid)).read().split('IJK=')[1].split('/')[0].split(',')[:3])
            ref = cr.restart_fields(rs, n)
            worst = 0.0
            print('final fields vs the reference rerun restart (valid cells): max|d|, max|d|/max|field|; gate %.1e' % tol_f)
            for nm in ['U', 'V', 'W', 'H', 'HS', 'US', 'VS', 'WS', 'D', 'DS', 'RHO', 'TMP']:
                a = cr.load_bin(out, chid, nm, n)
                dd = np.abs(a - ref[nm]); sc = max(np.abs(ref[nm]).max(), 1e-300)
                lim = max(tol_f * sc, 1e-12)   # absolute floor 1e-12 for fields that are zero up to round-off (D, DS)
                print('  %-4s %.3e  %.3e  limit %.1e %s' % (nm, dd.max(), dd.max() / sc, lim, 'ok' if dd.max() <= lim else 'FAIL')); worst = max(worst, dd.max() / lim)
            ok = ok and worst <= 1.0
    elif mode == 'glmat':
        # T/DT of the multi-mesh baseline differ from a single level solve from step 1 on (different pressure coupling at the mesh interfaces, ruling: T2-only): informational
        ok = True
        brow = {}
        for l in list(open('%s/%s_steps.csv' % (base, chid)))[2:]:
            q = [x.strip() for x in l.split(',')]; brow[int(q[0])] = (float(q[3]), float(q[2]))
        dr = {int(l.split(',')[0]): (float(l.split(',')[1]), float(l.split(',')[2])) for l in list(open('%s/%s_driver_steps.csv' % (out, chid)))[1:]}
        rt = max(abs(dr[k][0] - v[0]) / v[0] for k, v in brow.items()); rd = max(abs(dr[k][1] - v[1]) / v[1] for k, v in brow.items())
        print('T/DT is informational for the 4-mesh baseline (T2-only by ruling): max relative T difference %.2e, DT %.2e over the %d listed steps; driver took %d steps, baseline %d' % (rt, rd, len(brow), len(dr), max(brow)))
        # mesh k (1-based) lower corner in cells: (ix, iz)
        corner = {1: (0, 0), 2: (16, 0), 3: (16, 16), 4: (0, 16)}
        n = (32, 1, 32); nm_ = (16, 1, 16)
        # MMS of mesh 1
        bl = open('%s/%s_mms.csv' % (base, chid)).read().split('\n')
        tb = float(bl[1]); bd = np.array([[float(x) for x in l.split(',')] for l in bl[2:] if l.strip()])
        solb, _ = cr.shunn(32, tb)
        eb = {}
        for c, key in enumerate(['rho', 'z', 'u', 'w', 'H']):
            eb[key] = np.linalg.norm((bd[:, c].reshape(16, 16) - solb[key].T[:16, :16]).ravel()) / 16   # mesh 1 = cells 1..16 in x and z
        lines = open('%s/%s_mms.csv' % (out, chid)).read().split('\n')
        td = float(lines[1]); dd = np.array([[float(x) for x in l.split(',')] for l in lines[2:] if l.strip()])
        sol, _ = cr.shunn(32, td)
        ed = {}
        for c, key in enumerate(['rho', 'z', 'u', 'w', 'H']):
            a = dd[:, c].reshape(32, 32)[:16, :16]; ref = sol[key].T[:16, :16]
            ed[key] = np.linalg.norm((a - ref).ravel()) / 16
        print('MMS (mesh 1 quadrant) at T=%.9f (baseline) / %.9f (driver)' % (tb, td))
        for key, nme in (('rho', 'e_rho'), ('z', 'e_Z'), ('u', 'e_u'), ('H', 'e_H')):
            good = ed[key] <= 1.05 * eb[key]; ok = ok and good
            print('  %-6s baseline %.4e driver %.4e ratio %.4f  %s' % (nme, eb[key], ed[key], ed[key] / eb[key], 'T2 pass (e <= 1.05 e_base)' if good else 'T2 FAIL (e > 1.05 e_base)'))
        ref = {}
        for m in (1, 2, 3, 4):
            f = cr.restart_fields('%s/%s_%d.restart' % (base, chid, m), nm_)
            for nme, a in f.items():
                if nme not in ref: ref[nme] = np.zeros(n)
                ix, iz = corner[m]
                ref[nme][ix:ix + 16, :, iz:iz + 16] = a
        print('final fields vs the assembled 4 baseline restart files (valid cells): max|d|, max|d|/max|field|')
        for nme in ['U', 'W', 'H', 'HS', 'D', 'DS', 'RHO', 'TMP']:
            a = cr.load_bin(out, chid, nme, n)
            x = np.abs(a - ref[nme]); sc = max(np.abs(ref[nme]).max(), 1e-300)
            print('  %-4s %.3e  %.3e' % (nme, x.max(), x.max() / sc))
        bm = [l.strip().split(',') for l in open('%s/%s_mass.csv' % (base, chid))][2:]
        dm = [l.strip().split(',') for l in open('%s/%s_mass.csv' % (out, chid))][2:]
        mb = np.array([float(r[1]) for r in bm]); md = np.array([float(r[1]) for r in dm])
        print('mass: baseline drift max|M-M0| %.3e, driver drift %.3e' % (np.abs(mb - mb[0]).max(), np.abs(md - md[0]).max()))
    print('M2A COMPARE %s %s' % (mode, 'PASS' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
