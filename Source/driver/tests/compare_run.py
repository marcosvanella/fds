#!/usr/bin/env python3
"""Compare a driver run (fds_amr --run) with the single-mesh FDS baseline of a 2-D MMS case (M2a, S5).

Kernel-facing rules: (a) passive scalars are not involved here; (b) uniform Cartesian metrics only.

usage: compare_run.py <driver_outdir> <baseline_dir> <chid> [--dump ref.dump]   (single-mesh baselines only; for a multi-mesh case use its __1mesh baseline)
Reports: T/DT against <chid>_steps.csv (DT 3 digits, T 7 digits), the dump STEP records (bitwise), the final fields against the restart file, the mass file,
and the MMS errors of the shunn3 family (analytic solution) for driver and baseline with the T2 statement.
"""
import struct, sys, os
import numpy as np

def restart_fields(fn, n):
    """valid-cell arrays of the FDS restart file (single mesh IBAR=JBAR=1.., order of dump.f90): U V W D H US VS WS DS HS RHO TMP ... ZZ"""
    f = open(fn, 'rb').read()
    recs = []; p = 0
    while p < len(f):
        m = struct.unpack('<i', f[p:p + 4])[0]
        recs.append(f[p + 4:p + 4 + m]); p += 8 + m
    nx, ny, nz = n
    def arr(r, lb, ext):
        a = np.frombuffer(r, dtype='<f8').reshape(ext[::-1]).transpose(2, 1, 0)   # a[i-lb0, j-lb1, k-lb2]
        return a
    ib, jb, kb = nx + 1, ny + 1, nz + 1
    out = {}
    names = ['U', 'V', 'W', 'D', 'H', 'US', 'VS', 'WS', 'DS', 'HS']
    lbs = {'U': (-1, 0, 0), 'V': (0, -1, 0), 'W': (0, 0, -1)}
    exts = {'U': (nx + 3, ny + 2, nz + 2), 'V': (nx + 2, ny + 3, nz + 2), 'W': (nx + 2, ny + 2, nz + 3)}
    for i, nm in enumerate(names):
        b = nm[0] if nm[0] in 'UVW' else None
        key = nm.rstrip('S') if nm in ('US', 'VS', 'WS') else nm
        lb = lbs.get(key, (0, 0, 0)); ext = exts.get(key, (nx + 2, ny + 2, nz + 2))
        a = arr(recs[i], lb, ext)
        # valid data: cells 1..n (faces U(I), I=1..nx = high face of cell I)
        out[nm] = a[1 - lb[0]:nx + 1 - lb[0], 1 - lb[1]:ny + 1 - lb[1], 1 - lb[2]:nz + 1 - lb[2]]
    ext = (nx + 4, ny + 4, nz + 4)
    out['RHO'] = arr(recs[10], (-1, -1, -1), ext)[2:nx + 2, 2:ny + 2, 2:nz + 2]
    out['TMP'] = arr(recs[11], (-1, -1, -1), ext)[2:nx + 2, 2:ny + 2, 2:nz + 2]
    return out

def load_bin(d, chid, nm, n):
    a = np.fromfile('%s/%s_final_%s.bin' % (d, chid, nm), dtype='<f8')
    return a.reshape((n[2], n[1], n[0])).transpose(2, 1, 0)

def shunn(nx, t):
    r0, r1, uf, vf, k, w, L = 5.0, 1.0, 0.5, 0.5, 2.0, 2.0, 2.0
    dx = L / nx
    xc = -1 + (np.arange(nx) + 0.5) * dx
    def fields(x, y):
        s = np.sin(np.pi * k * (x - uf * t)) * np.sin(np.pi * k * (y - vf * t)) * np.cos(np.pi * w * t)
        z = (1 + s) / ((1 + r0 / r1) + (1 - r0 / r1) * s)
        rho = 1 / (z / r1 + (1 - z) / r0)
        u = uf + (r1 - r0) / rho * (-w / (4 * k)) * np.cos(np.pi * k * (x - uf * t)) * np.sin(np.pi * k * (y - vf * t)) * np.sin(np.pi * w * t)
        v = vf + (r1 - r0) / rho * (-w / (4 * k)) * np.sin(np.pi * k * (x - uf * t)) * np.cos(np.pi * k * (y - vf * t)) * np.sin(np.pi * w * t)
        return z, rho, u, v
    X, Y = np.meshgrid(xc, xc, indexing='ij')
    z, rho, u, v = fields(X, Y)
    _, _, uu, _ = fields(X + dx / 2, Y)      # U is compared at the face x[i+1]
    _, _, _, vv = fields(X, Y)               # W (second coordinate) at the cell centre
    H = 0.5 * (u - uf) * (v - vf)
    return dict(rho=rho, z=z, u=uu, w=vv, H=H), dx

def mms_errors(fn, nx, t_expected=None):
    lines = open(fn).read().split('\n')
    t = float(lines[1])
    d = np.array([[float(x) for x in l.split(',')] for l in lines[2:] if l.strip()])
    sol, dx = shunn(nx, t)
    # rows: K outer, J, I inner; ny = 1: row index = k*nx + i ; FDS z <-> second coordinate
    e = {}
    for c, key in enumerate(['rho', 'z', 'u', 'w', 'H']):
        a = d[:, c].reshape(nx, nx)   # [k, i]
        ref = sol[key].T              # sol[i,k] -> [k,i]
        e[key] = np.linalg.norm((a - ref).ravel()) / nx
    return t, e

def main():
    out, base, chid = sys.argv[1:4]
    dump = sys.argv[sys.argv.index('--dump') + 1] if '--dump' in sys.argv else None
    drows = {int(l.split(',')[0]): l.strip().split(',') for l in list(open('%s/%s_driver_steps.csv' % (out, chid)))[1:]}
    brows = {}
    for l in list(open('%s/%s_steps.csv' % (base, chid)))[2:]:
        p = [x.strip() for x in l.split(',')]; brows[int(p[0])] = p
    bad = []
    for k, p in sorted(brows.items()):
        if k not in drows: continue
        dt = float(drows[k][2]); t = float(drows[k][1])
        m, e = ('%.2E' % dt).split('E'); fm = '0.%s%sE%+03d' % (m[0], m[2:4], int(e) + 1)
        if fm != p[2]: bad.append(('DT', k, dt, p[2]))
        if '%.7f' % t != p[3]: bad.append(('T', k, t, p[3]))
    print('steps.csv: %d rows compared, %d mismatches (DT 3 digits, T 7 digits) %s' % (len(brows), len(bad), bad[:4]))
    print('driver steps %d, baseline steps %s, final T %s' % (len(drows), 'n/a', drows[max(drows)][1]))
    if dump:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'refdump')); sys.dont_write_bytecode = True
        import rd
        recs = rd.read(dump, False)
        st = [r for r in recs if r['name'] == 'STEP']
        nb = 0; mx = 0.0
        for r in st:
            k = r['icyc']
            if k not in drows: continue
            t0 = float(drows[k - 1][1]) if k > 1 else 0.0
            dt = float(drows[k][2])
            if not (dt == r['DT'] and t0 == r['T']): nb += 1
            mx = max(mx, abs(dt - r['DT']) / r['DT'])
        print('dump STEP records: %d of %d bitwise equal in (T at start, DT); max relative DT difference %.3e' % (len(st) - nb, len(st), mx))
    # mms
    bmms = '%s/%s_mms.csv' % (base, chid); dmms = '%s/%s_mms.csv' % (out, chid)
    if os.path.exists(bmms) and os.path.exists(dmms):
        nx = int(open(bmms).readline().split(',')[1])
        tb, eb = mms_errors(bmms, nx); td, ed = mms_errors(dmms, nx)
        print('MMS at T=%.9f (baseline) / %.9f (driver)' % (tb, td))
        for key, nm in (('rho', 'e_rho'), ('z', 'e_Z'), ('u', 'e_u'), ('H', 'e_H')):
            ok = ed[key] <= 1.05 * eb[key]
            print('  %-6s baseline %.4e driver %.4e ratio %.4f  %s' % (nm, eb[key], ed[key], ed[key] / eb[key], 'T2 pass (e <= 1.05 e_base)' if ok else 'T2 FAIL (e > 1.05 e_base)'))
    # final fields against the restart file
    rs = '%s/%s_1.restart' % (base, chid)
    if os.path.exists(rs) and os.path.exists('%s/%s_final_U.bin' % (out, chid)):
        nxx = int(open(bmms).readline().split(',')[1]); n = (nxx, 1, nxx)
        ref = restart_fields(rs, n)
        print('final fields vs baseline restart (valid cells): max|d| and max|d|/max|field|')
        for nm in ['U', 'V', 'W', 'H', 'HS', 'US', 'VS', 'WS', 'D', 'DS', 'RHO', 'TMP']:
            a = load_bin(out, chid, nm, n)
            d = np.abs(a - ref[nm]); sc = max(np.abs(ref[nm]).max(), 1e-300)
            print('  %-4s %.3e  %.3e' % (nm, d.max(), d.max() / sc))
    # mass
    bm = [l.strip().split(',') for l in open('%s/%s_mass.csv' % (base, chid))][2:]
    dm = [l.strip().split(',') for l in open('%s/%s_driver_mass.csv' % (out, chid))][1:]
    print('mass: baseline rows %d driver rows %d; total mass first/last baseline %s %s, driver %s %s' % (len(bm), len(dm), bm[0][1], bm[-1][1], dm[0][1], dm[-1][1]))
    mb = np.array([float(r[1]) for r in bm]); md = np.array([float(r[1]) for r in dm])
    print('  max |total mass driver - baseline final row| %.3e ; driver drift max|M-M0| %.3e ; baseline drift %.3e' % (abs(md[-1] - mb[-1]), np.abs(md - md[0]).max(), np.abs(mb - mb[0]).max()))

if __name__ == '__main__':
    main()
