#!/usr/bin/env python3
"""UII (INTEGRATED INTENSITY, kW/m2) slice y=Y: FDS .sf frames vs fds_amr RADSLC log lines (FDS_RAD_SLICE=y,t; W/m2 -> kW/m2), grid outlines for AMR.
usage: val_plot.py OUT.png T --fds c=DIR f1=DIR [f2=DIR] --amr 1=DIR 2=DIR [3=DIR]   (DIR with in_*_1.sf resp. run.log; Python from /workspace/venv_plot)"""
import sys, glob, struct, argparse
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def read_sf(fn):
    b = open(fn, 'rb').read(); p = 0; lab = []
    for _ in range(3): n = struct.unpack('i', b[p:p+4])[0]; lab.append(b[p+4:p+4+n].decode().strip()); p += 8 + n
    p += 4; I = struct.unpack('6i', b[p:p+24]); p += 28
    ni, nj, nk = I[1]-I[0]+1, I[3]-I[2]+1, I[5]-I[4]+1; n = ni*nj*nk; fr = []
    while p + 16 + 4*n <= len(b):
        t = struct.unpack('f', b[p+4:p+8])[0]; p += 16; fr.append((t, np.frombuffer(b[p:p+4*n], '<f4').reshape(ni, nj, nk, order='F'))); p += 4*n + 4
    return I, fr

def fds_slice(d, T):
    """cell-centered UII frame nearest T; returns list of (x0,z0,dx,dz,val); the slice is a single y cell plane"""
    out = []
    for fn in sorted(glob.glob(f'{d}/in_*_1.sf')):
        I, fr = read_sf(fn); t, a = min(fr, key=lambda f: abs(f[0]-T)); out.append((t, a[:, 0, :] if a.shape[1] == 1 else a[:, a.shape[1]//2, :]))
    return out

def main():
    p = argparse.ArgumentParser(); p.add_argument('out'); p.add_argument('T', type=float); p.add_argument('--fds', nargs='+'); p.add_argument('--amr', nargs='+')
    p.add_argument('--L', type=float, nargs=3, default=[3.6, 2.4, 2.4]); p.add_argument('--n0', type=int, nargs=2, default=[36, 24])
    o = p.parse_args(); Lx, Ly, Lz = o.L; panels = []
    for kv in o.fds or []:
        k, d = kv.split('='); parts = fds_slice(d, o.T); a = parts[0][1]; nx = a.shape[0]
        if len(parts) > 1:     # 4 meshes (2 x 1 x 2 in x,z) ordered mesh1:(x0,z0) mesh2:(x0,z1) mesh3:(x1,z0) mesh4:(x1,z1)
            nz = a.shape[1]; A = np.zeros((2*nx, 2*nz))
            for i, (t, m) in enumerate(parts): A[(i//2)*nx:(i//2+1)*nx, (i % 2)*nz:(i % 2+1)*nz] = m
            a = A
        panels.append((f'FDS {k} ({a.shape[0]}x{a.shape[1]}) t={parts[0][0]:.2f}', [(0, 0, Lx/a.shape[0], Lz/a.shape[1], a)], None))
    for kv in o.amr or []:
        k, d = kv.split('='); cells = []
        for l in open(f'{d}/run.log', errors='replace'):
            if l.startswith('RADSLC'):
                f = l.split(); cells.append((int(f[3]), float(f[4]), float(f[5]), float(f[6]), float(f[7]), float(f[8])*1e-3, int(f[9]), float(f[1])))
        if cells: panels.append((f'AMR {k}-lev t={cells[0][7]:.2f}', None, cells))
    vals = [v for _, g, c in panels for v in ((np.concatenate([x[4].ravel() for x in g]),) if g else (np.array([x[5] for x in c if not x[6]]),))]
    vmin, vmax = min(v.min() for v in vals), max(v.max() for v in vals)
    n = len(panels); fig, ax = plt.subplots(1, n, figsize=(3.6*n, 3.2), squeeze=False)
    for a_, (title, g, c) in zip(ax[0], panels):
        if g:
            x0, z0, dx, dz, A = g[0]; im = a_.imshow(A.T, origin='lower', extent=[0, Lx, 0, Lz], vmin=vmin, vmax=vmax, cmap='inferno', aspect='equal')
        else:
            for lev in sorted({x[0] for x in c}):
                for (L_, xc, zc, dx, dz, u, cov, t) in c:
                    if L_ != lev or cov: continue
                    a_.add_patch(Rectangle((xc-dx/2, zc-dz/2), dx, dz, facecolor=plt.cm.inferno((u-vmin)/(vmax-vmin)), edgecolor='none'))
            for L_ in sorted({x[0] for x in c}):   # outline: one rectangle edge per leaf cell of level >= 1
                if L_ == 0: continue
                for (l2, xc, zc, dx, dz, u, cov, t) in c:
                    if l2 == L_: a_.add_patch(Rectangle((xc-dx/2, zc-dz/2), dx, dz, fill=False, edgecolor='c' if L_ == 1 else 'lime', lw=0.15, alpha=0.6))
            a_.set_xlim(0, Lx); a_.set_ylim(0, Lz); a_.set_aspect('equal'); im = plt.cm.ScalarMappable(norm=plt.Normalize(vmin, vmax), cmap='inferno')
        a_.set_title(title, fontsize=8); a_.set_xlabel('x (m)', fontsize=7); a_.tick_params(labelsize=6)
    fig.colorbar(im, ax=ax[0].tolist(), shrink=0.8, label='UII (kW/m2)'); fig.suptitle('UII slice y=1.2 m', fontsize=9); fig.savefig(o.out, dpi=130, bbox_inches='tight')

if __name__ == '__main__': main()
