#!/usr/bin/env python3
"""Generates a fully periodic 3-D LES case of N^3 cells that the driver accepts (the csmag_32 recipe at another resolution, with a Taylor-Green initial velocity).
usage: make_periodic_case.py <N> <outdir> [t_end] [sx sy sz]     writes <outdir>/tg<N>.fds and <outdir>/tg<N>_uvw.csv (N^3 lines u,v,w, I fastest, then J, then K: UVW_INIT, PERIODIC_TEST=2)
The case is csmag_32 with IJK=N,N,N (same box 0.56549 m, six PERIODIC vents (PB form, valid for several meshes), constant Smagorinsky), no slices (they dominate the output size at large N), one KE device.
The velocity is divergence-free and smooth: u = A sin(kx) cos(ky) cos(kz), v = -A cos(kx) sin(ky) cos(kz), w = 0, A = 0.3 m/s, k = 2 pi / L, sampled at the cell centres.
Run (1 rank): mpirun -np 1 fds_amr tg<N>.fds --run --outdir . --chid tg<N> [--steps n]; the run writes tg<N>_driver_perf.csv (t_loop_s, f_pres, peak RSS).
With sx sy sz (divisors of N): the box is split into sx*sy*sz &MESH blocks (mesh order x fastest); FDS reads the UVW file from its start for every mesh, so the file holds one block and
the field is the Taylor-Green with the wavelength of a block (periodic per block, hence continuous across the block interfaces). Run with mpirun -np sx*sy*sz.
The csv of 128^3 is about 70 MB: generate it where it is used, do not commit it."""
import sys, os, numpy as np
n = int(sys.argv[1]); out = sys.argv[2]; t_end = sys.argv[3] if len(sys.argv) > 3 else '0.67'
sx, sy, sz = (int(a) for a in sys.argv[4:7]) if len(sys.argv) > 6 else (1, 1, 1)
L = 0.56549; os.makedirs(out, exist_ok=True)
assert n % sx == 0 and n % sy == 0 and n % sz == 0
nx_, ny_, nz_ = n // sx, n // sy, n // sz            # cells of one block
dx = L / n
X = ((np.arange(nx_) + 0.5) * dx)[None, None, :]; Y = ((np.arange(ny_) + 0.5) * dx)[None, :, None]; Z = ((np.arange(nz_) + 0.5) * dx)[:, None, None]
kx, ky, kz_ = 2 * np.pi / (L / sx), 2 * np.pi / (L / sy), 2 * np.pi / (L / sz); A = 0.3   # wavelength of one block
u = A * np.sin(kx * X) * np.cos(ky * Y) * np.cos(kz_ * Z)
v = -A * np.cos(kx * X) * np.sin(ky * Y) * np.cos(kz_ * Z)
w = np.zeros_like(u)
with open('%s/tg%d_uvw.csv' % (out, n), 'w') as f:   # C order of (z, y, x) = I fastest
    for kk in range(nz_):
        np.savetxt(f, np.stack([u[kk].ravel(), v[kk].ravel(), w[kk].ravel()], axis=1), fmt='%.8g', delimiter=',')
open('%s/tg%d.fds' % (out, n), 'w').write("""&HEAD CHID='tg%(n)d', TITLE='periodic Taylor-Green velocity, constant Smagorinsky, %(n)d^3 (driver size case)' /
%(meshes)s
&TIME T_END=%(t)s /
&MISC STRATIFICATION=.FALSE.
      PERIODIC_TEST=2
      UVW_FILE='tg%(n)d_uvw.csv'
      TURBULENCE_MODEL='CONSTANT SMAGORINSKY',C_SMAGORINSKY=0.2/
&DUMP RAMP_UVW='times' /
&RAMP ID='times', T=0.00 /
&RAMP ID='times', T=%(t)s /
&RADI RADIATION=.FALSE. /
&VENT PBX=0, SURF_ID='PERIODIC'/
&VENT PBX=%(L)g, SURF_ID='PERIODIC'/
&VENT PBY=0, SURF_ID='PERIODIC'/
&VENT PBY=%(L)g, SURF_ID='PERIODIC'/
&VENT PBZ=0, SURF_ID='PERIODIC'/
&VENT PBZ=%(L)g, SURF_ID='PERIODIC'/
&DEVC XB=0,%(L)g,0,%(L)g,0,%(L)g,QUANTITY='KINETIC ENERGY',ID='KE',SPATIAL_STATISTIC='MEAN' /
&TAIL /
""" % {'n': n, 'L': L, 't': t_end, 'meshes': '\n'.join(
    "&MESH IJK=%d,%d,%d, XB=%.8g,%.8g,%.8g,%.8g,%.8g,%.8g /" % (nx_, ny_, nz_, i * nx_ * dx, (i + 1) * nx_ * dx if i < sx - 1 else L, j * ny_ * dx, (j + 1) * ny_ * dx if j < sy - 1 else L, k * nz_ * dx, (k + 1) * nz_ * dx if k < sz - 1 else L)
    for k in range(sz) for j in range(sy) for i in range(sx))})
print('wrote %s/tg%d.fds and tg%d_uvw.csv' % (out, n, n))
