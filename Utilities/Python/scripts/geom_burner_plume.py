# geom_burner_plume.py
# Burner turned against the grid (Verification/Complex_Geometry/geom_burner_*.fds).
#
# A 0.6 m x 0.6 m x 0.1 m burner of 360 kW burns propane in open air. The burner is aligned with the grid (OBST and GEOM) or turned
# by 60 degrees about the axis (-1,1,1) (GEOM, and a stair-stepped OBST version), with gravity turned with it, so that the plume
# rises along a = (1,2,2)/3. Measures along the axis through the center of the burner:
#   - the mean temperature and the mean velocity along the axis, a.(u,v,w), at 18 points (line devices, t > 5 s),
#   - the flame height: the distance above the burner base below which 95 % of the heat release of the time-mean 3D slice of
#     HRRPUV lies, and the heat release of the slice (the burner releases 360 kW).
# The run is turbulent and one seed is used, so differences between decks have to be compared with the seed-to-seed spread
# (see the Verification Guide).

import os
import struct
import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

datadir = '../../Verification/Complex_Geometry/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

TSTAT = 5.
ZP = 0.125 + 1.85/17*np.arange(18)           # distance of the 18 points above the burner base (m)
RT = np.array(((2, -2, 1), (1, 2, 2), (-2, -1, 2)))/3.      # rotation of 60 degrees about (-1,1,1)

# deck: (label, axis of the plume, marker style)
decks = {
   'geom_burner_obst':      ('OBST, aligned', np.array((0., 0., 1.)), 'ko-'),
   'geom_burner_0deg_ref':  ('GEOM, aligned', np.array((0., 0., 1.)), 'rs-'),
   'geom_burner_3d_ref':    ('GEOM, 60 deg about (-1,1,1)', RT[:, 2], 'bs--'),
   'geom_burner_3d_stair':  ('stair-stepped OBST, 60 deg about (-1,1,1)', RT[:, 2], 'go--'),
}


def read_line_file(path):
   with open(path) as f:
      lines = f.read().split('\n')
   names = [x.strip() for x in lines[1].split(',')]
   cols = {n: [] for n in names}
   for ln in lines[2:]:
      for n, v in zip(names, ln.split(',')):
         if v.strip() and v.strip() != 'NaN':
            cols[n].append(float(v))
   return {n: np.array(v) for n, v in cols.items()}


def fortran_records(path):
   """Records of an unformatted FDS file: a 4 byte length, the data, the 4 byte length."""
   b = open(path, 'rb').read()
   p, out = 0, []
   while p + 4 <= len(b):
      n = struct.unpack('<i', b[p:p + 4])[0]
      if p + 8 + n > len(b):
         break
      out.append(b[p + 4:p + 4 + n])
      p += 8 + n
   return out


def mesh_grids(smv):
   """(IBAR, JBAR, KBAR, XB) of every mesh of the .smv file."""
   t = open(smv, errors='replace').read().split('\n')
   out = []
   for i, ln in enumerate(t):
      if ln.startswith('GRID'):
         ib, jb, kb = (int(v) for v in t[i + 1].split()[:3])
         j = i + 1
         while not t[j].startswith('PDIM'):
            j += 1
         out.append((ib, jb, kb, [float(v) for v in t[j + 1].split()[:6]]))
   return out


def flame(chid, axis):
   """Flame height (m) and heat release (kW) from the time mean (t > TSTAT) of the 3D HRRPUV slice."""
   s_all, q_all = [], []
   for m, (ib, jb, kb, xb) in enumerate(mesh_grids(datadir + chid + '.smv'), 1):
      f = datadir + f'{chid}_{m}_1.sf'
      if not os.path.exists(f):
         continue
      rec = fortran_records(f)
      i1, i2, j1, j2, k1, k2 = struct.unpack('<6i', rec[3])
      nx, ny, nz = i2 - i1 + 1, j2 - j1 + 1, k2 - k1 + 1
      acc, nfr = 0., 0
      for t_, dat in zip(rec[4::2], rec[5::2]):
         if struct.unpack('<f', t_)[0] >= TSTAT - 1e-6 and len(dat) == 4*nx*ny*nz:
            acc = acc + np.frombuffer(dat, dtype='<f4').reshape(nz, ny, nx).astype(float)
            nfr += 1
      if nfr == 0:
         return None
      acc /= nfr
      dx, dy, dz = (xb[1] - xb[0])/ib, (xb[3] - xb[2])/jb, (xb[5] - xb[4])/kb
      ii, jj, kk = np.arange(i1, i2 + 1), np.arange(j1, j2 + 1), np.arange(k1, k2 + 1)
      xc, yc, zc = xb[0] + (ii - 0.5)*dx, xb[2] + (jj - 0.5)*dy, xb[4] + (kk - 0.5)*dz
      inside = (((ii >= 1) & (ii <= ib))[None, None, :] & ((jj >= 1) & (jj <= jb))[None, :, None]
                & ((kk >= 1) & (kk <= kb))[:, None, None])
      Z, Y, X = np.meshgrid(zc, yc, xc, indexing='ij')
      s = X*axis[0] + Y*axis[1] + Z*axis[2] + 0.05          # distance above the burner base plane
      s_all.append(s[inside])
      q_all.append((acc*dx*dy*dz)[inside])
   if not s_all:
      return None
   s, q = np.concatenate(s_all), np.concatenate(q_all)
   edges = np.arange(-1.5, 4.0, 0.025)
   h, _ = np.histogram(s, bins=edges, weights=q)
   cum = np.cumsum(h)
   return float(np.interp(0.95, cum/cum[-1], edges[1:])), float(q.sum())


res = {}
for chid, (label, axis, style) in decks.items():
   path = datadir + chid + '_line.csv'
   if not os.path.exists(path):
      print(f"[WARN] {path} not found. Skipping.")
      continue
   L = read_line_file(path)
   T = L['T']
   Wa = axis[0]*L['U'] + axis[1]*L['V'] + axis[2]*L['W']
   fl = flame(chid, axis)
   res[chid] = (T, Wa, fl)
   print(f"  {chid:24s} T on the axis {T.min():.0f} to {T.max():.0f} C, flame height and heat release: {fl}")
if not res:
   print("[WARN] No data. No figures.")
   raise SystemExit(0)

version_string = ''
for chid in decks:
   try:
      version_string = fdsplotlib.get_version_string(datadir + chid + '_git.txt')
      break
   except OSError:
      pass

for k, (fname, ylabel, ymax) in enumerate((('geom_burner_temperature', 'Mean temperature on the axis ($^\\circ$C)', 1500),
                                            ('geom_burner_velocity', 'Mean velocity along the axis (m/s)', 11))):
   fig = None
   for chid, (label, axis, style) in decks.items():
      if chid not in res:
         continue
      kw = dict(x_data=ZP, y_data=res[chid][k], marker_style=style, data_label=label, figure_handle=fig)
      if fig is None:
         kw.update(revision_label=version_string, x_min=0, x_max=2, y_min=0, y_max=ymax,
                   x_label='Distance above the burner (m)', y_label=ylabel, legend_location='upper left')
      fig = fdsplotlib.plot_to_fig(**kw)
   plt.savefig(plotdir + fname + '.pdf', format='pdf')
   plt.close()

# Flame height and heat release in the slice, one point per deck
names = [c for c in decks if c in res and res[c][2] is not None]
if names:
   for fname, idx, ylabel, ymax in (('geom_burner_flame_height', 0, 'Flame height (m)', 4), ('geom_burner_heat_release', 1, 'Heat release in the slice (kW)', 600)):
      fig = None
      for k, c in enumerate(names):
         kw = dict(x_data=np.array([k + 1.]), y_data=np.array([res[c][2][idx]]), marker_style=decks[c][2][:2], data_label=decks[c][0],
                   figure_handle=fig)
         if fig is None:
            kw.update(revision_label=version_string, x_min=0.5, x_max=len(names) + 0.5, y_min=0, y_max=ymax, y_label=ylabel,
                      xticks=list(range(1, len(names) + 1)), legend_location='upper left')
         fig = fdsplotlib.plot_to_fig(**kw)
      fig.axes[0].set_xticklabels(['OBST', 'GEOM', 'GEOM 3D', 'stair 3D'][:len(names)])
      plt.savefig(plotdir + fname + '.pdf', format='pdf')
      plt.close()
