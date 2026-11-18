# geom_exact_flows_error.py
# Grid convergence of the cut-cell momentum on exact flow solutions (&MISC PERIODIC_TEST=31 to 37), for the per-piece
# scheme (decks geom_<case>_N<n>.fds, &MISC CC_O3E=.TRUE.) and the reference scheme (decks geom_<case>_N<n>_ref.fds).
#
# FDS writes <CHID>_mms_m<mesh>_s<step>.csv at MMS_TIMER, one row per gas velocity face and per gas piece of a cut face:
#   MESH,AXIS,CLASS,X,Y,Z,AREA,VEL,VEXACT      (CLASS 0 whole Cartesian face, 1 piece of a cut face)
# The error is the area-weighted L1 norm of (VEL-VEXACT)/U_wall over all rows, U_wall the wall speed of the case.

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

datadir = '../../Verification/Complex_Geometry/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

# name: (PERIODIC_TEST, length such that the grid spacing is length/N, wall speed U_wall (m/s), grids N)
# N is the number of cells across the domain, for the channel decks (34, 36) the number of cells per metre.
cases = {
   'rigid_rotation_3d':           (31, 2.1, 0.1, [21, 42]),
   'couette_3d':                  (32, 2.1, 0.1, [21, 42]),
   'couette_2d':                  (33, 3.0, 0.2, [24, 48, 96]),
   'poiseuille_2d_30deg':         (34, 1.0, 1.0, [4, 8, 16]),
   'uniform_flow_2d':             (35, 3.0, 0.2, [24, 48, 96]),
   'blowing_channel_2d_30deg':    (36, 1.0, 1.0, [4, 8, 16]),
   'oscillating_rotation_2d':     (37, 3.0, 0.2, [24, 48, 96]),
}

# the reference scheme is drawn first, with the larger markers, so that the two curves stay visible where they coincide
variants = [('_ref', 'FDS (reference scheme)', 'rs-', 8), ('', 'FDS', 'bo--', 4)]


def l1_error(chid, uwall):
   """Area-weighted L1 of (VEL-VEXACT)/uwall over the faces of the CSV file(s) of the run chid, or None."""
   files = sorted(glob.glob(os.path.join(datadir, chid + '_mms_m*_s*.csv')))
   if not files:
      print(f"[WARN] No {chid}_mms_m*_s*.csv file in {datadir}. Skipping.")
      return None
   rows = np.vstack([np.genfromtxt(f, delimiter=',', skip_header=2) for f in files])
   area, vel, vex = rows[:, 6], rows[:, 7], rows[:, 8]
   return np.sum(area * np.abs(vel - vex)) / np.sum(area) / uwall


for name, (ptest, length, uwall, grids) in cases.items():

   dx = {}
   err = {}
   for suffix, label, style, msize in variants:
      dx[suffix], err[suffix] = [], []
      for n in grids:
         e = l1_error(f'geom_{name}_N{n}{suffix}', uwall)
         if e is not None:
            dx[suffix].append(length / n)
            err[suffix].append(e)
         print(f"  {name:26s} PERIODIC_TEST={ptest} N={n:3d}{suffix:5s} L1/U_wall = {e}")
   if not err[''] and not err['_ref']:
      print(f"[WARN] No data for case {name}. No figure.")
      continue

   allerr = np.array(err[''] + err['_ref'])
   alldx = np.array(dx[''] + dx['_ref'])
   decade = lambda v, f: 10.**f(np.log10(v))
   axisval = [alldx.min() / 1.5, alldx.max() * 1.5,
              decade(allerr.min() / 4., np.floor), decade(allerr.max() * 4., np.ceil)]

   try:
      version_string = fdsplotlib.get_version_string(datadir + f'geom_{name}_N{grids[-1]}_ref_git.txt')
   except OSError:
      version_string = ''

   fig = None
   for suffix, label, style, msize in variants:
      if not err[suffix]:
         continue
      kw = dict(x_data=np.array(dx[suffix]), y_data=np.array(err[suffix]), marker_style=style, markersize=msize,
                data_label=label, figure_handle=fig)
      if fig is None:
         kw.update(plot_type='loglog', revision_label=version_string, x_min=axisval[0], x_max=axisval[1],
                   y_min=axisval[2], y_max=axisval[3], x_label='$\\Delta x$ (m)', y_label='$L_1$ error of $u_n$ / $U_{wall}$')
      fig = fdsplotlib.plot_to_fig(**kw)

   # first and second order slopes, started below the largest error of the reference scheme, which has every grid
   dxr = np.array(dx['_ref'])
   e0 = 0.3 * max(err['_ref'])
   fdsplotlib.plot_to_fig(x_data=dxr, y_data=e0 * dxr / dxr[0], marker_style='k--', figure_handle=fig, data_label='$O(\\Delta x)$')
   fdsplotlib.plot_to_fig(x_data=dxr, y_data=e0 * (dxr / dxr[0])**2, marker_style='k-', figure_handle=fig,
                          data_label='$O(\\Delta x^2)$')

   plt.savefig(plotdir + f'geom_{name}_convergence.pdf', format='pdf')
   plt.close()
