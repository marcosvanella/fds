# geom_hvac_supply_jet.py
# Supply box of an HVAC system turned against the grid (Verification/Complex_Geometry/geom_hvac_supply_*.fds).
#
# A box blows 5 m3/s of air through one face into a room. The box is turned by 15, 30 and 45 degrees about the vertical axis and
# by 60 degrees about the axis (-1,1,1). The jet should leave along the normal n of the supply face, with the momentum
# flux of the aligned box. Line devices (four rows, 0.2 m apart, along the face's short side) record the mean velocity (t > 10 s)
# on planes normal to n at 0.5, 1 and 2 m from the center of the face. On each plane
#    J_n = rho sum(u_n^2) dA   momentum flux along n, given as J_n/J_0 with J_0 = rho U_0^2 A_0 of the supplied air,
#    J_t = rho sum(u_n u_t) dA  momentum flux along a direction t in the plane of the face,
# and the angle between the jet and n is atan(sqrt(J_t1^2 + J_t2^2)/J_n).

import os
import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

datadir = '../../Verification/Complex_Geometry/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

RHO, U0, A0 = 1.2, 5.0/0.96, 0.96
ROOM = np.array((10., 10., 6.))
C = np.array((7.1, 5.0, 2.0))          # center of the box
HALF_LENGTH = 1.5                      # half length of the box along its own x axis
DIST = (0.5, 1.0, 2.0)
ROWS = (-0.3, -0.1, 0.1, 0.3)
S1 = np.linspace(-2.5, 2.5, 51)


def rotation(axis, deg):
   a = np.array(axis, float)
   a /= np.linalg.norm(a)
   t = np.radians(deg)
   K = np.array(((0, -a[2], a[1]), (a[2], 0, -a[0]), (-a[1], a[0], 0)))
   return np.eye(3) + np.sin(t)*K + (1. - np.cos(t))*K@K


# deck: (label, turn angle (deg), rotation matrix, kind); kind: obst, geom or stair
decks = {
   'geom_hvac_supply_obst':        ('OBST, aligned', 0, rotation((0, 0, 1), 0), 'obst'),
   'geom_hvac_supply_0deg_ref':    ('GEOM, 0 deg', 0, rotation((0, 0, 1), 0), 'geom'),
   'geom_hvac_supply_15deg_ref':   ('GEOM, 15 deg', 15, rotation((0, 0, 1), 15), 'geom'),
   'geom_hvac_supply_30deg_ref':   ('GEOM, 30 deg', 30, rotation((0, 0, 1), 30), 'geom'),
   'geom_hvac_supply_30deg_stair': ('stair-stepped OBST, 30 deg', 30, rotation((0, 0, 1), 30), 'stair'),
   'geom_hvac_supply_45deg_ref':   ('GEOM, 45 deg', 45, rotation((0, 0, 1), 45), 'geom'),
   'geom_hvac_supply_3d_ref':      ('GEOM, 60 deg about (-1,1,1)', 60, rotation((-1, 1, 1), 60), 'geom'),
   'geom_hvac_supply_3d_stair':    ('stair-stepped OBST, 60 deg about (-1,1,1)', 60, rotation((-1, 1, 1), 60), 'stair'),
}


def read_line_file(path):
   """Columns of a FDS _line.csv file: a units row, a names row, then the data rows."""
   with open(path) as f:
      lines = f.read().split('\n')
   names = [x.strip() for x in lines[1].split(',')]
   cols = {n: [] for n in names}
   for ln in lines[2:]:
      for n, v in zip(names, ln.split(',')):
         if v.strip() and v.strip() != 'NaN':
            cols[n].append(float(v))
   return {n: np.array(v) for n, v in cols.items()}


def jet(chid, R):
   """J_n/J_0 and the angle (deg) between the jet and the face normal on each plane, or None if the run is missing."""
   path = datadir + chid + '_line.csv'
   if not os.path.exists(path):
      print(f"[WARN] {path} not found. Skipping.")
      return None
   T = read_line_file(path)
   n, t1, t2 = R @ (-1., 0., 0.), R @ (0., 1., 0.), R @ (0., 0., 1.)
   fc = C + R @ (-HALF_LENGTH, 0., 0.)
   jn, ang = [], []
   for d in DIST:
      Jn = Jt1 = Jt2 = 0.
      for ir, s2 in enumerate(ROWS):
         key = f'_d{int(round(d*100)):03d}_r{ir}'
         u, v, w = T['U' + key], T['V' + key], T['W' + key]
         un = u*n[0] + v*n[1] + w*n[2]
         ut1 = u*t1[0] + v*t1[1] + w*t1[2]
         ut2 = u*t2[0] + v*t2[1] + w*t2[2]
         dA = 0.1*0.2                  # point spacing along the line times row spacing
         Jn += RHO*np.sum(un*un)*dA
         Jt1 += RHO*np.sum(un*ut1)*dA
         Jt2 += RHO*np.sum(un*ut2)*dA
      jn.append(Jn/(RHO*U0*U0*A0))
      ang.append(np.degrees(np.arctan2(np.hypot(Jt1, Jt2), Jn)))
   return np.array(jn), np.array(ang)


res = {}
for chid, (label, turn, R, kind) in decks.items():
   r = jet(chid, R)
   if r is not None:
      res[chid] = r
      print(f"  {chid:32s} J/J0 at 0.5, 1, 2 m: {r[0][0]:.3f} {r[0][1]:.3f} {r[0][2]:.3f}   angle off the normal (deg): "
            f"{r[1][0]:.1f} {r[1][1]:.1f} {r[1][2]:.1f}")
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

# Figure 1: angle between the jet and the face normal (mean over the three planes) against the turn angle.
series = [('OBST, aligned', ['geom_hvac_supply_obst'], 'ko', None),
          ('GEOM, turned about vertical axis', ['geom_hvac_supply_0deg_ref', 'geom_hvac_supply_15deg_ref', 'geom_hvac_supply_30deg_ref',
                                                'geom_hvac_supply_45deg_ref'], 'rs-', None),
          ('GEOM, 60 deg about (-1,1,1)', ['geom_hvac_supply_3d_ref'], 'rs', None),
          ('stair-stepped OBST', ['geom_hvac_supply_30deg_stair', 'geom_hvac_supply_3d_stair'], 'bo--', None)]
fig = None
for label, names, style, _ in series:
   names = [c for c in names if c in res]
   if not names:
      continue
   kw = dict(x_data=np.array([decks[c][1] for c in names]), y_data=np.array([np.mean(res[c][1]) for c in names]), marker_style=style,
             data_label=label, figure_handle=fig)
   if fig is None:
      kw.update(revision_label=version_string, x_min=-5, x_max=65, y_min=0, y_max=45,
                x_label='Turn angle of the box (deg)', y_label='Angle between jet and face normal (deg)')
   fig = fdsplotlib.plot_to_fig(**kw)
plt.savefig(plotdir + 'geom_hvac_supply_angle.pdf', format='pdf')
plt.close()

# Figure 2: momentum flux of the jet along the face normal against the distance from the face.
styles = {'geom_hvac_supply_obst': 'ko-', 'geom_hvac_supply_0deg_ref': 'rs-', 'geom_hvac_supply_15deg_ref': 'gs-',
          'geom_hvac_supply_30deg_ref': 'bs-', 'geom_hvac_supply_45deg_ref': 'ms-', 'geom_hvac_supply_3d_ref': 'cs-',
          'geom_hvac_supply_30deg_stair': 'bo--', 'geom_hvac_supply_3d_stair': 'co--'}
fig = None
for chid, (label, turn, R, kind) in decks.items():
   if chid not in res:
      continue
   kw = dict(x_data=np.array(DIST), y_data=res[chid][0], marker_style=styles[chid], data_label=label, figure_handle=fig)
   if fig is None:
      kw.update(revision_label=version_string, x_min=0.25, x_max=2.25, y_min=0, y_max=1,
                x_label='Distance from the face (m)', y_label='Momentum flux $J_n/J_0$', legend_location='lower left')
   fig = fdsplotlib.plot_to_fig(**kw)
plt.savefig(plotdir + 'geom_hvac_supply_momentum.pdf', format='pdf')
plt.close()
