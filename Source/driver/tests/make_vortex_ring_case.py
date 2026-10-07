#!/usr/bin/env python3
"""Generate a 3-D vortex-ring cold-flow case for the owner-show adaptive TAG_VORT demo.

Writes <outdir>/vortex_ring_<N>_l0.fds and vortex_ring_<N>_uvw.csv
(UVW_INIT / PERIODIC_TEST=2: N^3 lines u,v,w, I fastest, then J, then K).

IC: Rosenhead-regularised Biot-Savart circular filament in the x-y plane
(radius R, core a, circulation Gamma), centred near mid-box; self-advects in +z.
FREEZE_VELOCITY is off. No dye / TAG_SPEC — tagging is TAG_VORT on undivided |ω|.

usage: make_vortex_ring_case.py [N] [outdir]
"""
from __future__ import annotations
import os, sys
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 32
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "cases")
os.makedirs(out, exist_ok=True)

L = 1.0
dx = L / N
R = 0.20
a = 0.045
Gamma = 0.8
cx, cy, cz = 0.5, 0.5, 0.40
nseg = 256
theta = np.linspace(0.0, 2.0 * np.pi, nseg, endpoint=False)
dth = 2.0 * np.pi / nseg
rx = cx + R * np.cos(theta)
ry = cy + R * np.sin(theta)
rz = np.full(nseg, cz)
dlx = -R * np.sin(theta) * dth
dly = R * np.cos(theta) * dth
dlz = np.zeros(nseg)

xc = (np.arange(N) + 0.5) * dx
yc = (np.arange(N) + 0.5) * dx
zc = (np.arange(N) + 0.5) * dx
# (K, J, I) so I is fastest when ravelled C-order within each K slab written below
Z, Y, X = np.meshgrid(zc, yc, xc, indexing="ij")
U = np.zeros_like(X)
V = np.zeros_like(X)
W = np.zeros_like(X)
for s in range(nseg):
    rxv = X - rx[s]
    ryv = Y - ry[s]
    rzv = Z - rz[s]
    r2 = rxv * rxv + ryv * ryv + rzv * rzv + a * a
    r3 = r2 * np.sqrt(r2)
    cxv = dly[s] * rzv - dlz[s] * ryv
    cyv = dlz[s] * rxv - dlx[s] * rzv
    czv = dlx[s] * ryv - dly[s] * rxv
    fac = (Gamma / (4.0 * np.pi)) / r3
    U += fac * cxv
    V += fac * cyv
    W += fac * czv

# Undivided |ω| estimate on collocated samples (for the note / threshold hint only).
# Physical ω_φ peak ~ Gamma/(pi a^2); undivided scale ~ that * dx.
omega_phys_peak = Gamma / (np.pi * a * a)
omega_undiv_hint = omega_phys_peak * dx
thr_lo, thr_hi = 0.25 * omega_undiv_hint, 0.5 * omega_undiv_hint

uvw = os.path.join(out, "vortex_ring_%d_uvw.csv" % N)
with open(uvw, "w") as f:
    for kk in range(N):
        np.savetxt(
            f,
            np.stack([U[kk].ravel(), V[kk].ravel(), W[kk].ravel()], axis=1),
            fmt="%.10g",
            delimiter=",",
        )

fds = os.path.join(out, "vortex_ring_%d_l0.fds" % N)
open(fds, "w").write(
    """! 3-D vortex ring, level-0 only (driver --two-level-run grows fine via TAG_VORT).
! IC: Rosenhead Biot-Savart ring (see make_vortex_ring_case.py). Self-advecting (no FREEZE_VELOCITY).
! UVW_FILE + PERIODIC_TEST=2 (same path as csmag_32 / make_periodic_case.py).
! Hint: undivided |ω| peak ~ %(w)g; try --tag-vort auto or --tag-vort in [%(lo).4g, %(hi).4g].
! See Source/driver/notes/vortex-ring-owner-show.md.

&HEAD CHID='vortex_ring_%(n)d_l0', TITLE='3-D vortex ring cold-flow, %(n)d^3, owner-show TAG_VORT adaptive' /

&MESH IJK=%(n)d,%(n)d,%(n)d, XB=0.0,%(L)g,0.0,%(L)g,0.0,%(L)g /

&TIME T_END=2.0 /

&MISC SIMULATION_MODE='DNS', STRATIFICATION=.FALSE., NOISE=.FALSE., CFL_MAX=0.25,
      PERIODIC_TEST=2, UVW_FILE='vortex_ring_%(n)d_uvw.csv' /

&PRES VELOCITY_TOLERANCE=1.E-6, MAX_PRESSURE_ITERATIONS=100 /

&RADI RADIATION=.FALSE. /

&SPEC ID='BACKGROUND', BACKGROUND=.TRUE., MW=28.0, VISCOSITY=1.0E-4, DIFFUSIVITY=1.0E-4, SPECIFIC_HEAT=1.0 /

&VENT PBX=0, SURF_ID='PERIODIC' /
&VENT PBX=%(L)g, SURF_ID='PERIODIC' /
&VENT PBY=0, SURF_ID='PERIODIC' /
&VENT PBY=%(L)g, SURF_ID='PERIODIC' /
&VENT PBZ=0, SURF_ID='PERIODIC' /
&VENT PBZ=%(L)g, SURF_ID='PERIODIC' /

&TAIL /
"""
    % {"n": N, "L": L, "w": omega_undiv_hint, "lo": thr_lo, "hi": thr_hi}
)
print("wrote", fds)
print("wrote", uvw)
print(
    "IC: R=%.3g a=%.3g Gamma=%.3g center=(%.2f,%.2f,%.2f) |U|max=%.4g |W|max=%.4g undiv_|ω|_hint=%.4g thr~[%.4g,%.4g]"
    % (R, a, Gamma, cx, cy, cz, float(np.max(np.abs(U))), float(np.max(np.abs(W))), omega_undiv_hint, thr_lo, thr_hi)
)
