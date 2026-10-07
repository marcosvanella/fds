#!/usr/bin/env python3
"""Generate a 3-D vortex-ring cold-flow case for the owner-show adaptive TAG_VORT demo.

Writes:
  vortex_ring_<N>_l0.fds + vortex_ring_<N>_uvw.csv          (1 mesh)
  vortex_ring_<N>_4m.fds + vortex_ring_<N>_uvw.csv          (2x2x1, np=4)
       + vortex_ring_<N>_4m_uvw_m{1..4}.csv                 (sliced UVW for CSVF)
  vortex_ring_<N>_8m.fds + ..._8m_uvw_m{1..8}.csv           (2x2x2, np=8)
  optional: --sx --sy --sz for custom splits (NMESHES=sx*sy*sz)

FDS UVW_INIT (PERIODIC_TEST=2) reads each mesh from the *start* of UVW_FILE, so a
global N^3 CSV cannot tile multi-mesh. Multi-mesh cases therefore emit per-mesh
slices and &CSVF UVWFILE='..._uvw_m1.csv' (UVW_RESTART overwrites after the
PERIODIC_TEST=2 read). UVW_FILE is set to the m1 slice (mesh-sized) so the first
read does not EOF.

IC: Rosenhead-regularised Biot-Savart circular filament (x–y plane), self-advects +z.
FREEZE_VELOCITY off. No dye / TAG_SPEC — tagging is TAG_VORT on undivided |ω|.

usage:
  make_vortex_ring_case.py [N] [outdir]
  make_vortex_ring_case.py [N] [outdir] --4m
  make_vortex_ring_case.py [N] [outdir] --8m
  make_vortex_ring_case.py [N] [outdir] --sx 2 --sy 2 --sz 1
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

def ring_field(N: int, L: float = 1.0):
    dx = L / N
    R, a, Gamma = 0.20, 0.045, 0.8
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
    Z, Y, X = np.meshgrid(zc, yc, xc, indexing="ij")
    U = np.zeros_like(X)
    V = np.zeros_like(X)
    W = np.zeros_like(X)
    for s in range(nseg):
        rxv, ryv, rzv = X - rx[s], Y - ry[s], Z - rz[s]
        r2 = rxv * rxv + ryv * ryv + rzv * rzv + a * a
        r3 = r2 * np.sqrt(r2)
        fac = (Gamma / (4.0 * np.pi)) / r3
        U += fac * (dly[s] * rzv - dlz[s] * ryv)
        V += fac * (dlz[s] * rxv - dlx[s] * rzv)
        W += fac * (dlx[s] * ryv - dly[s] * rxv)
    omega_phys_peak = Gamma / (np.pi * a * a)
    omega_undiv_hint = omega_phys_peak * dx
    meta = dict(R=R, a=a, Gamma=Gamma, cx=cx, cy=cy, cz=cz,
                w=omega_undiv_hint, lo=0.25 * omega_undiv_hint, hi=0.5 * omega_undiv_hint,
                umax=float(np.max(np.abs(U))), wmax=float(np.max(np.abs(W))))
    return U, V, W, dx, meta

def write_uvw_block(path, U, V, W):
    with open(path, "w") as f:
        for kk in range(U.shape[0]):
            np.savetxt(f, np.stack([U[kk].ravel(), V[kk].ravel(), W[kk].ravel()], axis=1),
                       fmt="%.10g", delimiter=",")

def vents(L):
    return f"""&VENT PBX=0, SURF_ID='PERIODIC' /
&VENT PBX={L:g}, SURF_ID='PERIODIC' /
&VENT PBY=0, SURF_ID='PERIODIC' /
&VENT PBY={L:g}, SURF_ID='PERIODIC' /
&VENT PBZ=0, SURF_ID='PERIODIC' /
&VENT PBZ={L:g}, SURF_ID='PERIODIC' /
"""

def write_l0(out, N, L, U, V, W, meta):
    uvw = os.path.join(out, f"vortex_ring_{N}_uvw.csv")
    write_uvw_block(uvw, U, V, W)
    fds = os.path.join(out, f"vortex_ring_{N}_l0.fds")
    open(fds, "w").write(f"""! 3-D vortex ring, level-0 only (driver --two-level-run grows fine via TAG_VORT).
! IC: Rosenhead Biot-Savart ring (see make_vortex_ring_case.py). Self-advecting (no FREEZE_VELOCITY).
! UVW_FILE + PERIODIC_TEST=2 (same path as csmag_32 / make_periodic_case.py).
! Hint: undivided |ω| peak ~ {meta['w']:g}; try --tag-vort auto or --tag-vort in [{meta['lo']:.4g}, {meta['hi']:.4g}].
! Owner overnight: prefer --maxsize 8 --blocking 4, TAG_KEEP~0.85–0.9 (tight core), np=NMESHES.
! See Source/driver/notes/vortex-ring-owner-show.md.

&HEAD CHID='vortex_ring_{N}_l0', TITLE='3-D vortex ring cold-flow, {N}^3, owner-show TAG_VORT adaptive' /

&MESH IJK={N},{N},{N}, XB=0.0,{L:g},0.0,{L:g},0.0,{L:g} /

&TIME T_END=2.0 /

&MISC SIMULATION_MODE='DNS', STRATIFICATION=.FALSE., NOISE=.FALSE., CFL_MAX=0.25,
      PERIODIC_TEST=2, UVW_FILE='vortex_ring_{N}_uvw.csv' /

&PRES VELOCITY_TOLERANCE=1.E-6, MAX_PRESSURE_ITERATIONS=100 /

&RADI RADIATION=.FALSE. /

&SPEC ID='BACKGROUND', BACKGROUND=.TRUE., MW=28.0, VISCOSITY=1.0E-4, DIFFUSIVITY=1.0E-4, SPECIFIC_HEAT=1.0 /

{vents(L)}
&TAIL /
""")
    print("wrote", fds)
    print("wrote", uvw)
    return fds, uvw

def write_multimesh(out, N, L, U, V, W, meta, sx, sy, sz, tag):
    """tag e.g. '4m' or '8m' or f'{sx}x{sy}x{sz}'."""
    assert N % sx == 0 and N % sy == 0 and N % sz == 0
    nx, ny, nz = N // sx, N // sy, N // sz
    dx = L / N
    # Always also write global UVW (single-mesh companion / option-b source).
    global_uvw = os.path.join(out, f"vortex_ring_{N}_uvw.csv")
    if not os.path.isfile(global_uvw):
        write_uvw_block(global_uvw, U, V, W)
        print("wrote", global_uvw)

    mesh_lines = []
    nm = 0
    stem = f"vortex_ring_{N}_{tag}"
    for k in range(sz):
        for j in range(sy):
            for i in range(sx):
                nm += 1
                i0, i1 = i * nx, (i + 1) * nx
                j0, j1 = j * ny, (j + 1) * ny
                k0, k1 = k * nz, (k + 1) * nz
                xb0, xb1 = i0 * dx, (i1 * dx if i < sx - 1 else L)
                yb0, yb1 = j0 * dx, (j1 * dx if j < sy - 1 else L)
                zb0, zb1 = k0 * dx, (k1 * dx if k < sz - 1 else L)
                mesh_lines.append(
                    f"&MESH IJK={nx},{ny},{nz}, XB={xb0:.8g},{xb1:.8g},{yb0:.8g},{yb1:.8g},{zb0:.8g},{zb1:.8g}, MPI_PROCESS={nm-1} /"
                )
                Us, Vs, Ws = U[k0:k1, j0:j1, i0:i1], V[k0:k1, j0:j1, i0:i1], W[k0:k1, j0:j1, i0:i1]
                mpath = os.path.join(out, f"{stem}_uvw_m{nm}.csv")
                write_uvw_block(mpath, Us, Vs, Ws)
                print("wrote", mpath, "cells", Us.size)

    nmesh = sx * sy * sz
    fds = os.path.join(out, f"{stem}.fds")
    # UVW_FILE = m1 slice (mesh-sized) so PERIODIC_TEST=2 first read does not EOF;
    # &CSVF then UVW_RESTART-overwrites each mesh from its own slice.
    open(fds, "w").write(f"""! 3-D vortex ring, {N}^3 as {nmesh} meshes ({sx}x{sy}x{sz} of {nx}x{ny}x{nz}) for np={nmesh}.
! IC: Rosenhead ring; per-mesh UVW slices + &CSVF (UVW_RESTART). PERIODIC_TEST=2 kept for face copy.
! FDS UVW_INIT reads each mesh from the start of UVW_FILE — global N^3 CSV does NOT tile; use slices.
! Owner overnight: --two-level-run --tag-vort auto --maxsize 8 --blocking 4 (tight L1); TAG_KEEP~0.85.
! Hint: undivided |ω| peak ~ {meta['w']:g}; try --tag-vort auto or [{meta['lo']:.4g}, {meta['hi']:.4g}] (raise thr for tighter).
! See Source/driver/notes/vortex-ring-owner-show.md.

&HEAD CHID='{stem}', TITLE='3-D vortex ring cold-flow, {N}^3 {tag}, owner-show TAG_VORT adaptive' /

{chr(10).join(mesh_lines)}

&TIME T_END=2.0 /

&MISC SIMULATION_MODE='DNS', STRATIFICATION=.FALSE., NOISE=.FALSE., CFL_MAX=0.25,
      PERIODIC_TEST=2, UVW_FILE='{stem}_uvw_m1.csv' /

&CSVF UVWFILE='{stem}_uvw_m1.csv' /

&PRES VELOCITY_TOLERANCE=1.E-6, MAX_PRESSURE_ITERATIONS=100 /

&RADI RADIATION=.FALSE. /

&SPEC ID='BACKGROUND', BACKGROUND=.TRUE., MW=28.0, VISCOSITY=1.0E-4, DIFFUSIVITY=1.0E-4, SPECIFIC_HEAT=1.0 /

{vents(L)}
&TAIL /
""")
    print("wrote", fds, f"(NMESHES={nmesh})")
    return fds

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("N", nargs="?", type=int, default=32)
    ap.add_argument("outdir", nargs="?", default=None)
    ap.add_argument("--4m", dest="do4", action="store_true", help="also/emit 2x2x1 (np=4)")
    ap.add_argument("--8m", dest="do8", action="store_true", help="also/emit 2x2x2 (np=8)")
    ap.add_argument("--sx", type=int, default=None)
    ap.add_argument("--sy", type=int, default=None)
    ap.add_argument("--sz", type=int, default=None)
    ap.add_argument("--l0-only", action="store_true", help="only single-mesh l0 (default if no multi flags)")
    args = ap.parse_args()
    N = args.N
    out = args.outdir or os.path.join(os.path.dirname(__file__), "cases")
    os.makedirs(out, exist_ok=True)
    L = 1.0
    U, V, W, dx, meta = ring_field(N, L)

    multi = args.do4 or args.do8 or (args.sx is not None)
    if not multi or args.l0_only or not multi:
        # Always write l0 unless user asked only for a custom split without l0 — keep l0 always.
        write_l0(out, N, L, U, V, W, meta)
    if args.do4:
        write_multimesh(out, N, L, U, V, W, meta, 2, 2, 1, "4m")
    if args.do8:
        write_multimesh(out, N, L, U, V, W, meta, 2, 2, 2, "8m")
    if args.sx is not None:
        sx, sy, sz = args.sx, args.sy or 1, args.sz or 1
        tag = f"{sx}x{sy}x{sz}"
        write_multimesh(out, N, L, U, V, W, meta, sx, sy, sz, tag)

    print(
        "IC: R=%.3g a=%.3g Gamma=%.3g center=(%.2f,%.2f,%.2f) |U|max=%.4g |W|max=%.4g undiv_|ω|_hint=%.4g thr~[%.4g,%.4g]"
        % (meta["R"], meta["a"], meta["Gamma"], meta["cx"], meta["cy"], meta["cz"],
           meta["umax"], meta["wmax"], meta["w"], meta["lo"], meta["hi"])
    )

if __name__ == "__main__":
    main()
