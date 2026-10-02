#!/usr/bin/env python3
"""Invariants of the periodic-only driver fixes, read from a FDSTL_STAGEG dump (raw_p1_prevflux_*_b0.bin, one box, fully periodic domain).
usage: check_periodic_raw.py <run-dir> <check>   check = face_match | mu_edge_corner | kres_edge_corner
  face_match       periodic domain-face velocity match (FDS MATCH_VELOCITY): U(0,j,k) == U(IBAR,j,k), V(i,0,k) == V(i,JBAR,k), W(i,j,0) == W(i,j,KBAR) to FACE_TOL (the dump point of a later cycle
                   is after a corrector update, which leaves 1e-16 round-off between the two faces; an unmatched run differs by 1e-3)
  mu_edge_corner   MU in the domain edge and corner cells is what COMPUTE_VISCOSITY of FDS ends with (clamped copy of the adjacent gas cell, corners from the interior corner cell)
  kres_edge_corner the same for KRES
exit status 0 = invariant holds, 1 = violated (the number of violating cells is printed)"""
import sys, glob, numpy as np
FACE_TOL = 1e-12

def load(d, name, b):
    p = '%s/raw_p1_prevflux_%s_b%d.bin' % (d, name, b)
    h = np.fromfile(p, dtype=np.int32, count=16)
    lo = h[0:3]; hi = h[3:6]; v0 = h[8:11]; v1 = h[11:14]
    shape = (hi - lo + 1)[::-1]
    a = np.fromfile(p, dtype=np.float64, offset=64).reshape(shape)  # a[k, j, i], AMReX index = lo + position
    return a, lo, v0, v1

def nboxes(d):
    return len(glob.glob('%s/raw_p1_prevflux_MU_b*.bin' % d))

def edges(A, XL, XH, YL, YH, ZL, ZH):
    """the edge/corner statements of FDS COMPUTE_VISCOSITY (as in fds_ghost_bc.f90, FDS_G_MU_EDGES_DOM); a statement applies when both sides it names are domain sides
    (the domain sides of this box are XL..ZH). A has FDS indices 0..IBP1 (position == FDS index with one ghost layer)."""
    B = A.copy()
    nz, ny, nx = B.shape; I1, J1, K1 = nx - 2, ny - 2, nz - 2; Ip, Jp, Kp = I1 + 1, J1 + 1, K1 + 1
    if XL and ZL: B[0, :, 0] = B[1, :, 1]
    if XH and ZL: B[0, :, Ip] = B[1, :, I1]
    if XH and ZH: B[Kp, :, Ip] = B[K1, :, I1]
    if XL and ZH: B[Kp, :, 0] = B[K1, :, 1]
    if YL and ZL: B[0, 0, :] = B[1, 1, :]
    if YH and ZL: B[0, Jp, :] = B[1, J1, :]
    if YH and ZH: B[Kp, Jp, :] = B[K1, J1, :]
    if YL and ZH: B[Kp, 0, :] = B[K1, 1, :]
    if XL and YL: B[:, 0, 0] = B[:, 1, 1]
    if XH and YL: B[:, 0, Ip] = B[:, 1, I1]
    if XH and YH: B[:, Jp, Ip] = B[:, J1, I1]
    if XL and YH: B[:, Jp, 0] = B[:, J1, 1]
    for kz, kk, zs in ((0, 1, ZL), (Kp, K1, ZH)):
        for jy, jj, ys in ((0, 1, YL), (Jp, J1, YH)):
            for ix, ii, xs in ((0, 1, XL), (Ip, I1, XH)):
                if xs and ys and zs: B[kz, jy, ix] = B[kk, jj, ii]
    return B

def main(d, chk):
    nb_ = nboxes(d)
    if chk == 'face_match':
        bad = 0
        for n, ax in (('U', 0), ('V', 1), ('W', 2)):   # ax: direction of the face normal (x,y,z)
            boxes = [load(d, n, b) for b in range(nb_)]
            N = max(int(v1[ax]) for _, _, _, v1 in boxes)          # the high domain face index
            lo_pl, hi_pl = {}, {}
            for a, lo, v0, v1 in boxes:
                for plane, face in ((lo_pl, 0), (hi_pl, N)):
                    if not (v0[ax] <= face <= v1[ax]): continue
                    if (plane is lo_pl and v0[ax] != 0) or (plane is hi_pl and v1[ax] != N): continue
                    sl = [slice(int(v0[q] - lo[q]), int(v1[q] - lo[q]) + 1) for q in (2, 1, 0)]   # a[k, j, i]
                    sl[2 - ax] = int(face - lo[ax])
                    t = [q for q in (0, 1, 2) if q != ax]
                    vals = a[tuple(sl)]       # remaining axes in (k, j, i) order without ax
                    other = [q for q in (2, 1, 0) if q != ax]
                    for idx in np.ndindex(vals.shape):
                        key = tuple(int(v0[q]) + idx[m] for m, q in enumerate(other))
                        plane[key] = vals[idx]
            keys = set(lo_pl) & set(hi_pl)
            diff = [abs(lo_pl[k] - hi_pl[k]) for k in keys]
            nbad = int(sum(1 for x in diff if x > FACE_TOL)); bad += nbad
            print('  %s: low and high domain faces differ by more than 1e-12 in %d of %d cells (max |diff| %.3e)' % (n, nbad, len(diff), max(diff)))
            assert len(keys) > 0
        return bad
    f = {'mu_edge_corner': 'MU', 'kres_edge_corner': 'KRES'}[chk]
    bad = 0
    boxes = [load(d, f, b) for b in range(nb_)]
    Nd = [max(int(v1[q]) for _, _, _, v1 in boxes) for q in range(3)]    # highest valid cell index per direction
    for a, lo, v0, v1 in boxes:
        assert tuple(lo) == tuple(v0 - 1), (lo, v0)
        sides = [bool(v0[q] == 0) for q in range(3)] + [bool(v1[q] == Nd[q]) for q in range(3)]
        XL, YL, ZL, XH, YH, ZH = sides
        bad += int((edges(a, XL, XH, YL, YH, ZL, ZH) != a).sum())
    print('  %s: %d edge/corner cells differ from the FDS clamped copies (%d box(es))' % (f, bad, nb_))
    return bad

if __name__ == '__main__':
    bad = main(sys.argv[1], sys.argv[2])
    print('  %s: %s' % (sys.argv[2], 'HOLDS' if bad == 0 else 'VIOLATED (%d)' % bad))
    sys.exit(1 if bad else 0)
