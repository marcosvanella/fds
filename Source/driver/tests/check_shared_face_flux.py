#!/usr/bin/env python3
"""Do two level-0 boxes that share a face hold bitwise the same face flux? (question of the flux-override note, section 6.2)
usage: check_shared_face_flux.py <run-dir 1 box> <run-dir 4 boxes> <cycle>
Both runs are of tests/cases/dec1.fds and dec2.fds (32x1x32 periodic in x and z; dec2 = 2x2 boxes of 16x16), 1 rank, produced with
   FDSTL_STAGES=<cycle> FDSTL_STAGEG=<cycle> fds_amr <case>.fds --run --outdir . --chid <case> --quiet --steps <cycle>
Compared (predictor, pass 1): FX and FZ (advective face value, MASS_FINITE_DIFFERENCES), SWORK1 and SWORK3 after DIVERGENCE_PART_1 (diffusive species flux RHO_D_DZDX, RHO_D_DZDZ),
and the face velocities U and W before the flux (raw stage dump), at the box-to-box x faces (b0|b1, b2|b3), the periodic-wrap x faces (b1|b0, b3|b2), and the box-to-box z faces.
Prints, per array: faces compared, differences between the two boxes that hold the face, differences from the 1-box run. Exit status 1 if any two boxes disagree."""
import sys, numpy as np
d1, d2, cyc = sys.argv[1], sys.argv[2], sys.argv[3]
def ld(p, rank_hdr=True):
    h = np.fromfile(p, dtype=np.int32, count=16)
    if rank_hdr: r = h[0]; lb = h[1:1 + r]; ub = h[5:5 + r]
    else: lb = h[0:3]; ub = h[3:6]
    return np.fromfile(p, dtype=np.float64, offset=64).reshape((ub - lb + 1)[::-1]), lb
bad_boxes = 0
# box layout of dec2: b0 (x 0-15, z 0-15), b1 (x 16-31, z 0-15), b2 (x 0-15, z 16-31), b3 (x 16-31, z 16-31); FDS box-local index I = 0..16 (face), J = 1, K or I = 1..16
def check(name, tag, raw, direction, comps):
    """direction 'x': faces normal to x (pairs b0|b1, b2|b3 and the periodic wrap b1|b0, b3|b2); 'z': pairs b0|b2, b1|b3. raw: AMReX-indexed 3-D stage dump, else FDS-indexed scratch dump."""
    global bad_boxes
    def get(d, b):
        return ld('%s/raw_%s_%s_b%d.bin' % (d, tag, name, b), False) if raw else ld('%s/scr_%s_%s_b%d.bin' % (d, tag, name, b))
    def val(t, face, tang, n, glob):
        """face: the face index along `direction` and tang: the cell index along the other direction, in FDS box-local indices; glob: (x offset, z offset) of the box, for raw files"""
        a, lb = t
        i, k = (face, tang) if direction == 'x' else (tang, face)
        if raw: return a[(k - 1 + glob[1] if direction == 'x' else face + glob[1]) - lb[2], 0 - lb[1], (face + glob[0] if direction == 'x' else tang - 1 + glob[0]) - lb[0]]
        return a[n - lb[3], k - lb[2], 1 - lb[1], i - lb[0]]
    one = get(d1, 0); bx = [get(d2, i) for i in range(4)]
    org = [(0, 0), (16, 0), (0, 16), (16, 16)]
    if direction == 'x': pairs = [(0, 1, 16, 0), (2, 3, 16, 0), (1, 0, 16, 0), (3, 2, 16, 0)]   # (box with the high face, box with the low face, face of the first, face of the second)
    else: pairs = [(0, 2, 16, 0), (1, 3, 16, 0)]
    tot = nb = nr = 0
    for a, b, fa, fb in pairs:
        wrap = (direction == 'x' and a in (1, 3))
        for n in comps:
            for t in range(1, 17):
                va = val(bx[a], fa, t, n, org[a]); vb = val(bx[b], fb, t, n, org[b])
                tot += 1; nb += (va != vb)
                if not wrap:   # compare with the 1-box run at the global face (x: 16, z: 16; second row of boxes: +16 along the tangential direction)
                    toff = org[a][1] if direction == 'x' else org[a][0]
                    vg = val(one, 16, t + toff, n, (0, 0)); nr += (va != vg or vb != vg)
    bad_boxes += nb
    print('%-7s %s faces: %2d compared, %d differ between the two boxes holding the face, %d differ from the 1-box run' % (name, direction, tot, nb, nr))
for name, tag, raw, direction, comps in (('FX', 'p1_mfd', False, 'x', (1, 2)), ('FZ', 'p1_mfd', False, 'z', (1, 2)), ('SWORK1', 'p1_div1', False, 'x', (1, 2)), ('SWORK3', 'p1_div1', False, 'z', (1, 2)),
                                         ('U', 'p1_prevflux', True, 'x', (0,)), ('W', 'p1_prevflux', True, 'z', (0,))):
    check(name, tag, raw, direction, comps)
print('SHARED FACES: the boxes agree bitwise' if bad_boxes == 0 else 'SHARED FACES: DIFFER (%d)' % bad_boxes)
sys.exit(1 if bad_boxes else 0)
