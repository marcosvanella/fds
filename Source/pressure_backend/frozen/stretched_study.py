#!/usr/bin/env python3
"""Stretched-mesh zero-mode study against FDS ULMAT/GLMAT dumps (numpy only).

Usage: stretched_study.py [--periodic-xy] DUMP [DUMP ...]   (dbg_ULMAT_<n>.txt / dbg_GLMAT_<n>.txt from the scratch hook)

Rebuilds the FDS volume-scaled 7-point finite-volume matrix K (rows scaled by the cell volume, entries
area/centre-to-centre distance, periodic faces wrapped) from the dumped cell widths and solves K x = F1 with the
identity-row pin on the last unknown (the ULMAT/UGLMAT HYPRE setup), for two ways of removing the zero mode of the
right-hand side F (the dumped FDS F_H, boundary terms and any scratch offset included):
  arith  F1 = F - mean(F)                      what FDS does (plain count mean of the volume-scaled RHS)
  vol    F1 = F - v * sum(F)/sum(v)            volume-weighted mean of F/v removed from F/v
  none   F1 = F                                 no removal (what whole-domain GLMAT does on a single MPI rank, see note)
then the same final gauge as FDS (sum(rho*v*(KRES + x)) / sum(rho*v) removed). Compares with the dumped FDS fields."""
import sys
import numpy as np

def load(path):
    hdr, rows = {}, []
    for line in open(path):
        t = line.split()
        if not t: continue
        if t[0] in ("DIMS", "MTYPE"): hdr[t[0]] = [float(x) for x in t[1:]]
        elif t[0] in ("DX", "DY", "DZ", "DXN", "DYN", "DZN"): hdr[t[0]] = np.array([float(x) for x in t[1:]])
        elif t[0] == "COLS": cols = t[1:]
        else: rows.append([float(x) for x in t])
    a = np.array(rows)
    nx, ny, nz = (int(x) for x in hdr["DIMS"][:3])
    f = {}
    for c, name in enumerate(cols):
        if name in "IJK": continue
        v = np.full((nx, ny, nz), np.nan)
        v[a[:, 0].astype(int) - 1, a[:, 1].astype(int) - 1, a[:, 2].astype(int) - 1] = a[:, c]
        f[name] = v
    return hdr, f, (nx, ny, nz)

def build_K(hdr, n, periodic):
    nx, ny, nz = n
    DX, DY, DZ = hdr["DX"], hdr["DY"], hdr["DZ"]
    DN = (hdr["DXN"], hdr["DYN"], hdr["DZN"])           # DXN(i): centre i to centre i+1, i = 0..N
    W = (DX, DY, DZ)
    N = nx * ny * nz
    K = np.zeros((N, N))
    idx = lambda i, j, k: i + nx * (j + ny * k)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                a = idx(i, j, k)
                for d in range(3):
                    ijk = [i, j, k]
                    ijk[d] += 1
                    if ijk[d] == n[d]:
                        if not periodic[d]: continue
                        ijk[d] = 0
                    # face area: product of the other two widths of cell a; distance DN[d][index of low cell + 1]
                    others = [W[e][(i, j, k)[e]] for e in range(3) if e != d]
                    c = others[0] * others[1] / DN[d][(i, j, k)[d] + 1]    # arrays are 0-based on 1..N cells; DN has N+1 entries
                    b = idx(*ijk)
                    K[a, a] += c; K[b, b] += c; K[a, b] -= c; K[b, a] -= c
    return K

def rel(a, b, mask=None):
    return float(np.linalg.norm((a - b).ravel()) / np.linalg.norm(b.ravel()))

def grad_diff(H1, H2, hdr):
    """Relative L2 difference of the interior face gradients (velocity correction ~ -dt dH/dn)."""
    num = den = 0.0
    for d in range(3):
        D = np.diff(H1, axis=d) - np.diff(H2, axis=d)
        G = np.diff(H2, axis=d)
        sh = [1, 1, 1]; sh[d] = -1
        dn = hdr[("DXN", "DYN", "DZN")[d]][1:-1].reshape(sh)       # interior faces: DN(1..N-1)
        num += np.sum((D / dn) ** 2); den += np.sum((G / dn) ** 2)
    return float(np.sqrt(num / den))

def study(path, periodic):
    hdr, f, n = load(path)
    vol = hdr["DX"][:, None, None] * hdr["DY"][None, :, None] * hdr["DZ"][None, None, :]
    F0 = f["FH0"]; FH1 = f["FH1"]; rho = f["RHOP"]; kres = f["KRES"]
    order = lambda a: a.ravel(order="F")
    print(f"== {path}: n={n} periodic={periodic} MTYPE={int(hdr['MTYPE'][0])} rhs_offset={hdr['MTYPE'][1]:g}")
    print(f"   volume ratio max/min = {vol.max() / vol.min():.3f}; rho range {rho.min():.4f}..{rho.max():.4f}")
    inner = np.zeros(n, bool); inner[1:-1, 1:-1, 1:-1] = True
    if hdr['MTYPE'][1] == 0:
        print(f"   F_H0 / (PRHS*volume) - 1 on interior cells: max {np.abs(f['FH0'][inner] / (f['PRHS'][inner] * vol[inner]) - 1).max():.2e}")
    F = order(F0); v = order(vol)
    meanF = F.mean(); c_vol = F.sum() / v.sum()
    rmsF = np.sqrt(np.mean(F ** 2)); rms_b = np.sqrt(np.mean((F / v) ** 2))
    print(f"   FDS removed mean of F: mean(F) = {order(F0).mean() - order(FH1).mean():.6e} (dumped F_H0 - F_H1), numpy mean(F) = {meanF:.6e}, mean(F)/rms(F) = {abs(meanF) / rmsF:.3e}")
    print(f"   volume-weighted alternative: c = sum(F)/sum(v) = {c_vol:.6e}, |c|/rms(F/v) = {abs(c_vol) / rms_b:.3e}")
    K = build_K(hdr, n, periodic)
    out = {}
    for name, F1 in (("arith", F - meanF), ("vol", F - v * c_vol), ("none", F.copy())):
        x = np.zeros_like(F1)
        x[:-1] = np.linalg.solve(K[:-1, :-1], F1[:-1])
        res_full = np.linalg.norm(K @ x - F1) / np.linalg.norm(F1)
        sumF1 = F1.sum() / np.abs(F1).sum()
        if name == "none": F1 = F1.copy()
        # solution zero modes: arithmetic and volume-weighted (irrelevant after the final gauge, shown pre-gauge)
        xa = x - x.mean()
        xv = x - np.sum(v * x) / v.sum()
        r = order(rho); k = order(kres)
        gauge = lambda y: y - np.sum(v * r * (k + y)) / np.sum(v * r)
        H = (-gauge(xa)).reshape(n, order="F")
        Hm = H - H.mean()
        Hf = f["HP"]; Hfm = Hf - Hf.mean()
        out[name] = H
        print(f"   [{name}] true residual of the scaled system ||K x - F1||/||F1|| = {res_full:.2e}; sum(F1)/sum|F1| = {sumF1:.2e}")
        print(f"   [{name}] x (pinned solve) mean/rms = {x.mean() / np.sqrt(np.mean(x ** 2)):.3e}")
        print(f"   [{name}] vs FDS dumped F_H1 (after its mean removal): rel L2 = {rel(F1.reshape(n, order='F'), FH1):.3e}")
        print(f"   [{name}] vs FDS pre-gauge X after its arithmetic mean removal (XH1): arith-mean removal rel L2 = "
              f"{rel(xa.reshape(n, order='F'), f['XH1']):.3e}, volume-mean removal rel L2 = {rel(xv.reshape(n, order='F'), f['XH1']):.3e}")
        print(f"   [{name}] final H vs FDS: mean-removed rel L2 = {rel(Hm, Hfm):.3e}, raw rel L2 = {rel(H, Hf):.3e}, "
              f"max abs raw = {np.abs(H - Hf).max():.3e}, gradient (velocity correction) rel L2 = {grad_diff(H, Hf, hdr):.3e}")
    print(f"   [arith vs vol] difference of the two solutions (mean-removed rel L2) = {rel(out['arith'] - out['arith'].mean(), out['vol'] - out['vol'].mean()):.3e}")

if __name__ == "__main__":
    args = sys.argv[1:]
    per = [False, False, False]
    if args and args[0] == "--periodic-xy":
        per = [True, True, False]; args = args[1:]
    for p in args:
        study(p, per)
