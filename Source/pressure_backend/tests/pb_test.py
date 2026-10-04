#!/usr/bin/env python3
"""CTest driver for the pressure backend harness. Parses the harness RESULT/CMP/CHECK lines.
Subcommands: selector, exactsum, fftmlmg, frozen, decomp, repeat, singular, ulmat, ulmatgauge, meankind,
residual check: hypresid, rescheck, reslimit, hypresidcomp;
composite: compconv, compfull, compdecomp, comp3, compgrad, compshape, compmixed, compns2d, compsel, compws (composite two-/three-level
solves; manufactured solutions, see harness/composite_modes.cpp and frozen/composite-notes.md).
Exit 0 = pass. ulmat/ulmatgauge/meankind compare against an independent numpy computation that follows FDS
ULMAT (volume-scaled rows, arithmetic mean removal of F and of X, identity-pin reduced system, rho*volume gauge)."""
import argparse, math, os, re, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("cmd")
ap.add_argument("--harness", required=True)
ap.add_argument("--mpiexec", required=True)
ap.add_argument("--work", required=True)
ap.add_argument("--bc", default="neumann")
ap.add_argument("--n", default="64 64 64")
ap.add_argument("--ratio", type=int, default=2)
ap.add_argument("--ns", default="32 64")       # coarse sizes for the composite convergence runs
ap.add_argument("--plane", type=int, default=0)
ap.add_argument("--be", default="fft mlmg")  # backends for decomp/repeat/mixedfaces/bcdata (space separated)
ap.add_argument("--faces", default="")        # per-face types, e.g. ND,NN,DD (letters N D P; low then high face per direction)
A = ap.parse_args()
os.makedirs(A.work, exist_ok=True)
fails = []

def kv(line): return dict(re.findall(r"(\w+)=(\S+)", line))
def eps_H(n): m = max(int(x) for x in n.split()); return max(1e-8, 2.4e-12 * m * m)

def run(np_, **kw):
    # Each key=value is one argv entry (values may contain spaces, e.g. n_cell="64 64 64").
    cmd = [A.mpiexec, "--oversubscribe", "--bind-to", "none", "-np", str(np_), A.harness]
    cmd += [f"{k}={v}" for k, v in kw.items()]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=500)
    return p.returncode, p.stdout + p.stderr

def lines(out, tag): return [l for l in out.splitlines() if l.startswith(tag + " ")]
def ok(cond, what):
    print(("PASS " if cond else "FAIL ") + what)
    if not cond: fails.append(what)

def base(**extra):
    d = dict(mode="solve", n_cell=A.n, bc=A.bc, verbose=1)
    d.update(extra); return d

def fftmlmg():
    rc, out = run(4, **base(max_grid_size=16 if A.n != "64 64 64" else 32, backends="fft mlmg"))
    ok(rc == 0, "harness exit 0"); print(out)
    cmp_ = lines(out, "CMP"); ok(len(cmp_) == 1, "one CMP line")
    d = kv(cmp_[0]); e = eps_H(A.n)
    print(f"RESULT rel_l2={d['rel_l2']} eps_H={e:.3e}")
    ok(float(d["rel_l2"]) <= e, "FFT vs MLMG rel L2 <= eps_H")
    for l in lines(out, "RESULT"):
        r = kv(l)
        ok(r["status"] == "Ok" and r["residual_ok"] == "1" and r["nwarn"] == "0", f"{l.split()[1]} status Ok, true residual within tolerance, no warning")
        ok(float(r["true_rel2"]) <= 1e-12, f"{l.split()[1]} true rel2 {r['true_rel2']} <= 1e-12")

def frozen():
    pre = os.path.join(A.work, "case")
    rc, out = run(2, mode="gen", n_cell=A.n, bc=A.bc, max_grid_size=16, out=pre); print(out)
    ok(rc == 0, "generator exit 0")
    rc, out = run(4, **base(max_grid_size=8, backends="mlmg fft", rhs_file=pre + "_rhs.bin", ref_file=pre + "_H.bin"))
    print(out); ok(rc == 0, "harness exit 0")
    e = eps_H(A.n)
    for l in lines(out, "CMP"):
        d = kv(l); ok(float(d["rel_l2"]) <= e, f"{l.split()[1]} rel L2 {d['rel_l2']} <= {e:.2e}")
    ok(len(lines(out, "CMP")) == 2, "two CMP lines")

def solve_files(np_, mgs, tag):
    pre = os.path.join(A.work, tag)
    rc, out = run(np_, **base(max_grid_size=mgs, backends=A.be, out=pre))
    ok(rc == 0, f"run {tag} exit 0")
    hashes = {l.split()[1]: kv(l)["hash"] for l in lines(out, "RESULT")}
    pins = {l.split()[1]: (kv(l).get("comp0_pin"), kv(l).get("comp0_removed_mean")) for l in lines(out, "RESULT")}
    return pre, hashes, pins

def diffrel(a, b):
    rc, out = run(1, mode="diff", n_cell=A.n, bc=A.bc, a=a, b=b)
    d = kv(lines(out, "CMP")[0]); return float(d["rel_l2"])

def decomp():
    e = eps_H(A.n)
    ref_pre, ref_h, ref_pins = solve_files(1, 32, "np1_mgs32")
    for np_ in (1, 2, 4):
        for mgs in (32, 16):
            if (np_, mgs) == (1, 32): continue
            pre, h, pins = solve_files(np_, mgs, f"np{np_}_mgs{mgs}")
            for b in A.be.split():
                r = diffrel(f"{pre}_{b}.bin", f"{ref_pre}_{b}.bin")
                bit = h[b] == ref_h[b]
                print(f"DECOMP {A.bc} np={np_} mgs={mgs} {b}: rel_l2 vs np1/mgs32 = {r:.3e} bitwise={'yes' if bit else 'no'}")
                ok(r <= e, f"{b} np={np_} mgs={mgs} within eps_H")
            ok(pins == ref_pins, f"pin and removed mean bitwise identical np={np_} mgs={mgs}")

def repeat():
    hs = []
    for i in range(3):
        _, h, _ = solve_files(4, 16, f"rep{i}"); hs.append(h)
    for b in A.be.split():
        same = len({h[b] for h in hs}) == 1
        print(f"REPEAT {A.bc} {b}: hashes {[h[b] for h in hs]}")
        ok(same, f"{b} 3 repeats bitwise identical")

def selector():
    rc, out = run(1, mode="selector"); print(out); ok(rc == 0, "selector checks")
    rc, out = run(2, mode="selector"); print(out); ok(rc == 0, "selector checks, 2 ranks")

def exactsum():
    res = {}
    for np_, mgs in ((1, 32), (2, 16), (4, 8), (3, 12)):
        rc, out = run(np_, mode="exactsum", n_cell="40 24 36", max_grid_size=mgs, bc="neumann")
        ok(rc == 0, f"exactsum np={np_} mgs={mgs} checks")
        es = kv(lines(out, "EXACTSUM")[0]); mr = lines(out, "MEANREMOVAL")
        res[(np_, mgs)] = (es["sum0"], es["sum1"], es["n0"], es["n1"], kv(mr[0])["hash1"], kv(mr[0])["removed0"], kv(mr[0])["removed1"])
    print(res)
    ok(len(set(res.values())) == 1, "exact sums, removed means and mean-removed fields bitwise identical across np/max_grid_size")

def singular():
    n = "16 16 16"
    rc, out = run(1, mode="solve", n_cell=n, bc="neumann", max_grid_size=16, backends="mlmg", remove_mean=1, rhs_offset=1.0, verbose=1)
    print(out); r = kv(lines(out, "RESULT")[0])
    ok(r["residual_ok"] == "1" and float(r["true_rel2"]) <= 1e-12, "mean removal on: true residual passes with offset RHS")
    ok(float(r["comp0_removed_rel"]) > 1e-10 and r["nwarn"] == "1", "removed-mean diagnostic warns on injected offset")
    rc, out = run(1, mode="solve", n_cell=n, bc="neumann", max_grid_size=16, backends="mlmg", remove_mean=0, rhs_offset=1.0, verbose=1)
    print(out); r = kv(lines(out, "RESULT")[0])
    ok(r["residual_ok"] == "0" and float(r["true_rel2"]) > 1e-6, f"hook off: true residual {r['true_rel2']} exceeds tolerance and warns")
    rc, out = run(1, mode="solve", n_cell=n, bc="neumann", max_grid_size=16, backends="mlmg", remove_mean=1, verbose=1)
    r = kv(lines(out, "RESULT")[0])
    ok(r["nwarn"] == "0" and float(r["comp0_removed_rel"]) < 1e-10, "quiet without offset (removed mean at round-off)")

# ---- HYPRE assembled-matrix backend (frozen/hypre-notes.md) ----------------------------------------------
def hypsingle():
    """Single level: FFT, MLMG and HYPRE on the same synthetic problem; HYPRE agrees with both within eps_H."""
    e = eps_H(A.n)
    pre = os.path.join(A.work, "s")
    rc, out = run(4, **base(max_grid_size=16, backends="fft mlmg hypre", out=pre)); print(out)
    ok(rc == 0, "harness exit 0")
    for l in lines(out, "RESULT"):
        r = kv(l)
        ok(r["status"] == "Ok" and r["residual_ok"] == "1" and r["nwarn"] == "0", f"{l.split()[1]} status Ok, true residual within tolerance, no warning")
    h = [kv(l) for l in lines(out, "RESULT") if l.split()[1] == "hypre"][0]
    ok(h["backend"] == "HYPRE", "the HYPRE result was produced by the HYPRE backend")
    ok(float(h["true_rel2"]) <= 1e-11, f"HYPRE true rel2 {h['true_rel2']} <= 1e-11 (independent residual)")
    r_hf = diffrel(f"{pre}_hypre.bin", f"{pre}_fft.bin"); r_hm = diffrel(f"{pre}_hypre.bin", f"{pre}_mlmg.bin"); r_mf = diffrel(f"{pre}_mlmg.bin", f"{pre}_fft.bin")
    print(f"HYPSINGLE {A.bc} n=[{A.n}] hypre_iters={h['iters']} rel_l2 hypre-fft={r_hf:.3e} hypre-mlmg={r_hm:.3e} mlmg-fft={r_mf:.3e} eps_H={e:.1e}")
    ok(r_hf <= e, f"HYPRE vs FFT rel L2 {r_hf:.2e} <= eps_H"); ok(r_hm <= e, f"HYPRE vs MLMG rel L2 {r_hm:.2e} <= eps_H")

def hypci():
    """CI check: the same frozen (generated, reference H from the FFT) solves through MLMG and HYPRE must agree within eps_H."""
    for bc in ("neumann", "periodic", "dirichlet"):
        pre = os.path.join(A.work, f"case_{bc}")
        rc, out = run(2, mode="gen", n_cell=A.n, bc=bc, max_grid_size=16, out=pre); ok(rc == 0, f"{bc}: generator exit 0")
        rc, out = run(4, mode="solve", n_cell=A.n, bc=bc, verbose=1, max_grid_size=8, backends="mlmg hypre", rhs_file=pre + "_rhs.bin", ref_file=pre + "_H.bin", out=pre + "_o")
        print(out); ok(rc == 0, f"{bc}: harness exit 0")
        e = eps_H(A.n)
        for l in lines(out, "CMP"):
            d = kv(l); ok(float(d["rel_l2"]) <= e, f"{bc}: {l.split()[1]} vs the frozen reference rel L2 {d['rel_l2']} <= {e:.2e}")
        r = subprocess.run([A.mpiexec, "--oversubscribe", "-np", "1", A.harness, "mode=diff", f"a={pre}_o_hypre.bin", f"b={pre}_o_mlmg.bin", f"n_cell={A.n}", f"bc={bc}"], capture_output=True, text=True)
        d = kv(lines(r.stdout + r.stderr, "CMP")[0]); rl = float(d["rel_l2"])
        print(f"HYPCI {bc} n=[{A.n}] hypre-vs-mlmg rel_l2={rl:.3e} max_abs={d.get('max_abs')} eps_H={e:.1e}")
        ok(rl <= e, f"{bc}: HYPRE vs MLMG rel L2 {rl:.2e} <= eps_H {e:.1e} (CI gate)")
    # negative control: a loosely converged pair must be flagged by the same comparison (the gate can fail)
    pre = os.path.join(A.work, "case_neumann")
    rc, out = run(4, mode="solve", n_cell=A.n, bc="neumann", verbose=0, max_grid_size=8, backends="mlmg hypre", tol_rel=1e-3, out=pre + "_loose")
    r = subprocess.run([A.mpiexec, "--oversubscribe", "-np", "1", A.harness, "mode=diff", f"a={pre}_loose_hypre.bin", f"b={pre}_loose_mlmg.bin", f"n_cell={A.n}", "bc=neumann"], capture_output=True, text=True)
    d = kv(lines(r.stdout + r.stderr, "CMP")[0]); rl = float(d["rel_l2"])
    print(f"HYPCI control (tol_rel=1e-3) hypre-vs-mlmg rel_l2={rl:.3e}")
    ok(rl > eps_H(A.n) and "FAIL" in (r.stdout + r.stderr), f"negative control: loosely converged HYPRE and MLMG differ by {rl:.1e} > eps_H and the comparison reports FAIL")

def hypcomp():
    """Composite: HYPRE vs MLMG, independent true residual, bitwise repeatability, then 1/2/4-rank independence."""
    cases = [("neumann", 2, 2, "layout=0"), ("neumann", 2, 4, "layout=0"), ("periodic", 2, 2, "layout=1"), ("periodic", 2, 4, "layout=0"),
             ("dirichlet", 2, 2, "layout=0"), ("dirichlet", 2, 4, "layout=2"), ("neumann", 3, 2, "layout=0"), ("periodic", 3, 4, "layout=0"),
             ("neumann", 2, 2, "full=1"), ("neumann", 2, 2, "plane2d=1"), ("periodic", 2, 2, "plane2d=1"),
             ("neumann", 2, 2, "bcfaces=ND,DN,NN"), ("neumann", 2, 2, "bcfaces=PP,NN,DD")]
    for bc, nlev, ratio, extra in cases:
        kw = dict(mode="hypre_cmp", n=16, nlev=nlev, ratio=ratio, mgs=8)
        if extra.startswith("bcfaces"): kw["bcfaces"] = extra.split("=")[1]
        else: kw["bc"] = bc; kw[extra.split("=")[0]] = extra.split("=")[1]
        rc, out = run(3 if nlev == 3 else 2, **kw)
        ok(rc == 0, f"hypre_cmp {bc} nlev={nlev} ratio={ratio} {extra}: harness checks")
        if rc != 0: print(out)
        l = lines(out, "HYPCMP")[0]; print(l); d = kv(l)
        ok(float(d["rel2_diff"]) <= 1e-8 and float(d["hypre_true_rel2"]) <= 1e-8 and d["repeat_bitwise"] == "1",
           f"{bc} nlev={nlev} ratio={ratio} {extra}: HYPRE vs MLMG rel L2 {d['rel2_diff']}, true residual {d['hypre_true_rel2']}, iters hypre {d['hypre_iters']} / mlmg {d['mlmg_iters']} ({d['method']})")
    # rank independence: the composite solution on 1, 2, 4 ranks (different box distribution) within eps_H of the 1-rank result
    for bc, ratio in (("neumann", 2), ("periodic", 4)):
        ref = None
        for np_ in (1, 2, 4):
            pre = os.path.join(A.work, f"cd_{bc}_{np_}")
            rc, out = run(np_, mode="hypre_cmp", n=16, nlev=2, ratio=ratio, bc=bc, mgs=8, out=pre)
            ok(rc == 0, f"hypre_cmp {bc} r{ratio} np={np_}: harness checks")
            x = np.fromfile(pre + "_phi.bin")
            if ref is None: ref = x; continue
            r = float(np.linalg.norm(x - ref) / np.linalg.norm(ref))
            print(f"HYPDECOMP composite {bc} ratio={ratio} np={np_} vs np=1: rel_l2 {r:.3e}")
            ok(r <= 1e-8, f"composite HYPRE {bc} ratio={ratio}: np={np_} equals np=1 within eps_H (rel L2 {r:.2e})")

def hypcache():
    for np_ in (1, 2):
        rc, out = run(np_, mode="hypcache", n_cell=A.n, mgs=16, nsolve=5, bcpairs=A.bc); print(out)
        ok(rc == 0, f"HYPRE set-up cache checks (reuse, bitwise identical, rebuild invalidation), np={np_}")


# ---- independent numpy references (no AMReX) -------------------------------------------------------------
def read_field(path, n):
    nx, ny, nz = (int(x) for x in n.split())
    return np.fromfile(path, dtype="<f8").reshape(nz, ny, nx).transpose(2, 1, 0)   # raw is x fastest -> [i,j,k]

def build_K(n, h, bc):
    """FDS ULMAT matrix for uniform cells: sum over faces of AF/DXN * [[1,-1],[-1,1]] (= -volume * 7-point L)."""
    nx, ny, nz = (int(x) for x in n.split())
    N = nx * ny * nz
    K = np.zeros((N, N))
    idx = lambda i, j, k: i + nx * (j + ny * k)
    c = h * h / h                                   # AF / DXN
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                a = idx(i, j, k)
                for d, (ii, jj, kk) in enumerate(((i + 1, j, k), (i, j + 1, k), (i, j, k + 1))):
                    if (ii, jj, kk)[d] == (nx, ny, nz)[d]:
                        if bc != "periodic": continue
                        ii, jj, kk = (ii % nx, jj % ny, kk % nz)
                    b = idx(ii, jj, kk)
                    K[a, a] += c; K[b, b] += c; K[a, b] -= c; K[b, a] -= c
    return K

def ulmat_like(rhs, n, h, bc, rho=None, kres=None):
    """FDS ULMAT HYPRE path on uniform cells. Returns (H, mean(F)/rms(F), F-mean before removal)."""
    vol = h ** 3
    F = (rhs * vol).ravel(order="F")                # x fastest
    meanF = F.mean()
    F1 = F - meanF                                  # arithmetic mean of the scaled RHS (pres.f90 1683-1733)
    K = build_K(n, h, bc)
    x = np.zeros_like(F1)
    x[:-1] = np.linalg.solve(K[:-1, :-1], F1[:-1])  # identity-row pin on the last unknown (1748, 2965-2982)
    x -= x.mean()                                   # arithmetic mean of X_H (1779-1828)
    H = -x                                          # HP = -X_H (1892)
    if rho is not None:                             # rho*volume gauge on KRES + X_H (1830-1878)
        r = rho.ravel(order="F"); kr = kres.ravel(order="F")
        shift = np.sum(vol * r * (kr + x)) / np.sum(vol * r)
        H = -(x - shift)
    return H.reshape(rhs.shape, order="F"), abs(meanF) / np.sqrt(np.mean(F ** 2)), meanF

def relL2(a, b): return float(np.linalg.norm((a - b).ravel()) / np.linalg.norm(b.ravel()))

def ulmat():
    n = A.n; h = 1.0 / max(int(x) for x in n.split())
    pre = os.path.join(A.work, "gen")
    rc, out = run(1, mode="gen", n_cell=n, bc=A.bc, max_grid_size=16, out=pre); ok(rc == 0, "generator exit 0")
    for offset in (0.0, 0.75):                      # compatible RHS, and RHS with a large nonzero mean
        tag = f"off{offset}"
        rc, out = run(2, **base(max_grid_size=4, backends="fft mlmg", rhs_file=pre + "_rhs.bin", rhs_offset=offset,
                                mean_kind="scaled", out=os.path.join(A.work, tag)))
        ok(rc == 0, f"harness exit 0 (rhs offset {offset})"); print(out)
        rhs = read_field(pre + "_rhs.bin", n) + offset
        H, rel, meanF = ulmat_like(rhs, n, h, A.bc)
        print(f"ULMAT-like: mean(F)={meanF:.3e} rel to rms(F)={rel:.3e}")
        for l in lines(out, "RESULT"):
            name = l.split()[1]; r = kv(l)
            phi = read_field(os.path.join(A.work, f"{tag}_{name}.bin"), n)
            d = relL2(phi, H)
            print(f"EQUIV {A.bc} offset={offset} {name}: rel_l2 vs ULMAT-like = {d:.3e}  (removed_rel ours {float(r['comp0_removed_rel']):.6e}, numpy {rel:.6e})")
            ok(d <= 1e-9, f"{name} equals ULMAT-like arithmetic-mean path (rel L2 {d:.2e} <= 1e-9), offset {offset}")
            ok(abs(float(r["comp0_removed_rel"]) - rel) <= 1e-9 * max(rel, 1e-3) + 1e-12, f"{name} removed mean (relative) matches numpy arithmetic mean of F")
            ok(abs(phi.mean()) <= 1e-12 * np.abs(phi).max(), f"{name} solution has zero arithmetic (= volume-weighted) mean after the gauge")

def ulmatgauge():
    n = A.n; h = 1.0 / max(int(x) for x in n.split())
    pre = os.path.join(A.work, "gen")
    rc, out = run(1, mode="gen", n_cell=n, bc=A.bc, max_grid_size=16, out=pre); ok(rc == 0, "generator exit 0")
    outp = os.path.join(A.work, "g")
    rc, out = run(2, **base(max_grid_size=4, backends="fft mlmg", rhs_file=pre + "_rhs.bin", rhs_offset=0.4, gauge_rho=1, mean_kind="scaled", out=outp))
    ok(rc == 0, "harness exit 0"); print(out)
    rhs = read_field(pre + "_rhs.bin", n) + 0.4
    rho = read_field(outp + "_rho.bin", n); kres = read_field(outp + "_kres.bin", n)
    Hplain, _, _ = ulmat_like(rhs, n, h, A.bc)
    H, _, _ = ulmat_like(rhs, n, h, A.bc, rho, kres)
    print(f"gauge constant (rho*vol, KRES) minus plain mean gauge: {float((H - Hplain).mean()):.6e}")
    ok(abs(float((H - Hplain).mean())) > 1e-4, "rho-weighted KRES gauge differs from the plain gauge (test is sensitive)")
    for l in lines(out, "RESULT"):
        name = l.split()[1]
        phi = read_field(f"{outp}_{name}.bin", n)
        d = relL2(phi, H)
        print(f"GAUGE {A.bc} {name}: rel_l2 vs ULMAT-like rho*vol gauge = {d:.3e}")
        ok(d <= 1e-9, f"{name} matches the rho*volume, KRES gauge of ULMAT (rel L2 {d:.2e})")
        w = rho
        ok(abs(np.sum(w * (phi - kres))) <= 1e-12 * np.sum(w * np.abs(phi - kres)), f"{name} sum(rho*(H - KRES)) = 0")

def meankind():
    n = A.n
    outs = {}
    for np_, mgs in ((1, 32), (3, 6)):
        pre = os.path.join(A.work, f"mk_np{np_}")
        rc, out = run(np_, mode="meankind", n_cell=n, bc=A.bc, max_grid_size=mgs, out=pre, rhs_offset=0.5)
        ok(rc == 0, f"meankind np={np_} harness checks"); print(out)
        outs[np_] = (pre, [l for l in lines(out, "MEANKIND")])
    pre = outs[1][0]
    b0 = read_field(pre + "_b0.bin", n); v = read_field(pre + "_v.bin", n)
    rho = read_field(pre + "_rho.bin", n); kres = read_field(pre + "_kres.bin", n); phi0 = read_field(pre + "_phi0.bin", n)
    F = v * b0
    S = b0 - F.mean() / v                                    # FDS: arithmetic mean of F = v*b removed from F
    V = b0 - F.sum() / v.sum()                               # volume-weighted mean of b removed from b
    bS = read_field(pre + "_bS.bin", n); bV = read_field(pre + "_bV.bin", n)
    sc = np.abs(F).sum()
    ok(np.abs(bS - S).max() <= 1e-14 * np.abs(S).max(), "ScaledArithmetic equals b - mean(v*b)/v")
    ok(np.abs(bV - V).max() <= 1e-14 * np.abs(V).max(), "Volume equals b - sum(v*b)/sum(v)")
    ok(abs((v * bS).sum()) <= 1e-14 * sc and abs((v * bV).sum()) <= 1e-14 * sc, "sum of the scaled RHS is zero after either removal")
    dSV = relL2(bS, bV)
    print(f"stretched volumes: |b_S - b_V| / |b_V| = {dSV:.3e}")
    ok(dSV > 1e-3, "the two zero modes differ when volumes vary (test is sensitive)")
    g = phi0 - (v * phi0).sum() / v.sum()
    ok(np.abs(read_field(pre + "_gV.bin", n) - g).max() <= 1e-14 * np.abs(phi0).max(), "gauge with a volume field = volume-weighted mean removal")
    w = v * rho
    gr = phi0 - (w * (phi0 - kres)).sum() / w.sum()
    ok(np.abs(read_field(pre + "_gRK.bin", n) - gr).max() <= 1e-14 * np.abs(phi0).max(), "gauge with rho*volume weight and KRES offset")
    fl = ["b0", "v", "bS", "bV", "gV", "gRK"]
    for f in fl:
        a = open(outs[1][0] + f"_{f}.bin", "rb").read(); b = open(outs[3][0] + f"_{f}.bin", "rb").read()
        ok(a == b, f"{f} bitwise identical for np=1/mgs32 and np=3/mgs6")

def meankind_uniform():
    """Single level, uniform cells: Volume (default) and the FDS parity switch ScaledArithmetic are the same operation
    up to rounding, for FFT and MLMG, with a nonzero removed mean; the gauge (default, rho=1/KRES=0) is applied to both."""
    n = A.n; vol = (1.0 / max(int(x) for x in n.split())) ** 3
    sols = {}
    for mk in ("volume", "scaled"):
        pre = os.path.join(A.work, mk)
        rc, out = run(2, **base(max_grid_size=4, backends="fft mlmg", rhs_offset=0.75, mean_kind=mk, out=pre))
        ok(rc == 0, f"mean_kind={mk}: harness exit 0"); print(out)
        for l in lines(out, "RESULT"):
            name = l.split()[1]; r = kv(l)
            sols[(mk, name)] = (read_field(f"{pre}_{name}.bin", n), float(np.frombuffer(bytes.fromhex(r["comp0_removed_mean"]), ">f8")[0]), float(r["comp0_removed_rel"]))
            ok(r["status"] == "Ok" and r["residual_ok"] == "1", f"{mk} {name}: status Ok, true residual within tolerance")
    for name in ("fft", "mlmg"):
        (pv, mv, rv), (ps, ms, rs) = sols[("volume", name)], sols[("scaled", name)]
        d = relL2(ps, pv)
        print(f"MEANKIND uniform {A.bc} {name}: rel L2 scaled vs volume {d:.2e}; removed mean volume {mv:.15g}, scaled (mean of F) {ms:.15g}, ratio {ms / mv / vol:.15f}")
        ok(d <= 1e-12, f"{name}: ScaledArithmetic and Volume give the same solution on uniform cells (rel L2 {d:.1e})")
        ok(abs(ms - mv * vol) <= 1e-13 * abs(mv * vol), f"{name}: removed mean of ScaledArithmetic = Volume value times the cell volume")
        ok(abs(rv - rs) <= 1e-12 * abs(rv), f"{name}: relative removed mean identical for both kinds")
        ok(abs(pv.mean()) <= 1e-12 * np.abs(pv).max(), f"{name}: default gauge leaves a zero volume-weighted mean")


# ---------------------------------------------------------------------------------------------------------------
# Composite (multi-level) tests. The harness prints COMP / UNIFORM / DECOMP / GRAD / REPEAT / NS2D lines and CHECK
# PASS|FAIL lines; its exit code is nonzero when a CHECK fails.
EPS_COMP = 1e-8          # composite eps_H: true residual (relative)

def comp_run(np_, **kw):
    if A.faces: kw["bcfaces"] = A.faces
    rc, out = run(np_, mode="comp", bc=A.bc, **kw)
    ok(rc == 0, f"composite harness exit 0 ({' '.join(f'{k}={v}' for k, v in kw.items())})")
    if rc != 0: print(out)
    return out

def comp_lines(out):
    c = kv(lines(out, "COMP")[0]); u = kv(lines(out, "UNIFORM")[0]) if lines(out, "UNIFORM") else {}
    return c, u

def compconv():
    # Two-level refined patch (middle half of the box), error against the manufactured solution at two resolutions.
    ns = [int(x) for x in A.ns.split()]
    res = []
    for n in ns:
        out = comp_run(2, n=n, nlev=2, ratio=A.ratio, plane2d=A.plane, mgs=16 if n <= 32 else 32, dmkind=0)
        c, u = comp_lines(out); res.append((n, c, u)); print(lines(out, "COMP")[0]); print(lines(out, "UNIFORM")[0])
        ok(c["status"] == "Ok" and c["residual_ok"] == "1", f"n={n} status Ok, true composite residual within eps_H")
        ok(float(c["true_rel2"]) <= EPS_COMP, f"n={n} composite true residual {c['true_rel2']} <= eps_H {EPS_COMP:g}")
        ok(int(c["iters"]) <= 40, f"n={n} MLMG iterations {c['iters']} <= 40")
    (n0, c0, u0), (n1, c1, u1) = res[0], res[1]
    for key, what in (("err_l2", "composite error (all uncovered cells)"), ("comp_fine_err_l2", "composite error on the fine level")):
        src0 = c0 if key in c0 else u0; src1 = c1 if key in c1 else u1
        e0, e1 = float(src0[key]), float(src1[key]); order = math.log2(e0 / e1)
        print(f"ORDER {A.bc} ratio={A.ratio} plane={A.plane} {what}: {e0:.4e} -> {e1:.4e} order {order:.3f}")
        ok(order >= 1.8, f"{what}: order {order:.2f} >= 1.8")
    d0, d1 = float(u0["diff_fine_l2_abs"]), float(u1["diff_fine_l2_abs"]); od = math.log2(d0 / d1)
    print(f"ORDER {A.bc} composite vs uniform-fine solve on the fine region (abs L2): {d0:.4e} -> {d1:.4e} order {od:.3f}")
    ok(od >= 1.8, f"difference to the uniform fine solve converges: order {od:.2f} >= 1.8")

def compfull():
    # Fine level over the whole domain equals the single-level fine solve; coarse = average-down.
    for n, mgs in ((16, 8), (32, 16)):
        out = comp_run(2, n=n, nlev=2, ratio=A.ratio, full=1, plane2d=A.plane, mgs=mgs)
        c, u = comp_lines(out); print(lines(out, "COMP")[0]); print(lines(out, "UNIFORM")[0])
        ok(c["status"] == "Ok" and float(c["true_rel2"]) <= EPS_COMP, f"n={n} Ok, true residual {c['true_rel2']}")
        ok(float(u["diff_fine_l2_rel"]) <= 1e-10, f"n={n} full fine level vs single-level fine solve: rel L2 {u['diff_fine_l2_rel']} <= 1e-10")
        ok(float(u["diff_fine_linf"]) <= 1e-10, f"n={n} max abs difference {u['diff_fine_linf']} <= 1e-10")

def compdecomp():
    # Decomposition independence in process (box split + mapping), then across rank counts via raw dumps, repeatability.
    base = dict(n=32, nlev=2, ratio=A.ratio, plane2d=A.plane, uniform=0)
    hashes = {}
    dumps = {}
    for np_, mgs, dmk in ((1, 32, 0), (2, 16, 1), (2, 8, 2), (3, 16, 0)):
        pre = os.path.join(A.work, f"np{np_}_mgs{mgs}_dm{dmk}")
        out = comp_run(np_, mgs=mgs, dmkind=dmk, out=pre, **base)
        c, _ = comp_lines(out); print(lines(out, "COMP")[0])
        dumps[(np_, mgs, dmk)] = pre + "_phi.bin"; hashes[(np_, mgs, dmk)] = c["hash"]
        ok(float(c["true_rel2"]) <= EPS_COMP, f"np={np_} mgs={mgs} dmkind={dmk} true residual {c['true_rel2']}")
    ref = np.fromfile(dumps[(1, 32, 0)])
    for k, f in dumps.items():
        a = np.fromfile(f); rel = np.linalg.norm(a - ref) / np.linalg.norm(ref)
        print(f"DECOMP {k} vs np=1 mgs=32: rel L2 {rel:.3e}")
        ok(rel <= EPS_COMP, f"{k}: rel L2 {rel:.2e} <= eps_H {EPS_COMP:g}")
    out = comp_run(2, mgs=16, mgs2=8, dmkind=0, dmkind2=2, repeat=1, **base)
    print(lines(out, "DECOMP")[0]); print(lines(out, "REPEAT")[0])
    # run-to-run bitwise repeatability across processes
    h1 = kv(lines(comp_run(2, mgs=16, dmkind=0, **base), "COMP")[0])["hash"]
    h2 = kv(lines(comp_run(2, mgs=16, dmkind=0, **base), "COMP")[0])["hash"]
    ok(h1 == h2, f"run-to-run bitwise repeatability across processes ({h1} == {h2})")

def comp3():
    # Three levels (2 then 2, and 2 then 4 via ratio), basic residual and error checks, grad checks.
    out = comp_run(2, n=32, nlev=3, ratio=A.ratio, mgs=16, plane2d=A.plane, grad=1, mgs2=8, dmkind2=2)
    print(out)
    c, u = comp_lines(out)
    ok(c["status"] == "Ok" and float(c["true_rel2"]) <= EPS_COMP, f"three levels: Ok, true residual {c['true_rel2']} <= eps_H")
    ok(float(u["comp_fine_err_l2"]) < float(u["uni_coarse_err_l2"]), "three levels: finest-level error below the uniform coarse error")

def compgrad():
    out = comp_run(2, n=32, nlev=2, ratio=A.ratio, mgs=16, plane2d=A.plane, grad=1, uniform=0)
    print(out)
    out = comp_run(3, n=32, nlev=3, ratio=2, mgs=8, plane2d=A.plane, grad=1, uniform=0, dmkind=2)

def compshape():
    # Patch shapes: corner patch touching (and, if periodic, wrapping over) the domain faces; two separate patches;
    # a corner patch with a third level; ratio 2 and 4; odd rank count.
    for layout, nlev, ratio, np_ in ((1, 2, 2, 3), (1, 2, 4, 2), (2, 2, 2, 3), (2, 2, 4, 2), (1, 3, 2, 3)):
        out = comp_run(np_, n=32 if ratio == 2 else 16, nlev=nlev, ratio=ratio, mgs=8, layout=layout, grad=1, plane2d=A.plane, dmkind=2)
        c, u = comp_lines(out); print(lines(out, "COMP")[0]); print(lines(out, "UNIFORM")[0])
        ok(c["status"] == "Ok" and float(c["true_rel2"]) <= EPS_COMP, f"layout={layout} nlev={nlev} ratio={ratio} np={np_}: Ok, true residual {c['true_rel2']}")
        ok(float(u["comp_fine_err_l2"]) < float(u["uni_coarse_err_l2"]), f"layout={layout} nlev={nlev} ratio={ratio}: fine-level error below uniform coarse error")
        for l in lines(out, "GRAD"):
            if "div_consistency " in l: ok(float(kv(l)["rel2"]) <= EPS_COMP, f"layout={layout} nlev={nlev} ratio={ratio}: gradient divergence consistency {kv(l)['rel2']}")

def compmixed():
    # Per-direction mix of periodic and Neumann faces (all closed): 2-D plane periodic in x,z with a Neumann one-cell y
    # (the ns2d_16 pattern with a physical y wall), and 3-D mixes; patch in the middle and at the corner.
    for plane, b3, layout, ratio, np_ in ((1, "periodic neumann periodic", 0, 2, 2), (1, "neumann periodic neumann", 1, 2, 3),
                                          (0, "periodic neumann periodic", 0, 2, 2), (0, "neumann periodic neumann", 1, 4, 2),
                                          (0, "periodic periodic neumann", 2, 2, 3)):
        kw = dict(n=32 if ratio == 2 else 16, nlev=2, ratio=ratio, mgs=8, layout=layout, grad=1, plane2d=plane, dmkind=2, bcs3=f'"{b3}"')
        cmd = [A.mpiexec, "--oversubscribe", "--bind-to", "none", "-np", str(np_), A.harness, "mode=comp"] + [f"{k}={v}" for k, v in kw.items() if k != "bcs3"] + ["bcs3=" + b3]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=500); out = p.stdout + p.stderr
        ok(p.returncode == 0, f"mixed bc [{b3}] plane={plane} layout={layout} ratio={ratio} harness checks")
        if p.returncode != 0: print(out)
        c, u = comp_lines(out); print(lines(out, "COMP")[0]); print(lines(out, "UNIFORM")[0])
        ok(c["status"] == "Ok" and float(c["true_rel2"]) <= EPS_COMP, f"mixed bc [{b3}]: Ok, true residual {c['true_rel2']}")

def compns2d():
    for np_, mgs in ((1, 16), (2, 8), (4, 4)):
        rc, out = run(np_, **({"mode": "comp_ns2d", "mgs": mgs} | ({"backend": A.be} if A.be != "fft mlmg" else {})))
        print(out); ok(rc == 0, f"ns2d_16 two-level checks, np={np_} mgs={mgs}")
        l = kv(lines(out, "NS2D")[0]); ok(float(l["true_rel2"]) <= EPS_COMP, f"np={np_} true residual {l['true_rel2']}")

def read_levels(path, n, nlev, ratio, plane):
    """Level-concatenated full-domain dumps written by mode=comp_gauge (x fastest, zeros outside the BoxArray)."""
    raw = np.fromfile(path, dtype="<f8"); off = 0; out = []
    for l in range(nlev):
        r = ratio ** l
        nx, nz = n * r, n * r; ny = 1 if plane else n * r
        cnt = nx * ny * nz
        out.append(raw[off:off + cnt].reshape(nz, ny, nx).transpose(2, 1, 0)); off += cnt
    assert off == raw.size
    return out

def compgauge():
    """D-067 composite mean removal and gauge: independent numpy formulas on level dumps of a two-level hierarchy."""
    n, nlev, ratio = (32 if A.ratio == 2 else 16), 2, A.ratio
    pre = os.path.join(A.work, "g")
    rc, out = run(2, mode="comp_gauge", n=n, nlev=nlev, ratio=ratio, mgs=8 if A.ratio == 4 else 16, dmkind=1, bc=A.bc, plane2d=A.plane, out=pre)
    ok(rc == 0, "comp_gauge harness checks"); print(out)
    g = kv(lines(out, "GAUGE")[0])
    vol = [float(g[f"vol{l}"]) for l in range(nlev)]
    rhs, rho, kres, unc, pA, pB, pC = (read_levels(f"{pre}_{f}.bin", n, nlev, ratio, A.plane) for f in ("rhs", "rho", "kres", "unc", "phiA", "phiB", "phiC"))
    inside = [r != 0 for r in rho]                       # rho > 0 inside the BoxArray, 0 outside
    sel = [(u != 0) & i for u, i in zip(unc, inside)]
    def S(fields, weights=None):
        tot = np.longdouble(0)
        for l in range(nlev):
            t = fields[l][sel[l]].astype(np.longdouble) * np.longdouble(vol[l])
            if weights is not None: t = t * weights[l][sel[l]].astype(np.longdouble)
            tot += t.sum()
        return float(tot)
    nunc = sum(int(s_.sum()) for s_ in sel)
    ok(nunc == int(g["nunc"]), f"uncovered cell count {nunc} matches the solver")
    if A.bc == "dirichlet":
        ok(float(g["shiftA"]) == 0 and float(g["shiftB"]) == 0, "Dirichlet: no gauge shift (unique solution)")
        return
    Vtot = S([np.ones_like(r) for r in rho])
    meanV = S(rhs) / Vtot
    meanF = sum(float(vol[l]) * float(rhs[l][sel[l]].astype(np.longdouble).sum()) for l in range(nlev)) / nunc
    print(f"numpy: volume-weighted mean of b {meanV:.15g} (library {float(g['removedV_val']):.15g}); mean of F {meanF:.6e} (library {float(g['removedS_val']):.6e})")
    ok(abs(meanV - float(g["removedV_val"])) <= 1e-13 * abs(meanV), "removed mean (Volume) equals the independent exact volume-weighted mean of b")
    ok(abs(meanF - float(g["removedS_val"])) <= 1e-13 * abs(meanF), "removed mean (ScaledArithmetic) equals the independent arithmetic mean of v*b")
    ok(abs(meanV - meanF / vol[0]) > 1e-3 * abs(meanV), "the two zero modes differ on the hierarchy (test is sensitive)")
    # default gauge: zero volume-weighted mean
    ok(abs(S(pA)) <= 1e-13 * S([np.abs(f) for f in pA]), "default gauge: volume-weighted mean of phi is zero over the uncovered cells")
    # weighted gauge
    W = [r for r in rho]
    c = S([pA[l] - kres[l] for l in range(nlev)], W) / S(W)
    print(f"numpy gauge constant c = {c:.12e}; library shiftB - shiftA = {float(g['shiftB']) - float(g['shiftA']):.12e}")
    ok(abs(c) > 1e-3, f"weighted gauge constant {c:.3e} differs from the plain gauge (test is sensitive)")
    ok(abs((float(g["shiftB"]) - float(g["shiftA"])) - c) <= 1e-12 * max(abs(c), 1e-3), "library gauge constant equals the numpy sum(V rho (phi-KRES))/sum(V rho)")
    dev = max(np.abs((pB[l] - (pA[l] - c))[inside[l]]).max() for l in range(nlev))
    ok(dev <= 1e-13 * max(np.abs(f[i]).max() for f, i in zip(pA, inside)), f"phi with rho/KRES gauge = phi with plain gauge - c on all cells (max dev {dev:.2e})")
    res = S([kres[l] - pB[l] for l in range(nlev)], W)
    scale = S([np.abs(kres[l] - pB[l]) for l in range(nlev)], W)
    ok(abs(res) <= 1e-13 * scale, f"sum(rho V (KRES - H)) = {res:.2e} (scale {scale:.2e}) is zero over the uncovered cells")
    kscale = S([np.abs(kres[l]) for l in range(nlev)], W)
    ok(abs(res) <= 1e-14 * kscale, f"|sum(rho V (H-KRES))| = {abs(res):.2e} <= 1e-14 * sum(rho V |KRES|) = {1e-14*kscale:.2e}")

def compgaugedecomp():
    """Removed mean bitwise and gauge constants to round-off across rank counts, box sizes and distribution mappings."""
    runs = {}
    for np_, mgs, dmk, layout in ((1, 32, 0, 0), (2, 16, 1, 0), (3, 8, 2, 0), (2, 16, 0, 2), (1, 16, 0, 2), (4, 8, 2, 2), (1, 8, 1, 0), (2, 8, 1, 0), (4, 8, 1, 0)):
        rc, out = run(np_, mode="comp_gauge", n=32, nlev=2, ratio=2, mgs=mgs, dmkind=dmk, layout=layout, bc=A.bc, plane2d=A.plane)
        ok(rc == 0, f"np={np_} mgs={mgs} dmkind={dmk} layout={layout}: harness checks")
        if rc != 0: print(out)
        runs[(np_, mgs, dmk, layout)] = (kv(lines(out, "GAUGE")[0]), kv(lines(out, "GAUGE_EXACT")[0]))
    for layout in (0, 2):
        keys = [k for k in runs if k[3] == layout]
        r0 = runs[keys[0]]
        for k in keys:
            g, e = runs[k]
            print(f"DECOMP layout={layout} {k}: removedV={g['removedV']} removedS={g['removedS']} shiftA={g['shiftA']} shiftB={g['shiftB']} exact_shift={e['shift']}")
            ok(g["removedV"] == r0[0]["removedV"] and g["removedS"] == r0[0]["removedS"], f"{k}: removed means (Volume and ScaledArithmetic) bitwise identical")
            ok(e == r0[1], f"{k}: exact gauge sums of a fixed analytic field bitwise identical")
            for key in ("shiftA", "shiftB", "shiftC"):
                a, b = float(g[key]), float(r0[0][key])
                d = abs(a - b); sc = max(abs(a), abs(b), 1e-2)
                print(f"   {key}: |diff| = {d:.2e}")
                ok(d <= 1e-12 * sc, f"{k}: gauge constant {key} agrees with {keys[0]} to {d:.1e} (round-off of the solve)")


def compsel():
    for np_ in (1, 2):
        rc, out = run(np_, **({"mode": "comp_sel"} | ({"backend": A.be} if A.be != "fft mlmg" else {}))); print(out); ok(rc == 0, f"composite selector checks, np={np_}")

def compws():
    for np_ in (1, 2):
        rc, out = run(np_, **({"mode": "comp_ws"} | ({"backend": A.be} if A.be != "fft mlmg" else {}))); print(out); ok(rc == 0, f"workspace (D-058) checks, np={np_}")

def trigger1():
    for np_ in (1, 2):
        rc, out = run(np_, mode="trigger1", n_cell=A.n, bcpairs=A.bc, mgs=4 if np_ > 1 else 8)
        print(out); ok(rc == 0, f"FR-039 single-level trigger checks, bcpairs={A.bc}, np={np_}")

def comptrigger():
    for np_ in (1, 2):
        rc, out = run(np_, mode="comp_trigger"); print(out); ok(rc == 0, f"FR-039 composite trigger checks, np={np_}")

def fftcache():
    for np_ in (1, 2):
        rc, out = run(np_, mode="fftcache", n_cell=A.n, mgs=16, nsolve=6, bcpairs=A.bc); print(out)
        ok(rc == 0, f"FFT plan cache checks (bitwise identical, counters, invalidation), np={np_}")
        l = lines(out, "FFTCACHE")
        ok(len(l) == 1, "timing line present")

def dense_op(n, bc, skip=()):
    """Dense 7-point operator (unit cell width h = 1/max(n)); ghost of a Neumann face = +phi, Dirichlet = -phi, periodic wraps.
    A one-cell direction y drops its term (D-057 / FDS TWO_D); this reference is used for n >= 2 in every direction."""
    nx, ny, nz = n; N = nx*ny*nz; h = 1.0/max(n)
    ix = lambda i, j, k: i + nx*(j + ny*k)
    M = np.zeros((N, N))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                r = ix(i, j, k); c = [i, j, k]
                for d in range(3):
                    if d in skip: continue
                    for side in (0, 1):
                        t = c[:]; t[d] += -1 if side == 0 else 1
                        b = bc[d][side]
                        M[r, r] -= 1/h**2
                        if 0 <= t[d] < n[d]: M[r, ix(*t)] += 1/h**2
                        elif b == "P": t[d] %= n[d]; M[r, ix(*t)] += 1/h**2
                        elif b == "N": M[r, r] += 1/h**2
                        elif b == "D": M[r, r] -= 1/h**2
    return M

def mixedfaces():
    """Every per-direction combination of {PP,NN,DD,ND,DN} (125 sets of faces) against the dense reference, FFT and MLMG."""
    n = [6, 5, 4]; letters = ["PP", "NN", "DD", "ND", "DN"]
    combos = [(a, b, c) for a in letters for b in letters for c in letters]
    mats = {}
    for np_, mgs in ((1, 8), (2, 3)):
        for be in A.be.split():
            pre = os.path.join(A.work, f"mx_{be}_{np_}")
            rc, out = run(np_, mode="mixed1", n_cell="6 5 4", mgs=mgs, backend=be, sweep=1, out=pre)
            ok(rc == 0, f"mixed1 sweep {be} np={np_}")
            ml = lines(out, "MIXED"); ok(len(ml) == 125, f"{be} np={np_}: 125 sweep lines")
            worst = 0.0; nbad = 0
            for idx, cb in enumerate(combos):
                d = kv(ml[idx]); name = ",".join(cb)
                if d["status"] != "Ok": nbad += 1; print("NOT OK", ml[idx]); continue
                b = np.fromfile(f"{pre}_{idx}_rhs.bin"); x = np.fromfile(f"{pre}_{idx}_phi.bin")
                if name not in mats: mats[name] = dense_op(n, cb)
                M = mats[name]
                sing = all(ch in "NP" for pair in cb for ch in pair)
                if int(d["singular"]) != int(sing): nbad += 1; print("BAD singular flag", name, d["singular"])
                if sing:
                    xr = np.linalg.lstsq(M, b - b.mean(), rcond=None)[0]; xr -= xr.mean()
                else:
                    xr = np.linalg.solve(M, b)
                e = np.abs(x - xr).max()/np.abs(xr).max(); worst = max(worst, e)
                if e > 1e-10: nbad += 1; print("BAD", name, e)
            print(f"{be} np={np_}: worst relative error over 125 combinations {worst:.2e}")
            ok(nbad == 0 and worst <= 1e-10, f"{be} np={np_}: all 125 combinations agree with the dense operator (worst {worst:.2e})")

def bcdata_exact():
    """fold_boundary_data: discrete check (rhs built with explicit data ghosts; folded homogeneous solve reproduces the field) + negative control."""
    for bc in ("ND,DN,NN", "NN,NN,NN", "DD,DD,DD", "PP,ND,DN", "NN,DD,PP"):
        for be in (A.be.split() if A.be != "fft mlmg" else ("fft", "mlmg")):
            for np_ in (1, 2):
                rc, out = run(np_, mode="bcdata", part="exact", bcpairs=bc, backend=be, mgs=4)
                ok(rc == 0, f"bcdata exact bcpairs={bc} {be} np={np_}: harness checks")
                if rc != 0: print(out)
                d = kv(lines(out, "BCDATA")[0]); ok(float(d["max_err_rel"]) <= 1e-10, f"{bc} {be} np={np_}: max error {d['max_err_rel']} <= 1e-10")

def bcdata_mms():
    """Manufactured solution with exact Dirichlet / Neumann wall data folded into the rhs: second order."""
    for bc in ("DD,DD,DD", "NN,NN,NN", "ND,DN,NN", "DN,ND,DD"):
        errs = []
        for n in (16, 32, 64):
            rc, out = run(2, mode="bcdata", part="mms", bcpairs=bc, backend="fft", n_cell=f"{n} {n} {n}", mgs=32)
            ok(rc == 0, f"mms {bc} n={n}: harness checks")
            errs.append(float(kv(lines(out, "BCDATA")[0])["err_l2"]))
        orders = [math.log2(errs[i] / errs[i+1]) for i in range(2)]
        print(f"MMS {bc}: errors {errs} orders {orders}")
        ok(all(o >= 1.9 for o in orders), f"{bc}: second-order convergence {orders[0]:.2f}, {orders[1]:.2f}")
    # MLMG on the coarser grids gives the same answer as FFT
    rc, o1 = run(2, mode="bcdata", part="mms", bcpairs="ND,DN,NN", backend="fft", n_cell="32 32 32", mgs=16)
    rc2, o2 = run(2, mode="bcdata", part="mms", bcpairs="ND,DN,NN", backend="mlmg", n_cell="32 32 32", mgs=16)
    e1, e2 = float(kv(lines(o1, "BCDATA")[0])["err_l2"]), float(kv(lines(o2, "BCDATA")[0])["err_l2"])
    ok(rc == 0 and rc2 == 0 and abs(e1 - e2) <= 1e-6 * e1, f"FFT and MLMG give the same error {e1:.6e} {e2:.6e}")

def solve_dense(M, b, sing):
    if sing:
        xr = np.linalg.lstsq(M, b - b.mean(), rcond=None)[0]; return xr - xr.mean()
    return np.linalg.solve(M, b)

def rel(a, b): return np.abs(a - b).max() / max(np.abs(b).max(), 1e-300)

def pressure_bc_map():
    """D-057: FDS boundary strings of the driver sweep (one-cell y) -> backend / BC, 2-D operator, singular handling, refusals."""
    for np_ in (1, 2):
        rc, out = run(np_, mode="bcmap"); print(out)
        ok(rc == 0, f"bcmap status / selector / guard checks, np={np_}")
    n = [16, 1, 16]
    strings = [("PP,NN,PP", 91), ("NN,NN,NN", 86), ("PP,NN,NN", 25), ("DD,DD,DD", 1), ("ND,NN,NN", 17), ("DD,NN,NN", 4), ("NN,NN,ND", 4),
               ("DD,NN,ND", 3), ("ND,NN,ND", 2), ("DN,NN,DD", 1), ("NN,DD,NN", 0), ("NN,ND,PP", 0), ("DD,DN,NN", 0)]
    for s_, cnt in strings:
        bc = s_.split(","); sing = all(ch in "NP" for pair in [bc[0], bc[2]] for ch in pair)   # y never decides singularity
        M2 = dense_op(n, bc, skip=(1,)); M3 = dense_op(n, bc)
        for be in ("fft", "mlmg"):
            for np_ in (1, 2):
                pre = os.path.join(A.work, f"map_{s_.replace(',', '')}_{be}_{np_}")
                rc, out = run(np_, mode="mixed1", n_cell="16 1 16", mgs=8, backend=be, bcpairs=s_, out=pre)
                ok(rc == 0, f"{s_} ({cnt} inputs) {be} np={np_}: solved")
                d = kv(lines(out, "MIXED")[0])
                ok(d["status"] == "Ok" and d["backend"] == ("FFT" if be == "fft" else "MLMG"), f"{s_} {be}: status Ok, backend {d['backend']}")
                ok(int(d["singular"]) == int(sing), f"{s_} {be}: singular={d['singular']} (x and z faces decide)")
                ok(float(d["true_rel2"]) <= 1e-9, f"{s_} {be}: library true residual (2-D operator) {d['true_rel2']} <= 1e-9")
                b = np.fromfile(pre + "_rhs.bin"); x = np.fromfile(pre + "_phi.bin")
                e2 = rel(x, solve_dense(M2, b, sing))
                ok(e2 <= 1e-9, f"{s_} {be} np={np_}: equals the dense 2-D operator (rel err {e2:.1e})")
                # only a Dirichlet face in y distinguishes the 2-D from the full operator
                if "D" in bc[1]:
                    e3 = rel(x, solve_dense(M3, b, all(ch in 'NP' for pair in bc for ch in pair)))
                    ok(e3 > 1e-3, f"{s_} {be}: differs from the full 3-D operator with the y Dirichlet term (rel diff {e3:.2e}): the y term is dropped")
    # raw amrex::FFT::Poisson: the one-cell direction is ignored whatever its BC
    for ytype in ("PP", "NN", "DD", "ND", "DN"):
        s_ = f"NN,{ytype},NN" if ytype != "PP" else "NN,PP,NN"
        bc = s_.split(",")
        pre = os.path.join(A.work, f"raw_{ytype}")
        rc, out = run(1, mode="fftraw", n_cell="16 1 16", bcpairs=s_, out=pre); ok(rc == 0, f"raw FFT::Poisson {s_}: ran")
        b = np.fromfile(pre + "_rhs.bin"); x = np.fromfile(pre + "_phi.bin")
        sing = True   # x,z closed
        e2 = rel(x - x.mean(), solve_dense(dense_op(n, bc, skip=(1,)), b, True))
        ok(e2 <= 1e-9, f"raw FFT::Poisson y={ytype}: equals the 2-D operator, y term ignored whatever the BC (rel {e2:.1e})")
        if "D" in ytype:
            Mf = dense_op(n, bc)          # full operator: y Dirichlet term -2/h^2 keeps the matrix non-singular
            e3 = rel(x, np.linalg.solve(Mf, b))
            ok(e3 > 1e-3, f"raw FFT::Poisson y={ytype}: DISAGREES with the full 3-D operator (rel diff {e3:.2e})")
    # one-cell z / x with a Dirichlet face: raw FFT drops a term FDS keeps -> the interface must refuse
    for s_, nn in (("NN,NN,DD", [16, 16, 1]), ("DD,NN,NN", [1, 16, 16])):
        pre = os.path.join(A.work, f"raw_thin_{s_.replace(',', '')}")
        rc, out = run(1, mode="fftraw", n_cell=" ".join(map(str, nn)), bcpairs=s_, out=pre); ok(rc == 0, f"raw FFT::Poisson {s_} {nn}")
        b = np.fromfile(pre + "_rhs.bin"); x = np.fromfile(pre + "_phi.bin")
        bc = s_.split(",")
        e = rel(x, np.linalg.solve(dense_op(nn, bc), b))
        ok(e > 1e-3, f"raw FFT::Poisson {s_} on {nn}: wrong vs the full operator (rel diff {e:.2e}); the interface refuses it (bcmap)")
    for s_, nn in (("NN,NN,PP", [16, 16, 1]), ("PP,NN,NN", [1, 16, 16]), ("NN,NN,NN", [1, 16, 16])):
        pre = os.path.join(A.work, f"lib_thin_{s_.replace(',', '')}")
        rc, out = run(1, mode="mixed1", n_cell=" ".join(map(str, nn)), bcpairs=s_, backend="fft", mgs=8, out=pre); ok(rc == 0, f"library {s_} {nn}")
        b = np.fromfile(pre + "_rhs.bin"); x = np.fromfile(pre + "_phi.bin")
        e = rel(x, solve_dense(dense_op(nn, s_.split(",")), b, True))
        ok(e <= 1e-9, f"library {s_} on {nn}: one-cell x/z with N or P is exact (rel {e:.1e})")
    # singular problems use FFT::Poisson, never PoissonHybrid
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "FFTBackend.cpp")).read()
    code = "\n".join(l.split("//")[0] for l in src.splitlines())
    ok("FFT::Poisson<MultiFab>" in code and "PoissonHybrid" not in code, "FFTBackend uses FFT::Poisson only (PoissonHybrid is never selected)")

# ---- pin-aware residual check (frozen/hypre-notes.md, "Residual check") ------------------------------
def hypresid():
    """HYPRE on large singular single-level problems: the default residual_tol no longer warns although the raw residual
    exceeds it; the loose solves (negative controls) still warn; the pin row holds at most sqrt(N) times the other rows."""
    for ns in A.ns.split():
        n = f"{ns} {ns} {ns}"
        rc, out = run(2, mode="solve", n_cell=n, bc=A.bc, max_grid_size=32, backends="hypre", verbose=1); print(out)
        ok(rc == 0, f"n={ns} {A.bc}: harness exit 0")
        r = kv(lines(out, "RESULT")[0])
        raw, nop, chk, lim = (float(r[k]) for k in ("true_rel2", "rel2_nopin", "check", "limit"))
        ok(r["status"] == "Ok" and r["residual_ok"] == "1" and r["nwarn"] == "0", f"n={ns} {A.bc}: no warning at the default residual_tol 1e-12 (full {raw:.2e}, pin-excluded {nop:.2e})")
        fl = float(r["floor"])
        ok(lim == max(1e-12, fl), f"n={ns} {A.bc}: limit = max(residual_tol 1e-12, floor {fl:.2e}) = {lim:g}")
        ok(chk == nop and chk <= lim, f"n={ns} {A.bc}: checked value is the pin-excluded residual {nop:.2e} <= limit {lim:.2e}")
        N = int(ns) ** 3
        ok(float(r["pin_rel"]) <= math.sqrt(N) * nop * 1.0001 + 1e-30, f"n={ns} {A.bc}: pin-row residual {float(r['pin_rel']):.2e} <= sqrt(N)*||r_nopin|| (Cauchy-Schwarz on the zero-sum residual)")
        ok(abs(float(r["rel2_mr"]) - raw) <= 1e-3 * raw, f"n={ns} {A.bc}: mean removal of the residual changes nothing (mr {float(r['rel2_mr']):.3e} vs raw {raw:.3e})")
    n = "48 48 48"
    for tol in (1e-3, 1e-6, 1e-9):
        rc, out = run(2, mode="solve", n_cell=n, bc=A.bc, max_grid_size=32, backends="hypre", verbose=1, tol_rel=tol)
        r = kv(lines(out, "RESULT")[0])
        ok(r["residual_ok"] == "0" and r["nwarn"] == "1", f"negative control tol_rel={tol:g}: warns (checked {float(r['check']):.2e} > limit {float(r['limit']):.2e})")
    rc, out = run(2, mode="solve", n_cell=n, bc=A.bc, max_grid_size=32, backends="hypre", verbose=1, max_iter=4)
    r = kv(lines(out, "RESULT")[0])
    ok(r["residual_ok"] == "0" and r["status"] == "NotConverged", f"negative control max_iter=4: NotConverged and the residual check warns ({float(r['check']):.2e})")

def reslimit():
    """Size-scaled residual limit (FR-031, A-67) for every backend, singular (Neumann) and non-singular (Dirichlet) problems:
    limit = max(1e-12, floor), floor > 0 for both kinds and growing like N^2, unperturbed solves do not warn, loose solves do."""
    flo = {}
    for bc in ("neumann", "dirichlet"):
        for ns in A.ns.split():
            rc, out = run(2, mode="solve", n_cell=f"{ns} {ns} {ns}", bc=bc, max_grid_size=32, backends="fft mlmg hypre", verbose=1); print(out)
            ok(rc == 0, f"{bc} n={ns}: harness exit 0")
            for l in lines(out, "RESULT"):
                r = kv(l); be = l.split()[1]
                fl, lim, chk = float(r["floor"]), float(r["limit"]), float(r["check"])
                ok(r["status"] == "Ok" and r["residual_ok"] == "1" and r["nwarn"] == "0", f"{bc} n={ns} {be}: Ok, no warning (check {chk:.2e}, limit {lim:.2e})")
                ok(fl > 0.0 and lim == max(1e-12, fl), f"{bc} n={ns} {be}: floor {fl:.2e} > 0 and limit {lim:.3e} = max(1e-12, floor)")
                ok(chk <= lim, f"{bc} n={ns} {be}: checked value {chk:.2e} <= limit")
                flo[(bc, be, int(ns))] = fl
    ns_ = sorted(int(x) for x in A.ns.split())
    for bc in ("neumann", "dirichlet"):
        for be in ("fft", "mlmg", "hypre"):
            if (bc, be, ns_[0]) not in flo or (bc, be, ns_[-1]) not in flo: continue
            ratio = flo[(bc, be, ns_[-1])] / flo[(bc, be, ns_[0])]; want = (ns_[-1] / ns_[0]) ** 2
            ok(0.8 * want <= ratio <= 1.25 * want, f"{bc} {be}: floor grows like N^2 from n={ns_[0]} to {ns_[-1]} (ratio {ratio:.2f}, N^2 ratio {want:.2f})")
    # Negative controls on the non-singular (Dirichlet) problem: the scaled limit does not hide an under-converged solve.
    for be in ("mlmg", "hypre"):
        rc, out = run(2, mode="solve", n_cell="48 48 48", bc="dirichlet", max_grid_size=32, backends=be, verbose=1, tol_rel=1e-6)
        r = kv(lines(out, "RESULT")[0])
        ok(r["residual_ok"] == "0" and r["nwarn"] == "1" and float(r["check"]) > float(r["limit"]), f"dirichlet {be} tol_rel=1e-6: warns (check {float(r['check']):.2e} > limit {float(r['limit']):.2e})")
        rc, out = run(2, mode="solve", n_cell="48 48 48", bc="dirichlet", max_grid_size=32, backends=be, verbose=1, max_iter=4)
        r = kv(lines(out, "RESULT")[0])
        ok(r["residual_ok"] == "0" and float(r["check"]) > float(r["limit"]), f"dirichlet {be} max_iter=4: residual check warns (status {r['status']}, check {float(r['check']):.2e})")

def rescheck():
    """Residual check on the final H of a solve and on perturbed copies; synthetic sums (harness mode rescheck)."""
    for be, nn in (("hypre", A.n), ("mlmg", "32 32 32"), ("fft", "32 32 32")):
        rc, out = run(2, mode="rescheck", n_cell=nn, bc=A.bc, max_grid_size=16, backend=be); print(out)
        ok(rc == 0, f"rescheck {be} {A.bc} n=[{nn}]: all harness checks pass")
        base = kv(lines(out, "RESCHECK")[0])
        if be != "hypre":
            ok(base["pin_rel"] == "0" and base["nopin"] == base["raw"], f"{be}: no pin applied, pin-excluded residual equals the raw one")

def hypresidcomp():
    """Composite: default residual_tol, HYPRE and MLMG do not warn; loose solves warn."""
    for bc, nlev, ratio, extra in (("neumann", 2, 2, {}), ("neumann", 3, 2, {}), ("periodic", 2, 4, {}), ("periodic", 3, 2, {}), ("neumann", 2, 2, {"plane2d": 1})):
        kw = dict(mode="hypre_resid", n=16, nlev=nlev, ratio=ratio, mgs=8, bc=bc, **extra)
        rc, out = run(2, **kw); print(out)
        ok(rc == 0, f"hypre_resid {bc} nlev={nlev} ratio={ratio} {extra}: no warning at the default residual_tol (both backends)")
        rc, out = run(2, tol_rel=1e-4, expect_warn=1, **kw)
        ok(rc == 0, f"hypre_resid {bc} nlev={nlev} ratio={ratio} {extra}: negative control tol_rel=1e-4 warns (both backends)")


{"selector": selector, "exactsum": exactsum, "fftmlmg": fftmlmg, "frozen": frozen, "decomp": decomp,
 "repeat": repeat, "singular": singular, "ulmat": ulmat, "ulmatgauge": ulmatgauge, "meankind": meankind,
 "compconv": compconv, "compfull": compfull, "compdecomp": compdecomp, "comp3": comp3, "compgrad": compgrad, "compshape": compshape, "compmixed": compmixed,
 "compns2d": compns2d, "compsel": compsel, "compws": compws, "compgauge": compgauge, "meankind_uniform": meankind_uniform, "compgaugedecomp": compgaugedecomp, "trigger1": trigger1, "comptrigger": comptrigger, "fftcache": fftcache, "mixedfaces": mixedfaces, "bcdata_exact": bcdata_exact, "bcdata_mms": bcdata_mms, "pressure_bc_map": pressure_bc_map,
 "hypsingle": hypsingle, "hypci": hypci, "hypcomp": hypcomp, "hypcache": hypcache, "hypresid": hypresid, "rescheck": rescheck, "reslimit": reslimit, "hypresidcomp": hypresidcomp}[A.cmd]()
if fails:
    print("FAILED:", *fails, sep="\n  "); sys.exit(1)
print("ALL PASS")
