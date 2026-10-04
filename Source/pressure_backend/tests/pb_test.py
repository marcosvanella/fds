#!/usr/bin/env python3
"""CTest driver for the pressure backend harness. Parses the harness RESULT/CMP/CHECK lines.
Subcommands: selector, exactsum, fftmlmg, frozen, decomp, repeat, singular, ulmat, ulmatgauge, meankind.
Exit 0 = pass. ulmat/ulmatgauge/meankind compare against an independent numpy computation that follows FDS
ULMAT (volume-scaled rows, arithmetic mean removal of F and of X, identity-pin reduced system, rho*volume gauge)."""
import argparse, os, re, subprocess, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("cmd")
ap.add_argument("--harness", required=True)
ap.add_argument("--mpiexec", required=True)
ap.add_argument("--work", required=True)
ap.add_argument("--bc", default="neumann")
ap.add_argument("--n", default="64 64 64")
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
    rc, out = run(np_, **base(max_grid_size=mgs, backends="fft mlmg", out=pre))
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
            for b in ("fft", "mlmg"):
                r = diffrel(f"{pre}_{b}.bin", f"{ref_pre}_{b}.bin")
                bit = h[b] == ref_h[b]
                print(f"DECOMP {A.bc} np={np_} mgs={mgs} {b}: rel_l2 vs np1/mgs32 = {r:.3e} bitwise={'yes' if bit else 'no'}")
                ok(r <= e, f"{b} np={np_} mgs={mgs} within eps_H")
            ok(pins == ref_pins, f"pin and removed mean bitwise identical np={np_} mgs={mgs}")

def repeat():
    hs = []
    for i in range(3):
        _, h, _ = solve_files(4, 16, f"rep{i}"); hs.append(h)
    for b in ("fft", "mlmg"):
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
                                out=os.path.join(A.work, tag)))
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
    rc, out = run(2, **base(max_grid_size=4, backends="fft mlmg", rhs_file=pre + "_rhs.bin", rhs_offset=0.4, gauge_rho=1, out=outp))
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

{"selector": selector, "exactsum": exactsum, "fftmlmg": fftmlmg, "frozen": frozen, "decomp": decomp,
 "repeat": repeat, "singular": singular, "ulmat": ulmat, "ulmatgauge": ulmatgauge, "meankind": meankind}[A.cmd]()
if fails:
    print("FAILED:", *fails, sep="\n  "); sys.exit(1)
print("ALL PASS")
