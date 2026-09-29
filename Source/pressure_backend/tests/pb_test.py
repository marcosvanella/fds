#!/usr/bin/env python3
"""CTest driver for the pressure backend harness. Parses the harness RESULT/CMP/CHECK lines.
Subcommands: selector, exactsum, fftmlmg, frozen, decomp, repeat, singular. Exit 0 = pass."""
import argparse, os, re, subprocess, sys

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

{"selector": selector, "exactsum": exactsum, "fftmlmg": fftmlmg, "frozen": frozen, "decomp": decomp,
 "repeat": repeat, "singular": singular}[A.cmd]()
if fails:
    print("FAILED:", *fails, sep="\n  "); sys.exit(1)
print("ALL PASS")
