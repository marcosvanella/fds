#!/usr/bin/env python3
"""Driver of the bitwise tests of the host translations in FdsPressureLoops.cpp (pressure-area loops of FDS pres.f90).

  make-ref : assemble the Fortran reference (template + verbatim upstream loops, hash-pinned), compile it
             (gfortran: -O2 -ffp-contract=off -fno-fast-math; ifx/ifort: -O2 -fp-model=precise), run it; it writes the case files into <work>/cases.
  check    : run pb_fds_loops (or a mutant build of it) on the case files; exit 0 iff every element is bitwise equal.
             --expect-fail inverts the result (mutation check: the test must notice the mutant).
  check-fds: unpack the archive of case files written by the real FDS (frozen/fds_loops_cases, made with
             frozen/fds_loops_dump_hook.py) and run pb_fds_loops on them: the C++ functions against FDS's own arrays.
  check-fds-sets: unpack an archive with one sub-directory of FDS-written case files per input (frozen/fds_loops_cases/fds_loops_cases_r5.tar.gz,
             made with the same hook from frozen/fds_cases/ldump_*.fds) and run pb_fds_loops on every set; all must be bitwise equal.
  census   : read the same archive and check that every input really exercises what it was made for (the boundary codes, wall kinds, wind,
             ramp and eddy classes listed in CENSUS below); a changed input or hook that silently lost the feature fails here.
  drift    : the generator step on a copy of pres.f90 with one inserted line must REFUSE (the pinned line ranges moved).
"""
import argparse, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, **kw):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, **kw)
    return r.returncode, r.stdout


def read_records(path):
    """Stream records of tests/fds_loops/ref_loops.f90.in: name(16) kind(i4) rank(i4) dims(i4*rank) data; kind 2 = real(8), else int(4)."""
    import struct
    b = open(path, "rb").read()
    p, c = 0, {}
    while p < len(b):
        nm = b[p:p + 16].decode().strip(); p += 16
        kind, rank = struct.unpack("<ii", b[p:p + 8]); p += 8
        dims = struct.unpack("<%di" % rank, b[p:p + 4 * rank]); p += 4 * rank
        n = 1
        for d in dims:
            n *= d
        if kind == 2:
            v = struct.unpack("<%dd" % n, b[p:p + 8 * n]); p += 8 * n
        else:
            v = struct.unpack("<%di" % n, b[p:p + 4 * n]); p += 4 * n
        c[nm] = (dims, v)
    return c


TWENTY_EPS = 20.0 * 2.220446049250313e-16
OPEN, SOLID, DIRICHLET = 2, 1, 1


def census_set(d):
    """Facts about one set of FDS-written files: union over its files."""
    f = {"codes": set(), "open_ior": set(), "solid_dir_ww1": set(), "wind": False, "jbar": set(), "tbegin": set(), "prhs_i": set(),
         "ramp_t_eq": 0, "ramp_t_ne": 0, "noramp": 0, "eddy_walls": 0, "eddy_nonzero": 0, "z_ghost_k": set(), "bt": set(), "kbar": set()}
    for n in sorted(os.listdir(d)):
        if not n.endswith(".bin"):
            continue
        c = read_records(os.path.join(d, n))
        g = lambda k: c[k][1]
        if n.startswith("HBC"):
            f["codes"].add((g("LBC")[0], g("MBC")[0], g("NBC")[0]))
        if n.startswith("L1211"):
            f["prhs_i"].add((c["PRHS_IN"][0][0], g("IBAR")[0]))
        if n.startswith("L1209"):
            ior, bt, pbc, w1 = g("W_IOR"), g("W_BT"), g("W_PBC"), g("WALL_WORK1")
            f["wind"] = f["wind"] or g("OPEN_WIND")[0] == 1
            f["jbar"].add(g("JBAR")[0]); f["kbar"].add(g("KBAR")[0]); f["tbegin"].add(g("T_BEGIN")[0])
            f["bt"].update(bt)
            for i in range(len(ior)):
                if bt[i] == OPEN and pbc[i] == DIRICHLET:
                    f["open_ior"].add(ior[i])
                    if abs(ior[i]) == 3:
                        f["z_ghost_k"].add(g("W_K")[i] - (g("KBAR")[0] if ior[i] < 0 else 0))   # 0 for KK=0, 1 for KK=KBP1
                    ramp = g("WV_RAMP")[i]
                    eq = abs(g("W_TIGN")[i] - g("T_BEGIN")[0]) <= TWENTY_EPS
                    if ramp >= 1 and eq: f["ramp_t_eq"] += 1
                    elif ramp >= 1: f["ramp_t_ne"] += 1
                    else: f["noramp"] += 1
                    if g("WV_NEDDY")[i] > 0:
                        f["eddy_walls"] += 1
                        if "WV_EDDYV" in c and g("WV_EDDYV")[i] != 0.0: f["eddy_nonzero"] += 1
                if bt[i] not in (OPEN, 6) and pbc[i] == DIRICHLET and w1[i] != 0.0:
                    f["solid_dir_ww1"].add(ior[i])
    return f


# name of the set -> (description, predicate on the facts)
CENSUS = {
    "2d":     ("TWO_D, JBAR=1, mirror-type y walls",           lambda f: f["jbar"] == {1} and f["codes"] == {(4, 3, 4)}),
    "perx":   ("periodic in x: LBC=0, PRHS leading dimension IBAR",  lambda f: f["codes"] == {(0, 3, 4)} and all(a == b for a, b in f["prhs_i"])),
    "peryz":  ("periodic in y and z: MBC=NBC=0",               lambda f: f["codes"] == {(3, 0, 0)}),
    "dirx":   ("Dirichlet low-x face: LBC=2, SOLID-Dirichlet walls with WALL_WORK1 /= 0 at IOR=1, open IOR=1",
               lambda f: f["codes"] == {(2, 3, 3)} and f["solid_dir_ww1"] == {1} and f["open_ior"] == {1}),
    "dirxhi": ("Dirichlet high-x face: LBC=4, SOLID-Dirichlet with WALL_WORK1 /= 0 at IOR=-1, open IOR=-1",
               lambda f: f["codes"] == {(4, 3, 3)} and f["solid_dir_ww1"] == {-1} and f["open_ior"] == {-1}),
    "ylow":   ("Dirichlet low-y and low-z: MBC=2, NBC=2, open IOR=2 and 3",
               lambda f: f["codes"] == {(3, 2, 2)} and f["open_ior"] == {2, 3}),
    "openyz": ("open y faces, both z faces and +x: IOR +-2, 3, -3, -1 with wind; z walls at KK=0 and KK=KBP1",
               lambda f: f["wind"] and f["open_ior"] == {2, -2, 3, -3, -1} and f["z_ghost_k"] == {0, 1}),
    "ramp":   ("T_BEGIN /= 0; ramped open walls with T_IGN = T_BEGIN and with T_IGN /= T_BEGIN, and open walls without a ramp",
               lambda f: f["tbegin"] == {0.3} and f["ramp_t_eq"] > 0 and f["ramp_t_ne"] > 0 and f["noramp"] > 0),
    "eddy":   ("synthetic eddy inflow with wind: N_EDDY > 0 walls with non-zero eddy values",
               lambda f: f["wind"] and f["eddy_walls"] > 0 and f["eddy_nonzero"] > 0),
}


def unpack_sets(archive, work):
    import tarfile
    d = os.path.join(work, "fds_sets")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    with tarfile.open(archive) as t:
        t.extractall(d)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["make-ref", "check", "check-fds", "check-fds-sets", "census", "drift"])
    ap.add_argument("--pres", default=os.path.join(HERE, "..", "..", "pres.f90"))
    ap.add_argument("--git-repo", default=None)
    ap.add_argument("--fortran", default=os.environ.get("FC", "gfortran"))
    ap.add_argument("--work", required=True)
    ap.add_argument("--exe", default=None)
    ap.add_argument("--archive", default=None, help="tar of FDS-written case files (check-fds)")
    ap.add_argument("--exe-arg", action="append", default=None, help="extra argument for the executable (check-fds)")
    ap.add_argument("--expect-fail", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    gen = [sys.executable, os.path.join(HERE, "fds_loops", "fds_loops_gen.py"), "--template", os.path.join(HERE, "fds_loops", "ref_loops.f90.in")]
    if a.action == "make-ref":
        f90 = os.path.join(a.work, "ref_loops.f90")
        cmd = gen + ["--pres", a.pres, "--out", f90] + (["--git-repo", a.git_repo] if a.git_repo else [])
        rc, out = run(cmd)
        print(out, end="")
        if rc:
            return 1
        exe = os.path.join(a.work, "ref_loops")
        # D-079: Intel ifx/ifort reject GNU-only -ffp-contract=off/-fno-fast-math; use -fp-model=precise there.
        fc_base = os.path.basename(a.fortran)
        if "ifx" in fc_base or "ifort" in fc_base:
            fcflags = ["-O2", "-fp-model=precise"]
        else:
            fcflags = ["-O2", "-ffp-contract=off", "-fno-fast-math"]
        rc, out = run([a.fortran] + fcflags + ["-ffree-line-length-none", "-o", exe, f90], cwd=a.work)
        print(out, end="")
        if rc:
            return 1
        cases = os.path.join(a.work, "cases")
        shutil.rmtree(cases, ignore_errors=True)
        os.makedirs(cases)
        rc, out = run([exe, cases])
        print(out, end="")
        n = len([f for f in os.listdir(cases) if f.endswith(".bin")])
        print("reference wrote %d case files" % n)
        return 0 if rc == 0 and n > 0 else 1
    if a.action == "check":
        if not a.exe:
            sys.exit("--exe required")
        rc, out = run([a.exe, "dir=" + os.path.join(a.work, "cases")])
        print(out, end="")
        ok = rc == 0 and "SUMMARY" in out and "PASS" in out
        if a.expect_fail:
            if ok:
                print("MUTANT NOT CAUGHT: the bitwise test passed on a mutant")
                return 1
            print("MUTANT CAUGHT")
            return 0
        return 0 if ok else 1
    if a.action == "check-fds":
        import tarfile
        if not a.exe or not a.archive:
            sys.exit("--exe and --archive required")
        d = os.path.join(a.work, "fds_cases")
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        with tarfile.open(a.archive) as t:
            t.extractall(d)
        names = sorted(os.listdir(d))
        n1209 = [n for n in names if n.startswith("L1209_")]
        print("L1209 files in the archive:", len(n1209))
        if len(n1209) < 5:
            print("FAIL: the archive must hold the L1209 files (open-wall branches against real FDS)")
            return 1
        rc, out = run([a.exe, "dir=" + d] + (a.exe_arg or []))
        print(out, end="")
        good = rc == 0 and "SUMMARY" in out and "PASS" in out and "LOOP L1209 " in out and "LOOP L1209              cases   0" not in out
        if a.expect_fail:
            if good:
                print("NEGATIVE CONTROL NOT CAUGHT: the FDS-file comparison passed with the wrong velocity lower bound")
                return 1
            print("NEGATIVE CONTROL CAUGHT")
            return 0
        return 0 if good else 1
    if a.action == "check-fds-sets":
        if not a.exe or not a.archive:
            sys.exit("--exe and --archive required")
        d = unpack_sets(a.archive, a.work)
        sets = sorted(n for n in os.listdir(d) if os.path.isdir(os.path.join(d, n)))
        print("sets:", " ".join(sets))
        missing = [n for n in CENSUS if n not in sets]
        if missing:
            print("FAIL: sets missing from the archive:", missing)
            return 1
        allgood, anybad = True, False
        for n in sets:
            rc, out = run([a.exe, "dir=" + os.path.join(d, n)] + (a.exe_arg or []))
            ok = rc == 0 and "SUMMARY" in out and "PASS" in out
            print("%-8s %s  %s" % (n, "bitwise equal" if ok else "DIFFERENT", [l for l in out.splitlines() if l.startswith("SUMMARY")][:1]))
            if not ok:
                print(out)
            allgood = allgood and ok
            anybad = anybad or not ok
        if a.expect_fail:
            if not anybad:
                print("NEGATIVE CONTROL / MUTANT NOT CAUGHT on the FDS-written sets")
                return 1
            print("CAUGHT by the FDS-written sets")
            return 0
        return 0 if allgood else 1
    if a.action == "census":
        if not a.archive:
            sys.exit("--archive required")
        d = unpack_sets(a.archive, a.work)
        bad = 0
        for n, (what, pred) in CENSUS.items():
            f = census_set(os.path.join(d, n))
            ok = pred(f)
            print("%-8s %s: %s" % (n, "OK  " if ok else "FAIL", what))
            if not ok:
                print("         facts:", {k: (sorted(v) if isinstance(v, set) else v) for k, v in f.items()})
            bad += 0 if ok else 1
        return 1 if bad else 0
    if a.action == "drift":
        bad = os.path.join(a.work, "pres_drifted.f90")
        with open(a.pres, errors="replace") as f:
            txt = f.read()
        with open(bad, "w") as f:
            f.write("! one inserted line\n" + txt)
        rc, out = run(gen + ["--pres", bad, "--out", os.path.join(a.work, "drift_ref.f90")])   # no git fallback on purpose
        print(out, end="")
        if rc != 0 and "drifted" in out:
            print("DRIFT DETECTED (expected)")
            return 0
        print("DRIFT NOT DETECTED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
