#!/usr/bin/env python3
"""Driver of the bitwise tests of the host translations in FdsPressureLoops.cpp (pressure-area loops of FDS pres.f90).

  make-ref : assemble the Fortran reference (template + verbatim upstream loops, hash-pinned), compile it with gfortran
             (-O2 -ffp-contract=off, no fast-math), run it; it writes the case files into <work>/cases.
  check    : run pb_fds_loops (or a mutant build of it) on the case files; exit 0 iff every element is bitwise equal.
             --expect-fail inverts the result (mutation check: the test must notice the mutant).
  check-fds: unpack the archive of case files written by the real FDS (frozen/fds_loops_cases, made with
             frozen/fds_loops_dump_hook.py) and run pb_fds_loops on them: the C++ functions against FDS's own arrays.
  drift    : the generator step on a copy of pres.f90 with one inserted line must REFUSE (the pinned line ranges moved).
"""
import argparse, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, **kw):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, **kw)
    return r.returncode, r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["make-ref", "check", "check-fds", "drift"])
    ap.add_argument("--pres", default=os.path.join(HERE, "..", "..", "pres.f90"))
    ap.add_argument("--git-repo", default=None)
    ap.add_argument("--fortran", default=os.environ.get("FC", "gfortran"))
    ap.add_argument("--work", required=True)
    ap.add_argument("--exe", default=None)
    ap.add_argument("--archive", default=None, help="tar of FDS-written case files (check-fds)")
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
        rc, out = run([a.fortran, "-O2", "-ffp-contract=off", "-fno-fast-math", "-ffree-line-length-none", "-o", exe, f90], cwd=a.work)
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
        rc, out = run([a.exe, "dir=" + d])
        print(out, end="")
        return 0 if rc == 0 and "SUMMARY" in out and "PASS" in out else 1
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
