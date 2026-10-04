#!/usr/bin/env python3
"""Assemble the Fortran reference of the pressure-loop tests: the template with the VERBATIM upstream loop text.

The loop texts are copied by line range from FireX 36975d7 Source/pres.f90 (the survey source of docs/inventory). Each range is pinned by the
text of its first and last statement and by a SHA-256 of the whole range, so a drifted source stops the test instead of testing something else.
If the given file differs from the pinned text the script tries `git show 36975d7:Source/pres.f90` in --git-repo.
"""
import argparse, hashlib, subprocess, sys

# id: (first line, last line, first statement, last statement, sha256 of the lines joined with \n plus a final \n, first 16 hex digits)
RANGES = {
    "L1211": (250, 260, "DO K=1,KBAR", "ENDDO", "cbde46e0b48c72a3"),
    "L1207": (758, 764, "DO K=0,KBP1", "ENDDO", "59c5804f1242b021"),
    "L1220": (450, 462, "DO K=1,KBAR", "ENDDO", "b9603f38f0d679e1"),
    "L1221": (466, 477, "DO K=1,KBAR", "ENDDO", "0ac836d000b5581b"),
    "L1222": (481, 492, "DO J=1,JBAR", "ENDDO", "cfbc50c0b2fc4f55"),
    "L1209": (65, 228, "WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS", "ENDDO WALL_CELL_LOOP", "209508972d0bd20f"),
}


def extract(lines, lid):
    a, b, first, last, sha = RANGES[lid]
    seg = lines[a - 1:b]
    txt = "\n".join(seg) + "\n"
    ok = (len(seg) == b - a + 1 and seg[0].strip() == first and seg[-1].strip() == last and hashlib.sha256(txt.encode()).hexdigest().startswith(sha))
    return ok, seg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pres", required=True, help="Source/pres.f90 of the worktree")
    ap.add_argument("--template", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--git-repo", default=None, help="repository to try 'git show 36975d7:Source/pres.f90' in when the file has drifted")
    a = ap.parse_args()
    lines = open(a.pres, errors="replace").read().split("\n")
    git_lines = None
    texts = {}
    for lid in RANGES:
        ok, seg = extract(lines, lid)
        if not ok and a.git_repo:
            if git_lines is None:
                git_lines = subprocess.run(["git", "-C", a.git_repo, "show", "36975d7:Source/pres.f90"], capture_output=True, text=True, check=True).stdout.split("\n")
            ok, seg = extract(git_lines, lid)
        if not ok:
            sys.exit("fds_loops_gen: %s: lines %d-%d of %s are not the pinned survey text (FireX 36975d7); the loop drifted, review the translation" % (lid, RANGES[lid][0], RANGES[lid][1], a.pres))
        texts[lid] = "\n".join(seg)
    out = []
    used = set()
    for ln in open(a.template).read().split("\n"):
        s = ln.strip()
        if s.startswith("!@@") and s.endswith("@@"):
            lid = s[3:-2]
            out.append("! ---- verbatim %s: pres.f90:%d-%d (FireX 36975d7) ----" % (lid, RANGES[lid][0], RANGES[lid][1]))
            out.append(texts[lid])
            out.append("! ---- end verbatim %s ----" % lid)
            used.add(lid)
        else:
            out.append(ln)
    if used != set(RANGES):
        sys.exit("fds_loops_gen: template does not use %s" % sorted(set(RANGES) - used))
    open(a.out, "w").write("\n".join(out))


if __name__ == "__main__":
    main()
