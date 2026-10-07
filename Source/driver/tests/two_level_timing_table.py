#!/usr/bin/env python3
# Tables of notes/two-level-timing.md from timing.csv (run_two_level_timing.sh). Times in ms per step. "host stages" = t_loop - t_pressure - t_fds_outputs - avg-down - cf ghost
# fill - flux overwrite: the FDS stages outside the pressure scheme, same-level ghost exchange included. FillBoundary, fill_omesh and the FDS bc routines run inside and outside the
# pressure scheme and are listed, not subtracted. Model (ESTIMATE): two-level = single coarse + (fine patch cells /
# uniform fine cells) x uniform fine, i.e. one dt for both levels and no AMR machinery cost.
import csv, re, sys
rows = {r["tag"]: r for r in csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "timing.csv"))}
ms = lambda r, k: 1e3 * float(r[k]) / float(r["steps"])
slots = ["ghost_fill_boundary", "fill_omesh", "match_periodic_faces", "fds_bc_routines", "avg_down", "cf_ghost_fill", "flux_overwrite"]
print("| run | np | step | pressure scheme | of which solve | avg-down | cf ghost fill | flux overwrite | FillBoundary | fill_omesh | FDS bc routines | outputs | host stages |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t, r in rows.items():
    v = {k: ms(r, k) for k in ["t_loop", "t_pressure", "t_poisson_solve", "t_fds_outputs"] + slots}
    rest = v["t_loop"] - v["t_pressure"] - v["t_fds_outputs"] - v["avg_down"] - v["cf_ghost_fill"] - v["flux_overwrite"]
    print("| %s | %s | %.2f | %.2f | %.2f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f |" % (t, r["np"], v["t_loop"], v["t_pressure"], v["t_poisson_solve"], v["avg_down"],
          v["cf_ghost_fill"], v["flux_overwrite"], v["ghost_fill_boundary"], v["fill_omesh"], v["fds_bc_routines"], v["t_fds_outputs"], rest))
print()
print("| case | np | coarse cells | fine patch cells | single coarse | two-level run | model ESTIMATE | run / model | REF (patch = domain) | uniform fine single-level | AMR overhead per fine-patch cell | uniform fine cost per cell | two-level / uniform fine |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in rows:
    if not t.endswith("_two"): continue
    b = t[:-4]; case, np_ = b.rsplit("_np", 1); n = int(re.search(r"(\d+)$", case).group(1))
    s, w, f, rf = (ms(rows[b + x], "t_loop") for x in ("_single", "_two", "_fine", "_ref"))
    nf = (n // 2) ** 2 * 4; nu = (2 * n) ** 2
    model = s + f * nf / nu
    print("| %s | %s | %d | %d | %.2f | %.2f | %.2f | %.1fx | %.2f | %.2f | %.2f us | %.2f us | %.2fx |" % (case, np_, n * n, nf, s, w, model, w / model, rf, f, 1e3 * (w - s) / nf - 1e3 * f / nu, 1e3 * f / nu, w / f))
