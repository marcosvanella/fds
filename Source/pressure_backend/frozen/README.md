# frozen/

M1 uses synthetic frozen cases generated on the fly: `pb_harness mode=gen n_cell="..." bc=... out=<prefix>`
writes the RHS (`<prefix>_rhs.bin`) and the FFT reference H (`<prefix>_H.bin`) as raw float64, x fastest, and
the `frozen` CTest solves the stored RHS with MLMG and FFT and compares against the stored H. The RHS is the P2
smooth non-separable field (`rhs_func` in `harness/main.cpp`), made compatible once with the common-layer mean
removal. The committed store (manifest, small stored cases, reviewed generator script) is milestone M3;
Dirichlet and FDS-derived cases are M2 and M4.

`mean-removal-vs-fds.md` is the note on FDS mean removal and gauge. It is backed by `stretched_study.py` (numpy
reference against FDS dumps), `stretched_study_results.txt` (its output for the cases in the note),
`fds_dump_hook.py` (scratch-only edit of a copy of `pres.f90`; never applied to the reference tree) and the two small
FDS inputs in `fds_cases/`.

`hypre-notes.md` documents the HYPRE assembled-matrix backend (operator, C/F treatment, pin, tolerance, agreement and iteration
numbers, timing against MLMG, limitations). `mixed-nd-hierarchy-note.md` reports where the fine-level excess error of mixed
Neumann/Dirichlet directions sits on a hierarchy (measured with `pb_harness mode=err_map`).
