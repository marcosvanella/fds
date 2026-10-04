# Generator kernel entries for L1211 and L1207 (proposal for the generator owners)

Status: proposal only. No file in the generator worktree was edited. The probe ran against a scratch copy of the generator tree with the
command `s5gen.py --repo <src> --commit 36975d7 --markers <probe sidecar> --out <dir> --golden <file>`.

The independent CPU reference implementations are in `Source/pressure_backend/FdsPressureLoops.cpp` (`pres_p_from_h`,
`pres_compute_rhs_div`), tested bitwise against the verbatim upstream lines (see `Source/pressure_backend/frozen/fds-loops-notes.md`). A
generated kernel can be compared with the same case files.

## L1207 (pres.f90:758-764), accepted as is

```toml
[[kernel]]
name = "pres_p_from_h"
file = "pres.f90"
routine = "PRESSURE_SOLVER_CHECK_RESIDUALS"
lines = [758, 764]
anchor_line = 758
anchor = "DOK=0,KBP1"
layout = "exact"
```

The generator accepted the entry and
produced the correct bounds, including `RHOP(-1:IBAR+2, ...)`. Trap for the review: the box includes ghost cells in every direction.

## L1211 (pres.f90:250-260), needs a generator feature

```toml
[[kernel]]
name = "pres_rhs_div"
file = "pres.f90"
routine = "PRESSURE_SOLVER_COMPUTE_RHS"
lines = [250, 260]
anchor_line = 250
anchor = "DOK=1,KBAR"
private = ["TRM1", "TRM2", "TRM3", "TRM4"]
```

With the default layout the generated PRHS is `PRHS(0:IBAR+1, ...)`. FDS allocates `PRHS(ITRN,JTRN,KTRN)` where ITRN is IBP1 (IBAR when x
is periodic), JTRN is 1 when JBAR=1, so the strides and the index origin differ from the FDS array. Two readings:

- For the AMR route the natural ghosted layout (0:IBAR+1) is the right one and the default output is usable. The sum order generated,
  `((TRM1+TRM2)+TRM3)+TRM4`, matches the Fortran.
- To run the kernel on the FDS array itself (needed to compare against an FDS dump), `layout = "exact"` gives `PRHS(ITRN,JTRN,KTRN)`, but
  ITRN, JTRN and KTRN are then neither declared nor passed. Adding them to `args` fails with "derived argument set differs from the pinned
  set: only pinned ITRN JTRN KTRN". A `policy.arrays` entry with `dims = "ITRN,JTRN,KTRN"` fails the same way: only the rank-2 wall-kernel
  path (CONNECTED_ZONES) turns extent symbols into scalar arguments.

What the generator owners need to add: extent symbols of an `exact` allocation (here ITRN, JTRN, KTRN) become integer arguments for
cell-loop kernels too, and PRHS gets an index origin of 1. The same "layout contract" is what L1212-L1214 (transposed PRHS) need.

## What the CPU reference already covers

L1211 (IPS 1/4/7, non-cylindrical) with ITRN either IBAR+1 or IBAR; L1207. L1210 (cylindrical branch) and L1212-L1214 are not covered.
