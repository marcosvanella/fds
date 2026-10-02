# Sign-off request: patch 0003 (CHECK_MASS_DENSITY species clipping loop, two-pass split)

Reviewer: Species & Combustion Lead. Decision D-051, item 1 (family "`DELTA_RHO_ZZ` scatter, species and combustion loops").

**What changes.** `Source/mass.f90`, `CHECK_MASS_DENSITY`, the species loop (mass.f90:868-925 at FireX 36975d765f). The per-cell
computation of the neighbour amounts stays as it is. The seven additions into `DELTA_RHO_ZZ` (mass.f90:916-922) move into a second loop
over the cells in the same order (K, J, I) and read the amounts from two local scratch arrays filled by the first loop. Text of the change:
`0003-mass-check-density-two-pass.patch` (35 lines added, 8 removed, one routine).

**Why.** The GPU generator can run the first loop on the device only if nothing in it writes outside its own cell. The scatter is
inherently serial (seven cells write to one element, and the order of floating-point additions is part of the result), so it stays on
the host in the original order. Without the split this loop cannot be offloaded and stays on the CPU.

**Risk.**
- Results: none intended. Every added amount is the same expression, rounded the same way, added in the same order. Verified bitwise
  (patch description, checks 1 and 2).
- Memory: about 60 bytes per cell of the mesh, allocated per call and released on return (seven REAL(EB) amounts and one INTEGER flag).
  If this is too much, the arrays can move to the mesh work arrays; that changes the patch, not the idea.
- Time: one extra pass over the cells per tracked species (the flag store is made for every cell). Not measured on a large case.
- Not covered: MPI runs with more than one process, the full V&V baseline, and any case where clipping fires in many cells at once.

**What to check.**
1. The scatter order argument: for any element of `DELTA_RHO_ZZ`, the additions arrive in the same order as before (pass 2 visits the
   cells in the order of the original loop and applies the seven statements in the original order).
2. `CLIP_RHO_ZZ(N)` is now set after pass 1 from the flag array (flag 1 means "clipped", including cells with `SUM_MASS_N<=TWO_EPSILON_EB`
   that do not scatter), as the old line `CLIP_RHO_ZZ(N) = .TRUE.` did for every clipped, non-solid cell.
3. That the density loop above (mass.f90:799-846, same structure) can stay as it is for now; if the same split is wanted there, it is a
   second patch.
4. That two local allocatable arrays in this routine are acceptable in the FDS coding style, or which work arrays to use instead.
