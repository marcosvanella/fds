# Making FDS friendlier to GPUs: a primer for FDS developers

Audience: a Fortran CFD developer who knows FDS well but not AMReX or GPU programming. Status: draft for the owner; it informs the K1/K2 choice in ADR-001 and decides nothing. Numbers come from `docs/inventory/gpu_callgraph_survey.md` (FireX `36975d7`). Line numbers are for that commit.

**Terms used below.** *Kernel*: one loop nest that runs on the GPU, with thousands of threads each doing one cell. *Device / host*: the GPU and its memory / the CPU and its memory; they are separate address spaces. *MultiFab*: AMReX's container for one field (say `U`) across all the boxes of a mesh; each box's piece is a plain contiguous array. *K1*: kernels written as restricted C++ loop bodies called through AMReX's `ParallelFor`. *K2*: kernels kept in Fortran with OpenMP `target` directives. *Shim*: the small wrapper at the call site that pulls arrays out of FDS or AMReX data and passes them to a kernel.

## 1. Why GPUs choke on what FDS does inside loops

A GPU thread can only use memory that has been placed on the device, and the compiler must know where each array is. A flat array passed as an argument is one address and a set of bounds, so it is trivial to place. FDS reaches most data another way:

- **Derived types with allocatable components** (`WALL(IW)%...`, `M%U`, `SF%...`). Copying the outer variable to the device does not copy what its allocatable components point to; there is no automatic deep copy in OpenMP offload, CUDA Fortran, or C++. Each component has to be mapped by hand, and a table of them (like `WALL`) is a structure whose parts live in many places.
- **Pointers re-aimed at run time** (`UU => U`, `RHOP => RHOS`). The compiler cannot tell what a pointer refers to, whether two pointers overlap, or which array must be on the device. It serialises or refuses.
- **Module variables** used inside loops. They are hidden inputs; each must be mapped separately and are easily missed.
- **I/O, `ALLOCATE`, `STOP`** inside a loop. None of them exist on the device.

Survey numbers: of 829 loops that could run on the device, **542 have at least one "layout" blocker of this kind** (layout meaning pointers, derived-type components or allocatable data) and only 178 have none; about 294 routines (46,000 lines) would need rewrites; 275 loops re-point a pointer with `=>`, and 339 read allocatable or pointer components of derived types.

**Example A: a loop that is nearly clean (`mass.f90:442-455`, predictor step for density).**
Before:
```fortran
DO K=1,KBAR; DO J=1,JBAR; DO I=1,IBAR
   IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE        ! two-level indirection into a derived type
   RHS = - DEL_RHO_D_DEL_Z__0(I,J,K,N) + (FX(I,J,K,N)*UU(I,J,K)*R(I) - ...)*RDX(I)*RRN(I) + ...
   ZZS(I,J,K,N) = RHO(I,J,K)*ZZ(I,J,K,N) - DT*RHS
```
The arithmetic is already GPU-ready. The arrays are module pointers set by `POINT_TO_MESH` (and `UU => U`), and `SOLID` sits behind `CELL_INDEX`. After: the same body, but `UU, FX, ..., ZZS, RDX, R` arrive as dummy arguments with `INTENT`, and the solid flag arrives as an integer array `SOLID_MASK(I,J,K)` built once per mesh. Only the header and one line change.

**Example B: a wall loop (`mass.f90:424-436`).**
Before: `WC=>WALL(IW)`, `BC=>BOUNDARY_COORD(WC%BC_INDEX)`, then `UU(BC%IIG-1,BC%JJG,BC%KKG) = UVW_SAVE(IW)` inside a `SELECT CASE (BC%IOR)`. Two pointer hops through derived types, and a write to an index computed from data.
After: at setup, flatten the wall records used here into integer arrays `WALL_TYPE(IW)`, `WALL_IIG(IW)`, `WALL_JJG(IW)`, `WALL_KKG(IW)`, `WALL_IOR(IW)`; the loop reads those and writes `UU(...)`. Each wall cell writes its own face value, so threads do not collide; that is what the GPU needs.

**Example C: the hard case (`wall.f90:3641-3683`, the 3-D heat-conduction sweep over wall cells).** Each iteration follows `M%WALL(IW)`, `M%BOUNDARY_ONE_D(WC%OD_INDEX)`, `THR_D%NODE(I)%ALTERNATE_WALL_INDEX(II)` and `MESHES(NM2)%...`: derived types containing allocatable arrays of derived types, reached through other meshes, plus an `!$OMP CRITICAL` section. No amount of argument-passing tidies this; it needs a structure-of-arrays rewrite (one flat array per field, plus an integer table saying where each wall cell's data starts) or it stays on the host. This is the shape of most of the 294-routine rewrite list, and why ADR-001 allows host-only routines.

## 2. Code-hygiene changes, ranked by payoff against merge friction with upstream FireX

| # | Change | Payoff for GPU | Friction with upstream merges |
|---|---|---|---|
| 1 | Hot loops in their own routines that take **flat arrays as arguments** with `INTENT(IN/OUT)` and explicit bounds; the caller does `POINT_TO_MESH` and passes `M%U` | Highest: it is the boundary the shim needs, in either style | Low to moderate: a mechanical edit per routine, and loop bodies stay as they are (arithmetic-only hunks still merge). |
| 2 | **No `ALLOCATE`, I/O, `STOP` or `RANDOM_NUMBER` inside loops**; pre-size scratch arrays outside | Medium: 12 device-eligible loops allocate and 9 do I/O, so the count is small, but each one blocks its whole nest | Low: local edits. |
| 3 | **Small `PURE` leaf functions** for repeated pieces (species lookups, table interpolation); `ELEMENTAL` where scalar-in scalar-out | High: callees are inlined on the device only if the compiler can see them are side-effect free; 223 of the 224 reachable routines are not `PURE` | Low: adds a keyword; can break if a routine writes a module variable. |
| 4 | **Explicit `INTENT`** everywhere in hot routines | Medium: also documents which arrays a kernel reads and writes, which the merge checker uses | Low. |
| 5 | **Replace a derived-type indirection in a hot loop** by a plain array (`CELL(CELL_INDEX(I,J,K))%SOLID` to `SOLID_MASK(I,J,K)`) | Medium to high: 339 loops read such components | Moderate: keeps two representations in sync, or replaces one in every routine that reads it. |
| 6 | **Avoid module state in hot loops**: pass scalars as arguments; group read-only constants and lookup tables in one place | Medium: 4 loops write module variables (many more read them, which the survey does not count) | Moderate: touches call signatures upstream also changes. |
| 7 | **Split wall loops from cell loops** and give each wall loop flat index arrays (Example B), keep the physics call on the host or on a separate kernel | Medium: wall loops interleave with cell loops on the same arrays, so they force syncs | Higher: reshapes routines upstream edits often (`wall.f90` is among the most-changed candidate files, after `pres.f90` and `ccib.f90`). |
| 8 | Avoid pointer re-aiming (`UU => US`) inside a routine: pass the chosen array instead | Medium: 275 loops | Low to moderate. |

## 3. What could be proposed upstream to FireX

Items 1, 2, 3, 4 and 8 are ordinary code-quality changes with no GPU vocabulary in them: smaller routines, explicit `INTENT`, no allocation in loops, `PURE` leaves. They can be offered as small independent pull requests that change no results (verified bitwise by the Verification suite). Each one moves loops from the 542 with a layout blocker towards the 178 that have none, and lowers the cost of **both** K1 (the C++ translator has flat arguments to read) and K2 (the same Fortran body can take the `target` directive). Items 5 to 7 are more invasive and are best done in our branch first, one routine at a time, with the merge protocol in ADR-001 (kept names, rename script, mapping table) as the safety net. Upstream churn to keep in mind: 64% of changed loop hunks are body-only arithmetic, and about 5 to 6.6 commits a month touch a candidate routine.

Caution: this is a proposal to consider, not a plan. Upstream maintainers have their own priorities, and each request has a review cost.

## 4. How this bears on K1 versus K2 (neutral)

- **Both need the same data work.** Flat arguments, flattened wall records and no allocation in loops are required whichever style is chosen (ADR-001 "GPU kernel data layout"). Hygiene changes lower the cost of either.
- **Where K1 gains more.** K1 rewrites each loop body in C++. If upstream keeps bodies as pure arithmetic on flat arguments, each re-port is a translation of a small loop; if not, it is a translation of the loop and its data access. Hygiene shrinks the K1 re-port list (about 8 to 12.5 routine re-ports a month).
- **Where K2 gains more.** K2 keeps the Fortran body, so a hygienic routine can carry the directive in place and merge as ordinary Fortran. `PURE` leaves and explicit `INTENT` matter for the compiler's device-side analysis; K2's known weak points (nvfortran quirks with `private` scalars and callee-heavy loops) are lessened by simple, call-free bodies.
- **What does not change.** The toolchain question (nvfortran and OpenMP offload versus C++ and CUDA), the readability split among reviewers, and the need for bitwise tests on the GPU. The hygiene work makes the choice cheaper to reverse per kernel, which supports keeping K1 available as a per-kernel fallback if K2 is chosen.
