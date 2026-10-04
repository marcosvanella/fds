# M2 notes: solve triggers (FR-039), FFT plan cache, mixed faces, boundary data, D-057 mapping

This file is appended to as the M2 items land. Sections are independent.

## 1. Solve-trigger points (FR-039)

`PressureOptions::trigger` is a bit set of `SolveTrigger` values describing why this solve is being made:

| flag | value | meaning |
|---|---|---|
| `SolveRoutine` | 0 | ordinary time-step solve |
| `SolveFirst` | 1 | first solve with a workspace (detected automatically from `PressureWorkspace::num_solves() == 0`) |
| `SolveFirstAfterRegrid` | 2 | first solve after a regrid (set automatically by `PressureWorkspace::rebuild()` and when a solve sees a layout that differs from the workspace; the caller may also set it) |
| `SolveDebug` | 4 | debug run, check everything |

`PressureOptions::full_checks_on` (default `SolveFirst | SolveFirstAfterRegrid | SolveDebug`) selects which flags switch the full checks on.
`trigger` defaults to `SolveDebug`, so a caller that never touches the option gets the full checks on every solve, as before.
The effective flags (caller flags plus the ones the workspace detected) are returned in `PressureResult::triggers`;
`PressureResult::full_checks` says which path ran.

Full checks: structural layout validation, pin scan, rms / removed-mean compatibility diagnostic and its warning,
true residual (`residual_checked`).
Cheap path: selector, rhs/phi layout checks, mean removal, backend solve, gauge, convergence status.
`NotBuilt`, `InvalidInput` and `NotConverged` are reported on both paths.
The solution, the removed mean and the gauge constants are bitwise identical on both paths (tests below).

Tests: `pb_trigger_single_NNNNNN`, `pb_trigger_single_PPPPPP` (FFT and MLMG, np 1 and 2), `pb_comp_trigger` (hierarchy, np 1 and 2).

## 2. FFT plan caching

`PressureWorkspace` owns an `FFTBackend` that caches the `amrex::FFT::Poisson` object, keyed by BoxArray, DistributionMapping,
Geometry and the six boundary types. A repeated solve with the same key reuses the plan (`BackendStatus::plan_reused`,
`PressureWorkspace::fft_plan_reuses()`); any key difference, `rebuild()`, `discard()` or a composite rebuild builds a new plan
(`fft_plan_builds()`). A workspace result is bitwise identical to a no-workspace result (test `pb_fftcache_*`).
Single-level MLMG has no cache (the MLMG object is cheap to build; not requested).

Timing, 20 solves each, 1 core machine with oversubscribed ranks (so noisy; ratios are what matter), Neumann FFT, Routine trigger:

| grid | ranks | fresh plan per solve | cached plan per solve | ratio |
|---|---|---|---|---|
| 32^3 | 1 | 5.5 ms | 4.4 ms | 1.26 |
| 32^3 | 4 | 5.2 ms | 4.9 ms | 1.07 |
| 64^3 | 1 | 45.5 ms | 36.0 ms | 1.26 |
| 64^3 | 4 | 14.9 ms | 12.2 ms | 1.22 |
| 128^3 | 1 | 340 ms | 324 ms | 1.05 |
| 128^3 | 4 | 104 ms | 97 ms | 1.07 |

The saving is the plan construction, a few percent to about a quarter of a solve; it is not large because the transforms
dominate. `plan_build_ms` in the harness line is the time of an explicit `rebuild()`.

## 3. Mixed open/closed faces

Every per-direction combination of low/high face types in {periodic-periodic, Neumann-Neumann, Dirichlet-Dirichlet,
Neumann(lo)-Dirichlet(hi), Dirichlet(lo)-Neumann(hi)} is built on a single unmasked level, for the FFT backend
(`FFT::Poisson` with even/odd/periodic pairs; mixed even/odd pairs use the half-cell offset) and for MLMG, and on the
composite (MLMG) path. The old restriction "all six faces closed or all six open" is removed.

Verification (`pb_mixed_faces_sweep`): all 125 combinations of {PP, NN, DD, ND, DN} per direction on a 6x5x4 grid with a
non-smooth, incompatible right-hand side, FFT and MLMG, 1 and 2 ranks, against a dense numpy operator built independently
from the ghost rules (Neumann ghost = +phi, Dirichlet ghost = -phi, periodic wrap). Singular systems: mean-removed rhs, least
squares solution, zero-mean gauge. Worst relative difference: FFT 5e-15, MLMG 1e-13. The singular flag of every combination
is also checked (a problem is singular exactly when no face is Dirichlet).

Composite (`pb_comp_faces_*`): two-level hierarchy, manufactured solutions cos((k+1/2) pi x) for N(lo)-D(hi) and
sin((k+1/2) pi x) for D(lo)-N(hi), error against the exact solution at n = 32 and 64: second order (1.89 to 2.00), true residual
at round-off, vs the uniform fine solve converging at second order. For mixed directions the harness does not require the
fine-level error to be below the coarse uniform error at one resolution (the solution need not vanish at the patch edges); the order is
what is checked.

### One-cell directions (D-057), `effective_bc()`

A one-cell y direction is the FDS TWO_D case: its term is dropped by FFT::Poisson (factor 0 for any length-1 direction, whatever the BC)
and by MLPoisson for Neumann or periodic faces. Dirichlet faces of a one-cell y therefore act as Neumann:
`effective_bc()` (PressureIface.H) is applied in the FFT and MLMG backends, the composite solve, the reference operator
(`apply_operator`, so the true-residual check agrees) and the singular-component test (an all-open 2-D box with a one-cell y is
non-singular because of its x and z faces; a box with Dirichlet only in y is singular). A Dirichlet face in a one-cell x or z
direction keeps the term -2 phi/dx^2 and is NotBuilt (single level and composite), with a message. See d057-mapping.md.

## 4. Inhomogeneous boundary data: fold into the right-hand side

`fold_boundary_data(p, BoundaryData, rhs)` (PressureIface.H, implemented in CommonLayer.cpp) adds the known ghost-cell terms
of non-zero wall data to the right-hand side of the homogeneous problem the backends solve. Call it once on the rhs
before `solve_pressure`; the solver's mean removal then acts on the folded right-hand side (for a pure-Neumann problem with
compatible data the folded rhs sums to zero).

Sign convention (FDS pres.f90: header comment lines 57-58, boundary application 452-454, ULMAT F_H terms 1620-1679):

- Neumann data g is dH/dx_d along the increasing coordinate at both the low and the high face.
- Dirichlet data is H at the wall.
- Ghost values: low Neumann phi_g = phi_1 - h g; high Neumann phi_g = phi_n + h g; Dirichlet phi_g = 2 H_b - phi_1.
- Fold (h = cell width of the direction, boundary-adjacent layer of cells only, one term per face, nothing else multiplies it):

| face | rhs change |
|---|---|
| low Neumann | `+ g / h` |
| high Neumann | `- g / h` |
| Dirichlet, low or high | `- 2 H_b / h^2` |

Data per face is a cell-centred MultiFab (any BoxArray / DistributionMapping, read in the face-adjacent layer; uncovered
layer cells read as 0) or a uniform constant. Periodic faces take no data (InvalidInput, rhs untouched). A one-cell y
drops its term, so its data are ignored; a Dirichlet face in a one-cell x or z is NotBuilt. Uniform cell width only.

Tests: `pb_bcdata_exact` builds the rhs with an independent explicit ghost-value operator from slab and constant data (five face
mixes, FFT and MLMG, 1 and 2 ranks); the folded homogeneous solve reproduces the field to 1e-15 relative, and the same solve without the
fold is wrong by 0.3-0.5 (negative control). `pb_bcdata_mms` uses u = cos(1.3x+0.2) cosh(0.7y) sin(0.9z+0.4) with exact wall values and
exact increasing-direction derivatives: second order for DD,DD,DD; NN,NN,NN; ND,DN,NN; DN,ND,DD (error 1.4e-5 at n=64 for DD).
