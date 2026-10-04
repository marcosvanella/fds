# HYPRE assembled-matrix backend (ADR-002 backend (H))

Files: `HypreBackend.H/.cpp`, selection in `PressureIface.cpp` and `CompositeSolve.cpp`, harness modes in `harness/`, tests `pb_hypre_*`.

## What it is

`BackendKind::HYPRE`, requested through `PressureOptions::backend` (the option is on `PressureOptions`, not on `PressureProblem`).
Auto never selects it. It assembles the discrete operator of the MLMG backend into a HYPRE IJ matrix from host buffers
(`HYPRE_MEMORY_HOST`, no managed memory) and solves with a Krylov method preconditioned by one BoomerAMG V-cycle. Single level and
composite hierarchies use one code path (a single level is a hierarchy of one). The common layer does what it does for every
backend: validation, components, mean removal before the solve (D-067: volume-weighted default, `MeanKind::ScaledArithmetic` FDS parity
switch), the gauge after it, the independent true-residual check (MLMG's operator, not HYPRE's matrix).

## Operator (identical to MLMG's; verified entry by entry)

Unknowns are the uncovered cells of all levels, numbered rank-major, then level, then box, then x fastest (no communication is
needed to assemble). Covered cells keep a decoupled diagonal row. The matrix is the volume-scaled operator A = -V L (positive
definite, FDS ULMAT scaling), the right-hand side F = -V b.

- Same level, both cells unknowns: (phi_n - phi_c)/h^2 per face. A direction with one cell has no term (D-057).
- Physical faces: Neumann ghost = phi_c; homogeneous Dirichlet ghost = -phi_c; periodic wraps. Inhomogeneous boundary data is folded
  into the right-hand side by the common layer (`fold_boundary_data`), not here.
- Fine cell next to a C/F face (ratio r normal to the face): MLMG's ghost is the linear extrapolation through the boundary datum b
  at the centre of the coarse cell (r/2 fine cells outside the face) and the first interior cell:
  ghost = a b + (1 - a) phi_c, a = 2/(r+1) (2/3 for r = 2, 2/5 for r = 4). b is AMReX's order-3 sliding-parabola interpolation
  (`interpbndrydata_*_o3`) of the coarser level: centre, first and second difference along the tangential neighbours that are in
  the domain and not covered by the fine level, and the mixed term when all four diagonal neighbours are available. Each ghost is a
  fixed linear combination of coarse unknowns, so the matrix has no hidden state.
- Coarse cell next to a covered cell: the coarse face flux is replaced by the average of the fine face fluxes a (phi_f - b_f)/h_f
  over the r^(active tangential directions) fine faces (MLMG reflux, fine-flux consistency).
- Singular problems (no Dirichlet face in a direction that carries a term): the lowest-index uncovered cell of the coarsest level that
  has one is pinned, as FDS ULMAT/GLMAT does: its row is a diagonal-scaled identity with right-hand side 0 and its column is dropped
  from all other rows (this keeps a single level symmetric). The common layer's gauge then replaces whatever constant the pin gave.

Oracle (`mode=hypre_op`): HYPRE's matrix applied to a random field equals MLMG's composite operator (`compResidual` of b = 0) to
round-off: relative difference 2.5e-16 to 3.7e-16 for single level (periodic, Neumann, Dirichlet), two levels with ratio 2 and 4
(periodic, Neumann, Dirichlet), layouts middle/corner/two patches/split, and three levels with ratio 4 (3.2e-16).

## Krylov method and AMG settings

`HypreOptions::krylov = Auto`: PCG on one level, GMRES(30) on a hierarchy. The composite matrix is not symmetric (the C/F
interpolation is one-sided); PCG still converges on the cases tested (22 iterations against 22 for GMRES on a 32-cell two-level case,
BiCGSTAB 13 iterations of two preconditioner applications), but is not guaranteed to, so Auto uses GMRES there. PCG, GMRES and
BiCGSTAB can be requested explicitly. Defaults follow FDS (coarsen type 8 PMIS, relax type 18 L1-Jacobi, one sweep, HYPRE default
strength and interpolation). A tuned set (`coarsen=10 relax=6 strong=0.5 agg=1`) gives fewer iterations and a cheaper set-up
(table below) but is not the default: relax 6 is hybrid Gauss-Seidel, whose iterates depend on the rank count. The tolerance is
`PressureOptions::tol_rel` (the solver is run to 0.1 tol_rel because its recurrence residual is slightly optimistic) and the cap is
`max_iter`; the iteration count and HYPRE's own final relative residual are in `BackendStatus`.

## Workspace

The matrix, the AMG hierarchy and the solver are built on the first solve and cached in `PressureWorkspace` (single level and
composite), keyed by BoxArray, DistributionMapping, geometry, effective BC and `HypreOptions`; a second solve reuses them
(`plan_reused`, `hypre_built()`). `rebuild()`, a different layout or different HYPRE options drop and rebuild them (reported as
`workspace_rebuilt`). Nothing in the set-up points into the caller's MultiFabs (only index and BoxArray data), so the old hierarchy can be
destroyed after `rebuild` (tested in `pb_hypre_comp_workspace` and `pb_hypre_cache_*`).

## Not built

Same rules as MLMG, same messages: masked cells, covered cells on the single-level API, non-uniform widths, cylindrical, variable
coefficients, anisotropic ratios, ratios other than 2 and 4, a Dirichlet face in a one-cell x or z direction, an invalid hierarchy.
These are checked before any backend runs; `phi` is untouched. HYPRE is not cell-count limited beyond HYPRE's integer type (the
numbering is `HYPRE_BigInt`).

## Agreement (eps_H = 1e-8; every number is a relative L2 difference of the final H, tol_rel 1e-12, 4 ranks unless noted)

Single level, synthetic right-hand side (`pb_hypre_single_*`):

| case | HYPRE vs FFT | HYPRE vs MLMG | MLMG vs FFT | PCG iterations (MLMG) |
|---|---|---|---|---|
| Neumann 34x18x32 | 1.2e-13 | 1.2e-13 | 1.3e-14 | 22 (9-12) |
| periodic 34x18x32 | 3.4e-14 | 3.5e-14 | 9.7e-15 | 22 |
| Dirichlet 34x18x32 | 7.0e-15 | 4.9e-14 | 4.7e-14 | 20 |
| Neumann 64^3 | 9.8e-15 | 4.2e-14 | 4.9e-14 | 27 (12) |
| periodic 64^3 | 9.2e-15 | 4.4e-14 | 4.0e-14 | 26 (11) |
| Dirichlet 64^3 | 9.7e-15 | 4.0e-15 | 8.6e-15 | 23 (12) |
| one-cell y 16x1x16, Neumann | 8.3e-15 | 1.3e-12 | 1.3e-12 | 18 |
| one-cell y 16x1x16, periodic | 4.6e-15 | 7.2e-14 | 7.1e-14 | 19 |

Mixed faces (`pb_hypre_mixed_faces_sweep`): all 125 per-direction combinations of PP, NN, DD, ND, DN on 6x5x4, against a dense reference
operator, worst relative error 9.3e-15 (1 rank) and 1.1e-14 (2 ranks, 3-cell boxes). Inhomogeneous wall data (`pb_hypre_bcdata_exact`,
five face sets): max error 3.7e-14 to 1.7e-13.

CI gate (`pb_hypre_ci_vs_mlmg`): the same frozen solves (generated right-hand side and FFT reference H, 34x18x32) through MLMG and HYPRE
fail the test if they differ by more than eps_H: Neumann 1.3e-13, periodic 3.6e-14, Dirichlet 4.8e-14 (max abs 1.5e-14, 2.9e-15, 3.0e-16). The
test also runs a negative control with tol_rel = 1e-3, where the pair differs by 9.8e-6 and the comparison must report FAIL.

Composite, 16 coarse cells (`pb_hypre_comp_vs_mlmg`, GMRES; true composite residual by MLMG's operator, bound 1e-8):

| hierarchy | HYPRE vs MLMG | HYPRE true residual | iterations HYPRE (MLMG) |
|---|---|---|---|
| Neumann, 2 levels, ratio 2 | 1.5e-13 | 6.9e-14 | 22 (9) |
| Neumann, 2 levels, ratio 4 | 9.2e-14 | 2.9e-13 | 24 (11) |
| periodic, ratio 2, patch at the corner | 1.1e-13 | 9.6e-14 | 21 (9) |
| periodic, ratio 4 | 4.1e-13 | 3.3e-13 | 24 (9) |
| Dirichlet, ratio 2 / ratio 4 (two patches) | 2.5e-14 / 2.3e-14 | 8.0e-14 / 2.7e-13 | 21 (9) / 24 (10) |
| Neumann, 3 levels, ratio 2 | 1.9e-13 | 2.1e-13 | 22 (9) |
| periodic, 3 levels, ratio 4 | 4.5e-13 | 6.9e-13 | 28 (10) |
| Neumann, fine level over the whole domain | 1.0e-14 | 4.6e-14 | 23 (16) |
| one-cell y Neumann / periodic (plane2d) | 5.6e-14 / 3.4e-14 | 5.2e-14 / 3.9e-14 | 20 (9) / 20 (9) |
| mixed faces ND,DN,NN / PP,NN,DD | 2.0e-13 / 6.8e-13 | 9.3e-14 / 5.1e-14 | 21 (9) / 21 (9) |
| ns2d_16 two-level (1, 2, 4 ranks) | 9.7e-14 / 9.1e-14 / 9.3e-14 | 8.6e-14 / 4.1e-14 / 6.1e-14 | 19 / 20 / 19 (9) |

The harness also checks that the covered
coarse cells equal the average-down of the fine solution and that every HYPRE composite solve is bitwise repeatable.

## Decomposition and repeatability

- Single level, 64^3 (`pb_hypre_decomp_*`): ranks 1, 2, 4 with boxes of 32 and 16 against 1 rank/32: relative L2 3.4e-15 to 5.7e-15 (Neumann
  and periodic); not bitwise (the summation order of HYPRE changes), within eps_H by seven orders. The pin and the removed mean are bitwise equal.
- Composite (`pb_hypre_comp_vs_mlmg`): Neumann ratio 2 and periodic ratio 4 on 1, 2, 4 ranks against 1 rank: 4.8e-14, 5.4e-14 and 2.9e-13, 2.1e-13.
- Run-to-run (`pb_hypre_repeat_*`, three runs on 4 ranks at fixed layout) and the composite harness (fresh workspace, reused set-up, three
  solves): bitwise identical.

## Timing (single-core box shared with other jobs, so noisy; seconds for one solve inside `solve_pressure`, includes the common layer's
checks, which run MLMG's residual)

| problem | FFT | MLMG | HYPRE first solve (set-up + solve) | HYPRE tuned | HYPRE with cached set-up |
|---|---|---|---|---|---|
| 32^3 Neumann | 0.009 | 0.027 | 0.29 (0.22 + 0.06; 22 it) | 0.28 (19 it) | about 0.05 |
| 32^3 periodic | 0.009 | 0.019 | 0.37 (23 it) | 0.24 (19 it) | |
| 64^3 Neumann | 0.08 | 0.43 | 4.0 (2.2 + 1.7; 26 it) | 1.7 (0.47 + 1.2; 22 it) | about 1.7 |
| 64^3 periodic | 0.05 | 0.32 | 4.4 (26 it) | 2.8 (21 it) | |
| composite 2 levels, ratio 2, 32 coarse cells, 2 ranks | | 0.045 (10 it) | 0.25 (0.18 + 0.06; 22 it) | | 0.076 |
| composite 3 levels, ratio 4, 16 coarse cells (299k rows), 3 ranks | | 1.1 (10 it) | 2.5 (28 it) | | about 1.2 |

HYPRE is about 5-10 times slower than MLMG on one core for the same accuracy (MLMG's geometric V-cycle is cheaper on a structured grid
than BoomerAMG with the PMIS coarsening that FDS uses), and about 50 times slower than FFT where FFT applies. The set-up is roughly half of
a first solve and is amortised by the workspace. Where the box has many cores the picture may differ; this is not measured.

## Limitations

1. Residual floor on large singular problems. The identity pin concentrates the round-off residual of all other rows in the pin row, so the
   independent true residual of a singular problem grows like sqrt(N) times round-off and does not fall with tighter tolerance: 6e-13 at
   48^3 Neumann for tol_rel 1e-12, 1e-13, 1e-14 and 1e-15 alike (the largest entry is the pin cell), 1.3e-13 at 64^3 and 6.8e-12 at 96^3. The
   solution is correct (rel 1e-13 against FFT/MLMG after the gauge); with the default `residual_tol = 1e-12` the common layer warns at about
   10^6 unknowns. Raising `residual_tol` or comparing the residual away from the pin is the remedy; removing the leftover mean once more
   and restarting the Krylov solver were tried and change nothing. Non-singular problems (any Dirichlet face) are not affected.
2. The composite matrix is not symmetric; GMRES is the default there, PCG works on the cases tested but is not guaranteed.
3. `face_gradient_composite` is MLMG-based whatever backend produced H.
4. The backend option is on `PressureOptions`, not `PressureProblem`.
5. Speed: slower than MLMG on this box (above). No GPU, no device memory (host buffers only, by design).
6. Defaults follow FDS (PMIS, L1-Jacobi); the tuned set above is about twice as fast but is a user option, not the default.
