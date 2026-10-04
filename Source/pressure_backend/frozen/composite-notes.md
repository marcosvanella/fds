# Composite (multi-level) MLMG pressure solve: notes

Scope: first composite solve, CPU, `amrex::MLMG` over `MLPoisson`. Interface in `PressureIface.H`
(`PressureLevel`, `PressureProblem::levels`, `PressureWorkspace`, `solve_pressure`, `face_gradient_composite`).

## What is solved

Poisson problem `L phi = b` on the composite grid: every level l has valid cells; a cell is *uncovered* if no finer level
refines it. Unknowns are phi on all cells; the equations are enforced on uncovered cells only. The operator is the
flux-form 7-point Laplacian with reflux: at a coarse/fine (C/F) face the coarse flux is the average of the fine fluxes.
Covered coarse cells carry no equation; after the solve they hold the average-down of the fine solution.

Uncovered masks are derived inside the backend from the next finer level's BoxArray coarsened by that level's `ref_ratio`;
the caller passes no mask. RHS values in covered cells are ignored (zeroed in the work copy).

## Boundary conditions

Domain faces (coarse domain): Dirichlet (homogeneous), Neumann, periodic, per problem; a per-direction mix of Neumann and
periodic is accepted (all closed faces). Mixing open and closed faces is NotBuilt. C/F interface: AMReX MLMG, `maxorder 2`,
fine ghost cells interpolated from the coarse level (coarse-fine BC from the coarse level, `setLevelBC(l, nullptr)` on every
level). The C/F interpolation is first order at the interface, so the fine-side gradient at C/F faces is less accurate than
elsewhere (up to about 12 % for ratio 4 in the corner-patch tests); the solution itself converges at second order.

## Compatibility condition and mean removal (singular problems)

For all-Neumann or all-periodic faces the composite problem is solvable only if
`sum over uncovered cells of v_l * b = 0`, with `v_l` the cell volume of the cell's level. This is the discrete
integral of the divergence over the composite grid and it is the condition that makes reflux consistent. A plain arithmetic
(ScaledArithmetic, equal weights) mean would be wrong here because the number of cells and their volumes differ between
levels, so the composite path uses the volume-weighted ("Volume") mean. The removed constant is computed with the exact
sum (`ExactSum.H`, `exact_sum_multi`: one fixed-point scale for the whole hierarchy, so the result does not depend on the
decomposition or summation order), subtracted from uncovered cells only, and repeated once when needed; the
idempotence floor is 2^-52 of the largest |b| over uncovered cells. `removed_mean` and `removed_rel` (volume-weighted rms)
are reported as before. No pin is used with MLMG.

## Gauge

phi is defined up to a constant. The returned phi has zero exact volume-weighted mean over uncovered cells (same exact sum,
same scale); the shift is applied to all levels before average-down. `gauge_weight/offset` are NotBuilt for composite.

## True residual

After the final phi, `MLMG::compResidual` (reflux + average-down, covered cells masked out) gives `b - L phi` on uncovered
cells, independent of the MLMG convergence number. `residual_rel2` is the volume-weighted L2 over uncovered cells and
`residual_relmax` the max, both relative to the mean-removed RHS. The common layer's check compares against eps_H = 1e-8.
Measured: 1e-14 to 4e-13 in all tests. The tests additionally check the divergence of the returned face gradients against
the RHS independently of MLMG.

## One-cell direction (the 2-D case)

If the level-0 domain has one cell in a direction d (e.g. ns2d_16: y), the solve runs on an extruded copy: `ext_n` = 4
isotropic cells times the product of ratios per level in that direction, periodic, same boxes and DistributionMapping; data
is replicated in and the plane copied out. The gradient in d is exactly 0. Reason: AMReX `setHiddenDirection` converged
slowly with more than one level (73 to 200 iterations on a case that needs 9 in 3-D) and diverged for a fine level over the
whole domain; the extruded version converges in 7 to 10 iterations and matches a real 3-D 4-cell problem to 2.4e-13.
Cost: 4x the cells per level in that direction. The ratio in that direction must be 1; Dirichlet in a one-cell direction
and more than one one-cell direction are NotBuilt.

## Gradient helper

`face_gradient_composite` builds a local extruded workspace, then `prepareForSolve`/`prepareForFluxes` and
`MLPoisson::getFluxes(FaceCenter)` per level; gradient = -flux. Then `average_down_faces` from finest to coarsest so coarse
C/F faces hold the average of the fine fluxes. A coarse face on a periodic domain boundary exists twice (lo index and
hi+1 index); if one copy lies under a finer level the other copy takes the same fine average (ParallelAdd with periodicity).

## Workspace

`PressureWorkspace::rebuild(p)` stores layout copies (BoxArray, DistributionMapping, Geometry, ratios, flags) and no
pointers into any MultiFab. After a regrid call `rebuild` with the new hierarchy; a stale workspace is also rebuilt
automatically (`PressureResult::workspace_rebuilt`). With a null workspace each call builds and discards its setup. A
failed rebuild leaves the workspace empty.

## Nesting and layout rules (InvalidInput)

Level 0 covers its domain exactly; each level's domain equals the coarser domain refined by `ref_ratio`; boxes are
coarsenable by `ref_ratio`; properly nested with one coarse-cell buffer except at non-periodic domain faces (periodic
images considered); rhs/phi layouts match; Geometry periodicity matches the BC.

## NotBuilt (clear message, phi untouched)

Composite on the FFT backend; cylindrical or non-Cartesian geometry; non-uniform `cell_width`; coefficients;
`component_id`, `cell_class`, `uncovered` masks; `gauge_weight/offset`; mixed open/closed faces; Dirichlet in a one-cell
direction; anisotropic ratios (except 1 in the one-cell direction); ratios other than 2 and 4; non-homogeneous Dirichlet.
HYPRE is not used (default bottom solver); if a HYPRE bottom solver is added later, set `hypre.adjust_singular_matrix=1`.

## Tests (ctest, `pb_comp_*`)

Two-level centered patch with ratio 2 and 4 for Neumann/periodic/Dirichlet (order, residual); fine level over the whole
domain equals the single-level fine solve (1e-14 to 1.5e-13) with coarse = average-down; decomposition independence (1, 2,
3 ranks, box sizes 32/16/8, SFC/round-robin/knapsack) and bitwise repeatability; three levels; gradient checks (coarse =
fine average at C/F faces, divergence consistency, negative control); corner/wrapping/two-patch layouts; mixed
Neumann/periodic per direction; selector/NotBuilt messages; workspace reuse and rebuild; ns2d_16 hierarchy (16x1x16
level 0, level 1 refines coarse cells 4..11 in x and z with ratio (2,1,2)).

Convergence (L2 error vs the manufactured solution, coarse n to 2n): ratio 2 Neumann 1.045e-3 to 2.608e-4 (order 2.00),
periodic 4.18e-3 to 1.038e-3 (2.01), Dirichlet 9.65e-4 to 2.417e-4 (2.00). Ratio 4 (16 to 32): 2.01, 2.03, 2.00. Difference
against a uniform-fine solve on the fine region converges at 1.9 to 2.1.
