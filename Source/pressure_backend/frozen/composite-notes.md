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

## Compatibility condition and mean removal (singular problems, D-067)

For all-Neumann or all-periodic faces the composite problem is solvable only if
`sum over uncovered cells of v_l * b = 0`, with `v_l` the cell volume of the cell's level. This is the discrete
integral of the divergence over the composite grid and it is the condition that makes reflux consistent. A plain arithmetic
(equal weights) mean of `b` would be wrong here because the number of cells and their volumes differ between
levels.

`PressureProblem::mean_kind` selects the zero mode removed from `b` (default `Volume`, D-067):

* `MeanKind::Volume` (default): the composite volume-weighted mean `sum(v*b)/sum(v)` over the uncovered cells is subtracted
  from every uncovered cell of every level (D-032).
* `MeanKind::ScaledArithmetic` (FDS parity switch): the arithmetic mean of `F = v*b` over the uncovered cells is removed from
  `F`, i.e. `b_k -= mean(F)/v_l(k)` with a different constant on each level. This is what FDS ULMAT/UGLMAT do on a
  single mesh. On a hierarchy it attributes an incompatible part of `b` as a constant per cell instead of per unit volume;
  it is correct only when the part is truly per cell. It is for parity studies, not for production AMR runs.

Both make `sum(v*b)` zero, and they agree to round-off whenever `b` is already compatible. `removed_mean` is the constant
subtracted from `b` (Volume) or the mean of `F` (ScaledArithmetic); `removed_rel` is its size relative to the rms of `b`
(Volume) or of `F` (ScaledArithmetic). The sums are exact (`ExactSum.H`, `exact_sum_multi`: one fixed-point scale for the
whole hierarchy), so the removed constant is bitwise independent of the decomposition. Application is repeated once when
needed; the idempotence floor is 2^-52 of the largest |b| (or |F|) over uncovered cells. No pin is used with MLMG.
AMReX `makeSolvable` is not relied on.

## Gauge (D-067)

phi is defined up to a constant. For a singular component the returned phi satisfies
`sum over uncovered cells of V_l * rho * (KRES - phi) = 0`: the constant `sum(V_l*rho*(phi - KRES)) / sum(V_l*rho)` is
removed from all levels (exact sums, numerator and denominator each with one fixed-point scale over the whole hierarchy),
then the covered coarse cells take the average-down of the fine solution. `rho` is `PressureLevel::gauge_weight` (FDS RHOP:
RHO in the predictor, RHOS in the corrector) and `KRES` is `PressureLevel::gauge_offset`; both are per level, on the level's
own `BoxArray`/`DistributionMapping`, one component, read on uncovered cells only (covered values are ignored). They are
optional: null on every level means rho 1 and KRES 0, which is the plain exact volume-weighted mean zero. Giving a field on
some levels only, or with a different layout, is `InvalidInput`. The single-level `PressureProblem::gauge_weight` /
`gauge_offset` together with `levels` are `NotBuilt` (the message points to the level fields). A component that is not
singular (some Dirichlet face) has a unique solution and is not shifted, whatever the gauge fields. The same gauge is applied
by the single-level path for both the FFT and the MLMG backends (rho default 1, KRES default 0 gives the plain volume mean,
the earlier behaviour). FDS adds `20 eps` to the denominator; we do not.

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
`component_id`, `cell_class`, `uncovered` masks; the single-level `gauge_weight/offset` fields together with `levels`; mixed open/closed faces; Dirichlet in a one-cell
direction; anisotropic ratios (except 1 in the one-cell direction); ratios other than 2 and 4; non-homogeneous Dirichlet.
HYPRE is not used (default bottom solver); if a HYPRE bottom solver is added later, set `hypre.adjust_singular_matrix=1`.

## Tests (ctest, `pb_comp_*`)

D-067: `pb_comp_gauge_<bc>_r2`, `pb_comp_gauge_neumann_r4`, `pb_comp_gauge_plane2d_<bc>` (two-level hierarchy with per-level rho and
KRES and an offset RHS: removed means, gauge constant, the `sum(V rho (KRES-H)) = 0` condition and the shifted field against independent
numpy formulas on level dumps; Dirichlet: nothing shifted; divergence of the face gradients against `b - mean(v*b)/v` for the
parity switch, with the Volume zero mode as negative control), and `pb_comp_gauge_decomp_<bc>` (1 to 4 ranks, box sizes 8 to 32, three
mappings, two patch layouts: removed means bitwise identical, exact gauge sums of an analytic field bitwise identical, gauge
constants of the solve equal to 1e-16 absolute).


Two-level centered patch with ratio 2 and 4 for Neumann/periodic/Dirichlet (order, residual); fine level over the whole
domain equals the single-level fine solve (1e-14 to 1.5e-13) with coarse = average-down; decomposition independence (1, 2,
3 ranks, box sizes 32/16/8, SFC/round-robin/knapsack) and bitwise repeatability; three levels; gradient checks (coarse =
fine average at C/F faces, divergence consistency, negative control); corner/wrapping/two-patch layouts; mixed
Neumann/periodic per direction; selector/NotBuilt messages; workspace reuse and rebuild; ns2d_16 hierarchy (16x1x16
level 0, level 1 refines coarse cells 4..11 in x and z with ratio (2,1,2)).

Convergence (L2 error vs the manufactured solution, coarse n to 2n): ratio 2 Neumann 1.045e-3 to 2.608e-4 (order 2.00),
periodic 4.18e-3 to 1.038e-3 (2.01), Dirichlet 9.65e-4 to 2.417e-4 (2.00). Ratio 4 (16 to 32): 2.01, 2.03, 2.00. Difference
against a uniform-fine solve on the fine region converges at 1.9 to 2.1.
