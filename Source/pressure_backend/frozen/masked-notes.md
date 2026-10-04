# Masked branch of the common layer (A-68, FR-037 masked branch): design and status

Scope of this version: **single level, MLMG backend** (built and tested, see the measured results below). The FFT backend cannot take masks (the selector says so), HYPRE and the composite
paths still return `NotBuilt` for masks (status table at the end).

## Why it is needed

docs/pressure/07 section 14.3: masked cases (obstructions treated as solid, gap cells, known-value cells) are not legal for any product
backend, so the masked product-level comparison cannot be run. The masked solve itself exists only in study harnesses (assembled HYPRE and
MLMG on the FDS-style operator). The missing piece is the masked branch in the common layer, not a new solver.

## Cell classes (`PressureProblem::cell_class`, one `int` per cell, same BoxArray/DistributionMapping as `ba`)

| code | name | meaning |
|---|---|---|
| 0 | Gas | unknown of the system; has an equation |
| 1 | Solid | excluded: OBST or gap cell. No unknown, no equation. Every face of the cell has coefficient 0 (no flux, exact Neumann wall). |
| 2 | Known | known-value cell behind a face-Dirichlet condition (FDS OPEN value on a vent face). Not an unknown. The cell holds the value `g` on the face it shares with the gas cell (`PressureProblem::known_value`, null = 0). The gas neighbour gets the flux `2 (g - phi) / dx^2` through that face (the linear face-Dirichlet flux with ghost `2g - phi`, which is what the domain Dirichlet faces and FDS give). |

Any other code is `InvalidInput`. A pinned class is not offered: pins belong to the common layer and are not used by MLMG. An all-zero
`cell_class` (or null) is the unmasked problem and takes the unchanged unmasked code (bitwise).

## Discrete operator (defined in the common layer, `apply_operator`)

For a gas cell `i`, summed over the six neighbours `n` (periodic wrap included):

    (L phi)_i = sum_n c_n (phi_n - phi_i) / dx_d^2           c_n = 1 gas neighbour, 0 solid neighbour, 2 known neighbour with phi_n := g_n,
                                                              domain faces as in the unmasked operator (Neumann: 0, Dirichlet: -2 phi_i / dx^2)

`(L phi)` is zero (and phi is not read) on Solid and Known cells. A one-cell y direction drops its terms as before (D-057).

## Components (flood fill over faces with coefficient != 0)

Computed by the common layer (`label_components`): gas cells connected through gas-gas faces, periodic wrap included, are one component.
Algorithm: label = global x-fastest index of the cell, repeated sweeps of "label := min over gas neighbours" inside each box, then a ghost
exchange, until no label changes (global reduction); the surviving minimum is the lowest global index of the component, which is also its
**pin cell** (recorded, not applied by MLMG). Ids are the ranks of these minima in increasing order, so they do not depend on the
decomposition. A component is **singular** unless one of its gas cells touches a domain Dirichlet face (after `effective_bc`) or a Known
cell. Non-gas cells get label -1 (excluded from every sum). Driver-supplied `component_id` stays `NotBuilt` (only the computed partition
exists; the driver's zone ids can be added by mapping them onto the computed components when a case needs it).

## Mean removal, gauge, residual (D-067, per component, gas cells only)

Unchanged code, driven by the labels: `remove_mean` removes the volume-weighted mean per singular component over its gas cells (exact sum);
`apply_gauge` shifts each singular component by the `rho V` weighted mean of `H - KRES` over its gas cells (exact sum); non-singular
components are not shifted. Solid cells are 0 and Known cells hold `g` in the returned `phi`. The true-residual check zeroes the residual
and the right-hand side on non-gas cells, so norms, `||b||`, the floor of the size-scaled limit and the cell count are over gas cells.

## MLMG realisation: no overset mask, identity rows for excluded cells

MLMG's overset-mask route (docs/pressure/04 section 2(a2)) needs a pin per sealed component, a positive `a` on every masked cell and a
HYPRE bottom solver. Used here instead: a symmetric positive semi-definite `MLABecLaplacian` with **face coefficients and a cell
coefficient array**, no mask, so MLMG sees an ordinary (well-posed on every excluded cell) problem:

    A phi = s,    A = alpha - div(beta grad)        (a = b = 1)
    s = -rhs + sum_known 2 g / dx_d^2    on gas cells,       s = 0 on Solid and Known cells
    beta_face = 1 between two gas cells and on a domain face of a gas cell, 0 on every face that touches a Solid or Known cell
    alpha = sum over the Known faces of the cell of 2 / dx_d^2   on gas cells
    alpha = sum_d 1 / dx_d^2                                     on Solid and Known cells (identity row: phi = 0, decoupled, A positive)

Multiplying `L phi = rhs` by -1 and moving the known values to the right-hand side gives exactly this system on the gas cells, so the
solution on gas cells is the solution of the masked `L`. Excluded cells are decoupled by beta = 0 and solve to 0 (the common layer then
writes `g` into Known cells). The system is positive semi-definite; a sealed component keeps its null space (constants), which is why the
mean removal before and the gauge after are essential, exactly as for the unmasked MLMG. MLMG does not flag it singular (alpha > 0
elsewhere), so it does not subtract anything itself; its bottom solver sees a consistent system.

Known cost: MLMG coarsens `alpha` and `beta` by arithmetic averaging, so a coarse cell that mixes excluded and gas cells leaks on the
coarse levels. That changes the convergence rate, never the answer (the fine residual is exact); the iteration count is reported and tested.

## Dense reference test

`mode=masked` of the harness builds random solid cells (about 10 %), a solid slab that splits the box into two components, one component
with a Dirichlet domain face and one sealed, optionally Known cells, writes the right-hand side, the class, the known values and the
returned `H` to raw files; `tests/pb_test.py masked` builds the dense matrix with numpy from its own flood fill, solves each component
(least squares with the constant removed for sealed ones) and compares per component: solution to 1e-9 (up to the gauge), component count,
singular flags, per-component mean of the mean-removed right-hand side and of `H`.

## Measured results (ctest `pb_masked`, 12x10x8 cells, 1/2/3 ranks)

Six cases (sealed slab, Dirichlet face with a slab, Known cells, periodic wrap with the rho-weighted gauge, random solids with one periodic
direction, Dirichlet on x) all pass against the independent dense reference:

| case | components | iterations | true residual / limit | solution error vs dense |
|---|---|---|---|---|
| sealed, solid slab | 2 sealed | 31 | 4.3e-13 / 1e-12 | 3.4e-12 |
| slab, Dirichlet z-high | 1 sealed + 1 open | 38 | 6.9e-13 / 1e-12 | 5.2e-12 |
| slab + Known cells | 5 | 39 | 5.2e-13 / 1.1e-12 | 3.9e-12 |
| periodic x,y, Known cells, rho gauge | 3 | 13 | 1.2e-13 / 1.5e-12 | 5.2e-13 |
| random solids, periodic y, rho gauge | 1 | 41 | 5.6e-13 / 1.9e-12 | 5.8e-12 |
| Dirichlet x-low, slab | 2 | 199 | 7.7e-13 / 1e-12 | 1.3e-12 |

Labels, singular flags, cell counts and the removed mean agree exactly with the reference; the returned `H` and the labels are identical on
1, 2 and 3 ranks (to 1e-9 for `H`); the gauge sum is below 1e-12 (volume-weighted `H` for the default, `rho (H - KRES)` for the weighted gauge);
`H` is exactly 0 on Solid cells and exactly `g` on Known cells; a loose-tolerance negative control warns; masked + FFT is `NotBuilt`.

Multigrid quality (honest): with 10 % random solids the iteration count is 13 to 40 against about 10 without a mask; a Dirichlet face with a slab can need 100 to 200
iterations (the last row hits 199 of the default 200). Averaged coefficients leak on the
coarse levels, so the masked solve uses 8 pre/post smoothing sweeps and a BiCGStab+CG bottom solver (plain BiCGStab stalled on two sealed
components). The relative tolerance is scaled by `|rhs|_inf / |s|_inf` so that MLMG's own norm and the common-layer residual mean the same
when Known terms dominate `s`. Isolated sealed cells (no gas neighbour, no Dirichlet face) get an identity row and `H = 0` before the gauge.
A 50 % random-solid case (27 components) converges in 57 iterations; 30 % solids with Known cells in 111 (both checked by hand, not in ctest).
The slowest tested configuration is the Dirichlet x face with a slab (199 of 200 iterations): a masked problem with a Dirichlet face and a thin
slab may need a larger `max_iter`; a Krylov outer loop around the V-cycle or the HYPRE masked path would be the remedy.

## Status (what is built and what is not)

| Item | Status |
|---|---|
| Single level, MLMG, Gas/Solid/Known classes, computed components, per-component mean removal and gauge, residual check | built and tested (`pb_masked`) |
| FFT with masks | `NotBuilt` by design (no mask in `FFT::Poisson`); `Auto` picks MLMG for a masked single level |
| HYPRE assembled backend with masks | NotBuilt (next: rows of Solid and Known cells replaced by identity, Known neighbours moved to the right-hand side, the same system as above) |
| Composite (two or more levels) with masks | NotBuilt (components would have to be found on the composite graph) |
| Driver-supplied `component_id` | NotBuilt (computed components only) |
| Pinned class, covered cells on the single-level API | NotBuilt |
| Component labelling cost | one flood fill per solve (no cache yet); fine for tests, a workspace cache keyed on the mask is the obvious next step for runs with many solves per step |
