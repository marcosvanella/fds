# Masked branch: design decisions D1 to D8 (A-71; review docs/pressure/08 section 3.2)

Companion to `masked-notes.md` (operator, classes, MLMG and assembled-HYPRE realisation). For each request of the review: what is decided, why, which
code it binds, which test proves it, and its status. Single level + MLMG is built (D1, D2, D5, D7, D8). Single-level assembled HYPRE masked is built (D6).
Composite masked is built on assembled HYPRE (D3, D4, D6; Q1). MLMG composite masked stays NotBuilt. Closed decisions Q1–Q6 are recorded below; Q7 remains for the Architect.

## D1. Component labels are cached in the workspace (built: a2a4906647, copy-free A-71)

* **Decision.** Labels, singular flags, cell counts and pins (`ComponentMap`) are kept in `PressureWorkspace` and reused while the mask
  is unchanged. Key: content hash of `cell_class` and the layout (domain, face BCs, BoxArray, DistributionMapping). A regrid or any
  class change is a miss and relabels. Only the per-solve fields (`removed_mean`, `removed_rel`, `gauge_shift`) are reset on a hit. The
  labels are used in place (no per-solve copy of the label `iMultiFab`; A-71 removed one). The 2^31 cell limit is a selector refusal
  with a message (`PressureIface.cpp`), the assertion in `label_components` stays behind it. Layout: labels live on the problem's BoxArray/DM, so every sum traverses them in the same contiguous
  order as the field; each relaxation round exchanges labels only (one `FillBoundary`) plus one integer reduction.
* **Why.** Rounds grow with box-to-box hops along the longest path. Measured (4 ranks, box under load ~40):

  | geometry | cells | boxes | components | rounds | label time | study solve time (doc 06 section 8) |
  |---|---|---|---|---|---|---|
  | stairwell, `cutk` cut | 18.76 M | 648 (mgs 32) | 3 (2 sealed) | 18 | 1.76 s | Mf 2.78 s, H 1.09 s (4 CPU) |
  | stairwell, `cutk` cut | 18.76 M | 81 (mgs 64) | 3 (2 sealed) | 10 | 1.05 s | |
  | stairwell, uncut | 18.76 M | 81 (mgs 64) | 1 | 14 | 1.30 s | |
  | hallways | 0.52 M | 128 / 16 / 2 (mgs 16/32/64) | 1 | 19 / 9 / 4 | 0.71 / 0.26 / 0.14 s | |

  Labelling is 40 % to 160 % of a solve on the stairwell, so it must not run per solve. A cache hit costs one hash of `cell_class`
  (measured 4e-6 s to 1e-2 s on the test sizes).
* **Binds.** `PressureIface.cpp` (`solve`, cache lookup and reset), `CommonLayer.cpp` `label_components`, `PressureWorkspace`.
* **Tests.** `pb_masked` cache runs (np 1 and 3): second solve is a hit with a bitwise identical answer; a changed known value is a hit
  (D8); toggling one cell Gas to Solid is a miss and relabels. `pb_masked_hallways` (label `long`): every solve after the first is a hit.
  Harness `mode=masked_geo case=stairwell solve=0` reproduces the table.

## D2. Gauge key is (zone, component); a component spanning two zones is an error (built, A-71)

* **Decision.** `PressureProblem::component_id` now carries the driver's zone ids (one `int` per cell, read on gas cells). Components are
  still computed from the mask; each component takes its zone (`ComponentInfo::zone`). A component whose gas cells carry two zone ids is
  `InvalidInput` ("component c spans zones a and b ..."). Since a component lies in one zone, the key (zone, component) is the component
  itself: mean removal, gauge and pin stay per component, and the zone is a label used by the driver (D_PBAR_DT per zone). Two components
  with the same zone are legal (FDS zone split by a solid; each sealed part has its own null space).
* **Why.** Doc 06 section 5.1: one gauge per zone and connected component. One connected gas region is one Poisson problem and cannot have
  two background pressures.
* **Binds.** `CommonLayer.cpp` `assign_zones` (per component min/max zone, one integer min/max reduction), called from
  `PressureIface.cpp` `solve` after labelling; `PressureIface.H` `component_id`, `ComponentInfo::zone`. Single level only; the composite
  path keeps `NotBuilt` for `component_id` (`CompositeSolve.cpp`) until D4.
* **Tests.** `pb_masked`: `zones=1` (two zones separated by the one-cell solid slab) gives components with zones 0 and 1;
  `zones=2` (zone boundary cutting a component) gives `InvalidInput` with the message.

## D5. `Solid` is a gap cell; an OBST cell is `Solid` only under E-2 (documented)

* **Decision.** Gap cells (bounding box minus the meshes) are `Solid` in every configuration; this is the case the 40 FR-006 inputs need.
  An obstruction becomes `Solid` only when the case runs under E-2 (doc 01 REC-E1, doc 04 A4). Then for those faces the
  `NO_FLUX`/`WALL_WORK1` iteration is off (doc 01 section D: the wall is exact Neumann), `D = 0` inside them is irrelevant, and the average
  down over a masked fine patch uses gas cells only (doc 01 section E). Under E-1 (default, FDS fidelity) obstructions are Gas cells of
  the operator and FDS's forcing and pressure iteration act on them as before. The choice E-1/E-2 is the spec's, not the backend's.
* **Binds.** The driver's construction of `cell_class` (Role 1); the backend treats every `Solid` alike (`apply_operator`,
  `MLMGBackend::solve_masked`). The average-down rule binds D4.
* **Tests.** None in the backend (no backend code depends on the reason a cell is `Solid`). The E-2 run is a driver test.

## D7. Krylov outer loop: GMRES with one MLMG V-cycle as preconditioner (built, A-71)

* **Decision.** The masked single-level solve is `GMRESMLMG` (AMReX) preconditioned by one V-cycle with 8 pre and 8 post sweeps and the
  smoother as bottom solver (a fixed linear preconditioner, which GMRES requires). `max_iter` stays 200. The tolerance is mapped once
  from the common-layer norm to GMRES's own L2 norm of `s` (`rel = tol_rel |rhs|_2 / |s|_2` when Known terms dominate `s`); the true
  residual check is unchanged and decides. Replaces the V-cycle solve with BiCGStab+CG bottom and the max-norm retarget (fewer lines).
* **Why.** The V-cycle alone needed 13 to 199 iterations on the 12x10x8 cases and did not reach its own norm in 1000 V-cycles on hallways
  (the true residual was inside the limit; its own norm was dominated by Known terms). Masked coarse levels leak (averaged coefficients),
  which GMRES absorbs. Measured with GMRES:

  | case | old iterations | GMRES iterations | true residual / limit |
  |---|---|---|---|
  | sealed slab / Dirichlet slab / Known / periodic / random / Dirichlet x | 31 / 38 / 39 / 13 / 41 / 199 | 9 / 10 / 9 / 7 / 9 / 12 | 1.7e-13 to 7.1e-13 / 1e-12 to 1.1e-12 |
  | hallways 128x64x64 (0.52 M cells, 73,728 gas), mgs 16 / 32 / 64, 2 ranks | NotConverged at 1000 (mgs 32, 64) | 11 / 11 / 11 | 3.7e-12 / 6.8e-11 |

  Iterations are identical on 1, 2 and 3 ranks and over `max_grid_size` 16 to 64 (FR-005). Memory: GMRES keeps one Krylov vector of
  the level size per iteration (restart 30), so about 12 extra fields at 11 iterations. Times: hallways 6 to 17 s per solve on a box with
  load average about 40; not comparable with the study's quiet-machine rows (open hallways Mf 20 iterations, 0.28 s on 1 rank). The study
  used a HYPRE bottom solver and a dropped-box layout; see Q4 and Q5.
* **Binds.** `MLMGBackend.cpp` `solve_masked`.
* **Tests.** `pb_masked` fails when any run uses more than 95 % of the 200-iteration cap (18 runs); `pb_masked_hallways` (label `long`,
  hallways at mgs 32 and 64 on 2 ranks) checks status, true residual, the 5 % margin and the D1 cache. Stairwell solve: harness
  `mode=masked_geo case=stairwell` exists but was not run (Q4).

## D8. `known_value` is rebuilt every pressure iteration from the OPEN branch of L1209 (documented; backend side tested)

* **Decision.** For a Known cell behind an OPEN vent face, the driver writes `known_value = BXS/BXF/... ` of that vent wall from
  `pres_poisson_boundary_arrays` (L1209 OPEN branch: `KRES` of the gas cell when the normal velocity points out, `H0` otherwise,
  pres.f90:176-210) **before every solve**, including every pressure iteration. The value is the face value (the operator puts `g` on the
  face, ghost `2g - phi`). Under E-2 the other Dirichlet branches of L1209 (solid pseudo-Dirichlet plus `WALL_WORK1`, interpolated mesh
  interfaces) do not arise. The backend never caches `known_value` (only the class is in the cache key).
* **Binds.** Driver pressure iteration (Role 1) for the copy; backend: `PressureIface.cpp` reads `known_value` each solve; the OPEN branch
  is `FdsPressureLoops.H` `pres_poisson_boundary_arrays`.
* **Tests.** Backend: `pb_masked` cache run changes the known values between two solves with the same mask: cache hit and the solution
  changes (dev = 1.19). L1209 sign dependence: `pb_fdsloops_*` compare the OPEN branch with FDS bitwise, inflow and outflow. The
  end-to-end chain (flip the normal velocity, `g` switches between `KRES` and `H0`, `H` changes) is a driver integration test (Role 1).
  UGLMAT face placement stays `[VERIFY]` (Q3). Partial OPEN vents on a bounding-box face: Q6.

## D3. Composite masked: Solid and Known cells never enter a C/F interpolation stencil (built)

* **Rule.** (a) A coarse cell is usable in a C/F interpolation stencil only if it is uncovered `Gas`. The sliding-parabola stencil
  (`HypreBackend.cpp` `bdry_stencil`, the AMReX `interpbndrydata_o3` form) already drops to a one-sided or constant form when a tangential
  neighbour is not available; the change is one more condition in the availability predicate `avail()` (class == Gas). The mixed
  (cross-derivative) term is used only when all four diagonal neighbours are usable. (b) The C/F face coefficient is 0 when the fine-side
  cell is `Solid`, and when the coarse base cell of the fine ghost is `Solid`. (c) A coarse base cell that is `Known` is not allowed in the
  first version: a Known cell at a C/F face is `InvalidInput` (the regrid keeps vents inside one level; Q2). (d) Reflux: the coarse face
  flux is the average of the fine face fluxes, and a fine face with a `Solid` side contributes 0. (e) Average-down over a masked fine
  patch uses gas fine cells only (D5); a covered coarse cell with no gas fine cell is 0.
* **Why.** Solid and Known cells hold dummy `phi = 0` (identity rows); any stencil reading them corrupts the gas flux next to a corner.
* **Binds.** `HypreBackend.cpp` assembled composite (`avail_gas`, `bdry_stencil_gas`, C/F/reflux; same-level only on-level). MLMG composite with masks stays `NotBuilt`:
  MLMG's C/F interpolation cannot exclude cells per stencil (Q1).
* **Test.** Two-level case whose patch boundary touches a solid corner and a solid edge: (1) dense numpy reference of the composite
  masked matrix built from the same rule, solution to 1e-9; (2) changing the dummy value in Solid cells (0 versus 1e30 in the input
  `phi`) changes nothing bitwise; (3) true composite residual inside the limit.

## D4. Composite components (built)

* **Rule.** A component is a set of uncovered gas cells of all levels connected by same-level faces with coefficient /= 0 and by C/F
  faces between an uncovered coarse gas cell and an uncovered fine gas cell (fine side Gas and coarse side Gas, D3 b). Covered cells are
  not part of any component. One mean removal, one gauge and one pin per component, with the composite exact sum (level volumes as
  weights, the same `exact_sum_multi` as the composite residual). Pin: the lowest global index on the coarsest level that holds a cell
  of the component (matches the unmasked composite pin). Zone (D2): read per level on uncovered gas cells; a component spanning zones is
  `InvalidInput`.
* **Algorithm.** The same min-label relaxation, with labels unique across levels (level offset added to the global index) and two
  single-direction transfers per round besides the same-level sweep: fine to coarse, each covered coarse cell next to an uncovered gas
  cell takes the minimum label of the gas fine cells touching that shared face (a face-min, not a cell average); coarse to fine, fine
  ghost cells outside the patch take the label of their coarse cell when it is uncovered gas. Converges when no label changes on any
  level (one reduction per round). Cached as D1, keyed on the classes of all levels and the level layouts.
* **Why.** Per-level labelling would split a corridor whose middle is refined into two or three components with different gauges.
* **Binds.** `CommonLayer.cpp` (`label_components` on the composite), `CompositeSolve.cpp` (mean removal, gauge, residual, which today
  refuse masks and `component_id` with `NotBuilt`).
* **Tests.** (1) Corridor whose middle third is covered by a fine patch: one composite component, one gauge,
  `sum rho V (KRES - H) = 0` over the composite set to 1e-12. (2) Sealed pocket inside a fine patch: its own component, own pin, own
  mean. (3) Labels identical on 1 and 3 ranks and over `max_grid_size`.

## D6. Assembled HYPRE masked: one pin row per sealed component (built)

* **Rule.** (a) The assembled system has unknowns on gas cells only: Solid and Known cells get no row (row numbering is a per-box prefix
  count over gas cells, one pass), Known neighbours go to the right-hand side as `2 g / dx^2` (the same system as MLMG). This replaces the identity rows of `masked-notes.md` (fewer rows and
  FLOPs) and keeps the matrix conditioning of the gas problem. (b) Each sealed component gets one pin row at its
  `ComponentInfo::pin` (identity, value 0 after mean removal; the gauge sets the constant after). (c) `BackendStatus` reports a list of
  pin cells (level, cell); `residual_sums` and `evaluate_residual` (single level and composite) exclude every pin row.
* **Why.** A sealed component of an assembled system is singular; BoomerAMG needs it fixed. A single pin cell cannot handle two sealed
  components.
* **Binds.** `HypreBackend.cpp` (`MaskedHypre`: gas-only assembly, pin rows), `PressureIface.H` `BackendStatus::pins`, `CommonLayer.cpp`
  `residual_sums` / `evaluate_residual`, `CompositeSolve.cpp` residual.
* **Tests.** `pb_masked` HYPRE rows on the six small cases vs MLMG; sealed2: two pins reported, residual excludes both; `max_pins=1` warns (singular).
* **Status.** Built for single level and composite (multi-pin constructor on `HypreSystem`).
* **D-091.** Masked single-level `Auto` → HYPRE (Architect accept; A-56). Explicit MLMG override kept.

## Closed decisions (Q1–Q6) and remaining open

* **Q1 (closed).** Composite masked → assembled HYPRE only; MLMG composite masked stays `NotBuilt` (clear selector error).
* **Q2 (closed).** Known (OPEN) on a C/F face is `InvalidInput` in v1 (vents stay inside one level).
* **Q3 (closed, doc).** UGLMAT OPEN face-linear β=2; backend Known handling already matches. Lead reading recorded as accept.
* **Q4.** Stairwell product solve still not run on the shared box (memory); not blocking this path after D-091 HYPRE default.
* **Q5 (open, Role 1).** Covering↔sparse BoxArray API; do not block backend. Sparse ba/dm OK later.
* **Q6 (closed).** Driver pads Known/Solid halo for partial OPEN on a domain face; backend assumes that.
* **Q7 (open, Architect).** Zone merge / `component_id` field name when an obstruction opens (CONNECTED_ZONES); Role 1 confirms.

## Partial-fine C/F + mask (fixed)

* **Bug.** Fine C/F halo ghosts are not filled by `FillBoundary` (no fine box there). Trusting `cell_class`/`row_id` in that halo treated the neighbour as same-level and dropped the C/F coupling, so partial-fine masked HYPRE did not converge.
* **Fix.** Same-level ghost path only when the wrapped neighbour lies on this level's `BoxArray`; otherwise C/F / reflux / physical (D3). Solid-adjacent faces stay coeff 0; Covered coarse gas still refluxes; Known-on-C/F remains InvalidInput (Q2).
* **Evidence (`pb_comp_masked`, n=16, ratio=2, Dirichlet):**
  * full-fine + one fine Solid: Ok, 19 iters, residual_check 4.0e-12 / limit 1e-8
  * partial-fine corridor + Solid at uncovered coarse (1,1,1): Ok, 18 iters, residual_check 8.2e-12 / limit 1e-8
* **Product.** Composite masked remains assembled-HYPRE only (Q1); MLMG composite+mask stays NotBuilt. Pressure Lead sign-off still required before treating partial-fine as the product default path.
