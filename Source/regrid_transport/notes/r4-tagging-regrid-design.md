# R4 design note: tagging and dynamic regrid

Status: draft for the Chief Architect. No code yet. Covers milestone R4 of `plan.md`, on top of `docs/role3-numeric-scheme-note.md` (with addendum) and `docs/role3-fds-mesh-interface-note.md`. Spec FR-011..013, FR-015, FR-024, FR-025, IR-003, IR-008; ADR-001 (K2 default), ADR-002. Parameter names below are **working names**, final only with the User Guide draft (IR-003). Items not read in source or measured are marked **(assumption)** or **(to verify)**.

**Words used.** *Level*: a grid resolution; level 1 is finer than level 0 by the ratio 2 or 4. *Patch*: one rectangular box of cells on a level. *Tag*: a flag on a cell meaning "refine here". *Regrid*: rebuilding the finer levels from the current tags. *Tag buffer*: tags are grown by a few cells so a moving feature stays inside its patch between regrids. *Proper nesting*: a fine patch lies inside the coarse level with a margin of at least `N_PROPER` coarse cells, so the coarse-fine interface never touches a coarser interface. *Blocking factor*: every patch edge falls on a multiple of this many cells. *Clustering*: AMReX covers the tagged cells with few, mostly-full boxes (efficiency `GRID_EFF`). *Covered* coarse cell: a coarse cell under a fine patch. *Average-down*: replace covered coarse values with the volume average of their fine cells. *Face-centred*: stored on cell faces (velocity U, V, W). *Refinable region*: where IR-008 allows refinement. *Registry*: the driver's per-level field storage (Role 1 `LevelRegistry`).

## 1. Tagging

**Principle.** Each level is tagged from its own data at the end of a full step. Covered coarse cells already hold the fine average, so a criterion gives the same answer whether or not a patch exists below it, and a patch does not "tag itself away" (a difference of neighbouring cells shrinks on finer cells, so evaluating it on fine data would flicker). Tags are the OR of all active criteria, then (1) forced tags on the static boxes of the finer `&MESH` entries so they are never dropped (IR-002), (2) tags outside the refinable region discarded and counted (existing `TagClipper`, IR-008), (3) no tags at `MAX_LEVEL`.

| Criterion (FR-011) | Working names in `&AMR` | Test per cell |
|---|---|---|
| Temperature | `TAG_TMP_DIFF` (K), `TAG_TMP_ABOVE` (K above ambient) | max over 6 neighbours of abs(T - T_nb) > threshold; or T - T_ambient > threshold |
| Density | `TAG_RHO_DIFF` (relative) | max abs(rho - rho_nb) / rho > threshold |
| Species or mixture fraction | `TAG_SPEC_ID`, `TAG_SPEC_DIFF` | same difference test on Y of the named species |
| Heat release rate per unit volume | `TAG_HRRPUV` (kW/m3) | HRRPUV > threshold (**assumption**: the cell array is kept at end of step) |
| Vorticity | `TAG_VORT` (1/s) | magnitude of curl u from centred differences of face velocities > threshold |
| User boxes | existing `&AMR_REGION` (+ `LEVEL`) | forced inside the box (already R1) |
| Distance to OBST/VENT | `TAG_DIST`, `TAG_DIST_ID` | needs OBST on levels: **deferred to Phase 5** |

Differences are *undivided* (no division by dx), so a feature is refined until its jump is spread over cells and refinement then stops by itself. A threshold may be one value or one per level. Phase 3 implements temperature, density, species, HRRPUV and user boxes first; vorticity is cheap and follows, since Phase 3 velocity is prescribed (question 1).

**Buffer, nesting, interval.** `N_ERROR_BUF` (default 1, working) grows every tag; with a regrid every `REGRID_INTERVAL` steps (default 4) the buffer must cover the distance a feature travels in that many steps (advection CFL below about 0.5 per step gives at most 2 cells at interval 4, so buffer 2 is the safe value for moving features; the test checks it). `N_PROPER`, `BLOCKING_FACTOR`, `MAX_GRID_SIZE` are enforced by AMReX when it builds the grids and re-checked by the FR-013 debug assertions. `REGRID_INTERVAL = 0` keeps the static hierarchy. Regrid runs between full steps only, before the next dt is computed, so the minimum-dt rule sees the new hierarchy.

**Hysteresis.** A cell that is already refined stays tagged while its criterion is above `TAG_KEEP` times the threshold (working name, default 0.8); a new tag needs the full threshold. Together with the buffer this removes patches that appear and disappear at alternate regrids. A level whose new grids equal the old ones is left untouched (AMReX remakes a level only if its grids changed, or the level below changed; a remake with equal grids is then a plain copy, **to verify** as bitwise unchanged).

## 2. Kernels and parallel behaviour

**Kernels (ADR-001).** The only kernels I write are tagging kernels. Default is K2: one Fortran routine per criterion family on flat arrays (cell values with one ghost layer, a tag array, thresholds), OpenMP `target` loops, called through the driver shim, no AMReX types inside. K1 (C++ `ParallelFor`) is justified per kernel only where K2 needs data it cannot take as flat arrays: (a) counting tags per criterion for the log (a reduction, easy with AMReX reduction helpers); (b) the AMReX tag array (one byte per cell, owned by AMReX) if the shim cannot pass it flat. Each K1 use gets a recorded reason in `OWNERS.md`. Tags are cell-local; the only communication is one ghost-layer exchange of the inputs before tagging.

**Ranks and boxes.** Each rank tags its own patches; AMReX gathers the tag set, clusters it, cuts to `MAX_GRID_SIZE`, and every rank gets the same BoxArray (list of boxes). Level 0 is fixed by `&MESH`, not by the rank count, so only the *mapping of boxes to ranks* changes with ranks. The tag set depends only on cell data, which the driver already keeps independent of box split and ranks (exact sums, FR-005); hence the hierarchy at 1, 2 and 4 ranks should be identical (FR-015) and is tested. Risk: in Phase 4 the velocity comes from a pressure solve that need not be bitwise rank-independent, and a vorticity threshold could then flip one cell. FR-015 only demands repeatability at a fixed rank count; the 1-versus-4 test uses the Phase 3 prescribed velocity.

## 3. Data transfer at a regrid

AMReX regrids top-down from level 1: for each level it calls `MakeNewLevelFromCoarse` (level is new), `RemakeLevel` (grids changed) or `ClearLevel` (level gone). The coarser level is already final when the finer one is processed. `RegridAmrCore` forwards each call to the registry (`make_level`, `remake_level`, `clear_level`) and then fills the new level.

**Order for a new or remade level:** (1) average-down all levels (debug: check covered coarse = fine average) and record composite mass and species with the exact sum; (2) registry allocates the new fields, old fields stay readable; (3) fill coarse-level ghosts (same-level exchange, coarse-fine hook, physical boundary); (4) copy old fine valid data where the new patch overlaps the old one (bitwise); (5) interpolate the rest from the level below; (6) refill ghosts; (7) recompute derived fields; (8) run `DIVERGENCE_PART_1` on changed levels and print max abs(div u - D); (9) registry rebuilds side data (per-level lists for ghost exchange, wall/boundary records, exact-sum masks) and Role 3 rebuilds flux-override face lists; (10) old fields freed; then dt is computed.

| Class | New fine cells or faces | Dropped or covered | Properties |
|---|---|---|---|
| rho*Z_i (all tracked species), rho | Limited conservative linear interpolation of each rho*Z_i (AMReX conservative interpolators with slope limiter); rho = sum | Fine patch dropped: nothing to hand back, coarse already holds the average-down value | Children average equals parent to round-off; limiter gives no new extremes, so rho*Z_i >= 0, Y in [0,1], sum Y = 1; clips counted (FR-025) |
| Derived cell fields (T, MU, KRES, D, cp ...) | Recomputed by the kernels from rho, Z, background pressure | same | not interpolated, so no inconsistent T |
| D and lagged diffusion cell array | One divergence pass on changed levels | same | not conserved quantities |
| Face velocity | Interface faces (fine face on a coarse face): coarse value. Interior faces: `FaceDivFree` interpolation (**to verify**: ratio 4, 2-D, single-cell direction) | Coarse face = area average of fine faces (average-down of faces) | Each fine cell keeps its parent's divergence; leftover is D's variation inside the parent, absorbed by DDDT at the next solve; momentum not conserved (as FDS) |
| Stage arrays (US, RHOS, ZZS, HS, FV) | Not transferred, rebuilt in the next predictor | | |
| Pressure H | Not transferred; optional linear interpolation as solver first guess | | re-solved each stage |
| Ghost cells | Refilled, never transferred; both coarse-fine layers piecewise-constant (ruling b) | | |
| Wall/boundary data | None in Phase 3 (no walls); Phase 5: per-area values copied from the parent face, totals split by area | area-weighted mean | |

The post-regrid projection stays off by default; the printed max abs(div u - D) after each regrid is the measurement that decides whether it is turned on. In Phase 3 the velocity is a prescribed analytic field re-evaluated on the new faces, so only the cell transfer is exercised. At t = 0 the initial hierarchy is built by repeating tag, make level, fill until `MAX_LEVEL` (question 2).

**Pressure backend (Role 2).** After the grids change the composite operator and solver setup are rebuilt once (measured 0.016 to 0.11 s for the MAC-projection backend at about 1M unknowns, `docs/pressure/05`). I need one entry point "hierarchy changed, rebuild" on the backend interface, called after step (9); the backend must not keep pointers to old MultiFabs.

**Registry (Role 1).** Requests: (a) bracketing calls `begin_regrid()` / `end_regrid()` so old fields are freed only after all copies; (b) a `Fields` that can be built on an arbitrary BoxArray and DistributionMapping at run time (the Architect's choice of mesh objects for level above 0 is pending); (c) `run_divergence_part1(level)` and `max_divergence_error(level)` as already declared in `RegridInterface.H`.

## 4. Tests (`tests/`, 1 and 4 ranks, one thread)

1. Per criterion: tagged cells equal an independent naive evaluation, cell by cell, including box edges and rank boundaries (FR-011).
2. Hierarchy: fixed synthetic tag pattern gives the same boxes at 1, 2 and 4 ranks and on repeats; nesting, blocking factor, max grid size, forced static boxes, discarded-tag count (FR-013, FR-015).
3. Conservation: random positive multi-species field, several regrids that add, move and drop patches; composite mass and each species change <= 1e-12 relative per regrid, expected near 1e-15; overlapping old fine data bitwise unchanged; uniform field stays bitwise uniform; a linear field is reproduced exactly; Y in [0,1], sum 1.
4. Velocity: divergence-free coarse field gives fine interior divergence at round-off; a field with divergence keeps the parent value; interface faces equal the coarse value.
5. Moving blob: prescribed translation plus divergence-free periodic field, 2-D and 3-D, ratios 2 and 4, regrid every few steps; conservation at every regrid and step; every cell above the threshold is on the finest allowed level after each regrid; result against the uniform-fine run reported.
6. Negative controls (each must fail the test): non-conservative interpolation (point sampling) breaks test 3; buffer 0 with fast blob loses the feature in test 5; hysteresis off raises the patch change count; one rank's tags altered breaks test 2; skipping average-down before a drop trips the debug check.

## 5. Effort and open questions

Estimate (unmeasured): about 9 to 11 working days. Roughly 4 days need nothing from Role 1 (tagging kernels, transfer operators and tests 1 to 4 on bare MultiFabs); about 5 to 7 days follow Role 1's level-above-0 kernel loop (blob test, registry integration, timing).

Questions for the Architect:
1. Phase 3 criterion set: temperature, density, species, HRRPUV and user boxes first, vorticity next, OBST distance in Phase 5?
2. Initial hierarchy at t = 0: interpolate from coarse (same path as regrid) or evaluate the input's initial conditions directly on new levels (better for features smaller than a coarse cell)?
3. Hysteresis and difference-based tags: accept `TAG_KEEP` and undivided differences, or prefer a scaled gradient?
4. Approve the `begin_regrid`/`end_regrid` bracket and the backend "rebuild" entry point as asks to Roles 1 and 2.
5. If `FaceDivFree` does not support ratio 4 or the single-cell direction, accept coarse-injection plus a post-regrid projection for those cases?
6. Species realizability (Role 4): clip after interpolation counts against the 1e-12 budget; is that the intended accounting?
