# Role 3 plan: regrid and multi-level transport (Batch 2, Phase 3)

Status: draft for Architect approval. No code is written until approved. Spec v0.4.28; ADR-001 v0.7 (K2 default, K1 per-kernel fallback), ADR-002 v1.1, ADR-003 v1.1. Base: branch FDS-AMReX at the merge of FireX 835588bb54 (bee11f0329). Effort figures are rough (no measured throughput exists).

## 1. Scope and proposed directory
In: `&AMR` input and input-time checks (IR-003, FR-010 ratio {2,4}, blocking factor, single-cell directions, FR-013); hierarchy from finer `&MESH` boxes and declared refinable region (IR-002, IR-008); tagging (FR-011); regrid with conservative transfer (FR-012); coarse-to-fine ghost fill, average-down and refluxing for mass and species (FR-016 data path, FR-024, FR-020/021); determinism (FR-015). Out: pressure (Role 2; Phase 3 uses a prescribed or simplified velocity per roadmap), species chemistry and realizability policy (Role 4, FR-025; I provide the hook), walls/OBST on levels (Phase 5), subcycling.

Directory: **`Source/regrid_transport/`** (disjoint from `Source/driver` and `Source/pressure_backend`). Same ownership pattern: `OWNERS.md` (only Role 3 edits it), `README.md`, `notes/`, `tests/`, C++ glue only (D-043, NFR-049). Edits to existing FDS files (read.f90 for `&AMR`, CMakeLists) go to the Architect as a patch list behind `#ifdef WITH_AMREX`. The Architect makes the first commit (OWNERS.md, README stub, empty CMake target); this plan then moves to `Source/regrid_transport/notes/plan.md`. Until then it stays uncommitted under `docs/`.

## 2. Design in one paragraph
Today `Level0` is one box per `&MESH`, the driver does not use AmrCore, and `Fields` holds a const reference to one `Level0` (m2a-gate-report 2.2). I propose `regrid_transport` owns an `AmrCore` subclass (AMReX class that holds the hierarchy: grids, distribution of grids to ranks, and the callbacks "make level / remake level / clear level / tag cells") and the level operators. Role 1 keeps the field registry and the kernel shim, generalised to a level index. Kernels stay Fortran (K2); I add no physics kernels except tagging, written in C++ unless a Fortran tagging kernel is cheaper (per-kernel K1 reason recorded if so).

## 3. Milestones (each has a checkable exit)
| M | Content | Exit test | Weeks |
|---|---|---|---|
| R0 | Interface agreement with Role 1 and Architect (section 4); skeleton, CMake target; `&AMR` parser (patch list); ratio, blocking-factor and single-cell-direction input checks | Unit tests: 3:1 and direction-dependent input rejected naming mesh pair; `race_test_1` remesh stops with `bf=8`, starts with `2 8`; `int_1to2` passes | 0.5 |
| R1 | Static hierarchy: AmrCore subclass, grids from finer `&MESH` boxes and region boxes, proper nesting, hierarchy dump | Dump equals input boxes after nesting/blocking (FR-010); debug nesting assertions (FR-013); tags outside region discarded and counted once per level (IR-008 R1-R5) | 1 |
| R2 | Multi-level fields and static two-level transport data path: per-level allocation, `FillPatch` (AMReX call that fills ghost cells from same-level and coarser data) with piecewise-constant interpolation for coarse-to-fine, volume-weighted `average_down` (fine over coarse) for cells and `average_down_faces` for faces; prescribed-velocity advection of rho and rho*Z across the interface | FR-016: coarse-to-fine ghost values of rho, rho*Z, TMP, MU, KRES, D, DS equal baseline on frozen dumps of `ns2d_16_int_1to2_refinement` (instrumented build, A-09b); decomposition independence (boxes, 1 and 4 ranks) | 1.5 |
| R3 | Refluxing (a flux register, `YAFluxRegister`, stores the mismatch between coarse-face flux and the sum of fine-face fluxes, then corrects the coarse cells next to the interface) for advective and diffusive scalar fluxes; global dt, all levels advance together | FR-020/021: composite mass and each species, covered coarse cells excluded, change only by boundary flux, to round-off, using the D-028 exact sum; run on a two-level periodic case with a blob crossing the interface | 1.5-2 |
| R4 | Dynamic regrid: tagging (FR-011 set), buffering, region clipping, regrid interval, transfer by field class (section 5), side-data rebuild hook | Per-criterion unit case: tagged cells = independent evaluation, cell by cell; per regrid composite mass and species change <= 1e-12 relative (FR-012); identical hierarchy on repeated runs and at 1/4 ranks (FR-015) | 2 |
| R5 | Regrid under flow: moving blob across several regrids, with and without refined-to-coarse transitions; divergence of the transferred velocity measured | Mass/species conservation across all regrids; velocity-divergence error after next projection <= FR-031 tolerance is checked when Role 2's composite path exists (Phase 4), reported only until then | 0.5-1 |
| R6 | Hardening: debug assertions, regrid cost versus step time (target <= 10% of step, roadmap Phase 4 criterion), load-balance weight hook (NFR-035), User Guide text for `&AMR` (IR-003) | Every `&AMR` parameter has a test (IR-003); out-of-range values rejected; timing note | 0.5-1 |

Total about 7-10 implementer-weeks (Batch 2 proposal says 6-10). Widest: R3 (diffusive flux extraction) and R4 (velocity transfer). Gate to Role 4: R3 done (species transport runs on levels).

## 4. Interfaces needed from Role 1 (Source/driver)
1. **Hierarchy-aware registry.** `Fields` and `TimeLoop` take a level index and a per-level BoxArray/DistributionMapping, not a fixed `Level0`; a `Fields` instance per level, constructible and destroyable at regrid. Keep ghost widths (D-031) per level.
2. **Per-level kernel loop.** `TimeLoop` stage functions callable per level on that level's MultiFabs, with global dt (MIN over levels) and the D-028 exact sums over uncovered cells only (fine-mask input).
3. **Face fluxes out of the mass/species kernel.** Expose the advective flux arrays `FX/FY/FZ` (per scalar; today "not registered", `is_per_box_scratch`) as per-box output views the transport layer can read to fill flux registers. Diffusive face fluxes are not exposed today (DEL_RHO_D_DEL_Z is a cell array built in divg.f90): need a kernel-side extraction, owner to be agreed (Role 1 shim, my consumer).
4. **Coarse-fine ghost hook.** A way to supply ghost values that did not come from same-level `FillBoundary`: the EXTERNAL_GHOSTS_FILLED path (patches 0003/0004) is the natural one; confirm it also covers cells at a coarse-fine face, and that `BcStep`'s OMESH fill is bypassed there.
5. **Side-data rebuild per level.** `SideData` rebuild callable from my make-level/remake-level callbacks (host is fine until Phase 11, D-047). Masks stay NotBuilt in this phase (no OBST).
6. **Field transfer classes.** Role 1 and I agree a table (extension of `inventory/mesh_fields.csv`): conserved state (RHO, ZZ as rho*Z), face velocity (U,V,W and US,VS,WS), pressure (H, HS), derived (TMP, D, DS, MU, KRES, ...). Derived fields are recomputed by the kernels after a regrid, not interpolated.
7. Test fixtures: reuse `tests/` scripts, `ExactSum.H`, field-dump format and `compare_run.py`.
Who owns the AmrCore subclass is an Architect decision; my proposal is Role 3 (section 2).

## 5. Regrid transfer rules (proposal)
- rho and rho*Z: new fine cells from the coarse level by limited conservative interpolation (AMReX `CellConservativeProtected` or `CellConservativeLinear` with limiting), where fine data already exists it is copied; coarse under fine is replaced by `average_down`. rho is recomputed as the species sum; any clip that fires is counted and logged (FR-025).
- Face velocity: `FaceDivFree` (AMReX interpolator that keeps the divergence of the coarse field on fine faces) as the starting choice, with `average_down_faces` for coarse faces under fine. FR-012 leaves the method TBD(Pressure Lead): open question 3.
- Scalar diagnostics: composite budgets with the fixed-point exact sum before and after every regrid, as a debug check on every regrid and a log line otherwise.

## 6. Tests
- **Single-level equivalence.** (a) `max_level=0` through the new path gives byte-identical field dumps to the M2a driver on `shunn3_32`, `csmag_32`, `shunn3_4mesh_32` (`tests/run_driver_tests.sh`, `run_decomp_check.sh`). (b) A level 1 covering the whole periodic domain, ratio 2, prescribed velocity: the level-1 result equals a uniform fine single-level run (`average_down` makes the coarse level follow it); level 0 equals the average-down of level 1 to round-off. (c) A refined patch with an identical-value field (uniform state) stays uniform to bit level across the interface.
- **Conservation.** Composite mass and species (exact sum, covered coarse cells excluded): per step with refluxing (FR-020/021 to round-off), per regrid (<= 1e-12 relative, FR-012). Negative control: switch refluxing off and show the budget error is no longer round-off (proves the test can fail).
- **Static two-level equivalence (FR-016, scalars only).** Ghost fill of transported fields against baseline on frozen data, as in R2. Velocity and H are excluded by spec.
- **Hierarchy.** Dump versus input; nesting, blocking factor, region clipping; decomposition independence (FR-015); tag criterion unit cases (FR-011); input rejection cases (FR-010, FR-004, IR-003).
- **Realizability** (FR-025, with Role 4): bounds and sum of Y after fill, regrid, average-down and reflux, in debug builds.
All tests are scripts under `Source/regrid_transport/tests/`, run on the 8-core development machine at 1 and 4 ranks, one thread.

## 7. Risks
Diffusive flux extraction needs a kernel-side change (R3). AMReX `FillPatch` and `average_down` assume a uniform grid index map; FDS face indexing is offset by one (Fields.H index map), so every call goes through a thin, unit-tested adapter. 2-D inputs need `ref_ratio_vect` with 1 in the single-cell direction and a y blocking factor of 1 (FR-010, provisional). Regrid rebuild cost is unmeasured (R-26). Prescribed velocity in Phase 3 means the transport tests do not exercise a projection at the interface: that is Phase 4.

## 8. Open questions for the Architect
1. Phase 3 entry is M2 in roadmap.md, the hiring proposal releases Batch 2 after M2a. Is R0-R3 (static, periodic, no walls) acceptable on M2a, with anything touching non-box level 0 or OBST waiting for M2?
2. Who owns the AmrCore subclass and `Level0` generalisation (me or Role 1)? IR-002 hierarchy-from-`&MESH` parsing: Role 1's `assemble_level0` currently rejects unequal cell sizes.
3. FR-012 velocity transfer method is TBD(Pressure Lead): accept `FaceDivFree` as the Phase 3 default?
4. Four M2a items are NOT YET (gate report 4.1-4.5, e.g. no official `csmag_32` baseline). Do my equivalence tests use `shunn3_32` and `shunn3_4mesh_32` only until the baseline exists?
5. Phase 3 velocity: prescribed field (my proposal) or the simplified projection the roadmap also allows?
