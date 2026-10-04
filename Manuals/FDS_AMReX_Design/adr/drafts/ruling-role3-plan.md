# Rulings on the Role 3 plan (regrid and multi-level transport), 2026-10-02

Plan: `Source/regrid_transport/notes/plan.md`. Approved with the rulings below; the first commit of the directory is `ba3b5515e5`.

1. **AmrCore subclass:** Role 3 owns it and the level operators (FillPatch, average-down, reflux, regrid). Role 1 keeps the field registry and the kernel shim and generalises `Fields` and `TimeLoop` to a level index and a per-level BoxArray and DistributionMapping, after M2a closes. Role 1 and Role 3 agree the interface header in R0.
2. **Hierarchy parsing (IR-002):** Role 3 owns it. It provides a function that groups the `&MESH` entries by cell size and returns the level-0 meshes and the finer-level box lists. Role 1's `assemble_level0` receives only the level-0 meshes, so its rejection of unequal cell sizes stays as a guard.
3. **Sequencing:** R0 and R1 now (they do not depend on Role 1). R2 and R3 start once Role 1 delivers the per-level interface (items 1 to 5 of plan section 4), on periodic cases with no walls. Non-box level 0 and OBST work waits for M2.
4. **Face-velocity transfer (FR-012):** `FaceDivFree` is the Phase 3 default, provisional, behind a setting. The Pressure Solver Lead confirms or replaces it; R5 measures the divergence of the transferred velocity.
5. **Phase 3 velocity:** prescribed, analytic and constant in time (translation plus a divergence-free periodic field). No simplified projection in Phase 3; the projection at the interface is Phase 4.
6. **Diffusive face fluxes:** Role 1 implements the kernel-side extraction (an extra face-flux output in the generated kernel; the cell arrays stay bitwise unchanged, no change to FDS arithmetic). Role 3 specifies the array shapes and consumes them. The AMReX Integration Lead reviews. R3 may deliver advective-flux reflux first and diffusive second; if the extraction needs more than an added output, escalate to the Chief Architect.
7. **Baselines:** the V&V Lead's M2a baselines now exist (`vv-runs/baseline/gnu_ompi_firex-36975d7/`: `shunn3_32`, `csmag_32`, `shunn3_4mesh_32` and variants), so the equivalence tests may use them. Frozen two-level dumps still need the instrumented build (A-09b).
8. **Plan wording:** the tests run on the shared test machine, not "the 8-core development machine", in any document that is synced.

## Update 2026-10-02 (owner ruling D-050: single global dt, no subcycling)
Role 3 note: `docs/role3-numeric-scheme-note.md`. Reflux is dropped from R3; it becomes the **interface flux overwrite** of ADR-002 v1.2 (coarse face flux = area-sum of the fine face fluxes, per stage, advective and diffusive; average-down afterwards; overwrite-off negative control). Approved by the Chief Architect.
- **Ruling 6 is replaced.** Role 1 provides (a) the stage face fluxes (advective product face value times velocity, and the diffusive face fluxes of `divg.f90`) as readable arrays and (b) an override input for the coarse faces, applied before the flux divergence is formed. With an empty override set the result must be bitwise unchanged. The points are between `MASS_FINITE_DIFFERENCES` and `DENSITY`, and inside `divg.f90`. Role 3 specifies the shapes and the face lists; the Integration Lead reviews. Interim fallback for prototyping (post-kernel correction of covered-adjacent coarse cells) is allowed, not as the deliverable.
- **Second ghost layer at level jumps (Role 3 question b).** Phase 3 reproduces the FDS behaviour: both ghost layers by piecewise-constant coarse-to-fine fill (layer 2 equal to the layer-1 value, which is FDS's zero-gradient copy). FR-016 compares layer 1 against the baseline and the second layer is covered by the same fill, so the baseline comparison stays meaningful. Limited linear interpolation of layer 2 is an option behind a setting, evaluated with a convergence test after Phase 3, not before.
- **Time step.** The driver's `dt` is the minimum over all levels and ranks, the same in both stages (ADR-002 v1.2). Role 1's TimeLoop already takes the minimum over levels (plan section 4, item 2).
- **`&AMR` absent but finer `&MESH` entries present (patch list item 4).** It is an error that names the missing `&AMR` line and the mesh pair with unequal cell sizes. AMR mode is never inferred. Spec Lead to confirm against IR-003.
- **Patch list (`notes/readf90-cmake-patch-list.md`).** Item 2: no Fortran change; the C++ driver reads the input file itself. Items 1 and 3 (`read.f90` skip of the `&AMR` groups behind `WITH_AMREX`, CMake subdirectory): I apply them when R1 needs them in the driver build; Role 3 tells me then.
- **Spec wording.** The Spec & Program Lead changes "flux register / refluxing" to "interface flux overwrite" in D-023, FR-016, FR-020/021, FR-024, IR-003, roadmap Phase 3 and the related spec text.

## Update 2026-10-02 (b): rulings on the FDS mesh-interface note (`docs/role3-fds-mesh-interface-note.md`) and on R1 decisions
1. **Advective override velocity.** Approved. Use the single interface-face velocity (`UVW_SAVE` equal to it) in the multi-level path; the coarse face value is the area average of the fine faces. `MATCH_VELOCITY` is not called at level interfaces; same-level box faces are single-valued in AMReX. The override replaces the product (face value times velocity), never a product of averages. Overwrite off must reproduce FDS interface behaviour; that is the FR-016 baseline.
2. **Diffusive override.** Approved. Override = mean of the fine post-species-sum-correction fluxes, which equals FDS's area-weighted sum (`wall.f90:891-958`) and sums to zero over species. The `COARSE_MESH_IF`/`OMESH` branch is bypassed at level interfaces (Role 1 confirms in the hook design). The override needs the fine-level fluxes before the coarse divergence is formed, so the hook is a two-phase split (all-level fluxes, then override, then divergence), finest level first; the post-kernel correction of coarse cells stays prototyping only.
3. **Heat-flux override.** Dropped for Phase 3. Reproduce FDS (conduction unmatched; enthalpy diffusion from the overwritten species flux). Spec: mass and species conserved to round-off across levels, energy as in FDS. ADR-002 v1.2 text corrected accordingly.
4. **New fine interface faces at regrid.** Take the coarse face value (area average then equals the coarse value); interior fine faces use `FaceDivFree`. No 0.5-blend.
5. **Post-regrid projection.** Optional, default off in Phase 3. `DIVERGENCE_PART_1` runs on new levels, and a diagnostic reports max |div u - D| after each regrid.
6. **Level-0 meshes.** Approved: one box per `&MESH` in input order; the level-0 domain extent is subject to the non-box ruling N1a (padded to a multiple of bf0 for the domain, unpadded extents for `P_0`, snapping and output). Level-0 boxes need not be multiples of bf. Levels 1 and up are cut to `MAX_GRID_SIZE` and must be bf-aligned.
7. **Alignment.** A user-given fine `&MESH` that is not bf-aligned is an error and is never enlarged. Refinable-region boxes from `&AMR_REGION` are snapped outward to the blocking factor with a warning; the snapped region is the effective refinable region (it is what the dump, IR-008 and the ADR-004 output level use).
8. **`&AMR_REGION LEVEL`** = finest level allowed in that box (default `MAX_LEVEL`). Approved.
9. **Names.** `FACE_LINEAR`/`FACE_CONSERVATIVE` stay placeholders. Send the list of new namelist parameters and values to the Spec Lead, who owns the input names (IR-003).
10. **Patch list items 1 and 3** stay unapplied until R1 needs them in the driver build.


## Update 2026-10-02 (c): R4 design note rulings (D-058)

The R4 design note (tagging, dynamic regrid, data transfer) is accepted as the plan for R4, with these rulings on its six questions.

1. **Criteria.** Phase 3: temperature, density, species, HRRPUV and user boxes; vorticity next; OBST/VENT distance in Phase 5. Parameter names stay working names until the User Guide draft (IR-003).
2. **Initial hierarchy at t = 0.** Evaluate the input's initial conditions directly on each new level (repeat tag, make level, fill, until the finest allowed level), not interpolation from the coarse level, so features smaller than a coarse cell are resolved. After the last level is built, average down once so every covered coarse cell holds the fine average. The conservation budget (1e-12 per regrid) starts after that average-down. Role 1 provides the "initial fields on level L" entry point; in Phase 3 the prescribed analytic fields serve it. Interpolation from the coarse level remains the path for every regrid after t = 0.
3. **Tag criterion.** Accept undivided differences with `TAG_KEEP` hysteresis. No scaled gradient in Phase 3; it can be added later as an option if a case needs it.
4. **Interface asks.** Approved: `begin_regrid()` / `end_regrid()` on the registry (Role 1), the registry building `Fields` on an arbitrary BoxArray and DistributionMapping at run time (fine-level mesh objects: D-056, option B), and a "hierarchy changed, rebuild" entry point on the pressure backend that keeps no pointers to old MultiFabs (Role 2).
5. **FaceDivFree.** Verify ratio 4, 2-D and the single-cell direction first (it is test 4 and needs no Role 1 work). If it does not support a case, that case uses coarse injection plus the post-regrid projection, recorded in the design note as a limitation, with max abs(div u - D) printed after each regrid. In a single-cell direction the refinement ratio is 1 in that direction (2-D refinement is in-plane only).
6. **Species realizability.** Clipping is conservative: if the limited interpolation would give a negative or out-of-range value, the children of that parent are rescaled so their sum still equals the parent. The conservation budget therefore holds by construction. Every clip is counted (FR-025); the test on positive fields requires a count of zero, and a nonzero count in a case is reported as a warning with the level and the number of cells.

Also accepted from the note: finer `&MESH` boxes are force-tagged so a regrid never drops them; `RemakeLevel` called on a level whose grids are unchanged must be a bitwise copy (test it).

## Update 2026-10-02 (d): shared coarse ghost cells at corners (D-059)

In AMReX a covered coarse cell that is the ghost of two fine-patch faces (a patch corner, or a patch thinner than 4 coarse cells) holds one value, while FDS gives each mesh its own array. In the FR-016 comparison 4 of 16 corner-zone KRES cells differ; everything else agrees (fine side bitwise, coarse side to 3e-15 relative).

Ruling: **accept as a documented Phase 3 limitation (option a)**; per-face ghost storage (option b) is not built. Conditions: (1) the driver counts the covered cells that serve two or more faces and prints the count at each regrid; (2) the FR-016 report lists the corner-zone differences as the accepted exception, with the cause; (3) the blob tests (2-D and 3-D, patches with corners) must show mass and species conservation at round-off and no difference above round-off in the result against the uniform-fine run attributable to corner cells; (4) if any case at the Phase 3 gate shows an effect beyond round-off, the limitation is reopened and option (b) is costed for the fine kernels. The limitation is listed in the Phase 3 acceptance notes.

## Update 2026-10-02 (e): face prolongation (D-060)

Finding from the FaceDivFree probe (AMReX 26.09, 3-D): ratios 2, 4 and mixed keep each parent's divergence; ratio 1 in any direction aborts; interface faces get a slope instead of the coarse value.

Ruling: the driver's own normal-linear face prolongation replaces the "coarse injection plus post-regrid projection" fallback of Update (c), item 5. It handles ratio 1 (single-cell directions), keeps each parent's divergence to round-off, and the residual (variation of D inside a parent) is reported by max abs(div u - D) after each regrid. Interface faces keep the ruled rule (they take the coarse value), implemented by pre-fill plus mask. FaceDivFree stays in use where it applies only if the tests show identical behavior; otherwise one code path (the driver's) is used for all cases. Post-regrid projection stays off by default. The device (GPU) version is a later item; the Phase 3 host loops are accepted, with the GPU path listed in the Phase 4 backlog.

## Update 2026-10-03 (f): species transfer, post-regrid projection, fine-level solid plan (D-062 to D-064)

**1. Species transfer is mass-weighted (D-062).** Restriction, average-down and prolongation of species act on rho and rho*Z. Z (or Y) is derived from them; it is never averaged or interpolated linearly. Reason: linear Z averaging breaks species mass conservation (about 1e-4 relative in the earlier run; the negative control with linear Z breaks it by 1e-1). Implemented in Role 3 commit 7c85f23539 (`average_down_registry` and the coarse-fine ghost hook use it; rho, ZZ and ZZS on covered cells agree bitwise with Role 1's RegistryTransfer). The test has a negative control with linear Z, and the test fails against it.

**2. Post-regrid composite projection (D-063).** Finding from the dynamic-regrid blob tests: where a regrid creates new fine faces next to retained old ones, max abs(div u - D) at the junction is of velocity-gradient size (1.0 in 3-D, 3.5 at ratio 4); all-new cells match the parent divergence to about 1e-14. This is not acceptable as a final state, because the next pressure RHS assumes div u = D.

Ruling:
- After such a regrid the driver runs a composite projection: solve Lap phi = div u - D on the whole hierarchy, u -= grad phi, and fine fluxes overwrite coarse.
- Acceptance: max abs(div u - D) <= pressure tolerance on uncovered cells, with a negative control (projection off fails).
- Re-prolonging a band of one parent cell around the retained faces is rejected.
- Interim: the hook `project_after_regrid(levels, D)` is built with a mock solver. Until the composite MLMG lands (Role 2; its entry accepts an arbitrary per-level RHS) the test reports the number; once it lands the test asserts.
- A regrid without projection support is not a gate case.
- This changes the last sentence of Update (e) ("Post-regrid projection stays off by default") for regrids that retain faces. Where no old faces are retained, the projection is not needed.

A GPU-enabled AMReX build is needed to use the offload tagging kernels inside the C++ library; this stays a Phase 4 backlog item.

**3. Fine-level solid phase (D-064, plan `solid/05-fine-level-solid-plan.md`).**
1. The CALCULATE_ZZ_F mass-flux block (`wall.f90` 1185-1320) may be extracted verbatim as a record-local routine by a guarded WITH_AMREX patch. It stays a DRAFT until oneAPI and GNU Debug validation, and until the Legacy Mapper confirms that the block has no hidden dependence on owner loop state. A6 is extended to cover it.
2. Face key = (level, global integer index of the gas-side cell at that level's resolution, IOR in +-1..3), independent of mesh, box and rank. Sum order is deterministic: (level, k, j, i, IOR). The Spec Lead aligns FR-046 with this key.
3. One accessor returns the mesh object of (NM, level) via POINT_TO_BOX (D-056). It is introduced by a guarded local patch at HEAT_TRANSFER_COEFFICIENT (`func.f90` 3156) and at the back-side MESHES(NM) lookups. It is not an upstream patch.
4. Role 1 builds the fine-level OBST wall tables in Phase 5; Role 3 triggers the rebuild at regrid.
5. Area-mean TMP_F against the FR-022 budget is deferred until the WP6 measurement.
6. Thin-wall faces: records are never deeper than their owner.

## Update 2026-10-04 (g): solid wall/BC translation and GPU spike plan review (D-065, D-066)

**1. Solid wall/BC translation (D-065, plan `amrex/wall-bc-translation-plan.md`).**
- Q1: the `EWC%NIC>1` branch (INTERPOLATED_BC species flux match, L1485) is retired in the AMR route behind a host abort guard. This is conditional on a grep of the supported AMR input set showing no `NIC>1` case (no INTERPOLATED_BC between meshes of different resolution). Refused inputs are listed as FDS-only. Consistent with D-055.
- Q2: the back-wall heat transfer coefficient in the thick-wall and thin-wall passes uses a snapshot of the other side in the AMR route. FDS-only mode is untouched. The host AMR route also offers snapshot mode, so the difference is measured as algorithm and not as port. V&V signs the tolerances (`back_wall_test`, `heat_conduction_a`).
- Q4: running the whole wall pass on the host when `HVAC_SOLVE` is on is acceptable in Phase 4. HVAC cases are correctness gates only, not GPU performance gates. This is a known limitation; its cost is measured, and it is reopened if an HVAC case becomes a performance target.
- Q6: the neighbour obstruction-mass read returns the value at the start of WALL_BC of that stage (snapshot before the pass, read-only during it), identical on host and device. The difference from FDS order is recorded as for Q2. Burn-away is deferred (FR-042).

**2. GPU spike plan review (D-066, plan `amrex/stage1-gpu-spike-plan.md`).**
- Wall-seam work packages WP1b, WP9, WP9b, WP11 and WP12 are approved with changes; WP10 is approved.
- Flux-hook findings are approved with changes. Finding 1: gather straight from `ADV_F*`, no pack kernel.
- Data-movement infrastructure kernels (pack, gather, scatter, checksum) in the C++ driver layer may be C++ `ParallelFor` (K1) because they contain no physics. All physics kernels stay K2. This clarifies D-049; ADR-001 gets one sentence.
- The per-box wrapper for the 1-box and 4-box device tests lives in test code only (D-052).
- A-57 tooling (K2 CI check, kernel lint, zone-sum order script, `port_kernel_map`) must exist before the first new kernel is accepted. The owner is the Mesh Data Loops Engineer.
- The plan total is relabelled 55-68 work-days with ESTIMATE labels.
- The csmag_32 RHOS difference is not a blocker for the work packages, but it blocks FDS-agreement claims at edges and corners.

## Update 2026-10-04 (h): pressure gauge and mean removal (D-067)

D-067 pressure gauge and mean removal (source: `pressure/06-meanremoval-masked-stretched.md` sections 5-6).
- Default gauge: `sum(rho*V*(KRES-H)) = 0` per pressure zone and per connected component, with the exact (decomposition-independent) sum over uncovered cells (ULMAT/UGLMAT convention). It is always applied, because `p = rho*(H-KRES)` carries the constant into the baroclinic pass-2 RHS. It is also applied to FFT-solved cases, with PRES compared after removing a constant.
- Default mean removal: the composite volume-weighted mean (D-032). The FDS arithmetic removal is kept as a runtime parity switch (AMR mode is uniform-per-level only, so the stretched-grid advantage matters for parity studies against FDS-only runs). AMReX native `makeSolvable` is not relied on.
- The FDS single-rank GLMAT gauge defect (the one-rank periodic test 7 leaves a constant of +0.0304 and moves the mass-fraction slice by up to 1.5e-5 in 0.05 s) goes to the owner as an upstream patch file under D-051, not as a local change.
- The composite `gauge_weight`/offset is an open implementation item for the Pressure Backend role.
- Intel validation: draft patches 0005-0009, including the 0007 `DT_NEW(*)` amendment (`a72d491fb9`), pass oneAPI validation (0 warnings, OFF bitwise, driver tests, decomposition check with a max(2e-15, 1 ulp) gate, outputs check against a same-compiler reference). GNU Debug validation is pending. The kernelcheck driver-fill csmag DIV1 difference is an open item under bisection.

## Update 2026-10-04 (i): A-57 tooling, GPU spike plan follow-ups, patch validation, projection setting (D-068)

D-068 A-57 tooling and GPU spike plan follow-ups.
- Accepted: the tools in `docs/tools/` (`k2_ci_check.py`, `kernel_lint.py`, `zone_sum_order.py`, `port_kernel_map.py`, `kernel_registry.toml`).
- `AMReX_CUDA_FASTMATH=OFF` is pinned in the driver CUDA configure with FORCE, and the configure fails if it is ON. The lint checks the pin.
- `reduction(max:)` and `reduction(min:)` are added to the D-029 clause list (exact, order independent; a kernel using them counts as ported only after a device run). `reduction(+)` is forbidden except through the D-053 zone-sum order.
- Waivers: the hand-written S4 prototypes, the taskwait and do-concurrent variants and the 19 K1 sites are category `prototype` (evaluation code, expiry when replaced by generator output, not in production builds). The waivers for the rule-7 violations `s4k2_face_values` and `s4k2_clip_terms` are rejected: fix them, or keep them prototype-only with an expiry. `--strict` must pass for production kernels.
- The kernel header-note definition is accepted: a provenance comment `file.f90:a-b` above the kernel, plus passive-scalar and cylindrical notes for hand-written files.
- Design constraint: "streams per box" means host threads per box for K2 launches. K2 launches block and do not wait for AMReX streams.
- Patches 0005 to 0009 pass oneAPI validation and GNU Debug validation (amended 0007 is `a72d491fb9`; the fine-b shadow check runs on a single level-0 mesh and 1 rank only). The gate on level>0 physics is lifted.
- Setting for D-063: `POST_REGRID_PROJECTION = AUTO (default) | ON | OFF`. AUTO projects when some fine level both keeps old cells and gains new ones. ON always projects. OFF is diagnostics-only: it prints a warning in the output and the run is not a gate run.

## Update 2026-10-04 (j): Intel build flags for bitwise comparison with the driver (D-069)

D-069 Intel build flags for bitwise comparison with the driver.
- FDS ghost values of TMP, RSUM, ZZ and ZZS are computed in ASSIGN_GHOST_VALUE (`RHO_ZZ_OTHER_2/RHO_OTHER_2` and `PBAR_P_2/(RSUM_TMP*RHOP)`) and are not bit-copies of the interior values. The ifx default approximate division makes DENS, DENSCLIP and DIV1 differ against the driver's exact-copy ghost fill.
- The fix is flags only: build the Fortran of both the reference FDS and the driver with `-O2 -prec-div` (or `-fp-model=precise`), with identical Fortran flags on both sides. The C++ keeps `-fp-model=precise` (compiler `mpiicpx`; the icx default reassociation breaks `tile_race`).
- Kernelcheck references must come from a `-prec-div` ifx FDS build.
- The earlier csmag DIV1 difference was a reference-dump reconstruction error (`FDSREF_STEPS` must be 2,3 to match the recorded dump), not a compiler effect.
- Remaining known: the plain (no-BC) full/face VISC/VFLUX snapshot effect (the dump is not the exact pre-boundary-step state; `+strips` removes it), the same as with gfortran. The decomposition check `p1_div_DDDT` differs by 1 ulp, covered by the `max(2e-15, 1 ulp)` gate.
- D-063 context: Role 3 driver-level regrid conservation (R2b): advection with interface flux overwrite closes to 2e-14 (off: 4.3e-5); the diffusion corrector residual of 4.4e-7 (scales with dt^2) is open and not acceptable for the Phase 3 gate.

## Update 2026-10-04 (k): libm tolerance and pressure-backend/gate-threshold rulings (D-070, D-071)

D-070 libm tolerance for kernels that call transcendental functions.
- Kernels whose arithmetic is only `+ - * /` and `sqrt` stay bitwise (with the +0/-0 rule).
- Kernels that call libm transcendentals (`**` with a non-integer exponent, `exp`, `log`, trigonometric functions) get the category `libm` in the kernel registry, with a per-value tolerance of 2 ulp (measured: 1). First case: `cfl_wall_max` / UVWMAX with `(ABS(Q)/RHO)**(1/3)`; the device result differs from the host by 1 ulp for about 13% of the arguments (3 of 13 scenarios differ in the last bit).
- No own device power routine: the host libraries (gfortran libm, Intel libm) already differ from each other, so reproducing one host library bit for bit buys nothing.
- Dt-coupled quantities (UVWMAX feeds dt) use the run-level tolerance already used for cross-compiler runs.
- The K2 CI check flags every libm call, so none enters unlisted.

D-071 Pressure-backend and gate-threshold rulings.
- (a) TWO_D contract: a Dirichlet face in a one-cell direction is treated as Neumann and the term is dropped (as FDS TWO_D does).
- (b) Fold sign convention: low Neumann `rhs += g/h`, high Neumann `rhs -= g/h`, Dirichlet `rhs -= 2*H_b/h^2`. The Pressure Lead confirms it with a nonzero-data test against FDS H.
- (c) Hierarchy fold helper: wait for a consumer.
- (d) The 23 inputs with Dirichlet in a one-cell x or z direction stay refused in AMR mode (FDS-only, D-057).
- (e) Mixed Neumann/Dirichlet on a hierarchy is accepted if the order (about 2) and the true residual are right; where the fine-level excess error sits goes to the Pressure Lead.
- (f) A-58 provisional thresholds confirmed. P3-F02 corner limits are a gate: fine side none; KRES coarse differing cells only in the declared D-059 set; RHO/TMP coarse <= 3e-15; if exceeded, reopen D-059, do not widen. P3-B09 (2): corner max <= edge max + 4 ulp. P3-B09 (1): the propagation distance is derived from stages times stencil half-width per step instead of n+1. The 0.5 discrimination rule `||AMR-F|| <= 0.5||C-F||` on the tracer slice is a gate for P3-B02 and P3-B07 and report-only for P3-B08 (ratio 4). P3-R02 pressure tolerance: working bound 1e-9*U/dx_fine at solver tolerance 1e-12 (the Pressure Lead confirms the mapping).
- (g) Fine-box domain-edge defect: a fine box with faces on the physical domain edge (for example a one-cell y direction in 2-D) had unfilled ghost layers, giving a spurious diffusive flux (about -5e-5), which explains the E3b residual 4.4e-7 and the E1 ZZ gap 3.9e-5. Ruling: fine boxes get the same physical-domain boundary treatment as level 0 on every domain-edge face (wall cells or an identical mirror/periodic ghost rule). Phase 3 scope is periodic and Neumann/mirror domain faces; other boundary types abort with a clear message. Role 1 owns it; Role 3 tightens MASS_TOL and E1_ZZ_TOL to round-off afterwards.
- Cause of the earlier vcorr 1-ulp difference: FVX/FVZ change inside the pressure step between the dump and the corrector, so dump-based comparisons use corrector-time inputs.

## Update 2026-10-04 (l): kernel tooling, solid-phase and radiation sweep rulings (D-072, D-073)

D-072 kernel tooling and solid-phase rulings.
- (a) libm function list for the K2 check confirmed wide: `**` with a non-integer exponent, `exp`, `log`, `log10`, trigonometric, hyperbolic, inverse trigonometric, `erf`, `gamma`, Bessel, `hypot`. An integer exponent is judged by the value as written in the source (a literal or an integer-typed exponent, such as `x**2`, `x**2._EB`, `x**3`). The source should use explicit multiplication or an integer exponent so that no `pow` is emitted; `pow` emission at low optimisation is a build-flag check.
- (b) The 2 ulp comparison gate belongs to the V&V device-tier harness; the registry only records the category.
- (c) A kernel is "ported" when a device run is on record on real fields within its class tolerance (bitwise, or 2 ulp for libm), the kernel owner has reviewed it, and it passes `ci_checks --strict`. The Integration Lead sets the status in the map.
- (d) ZS-06 (wall-list kernels accumulating into the gas cell in list order, one thread per gas cell) is accepted as a NOTE; see `docs/solid/07-sp2-sp3-kernel-review.md`. The rule-7 pair `s4k2_face_values` / `s4k2_clip_terms` stays prototype-only, with an expiry.
- (e) Solid-phase 1-D wall solve: the sub-step cap is 10**6, with status 301 for non-finite input only (the harness asserts that it never triggers on finite gate cases). The D-065 Q2 back-side snapshot covers only the four fields read at `wall.f90:2163-2183`. The flip budget between host and device decisions is measured and gated: frozen inputs per call, decisions exact, then values at 2 ulp; flipped calls are compared at 1e-10 relative on temperature and heat flux; the working hypothesis is 1 flip per 10**4 record-calls.
- (f) D-065 Q1 condition is a design requirement: the AMR route never builds FDS EXTERNAL_WALLs at a coarse/fine level jump (the driver handles the jump by ghost fill and the interface flux overwrite), so the `NIC>1` branches (`wall.f90:885`, `divg.f90:210`) are not reached. The seven multi-mesh level-0 inputs with different-resolution meshes that could reach them stay FDS-only, or are refused in AMR mode, unless converted to refinement levels. The Legacy Mapper verifies and lists their AMR-mode status.
- (g) Kernelcheck: the csmag `DIV1` mismatches (Intel and GNU, `csmag_32` and `csmag_32_fishpak`) were reference-dump errors: the dumps were generated with steps 1,2, and the correct setting is `FDSREF_STEPS=2,3`. `DIV1` is bitwise in all modes after regeneration. The plain full/face `VISC`/`VFLUX` snapshot effect remains documented.

D-073 radiation sweep kernel design (`docs/radiation/06-fr062-sweep-kernel-design.md`).
- (a) Plane-parallel across boxes, serial over angles, is the baseline. Angle batching comes only after measurement, with ordered chains and B=1 for mirror and cylindrical.
- (b) Ghost-cell `UII`/`UIID` may be not bit-reproduced if V&V confirms that nothing reads them.
- (c) Box split versus multi-mesh FDS byte-identity (2 and 8 boxes, 1/2/4 ranks) is an accepted T0 test.
- (d) A hand-written K2 sweep kernel is acceptable if the cell body is copied verbatim by line range, a drift check against `radi.f90` runs in CI, the K2 check and lint pass, the registry entry is tagged hand-written/FR-062, and there is no FMA contraction. The sweep is registered as claim L1242 with the hand-written tag.
- (e) Cylindrical geometry is out of the first release (FDS-only in AMR mode).
- (f) Box-level launch batching is a driver-layer item (Integration Lead).
- (g) The Radiation Lead R1 sign-off is applied in `blocked-loop-families.md`.

## Update 2026-10-04 (m): FR-003 scope and pressure results (D-074, D-075)

D-074 FR-003 scope for the pressure-code-0 inputs (A-62).
- (a) The 23 inputs with pressure code 0 (21 `soborot_*`, `bound_test_1`, `bound_test_2`; one-cell Dirichlet, no pressure solve) are FDS-only in every mode.
- (b) `fds_amr` aborts on them with a clear message that names the FDS executable. No no-pressure-solve mode is added to the driver; revisit only if the owner wants them in the AMR executable.
- (c) FR-003 coverage is the verification set minus these 23, listed in the spec. No verification input has a one-cell x or z direction (wording confirmed with the V&V Lead).

D-075 pressure results and open answers (`docs/pressure/07` sections 11 and 12).
- (a) Fold sign convention confirmed (low Neumann `rhs += g/h`, high Neumann `rhs -= g/h`, Dirichlet `rhs -= 2*H_b/h^2`, `phi = H`): 24 of 24 cases against an FDS-style dense matrix, worst 4.4e-14, each sign flip gives 0.31 to 3.6 relative difference. No FDS run with nonzero wall data and an H dump exists, so A-09b stays open.
- (b) Mixed Neumann/Dirichlet fine-level excess is a smooth domain-wide discretisation effect (about 85% from lost coarse truncation-error cancellation inside the patch, 13 to 15% from the coarse/fine interface). The harness checks order per region and a constant-in-n ratio to the uniform fine error. Non-matching Dirichlet wall data is untested.
- (c) P3-R02 stays gated at `1e-9*U/dx_fine` at solver tolerance 1e-12; measured values (7,000 to 62,000 times below) are reported alongside; no tightening until a second case family confirms.
- (d) Slow stretched-grid HYPRE bottom solve: AMReX hands HYPRE a row-scaled non-symmetric matrix while PCG is selected. Workaround `hypre_solver=BiCGSTAB`, bottom tolerance 1e-11. A symmetric-scaling change in `habec_ijmat` is a candidate upstream patch, delivered as a patch file for the owner after AMReX tests.
- (e) Role 2 HYPRE backend: Krylov defaults accepted (PCG with BoomerAMG for one level, GMRES(30) for a hierarchy); pin row excluded from the residual check and reported separately; `residual_tol` not relaxed; backend option stays on `PressureOptions`.
- (f) Solid-phase face-write table: a table refused by `face_write_check` makes the driver abort with the report (owner: Chief Architect with the Wall Loops Engineer).
- (g) Radiation batched launch layer (one launch per wavefront plane across all boxes of a level): Integration Lead (D-073 (f)).

## Update 2026-10-04 (n): level-jump EXTERNAL_WALL check (D-076)

- (a) Fine boxes have no wall cells at a level jump; the jump is handled by ghost fill and flux overwrite. The NIC>1 loops (`wall.f90:897`, `divg.f90:219`) run over zero walls on a fine box.
- (b) The whole route is not provable by reading. It is proved at run time by the D-065 Q1 abort guard at the end of `INITIALIZE_MESH_EXCHANGE_1` (new numbered driver patch, Legacy Mapper, validated on both compilers) and a level-0 NIC=1 assertion test (V&V Lead).
- (c) Finer `&MESH` lines are not filtered in `read.f90`. The driver gets an input pre-pass (converter) that writes a level-0-only input and builds the `&AMR` hierarchy from the removed meshes. `main.cpp` calls `parse_amr_params` and the hierarchy builder; `assemble_level0` gets only level-0 meshes. Owner: Role 3 with Role 1, after the domain-edge fix. Until then, inputs with finer meshes abort at `FdsAmr.cpp:30` with a message that names the converter.
- (d) Classification: race_test_1/4 and the four derived 2:1 cases run via `_r4` copies or converter output (originals FDS-only); three inputs convert to refinement levels; `duct_flow_uglmat_refine` waits for FR-040 R3; stretched-grid and embedded inputs are FDS-only. The Legacy Mapper's list is the reference.
- (e) FM_Burner 5 mm inputs have a coarse/fine face with 2 or more species and are FDS-only in AMR mode until converted. The stale header comment at `LevelRegistry.H:8-11` is for its owner to fix.

## Update 2026-10-04 (o): flip-budget gate, ported list, test lock, libm entry (D-077)

- (a) The flip-budget gate in `solid/09` is ratified: 95% confidence per run and 99% at phase exit; PASS, FAIL and INCONCLUSIVE, with INCONCLUSIVE not signing; the denominator is exposed distinct-input calls; Level A (device equals forced host within 2 ulp) and Level B (`E_class` measured in a host forced-flip survey, signed by the Solid Phase Lead and the V&V Lead); shortfall rule and overrun responses as in D-077 (1).
- (b) The 1-D wall solve cap is 10**6 sub-steps with the predicate form: at the cap, finite state and `DT_BC` continue as upstream, otherwise status 301. The AMR CPU path carries the same cap and text. The driver aborts with the report.
- (c) Patch 0010 (NIC>1 guard) validated on both compilers and tests P-3 and N-5 passing on the host build come before the L1485 retirement is merged.
- (d) The V&V Lead is reviewer of record for `ported` status; the Architect co-signs `ported.toml` with the evidence scope per entry; `ported` counts for the map only until `ci_checks --strict` is green.
- (e) Shared generator test lock: one full pass per engineer per hour, FIFO under `flock`, 30-minute hold limit; scoped runs need no full lock.
- (f) `CUNNINGHAM` is added to the libm registry class (2 ulp, not time-step coupled) by the Species and Combustion Lead after a recorded grep check.
- (g) Ghost `UII`/`UIID`/`QR` are compared interior-only in restart-file checks. The `-gpu=nofma` build and the 8-thread sweep run on the GPU test machine.

## Update 2026-10-04 (p): libm entry for the wall viscosity, gsfv co-sign, Role 1 milestone, patch 0010 on Intel (D-078)

- (a) `wall_visc_les` (`velo.f90:306-351`, `exp` and `pow` with 2.5, 1.25 and 1.5; 2 ulp; `dt_coupled` true) is approved for the libm registry. `dt_coupled` kernels are compared at class tolerance at run level (time-step sequences agree to round-off, not bitwise). The Wall Loops Engineer appends the entry and adds the check.
- (b) The six `gsfv_*` kernels are co-signed as ported on the V&V review (`vv/ported-review.md` section 9); they count for the map only until `ci_checks --strict` is green; the limits stay in each note.
- (c) Domain-edge fix accepted. The new fine level takes `RSUM`, `MU`, `KRES`, `D`, `DS`, `H` and `HS` from the parent by `RegistryTransfer::derive`. The skip of `fill_omesh` and `VELOCITY_BC` on level>0 is an interim; the cause is an open item before the Phase 3 gate. For Role 3: the tolerances `MASS_TOL` and `E1_ZZ_TOL` can now be set to round-off against the first two-level numbers (mass and rho*Z drift 1.3e-14 over 40 steps, max abs(div u - D) 1.8e-12). Open: overwrite-off control, multi-box fine level, 4-rank two-level run.
- (d) Patch 0010 passes on Intel. The N-5 control needs an input that passes the D-076 converter, or a scratch build that bypasses it. Tests key on messages, not exit status.
