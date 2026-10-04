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
