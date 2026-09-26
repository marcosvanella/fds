# Architecture rulings FR-040 R3 (thin OBST across a coarse/fine interface) and FR-034 (leak area across levels)

| Field | Value |
|---|---|
| Status | **ACCEPTED (final), AMR Chief Architect, 2026-09-25.** Incorporates the AMR Pressure Solver Lead's answers and corrections (`pressure/03-fr034-fr040r3-answers.md`). Changes from the draft are listed in §5. The requirements text is in §6 for the Spec Lead (same batch as the IR-002 and FR-041a wording fixes). Folding these rulings into ADR-002 / ADR-003 is a later editorial step with no change of content. Amends D-009 (band rule R3-T) and adds to D-032 (removed-mean diagnostic); D-028 unchanged. |
| Date | Draft 2026-09-25; final 2026-09-25 |
| Source pin | All `*.f90` citations are **FireX 36975d765f** (`Source`, read-only), per requirements.md:16 (D-034 citation rule). FireX 36975d7 is the base of branch `FDS-AMReX` (local, and on `marcosvanella/fds`). AMReX citations are commit `99ddfda` (`(local AMReX checkout)`). All Verification-input citations are `Verification/...`. Docs are cited relative to `Manuals/FDS_AMReX_Design/`. |
| Answers | requirements.md:226 (FR-040 "Open point R3"); requirements.md:180 (FR-034 "Open point: leak area across levels"); vv/scope_alignment.md:132-133, 142-143. |
| Owners | Decided: AMR Chief Architect. Reviewed: AMR Pressure Solver Lead (Ruling 2; E-2 and interpolation items of Ruling 1). Actions: AMR Spec & Program Lead (§6 text into requirements.md); AMReX Integration Lead (exact-sum scale, zone numbering and seeding, §2.2 (g)); AMR V&V Lead (acceptance tests §1.6, §2.5; derived variants under A-41); FDS Legacy Mapper (§3.3, §3.2 Q4 [VERIFY]). |

Every claim about FDS behaviour below was read in the FireX source at the cited lines. Statements marked **[ruling]** are the binding design decisions of this document. Statements marked **[derived]** are arithmetic from the cited input lines, not runs. AMReX claims were read at `99ddfda` by the Pressure Solver Lead. No runs were made for this ruling.

---

## 0. Background facts used by both rulings

### 0.1 How FDS snaps an OBST, per mesh (READ_OBST)

- The OBST is processed separately for every mesh it intersects (per-mesh loop, read.f90:10991). It is clipped to the mesh (read.f90:11127-11139) and snapped to the nearest face of *that* mesh with `NINT` (read.f90:11157-11162).
- `THICKEN` to one cell: read.f90:11167-11173 (and y/z analogues). Collapse below 0.25 cell: read.f90:11174-11216 (includes the A-20 slips at 11176, 11193, 11210). Discard if zero-thickness in ≥ 2 directions: read.f90:11220-11224.
- `OB%THIN` is set when the snapped OBST has `I1==I2` (etc.) and the undivided input length is < 0.5 cell (read.f90:11315-11317). A zero-thickness input (`XB=1.0,1.0,...`) is therefore THIN on every mesh.
- A zero-thickness OBST blocks no cells: `BLOCK_CELL` fills cells `I1+1..I2` (read.f90:11708; func.f90:5503-5526), so `I1==I2` blocks nothing and the obstruction exists only as faces.
- Near a mesh boundary, an OBST lying just outside the mesh (within half a ghost cell) is collapsed onto the mesh boundary (read.f90:11041-11123). The ghost width comes from the neighbour mesh found by `SEARCH_OTHER_MESHES` only `IF (NOM>0 .AND. PROCESS(NM)==MY_RANK)` (read.f90:11044-11048 and analogues); otherwise the mesh's own `DX(0)`/`DX(IBP1)` is used. See §3.3.
- An OBST entirely inside a lower-numbered (finer) mesh is only made invisible in Smokeview (`COLOR_INDICATOR=-2`, read.f90:11461-11465). It still exists on the coarse mesh.

### 0.2 How FDS represents a thin OBST as wall faces

- Ordinary OBST face loops in `INITIALIZE_MESH_VARIABLES_1` create wall cells on exposed faces (init.f90:195-300). For `I1==I2` the `IOR=-1` face gets a wall cell in cell `I1` and the `IOR=+1` face gets one in cell `I1+1`: **one geometric face, two wall cells, one per gas side**. A face that points out of the mesh gets no wall cell and sets `EXPOSED_FACE_INDEX` instead (init.f90:199, 221). If a face already carries an external wall cell, the OBST reuses it (init.f90:204-211).
- `WC%THIN = .TRUE.` when the cells on both sides are gas and not exterior (init.f90:3315). This, not the `THIN_WALL` array, is the general thin-face representation.
- The `THIN_WALL` cells of `INIT_THIN_WALL_CELL` are created only when `ANY(SURFACE(OB%SURF_INDEX)%HT_DIM>1)`, i.e. HT3D (init.f90:118-191). HT3D is deferred (D-033), so this path is out of scope.
- `B1%AREA` of a wall cell is the geometric face area (init.f90:3060, 3066, 3088, 3105, stored at 3380). Thin-OBST wall coordinates come from the OBST coordinates (init.f90:3368-3375).
- ERROR(421): a SURF with `U_NORMAL` or a `LEAK_PATH` on a face whose back cell is gas and that belongs to an OBST (a thin obstruction) stops the run. ERROR(422): the same for an HVAC VENT (init.f90:1102-1115, in `WALL_LOOP_0`, init.f90:1080). For external wall cells the back cell is looked up in the neighbour mesh (init.f90:1089-1096). **So LEAK_PATH and HVAC faces never sit on thin OBSTs.** This matters for Ruling 2.

### 0.3 How FDS joins meshes of different resolution

- `INIT_WALL_CELL` finds the neighbour mesh and its covering cell range `IIO..KKO` (init.f90:3145-3180) and stops with ERROR(431) if fine faces do not tile the coarse face within `ALIGNMENT_TOLERANCE` (init.f90:3184-3226). An interface face is `INTERPOLATED` if both sides are gas; if either side is solid it takes the solid OBST's SURF (init.f90:3251-3262). The `IIO..KKO` range is stored per external wall cell (init.f90:3317-3326).
- Ghost values across the interface are averages over the `IIO..KKO` range, e.g. `H` in `NO_FLUX` (velo.f90:1388-1399); scalars are area-averaged (wall.f90:315-347, cited at requirements.md:142).
- With `UGLMAT`, the matrix setup enforces "solid wins" across the interface (pres.f90:3662-3902): a coarse INTERPOLATED face with at least one SOLID and one INTERPOLATED fine subface becomes SOLID (pres.f90:3726-3766); for a SOLID coarse face the INTERPOLATED fine subfaces become SOLID (pres.f90:3780-3800); a fine INTERPOLATED face whose coarse face is SOLID becomes SOLID (pres.f90:3801-3829). This is the FDS precedent for C/F face consistency; FDS applies it silently and only to the pressure matrix.
- Pressure-zone flood fill (`ASSIGN_PRESSURE_ZONE`, func.f90:5558-5729) does not cross a thin OBST face (func.f90:5656-5676). A thin wall that leaks on one level therefore also merges zones on that level.
- Cells covered by a lower-numbered (finer) mesh are flagged in `INTERPOLATED_MESH` (init.f90:1564-1587). `DSUM`/`PSUM` skip them (divg.f90:727-753) and `USUM` skips wall cells whose gas cell is covered (divg.f90:760).

---

## 1. Ruling FR-040 R3-T: thin OBSTs that cross or touch a coarse/fine interface

### 1.1 The test case, read from the input (`Pressure_Solver/duct_flow_uglmat_refine.fds`)

- Meshes: 4 coarse meshes of 16³, dx = 0.2 m, as a 2×2 tiling in x/y over z ∈ [-0.2, 3.0] (lines 5-6, `MULT` with `K_LOWER=0, K_UPPER=0`); 4 fine meshes of 32³, dx = 0.1 m, over z ∈ [3.0, 6.2] (lines 8-9, `K_LOWER=1, K_UPPER=1`). Line 11 (`MESH IJK=32,32,32 ...`) has no `&` and is a comment. `SOLVER='UGLMAT HYPRE'`, `MAX_PRESSURE_ITERATIONS=1` (line 17). `STRATIFICATION=.FALSE.` (line 16). `T_END=60` (line 14). Run with `-p 8` (Verification/FDS_Cases.sh:450).
- `SURF ID='DUCT'` has no `MATL_ID` (line 19), so all duct walls are **stateless** (FR-041a definition, requirements.md:230-231).
- HVAC: a 0.1 m OBST at x 1.0-1.1, y 1-2, z 1-2 (line 21) carries `SUCK`/`BLOW` vents joined by a fan duct, `VOLUME_FLOW=1` (lines 22-26). All of this is in the coarse region.
- Thin walls that cross z = 3.0 (lines 43-46):
  - line 43: `XB=1.0,1.0,4.0,5.0,1.0,5.0` (x = 1.0, z 1-5)
  - line 44: `XB=2.0,2.0,4.0,5.0,2.0,4.0` (x = 2.0, z 2-4)
  - line 45: `XB=1.0,2.0,4.0,4.0,2.0,4.0` (y = 4.0, z 2-4)
  - line 46: `XB=1.0,2.0,5.0,5.0,2.0,4.0` (y = 5.0, z 2-4)

  Together they form a vertical riser x ∈ [1,2], y ∈ [4,5] that carries the whole duct flow up through the interface.
- Devices: `flow_in` and `flow_out` are area integrals of U at x = 1.0 (lines 84-85). The dataplot row expects 1 m³/s (reference `Pressure_Solver/duct_flow.csv`), relative error 0.05 at the end (Utilities/Python/FDS_verification_dataplot_inputs.csv:161).
- **[derived]** Every wall coordinate lies on a coarse face and on a fine face. With the level-0 origin at -0.2: x = 1.0 is coarse face 6 / fine face 12; x = 2.0 is 11 / 22; y = 4.0 and 5.0 are coarse faces 21 and 26 from the domain origin (5 and 10 in the FDS mesh starting at y = 3.0); z = 3.0 is coarse 16 / fine 32. The walls are **normal** to the interface plane z = 3.0, so no thin face lies in that plane. In the riser cross-section, 1 m² = 25 coarse or 100 fine z-faces at z = 3.0, all gas–gas. Under R1 both levels produce the same wall planes, so the case has **zero M1-M4 mismatches** and the walls are THIN on both levels (zero input length, read.f90:11315-11317).
- **[derived]** FDS itself runs this case leak-free for the same reason: each coarse z-face at z = 3.0 lies wholly on one side of each wall, so the `IIO..KKO` averaging (velo.f90:1388-1399) never mixes the two sides.
- AMR mapping (FR-010, requirements.md:102-103; FR-016, requirements.md:125; D-009; the mapping of finer `&MESH` entries to levels is an assumption, see §3.1 C10): level 0 = 32³ at dx 0.2 over [-0.2, 6.2]³; level 1 (ratio 2) over z ∈ [3.0, 6.2]; the only C/F plane is z = 3.0. The coarse walls continue above z = 3.0 on level 0 (covered), where they coincide with the fine walls.

### 1.2 Problem statement

D-009 (README.md:44; ADR-003-geometry.md:79; requirements.md:219-224) forbids mismatches **on interface faces**. A thin wall *normal* to the interface has no interface face at all: every C/F face in duct_flow_uglmat_refine is gas–gas. The risk sits on the faces next to the wall:

- **Straddle.** If the wall's line of intersection with the C/F plane does not lie on coarse-face edges, one coarse C/F face straddles the wall. `average_down_faces` (pressure/01-amr-mapping-spec.md:244-251; D-032 (3), README.md:67) and reflux (FR-024, requirements.md:154) then give that one coarse face, and the coarse cell behind it, flux from *both* sides of the wall. Mass, species and energy cross a zero-thickness wall through the coarse cell. This is the AMR analogue of FDS's `IIO..KKO` averaging (velo.f90:1388-1399) when mesh walls disagree.
- **Partial coplanar wall.** A thin wall lying *in* the C/F plane that covers only some of the r² fine subfaces of a coarse face would be half-open on one level and closed on the other. FDS silently closes the whole face for UGLMAT (pres.f90:3726-3829).
- **Zone topology.** A leak through either route also merges pressure zones on one level but not the other (func.f90:5656-5676).
- **Static vs dynamic.** R3 ("grow … or stop", requirements.md:222) conflicts with FR-010 for static FDS multi-resolution inputs, where the hierarchy must match the input boxes (requirements.md:103).

### 1.3 Rule R3-T **[ruling]**

**(a) R1 unchanged.** Each level snaps every OBST with the unchanged FDS rule, using its own dx and the level-0 origin (D-009). Level 0 stays bitwise equal to the baseline. Each level builds its face mask from the **global** OBST list on valid+2 (ADR-001-driver-architecture.md:300), not from the FDS mesh the OBST came from. So the near-boundary collapse of read.f90:11041-11123 does not apply inside the AMR domain; see §3.3.

**(b) Ownership.** Every face belongs to the finest level whose valid region contains it (both adjacent cells). This follows from D-032 (3) (README.md:67) and pressure/01-amr-mapping-spec.md:244-251:
- A coarse face on the C/F boundary is **owned by the fine level**. Its value is always `average_down_faces` of the r² fine subfaces (the fine fluxes go through the flux register).
- A covered coarse face is also fine-owned and is overwritten the same way.
- A coarse face between two uncovered coarse cells is coarse-owned.
- A wall face (gas on one side, solid or thin wall on the other) belongs to the level whose **valid, uncovered** cell is its gas side. This is the FDS test `INTERPOLATED_MESH(IIG,JJG,KKG)>0` of divg.f90:760, generalised. For a thin wall the two sides are two wall faces (init.f90:195-300), and each is owned separately.

**(c) Interface-band consistency (the new check).** Let the **interface band** of a level pair (ℓ, ℓ+1) be all coarse cells within one coarse cell of the C/F boundary, on both sides (covered and uncovered). This contains the fine ghost region valid+2 at ratio 2, and more at ratio 4. Inside the band, **the fine face mask must be the refinement of the coarse face mask**:
- **T1 (coplanar).** A coarse face is a wall face (thin or solid) in the level-ℓ mask exactly when all of its r² fine subfaces are wall faces in the level-(ℓ+1) mask.
- **T2 (no sub-coarse walls).** No fine wall face lies on a fine face plane that is interior to a coarse cell of the band.
- Consequence for normal crossings: a thin wall crossing the C/F plane must meet it along coarse-face edges, and on both levels at the same coordinate. That is what rules out straddles.
- The rule generalises D-009's "no mismatch on interface faces" from the interface faces to the band. It is the "solid wins" condition of pres.f90:3726-3829 turned into a check that must already hold, instead of a silent fix.
- M1-M4 (requirements.md:223) stay allowed only in covered cells outside the band.

**(d) Flux register and reflux.** With (c) satisfied:
- The fine flux through every fine wall face is zero.
- Every coarse C/F face lies on one side of every wall, so `average_down_faces` and reflux move no mass across a wall.
- A thin wall lying in the C/F plane gives a zero average-down flux on the coarse face, which matches the coarse wall flag.
- A coarse C/F face that is a wall on the coarse side needs no flux-register entry (the flux is zero on both sides). The reflux correction stays well-defined and zero there.

**(e) Solid-face mask for the velocity and pressure operators.**
- **Velocity (FDS IBM, E-1, the Phase 4 baseline).** The operator is constant-coefficient (`MLPoisson`, β ≡ 1), so no β rule is needed. Each level applies the wall-face forcing (`NO_FLUX`, velo.f90:1348-1563; no H gradient at wall faces for UGLMAT/ULMAT, velo.f90:1483-1487) and the obstruction iteration on its **own owned** wall faces, using its own R1 mask. The coarse R1 mask is sufficient; no averaged fine thin-wall mask is needed on coarse faces (resolves §3.2 Q3). Covered and C/F coarse faces take the `average_down_faces` velocity (b). Under (c) that value is zero on walls and consistent elsewhere.
- **Pressure, option E-2 (β = 0 on wall faces; pressure/01-amr-mapping-spec.md:330-345, REC-E1).**
  - Every level sets β = 0 on every wall face (thin or solid) of its own R1 mask and β = 1 elsewhere, on all its faces, covered ones included.
  - Uncovered coarse faces (including those in the band outside the patch): β from the level-ℓ R1 mask. T1/T2 guarantee it matches the fine mask where they overlap.
  - Fine faces: β from the fine mask. C/F fluxes in the composite operator come from the fine level (reflux).
  - **Covered and C/F coarse faces: β = arithmetic mean of the r² fine subfaces.** The driver does not choose this value: `MLABecLaplacian::prepareForSolve` always calls `averageDownCoeffs`, which overwrites them (AMReX_MLABecLaplacian.H:497-510, 799-815; kernel `amrex_avgdown_faces`, AMReX_MultiFabUtil_3D_C.H:95-120). Under T1 the mean is exactly 0 or 1 inside the band, so it equals the R1 coarse mask there. Fractional β appears only on covered M1-M4 faces outside the band and on coarse MG levels (`averageDownCoeffsSameAmrLevel`, AMReX_MLABecLaplacian.H:711-737). It affects convergence, never the composite answer. This resolves the pressure/01-amr-mapping-spec.md §E "[VERIFY]" and §3.2 Q1.
- **C/F interpolation stencils: standard MLMG, no masked interpolant (resolves §3.2 Q2).**
  - The normal fine-side ghost at maxorder 2 uses the face value and one interior fine cell and never crosses a wall normal to the interface (AMReX_MLLinOp_K.H:124-138).
  - The tangential coarse interpolation is hard-coded to order 3 (AMReX_InterpBndryData.H:126; AMReX_MLCellLinOp.H:784, 862, 962, 977). It uses coarse `j±1`, `k±1` and corners (AMReX_InterpBndryData_3D_K.H:23-51), and its mask knows only "covered" and "outside domain" (AMReX_BndryData.H:49-51). So it **does** cross a thin wall that meets the C/F plane on a coarse-face edge, and R3-T does not prevent that. β plays no part in the ghost fill, so β = 0 does not help.
  - Accepted as is. The error is local to C/F faces next to the line where the wall meets the interface, and it scales with the H jump J across the wall (≈ J/8 + J/32 at ratio 2) **[derived]**. J is at dynamic-pressure scale, because the zone pressure lives in `PBAR`, not H (≈ 0.5 m²/s² in the duct_flow_uglmat_refine riser) **[derived]**. No mass crosses the wall, and reflux still conserves.
  - FillPatch slopes are limited (AMReX_MFInterp_3D_C.H:113-137), so contamination there is O(h·|∇φ|); the projection does not use FillPatch for H.
  - Guard: the coarse/fine flux check of §1.6 item 7. **Fallback, only if that check fails:** an AMReX patch that makes the tangential order a per-operator setting and uses order 1 (the existing `interpbndrydata_o1`, AMReX_InterpBndryData.H:211-217) only on faces whose tangential stencil crosses a wall face. Order 1 everywhere is rejected (it would lower the P2-measured C/F order of ≥ 1.92 toward 1).
- **Zone flood fill** runs per level with each level's thin faces (func.f90:5656-5676). Because of (c), the band cells of both levels get the same zone. The composite zone map (D-032 mean removal per zone, requirements.md:208) takes each uncovered cell's zone from its own level. A zone-ID mismatch between a fine cell and its covering coarse parent in the band is a setup error.

**(f) Differing snap on the two levels.** Example **[derived]**: move line 43's wall to x = 1.07. Coarse `NINT((1.07+0.2)/0.2)=NINT(6.35)=6` gives x = 1.0. Fine `NINT(12.7)=13` gives x = 1.1. The wall then crosses z = 3.0 at different planes, and T1/T2 fail in the band.
- **Static FDS multi-resolution input** (hierarchy fixed by the input, FR-010): stop at setup with an error. The error names the OBST (input line / ordinal), the two levels, both snapped planes, and the C/F plane it crosses.
  - No automatic growth: it would contradict FR-010's "hierarchy matches the input boxes", and for a normal crossing it only moves the problem to the new boundary.
  - No silent fix: the FDS silent solidification of pres.f90:3726-3829 would change the geometry without telling the user.
  - The user remeshes deliberately (as D-030 does for ratios, requirements.md:108).
- **Dynamic refinement:** the regridder tags the **whole footprint** of any OBST whose R1 snap differs between ℓ and ℓ+1 (the union of both levels' snapped extents) plus one coarse cell, or excludes it entirely from level ℓ+1 if it carries no user tag.
  - Local one-cell growth (the current R3 text) is not enough for walls normal to the interface, because the crossing just moves with the boundary.
  - If the tagged region exceeds `MAX_LEVEL` or memory limits, stop with an error (same pattern as FR-041a's negative test, requirements.md:233).
- The mismatch list is computed once at setup from the global OBST list for every level pair (it does not depend on the hierarchy), so regrids never need to rediscover it.

**(g) Regrid behaviour.** After every regrid:
- rebuild each level's wall mask from the global OBST list (stateless faces are rebuilt from input, FR-041a, requirements.md:230);
- re-run the band check (c) and the zone check (e);
- recompute the FR-034 owned-face sums (Ruling 2).

A failed check after a regrid is a bug in the tagger, not a user error, and must stop the run. The FR-039 true-residual check and the removed-mean diagnostic of §2.2 (j) run on the first solve after the regrid. Stateful thin walls (thin OBST with `MATL_ID`) follow D-010/FR-041a. See the clarification in §3.1 C6 for static hierarchies.

**(h) Validity checks run on the owning level.** ERROR(421)/(422) (init.f90:1102-1115) and `THIN` (read.f90:11315-11317) depend on the snapped geometry.
- They are evaluated on the level that owns the face.
- A covered-only M2 case (thin on the coarse level only, thick on the fine) must not raise ERROR(421) from the covered coarse copy.
- The same face owned by an uncovered coarse level must raise it, exactly as single-mesh FDS at that resolution would.

**(i) Disconnected parts. [ruling] [VERIFY]** If thin walls split a level into disconnected gas parts, the fine operator has one null vector per part while MLMG's singular handling knows one constant per level, and coarse MG levels reconnect the parts. Under E-2 each disconnected part gets one pinned cell (mask = 0), in addition to D-032's per-zone mean removal. Parts are identified by the zone flood fill of (e). Convergence rate on such a case is **[VERIFY]**; it goes on the E-2 test list (e.g. zone_shape_2 with E-2). Under E-1 the all-cell operator is connected and nothing is needed.

### 1.4 Rationale

- Keeps D-009/R1 and level-0 bit parity.
- Makes the interface conservative for walls, which FR-020/FR-024 need (reflux at round-off, requirements.md:141-143, 154-156).
- Mirrors the only FDS precedent (pres.f90:3726-3829) but as a hard check.
- Covers all three topologies: coplanar walls, normal crossings, and walls ending at the plane.
- A wall ending exactly on the C/F plane from either side satisfies (c), because no coarse face in the band straddles it.

### 1.5 Rejected alternatives

1. **Accept FDS per-mesh snapping as given (no check).** This lets normal walls straddle coarse C/F faces and leak through reflux. It is silent, and it depends on the resolution.
2. **Snap fine to coarse for OBSTs that cross the interface.** This breaks R1 on the fine level, loses fine fidelity, and makes the fine geometry depend on the hierarchy (it changes at every regrid).
3. **FDS-style silent solidification (pres.f90:3726-3829).** It changes the geometry without a message, and FDS applies it to the pressure matrix only, not to transport.
4. **ADR-003 G1 ">50%" rule / finest-extent rule.** This is still Spike G1's job (requirements.md:224). It does not solve the normal-crossing straddle by itself.
5. **Auto-grow static inputs.** This contradicts FR-010 (requirements.md:103) and D-030's "remeshed deliberately, never auto-converted" (requirements.md:108).
6. **EB/cut-cell thin walls.** GEOM/CC_IBM is deferred (D-033, README.md:68), and EB would force the C++ operator path (FR-030 constraint, requirements.md:162).
7. **Forbid thin walls crossing C/F boundaries.** Impossible for static inputs: the duct_flow_uglmat_refine riser crosses z = 3.0 by construction. Kept only as a preference in the dynamic tagger.
8. **Masked or order-1 tangential C/F interpolation now.** Needs an AMReX patch for an error that is local and O(J). Kept as the fallback of §1.3 (e), triggered by §1.6 item 7.
9. **Covered-face β from the R1 coarse mask** (the draft). Not implementable without patching AMReX: MLMG overwrites covered and C/F coarse β with the fine average before every solve.

### 1.6 Acceptance test (duct_flow_uglmat_refine)

Status: **IN** for FR-006 with the static hierarchy as given (resolves vv/scope_alignment.md:143, 165; requirements.md:226).

1. **Hierarchy.** The hierarchy dump shows level 0 = 32³ at dx 0.2 over [-0.2, 6.2]³ and level 1 (ratio 2) over z ∈ [3.0, 6.2] (FR-010).
2. **Mask checker.** Zero mismatches outside M1-M4, zero on interface faces (FR-040), and **zero R3-T band violations**. The four walls of lines 43-46 are reported at x = 1.0/2.0 and y = 4.0/5.0 on both levels, with identical zone maps in the band.
3. **Physics.** `flow_in` and `flow_out` are within 5 % of 1 m³/s at the end (dataplot row 161). T2 against the FireX UGLMAT-HYPRE baseline of the same input (FR-016 reference, requirements.md:125-136).
4. **Conservation.** The composite mass budget closes to FR-020 tolerances (per step ≤ 1e-12, cumulative ≤ 1e-10; requirements.md:143). The reflux correction on the C/F faces of the riser cross-section (25 coarse faces) and across all thin-wall faces is at round-off.
5. **Negative test** (derived copy under `vv-runs/inputs/`): line 43 moved to x = 1.07 must stop at setup with the named error of §1.3 (f). A second copy with a 0.05 m-thick wall that is THIN on level 0 only (M2 inside the band) must stop the same way.
6. **Rank and layout.** Masks and the checker report are bitwise identical at 1/2/4/8 ranks and two `max_grid_size` values (FR-005 (i)).
7. **Coarse/fine flux check at the wall/interface junctions.** Report max |J| (H jump) across the riser walls on the z = 3.0 C/F faces, and the flux error on the C/F faces adjacent to the walls (difference between the fine-face flux and the flux computed with a piecewise-constant ghost, both at the converged solution). **Trigger for the §1.3 (e) order-1 fallback:** those faces carry more than 1 % of the riser flow error, or item 3 (T2) fails.

---

## 2. Ruling FR-034: leak area across levels

### 2.1 How FDS does it (read in source)

- **Zone-pair scalar.** `&ZONE LEAK_AREA(n)`, `LEAK_PRESSURE_EXPONENT`, `LEAK_REFERENCE_PRESSURE` and `DISCHARGE_COEFFICIENT` are per zone-pair arrays (read.f90:13611-13616, 13663-13670, 13712-13717). Specifying a pair from both sides is ERROR(871); otherwise the value is copied to the partner zone (read.f90:13719-13733).
- **One duct per pair.** `SURF LEAK_PATH` pairs mark zone pairs (hvac.f90:207-218). `LEAKAGE_HVAC` creates one duct per pair with `DU%AREA_INITIAL = DU%AREA = P_ZONE(NZ2)%LEAK_AREA(NZ1)` (hvac.f90:3370-3378). Zone 0 is ambient.
- **Localized leakage.** `&HVAC TYPE_ID='LEAK'` creates a duct with its own `AREA` and two vent nodes (hvac.f90:598-713; `DU%AREA_INITIAL = AREA`, 684-685).
- **Pressure dependence.** `ADJUST_LEAKAGE_AREA` sets `DU%AREA = C_D·AREA·(ΔP/ΔP_ref)^(n−0.5)` (hvac.f90:3585-3609). **The leak flow is set by the zone-pair scalar and the node pressures, never by summed face area.**
- **Face distribution.** `HVAC_BC_IN` (hvac.f90:2299-2539) loops over all external and internal wall cells:
  - skips `NULL_BOUNDARY` (2343) and faces whose gas cell is solid (2376);
  - for LEAK_PATH faces in their own zone (2494-2504), accumulates `NODE_AREA += B1%AREA` (2511-2513; raw geometric area, not `AREA_ADJUST`) plus area-weighted `NODE_RHO/ZZ/TMP/H/X/Y/Z/P`, with `P_AVE` from `PBAR(KK,zone)` (2522-2533);
  - HVAC vents take the same path (`NODE_AREA` at 2415).
  - It does **not** test `INTERPOLATED_MESH` (`rg INTERPOLATED_MESH hvac.f90` returns nothing), unlike `USUM` (divg.f90:760).
- **Reduction.** `EXCHANGE_HVAC_BC` is a floating-point `MPI_ALLREDUCE` SUM of the node properties plus a MAX of `NODE_ZONE` (main.f90:4970-4983). These are order-dependent sums (inventory/global_reductions.csv row 87 and rows 163-201, the per-process presums in hvac.f90:2415-2533). The network is solved on rank 0 and `NODE_AREA_EX` and `DUCT_MF` are broadcast (main.f90:4988-5017). This runs at init (main.f90:562-575) and at the first pass of every step (main.f90:811-830). `HVAC_CALC` in the corrector (main.f90:980) returns early (hvac.f90:1418-1425), so `DUCT_MF` is frozen for the step. The **face velocity is recomputed every stage**: `CALC_HVAC_BC` runs from `WALL_BC` (wall.f90:178; main.f90:838 predictor, 1009 corrector) with that stage's ρ_F and `PBAR_P`.
- **Node area gate.** `COLLAPSE_HVAC_BC` sets `NODE_AREA_EX = NODE_AREA` (hvac.f90:3134). If the area is ≤ 20ε it zeroes the duct area and velocity; otherwise it restores `AREA_INITIAL` (hvac.f90:3125-3133). `NODE_AREA_EX` enters the duct equations (hvac.f90:2053, 2066, 2128, 2162).
- **Per-face flux.** `CALC_HVAC_BC`: `MFT = −DIR·DUCT_MF/NODE_AREA_EX` (wall.f90:1761). For MFT ≥ 0, `U_NORMAL = MFT/RHO_F` (wall.f90:1784-1791). For MFT < 0 the code sets `M_DOT_G_PP_ADJUST = −NODE_ZZ_EX·MFT` and the face velocity comes from the mass-flux branch of `CALCULATE_ZZ_F` (wall.f90:1446-1447). In both cases Σ_faces MFT·A = DUCT_MF **only if the flux is applied on exactly the face set that was summed into NODE_AREA**.
- **Zone sums.** `USUM = Σ U_NORMAL·B1%AREA` over SOLID_BOUNDARY wall cells whose gas cell is not covered (divg.f90:757-767, test at 760). Leak and HVAC faces are SOLID_BOUNDARY wall cells, so **their USUM contribution already carries the covered mask in FDS**. The sums are reduced at main.f90:2028-2054, giving `D_PBAR_DT = (DSUM−USUM)/PSUM` (divg.f90:1519). Leak and HVAC faces are excluded from the U_NORMAL prediction (divg.f90:1351-1352, 1415-1416). `Q_LEAK` enters the divergence at divg.f90:543-544.
- **No leak faces on thin OBSTs** (ERROR(421), init.f90:1102-1109; §0.2).

**Consequence for AMR.** Summing face-distributed quantities over all levels makes `NODE_AREA` too large by the covered area A_c. The flux applied on owned faces then delivers only (A−A_c)/A of `DUCT_MF`, so mass is not conserved against the HVAC network (FR-020, requirements.md:141), and the steady leak velocity changes by the factor A/(A−A_c). FDS has the same defect for embedded meshes (the missing `INTERPOLATED_MESH` test in `HVAC_BC_IN`); embedded meshes are FDS-only references, so this is not reproduced. The gap is only in `HVAC_BC_IN` and the `MFT` application; `USUM` is already masked. The former pressure/01-amr-mapping-spec.md:363 ("HVAC and leakage … unchanged") has been replaced by the Pressure Lead with REC-F1 (§F), which matches this ruling.

### 2.2 Rule **[ruling]**

**(a) Prescribed areas are zone-level (pair) scalars, level-independent.** `LEAK_AREA`, the exponent, the reference pressure, `C_D`, and `&HVAC TYPE_ID='LEAK' AREA`:
- are read once;
- are never summed from faces and never scaled by resolution, covered fraction or level;
- live on the network (rank 0) exactly as in FDS (hvac.f90:3370-3378, 684-685, 3585-3609).

**(b) Face-distributed quantities are summed over owned faces only.** This covers `NODE_AREA`, the area-weighted node state (`NODE_RHO/ZZ/TMP/H/X/Y/Z/P`), and the leak/HVAC contributions to `USUM`.
- **Owned face:** a wall face whose gas-side cell is a valid, **uncovered** cell of its level (`makeFineMask`), i.e. the FDS tests hvac.f90:2376 (gas side not solid) and divg.f90:760 (gas side not covered) combined.
- **Area:** `B1%AREA` = geometric face area of the owning level (init.f90:3380), not `AREA_ADJUST` (FDS uses `B1%AREA` at hvac.f90:2511 and divg.f90:765).
- **One mask for all three.** One owned-face mask serves `NODE_*`, the per-face application of `MFT` and `USUM`, and it is the same finest-owner rule as REC-C1 (pressure/01-amr-mapping-spec.md:262-263) and REC-F1 (§F). For `USUM` this is what FDS already does (divg.f90:760); the mask is new for `NODE_*` (hvac.f90:2299-2539) and for `MFT` (wall.f90:1761).

**(c) Exact, decomposition-independent summation.** All sums in (b) use the FR-005 (ii) exact fixed-point accumulation, from the per-face term to the global reduction (requirements.md:78 names "HVAC nodes" among the sums it covers; D-028/FR-005 (v) for setup-time areas, requirements.md:81). The rank-0 network solve and broadcast stay as in FDS.
- **Scale requirement for the Integration Lead:** the fixed-point scale must represent the finest-level face area exactly.
- **[derived]** When `dx_ℓ = dx_0 / r^ℓ` with r ∈ {2,4}, `fl(dx_f·dy_f) = fl(dx_c·dy_c)/r²` exactly (scaling by a power of two commutes with rounding, barring underflow). So the r² fine subfaces sum *exactly* to the coarse face area.
- With exact accumulation, `NODE_AREA` is then **bitwise identical** whatever the box layout, rank count, level structure or regrid history, provided the leak surface's edges lie on coarse faces.
- Check for the V&V Lead: dx 0.2 → 0.1 gives 0.04000000000000001 = 4 × 0.010000000000000002 exactly (checked with Python `fractions`).
- Against single-mesh FDS (sequential floating-point sums) the target is T1, within 1e-10 relative (FR-005 (v) wording).

**(d) Faces on coarse/fine boundaries.** Ownership decides every case; nothing is counted twice:
- **Leak face in the plane of the patch edge** (an OBST face or domain face coincident with the C/F plane): exactly one side is gas. The face belongs to the level of its gas cell. FDS already does this across mesh interfaces: the interface face takes the solid OBST's SURF (init.f90:3251-3262), and the solid-side copy is skipped (hvac.f90:2376).
- **Leak surface cut by the C/F boundary line** (e.g. a floor VENT partly under a patch): the uncovered coarse part plus the fine part add up exactly to the single-level area under (c).
- **Covered coarse leak face:** not owned. It contributes nothing to `NODE_*` or `USUM`. Its normal velocity is overwritten by `average_down_faces` of the fine faces (D-032 (3); pressure/01-amr-mapping-spec.md:244-251). It never gets its own `CALC_HVAC_BC` flux.
- **Leak surface whose snapped extent differs between levels** (VENT or OBST edge not on a coarse face): allowed. The total delivered flow still equals `DUCT_MF` because (b) uses one face set for the sum and for the flux. Only the distribution and the small `NODE_AREA_EX` term in the duct matrix (hvac.f90:2053, 2128) change. The setup report must print it as an M1-type note.

**(e) Per-face application.** `MFT = −DIR·DUCT_MF/NODE_AREA_EX` (wall.f90:1761) is applied at every stage on owned faces only (both MFT branches, wall.f90:1784-1791 and 1446-1447), with the exactly summed `NODE_AREA_EX`. Then Σ_owned MFT·A = DUCT_MF to round-off, which closes the FR-020 budget including HVAC (requirements.md:141). The same rule holds for ordinary HVAC vents (hvac.f90:2415) and localized `TYPE_ID='LEAK'` vents.

**(f) Regrid.** FDS recomputes the node sums and solves the network once per step (main.f90:811-830; corrector early return, hvac.f90:1418-1425), and the face velocities every stage. AMR does the same with the current owned mask, so a regrid needs no special transfer. Under (c), `NODE_AREA` for coarse-aligned leak surfaces is bitwise unchanged across regrids, so the network sees no jump. If ownership leaves a node with zero area, the FDS ≤ 20ε gate zeroes the duct (hvac.f90:3125-3133). That can only happen through a mapping error, and it must be reported.

**(g) Zone identity.**
- `B1%PRESSURE_ZONE` comes from the gas-side cell (init.f90:1087) on its own level.
- Explicit zones (`&ZONE XYZ`) are seeded per level from the same point (main.f90:2630-2641).
- **Owner: AMReX Integration Lead** (setup/driver code and FR-005 layout independence; the Pressure Lead is a consumer). Zone IDs must be identical across layouts, ranks and regrids, and identical between a fine cell and its covering parent in the band.
- **Auto-numbered zones** (main.f90:2653-2682, numbered in mesh scan order; not re-read for the final) must get a deterministic, hierarchy-independent numbering (e.g. from the level-0 scan). LEAK_PATH and LEAK_AREA refer to zones by integer index, so a numbering that depends on the box layout would silently re-wire the leak network.
- A `&ZONE` without `XYZ` is seeded from "the first non-solid cell in the first mesh controlled by this MPI process" (read.f90:13695-13707). That point depends on the rank. AMR must seed it from level 0 deterministically.

**(h) Background pressure: one composite definition.** FDS holds one global `P_0(z)` (init.f90:452-484), initialises the per-mesh `PBAR` from it (init.f90:1044-1054), and shifts it uniformly by the zone scalar each stage (mass.f90:548, 730), with one globally reduced `D_PBAR_DT` per zone (divg.f90:1519).
- AMR stores **one scalar ΔP_zone(t) per zone**, advanced once per stage from the composite `D_PBAR_DT`.
- `P_0(z)` stays a function (the ramp), built once from global extents, not a level-0 array. Ramp spacing `RP%RDT` against the finest dz is **[VERIFY]** (init.f90:458-469).
- Each level (or box) evaluates `PBAR_ℓ(k,zone) = P_0(z_k^ℓ) + ΔP_zone`, with its own `R_PBAR` and `PBAR_S`, at its own cell centres, and at faces where FDS averages (`P_AVE`, hvac.f90:2522-2530). `NODE_P` takes the owning level's sample.
- Rejected: integrating `D_PBAR_DT` per level (breaks the single per-zone compatibility of divg.f90:1540); one level-0 z-array interpolated to fine levels (adds interpolation error in `R_PBAR`, `RHO_F`, `NODE_P`).

**(i) Order within each stage and the removed-mean diagnostic.**
1. owned-face `CALC_HVAC_BC`;
2. exact per-zone DSUM/PSUM/USUM over owned faces and uncovered cells;
3. `D_PBAR_DT = (DSUM−USUM)/PSUM` (divg.f90:1519);
4. the D correction (divg.f90:1540);
5. RHS assembly;
6. D-032 composite per-zone mean removal;
7. solve.

With the same stage's face velocities in USUM and in the Poisson BC, the mean removal strips only round-off. The frozen `DUCT_MF` is a first-order time lag that FDS also has; it does not break the balance.

**(j) Diagnostic (addition to FR-039, not a change to D-032).** Mean removal would silently hide a face-set mismatch (e.g. USUM on all levels, MFT on owned faces: the per-zone mean becomes ΔU·A_covered/V_zone). So at the FR-039 true-residual trigger points (first solve of the run, first solve after every regrid, every solve in debug) the driver logs the removed mean per zone relative to ‖b‖, and prints a WARNING when it exceeds 1e-10·‖b‖. Expected level with exact sums: about 1e-14.

**(k) Leak nodes.** One leak node per zone pair, with owned faces on any number of levels. FDS already builds each node from faces on many meshes with one global reduction (main.f90:828) and one network solve (main.f90:829); levels are just more face owners. No per-level node, no flow partitioning, no network change. ERROR(552) ("ductnode must lie within a single pressure zone", hvac.f90:2410) is checked on owned faces with each level's zone map; R3-T (§1.3 (e)) guarantees a fine cell and its band parent share a zone.

**(l) Phase placement.** Multi-zone with leakage gates **Phase 4**, with the refined variants of zone_shape_2 and HVAC_leak_exponent. The uniform variants gate Phase 2 through FR-006. Nothing in the pressure design needs Phase 6 for multi-zone: there is no subcycling, and the composite per-zone sums, `D_PBAR_DT` and mean removal work for any number of zones in one solve. Phase 6's dP0/dt sync exists only for subcycling.

### 2.3 Rationale

- FDS already separates the two quantities: the scalar decides *how much* flows (hvac.f90:3370-3378, 3585-3609), and the face sum decides only *where* (wall.f90:1761). Counting area twice therefore breaks conservation, not the leak law.
- One owned-face mask for `NODE_*`, `MFT` and `USUM` is the only choice that makes Σ MFT·A = DUCT_MF and the `USUM` term consistent with each other and with REC-C1.
- The exact sum plus the power-of-two argument makes the leak areas invariant to the hierarchy. This covers the setup-time area half of D-028 and the runtime half of FR-005 (ii).

### 2.4 Rejected alternatives

1. **Sum over all levels** (FDS's embedded-mesh behaviour, hvac.f90:2299-2539 without an `INTERPOLATED_MESH` test): double counts, loses mass against the network.
2. **Per-level leak area or per-level networks:** there is one physical leak per zone pair (hvac.f90:3370-3378). Splitting it by level would need arbitrary flow partitioning.
3. **Scale `LEAK_AREA` by the uncovered fraction:** changes the prescribed physics with the hierarchy.
4. **Floating-point `MPI_Allreduce` of node sums** (FDS main.f90:4979): order-dependent, which violates FR-005 (ii) (requirements.md:78).
5. **Use `AREA_ADJUST`:** FDS uses the raw `B1%AREA` for leak/HVAC nodes (hvac.f90:2511) and for `USUM` (divg.f90:765). Changing it would break T1 parity.

### 2.5 Acceptance tests

**zone_shape_2** (`Pressure_Effects/zone_shape_2.fds`; 8 meshes of 16³, dx 0.2, lines 3-4; `-p 8`, FDS_Cases.sh:431)
- **Input facts.** Zone 1 at (0,0,0), zone 2 at (1.5,1.5,1.5) with `LEAK_AREA(1)=0.007746` (lines 6-7). `SURF 'DUCT LEAK' LEAK_PATH=2,1` (line 16). Fan `VOLUME_FLOW=0.1` into the duct, which is zone 2 (lines 18-23), enclosed by thin duct walls (lines 25-61). The leak faces are:
  - the +x face of `OBST XB=1.0,1.2,1.0,2.0,4.0,5.0` (line 50, zone 2 side, 1.0 m²);
  - the floor VENT `XB=4.0,5.0,4.0,5.0,-0.2,-0.2` (line 51, zone 1, 1.0 m²).

  `VELOCITY_TOLERANCE=0.001` (line 13). `T_END=300` (line 10).
- **[derived]** Steady ΔP = ρ/2·(0.1/0.007746)² ≈ 100 Pa, which matches the expected value.
- **Variants.**
  - **U:** the input as given (a single uniform level in AMR mode).
  - **R1-leak (A-41 derived copy):** a level-1 patch (ratio 2) that covers the leak OBST, e.g. x ∈ [0.6, 2.2], y ∈ [0.6, 2.2], z ∈ [3.8, 5.4] (edges on 0.8 m blocking-factor lines). The thin walls at y = 1.0/2.0 and z = 4.0/5.0 cross its x = 2.2 face normally, on coarse faces, so R3-T passes.
  - **R1-split:** a patch x ∈ [4.6, 5.4], y ∈ [3.8, 5.4], z ∈ [-0.2, 0.6] that splits the floor VENT at x = 4.6 into 0.6 m² coarse + 0.4 m² fine.
  - The V&V Lead finalises the variants under A-41 (vv/scope_alignment.md:118).
- **Pass criteria.**
  1. The setup report prints `NODE_AREA_EX` of both leak nodes = 1.0 m² (the float sum of the owned faces). The values are **bitwise identical** across U, R1-leak and R1-split, 1/2/4/8 ranks and two layouts (T0), and within 1e-10 relative of single-mesh FDS (T1).
  2. The zone-2 flood fill is identical on both levels, with zero R3-T violations (Ruling 1).
  3. `Delta p` = 100 Pa as the mean over 200-300 s, within 1 %; total mass 313 kg at the end, within 1 % (dataplot rows 802-803).
  4. FR-020 budget including HVAC at round-off: Σ_owned MFT·A = DUCT_MF per step to ≤ 1e-12 relative.
  5. Negative check (debug option): summing all levels must fail criterion 4 in R1-leak.
  6. The removed-mean diagnostic of §2.2 (j) stays below 1e-10·‖b‖ at every trigger point, and in the negative check of criterion 5 it warns.

**HVAC_leak_exponent** (`HVAC/HVAC_leak_exponent.fds`; serial, FDS_Cases.sh:381)
- **Input facts.** Three disjoint 20³ meshes of dx 0.5 at x ∈ [0,10], [20,30], [40,50] (lines 9-11). Zones with `LEAK_AREA(0)=0.01`, exponents 0.5/0.6/0.6 and reference pressures 4/4/10 Pa (lines 13-15). `LEAK_PATH=n,0` floor and ceiling VENTs of 100 m² each (SURFs lines 21-23, VENTs lines 25-35), so 200 m² per zone node. `INFLOW VEL=-0.01` on a 4×4 m VENT gives 0.16 m³/s. `L1-L3` are `DUCT VELOCITY` of `LEAK 0 n` (lines 37-39). Expected −16, −11.79745, −12.7282 m/s, relative 0.01 at the end (dataplot row 378).
- **Precondition.** AMR mode needs a box level 0 over [0,50]×[0,10]×[0,10] with the gaps x ∈ [10,20] and [30,40] filled with solid. The requirements now record this non-box rule as *confirmed by the AMReX Integration Lead; final ruling: AMR Chief Architect* (requirements.md:307). Older docs still call it unconfirmed (vv/scope_alignment.md:84, 183; vv/test-plan.md:334). Also, the `INFLOW` VENTs at x = 20 and x = 40 (lines 31, 35) sit on FDS mesh faces inside the bounding box, so they become fluid–solid boundary conditions on the fill (requirements.md:307, the A-46 flag). This is not a leak face and does not change the ruling.
- **Variants.**
  - **U:** the input mapped to one uniform level.
  - **R (A-41 derived copy):** level 1 (ratio 2, dx 0.25) over zone 1's floor strip x ∈ [0, 4], z ∈ [0, 2] (40 m² fine + 60 m² coarse of the floor), and over zone 2's whole ceiling x ∈ [20, 30], z ∈ [8, 10].
  - dx 0.5 → 0.25 is dyadic, so every face area is exact.
- **Pass criteria.**
  1. `NODE_AREA_EX` = 200 m² exactly for each of the three zone leak nodes, bitwise identical across U/R, ranks and layouts (T0), and equal to single-mesh FDS (T1 ≤ 1e-10).
  2. L1/L2/L3 = −16, −11.79745, −12.7282 m/s within 1 % at the end.
  3. FR-020 budget including HVAC at round-off per zone.
  4. **[derived] sensitivity:** double counting in variant R would make zone 1's node area 240 m², deliver 200/240 of the duct flow, and drive |L1| toward 16 × 1.2 = 19.2 m/s, a 20 % failure. So criterion 2 detects the defect this ruling prevents.

### 2.6 Resolved questions (Pressure Solver Lead, pressure/03-fr034-fr040r3-answers.md)

- **P1 (USUM mask).** Yes, one owned-face mask; FDS already masks USUM (divg.f90:760). The spec line 363 is replaced by REC-F1. See §2.2 (b).
- **P2 (PBAR per level).** Per-level arrays sampled from one composite definition. See §2.2 (h).
- **P3 (phase).** Phase 4. See §2.2 (l).
- **P4 (per-step coupling vs per-stage D_PBAR_DT).** No clash with D-032; stage order and diagnostic in §2.2 (i)-(j).
- **P5 (zones spanning levels).** One node per zone pair. See §2.2 (k).
- **P6 (zone numbering and seeding).** AMReX Integration Lead. See §2.2 (g).

---

## 3. Cross-cutting items

### 3.1 Facts that contradict or outdate the existing ADRs and docs

- **C1.** ADR-003-geometry.md:100 and :115 still list thin walls under an open Q4, and :94 says "HT3D and thin walls … FR-004 rejection in AMR mode until assessed". D-033 (README.md:68; spec-responses.md:170) made thin OBSTs **required** with refinement. ADR-003 needs a v0.3 edit.
- **C2.** ADR-003-geometry.md:27 cites `init.f90:118-122` as the handling of zero-thickness OBSTs. Those lines are the HT3D-only `THIN_WALL` path (created only if `HT_DIM>1`, init.f90:118-191). General thin faces are ordinary wall cells with `WC%THIN` (init.f90:3315), on two cells per face (init.f90:195-300).
- **C3.** FR-040 R3 ("grow … or stop", requirements.md:222) conflicts with FR-010 ("hierarchy matches the input boxes", requirements.md:103) for static FDS multi-resolution inputs. This ruling resolves it: stop for static inputs, whole-footprint tagging for dynamic runs.
- **C4.** D-009's guarantee "no mismatch on interface faces" (README.md:44) does not cover thin walls **normal** to the interface. There the interface faces are all gas, but a straddle still leaks through reflux. The band rule R3-T is needed in addition.
- **C5.** vv/test-plan.md:382 lists duct_flow_uglmat_refine as "supported 2:1, no action", while vv/scope_alignment.md:143/165 and requirements.md:226 carry it as UNCLEAR pending this ruling. [derived] The case in fact has zero mismatches and passes R3-T, so it becomes IN with a real check behind it.
- **C6.** D-010/FR-041a ("the coarse-fine boundary never cuts a stateful OBST face", requirements.md:231) is worded absolutely. A static multi-resolution input with a stateful thin wall crossing the interface would be rejected, although no regrid (and so no state transfer) ever happens. This needs clarification. duct_flow's walls are stateless (line 19), so the test case is not affected.
- **C7. (Resolved.)** pressure/01-amr-mapping-spec.md:363 said HVAC and leakage were "unchanged". The Pressure Lead replaced it with REC-F1 (§F): the node sums and the `MFT` application need the owned-face mask, because FDS `HVAC_BC_IN` lacks the covered-cell test that `USUM` already has (hvac.f90:2299-2539 vs divg.f90:760).
- **C11.** requirements.md:179 ("TBD(Pressure Lead) whether multi-zone moves to Phase 6") is resolved: Phase 4 (§2.2 (l)).
- **C8.** The brief describes the FDS multi-mesh-to-level mapping as "D-006/D-007/D-009". In the decision log, D-006 is the M2a demo gate and D-007 is the Phase 6 scope rule (README.md:41-42). The mapping itself is FR-010 (requirements.md:102-108) and IR-002 (`&MESH` semantics, requirements.md:306) plus D-009 (README.md:44). The brief's references should be corrected.
- **C9.** vv/scope_alignment.md:84/183 and vv/test-plan.md:334 still treat the solid-filled non-box domain as unconfirmed. requirements.md:307 now records it as confirmed by the Integration Lead, with the final ruling left to the Chief Architect. HVAC_leak_exponent's IN status depends on that ruling.
- **C10.** IR-002 (requirements.md:306) says `&MESH` entries define the level-0 domain and that "all level-0 meshes must share one uniform cell size". It does not say what finer `&MESH` entries become. FR-016 (requirements.md:125) and vv/test-plan.md:382 assume that an FDS multi-resolution input such as duct_flow_uglmat_refine maps to a static two-level hierarchy. Ruling 1 rules, **[ruling]**: the level-0 domain is the bounding box of all `&MESH` entries at the coarsest cell size; each finer `&MESH` becomes a static box on level log_r(dx_0/dx) (FR-010). IR-002 should say this explicitly. Read literally, IR-002 would reject the input (0.2 m and 0.1 m meshes).

### 3.2 Former open questions

- **Q1 (β on covered and band faces). Resolved:** fine-face average set by MLMG; exactly 0 or 1 in the band under T1. See §1.3 (e).
- **Q2 (C/F interpolation across a thin wall). Resolved:** standard MLMG plus the §1.6 item 7 check; order-1 AMReX patch only as the fallback. See §1.3 (e).
- **Q3 (E-1 coarse mask). Resolved:** the coarse R1 mask suffices; forcing on each level's owned wall faces. See §1.3 (e).
- **Q4 (OBSTs created or removed at run time). [ruling]** The R3-T mismatch list is computed at setup for **every** OBST in the input, whatever its initial activation state, so activating one at run time cannot create an unchecked band mismatch. **[VERIFY]** (FDS Legacy Mapper, before Phase 5): that DEVC/CTRL creation and removal only toggle OBSTs present in the input list and never move them.
- **Still open:** the disconnected-parts pin of §1.3 (i) and the `P_0` ramp spacing of §2.2 (h), both **[VERIFY]**.

### 3.3 Observation for the FDS Legacy Mapper (reading only, not run)

The near-boundary snap in READ_OBST takes the neighbour's ghost width only `IF (NOM>0 .AND. PROCESS(NM)==MY_RANK)` (read.f90:11044-11048 and the analogous branches through 11123). Otherwise it uses the mesh's own `DX(0)`/`DX(IBP1)`. As written, the test is on the rank that owns NM, not on whether NOM's data is available. Whether that can change the snap with rank count or resolution is not established. It is irrelevant to the AMR level masks, which are built from the global OBST list (§1.3 (a)), but it could affect FDS baselines of multi-resolution inputs such as duct_flow_uglmat_refine. Please confirm or dismiss it (candidate for `inventory/uniform_grid_assumptions.csv`).

---

## 5. Changes from the draft

- **§1.3 (e), E-2 β (Pressure Lead correction a).** Covered and C/F coarse faces do not keep β from the R1 mask; MLMG overwrites them with the fine average before every solve. Rewritten; draft wording added to §1.5 as rejected alternative 9.
- **§2.1, U_NORMAL (correction b).** `U_NORMAL = MFT/RHO_F` holds only for MFT ≥ 0; MFT < 0 goes through wall.f90:1446-1447. Conclusion unchanged.
- **§2.1, §2.2 (f), timing (correction c).** Node sums and `DUCT_MF` update once per step, but the HVAC face velocity is recomputed every stage.
- **§2.1, §2.2 (b), C7, USUM (correction d).** FDS already masks USUM against covered cells, leak and HVAC faces included. Only `HVAC_BC_IN` and the `MFT` application lack the mask.
- **New:** §1.3 (i) disconnected-parts pin [VERIFY]; §1.3 (e) interpolation ruling and fallback; §1.6 item 7; §1.5 items 7-9; §2.2 (h) composite PBAR; §2.2 (i)-(j) stage order and removed-mean diagnostic; §2.2 (k) leak nodes and ERROR(552); §2.2 (l) Phase 4; §2.5 zone_shape_2 criterion 6; §3.2 Q4 ruling.
- **Moved:** the FR-034 wording (former §2.2 (h)) is now in §6, extended.

## 6. Requirements text for the Spec Lead (requirements.md)

Same batch as the IR-002 and FR-041a fixes. Line numbers are requirements.md as of that date.

**FR-034, line 179 (phase), replace the last two sentences with:** "Phase 4, including multi-zone with leakage (refined variants of `zone_shape_2` and `HVAC_leak_exponent`; uniform variants gate Phase 2 through FR-006). Phase 6 adds only the dP0/dt sync correction needed by subcycling (ADR-002). (Pressure Lead, resolved 2026-09-25.)"

**FR-034, line 180 (open point), replace with a rule (Chief Architect ruling, `adr/drafts/rulings-FR040R3-FR034.md` §2):**
- "Prescribed leakage areas (`&ZONE LEAK_AREA` and its exponent, reference pressure and discharge coefficient; `&HVAC TYPE_ID='LEAK' AREA`) are zone-pair scalars, read once, independent of the hierarchy, and never summed from faces."
- "Face-distributed leak/HVAC quantities (`NODE_AREA`, area-weighted node state, per-face mass flux `MFT`, the leak/HVAC part of `USUM`) use one owned-face mask: faces whose gas-side cell is valid, uncovered and not solid, with `B1%AREA` of the owning level and the FR-005 (ii) exact sum. Covered coarse leak/HVAC faces get no flux of their own and take `average_down_faces`. Check: Σ_owned MFT·A = DUCT_MF per step to ≤ 1e-12 relative (FR-020)."
- "One leak node per zone pair, with owned faces on any level; ERROR(552) is evaluated on owned faces."
- "Background pressure: one ΔP per zone, advanced each stage from the composite `D_PBAR_DT`; each level evaluates `P_0(z)` at its own cell centres."
- "Acceptance: `zone_shape_2` (`NODE_AREA_EX` = 1.0 m² bitwise across variants, ranks and layouts; Δp 100 Pa and mass 313 kg within 1 %) and `HVAC_leak_exponent` (200 m² per node; L1/L2/L3 = −16, −11.79745, −12.7282 m/s within 1 %). Both move from UNCLEAR to IN for FR-006 (`HVAC_leak_exponent` also depends on the non-box level-0 ruling, IR-002)."
- "Zone numbering and seeding: deterministic and hierarchy-independent. Owner: AMReX Integration Lead."

**FR-039, add a bullet after the true-residual check:** "Removed-mean diagnostic: at the same trigger points the driver logs the RHS mean removed per zone, relative to ‖b‖, and prints a WARNING when it exceeds 1e-10·‖b‖ (expected ≈ 1e-14 with exact sums). It detects a face-set mismatch between `USUM` and the leak/HVAC flux that mean removal would otherwise hide."

**FR-040, line 222 (R3), replace with:** "R3 (R3-T, D-009 amendment): no mismatch on coarse-fine interface faces, and inside the interface band (coarse cells within one coarse cell of the C/F boundary, covered and uncovered) the fine face mask is the exact refinement of the coarse face mask: a coarse face is a wall exactly when all its r² fine subfaces are walls (T1), and no fine wall face lies on a fine face plane interior to a coarse cell of the band (T2). Static FDS multi-resolution inputs that violate R3-T stop at setup with an error naming the OBST, both levels, both snapped planes and the C/F plane. Dynamic runs tag the whole footprint of every OBST whose R1 snap differs between levels plus one coarse cell, or keep it off the finer level; if that exceeds `MAX_LEVEL` or memory, stop. Masks, the band check and the zone check are rebuilt after every regrid; a failure after a regrid stops the run. Validity checks (ERROR(421)/(422), `THIN`) are evaluated on the owning level. MLMG keeps its standard C/F interpolation."

**FR-040, line 226 (open point R3), replace with:** "`Pressure_Solver/duct_flow_uglmat_refine` is IN for FR-006 with the static two-level hierarchy as given; it has zero M1-M4 mismatches and passes R3-T. Acceptance (ruling §1.6): hierarchy dump, zero band violations, `flow_in`/`flow_out` within 5 % of 1 m³/s and T2 vs the FireX UGLMAT baseline, composite mass budget at FR-020 tolerances, negative tests (wall at x = 1.07; 0.05 m wall THIN on level 0 only), bitwise masks across ranks and layouts, and a coarse/fine flux check at the wall/interface junctions whose failure (> 1 % of the riser flow error, or T2 fails) triggers the order-1 tangential interpolation fallback." Also: vv/test-plan.md:382 and vv/scope_alignment.md:143, 165 change from "no action"/UNCLEAR to IN with this acceptance (V&V Lead).

**FR-041a, line 231, append to the Rule:** "Clarification (D-010): the rule applies to hierarchies that change during the run. In a static hierarchy fixed by the input, a stateful wall may cross a coarse-fine interface that never moves, because no state transfer ever occurs; R3-T still applies."

**IR-002, line 308, append:** "Finer `&MESH` entries: the level-0 domain is the bounding box of all `&MESH` entries at the coarsest cell size, and each finer `&MESH` becomes a static box on level log_r(dx_0/dx), r ∈ {2, 4}, with the hierarchy fixed by the input (FR-010). A finer mesh whose cell size is not dx_0/r^ℓ, or whose edges do not lie on its parent level's faces, stops at setup (D-030)." The rule "all level-0 meshes must share one uniform cell size" then reads "all `&MESH` entries that map to the same level".

---

## 4. Citation index (FireX 36975d765f unless a doc)

| Topic | Location |
|---|---|
| Per-mesh OBST snap, THICKEN, collapse, discard, THIN | read.f90:10991, 11127-11139, 11157-11162, 11167-11173, 11174-11216, 11220-11224, 11315-11317 |
| Near-boundary collapse, rank test | read.f90:11041-11123 (11044-11048) |
| Finer-mesh OBST invisible only | read.f90:11461-11465 |
| Zero-thickness blocks no cells | read.f90:11708; func.f90:5503-5526 |
| OBST wall faces (two per thin face); exposed faces | init.f90:195-300 (199, 204-211, 221) |
| HT3D-only THIN_WALL | init.f90:118-191 |
| INIT_WALL_CELL C/F search, ERROR(431), solid wins at interface, WC%THIN, area | init.f90:3145-3180, 3184-3226, 3251-3262, 3315, 3317-3326, 3368-3375, 3380 |
| ERROR(421)/(422) thin OBST with LEAK/HVAC | init.f90:1080-1115 |
| Covered-cell mask | init.f90:1564-1587 |
| Ghost H averaging over IIO..KKO | velo.f90:1388-1399 |
| NO_FLUX wall forcing; DHFCT=0 for UGLMAT/ULMAT | velo.f90:1348-1563, 1483-1487 |
| UGLMAT C/F solid-wins stages | pres.f90:3662-3902 (3726-3766, 3780-3800, 3801-3829) |
| Zone flood fill blocked by thin OBST | func.f90:5558-5729 (5656-5676) |
| Zone setup, explicit/auto zones | main.f90:2606-2682 |
| ZONE LEAK_AREA read, symmetry, ERROR(871); no-XYZ seed | read.f90:13611-13739 (13695-13707, 13719-13733) |
| Leak ducts from LEAK_PATH; LEAKAGE_HVAC area | hvac.f90:207-218, 3370-3378 |
| Localized leakage AREA | hvac.f90:598-713 |
| HVAC_BC_IN face sums (no covered test) | hvac.f90:2299-2539 (2343, 2376, 2415, 2494-2533) |
| COLLAPSE_HVAC_BC, NODE_AREA_EX | hvac.f90:3101-3221 (3125-3134) |
| ADJUST_LEAKAGE_AREA | hvac.f90:3585-3609 |
| NODE_AREA_EX in duct equations | hvac.f90:2053, 2066, 2128, 2162 |
| Per-face HVAC flux | wall.f90:1745-1797 (1761, 1784-1791) |
| HVAC exchange (float Allreduce), broadcast, drivers | main.f90:4970-5017, 562-575, 811-830 |
| USUM, D_PBAR_DT, DSUM/PSUM mask | divg.f90:727-767, 1495-1507, 1519, 1540 |
| Test inputs | duct_flow_uglmat_refine.fds:5-9, 11, 14, 16, 17, 19, 21-26, 43-46, 84-85; zone_shape_2.fds:3-4, 6-7, 10, 13, 16, 18-23, 25-61; HVAC_leak_exponent.fds:9-15, 21-23, 25-35, 37-39 |
| HVAC corrector early return; per-stage WALL_BC | hvac.f90:1418-1425; wall.f90:178; main.f90:828-830, 838, 980, 1009 |
| MFT < 0 branch; ERROR(552) | wall.f90:1446-1447; hvac.f90:2410 |
| P_0, PBAR init and update | init.f90:452-484, 1044-1054; mass.f90:548, 730 |
| MLABecLaplacian β average-down (AMReX 99ddfda) | AMReX_MLABecLaplacian.H:497-510, 711-737, 799-815; AMReX_MultiFabUtil_3D_C.H:95-120 |
| C/F interpolation, order 3 and mask (AMReX 99ddfda) | AMReX_InterpBndryData.H:126, 186-260 (o1 at 211-217); AMReX_InterpBndryData_3D_K.H:23-51; AMReX_BndryData.H:49-51; AMReX_MLCellLinOp.H:784, 862, 962, 977 |
| Normal ghost; FillPatch limited slopes (AMReX 99ddfda) | AMReX_MLLinOp_K.H:124-138; AMReX_MFInterp_3D_C.H:113-137 |
| Dataplot expectations | Utilities/Python/FDS_verification_dataplot_inputs.csv:161, 378, 802-803 |
| Docs | README.md:41-44, 63, 67-68; requirements.md:16, 78, 81, 102-108, 306-307, 141-143, 154-156, 160-162, 171, 178-180, 125-136, 208, 219-233; ADR-001-driver-architecture.md:300; ADR-003-geometry.md:27, 79-82, 94, 100, 115; pressure/01-amr-mapping-spec.md:244-264, 330-345, 363; vv/scope_alignment.md:84, 108, 118-119, 132-133, 142-143, 164-165, 183; vv/test-plan.md:334, 382, 410 |
