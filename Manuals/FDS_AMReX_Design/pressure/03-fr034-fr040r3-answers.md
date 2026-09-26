# 03 — Answers to the AMR Chief Architect: FR-034 (leak area across levels) and FR-040 R3 (thin walls at C/F)

| Field | Value |
|---|---|
| Author | AMR Pressure Solver Lead |
| Date | 2026-09-25 |
| Answers | `docs/adr/drafts/rulings-FR040R3-FR034.md` (DRAFT): §2.6 P1–P5 (FR-034) and §3.2 Q1–Q2 (FR-040 R3, E-2 items) |
| Source pins | FireX `36975d765f` (`Source`, read-only). AMReX `99ddfda` (`(local AMReX checkout)`, read-only). Docs cited relative to `docs/`. |
| Status | Proposal from the Pressure Lead. Nothing here changes D-009, D-028, D-032 or a requirement. The Chief Architect decides. |

Every source claim below was read at the cited lines. Items not read are marked **[VERIFY]**. Nothing was built or run for this note
(reading and `rg` only; the exclusive timing window was respected).

## Assumptions
* **A1.** The decisions in force are those in `docs/README.md`: D-009 (README.md:44), D-021 (:56), D-028 (:63), D-032 (:67, ADR-002 v0.2.1), D-033 (:68).
  One global dt, no subcycling, composite MLMG at predictor and corrector. eps_H = max(1e-8, 2.4e-12·N²) (requirements.md:18-19).
* **A2.** Phase 4 baseline for solids is **E-1** (constant-coefficient operator, IBM forcing on wall faces, FDS FFT/GLMAT default).
  **E-2** (β = 0 on wall faces, `MLABecLaplacian`, overset mask) follows (01 §E, REC-E1). Answers are given for both where they differ.
* **A3.** The Architect's R3-T band rule (draft §1.3 (c)) is adopted as written: inside the interface band, the fine face mask is the
  refinement of the coarse face mask (T1, T2). Several answers below depend on it.
* **A4.** The "owned face" definition is the draft's §2.2 (b): a wall face whose gas-side cell is a valid, **uncovered** cell of its level
  (the FDS tests `hvac.f90:2376` and `divg.f90:760` combined).

---

## 0. Verification of the draft's FDS claims

| Draft claim | Verdict | Evidence |
|---|---|---|
| `LEAKAGE_HVAC` sets `DU%AREA_INITIAL = DU%AREA = P_ZONE(NZ2)%LEAK_AREA(NZ1)` (hvac.f90:3370-3378) | **Correct** | assignment at hvac.f90:3374-3375, inside 3366-3380 |
| `ADJUST_LEAKAGE_AREA` sets area from C_D and (ΔP/ΔP_ref)^(n−0.5) (hvac.f90:3585-3609) | **Correct** | hvac.f90:3585-3609; ΔP from `DUCTNODE%P` of the two nodes |
| Localized leak `AREA_INITIAL = AREA` (hvac.f90:684-685) | **Correct** | hvac.f90:684-685 |
| `HVAC_BC_IN` (hvac.f90:2299-2539) has no covered-cell test | **Correct** | skips `NULL_BOUNDARY` (2343) and solid gas cells (2376) only; `rg INTERPOLATED_MESH hvac.f90` returns nothing. Vent `NODE_AREA` at 2415, leak branch from 2494, `NODE_AREA` at 2513, `P_AVE` at 2522-2533 |
| So summing every level would inflate `NODE_AREA` | **Correct** | follows from the row above; FDS has the same defect for embedded meshes |
| `USUM` skips covered gas cells (divg.f90:760) | **Correct, and it matters more than the draft says** | `WALL_LOOP4` (divg.f90:757-767) skips `INTERPOLATED_MESH(IIG,JJG,KKG)>0` (760), `IPZ<1` and every non-`SOLID_BOUNDARY` wall cell. Leak and HVAC faces are `SOLID_BOUNDARY` wall cells, so **their USUM contribution already carries the covered mask in FDS**. Only `NODE_*` and the per-face `MFT` application lack it. |
| `CALC_HVAC_BC`: `MFT = −DIR·DUCT_MF/NODE_AREA_EX` (wall.f90:1761), `U_NORMAL = MFT/RHO_F` (wall.f90:1784-1791) | **Partly wrong** | 1761 correct. `U_NORMAL(_S) = MFT/RHO_F` with `RHO_F = PBAR_P(KK,zone)/(RSUM_F·TMP_G)` holds only for **MFT ≥ 0** (1784-1791). For MFT < 0 the code sets `M_DOT_G_PP_ADJUST = −NODE_ZZ_EX·MFT` and the face velocity comes from the mass-flux branch of `CALCULATE_ZZ_F` (wall.f90:1446-1447). Σ MFT·A = DUCT_MF still holds on the face set; the conclusion is unaffected. |
| "The ΔP-driven leak coupling updates at FIRST_PASS only (main.f90:811-830)" | **Incomplete** | Node sums (`HVAC_BC_IN`, main.f90:822, `IF (FIRST_PASS .AND. HVAC_SOLVE)`) and the network solve (`HVAC_CALC`, main.f90:829) are once per step; `HVAC_CALC` in the corrector (main.f90:980) returns early (hvac.f90:1418-1425), so `DUCT_MF` is frozen for the step. But `CALC_HVAC_BC` runs **every stage** from `WALL_BC` (wall.f90:178; main.f90:838 predictor, 1009 corrector) and recomputes the face velocity with that stage's ρ_F and `PBAR_P`. |
| Leak/HVAC faces are excluded from the U_NORMAL prediction (divg.f90:1351-1352, 1415-1416) | **Correct** | exit on `NODE_INDEX>0` or `LEAK_PATH>0` at 1351-1352; no copy of `U_NORMAL_S` in the corrector at 1415-1416 |
| `Q_LEAK` enters D (divg.f90:543-544) | **Correct** | set at wall.f90:1760/1773 |
| Ghost H across a mesh interface is the `IIO..KKO` average (velo.f90:1388-1399) | **Correct** | piecewise-constant average of the covering cells on the other mesh (1390-1399). FDS therefore uses **order-1** coarse data at multi-resolution interfaces; MLMG uses order 3 (Q6). |
| UGLMAT/ULMAT: no H gradient at wall faces (velo.f90:1483-1487) | **Correct** | `DHFCT=0` for UGLMAT/ULMAT, and for GLMAT on external wall cells (1483-1487) |
| E-2 covered coarse faces take β "from the R1 coarse mask" (draft §1.3 (e)) | **Does not survive MLMG setup** | `MLABecLaplacian::prepareForSolve` always calls `averageDownCoeffs` (AMReX_MLABecLaplacian.H:497-510), which overwrites the coarse AMR level's face β under the fine level with the arithmetic mean of the r² fine subfaces (`averageDownCoeffsToCoarseAmrLevel`, :799-815; kernel `amrex_avgdown_faces`, AMReX_MultiFabUtil_3D_C.H:95-120, factor `1/(facy·facz)` at :109). Whatever the driver puts on covered faces is replaced. See Q5. |

Not re-read for this note (cited by the draft, not needed for the answers): read.f90 snap lines, init.f90:195-300/3315 thin-wall
wall cells, func.f90:5656-5676 zone flood fill, main.f90:2653-2682 auto-zone numbering **[VERIFY, not re-read by me]**.

---

## 1. FR-034 answers

### Q1 (draft P1). Does the leak/HVAC part of USUM take the owned-face mask?

**Yes.** In FDS it already does: `USUM` sums `U_NORMAL(_S)·B1%AREA` only over `SOLID_BOUNDARY` wall cells whose gas cell is not
covered (divg.f90:757-767, test at 760), and leak/HVAC faces are such wall cells. What FDS lacks is the same mask on `NODE_AREA`
and the area-weighted node state (`HVAC_BC_IN`, hvac.f90:2299-2539) and on the faces where `CALC_HVAC_BC` applies `MFT`
(wall.f90:1745-1796). So the AMR rule is: **one owned-face mask for `NODE_*`, the `MFT` application and `USUM`**, each summed with the
FR-005 (ii) exact accumulation. Only then do Σ_owned MFT·A = DUCT_MF and the zone budget in `D_PBAR_DT` (divg.f90:1519) describe the
same faces. Otherwise D-032's mean removal would silently absorb the mismatch (see Q3).

Covered coarse leak/HVAC faces carry no flux of their own: their normal velocity is overwritten by `average_down_faces` (D-032 (3)).

**Amendment** (applied in `01-amr-mapping-spec.md` §F, replacing the "leakage unchanged" sentence, former line 363):

> **HVAC and leakage.** Leak/HVAC faces are prescribed-normal-velocity faces. `CALC_HVAC_BC` sets `U_NORMAL(_S)` at every stage
> (wall.f90:1745-1796, called from `WALL_BC`, wall.f90:178). So in MAC form they enter `Ũ` as boundary-face velocities (homogeneous
> Neumann, first row of the table). They also enter the zone sums (`USUM`, divg.f90:757-767, giving dP̄/dt, divg.f90:1519) and D through
> `Q_LEAK` (divg.f90:543-544). The zone-pair leak areas (`LEAK_AREA`, hvac.f90:3374-3375; localized `AREA`, hvac.f90:684-685;
> pressure scaling, hvac.f90:3585-3609) are hierarchy-independent scalars. They are never summed from faces.
>
> **[REC-F1] Owned-face mask.** Every face-distributed leak/HVAC quantity uses **owned faces only**: faces whose gas-side cell is valid,
> uncovered (`makeFineMask`) and not solid. This covers `NODE_AREA`, the area-weighted node state, the per-face `MFT`
> application and the leak/HVAC part of `USUM`. All use `B1%AREA` of the owning level and the FR-005 (ii) exact sum.
> `USUM` already has this mask in FDS (divg.f90:760); `HVAC_BC_IN` does not (hvac.f90:2299-2539, no `INTERPOLATED_MESH` test),
> so the mask is new for the node sums and the flux application.
>
> Covered coarse leak/HVAC faces get no `CALC_HVAC_BC` flux. Their velocity is `average_down_faces` of the fine faces
> (D-032 (3), §C.2). Check: Σ_owned MFT·A = DUCT_MF per step to round-off (FR-020). Details:
> `03-fr034-fr040r3-answers.md` Q1–Q4.

### Q2 (draft P2). Per-level PBAR(KK,zone) or one composite profile?

**Per-level arrays are acceptable, but they must be samples of one composite definition.** Neither per-level integration nor a
single stored z-array is acceptable. Here is what FDS does:
* `P_0(z)` is one global profile. It is either an analytic/ramp hydrostatic profile, evaluated at each mesh's own `ZC(K)` when
  `STRATIFICATION` is on (init.f90:452-478), or `P_INF` otherwise (init.f90:480-484).
* `PBAR(0:KBP1,0:N_ZONE)` is allocated per mesh and initialised to `P_0(K)` (init.f90:1044-1054).
* `PBAR` then advances **uniformly in z** by the zone scalar:
  * predictor: `PBAR_S(:,I) = PBAR(:,I) + D_PBAR_DT(I)·DT` (mass.f90:548);
  * corrector: `PBAR(:,I) = ½(PBAR + PBAR_S + D_PBAR_DT_S·DT)` (mass.f90:730).
* `D_PBAR_DT` is a globally reduced per-zone scalar (divg.f90:1519 from the reduced DSUM/USUM/PSUM, main.f90:2028-2054).
* So every mesh holds **the same function** `PBAR(z,zone,t) = P_0(z) + ΔP_zone(t)`, sampled at its own cell heights.

**[REC] Rule.**
* Store one scalar ΔP_zone(t) per zone, advanced once per stage from the composite `D_PBAR_DT` (FR-034).
* Keep `P_0(z)` as a function (the ramp), not as a level-0 array.
* Each level (or box) evaluates `PBAR_ℓ(k,zone) = P_0(z_k^ℓ) + ΔP_zone`, together with its own `R_PBAR` and `PBAR_S`, at its own cell
  centres (and at faces where FDS averages, e.g. `P_AVE`, hvac.f90:2522-2530).

**Why.**
* This is layout-independent by construction: a box split or regrid changes only where the function is sampled. It is identical to
  FDS multi-mesh semantics, and it gives the fine level the correct stratification at its own z (O(Δz²) sampling differences only).
* Per-level *integration* of `D_PBAR_DT` would break FR-034. The composite RHS compatibility of §C.4 (divg.f90:1540) needs **one**
  `D_PBAR_DT` per zone.
* A single level-0 z-array interpolated to fine levels would add an unnecessary interpolation error in `R_PBAR`, `RHO_F` and `NODE_P`.

**For HVAC nodes**, `P_AVE` comes from the owning level of each face (Q1). The resulting `NODE_P` differs from single-mesh FDS
only by that level's z-sampling.

**Open point.** P_0 is built on a ramp table whose range depends on `ZS_MIN`/`NODE_Z_MIN` (init.f90:458-469). It must be built once
from global extents, not per level **[VERIFY: `RP%RDT` spacing vs the finest dz]**.

### Q3 (draft P4). Per-step leak coupling vs per-stage D_PBAR_DT: clash with D-032 mean removal?

**No clash, provided the same stage's face velocities are used everywhere in that stage.**

**Per stage, FDS does the following:**
* `WALL_BC` → `CALC_HVAC_BC` sets the leak/HVAC `U_NORMAL(_S)` (main.f90:838/1009).
* `DIVERGENCE_PART_1` (main.f90:840/1056) forms `USUM` from the same `U_NORMAL_S` (predictor) or `U_NORMAL` (corrector)
  (divg.f90:757-767). It then forms `D_PBAR_DT_P = (DSUM−USUM)/PSUM` (1519) and corrects D cell by cell with
  `DP −= (R_PBAR − RTRM)·D_PBAR_DT` (1540).
* The predictor uses `D_PBAR_DT_S`, the corrector `D_PBAR_DT` (1513-1514).
* After that correction, per zone, ∫D dV = USUM **of the same stage**. The Poisson BC on those faces is the same stage's prescribed
  velocity, so the per-zone RHS is compatible up to round-off and discretisation mismatch. That residual is exactly what REC-C2 / D-032 (2) removes.

**Freezing `DUCT_MF` for the step is not a compatibility issue.** `HVAC_CALC` returns early in the corrector (hvac.f90:1418-1425), so the
corrector reuses the predictor's `DUCT_MF`. That is a first-order-in-time lag in the leak law, which FDS has too. It does not
break the balance, because `U_NORMAL` is re-evaluated per stage from that `DUCT_MF` with that stage's ρ_F.

**[REC] Order within each stage:**
1. owned-face `CALC_HVAC_BC`;
2. exact per-zone DSUM/PSUM/USUM over owned faces and uncovered cells;
3. `D_PBAR_DT`;
4. the D correction (divg.f90:1540);
5. RHS assembly;
6. D-032 composite per-zone mean removal;
7. solve.

**Risk to guard.** Mean removal hides a face-set mismatch. If `USUM` were summed on a different face set from the one carrying `MFT`
(e.g. all levels vs owned), the per-zone RHS mean would be ΔU·A_covered/V_zone instead of round-off, and D-032 would remove it
without complaint. The pressure would look fine while mass is lost against the network.

**[REC] Diagnostic:**
* Log the removed mean per zone, relative to the RHS norm, at the FR-039 true-residual trigger points (first solve of the run, first
  solve after a regrid, every solve in debug).
* WARN when it exceeds 1e-10·‖b‖ (estimate: it should sit at ~1e-14 with exact sums).
* This is a proposed addition to FR-039's check, not a change to D-032.

### Q4 (draft P5). One leak node per zone pair with faces on several levels?

**Yes.** The pressure spec does not assume leak faces sit on one level.
* A node is a network object. FDS already assembles each node from faces on many meshes by a global reduction
  (`EXCHANGE_HVAC_BC`, main.f90:828) and solves the network once (rank 0, main.f90:829). Levels are just more face owners.
* The only AMR-specific requirement is Q1's owned-face mask, so that no geometric face is counted twice.
* The consistency rule ERROR(552) ("Ductnode must lie with a single pressure zone" [sic], hvac.f90:2410) must be evaluated on
  owned faces with each level's zone map. This works because R3-T (draft §1.3 (e)) guarantees a fine cell and its covering coarse
  parent in the band have the same zone.
* Pressure-side needs:
  * the node's `NODE_P` uses the owning level's `PBAR` sample (Q2);
  * the zone's `USUM` gets each owned face once (Q1).
* No per-level node, no flow partitioning, no network change.

### P3 and P6 (not in the six questions, answered briefly)

**P3 (Phase placement).**
* Nothing in the pressure design needs Phase 6 for multi-zone. There is no subcycling, and the composite per-zone
  DSUM/PSUM/USUM, `D_PBAR_DT` and mean removal work for any number of zones in one composite solve.
* Phase 6's dP0/dt sync (FR-034) exists only for subcycling.
* **Recommendation:** multi-zone with leakage gates **Phase 4** with the refined variants of zone_shape_2 / HVAC_leak_exponent. The
  uniform variants gate Phase 2 via FR-006.
* This resolves the requirements.md:179 "TBD(Pressure Lead)" on the pressure side. The Spec Lead updates FR-034. I did not edit it.

**P6 (zone numbering and seeding).** The Integration Lead should own this: it is setup/driver code (flood fill, seeding, numbering) and
FR-005 layout-independence. The Pressure Lead is a consumer. What the pressure side needs is: zone IDs identical across layouts, ranks
and regrids, and identical between a fine cell and its covering parent in the band.

---

## 2. FR-040 R3 answers

### Q5 (draft §3.2 Q1). Which β on covered faces and interface-band faces?

**Concrete rule.**
* **Under E-1 (the Phase 4 baseline, A2):**
  * no rule is needed: the operator is constant-coefficient (`MLPoisson`, β ≡ 1 everywhere, H formulation);
  * thin walls act only through the wall-face velocity forcing (`NO_FLUX`, velo.f90:1348-1563) and the obstruction iteration (01 §D),
    applied by each level on its own **owned** wall faces;
  * covered coarse faces take `average_down_faces` velocities (D-032 (3)).
* **Under E-2:**
  1. Each level sets β_face = 0 on every wall face (thin or solid) of **its own R1 mask** (D-009), and β_face = 1 elsewhere
     (01 §E). Do this on all faces of the level, covered ones included, for simplicity.
  2. **Covered coarse faces and coarse faces on the C/F boundary:** AMReX replaces them with the arithmetic mean of the r² fine
     subfaces (read in AMReX, see §0 last row). The driver cannot choose anything else short of patching AMReX, so the rule simply
     *is* "average of the fine β". With R3-T T1 in the band, a coarse face in the band is a wall exactly when all its fine subfaces are
     walls, so the average is exactly 0 or exactly 1 there. The band has one consistent β, as the draft wants. Outside the band, the
     M1–M4 covered mismatches give fractional β (e.g. ½). That affects only the coarse AMR level's smoothing under the fine level,
     never the composite residual, whose C/F fluxes come from the fine level (reflux) and whose covered residual is replaced by the
     averaged fine residual.
  3. **Uncovered coarse faces in the band** (coarse side, outside the fine patch): β from the level-ℓ R1 mask. R3-T T1/T2
     guarantees these match the fine mask where it overlaps (valid+2 ghosts).
  4. **MG coarsening within a level:** `averageDownCoeffsSameAmrLevel` also averages face β arithmetically
     (AMReX_MLABecLaplacian.H:711-737, `average_down_faces` at 736).
     * A thin wall lying on a face plane *interior* to a coarsened MG cell disappears on that MG level.
     * A partial wall gives fractional β.
     * This is the "leak on coarse MG levels" of 01 §E. It changes the convergence rate only; the fine-level operator (and hence the
       answer at the solver tolerance) is unaffected.
     * Caveat: if thin walls split a level into **disconnected components**, coarse MG levels reconnect them. The fine operator
       then has one null vector per component, while MLMG's singular handling knows one constant per level. D-032 (2)'s per-zone
       mean removal makes the RHS compatible per zone, and a per-component pin (01 §E, "one mask = 0 cell per zone") makes each
       component non-singular. **[VERIFY]** convergence rate on such a case; it belongs in the E-2 test list (e.g. zone_shape_2
       with E-2).
* This resolves the 01 §E "[VERIFY] coarse-level β averaging behaviour" (applied in 01).

**Correction to the draft's §1.3 (e) bullet "Covered coarse faces: β from the R1 coarse mask":** the driver's covered-face value is
overwritten by MLMG. It should read "Covered and C/F coarse faces: β = arithmetic mean of the r² fine subfaces (set by
`MLABecLaplacian::averageDownCoeffs`); under R3-T this equals the R1 coarse mask inside the band."

### Q6 (draft §3.2 Q2). C/F interpolation across a thin wall normal to the interface: is a masked interpolant needed?

**What MLMG does at C/F (read in AMReX).**
* **Normal direction (fine side).** `applyBC` fills the fine ghost cell from the boundary value on the C/F face plus interior fine
  cells along the normal, using `maxorder`. With `setMaxOrder(2)` (D-032 (1)) that is a linear fit through the face value and **one**
  interior fine cell (AMReX_MLLinOp_K.H:124-138, `NX = min(blen+1, maxorder)`; call sites AMReX_MLCellLinOp.H:986-1092). This
  stencil never crosses a wall normal to the interface.
* **Tangential direction (coarse data onto the C/F face).**
  * `InterpBndryData::setBndryValues` with `IBD_max_order_DEF = 3` is hard-coded at every MLMG call site (AMReX_MLCellLinOp.H:784,
    862, 962, 977; AMReX_InterpBndryData.H:126). The order-3 path (AMReX_InterpBndryData.H:191-260) calls `interpbndrydata_{x,y,z}_o3`.
  * For a face normal to x, the kernel uses the coarse cell behind the face `(ic,jc,kc)`, its tangential neighbours `jc±1`, `kc±1`
    and the four diagonal corners. It forms central slopes, second differences and the cross term
    (AMReX_InterpBndryData_3D_K.H:23-51).
  * A neighbour is dropped (one-sided slope, zero curvature) **only** when the boundary mask says it is not `not_covered`, i.e. when
    it is covered by same-level fine grids or outside the domain (32-33, 38-39, 44-45; mask values AMReX_BndryData.H:49-51).
  * The mask knows nothing about thin walls or β.
  * **So yes: a thin wall normal to the interface between coarse cells `jc` and `jc+1` is crossed by the tangential slope and
    curvature.** R3-T does not prevent this: it keeps the wall on a coarse-face edge, which is exactly the position where the stencil
    straddles it.
* **β does not enter.** `setBndryValues`/`updateBndryValues` take the coarse data, ratio and order only (AMReX_MLCellLinOp.H:782-785,
  859-862, 961-963, 976-978). Zero face β in `MLABecLaplacian` changes the operator's fluxes, not the C/F ghost stencil.
* **FillPatch** (`mf_cell_cons_interp`, the default cell-conservative linear interpolater) also takes central slopes from i±1. But
  they are **limited** (minmod/MC with the one-sided differences, AMReX_MFInterp_3D_C.H:113-137 for the `llslope` variant). With a
  jump on one side and smooth data on the other, the limiter caps the slope at twice the smooth one-sided difference. The
  contamination is then O(h·|∇φ|_smooth), not O(jump). Also, the projection does not use FillPatch for H: the C/F fluxes come from
  `getFluxes` plus `average_down_faces` (D-032 (3)).

**How big is the error, and does it break anything?**
* It does not move mass across the wall. The contaminated values are ghost values on the fine side of C/F faces that lie wholly on
  one side of the wall (R3-T). The wall faces themselves carry zero flux (E-2: β = 0; E-1: forced U = 0), and reflux conserves.
* It is a local accuracy error in the fine fluxes through the C/F faces adjacent to the line where the wall meets the interface.
  The ghost error is ≈ |y·dy + y²·dy2| with dy ≈ J/2, i.e. ≈ J/8 + J/32 at ratio 2 (y = ±¼), where J is the H jump across the wall
  **[derived]**.
* **J is small in FDS's formulation.** Zone-to-zone pressure differences are carried by `PBAR(z,zone)`, not by H (Q2; divg.f90:1540,
  mass.f90:548/730), so H jumps across a thin wall are dynamic-pressure scale, O(ρu²)/ρ **[derived]**. Example: the
  duct_flow_uglmat_refine riser, 1 m³/s through 1 m², gives u ≈ 1 m/s and J ≈ 0.5 m²/s².
* Under E-1 the all-cell operator already couples H across thin walls exactly as FDS's FFT does, so the stencil is "as FDS" except for
  order (FDS uses the piecewise-constant `IIO..KKO` average, velo.f90:1390-1399).
* Reversing that: FDS's own C/F coupling is first order, while ours is higher order everywhere except O(J) near wall/interface
  junctions.

**Options considered.**

| Option | What | Verdict |
|---|---|---|
| (i) Tagging/blocking rule: thin walls must not cross C/F boundaries | forbid crossings for dynamic runs (tag the whole wall) | Not possible for static inputs: duct_flow_uglmat_refine's riser crosses z = 3.0 by construction (draft §1.1). Useful as a *preference* in the dynamic tagger, not as a rule. |
| (ii) R3-T only, keep MLMG as is | accept the O(J) local ghost error at wall/interface junctions | **Recommended** (with the acceptance check below). No AMReX change. Conservation and zone integrity are unaffected. |
| (iii) Masked or order-1 tangential interpolant in the band | pass a wall-aware mask, or use order 1 (`interpbndrydata_o1`, piecewise constant, which exists at AMReX_InterpBndryData.H:211-217) next to thin faces | **Fallback.** Needs an AMReX patch, because the order is hard-coded at MLCellLinOp:784/862/962/977 and the mask type has no wall state. An upstream-able change is to make the order a per-operator setting, like `setMaxOrder`, costing ~20–50 lines **[estimate]**. Order 1 everywhere would lower the P2-measured C/F order (≥ 1.92 at maxorder 2 with order-3 tangential) toward 1. So if needed, apply it only on faces whose tangential stencil crosses a wall face. |
| (iv) Zero-β faces | set β = 0 on the wall faces | Does **not** change the interpolation stencil (see above). Necessary under E-2 for the operator, irrelevant for this issue. |
| (v) Masked slope limiter in FillPatch | wall-aware slopes in `mf_cell_cons_interp` | Not needed for pressure (FillPatch is not in the projection path). The limiter already bounds the error to O(h·∇φ). Scalars and velocity are the Integration Lead's call. |

**Recommendation.**
* Adopt (ii): R3-T as the geometric rule, the standard MLMG C/F interpolation, and E-1 as the Phase 4 baseline. Do not require a
  masked interpolant now.
* **Add one acceptance item** to the draft's §1.6 for duct_flow_uglmat_refine:
  * report max |J| across the riser walls on the z = 3.0 C/F faces;
  * report the flux error on the C/F faces adjacent to the walls (difference between the fine-face flux and the one computed with
    a piecewise-constant ghost, both at the solution);
  * keep `flow_in`/`flow_out` within 5 % and T2 vs the FireX UGLMAT baseline.
* **Trigger for option (iii):** those C/F faces carry > 1 % of the riser flow error, or T2 fails.
* E-2 changes nothing here, because the stencil does not depend on β.

---

## 3. Changes made to other files in this folder
* `01-amr-mapping-spec.md`:
  * §E MG-robustness row: the [VERIFY] is resolved with the AMReX facts of Q5;
  * §F: the line-363 sentence is replaced with the HVAC/leakage paragraph and REC-F1 (Q1).
* `README.md`:
  * Q4 note (E-2 β rule);
  * new Q17 (FR-034 face mask and phase placement), Q18 (C/F tangential interpolation across thin walls) and Q19 (masked
    level 0, see `04-masked-domain-pressure.md`);
  * "Not verified" list updated (coarse-level β averaging is now verified).

## 4. Citation index
| Topic | Location |
|---|---|
| USUM covered mask, DSUM/PSUM | divg.f90:727-753 (734), 757-767 (760) |
| D_PBAR_DT, D correction, stage selection | divg.f90:1495-1506, 1513-1514, 1519-1520, 1540 |
| Leak/HVAC excluded from U_NORMAL prediction | divg.f90:1351-1352, 1415-1416 |
| Q_LEAK | divg.f90:543-544; wall.f90:1760, 1773 |
| CALC_HVAC_BC | wall.f90:178, 1745-1796 (1761, 1784-1791), 1446-1447 |
| HVAC_BC_IN | hvac.f90:2299-2539 (2343, 2376, 2410, 2415, 2452-2456, 2494, 2513, 2522-2533) |
| Leak duct area, pressure scaling, localized leak, collapse | hvac.f90:3366-3380 (3374-3375), 3585-3609, 684-685, 3125-3134 |
| HVAC_CALC corrector early return | hvac.f90:1418-1425 |
| Time loop | main.f90:748, 772, 810-813, 822, 828-830, 838, 840, 845, 855, 928, 980, 1009, 1056, 1065, 1090 |
| P_0, PBAR init and update | init.f90:452-484, 1044-1054; mass.f90:548, 730 |
| Ghost H over IIO..KKO; DHFCT | velo.f90:1388-1399, 1483-1487 |
| MLABecLaplacian β averaging | AMReX_MLABecLaplacian.H:474-510, 693-707, 711-737, 799-815 |
| average_down_faces kernel | AMReX_MultiFabUtil_3D_C.H:95-120, 173-188 |
| C/F tangential interpolation (order 3, mask) | AMReX_InterpBndryData.H:126, 186-260; AMReX_InterpBndryData_3D_K.H:23-51; AMReX_BndryData.H:49-51; AMReX_MLCellLinOp.H:775-790, 855-866, 955-980 |
| Normal ghost extrapolation | AMReX_MLLinOp_K.H:124-138; AMReX_MLCellLinOp.H:986-1092; AMReX_MLLinOp.H:310, 886 |
| FillPatch limited slopes | AMReX_MFInterp_3D_C.H:113-137 |
