# ADR-003: Solid geometry on AMR levels — FDS OBST/WALL staircase vs AMReX EB vs porting GEOM/CC_IBM

| Field | Value |
|---|---|
| Status | Proposed (framing memo; leaning recorded, not decided). v0.2 records accepted sub-decisions D-009 and D-010; v0.2.1 records charter Q4 as answered (D-033), the D-009 amendment R3-T (FR-040 R3) and a D-010 clarification (see "Accepted decisions") |
| Version | v0.2.1 (2026-09-25): Q4 answered by D-033 (thin OBSTs required with refinement; GEOM/CC_IBM and HT3D deferred); zero-thickness citation corrected; D-009 amendment R3-T (thin OBSTs at coarse/fine interfaces); D-010 clarified (static crossing allowed, regrid cut forbidden); see "v0.2.1 changes". v0.2 (2026-09-25): D-009 and D-010 accepted |
| Type | Framing memo |
| Date | 2026-09-25 |
| Owner | Chief Architect, with AMR Pressure Solver Lead and AMReX Integration Lead |
| Deciders | Project owner (Q4 scope); AMR Spec & Program Lead; AMR Pressure Solver Lead |
| Blocks | **ADR-002** (whether the pressure iteration survives for solids) and **ADR-001** partially (EB is C++-only; wall/particle state is the shim's hardest data) |
| **Code base** | **FireX `36975d765f` (2026-09-24), local branch `AMReX`, this repository (read-only)** — primary code base and acceptance baseline. All `Source/` line numbers are FireX. |
| Other evidence | AMReX `99ddfda` (`(local AMReX checkout)`); ERF, incflo, PeleLMeX clones under `ref/` |

Cross-references: risks R-06, R-07, R-15, R-23, R-26, R-33; requirements FR-004, FR-035, FR-040, FR-041a/b, FR-042..044, FR-022; decisions D-009 (amended v0.2.1, R3-T), D-010 (clarified v0.2.1), D-030, D-032, D-033; charter Q4 (answered, D-033); pressure spec `docs/pressure/01-amr-mapping-spec.md` §D, §E (E-1/E-2/E-3, REC-E1); AMReX mapping `docs/amrex/mapping.md` §5. Note: FR-042 cites `main.f90:1651` and FR-035 cites `read.f90:10045` (merge-base `ce1f659` numbering); FireX `CREATE_OR_REMOVE_OBSTRUCTIONS` is `main.f90:1782-1798`.

## v0.2.1 changes (2026-09-25)
- **Q4 answered (owner decision D-033):** thin OBSTs, HVAC, pressure ZONEs and level-set wildfire must work with refinement; GEOM/CC_IBM and HT3D are deferred, not permanently dropped. Thin walls are removed from every "open"/"FR-004 rejection" statement (risks, Spec Lead questions, "Decision needed from the project owner"). Under the Q4 framing below, "no" for GEOM confirms Option B with D deferred; formal acceptance of the ADR as a whole stays with the Chief Architect.
- **Citation fix (Context):** `init.f90:118-122` is the HT3D-only `INIT_THIN_WALL_CELL` path (`HT_DIM>1`), not general thin-OBST handling. The zero-thickness bullet now cites the snap and collapse code (`read.f90:11157-11162, 11167-11173, 11174-11216, 11220-11224, 11315-11317`), the no-blocked-cell rule (`read.f90:11708`, `func.f90:5503-5526`) and the thin-face wall cells (`init.f90:195-300, 3315`). Also corrected: the `NINT` snap is `read.f90:11157-11162` (was 11155-11158) and `OBSTRUCTION_TYPE%THIN` is `type.f90:1135` (was 1136).
- **New accepted-decision text: D-009 amendment R3-T** (FR-040 R3, thin OBSTs that cross or touch a coarse/fine interface): interface-band exact-refinement check, setup error for static inputs, footprint-plus-one-coarse-cell tagging for dynamic refinement, rebuild and re-check after every regrid. Two sub-items are **pending the AMR Pressure Solver Lead** (β on covered/band faces; masked C/F interpolant). Source analysis: `adr/drafts/rulings-FR040R3-FR034.md` Ruling 1 (draft, kept as-is).
- **D-010 clarified:** a stateful wall crossing a C/F interface in a static hierarchy is allowed (no wall-state transfer happens); only a regrid that moves a C/F boundary through a stateful face is forbidden.
- Follow-ups outside this ADR (Spec & Program Lead): align requirements.md FR-040 R3 and the FR-040 "Open point R3", and the FR-041a "Rule" wording, with v0.2.1.

## v0.2 changes (2026-09-25)
- New section "Accepted decisions (v0.2)": D-009 (FR-040 snap rule: unchanged FDS rule per level) and D-010 (Phase 5 stateful-wall coverage freeze, FR-041a), with the breadth of the freeze recorded (R-33/FR-041a widened); FR-041b's target is set after Spike G2's cost estimate.
- Existing text brought in line: the "never refine across burning OBST faces" alternative is now the accepted Phase 5 rule and covers every stateful wall, not only burning ones; the snap rule is no longer "choose in G1" (G1 can overturn D-009); G1/G2 pass signals and the Spec Lead questions are updated.
- The pressure-solver, refinement-ratio, reflux and decomposition items accepted the same day are recorded in ADR-002.

## Context

**FDS rectilinear geometry (OBST/WALL), verified.**
- `OBSTRUCTION_TYPE` (`type.f90:1082-1147`) stores integer extents `I1..K2` **per mesh**. They are snapped with `OB%I1 = NINT(GINV(XB1-XS,1,NM)*RDXI)` (`read.f90:11157-11162`, inside `READ_OBST` `read.f90:10742-11726`), so the same OBST can occupy different physical extents on meshes of different resolution.
  - Illustration (arithmetic, not a measured case): with mesh origin 0, x∈[0.375,1.125] m snaps to [0.5,1.0] at dx=0.5 but to [0.5,1.25] at dx=0.25 (values exact in binary; Fortran `NINT` rounds halves away from zero). A child's solid mask averaged down would not match the parent's → FR-040 must allow documented mismatches.
  - An OBST thicker than a fine cell but thinner than a coarse cell becomes **zero-thickness** on the coarse level only. After the per-mesh `NINT` snap (`read.f90:11157-11162`), the collapse branch for OBSTs thinner than 0.25 cell sets `I1==I2` (`read.f90:11174-11216`; the `THICKEN` alternative is `read.f90:11167-11173`; an OBST zero-thick in ≥ 2 directions is discarded, `read.f90:11220-11224`), and `OB%THIN` (`type.f90:1135`) is set when the input length is < 0.5 cell (`read.f90:11315-11317`). A zero-thickness OBST blocks no cells (`BLOCK_CELL` fills `I1+1..I2`, `read.f90:11708`; `func.f90:5503-5526`). It exists only as faces: two wall cells per face, one per gas side (OBST face loops, `init.f90:195-300`), flagged `WC%THIN` when both sides are gas and not exterior (`init.f90:3315`). `INIT_THIN_WALL_CELL` (`init.f90:118-191`) is created only for `HT_DIM>1` (`init.f90:118`), i.e. the HT3D path, deferred with HT3D (D-033). So topology, not just extent, is level-dependent.
- Cells are marked through `CELL_TYPE` (`type.f90:2175-2187`, incl. `THIN_WALL_INDEX/THIN_SURF_INDEX/THIN_OBST_INDEX`) via `BLOCK_CELL` (`func.f90:5503-5526`).
- Every gas/solid face is a `WALL_TYPE` (`type.f90:434-457`) pointing into `BOUNDARY_ONE_D_TYPE` (`type.f90:217-264`) and `BOUNDARY_PROP1` storage: 1-D conduction profiles, pyrolysis mass, surface temperature. Built once by `INIT_WALL_CELL` (`init.f90:2975-3386`); applied each step in `WALL_BC` (`wall.f90:26-279`). HT3D adds 3-D solid conduction (`type.f90:259, 288, 931`).
- Geometry changes at runtime: `CREATE_OR_REMOVE_OBSTRUCTIONS` (`main.f90:1782-1798`) → `CREATE_OR_REMOVE_OBST` (`init.f90:4852-4891`) → `REASSIGN_WALL_CELLS` (`init.f90:4898-5228`) → `GLOBAL_MATRIX_REASSIGN` (`main.f90:1806-1835`). **This is existing prior art for rebuilding wall cells and carrying state mid-run**, which is exactly what a regrid needs.
- Coarse-fine interfaces between meshes are themselves wall cells (`EXTERNAL_WALL_TYPE`, `type.f90:462-482`; ghosts averaged over `IIO..KKO` ranges, `wall.f90:321-330`). NIST's 2014 SAMR notes planned embedded-mesh interfaces the same way (wall-cell indexing `II`/`IIG` at the embedded boundary, `Manuals/FDS_Technical_Reference_Guide/samr_notes.tex:50-110`). They contain nothing on obstructions crossing refinement.
- Pressure: the solid treatment decides whether the pressure loop survives (spec §D reason 2). FDS default FFT/GLMAT = IBM forcing on an all-cell operator (`NO_FLUX`, `velo.f90:1348-1563`); ULMAT/UGLMAT = solids removed from the operator (UG `FDS_User_Guide.tex:9455`).

**FDS complex geometry (GEOM/CC_IBM).** `geom.f90` 27,735 lines + `ccib.f90` 24,046 lines; `CFACE_TYPE` `type.f90:1361-1382`. GEOM forces UGLMAT (`FDS_User_Guide.tex:1685`). Deferred, not permanently dropped (owner decision D-033; FR-004, FR-044, R-07).

**AMReX EB, verified.**
- `EB2` implicit functions `IF_Box`/`IF_Union`/`IF_Plane` and STL (`Src/EB/AMReX_EB2_IF_Box.H`, `AMReX_EB2_IndexSpace_STL.H`), `EBFluxRegister` (`Src/EB/AMReX_EBFluxRegister.H`), state redistribution (`AMReX_EB_Redistribution*`, `StateRedist*`), EB MLMG operators (`MLEBABecLap`, `MLEBTensorOp`).
- EB is built on the finest level and **coarsened**; "coarsening could create multi-valued cells even if the fine level does not have any ... multi-valued cells are not supported, it will cause a runtime error" (`Docs/sphinx_documentation/source/EB.rst:173-177`). Thin FDS walls and thin plates near other solids are exactly the shapes that go multi-valued when coarsened.
- EB is built once at init; runtime changes mean rebuilding index space and factories on all levels (mapping.md R5.1). No EB bindings in `F_Interfaces` (mapping.md §5), so EB implies ADR-001 Option A.
- Precedent: ERF represents buildings through EB (`ref/ERF/Source/EB/ERF_EBIFBuildings.H`); incflo and PeleLMeX are EB-capable. None has burn-away solids or 1-D in-depth solid conduction per wall face (not found in the clones; evidence thin).
- Non-EB solid removal in MLMG: `MLCellABecLap` overset mask (`Src/LinearSolvers/MLMG/AMReX_MLCellABecLap.H:70-71,128`) supports E-2; not exposed in `F_Interfaces`. The Integration Lead suggests emulating it from Fortran with β=0 on solid faces and `acoef`=1 in solid cells (`docs/amrex/driver-options.md` row 4, **unverified**). If that works, E-2 does not force C++, and **only EB (E-3) forces ADR-001 Option A** (driver-options §1.1).

## Decision drivers (ranked)
1. **Fidelity to FDS V&V** for rectilinear cases (FR-001..003, FR-022 walls): the bulk of the FDS validation suite is OBST-based.
2. **Wall state survives regrid** (FR-041) and runtime create/remove works (FR-042).
3. **Mask consistency across levels** (FR-040) and solid velocity error (FR-035).
4. **Pressure coupling**: which operator (E-1/E-2/E-3) and whether the loop survives (ADR-002).
5. Effort and schedule (Phase 5 per roadmap), and future complex-geometry scope (Q4).

## Options

### Option A — AMReX EB cut cells for all solids
- **Pros:** single-valued cut cells, conservative C/F with `EBFluxRegister`, EB-aware MLMG (no pressure loop for solids; FR-035 "zero" branch); regrid geometry for free; path toward GEOM-like fidelity.
- **Cons:** FDS walls are faces of blocked cells, not cut cells, so V&V comparability is lost for the OBST suite; multi-valued coarse cells from thin walls (`EB.rst:173-177`) → coarse-level topology must be simplified; burn-away needs EB rebuilds; `BOUNDARY_ONE_D` state has no home in EB data (per cut face `bcent/bnorm` only); C++ only; small-cell redistribution changes numerics everywhere near solids.

### Option B — Keep FDS OBST/WALL on each level (staircase), rebuilt at regrid (leaning)
Each level has its own solid mask (`iMultiFab`) and face masks from the global OBST list; wall records per box; pressure via E-1 (IBM on all-cell `MLPoisson`) then E-2 (overset mask / β=0 faces).
- **Pros:** bit-comparable to FDS on single level; reuses `BLOCK_CELL`, `INIT_WALL_CELL`, `WALL_BC`, `REASSIGN_WALL_CELLS` logic; create/remove stays a mask flip plus wall reassignment; matches pressure REC-E1 and mapping.md §5 "Phase 1: no EB".
- **Cons:** per-level snapping mismatches and level-dependent `THIN` topology (FR-040 needs a documented-mismatch policy); wall state must be re-homed at every regrid (R-06) — the dominant new data model; under E-1 the pressure loop survives for solids (R-15); staircase accuracy unchanged from FDS.
- **Sub-choice for wall state (FR-041):** key 1-D solid state to the OBST surface patch at the **finest resolution that has ever covered it**, not to the box (mapping.md §5 proposal). A coarse face then aggregates r² fine records (average for BCs, energy-conserving), and refinement copies/splits. Phase 5 instead uses the accepted coverage freeze (D-010, FR-041a): no coarse-fine boundary across any stateful wall face, which (breadth, above) means every wall with a material layer, not only burning faces. Cheaper, but a tagging constraint users will hit (R-33). The patch-keyed transfer is FR-041b, after Phase 6.

### Option C — Port FDS GEOM/CC_IBM onto AMReX levels
- **Pros:** keeps FDS's own cut-cell numerics and `CFACE` solid model.
- **Cons:** ~51.8k lines built around per-mesh unstructured cut-cell tables and their own exchange; forces UGLMAT-like operator; no AMReX precedent. Largest option by far; only justified if Q4 makes GEOM mandatory.

### Option D — Hybrid: OBST staircase (B) now; EB (A) later as the target for GEOM only
EB is introduced only if/when complex geometry is in scope, feeding EB2 from STL/GEOM input and dropping `geom.f90` intersection. OBSTs stay staircase so the rectilinear V&V suite remains comparable.

## Leaning
**Option B now, D as the long-term path.** Rectilinear OBSTs as per-level masks rebuilt at regrid, with pressure **E-1 first, then E-2** (pressure spec REC-E1); EB only together with a GEOM decision (Q4), and not for OBSTs. Snap rule: unchanged FDS rule per level (D-009). Phase 5: stateful walls keep their level (D-010, FR-041a); later, wall state keyed to OBST surface patches, not boxes (FR-041b). Thin/zero-thickness walls are a **blocker for EB** at coarse levels and a documented-mismatch case for B; they are required with refinement (D-033) and handled under B by D-009 plus R3-T.

Strongest evidence: FDS's V&V is OBST-based and snapping is already per-mesh (`read.f90:11157-11162`); FDS already rebuilds walls and matrices at runtime (`init.f90:4898-5228`, `main.f90:1806`); EB forbids the multi-valued cells that coarsened thin walls produce (`EB.rst:173-177`). Confidence: medium-high for "not EB for OBST"; medium for the wall-state data model (unmeasured; Spike G2).

**Consequences for ADR-002:** under B+E-1 the pressure loop survives for solids (plus the baroclinic lag under A2-a); under B+E-2 or A it survives only for the baroclinic lag. **For ADR-001:** B keeps the `POINT_TO_BOX` shim viable for walls (per-box `WALL` lists), but wall/particle state is the shim's hardest data (R-26); A would force Option A's C++ driver and new kernels.

## Accepted decisions (v0.2, 2026-09-25; amended v0.2.1)
These are accepted sub-decisions. The ADR as a whole stays Proposed; Q4 is answered (D-033), and formal acceptance of the ADR is for the Chief Architect.

- **D-009, FR-040 snap rule (accepted).** Each level snaps OBSTs with the unchanged FDS rule (`read.f90:11157-11224, 11315-11317`) using that level's dx and the level-0 origin. The fine level is authoritative under cover. No mismatch is allowed on coarse-fine interface faces; the regridder grows around mismatches. Allowed mismatches M1-M4 (FR-040) occur in covered regions only. This keeps level 0 identical to baseline, including the two upstream slips at `read.f90:11193` and `:11210` (A-20). Spike G1 can overturn it (replacement: "finest-level extent + >50 %") if G1 shows the no-mismatch rule on interface faces cannot be met.
- **D-009 amendment R3-T, thin OBSTs at coarse/fine interfaces (FR-040 R3; v0.2.1, Chief Architect ruling).** Source analysis and citations: `adr/drafts/rulings-FR040R3-FR034.md` Ruling 1. Applies to every OBST, and is needed in particular for thin walls normal to an interface: there every interface face is gas–gas, so D-009's "no mismatch on interface faces" holds, yet a coarse C/F face that straddles the wall would, through `average_down_faces` and reflux (D-032 (3), FR-024), move mass across a zero-thickness wall (FDS analogue: ghost averaging over `IIO..KKO`, `velo.f90:1388-1399`).
  - **Ownership.** Each level builds its face mask from the global OBST list (not per FDS mesh). Every face belongs to the finest level whose valid region contains it: C/F and covered coarse faces take `average_down_faces` of the fine faces. A wall face belongs to the level whose valid, uncovered cell is its gas side (the FDS covered-cell test, `divg.f90:760`).
  - **Interface-band exact-refinement check.** The band of a level pair (ℓ, ℓ+1) is every coarse cell within one coarse cell of the C/F boundary, on both sides. In the band the fine face mask must be the exact refinement of the coarse face mask: (T1) a coarse face is a wall face (thin or solid) if and only if all r² fine subfaces are wall faces; (T2) no fine wall face lies on a fine face plane interior to a coarse cell. A thin wall crossing the C/F plane therefore meets it along coarse-face edges at the same coordinate on both levels. This is FDS's UGLMAT "solid wins" rule (`pres.f90:3726-3829`) turned into a check instead of a silent fix. M1-M4 remain allowed only in covered cells outside the band. Fine band cells must carry the same pressure-zone ID as their coarse parents (the zone flood fill stops at thin faces, `func.f90:5656-5676`).
  - **Static inputs: setup error.** An FDS multi-resolution input whose hierarchy is fixed by the input (FR-010) and that fails the check stops at setup with an error naming the OBST (input line/ordinal), the two levels, both snapped planes and the C/F plane. No automatic growth (FR-010: the hierarchy matches the input boxes) and no silent solidification; the user remeshes deliberately, as under D-030.
  - **Dynamic refinement: footprint plus one coarse cell.** The mismatch list is computed once at setup from the global OBST list for every level pair. The tagger refines the whole footprint (the union of both levels' snapped extents) of any OBST whose snap differs between ℓ and ℓ+1, plus one coarse cell. Local one-cell growth is not enough for a wall normal to the interface, because the crossing moves with the boundary. If the tagged region exceeds `MAX_LEVEL` or memory limits, the run stops with an error.
  - **Rebuild after regrid.** After every regrid: rebuild each level's masks from the global OBST list, re-run the band check and the zone check, and recompute the owned-face sums. A failure after a regrid is a tagger bug and stops the run. The geometry-dependent FDS checks ERROR(421)/(422) (`init.f90:1102-1115`) and `THIN` (`read.f90:11315-11317`) are evaluated on the level that owns the face.
  - **Pending the AMR Pressure Solver Lead** (not yet decided): (P-a) which β applies on covered coarse faces and on band faces under E-2 (pressure spec §E; the `[VERIFY]` on coarse-level β averaging at `pressure/01-amr-mapping-spec.md:338`); (P-b) whether a masked (or piecewise-constant) coarse/fine interpolant is needed in band cells next to thin faces, where MLMG's tangential C/F interpolation or FillPatch slopes would otherwise take coarse neighbours from across the wall.
  - **Test case.** `Pressure_Solver/duct_flow_uglmat_refine` is in scope with its static two-level hierarchy (level 1 over z ∈ [3.0, 6.2]). Its four thin walls (input lines 43-46) lie at x = 1.0/2.0 and y = 4.0/5.0, on coarse faces and normal to z = 3.0, so the check reports zero mismatches and zero band violations (derived from the input, not run). A derived copy with one wall moved to x = 1.07 (coarse snap 1.0, fine snap 1.1) must stop at setup. Full acceptance: draft §1.6.
- **D-010, Phase 5 stateful-wall coverage freeze (accepted for Phase 5).** Wall faces with solid-phase history keep their level; the coarse-fine boundary never cuts them, so there is no cross-level wall-state transfer in Phase 5 (FR-041a).
  - **Breadth:** any `SURF` with a `MATL_ID` carries 1-D conduction history from step 0, so nearly every non-inert wall is stateful before the first regrid. The freeze therefore effectively means static refinement near every wall with a material layer, not only near burning faces (R-33 and FR-041a widened). Acceptable for the Phase 5 burner case; it limits flame-spread use cases until FR-041b.
  - **FR-041b** (cross-level wall-state transfer) is deferred until after Phase 6. Its target milestone is set after Spike G2's cost estimate.
  - **Clarification (v0.2.1): static crossing allowed, regrid cut forbidden.** A stateful wall face may cross a C/F interface whose position across that face never changes during the run (e.g. a static hierarchy fixed by the input): each of its wall cells lives on one level for the whole run, so no wall-state transfer happens. What the freeze forbids is a regrid that creates, moves or removes a C/F boundary through a stateful face; that would need FR-041b. The R3-T band check still applies to the geometry. requirements.md FR-041a "Rule" ("the coarse-fine boundary never cuts a stateful OBST face") is to be aligned by the Spec & Program Lead.

## Rejected alternatives
- **A for OBSTs now:** loses OBST V&V comparability, breaks on thin walls, and does not carry 1-D solid state or burn-away.
- **C now:** size (~51.8k lines) and no precedent; only with an explicit Q4 mandate.
- **Snap every OBST to the level-0 grid on all levels** (identical masks, trivially FR-040): throws away refinement exactly where fires attach to walls.

## Consequences and risks
- **+** Rectilinear suite stays comparable; runtime create/remove keeps working; no C++ lock-in from geometry alone.
- **−** Wall-state re-homing at regrid is new code with conservation obligations (R-06, FR-041b, FR-022); deferred until after Phase 6 by the D-010 freeze, which in turn forces static refinement near every wall with a material layer in Phase 5 (R-33).
- **Risk: mask mismatch/topology change between levels** (FR-040). Mitigation: mask checker; the FDS rule per level is accepted (D-009) with no mismatch on interface faces, the R3-T interface-band check (v0.2.1), and M1-M4 only under cover outside the band. The "finest covering level defines the extent, coarse cells solid if >50%" rule (mapping.md §5) is the G1 fallback.
- **Risk: solid velocity error semantics** under E-1 (R-15, FR-035). Mitigation: E-2 right after E-1.
- **Risk: thin walls cross boxes and levels.** Thin OBSTs are required with refinement (D-033), so FR-004 rejection is not available for them. Mitigation: D-009 per-level snap plus the R3-T interface-band check, setup error for static inputs and footprint tagging for dynamic runs (v0.2.1). **HT3D** (3-D conduction through thin walls, `init.f90:118-191`) is deferred (D-033) and rejected in AMR mode under FR-004.
- **Risk: GEOM demand arrives late** (R-07) and forces EB plus C++. Q4 is answered for now (D-033: GEOM/CC_IBM deferred, not dropped); the risk returns when GEOM is brought back. Mitigation: keep Option D as the path.

## Open questions for teammates
- **FDS Legacy Mapper:** deliver `obst_wall_cface_indexing.csv` (pending): every index from `WALL`/`THIN_WALL`/`CFACE` into `BOUNDARY_*` storage, and every place that assumes wall indices are stable within a run. Which parts of `REASSIGN_WALL_CELLS` are reusable for a regrid rebuild?
- **AMReX Integration Lead:** cost of rebuilding per-box wall lists in `LayoutData` at regrid; `MLCellABecLap` overset-mask behaviour across levels (E-2); can EB2 be restricted to some levels or avoid coarsening thin features?
- **AMR Spec & Program Lead:** *answered* (v0.2): FR-040 snap rule and M1-M4 (D-009); Phase 5 freeze across all stateful walls (D-010). *Answered* (v0.2.1): Q4 scope by D-033 (thin OBSTs required; GEOM/CC_IBM and HT3D deferred). Open: FR-041b milestone once G2 reports; align requirements.md FR-040 R3 (and its "Open point R3") and FR-041a "Rule" with the v0.2.1 R3-T amendment and D-010 clarification.
- **AMR Pressure Solver Lead:** confirm E-1 → E-2 sequencing and per-zone gauge under E-2 with level-dependent masks (spec §E, zones separated by solids). **Pending for R3-T (v0.2.1):** (P-a) β on covered coarse faces and on interface-band faces under E-2; (P-b) whether a masked coarse/fine interpolant is needed next to thin faces in the band.
- **AMR V&V Lead:** OBST-heavy anchor cases for FR-022/FR-041 (a pyrolysing slab crossing a refinement boundary; a burn-away case); energy budget through a regrid.
- **GNU/Intel build chiefs:** cost of enabling `AMReX_EB` (default OFF) on both toolchains, even if unused in Phase 5.

## Spike plan
| Spike | Scope | Pass / fail signal |
|---|---|---|
| **G1 Mask consistency** (3 days, no solver) | Take OBST inputs from the verification suite; compute per-level masks with FDS `NINT` snapping vs "finest-level extent + >50%" rule for r=2,4 | count of mismatches and `THIN` topology changes per case; confirms D-009 unless the no-mismatch rule on interface faces cannot be met (then the >50 % rule) |
| **G2 Wall-state regrid** (2 wk) | Static 2-level box with a heated OBST; refine then coarsen across its face; transfer `BOUNDARY_ONE_D` via OBST-patch keying | solid energy conserved to round-off across regrid; surface T continuous; cost estimate that sets the FR-041b target milestone (D-010) |
| **G3 Create/remove on levels** (1 wk) | Burn-away OBST inside a refined patch using mask flip + wall reassignment | matches single-level FireX baseline (T2) and FR-042 |
| **G4 E-1 vs E-2 solid error** (with Pressure Lead) | `duct_flow`-type case on 2 levels | wall normal velocity ≤ `VELOCITY_TOLERANCE` (E-1) / ≈0 (E-2); FR-035 |
| **G5 EB feasibility probe** (3 days, optional) | Build EB2 `IF_Box` union from a thin-wall OBST set; coarsen to 2–3 levels | whether multi-valued-cell errors occur; informs D only |

## Decision needed from the project owner
**Q4 (answered, D-033, the project owner, 2026-09-25):** thin OBSTs, HVAC, pressure ZONEs and level-set wildfire must work with refinement; GEOM/CC_IBM and HT3D are deferred, not permanently dropped. Under the original framing ("No" confirms B with D deferred; "yes" adds EB (C++), a GEOM→EB2 path and likely a large extra phase, R-07), this confirms Option B with D deferred. Nothing further is needed from the project owner for this ADR at present.

## References
FireX `36975d765f` `Source/{type,read,func,init,main,wall,velo,pres,geom,ccib}.f90`; `Manuals/FDS_User_Guide/FDS_User_Guide.tex` (1685, 9455); `Manuals/FDS_Technical_Reference_Guide/samr_notes.tex:50-110`. AMReX `Src/EB/*`, `Src/LinearSolvers/MLMG/AMReX_MLCellABecLap.H`, `Docs/sphinx_documentation/source/EB.rst`. `ref/ERF/Source/EB/`. Teammate: `docs/amrex/mapping.md` §5, `docs/pressure/01-amr-mapping-spec.md` §D–E, `docs/{requirements,risks,charter}.md`.
