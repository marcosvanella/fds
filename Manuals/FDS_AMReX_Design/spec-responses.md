# Spec & Program Lead responses to the ADR drafts

Owner: Spec & Program Lead · Status: draft v1.4 (2026-09-25; changelog in README.md v0.4.17) · Source pin: FireX 36975d765f on branch `FDS-AMReX` (this repository; renamed from `AMReX`, D-037)

This file answers the questions that `docs/adr/ADR-001..003` put to the Spec & Program Lead. The ADR files themselves are not edited. Labels:
- **Proposed**: my recommendation, not yet agreed by the team.
- **Pending the project owner's decision**: an owner decision (charter §9). My recommendation is given, but it is not a decision.

Each answer points to the doc changes that carry it: requirements.md, roadmap.md, and the README decision/action logs.

---

## ADR-001 (driver architecture), "Open questions: AMR Spec & Program Lead"

### (a) Minimum feature set for the Phase 2 demo (Accepted by Chief Architect 2026-09-25; D-006)

The demo is the **M2a** gate inside Phase 2. It extends ADR-001 spike S1 (single level, periodic, no walls; shim binds `RHO, ZZ, U, V, W, TMP`; runs `DENSITY`, `COMPUTE_VISCOSITY`) from two kernels to **full FDS time steps** on three existing inputs. That is the smallest set that tests FR-001 (single box) and FR-002 (multi-box, one level), in 2-D and 3-D, without walls.

**Demo anchor cases.** All exist in FireX `Verification/` and are in `FDS_Cases.sh` (lines 784, 801, 849).

| Case | What it proves | Features exercised | Reference / tolerance |
|---|---|---|---|
| `Scalar_Analytical_Solution/shunn3_32.fds` | FR-001, first case (S1's natural successor) | 2-D (J=1), 1 mesh, all-periodic, two species with variable density (ρ₀=5, ρ₁=1), DNS-like transport, `RADIATION=.FALSE.`, `CHECK_POISSON`, `MASS_FILE` | Bitwise (T0/T1 per Q1) on field dumps of the explicit kernels (mass, velocity predictor, divergence) on frozen input vs FireX baseline; whole run T2 vs FireX baseline (CSV and field dumps of `RHO, ZZ, U, V, W, TMP, H` at `T_END`), since bit-parity ends at the first pressure solve (D-022; field-dump tool TBD(V&V)) |
| `Turbulence/csmag_32.fds` | FR-001 in 3-D | 32³, 1 mesh, all-periodic, constant Smagorinsky LES, `UVW_FILE` initial field, `RADIATION=.FALSE.`, `SPATIAL_STATISTIC='MEAN'` DEVC | Bitwise explicit-kernel checks on frozen input vs FireX baseline; whole run T2 (D-022) |
| `Scalar_Analytical_Solution/shunn3_4mesh_32.fds` | FR-002 | Same problem as `shunn3_32` on 4 periodic 16×16 meshes, 4 ranks; not input-identical (`shunn3_32` differs in `&PRES` tolerance and `FISHPAK_BC`, `CONDUCTIVITY`, `MMS_TIMER`) | Bitwise (T0/T1 per Q1) explicit-kernel checks on frozen input vs a derived single-mesh copy run on the FireX baseline (A-24, D-022); full steps within tolerance (T2) vs that single-mesh copy and vs 4-mesh FireX baseline (FFT and `SOLVER='GLMAT'`); `H` from one frozen solve within eps_H vs GLMAT (FR-002) |

**In the demo**
- `USE_AMREX` CMake option. gfortran 14.2 + OpenMPI 5 only (oneAPI was restored on 2026-09-25, A-18 closed, but the M2a demo stays gfortran-only).
- C++ AmrCore driver with one level. One box per `&MESH` (no re-boxing); `max_grid_size` = mesh size.
- Shim (`POINT_TO_BOX`) binds the S1 fields plus every other field the full step touches for these three cases. The Legacy Mapper derives the list from the kernel footprints (`inventory/kernel_footprint_mass.*` and successors).
- Whole-step kernels through the shim: mass/species transport, viscosity and constant-Smagorinsky turbulence, velocity flux, divergence, predictor/corrector velocity update, stability/dt control, and `UVW_FILE` initialisation.
- Ghost fill with AMReX `FillBoundary`: periodic plus same-level box-box. This replaces `MESH_EXCHANGE` for these fields.
- Pressure outside the shim (ADR-001 decision 1):
  - single box: the FFT fast path behind the composite-solver interface (FR-037: `amrex::FFT::Poisson`, decided, `pois.f90` not ported, D-021; FFTW installed and AMReX configured with `AMReX_FFT=ON`, A-23 closed);
  - multi-box: a single-level composite MLMG solve (= prototype P2 / ADR-001 S2 / ADR-002 S2).
- Output: `CHID_devc.csv`, `CHID_hrr.csv`, `CHID_mass.csv`, the `.out` file, and field dumps for comparison. Global CFL/dt selection.
- `IR-006`: with `USE_AMREX=OFF`, the three cases are T0 vs baseline.
- Runs at 1 rank (single mesh) and 4 ranks (`shunn3_4mesh_32`), `OMP_NUM_THREADS=1`.
- NFR-030 measured (not required), including the pressure share of the step (A-22).

**Out of the demo** (moved to the M2 exit or later)
- OBST, VENT other than `PERIODIC`, OPEN boundaries and wall cells. `ns2d_16` moves to M2: it does not switch radiation off, so the radiation step runs by default.
- Combustion, radiation, particles, HVAC, pressure zones beyond one, restart.
- Smokeview/VTK output checks. Slice lines may write files, but they are not a demo criterion.
- Refinement and regridding, OpenMP > 1 thread, the oneAPI build, `MPI_PROCESS` re-mapping.

**Relation to the Phase 2 exit (M2).** M2 keeps the full uniform-mode list: FR-001/002 on all single- and multi-mesh anchors (walls, radiation, particles in uniform mode), FR-038, FR-061, IR-001, NFR-010/011/012, NFR-020, NFR-022. The demo is a proper subset reached first. If the demo fails, the shim approach is re-examined against the ADR-001 S1 "overturns if" criteria before the rest of Phase 2 starts.

### (b) R-26 exit milestone in the roadmap (done)

It is in `roadmap.md`, in two places:
- Phase 1 exit criteria: "If ADR-001 adopts the per-box `POINT_TO_MESH` shim, it states the shim's exit condition (R-26)…".
- The **M4** milestone row: "R-26 shim-exit review (kernel-extraction date set, or regrid rebuild ≤ 10% of step time) if ADR-001 adopts the shim".

v0.2.1 also adds the review to the Phase 4 exit criteria and records the Phase 2 regrid-cost measurement from prototypes P1/P3 as its input. The thresholds are ADR-001's: 10% of step time, and a kernel-extraction date by M4.

### (c) Stretched grids and Smokeview as hard requirements (Answered by the project owner: stretched, D-030; Smokeview, FR-072)

**Owner decision on Smokeview (Q5, the project owner, 2026-09-25):** AMR runs must write both Smokeview-format output and VTK; Smokeview cannot read VTK, so VTK does not replace it, and Smokeview output from refined levels is required, not just a fixed-mesh view. This overrides the Smokeview recommendation below (kept as history). The writing approach is chosen by an ADR after a small-case test (A-40, R-43). Correction to the reasons below: FDS writes the `GRID` list once at setup (dump.f90:2457-2541); whether Smokeview itself needs a static list is unverified.

**Stretched grids (`&TRNX/Y/Z`): recommend "not a hard requirement with refinement".**
- AMReX levels are uniform per level.
- The ADR-001 shim synthesises uniform per-box metrics.
- Local refinement is the modern replacement for stretching.
- Usage is small: 17 of 933 Verification inputs and 49 of 3,132 Validation inputs contain a `&TRN` line (grep count; heuristic).
- `CYLINDRICAL`: same answer (11 Verification inputs).

Proposal: reject both when refinement is requested (FR-004, IR-002). They keep running wherever uniform mode is served by the legacy path. Note that `amrex::FFT::PoissonHybrid` supports stretching in z only, and z cannot be periodic (Pressure Lead, 2026-09-25; `AMReX_FFT_Poisson.H:163`). FDS x/y stretching would therefore need a new solver path even in uniform mode under AMReX, which is a further cost argument. Whether the `USE_AMREX=ON` executable keeps a legacy fallback driver for such inputs, or they need a `USE_AMREX=OFF` build, is TBD(Chief Architect).
Supporting evidence (Pressure Lead, 2026-09-25): stretched runs stay outside the AMReX driver, because AMReX `Geometry` is uniform and `PoissonHybrid` stretches in z only. **Decided (owner decision, 2026-09-25, relayed by the Chief Architect):** AMR mode uses only uniform grids on each level, and stretched cases stay FDS-only (D-030; charter Q11 (c)).

**Smokeview: recommend "hard requirement for uniform mode and for a fixed-mesh view of AMR runs; not for full-resolution refined-level output".**
- Uniform mode: formats unchanged (FR-072), since users and firebot depend on them.
- AMR mode: Smokeview output resampled onto the user's static `&MESH` list at each mesh's own resolution, so the `.smv` `GRID` list stays static. Levels finer than the output meshes are averaged down for Smokeview.
- Full-resolution hierarchy: FireX VTKHDF (`vtkf.f90`, FR-074) as the supported secondary path, read with ParaView. AMReX plotfiles are optional (FR-073).

Reasons:
- Smokeview needs a static `GRID` list (dump.f90:2499-2511), and changing Smokeview is a non-goal.
- FireX already ships a VTK writer whose `UnstructuredGrid` type is not tied to a static mesh list (unverified for AMR).
- Box-as-`GRID` output explodes file counts and changes every regrid (R-08).

---

## ADR-002 (time stepping), "Open questions: AMR Spec & Program Lead"

### (d) Cut Phase 6 to "sync corrections only, no subcycling"? (Accepted by Chief Architect 2026-09-25; D-007, Pending the project owner's decision for the final scope)

**Recommend: yes, as the default, but conditional on the ADR-002 S1 global-dt penalty spike, run on the NFR-032 benefit case (answer e) before M1.**

| S1 result: estimated work of global dt (A) / subcycled (B), 1–2 levels, ratio 2 | Phase 6 scope |
|---|---|
| ≤ 1.5× | **Sync corrections only.** Subcycling moves to "out of roadmap". NFR-032 is judged under global dt. |
| > 1.5× and ≤ 3× | Sync corrections first (6a). Subcycling (6b) is optional, with a go/no-go at M5 on the re-measured NFR-032 projection. The ADR-002 Option C "subcycling-ready" data model (per-level old/new state, time-weighted flux registers) is kept from Phase 3. |
| > 3× | **Subcycling pulled earlier.** Option C data model mandatory from Phase 3. Phase 6 (subcycling + sync corrections) starts at M4 in parallel with Phase 5 and is a precondition for NFR-032 at M10. |

Reasoning:
- With one global dt and a composite solve, the FDS step is algebraically a MAC projection with potential H: β = dt in the predictor (velo.f90:1603-1630) and dt/2 in the corrector (velo.f90:1723-1751). No synchronisation *projection* is needed (Pressure Lead, 2026-09-25).
- What remains before subcycling is coarse-fine flux/divergence consistency. Flux registers (FR-024, Phase 3) and composite zone integrals (FR-034) already cover it.
- So R-05 is a subcycling-only risk. Its gate moves from M4 to the Phase 6 subcycling decision (risks.md R-05; roadmap Phase 4/6).
- Subcycling is only worth its cost (it breaks four same-time global couplings, per ADR-002) if S1 shows a large global-dt penalty on a realistic case.

### (e) NFR-032 benefit case (Proposed; D-008; level-0 grid accepted as D-024)

**Primary: `Validation/Heskestad_Flame_Height/FDS_Input_Files/Qs=1_RI=10.fds` as the level-0 input, with one refinement level (ratio 2) around the burner and flame.**
- **Input:** single mesh, IJK 33×33×80 (87,120 cells), XB −1.8..1.8 × −1.8..1.8 × −0.45..8.59 m, `T_END=54` s. **Level 0 for AMR (D-024, accepted 2026-09-25):** IJK 32×32×80 (81,920 cells), high side trimmed to XB −1.8..1.690909 in x and y, same cell size, because AMReX needs the level-0 domain divisible by `blocking_factor` (`AMReX_AmrMesh.cpp:1252-1259`). The burner stays exactly 9×9 cells; a centred 32-cell domain would make it 10×10 (+23% burner area, hence a different diameter and Q*; total HRR is unchanged because FDS rescales HRRPUA to the input area, init.f90:3574-3575). The burner is a 1 m² OBST with `HRRPUA=1512.7` kW/m² (Q* = 1, ≈1.5 MW propane), `THICKEN=.TRUE.`, open on all sides.
- **Why it fits a global dt:** there are no other fast flows. The fastest gas velocities are in the flame and near plume above the burner, which is where the refinement goes. So the CFL-limiting cells should lie inside level 1. ADR-002 S1 confirms this with `&DUMP CFL_FILE=T` before M1; if the limiting cells are outside, the case is replaced.
- **Refined region (proposed):** |x|,|y| ≤ 0.9 m, z ≤ ~4 m, about 1.5× the expected flame height (≈2.7 m from Heskestad's correlation, L/D = 3.7 Q*^0.4 − 1.02; *estimate*). That is about 12% of the domain volume.
- **Cost estimate** (arithmetic, not measured): about 86k level-1 cells + 82k level-0 cells ≈ 26% of the uniform-fine cell count (64×64×160 = 655,360, D-024). Under a global dt that leaves margin to the ≤ 50% target for regrid and ghost-fill overhead. Uniform-fine memory is well inside 12 GB (≈ 1.6 GB, V&V estimate). Uniform-fine runtime on 8 ranks ≈ 2.3 h (range ≈ 0.5–4.5 h, V&V estimate, not measured); a median-of-3 timing needs ≈ 7 h of the whole idle machine (R-37, A-26). This replaces the earlier "tens of minutes" guess.
- **Metric for T3:** the case's own outputs, flame height `Lf_95/97/99` (running averages after 27 s, from the `HRRPUV` profile DEVC) and the HRR. Proposed T3: |Lf − Lf_fine| ≤ one fine cell height (0.057 m) and HRR within 1%. V&V's provisional T3 values (vv/test-plan.md §9, R-2) are now in requirements.md NFR-032 and §2.2.
- **Reference:** a derived uniform-fine input at exactly twice the resolution, IJK = 64,64,160 over the trimmed XB (D-024), effectively RI ≈ 20. V&V generated it in `(local V&V run directory)/inputs/A-19/` (`Qs1_RI10_fine_64x64x160.fds`, with level 0 `Qs1_RI10_coarse_32x32x80.fds`; A-19 closed); burner 9×9 coarse / 18×18 fine, 0.963967 m², 1512.7 kW on both. The existing `Qs=1_RI=20.fds` uses IJK 65,65,160 (ratio 1.97 in x/y), so it serves only as an informational cross-check.
- **Compatible with (g):** the burner is stateless (fixed `HRRPUA`, `INERT` sides; no in-depth conduction or burn-away). With `THICKEN`, it is one cell thick on each level, which is allowed mismatch M4 in (f).

**Fallback: `Validation/McCaffrey_Plume/FDS_Input_Files/McCaffrey_14_kW_5.fds` as level 0 (27 meshes of 17×17×18, 140,454 cells, `T_END=30` s, 0.3 m burner OBST at `HRRPUA=160`), refined over the burner and lower plume column.**
- Its metric is centreline temperature/velocity versus height.
- Its reference is a derived IJK 34,34,36-per-mesh input (≈1.12M cells). The native `_11` file (37 cells per mesh) is not ratio-2 nested.
- It is used if Heskestad's CFL-limiting cells fall outside the plausible refined region, or if the uniform-fine reference takes too long on the development machine.

All four input files were checked present in `Validation` (read-only).

---

## ADR-003 (geometry), "Open questions: AMR Spec & Program Lead"

### (f) FR-040 snap rule and allowed-mismatch list (Accepted by Chief Architect 2026-09-25; D-009; Spike G1 can overturn)

**Actual FDS rule** (FireX `READ_OBST`, per mesh; each step uses that mesh's grid through `GINV`, func.f90:6722):
1. Bounds snap to the **nearest cell face**: `OB%I1 = NINT(GINV(XB1-XS,1,NM)*RDXI)`, likewise I2..K2 (read.f90:11157-11162).
2. **Thin collapse:** without `THICKEN`, an OBST thinner than 0.25 cell in a direction whose snapped indices differ collapses to the nearer face, giving zero thickness (x: read.f90:11174-11181; y: 11184-11198; z: 11201-11215).
3. **Thicken:** with `THICKEN` (per OBST, default from `&MISC THICKEN_OBSTRUCTIONS`, read.f90:1815, 10890), a zero-thickness OBST becomes one cell thick around its midpoint (read.f90:11167-11173 and y/z analogues). This is not done if the OBST lies mostly outside the mesh (read.f90:11034-11039).
4. **Discard:** an OBST of zero thickness in two or more directions is thrown out on that mesh (read.f90:11218-11224).
5. **THIN flag:** a zero-thickness OBST whose input thickness is < 0.5 cell is flagged `THIN` (read.f90:11315-11317).

*Observation, confirmed (Legacy Mapper, 2026-09-25; A-20 closed):* in the thin-OBST collapse branch (`THICKEN_DIR` off), the y snap uses `GINV(XB3-XS,2,NM)`, x origin `XS` instead of `YS`, while the guard at read.f90:11192 and the other term use `YS` (read.f90:11193); the z snap compares `REAL(OB%I1)`/`REAL(OB%I2)` but assigns `OB%K1`/`OB%K2` (read.f90:11210). Both are also in `ce1f659` (read.f90:11060, 11077). They are upstream slips that can pick the wrong face for sub-cell OBSTs collapsed in y or z. Level 0 must reproduce them exactly (FR-040 R1, D-009). A third slip is confirmed on FireX at read.f90:11176 (an integer index subtracted from a length in metres); all three are logged in `inventory/uniform_grid_assumptions.csv` (Legacy Mapper, inventory complete).

**Proposed AMR rule**
- **R1, snap per level with the unchanged FDS rule:** each level is treated like an FDS mesh of that resolution, with level-0 origin and level dx. This makes level 0 identical to baseline (FR-001), and each fine level matches a baseline mesh of the same resolution (useful for FR-016-style comparisons). It is the rule ADR-003 Spike G1 compares against "finest-level extent + >50%". G1 overturns R1 only if its mismatch counts break R3 on the verification OBST set.
- **R2, covered regions: the fine level is authoritative.** Coarse cells under fine cover may disagree with the averaged-down fine mask. `average_down` averages gas children only; a coarse cell whose children are all solid is solid.
- **R3, coarse-fine interface: no mismatch allowed.** A face on the coarse-fine boundary that is a solid/gas face on one level must be one on the other. The regridder grows the fine region by one coarse cell around any OBST whose snapped extent differs between levels, so the boundary never cuts through a mismatch. If that is impossible, the run stops with an error (FR-004 style); it never runs silently.
- **R4, thinner than a coarse cell:** apply R1 on each level. On the coarse level, the OBST becomes zero-thickness (and `THIN`) if thinner than 0.25 cell, 0 or 1 cell by `NINT` otherwise, or disappears (step 4). On the fine level it keeps whatever thickness the finer grid resolves.

**Allowed mismatches** (only in covered regions, per R2/R3; listed per case by the FR-040 mask checker)

| ID | Mismatch | Source |
|---|---|---|
| M1 | Extent rounding: an OBST face differs by up to half a coarse cell (one fine cell at ratio 2) because `NINT` runs on different grids | step 1; ADR-003 illustration |
| M2 | Thin topology: zero-thickness/`THIN` on the coarse level, finite thickness on the fine level | steps 2, 5 |
| M3 | Vanishing: discarded on the coarse level (zero thickness in ≥ 2 directions), present on the fine level | step 4 |
| M4 | Thickened thickness: one coarse cell vs one fine cell | step 3 |

Anything else is an error. OBSTs aligned with coarse faces (the usual practice) snap identically on all levels at ratio 2 and produce no mismatch.

### (g) "No refinement boundary across burning OBST faces" for Phase 5? (Accepted for Phase 5; D-010; breadth in ADR-003 v0.2)

**Recommend: accept for Phase 5, generalised to a "stateful-wall coverage freeze".**
- **Stateful wall face:** a wall face carrying history that cannot be re-derived from the input and current gas state. This covers in-depth conduction (`BOUNDARY_ONE_D`), pyrolysis/accumulated mass, burn-away, ignition time from `IGNITION_TEMPERATURE`, and deposition.
- **Stateless faces** (fixed-temperature `INERT`, `HRRPUA` burners, OPEN, MIRROR, PERIODIC) can be rebuilt at any regrid.
- **Rule:**
  - A stateful face keeps the level that covers it after the initial grid is built.
  - Regrids may not refine or coarsen across it, and the coarse-fine boundary may not cut a stateful OBST face; tagging adds a one-coarse-cell buffer.
  - Wall state is therefore never transferred between levels in Phase 5.
- **Consequences for FR-041:** it splits into:
  - **FR-041a (Phase 5):** under the freeze, stateful wall data is bitwise unchanged across every regrid, and solid mass/enthalpy is conserved to round-off.
  - **FR-041b (deferred to after Phase 6; target TBD(Chief Architect) until ADR-003 Spike G2 gives a cost estimate):** real cross-level transfer by OBST-patch keying, per ADR-003 Spike G2.
- **Tests:**
  - a pyrolysing-slab case with regrids elsewhere in the domain (FR-041a bitwise);
  - a tagging test showing the fine region grows rather than cutting a stateful face;
  - a negative test where the freeze would force refinement beyond `MAX_LEVEL`/memory limits, which must stop with an error.
  - FR-022 (walls) runs on the same cases.
- **Risk (R-33):** the freeze blocks the most attractive AMR use case, refinement that follows flame spread over burning surfaces. It can also inflate the fine region around large burning objects. Breadth (AMR Chief Architect, 2026-09-25): any `SURF` with a `MATL_ID` carries 1-D conduction history from step 0, so nearly every non-inert wall is stateful before the first regrid, and the freeze means static refinement near all such walls, not just burning faces. This is acceptable for Phase 5 because the NFR-032 case (e) has a stateless burner. It must be lifted (FR-041b) before any claim about flame-spread cases.

### (h) Q4 scope: GEOM, HT3D, thin walls (Answered by the project owner, D-033)

**Owner decision (2026-09-25; D-033):** GEOM/CC_IBM and HT3D are out for now, deferred rather than permanently dropped (matches the recommendations below); thin obstructions must work with refinement (matches); HVAC, pressure ZONEs and level-set wildfire must also work with refinement (beyond this answer's scope); the AMReX code shall run the FDS Verification inputs that cover in-scope features (FR-006). The recommendation text below is kept as history.

- **GEOM/CC_IBM: recommend not required with refinement for "done".**
  - Reject with refinement (FR-044). Future path is ADR-003 Option D (EB for GEOM only, C++).
  - Cost is high: geom.f90 + ccib.f90 ≈ 51.8k lines, GEOM forces UGLMAT, EB needs C++.
  - Usage is significant (128 Verification and 366 Validation inputs contain `&GEOM`, grep count), so those inputs stay uniform-mode only. The project owner should weigh that before answering.
- **HT3D: recommend not required with refinement.** Reject with refinement (FR-004). 3-D solid conduction crossing boxes and levels is unassessed (ADR-003). It must still work in uniform mode (FR-003); whether the shim covers it is TBD(Legacy Mapper). Usage: 31 Verification, 12 Validation inputs.
- **Thin (zero-thickness) OBSTs / thin walls: recommend required with refinement (Phase 5).** They are core OBST semantics (read.f90:11174-11224, 11315-11317; `CELL_TYPE` `THIN_*` indices, type.f90:2175-2187), and excluding them would reject much of the rectilinear suite. They are supported under the FR-040 rules (M2/M3). Thin OBSTs with solid-phase state follow the (g) freeze.
