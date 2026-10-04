# Radiation on the AMReX Hierarchy: Spec Outline

**Status:** Draft outline (not reviewed, not approved). Rev 2026-09-26: owner decision on Q12, S-B default (see §2).
**Author role:** Radiation Lead
**Depends on:** requirements.md (FR-005, FR-006, FR-022, FR-040, FR-041, FR-060, FR-061, NFR-030/031, NFR-043/044, NFR-047), ADR-001, ADR-002 (leaning), ADR-003 v0.2.1, risks.md (R-16, R-36), amrex/mapping.md:287
**Source pin:** FireX `36975d765f` (local branch `AMReX`, this repository (read-only)). Citations get re-pinned per milestone (D-034). *Re-pin note:* the `radi.f90` citations in this file now follow the local branch after the upstream merge of `afb5e31a48` (RTE_SOURCE precompute): lines up to 3797 are unchanged, and lines after 3797 moved by +1 to +8 inside the band loop and by +4 from the end of the band loop onwards (file length 5,323 to 5,327). `main.f90` and `wall.f90` citations were also moved to the local branch at `bee11f0329` (main.f90 is +73 lines from `main.f90:~1000` on, wall.f90 +6); `mesh.f90`, `dump.f90`, `func.f90`, `type.f90`, `read.f90` are unchanged between the two pins. See `02-radiation-gpu-candidates.md` Part B.
**Legend:** [REC] recommended proposal · [ALT] alternative · [OPEN] undecided · [VERIFY] not yet confirmed in source. Every decision in this document is a proposal.

## 1. Current solver (FireX)

**Structure and settings**
- radi.f90 has 5,327 lines. The RAD module spans :2754-5326 and MIEV :331-2749. RADCAL lives separately in rcal.f90.
- The namelist is READ_RADI (read.f90:10206-10371).
- Defaults: NUMBER_RADIATION_ANGLES=100 (104 after the angle set is built), TIME_STEP_INCREMENT=3, and ANGLE_INCREMENT=MAX(1,MIN(5,NRA/15))=5 for gray. The 2D defaults are 60 angles and TSI=2.
- The wide-band model uses NSB=6 and WSGG uses NSB=5. Both force ANGLE_INCREMENT=1.
- INIT_RADIATION (radi.f90:2788-3296) builds the RADCAL κ table: 44 temperatures × 50 concentrations (:2770-2775). The table does not depend on grid spacing, because PATH_LENGTH defaults to 0.1 m (:3136). The same κ can therefore be used on every level.
- Angles are set up in CALCULATE_FVM_ANGLES (:3301-3406) and CALCULATE_DIRECTION_COEFFICIENTS (:3431-3596). The optional random rotation is drawn on rank 0 and broadcast with MPI_BCAST (:3443-3447).
- INTERPOLATE_IL (:3599-3717) remaps stored boundary intensities to the rotated angles.

**Per step (RADIATION_FVM, :3772-5069)**
- An update happens when MOD(RAD_CALL_COUNTER,TSI)==0, during initialization, at ICYC==1, or when UPDATE_ALL_ANGLES is set (:3856-3862). The counter is per mesh (mesh.f90:339).
- Each update sweeps only the angle subset N = NRA-AIC+1, stepping down by ANGLE_INCREMENT (:4315-4321). That is about 21 of 104 angles, so a full angular cycle takes 15 steps at defaults.
- UII is the sum of UIIDIM slots, and those slots have different ages (:4299-4303, :4961-4976).
- κ comes from GET_KAPPA (:5174-5212). Particles add KAPPA_PART and KFST4_PART in the cell that contains them (:3980-3992).
- Gray RTE source correction: the partial sums RAD_Q_SUM/KFST4_SUM (:4088-4125) are combined by MPI_ALLREDUCE at main.f90:1833-1834, then damped and clipped at main.f90:1840.
- The 3D sweep is a hyperplane (i+j+k) wavefront with `!$OMP PARALLEL DO` on each plane (:4494-4718). Solid cells use the upwind override from CELL_ILW (:4542-4548).
- Differencing: STEP is at :4621-4628. The FireX-only DIAMOND and EXPONENTIAL schemes are at :4635-4714 (commit 7d8dcb2707). No Verification input uses them.
- RTE source (upstream `afb5e31a48`): the angle-independent term `RTE_SOURCE = KFST4_GAS + KFST4_PART + RSA_RAT*(SCAEFF+SCAEFF_G)*UIIOLD` is built once per band (:4220-4222, freed at :4942) and read by the cylindrical, 2D and 3D sweeps (:4465, :4489, :4628, :4672, :4708) instead of being recomputed inline for every angle. Same operation order as the inline form, whole-array on `(0:IBP1,...)` like `KFST4_GAS`, ghost values never read by the sweeps. It is a box-local elementwise kernel computed once per band per box before the sweeps, with no new obstacle for FR-062 (review: `02-radiation-gpu-candidates.md`, Part B).
- Walls: WALL_LOOP1 sets the incoming boundary intensity (:4327-4372). WALL_LOOP2 stores ILW per wall, per angle, per band and updates INRAD_W (:4747-4775). Q_RAD_IN is set at :4902-4920.
- QR = κ·UII − KFST4 (:4981-4988). QR enters the divergence at divg.f90:567/579.
- Memory: intensity is not stored for all angles per cell. IL is one scratch array per angle. Per-angle storage exists only for wall cells (NRA×NSB per wall) and mesh-interface cells (IL_S/IL_R/IL_R_OLD, type.f90:1035-1039).

**How sweeps cross meshes and MPI today**
- Each mesh sweeps on its own. At an INTERPOLATED boundary, the ghost intensity is the arithmetic mean of the neighbour-mesh IL_R over NIC_MIN..NIC_MAX (:4353-4361).
  - At 2:1 in 3D, a coarse ghost averages 4 fine cells and a fine ghost copies 1 coarse cell (NIC setup main.f90:2166-2233).
- After each angle, IL_S is packed for all other meshes (:4826-4843). The loop skips NM==NOM, so periodic self-coupling may never be filled [VERIFY].
- MESH_EXCHANGE code 2 packs only angles with DLN(IOR,N)>0 from the current subset (main.f90:3559-3583). It uses persistent MPI_STARTALL on REQ5 (:3671-3674) and unpacks at :3815-3834. Same-rank meshes copy directly (:3507-3508).
- With the default RADIATION_ITERATIONS=1, the exchange runs at the end of the step (main.f90:1183-1189). Intensity therefore lags by one angular cycle per mesh crossed. This is block-Jacobi in time and depends on the decomposition. The User Guide documents it (FDS_User_Guide.tex:5894).
- RADIATION_ITERATIONS>1 repeats the whole solve K times per step, with a MESH_EXCHANGE(2) after each pass (main.f90:1095-1117). Each pass re-sweeps the same angle subset, because RAD_CALL_COUNTER advances only on the last pass (radi.f90:3866). So each extra pass carries intensity across one more mesh interface within the step, and costs a full extra radiation solve. After the first cycle only one exchange runs per pass, whatever ANGLE_INCREMENT is (`IF (ICYC>1) EXIT`, main.f90:1111-1114). Exchange happens only on intensity-update steps (EXCHANGE_RADIATION, radi.f90:3856-3862). INITIAL_RADIATION_ITERATIONS defaults to 3 (read.f90:10224).
- Creating or removing an OBST sets UPDATE_ALL_ANGLES (main.f90:1868).
- Restart writes UIID (dump.f90:3921) and RAD_Q_SUM/KFST4_SUM/RTE_SOURCE_CORRECTION_FACTOR (dump.f90:3950). Wall ILW is packed by PACK_BOUNDARY_RADIA (func.f90:5128-5146).

**Interaction with in-scope features**
- **Thin OBSTs** (zero thickness, ordinary WALL cells on both sides, per ADR-003) block radiation through the CELL_ILW override (:4362-4365, :4545-4547). Behaviour on box faces and C/F faces is [VERIFY].
- **HT3D thin walls** are deferred (wall.f90:479-496).
- **Level set:** LEVEL_SET_MODE 1-3 turn radiation off. Modes 4 and 5 keep it on (read.f90:1986-2016), and all level-set modes set NO_PRESSURE_ZONES (:1983).
  - Boundary fuel absorbs through the wall Q_RAD_IN.
  - Vegetation particles absorb and emit through KAPPA_PART/KFST4_PART (:3980-3992), so particle deposition must go to the owning level.
- **Pressure zones:** only WSGG reads PBAR(K,PRESSURE_ZONE) (:4066). No Verification case uses WSGG.
- **HVAC:** no direct radiation coupling [VERIFY].
- **Output quantities affected** (Smokeview and VTK, data.f90):
  - Cell quantities: ABSORPTION COEFFICIENT :146, INTEGRATED INTENSITY :154, RADIATION LOSS :158, RADIATION EMISSION :437.
  - Wall and device quantities: RADIATIVE/NET/GAUGE/INCIDENT HEAT FLUX :1364/1410/1424/1451, RADIANCE :1445, RADIOMETER :1459.
  - HRR column Q_RADI (dump.f90:863).
  - Per-mesh angle-resolved RADF files (radi.f90:5046-5067) have no level rule yet [OPEN].
- **Parallelism today:** only OpenMP; radi.f90 has no GPU directives.
- **Global state:** the counters, UPDATE_ALL_ANGLES/EXCHANGE_RADIATION and WEIGH_CYL are host-global variables written at run time (inventory module_globals.csv). They block pure device kernels.

## 2. Sweeps on a load-balanced multi-level hierarchy

**Box coupling within one level** (owner decision Q12, 2026-09-26; FR-005(i) exemption for the radiation stage is being recorded by the Spec Lead; FR-005(iv) run-to-run reproducibility still holds):
- **S-B [DECIDED, default]: FDS-style lagged box-face exchange with RADIATION_ITERATIONS.** Each box sweeps its angle subset independently, using box-face ghost intensities from the previous exchange. A per-box face buffer replaces OMESH IL_S/IL_R. The exchange is an AMReX FillBoundary-style copy of the upwind face intensities per angle, and it follows the FDS cadence (main.f90:1095-1117, 1110-1116): after each pass when K>1, otherwise at the end of the step, and only on intensity-update steps.
  - The iteration loop keeps FDS semantics. K passes per step move intensity across K box faces within the step, each pass costing a full solve.
  - Reproducibility: the sweep has no atomics and the exchange is fixed, so results are run-to-run identical (FR-005(iv)). With the per-box default below they change only with the box split, which D-039 exempts.
- **How the lag differs from FDS meshes.** In FDS, intensity loses one exchange per mesh interface crossed, and users choose few, large meshes. AMReX boxes are much smaller (floor 16, D-024; typical max_grid_size 32-64), so a ray crosses several times more interfaces over the same distance. With K=1 information needs about one update cycle per box crossed. A 64-cell path at 16³ boxes crosses about 4 faces, against 0-1 for a typical FDS mesh layout. Regrid and load balancing also move the faces, so the lag pattern shifts over time, and new faces start with no lagged value.
  - **Default [DECIDED]: fully parallel per-box sweep.** Every box sweeps independently from lagged box-face intensities, like FDS meshes. There is no ordering between boxes, on the same rank or not. Faces are double-buffered: a box reads only face values from the previous exchange, never ones written in the current pass, and one uniform exchange serves same-rank and cross-rank faces (FR-062 (a), ruling in adr/drafts/rulings-IR008-FR062.md §2). Byte identity across ranks and threads also needs particle κ and source deposition to stay deterministic in the GPU port (§4). At a fixed box layout the result is therefore independent of rank count, ranks per GPU, and which rank owns which box (D-039 exempts box-split dependence only). RADIATION_ITERATIONS tightens the coupling: K passes carry intensity across K faces within a step.
  - New or moved faces after a regrid or rebalance are filled from the box interior on the upwind side, or from the coarse level where it covers them, never left at zero (FR-062).
- **Cost argument (estimates, to be measured).** Total sweep work is the same for any ordering: cells × angles in the current subset. What differs is the dependency chain, which dominates on GPU because each hyperplane is a dependent step (§4).
  - Per-box, K=1: one 16³ box has 46 planes per angle. All boxes and all angles in the subset run the same plane index together, so the chain is 46 steps per pass. It is the same on 1 rank or 100.
  - Ordered sweep over an n×n×n block of boxes: the chain grows to about 3·16·n planes, which is 190 for n=4. That is about 4× the latency of per-box K=1 for the same work, and it also idles the GPU at the start and end of each sweep.
  - Per-box, K=2 or 3: work is K× (a full extra solve per pass), and the chain is K×46 (92 or 138). In the latency-bound GPU regime that is still shorter than a 4³-box ordered sweep. On CPU, where work dominates, K=2 costs roughly 2× radiation time.
  - So per-box plus iterations is the cheaper route on GPU for comparable accuracy, and it is the default. The ordered sweep is justified only if K=1..3 cannot meet the lag-error tolerance.
- **Iteration default [REC]:** RADIATION_ITERATIONS=1 and INITIAL_RADIATION_ITERATIONS=3, matching FDS (read.f90:10224). The FR-062 lag-error check (`radiation_gas_panel` split into 16³ and 32³ boxes at K=1,2,3, against single mesh) sets the recommended K. If K=1 is not accurate enough, the options in order are: raise K; then, as an optional rank-independent mode, sweep radiation on a separate BoxArray with larger boxes (fewer lagged faces, at the cost of a ParallelCopy of κ, sources and QR each radiation step); and only last, the ordered sweep below.
- **S-A [ALT, conditional]: dependency-ordered sweep across boxes.** Build it only if the FR-062 lag-error check shows per-box plus K=1..3 is not accurate enough. It is not on the M2 critical path. If built, it must order over all boxes (the global upwind order), never rank-locally, so that results stay rank-count independent. The single-mesh FDS run is the lag-error reference, so the check does not need S-A.
- Both modes share one per-box sweep kernel. The driver chooses only the ordering and exchange.

**Coupling between levels:**
- **L-1 [REC, Phase 8a increment; the R-16 mitigation]: one radiation level** (default level 0, configurable to a fixed L_rad).
  - Before the sweep, finer data is averaged down conservatively: κ, KFST4 as κ·4σT⁴ (never T), and particle terms.
  - Afterwards, QR goes to finer levels by conservative prolongation. Fine emission stays local: QR_f = κ_f·U_c − KFST4_f.
  - Fine walls take incident flux from the covering coarse wall ILW or from a wall-only reconstruction.
  - This option is cheap, but it loses the resolution of the flame T⁴ peaks and of walls that exist only on fine levels.
- **L-2 [REC, Phase 8 target default]: per-level sweeps, coarse to fine.**
  - The coarse sweep uses averaged-down sources in covered cells.
  - For each angle, the fine level then sweeps with C/F ghost intensities prolonged from the coarse IL (the same one-way rule as the FDS fine ghost).
  - Interleaving per angle avoids storing all angles.
  - Composite QR/UII is the finest available. Each level's walls get their own Q_RAD_IN.
- **L-3 [ALT, rejected for now]: composite multi-level sweep with two-way coupling in one pass.** The DAG over all levels is complex and fits GPUs poorly.
- **L-4 [rejected]: P1 or MLMG radiation** (as in PeleMP, which uses composite MLMG). It changes the physics, which D-003 does not allow.

**Recommendation:** use the per-box S-B sweep on every level plus L-2, with L-1 kept as the first increment and as a user option (FR-060 already allows "level 0 only / fixed level / all").
- The RTE correction sums run over uncovered cells only, with the exact sum of FR-005(ii).
- Q_RADI and the FR-022 budget integrate over uncovered cells only.

## 3. Resolution policy and time-sharing (ADR-002)

- Angle set: the same NRA, rotation and ANGLE_INCREMENT on all levels [REC]. C/F exchange then needs no angular interpolation.
  - Angular refinement on fine levels [ALT, later] could reuse the INTERPOLATE_IL 4-neighbour weights (:3599-3717).
- Spatial: each level sweeps at its own dx under L-2. L-1 uses L_rad.
- Under the global dt with no subcycling (ADR-002 leaning), keep FDS cadence. Make RAD_CALL_COUNTER and ANGLE_INC_COUNTER global instead of per mesh, and have all levels sweep the same angle subset in each step.
- Under subcycling (only if D-007 pulls it in): update radiation on level-0 steps only, hold QR through fine substeps, and count level-0 steps.
- Regrid:
  - Fill all UIID slots on new fine cells from the coarse UIID.
  - Copy wall ILW from the covering coarse wall, or set it to OUTRAD, for FR-041.
  - Do not force UPDATE_ALL_ANGLES on every regrid (about 5× the cost). Whether to do so is [OPEN].
- A new random rotation needs INTERPOLATE_IL to also remap C/F and box-boundary buffers.

## 4. GPU execution (estimates)

- Kernel: hyperplane wavefront over precomputed cell lists. This can be written in either K1 (ParallelFor over a list) or K2 (`!$omp target teams loop`). It needs no atomics. UIID accumulation per angle is deterministic.
- Occupancy: one 16³ box has 46 planes per angle, which is too small. The plan is to batch (box, angle) pairs per plane index across all boxes and all angles of the current subset (about 21 angles in gray).
- Memory for angle-parallel sweeps: cells × angles in flight × 8 B. That is about 170 MB per 1M cells at 21 angles, and 832 MB per band for all 104 angles (wide band). Process band by band to stay within NFR-031.
- Per-box S-B (the default) has no serialization between boxes: all (box, angle) pairs at the same plane index run together. The conditional S-A would serialize along the sweep diagonal (cost argument in §2).
- Wall loops are unstructured per-wall kernels. RADCAL and table setup stay on the host at init. Particle κ deposition: current FDS has no OpenMP atomics in part.f90 or radi.f90 (grep, source pin). A GPU port must keep the scatter deterministic (for example a sorted segmented sum, not floating-point atomics) to preserve FR-062 byte identity.
- Several ranks per GPU (NFR-047) do not change per-box results or the lag, which depends only on the box layout.

## 5. Cost relative to the flow solve (estimates; TO MEASURE)

- Sweep work per step is about cells × NRA/(ANGLE_INCREMENT·TSI), roughly 7 angle-sweeps per cell per step in gray. Wide band costs about 30× more (104 angles × 6 bands / 3 steps).
- Measured on the development machine (GNU FireX baseline, RADI share of T_USED):
  - 2D cold cases: 1.6-10.5% (dancing_eddies_1mesh 10.5%).
  - 3D fm_burner 2 cm: 12.7%, but that run is dominated by initialization.
  - The User Guide cites about 20% (FDS_User_Guide.tex:218).
- Under L-2 with a global dt, radiation cost scales with total cells, like the flow. Expect a 10-20% share on CPU, with no pipeline penalty under S-B at K=1; each extra radiation iteration adds a full radiation solve (K=2 roughly doubles radiation time).
- On GPU the share will probably grow, because sweeps are latency- and dependency-bound.
- A 3D steady-fire measurement is needed (for example, the NFR-032 Heskestad case).

## 6. Verification cases the AMR code must run (FR-006) — pending coordination with the V&V Lead

- `plate_view_factor_cart_30` (anchor) and `_60`/`_100`: view factors; refined variant for FR-060.
- `radiating_polygon_square_20` (anchor) and `_40`/`_80`: flux from a hot polygon; refined variant for FR-060.
- `plate_view_factor_2D_30/60/100`: the 2D angle set.
- `radiation_box__20_*` and `radiation_box_100_*`: participating media, κ sweep; large-angle and large-cell cost.
- `radiation_plane_layer_1..6_1..5`: slab layer solutions (KAPPA0); 12 cases use the wide band.
- `emissivity`, `TC_heating`, `TC_view_factor`, `thermocouples`, `wall_internal_radiation`: wall and device coupling.
- `adiabatic_surface_temperature`, `net_rad_1000_obst`, `radiation_shield`: thin OBSTs (plus a ZONE in radiation_shield).
- `droplet_absorption_cart`, `part_attenuation`, `Sprinklers_and_Sprays/particle_(an)isotropic_radi`: particle κ and scattering; RADIATION_ITERATIONS>1.
- `radiation_gas_panel`: the only same-resolution multi-mesh interface-crossing case; used to measure S-B lag error against the single-mesh result (§2).
- `check_kappa`, `hot_spheres`, `target_test`: multi-mesh runs (completion or T2).
- `WUI/ground_vegetation_radi`, `radiation_gas-veg_consistency_gas/_veg`, `level_set_fuel_model_1`: vegetation particles and level-set mode 4.
- `Energy_Budget/*`, `Heat_Transfer/adiabatic_net_flux`, `Fires/ftp`, `Detectors/beam_detector`: energy closure (FR-022).
- `Thread_Check/race_test_1_r4`, `race_test_4_r4` (A-35 4:1 copies): the only fire + radiation + C/F candidates.
- Derived cases to create:
  - Static 2:1 patch variants of the two anchors, with the patch between emitter and target.
  - `radiation_gas_panel` split into multiple boxes and a refined variant.
  - A thin OBST at a C/F face.
  - An energy-budget case with a refined fire patch.
  - Cases for RANDOMIZE_RADIATION_DIRECTIONS and DIFFERENCING_SCHEME, which have zero coverage today.
- Deferred or unclear: GEOM cases (`geom_rad*`, `net_rad_1000_cgeom`, `plate_view_factor_ibm_*`) are deferred. Cylindrical cases (`plate_view_factor_cyl_*`, `droplet_absorption_cyl`) are UNCLEAR (Q4).

## 7. Open questions and risks

**Open questions**
- Q1 [DECIDED 2026-09-26]: S-B lagged exchange with RADIATION_ITERATIONS is the default; FR-005(i) exemption for the radiation stage. Remaining [OPEN]: the acceptance tolerance for S-B lag error versus single-mesh, and a minimum radiation box size if needed.
- Q2 [OPEN]: Is L-2 the default, with L-1 as the first increment? The ADR owner is the Chief Architect (FR-060).
- Q3 [OPEN]: Which level owns a wall's Q_RAD_IN, and what are the rules for thin OBSTs on box faces and C/F faces (ADR-003 R3-T)?
- Q4 [OPEN]: What happens at regrid: a full-angle update or UIID/ILW prolongation? Does ILW become part of the FR-041 wall state?
- Q5 [OPEN]: What is the RADF per-mesh output rule under AMR? Which radiation slice and boundary quantities go to refined-level Smokeview and VTK output?
- Q6 [OPEN]: Do FireX DIAMOND/EXPONENTIAL and random rotation stay in scope without Verification coverage (D-003)?
- Q7 [VERIFY]: Periodic self-coupling of IL_S (radi.f90:4826-4843).

**Risks (proposed)**
- (Retired with Q1 decision: S-A GPU and many-rank scaling applies only if the optional exact mode is built.)
- S-B lag grows with the number of lagged faces a ray crosses (small boxes, regrid moving faces), which degrades accuracy relative to FDS multi-mesh and makes results shift with load balancing. Mitigation: face refill at regrid, the FR-062 lag-error check, RADIATION_ITERATIONS>1, then an optional larger-box radiation BoxArray; the ordered sweep last.
- Angle-parallel GPU memory (cells × angles × bands) could exceed NFR-031.
- Regrid loses the time-shared UIID/ILW history and causes flux transients at new fine patches.
- Radiation host-global counters and flags block pure device kernels (R-35/R-46 class).
- R-16 as written ("level-0 plus QR interpolation") under-resolves T⁴ emission in refined flames.
