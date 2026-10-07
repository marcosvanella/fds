# Loop work list: device-eligible loops not yet translated

Owner: FDS Legacy Mapper (keeps this list; the claim protocol below is the only way to change an Owner cell).
Machine-generated tables sit between `GENERATED-BEGIN` and `GENERATED-END` markers. Everything outside the markers is hand-written. To refresh the tables: `python3 tools/inventory/loop_work_list.py` (see "Regenerating" at the end).

Companion files:
- `docs/amrex/loop_work_list.csv`: the same ranked list as data (one row per eligible loop, translated loops included).
- `docs/amrex/loop_claims.csv`: the claim register read by the script (key, owner, status, note).
- `docs/amrex/pressure-velocity-h-loops.md`: the detailed candidate list for the AMR Pressure Solver Lead.
- `docs/amrex/generator-howto.md`: how to add a loop to the generator and check it bitwise.
- `docs/amrex/blocked-loop-families.md`: the blocker families (P1 to P4, S1 to S4, SP1 to SP4, R1, V1, V2, O1 to O4) cited in the Blocker column.
- `docs/inventory/gpu_generator_loop_classes.csv`: the survey the list is built from (FireX 36975d7 line numbers).

## 1. What is on the list

The list covers every loop of the survey that is a device candidate and not yet in the generator at its committed state.

| Set | Loops | Share of total modelled time |
|---|---|---|
| Survey, non-geometry | 451 | 69.444 % |
| minus retired or replaced (`pois.f90` solver internals 205, exchange 2) | 207 | 18.332 % |
| minus host-side (set-up, I/O, MPI) | 13 | 1.935 % |
| **Translation-eligible (the denominator below)** | **231** | **49.177 %** |
| Geometry-deferred (cut-cell, CFACE, CCVAR), outside every denominator | 378 | 23.854 % |

Line numbers are those of the survey source (FireX 36975d7). The script checks that the first line of every listed loop is a `DO` statement in that source, so a line number cannot silently drift.

"Translated" means a kernel for the whole loop is in the committed generator sidecar (`amrex/s4_mass/s5_gen/s5_markers.toml`, `markers/pvf_kernels.toml`, `markers/gsfv_field_kernels.toml`) at the generator branch head the script reads with `git show`. Work that exists only as uncommitted edits in a shared worktree, or on another branch, is listed as in progress and is not counted.

### Coverage

- Committed at the generator branch head: **52 loops, 4.024 % of total modelled time (8.2 % of the eligible time)**. The survey CSV alone shows 45 loops (2.479 %); the sidecar adds PATCH_VELOCITY_FLUX L1390 (1.270 %) and the six GET_SCALAR_FACE_VALUE field kernels L0652 to L0657 (0.275 %).
- With the four finished mass-flux wall nests (`spec_wall_zz`, `spec_wall_rmw`, `mass_wall_zz`, `mass_wall_rmw`; finished, not yet merged to the generator branch): they complete L0880 (4.832 %) and L0882 (0.302 %), **+5.134 points, 9.158 % of total time (18.6 % of the eligible time), 54 loops**. Quote the committed figure until the merge is done.
- Shares are the survey model, not measurements. A loop's whole share is credited when a kernel covers its nests, so the figure is an upper bound for loops with several nests.

## 2. Columns and classes

| Column | Meaning |
|---|---|
| `#` | Rank among the not-yet-translated loops by modelled share (translated loops have no rank). |
| `Loop` | Survey loop id (also the id used in sidecar `[[kernel]]` entries and in `loop_claims.csv`). |
| `file:lines` | First and last line of the loop in the survey source. |
| `Share %` | Modelled share of total time from the survey. |
| `Domain` | pressure, velocity, species/combustion, radiation, solid phase, wall/BC, mass, other (the DIVERGENCE_PART_1 body and the like). Assigned by file and routine. |
| `Required generator features` | What the generator must handle that it does not handle for this loop today (or "front end accepts"). |
| `Blocker` | Blocker family or class: none, generator gap, table needed, callee flatten (T3), neighbour mesh (T4), race-uniqueness, pointer scratch (S2) and so on. |
| `Status` | translated (committed), in progress (work exists, not committed), claimed (an owner is named, no code yet), open, not planned (test-only or host). |
| `Class` | Work class, below. |
| `Owner` | Confirmed owner, or "open — available". |
| `Proposed owner` | Suggested lead for an unclaimed loop. A suggestion only: nobody owns the loop until a claim is recorded. |
| `Evid.` | `R`: the loop text was read for this list. `I`: the entry comes from the survey CSV and `blocked-loop-families.md` only; read the loop before starting. |

Work classes:
- **translated**: kernel committed.
- **claimed**: has an owner (in progress or claimed).
- **translatable-now**: the generator front end accepts the loop as it stands; what is missing is the kernel entry and the bitwise test.
- **needs-feature**: translatable once a named generator feature exists (a table kind, a layout contract, a loop shape).
- **blocked**: needs something outside the loop family itself (neighbour-mesh data, callee flattening, a uniqueness proof, a CSR table family).
- **host-side**: event-driven or file-output code that stays on the host (`REASSIGN_WALL_CELLS` L0823 and L0822, the RADF file write L1245). The modelled share is not credible for code that runs on an event only; the device deliverable, if any, is a trigger or a table refresh (`wall-bc-translation-plan.md`, section 11). Not counted as work to translate.
- **test-only**: guarded by a manufactured-solution or dead-code condition (`PERIODIC_TEST==7`, `.AND. .FALSE.`). Never executed in production runs. Their modelled shares (1.836 % in total) are overstated, so they are not planned and should not be ranked with real work.

## 3. Owners

| Owner | Holds |
|---|---|
| GPU Generator Engineer | Periodic-step cell loops, reductions and CYCLE handling, function callees, edge tables (L1376 to L1378 done; L1375 EVALUATE_RAMP open; L1392 and L1393 deferred), the DIVERGENCE_PART_1 nests (L0365, L0369, L0366, L0381 and the routine as a whole), the exact fixed-point sum switch (Decision B in `generator-decisions.md`). |
| Legacy Mapper | The partly translated nests (L0880, L0401), the S2 package (L0880, L0401, L0882, L0403, L0398), the pointer, array-constructor and alias machinery, golden signatures and the manifest, `port_merge_check`, and the cell-loop batch L1355 (Vreman), L1315 and L1316 (TEST_FILTER), L1321. (L1356 WALE released to GPU Wall Loops Engineer; L0876/L0877 released to GPU Mesh Data Loops Engineer under D-083.) |
| GPU Wall Loops Engineer | Wall hooks and tables, the CSR family, NEAR_SURFACE_GAS_VARIABLES per-surface tables, L1402, L1359, L0375, L1358; **L1356** WALE (`visc_wale`). |
| GPU Mesh Data Loops Engineer | Neighbour-mesh (T4): L1366 translated; L1363/L1364/L1367/L1399 not planned (host) under D-055 (same-mesh empty; `o2-edge-exchange-scope.md`); `s5_gsfv*` and `s5_pvf` builders (L1390 done); D-083 hand-written K2 `s5_clip_gather` for **L0876** and **L0877** (s5-gen `08f6016e28`; not markers-sidecar). |
| AMR Solid Phase Lead | SP4 scratch-sum reference and tests (see wall split), then L0394 and L0405 as CSR consumers. |
| AMR Pressure Solver Lead, AMR Species & Combustion Lead, AMR Radiation Lead | Proposed owners of the open loops of their domains. They own nothing until they claim. |

Interpretation used for the registers: where a routine is owned as a whole (DIVERGENCE_PART_1), its nests belong to that owner even if they are also partly translated. The Legacy Mapper holds the partly translated nests outside that routine.

### Wall split (agreed between the AMR Solid Phase Lead, the GPU Wall Loops Engineer and the Legacy Mapper)

| Who | Takes |
|---|---|
| AMR Solid Phase Lead | The SP4 scratch-sum reference and its tests (`wall_sums.py`, `test_wall_sums.py`): OBSTRUCTION%MASS and the gas-cell D_SOURCE and M_DOT_PPP sums, added in ascending wall index with no atomics. The WALL_BC loops **L1486**, **L1487**, **L1488**, **L1489** (see the plan in `wall-bc-translation-plan.md`). **L0394** and **L0405** were released: the GPU Wall Loops Engineer translated them (`wall_bc_dp`, `wall_spec_adv2`). |
| GPU Wall Loops Engineer | The CSR family itself, **L1402**, **L1359**, **L0375**, **L1358**, **L0394**, **L0405**. |
| Legacy Mapper | The S2 package (**L0880**, **L0401**, **L0882**, **L0403**, **L0398**). (**L0876** / **L0877** moved to GPU Mesh Data Loops Engineer under D-083.) |

L1488, L1489, L1486, L1487 (WALL_BC) are not translated; they depend on the SP4 reference and on callee flattening (T3). L1485 (`SURFACE_HEAT_TRANSFER` NIC>1) is **not planned (retired)** in the AMR route (host abort guard, patch 0010; D-065 Q1 / A-65). L1489 is the lateral thin-wall 1-D solve (wall.f90:192-194); `CALC_DEPOSITION` and `CALC_HVAC_BC` are called in L1488 (wall.f90:173, 178). L0823 and L0822 (`REASSIGN_WALL_CELLS`) are host-side.

### Work that is not a loop

These items gate several loops and have owners of their own; they are not rows of the table.

| Item | Owner | Gates |
|---|---|---|
| Per-surface NEAR_SURFACE_GAS_VARIABLES tables | GPU Wall Loops Engineer (in progress) | wall loops that read surface properties |
| CSR cell-to-wall table family (ascending wall index) | GPU Wall Loops Engineer | S3 (L0398, L0399, L0405, L0406), SP1 to SP4, P3 (L0394) |
| SP4 scratch-sum reference (`wall_sums.py`, `test_wall_sums.py`) | AMR Solid Phase Lead | L1488, L1489, L0394, L0405 and the WALL_BC family |
| Exact fixed-point sum switch (zone sums, Decision B) | GPU Generator Engineer | P1 zone sums, V1 reductions |
| Pointer, array-constructor and alias machinery | Legacy Mapper | S2 package |
| Wall index list per box and the `WALL_INDEX` cell table (on the s5-wall branch, being merged) | whoever completes that merge | L0399, L0406, L1366 |

## 4. Claim protocol

1. **Claim by message to the Legacy Mapper.** Name the loop ids (or one routine), the intended generator features and the tests you will write. Do not edit the Owner column yourself.
2. **One owner per loop.** A loop with an owner is not available to anyone else. If a claim overlaps an existing one, the Legacy Mapper answers with the current owner and the two owners settle it; the list does not change until they agree.
3. **The claim is recorded in the register.** The Legacy Mapper adds a row to `docs/amrex/loop_claims.csv` (key = loop id or `routine:NAME`, owner, status, note) and regenerates this file and `loop_work_list.csv`, so the Owner column always matches the register. A claim that shows no commit and no message for a long stretch is asked about, and released if the owner agrees.
4. **Status changes the same way.** Tell the Legacy Mapper when a loop moves from claimed to in progress, or when a commit lands. The status becomes `translated` only when the kernel entry is in the committed sidecar; the script reads that from the generator branch head and flips the row itself.
5. **Shared-tree rules** (the generator worktree is shared by several roles):
   - Commit only your own files, by path (`git add <path>` / `git commit <path>`), never `git add -A` or `git commit -a`.
   - A change to a function or a sidecar entry that another owner holds goes in a separate commit whose subject starts with `shared:`, and only after you have messaged that owner.
   - Put tests in your own test files, not in another owner's.
   - Do not touch another role's uncommitted files, and do not run the generator over the shared output directory while someone else holds the lock. Wrap any full run in the shared lock: `flock .s5gen.lock test/run_pvf.sh all` (lock file at the repository root of the generator worktree).
6. **Bitwise discipline.** A loop counts as translated only when its generated kernel is bitwise equal to the verbatim upstream loop text on the test data:
   - the reference is the upstream loop copied line for line from the survey source (use `test/make_r2_tests.py` for line-range references), not a rewrite;
   - all six flag sets (`O0 O2 O0omp O2omp O2omp_off O2omp_dpd`), both callee switches (`dpd` and `bind`), and 4 and 8 threads;
   - mutants that change one operand, one sign or one index must be caught (mutation check), and the negative checks of the generator must still fail on the cases they are meant to reject.
   Exact commands are in `generator-howto.md`.
7. **No unprovable assertions.** A wall-subscript store needs either a uniqueness proof (O3), a gather, or an idempotent marker justified by the source (the stored value depends on the target only). Do not assert uniqueness that the source does not guarantee (see S2 in `blocked-loop-families.md`).

### Recorded sign-offs (claim register notes)

The register (`loop_claims.csv`) carries the domain-lead verdicts below in the note column of the loop rows; each names its source document and section. Family-level sign-offs stay in `blocked-loop-families.md`.

| Loops | Lead | Verdict | Source |
|---|---|---|---|
| L0877, L0876 | Species & Combustion | clip gather prototype accept with three conditions; D-083 device ACCEPT-with-conditions then CLEARANCE YES (hand-written K2 on s5-gen `08f6016e28`) | `combustion/02` ("Clip gather prototype"); `combustion/05-signoff-check-mass-density-d083.md`; CLEARANCE docs `126662fc8f` |
| L1272, L0866, L0875, L0397 | Species & Combustion | host ACCEPT (generator path; gfortran bitwise) | `combustion/06-host-review-l1272-l0866-l0875-l0397.md` (docs `3e98931`) |
| L0401, L0398 | Species & Combustion | patches UP-0001 and UP-0002 (MP5 `Z_TEMP` pad) accept with a physics note | `combustion/02`, section "MP5 Z_TEMP pad" |
| L1239 | Radiation | accept | `radiation/00-r1-signoff.md`, section 2 (pasted on R1 in `blocked-loop-families.md`) |
| L1243, L1248 | Radiation | accept with change | same |
| L1245 | Radiation | accept (host code) | same |

Family-level Radiation R1 sign-off text from `radiation/00-r1-signoff.md` §2 is applied on the R1 **Sign-off** line in `blocked-loop-families.md`. The other radiation documents (`radiation/03`, `05-br-ilw-table-spec.md`, `06-fr062-sweep-kernel-design.md`) hold specifications and open review questions, not additional verdicts.

## 5. Ranking caveats

- Shares are the survey model. Nine test-only or dead-code loops (1.836 % in total, among them L1388, L1389, L1397, L1398, L1372, L1373, L0861, L0870, L0384) are in the table but will not run in production; they are listed so nobody picks them up as high-value work.
- Open does not mean easy. The largest open loops are blocked on shared infrastructure (callee flattening, neighbour-mesh data, the CSR family). The "Top 10 translatable now" table below lists the work that can start today.
- The shared generator worktree holds uncommitted work that is not in the register: a radiation sidecar and builder (`markers/rad_kernels.toml`, `s5_rad.py`) and more DIVERGENCE_PART_1 kernels. None of it is counted. Its owners should claim the loops by message so the Owner column and the status can follow.
- Rows with `Evid.` = `I` have not been read for this list. Read the loop before you claim it; if the text disagrees, tell the Legacy Mapper and the row is corrected.
- L1205 (COMPUTE_VELOCITY_ERROR) is classed geometry-deferred in the survey because its body is guarded by `IF (CC_IBM)` (see `pressure-velocity-h-loops.md`); it is outside the denominators.

<!-- GENERATED-BEGIN:worklist -->
### Counts by status (231 translation-eligible loops, 49.177 % modelled share)

| Status | Loops | Modelled share (% of total time) |
|---|---|---|
| claimed | 30 | 14.376 |
| in progress | 11 | 10.508 |
| not planned (host) | 11 | 2.957 |
| not planned (retired) | 9 | 1.838 |
| not planned (test-only) | 9 | 1.836 |
| open | 54 | 2.308 |
| translated | 104 | 12.787 |
| translated on a branch (merge pending) | 3 | 2.567 |

### Counts by work class

| Work class | Loops | Modelled share (% of total time) |
|---|---|---|
| blocked | 14 | 0.049 |
| claimed | 76 | 31.495 |
| host-side | 3 | 2.039 |
| needs-feature | 24 | 0.971 |
| test-only | 9 | 1.836 |
| translatable-now | 1 | 0.000 |
| translated | 104 | 12.787 |

### Counts by domain (not-yet-translated loops only)

| Domain | Open loops | Modelled share (%) | With a confirmed owner |
|---|---|---|---|
| other | 20 | 12.611 | 15 |
| wall/BC | 8 | 7.517 | 8 |
| solid phase | 14 | 4.851 | 10 |
| velocity | 37 | 3.885 | 13 |
| species/combustion | 4 | 3.524 | 4 |
| mass | 4 | 2.908 | 2 |
| pressure | 38 | 1.091 | 26 |
| radiation | 2 | 0.003 | 2 |

### Coverage figure
- Eligible denominator: 231 loops, 49.177 % of total modelled time (451 non-geometry loops, 69.444 %, minus 207 retired or replaced and 13 host-side). Geometry-deferred loops (378, 23.854 %) are outside every denominator.
- Inventory CSV status `translated+tested`: 46 loops, 2.553 %. Added by the sidecar at the generator ref (whole-loop kernels the CSV does not yet show: PATCH_VELOCITY_FLUX L1390 and the six GET_SCALAR_FACE_VALUE field kernels L0652-L0657): 58 loops, 10.234 %.
- **Translated at the generator ref: 104 loops, 12.787 % of total time (26.0 % of the eligible time).**
- The four finished mass-flux wall nests (spec_wall_zz, spec_wall_rmw, mass_wall_zz, mass_wall_rmw) complete L0880 (4.832 %, now partly translated through its cell nest) and L0882 (0.302 %, not translatable today). Counting them adds **5.134 points**, giving **17.921 % of total time (36.4 % of the eligible time)**, 106 loops. They are finished but not committed to s5-gen, so the committed figure above does not include them.
- The shares are inventory-model estimates; a loop's whole share is credited when a kernel covers its nests, so inner-nest time is an upper bound.

### Top 10 loops that are open and unclaimed (available)

| # | Loop | file:lines | Routine | Share % | Domain | Class | Blocker | Proposed owner |
|---|---|---|---|---|---|---|---|---|
| 16 | L1305 | turb.f90:1796-1861 | SYNTHETIC_TURBULENCE | 0.572 | velocity | needs-feature | table needed | - |
| 28 | L1317 | turb.f90:1027-1043 | TEST_FILTER | 0.284 | velocity | needs-feature | generator gap | - |
| 42 | L1275 | turb.f90:819-831 | CALC_VARDEN_LEONARD_TERM | 0.037 | velocity | needs-feature | generator gap | - |
| 43 | L1323 | turb.f90:633-645 | VARDEN_DYNSMAG | 0.037 | velocity | needs-feature | generator gap | - |
| 48 | L1325 | turb.f90:690-710 | VARDEN_DYNSMAG | 0.021 | velocity | needs-feature | generator gap | - |
| 51 | L0677 | func.f90:5183-5201 | PACK_CELL | 0.016 | other | blocked | neighbour mesh (T4) | - |
| 55 | L0911 | part.f90:4643-4763 | PARTICLE_MOMENTUM_TRANSFER | 0.015 | other | blocked | neighbour mesh (T4) | - |
| 62 | L1324 | turb.f90:666-678 | VARDEN_DYNSMAG | 0.011 | velocity | needs-feature | generator gap | - |
| 66 | L0890 | part.f90:605-607 | INSERT_VENT_PARTICLES | 0.008 | other | blocked | callee flatten (T3) | - |
| 69 | L1322 | turb.f90:570-581 | VARDEN_DYNSMAG | 0.005 | velocity | needs-feature | generator gap | - |

### Top 10 translatable now (front end accepts; the missing piece is the bitwise test)

| # | Loop | file:lines | Routine | Share % | Domain | What | Proposed owner |
|---|---|---|---|---|---|---|---|
| 92 | L0884 | part.f90:4606-4617 | CLIP_PARTICLE_DRAG | 0.000 | other | cell loop | - |

### Ranked list: 127 not-yet-translated loops

Columns: `Loop` is the survey id (sidecar and marker id) with `file:first-last` at FireX 36975d7; share is the inventory model, not a measurement; `Evid.` is `R` when the loop text was read for this list and `I` when the features come from the inventory CSV only.

| # | Loop | file:lines | Routine | Share % | Domain | Required generator features | Blocker | Status | Class | Owner | Proposed owner | Evid. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | L0365 | divg.f90:128-235 | DIVERGENCE_PART_1 | 7.679 | other | species loop N outermost; function callees (INTERPOLATE1D_UNIFORM, TENSOR_DIFFUSIVITY_MODEL); table column alias D_Z(:,N) | generator gap | in progress | claimed | GPU Generator Engineer | - | R |
| 2 | L0401 | divg.f90:995-1088 | SPECIES_ADVECTION_PART_1_NEW | 3.260 | species/combustion | wall gather; pointer scratch/array constructors; function callees (GET_SCALAR_FACE_VALUE_PT); Z_TEMP pad of the fourth element | MP5 wall element (upstream patch) | claimed | claimed | Legacy Mapper | - | I |
| 3 | L1488 | wall.f90:144-185 | WALL_BC | 2.688 | wall/BC | wall gather; function callees (SURFACE_HEAT_TRANSFER, CALCULATE_ZZ_F, ...); ragged per-wall table (BOUNDARY_ONE_D); scratch-sum reference (OBSTRUCTION%MASS, D_SOURCE) | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 4 | L0369 | divg.f90:287-421 | DIVERGENCE_PART_1 | 2.399 | other | species loop N outermost; cell loop nests | generator gap | in progress | claimed | GPU Generator Engineer | - | R |
| 5 | L0877 | mass.f90:868-939 | CHECK_MASS_DENSITY | 1.967 | mass | species loop N outermost; scatter into 7 cells (two-pass gather in source order); CYCLE; integer max-reduction CLIP_RHO_ZZ(N) | race-uniqueness | translated on a branch (merge pending) | claimed | GPU Mesh Data Loops Engineer | - | I |
| 6 | L1485 | wall.f90:888-947 | SURFACE_HEAT_TRANSFER | 1.832 | solid phase | wall gather; neighbour-mesh (OMESH copies through EWC%IIO_MIN..KKO_MAX); ragged per-wall table | neighbour mesh (T4) | not planned (retired) | claimed | AMR Solid Phase Lead | - | I |
| 7 | L0823 | init.f90:4940-5054 | REASSIGN_WALL_CELLS | 1.727 | solid phase | host-side event code (burn-away, control-driven creation or removal); allocation and derived-type copies | host-side event code | claimed | host-side | AMR Solid Phase Lead | - | I |
| 8 | L0381 | divg.f90:640-660 | DIVERGENCE_PART_1 | 1.491 | other | species loop N outermost; function callees (SPECIES_ADVECTION_PART_2, GET_SENSIBLE_ENTHALPY_Z); guarded call SET_EXIMRHOZZLIM_3D | generator gap | claimed | claimed | GPU Generator Engineer | - | R |
| 9 | L1355 | velo.f90:204-252 | COMPUTE_VISCOSITY | 1.467 | velocity | cell loop; private local arrays A_IJ(3,3), B_IJ(3,3) with inner 3x3 loops; live-out scalars (PRIVATE list); module scalar C_VREMAN as argument | generator gap | claimed | claimed | Legacy Mapper | - | R |
| 10 | L1489 | wall.f90:192-194 | WALL_BC | 1.370 | wall/BC | thin-wall loop (not a wall loop); function callees (SOLID_HEAT_TRANSFER); ragged per-wall table (BOUNDARY_ONE_D); one shared event (OBSTRUCTION%MASS=-1 for a burned-through thin wall) | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | R |
| 11 | L1363 | velo.f90:2891-3015 | MATCH_VELOCITY_FLUX | 0.994 | wall/BC | wall gather; neighbour-mesh (OMESH) | neighbour mesh (T4) | not planned (host) | claimed | GPU Mesh Data Loops Engineer | - | I |
| 12 | L1367 | velo.f90:1858-1897 | VELOCITY_BC | 0.795 | wall/BC | wall gather; neighbour-mesh (OMESH US/VS/WS) | neighbour mesh (T4) | not planned (host) | claimed | GPU Mesh Data Loops Engineer | - | I |
| 13 | L1399 | velo.f90:514-550 | VISCOSITY_BC | 0.710 | wall/BC | wall gather; neighbour-mesh (OMESH MU) | neighbour mesh (T4) | not planned (host) | claimed | GPU Mesh Data Loops Engineer | - | I |
| 14 | L0384 | divg.f90:696-716 | DIVERGENCE_PART_1 | 0.658 | other | species inner loop; source function (VD2D_MMS_Z_SRC); pointer SM | test-only | not planned (test-only) | test-only | GPU Generator Engineer | - | R |
| 15 | L0876 | mass.f90:799-849 | CHECK_MASS_DENSITY | 0.599 | mass | K,J nest; scatter into 7 cells (two-pass gather) | race-uniqueness | translated on a branch (merge pending) | claimed | GPU Mesh Data Loops Engineer | - | I |
| 16 | L1305 | turb.f90:1796-1861 | SYNTHETIC_TURBULENCE | 0.572 | velocity | cell/face loop; non-perfect nest (outer loop) | table needed | open | needs-feature | open — available | - | I |
| 17 | L1470 | wall.f90:3635-3677 | HT3D_TEMPERATURE_EXCHANGE | 0.500 | solid phase | wall gather; ragged per-wall table (BOUNDARY_THR_D); store into TMP | neighbour mesh (T4) | open | claimed | AMR Solid Phase Lead | - | I |
| 18 | L1471 | wall.f90:3681-3717 | HT3D_TEMPERATURE_EXCHANGE | 0.465 | solid phase | wall gather; ragged per-wall table (BOUNDARY_THR_D) | neighbour mesh (T4) | open | claimed | AMR Solid Phase Lead | - | I |
| 19 | L1364 | velo.f90:1376-1400 | NO_FLUX | 0.426 | pressure | wall gather; neighbour-mesh (OMESH(NOM)%H/HS); overlap sum in K,J,I order then divide | neighbour mesh (T4) | not planned (host) | claimed | GPU Mesh Data Loops Engineer | - | R |
| 20 | L1486 | wall.f90:109-119 | WALL_BC | 0.407 | wall/BC | wall gather; function callees (ASSIGN_GHOST_VALUE, NEAR_SURFACE_GAS_VARIABLES) | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 21 | L1487 | wall.f90:124-137 | WALL_BC | 0.354 | wall/BC | wall gather; function callees | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 22 | L1365 | velo.f90:1404-1459 | NO_FLUX | 0.347 | pressure | obstruction box loop (outer N over OBSTRUCTION); CELL%SOLID mask on both cells; pointer re-aim HP=>H/HS (hoist); PREDICTOR switch as argument | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 23 | L1389 | velo.f90:1048-1054 | MMS_VELOCITY_FLUX | 0.325 | velocity | face loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 24 | L1388 | velo.f90:1040-1046 | MMS_VELOCITY_FLUX | 0.315 | velocity | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 25 | L0822 | init.f90:4917-4936 | REASSIGN_WALL_CELLS | 0.310 | solid phase | host-side event code; neighbour-mesh reads | host-side event code | claimed | host-side | AMR Solid Phase Lead | - | I |
| 26 | L1315 | turb.f90:979-995 | TEST_FILTER | 0.284 | velocity | cell/face loop; live-out scalars (PRIVATE list) | generator gap | claimed | claimed | Legacy Mapper | - | I |
| 27 | L1316 | turb.f90:1003-1019 | TEST_FILTER | 0.284 | velocity | cell/face loop; live-out scalars (PRIVATE list) | generator gap | claimed | claimed | Legacy Mapper | - | I |
| 28 | L1317 | turb.f90:1027-1043 | TEST_FILTER | 0.284 | velocity | cell/face loop; rank/subscript contract | generator gap | open | needs-feature | open — available | - | I |
| 29 | L0398 | divg.f90:835-939 | ENTHALPY_ADVECTION_NEW | 0.260 | species/combustion | wall gather; pointer scratch/array constructors; CSR cell-to-wall gather (ascending IW); pointer remap U_TEMP=>U_WORK | pointer scratch (S2) | claimed | claimed | Legacy Mapper | - | I |
| 30 | L1452 | wall.f90:319-339 | ASSIGN_GHOST_VALUE | 0.199 | wall/BC | wall gather; neighbour-mesh | neighbour mesh (T4) | open | claimed | AMR Solid Phase Lead | - | I |
| 31 | L0861 | mass.f90:484-495 | DENSITY | 0.171 | mass | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 32 | L0870 | mass.f90:666-677 | DENSITY | 0.171 | mass | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 33 | L1321 | turb.f90:501-534 | VARDEN_DYNSMAG | 0.105 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | claimed | claimed | Legacy Mapper | - | I |
| 34 | L0364 | divg.f90:93-95 | DIVERGENCE_PART_1 | 0.099 | other | wall loop; wall gather; function callees; module scalar as argument | callee flatten (T3) | claimed | claimed | GPU Generator Engineer | - | I |
| 35 | L1121 | pres.f90:4012-4026 | CHECK_UNSUPPORTED_MESH | 0.086 | pressure | mesh/rank loop; MY_RANK test | blocked | open | claimed | AMR Pressure Solver Lead | - | R |
| 36 | L0371 | divg.f90:463-471 | DIVERGENCE_PART_1 | 0.068 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 37 | L1372 | velo.f90:1770-1778 | VELOCITY_CORRECTOR | 0.049 | pressure | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 38 | L1373 | velo.f90:1779-1787 | VELOCITY_CORRECTOR | 0.049 | pressure | face loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 39 | L1397 | velo.f90:1648-1656 | VELOCITY_PREDICTOR | 0.049 | pressure | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 40 | L1398 | velo.f90:1657-1665 | VELOCITY_PREDICTOR | 0.049 | pressure | face loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 41 | L0379 | divg.f90:608-616 | DIVERGENCE_PART_1 | 0.037 | other | cell/face loop; rank/subscript contract | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 42 | L1275 | turb.f90:819-831 | CALC_VARDEN_LEONARD_TERM | 0.037 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 43 | L1323 | turb.f90:633-645 | VARDEN_DYNSMAG | 0.037 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 44 | L1391 | velo.f90:1243-1253 | VELOCITY_FLUX_CYLINDRICAL | 0.026 | velocity | cell loop; pointer aliases to WORK2, WORK5 | none | open | claimed | AMR Pressure Solver Lead | - | R |
| 45 | L0377 | divg.f90:573-582 | DIVERGENCE_PART_1 | 0.025 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 46 | L0380 | divg.f90:622-629 | DIVERGENCE_PART_1 | 0.025 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 47 | L0386 | divg.f90:757-767 | DIVERGENCE_PART_1 | 0.022 | other | wall loop; wall gather; live-out scalars (PRIVATE list); zone table | table needed | claimed | claimed | GPU Generator Engineer | - | I |
| 48 | L1325 | turb.f90:690-710 | VARDEN_DYNSMAG | 0.021 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 49 | L0630 | func.f90:5509-5524 | BLOCK_CELL | 0.020 | other | cell/face loop; wall gather; derived-type designator table | neighbour mesh (T4) | not planned (host) | claimed | GPU Mesh Data Loops Engineer | - | I |
| 50 | L0817 | init.f90:4878-4889 | CREATE_OR_REMOVE_OBST | 0.017 | solid phase | cell/face loop; derived-type designator table | generator gap | claimed | claimed | AMR Solid Phase Lead | - | I |
| 51 | L0677 | func.f90:5183-5201 | PACK_CELL | 0.016 | other | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | - | I |
| 52 | L1381 | velo.f90:992-1000 | CORIOLIS_FORCE | 0.016 | velocity | cell loop; constant-subscript rank-1 array OVEC(n) as by-value scalars | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 53 | L1382 | velo.f90:1006-1014 | CORIOLIS_FORCE | 0.016 | velocity | cell loop; constant-subscript rank-1 array OVEC(n) | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 54 | L1383 | velo.f90:1020-1028 | CORIOLIS_FORCE | 0.016 | velocity | face loop; constant-subscript rank-1 array OVEC(n) | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 55 | L0911 | part.f90:4643-4763 | PARTICLE_MOMENTUM_TRANSFER | 0.015 | other | wall loop; wall gather; neighbour-mesh; pointer alias | neighbour mesh (T4) | open | blocked | open — available | - | I |
| 56 | L1392 | velo.f90:1270-1302 | VELOCITY_FLUX_CYLINDRICAL | 0.014 | velocity | edge tables; K,I nest with J fixed | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 57 | L1393 | velo.f90:1306-1337 | VELOCITY_FLUX_CYLINDRICAL | 0.014 | velocity | edge tables; K,I nest with J fixed | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 58 | L0373 | divg.f90:499-505 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 59 | L0378 | divg.f90:591-597 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 60 | L0382 | divg.f90:668-674 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 61 | L0383 | divg.f90:681-687 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop; rank/subscript contract | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 62 | L1324 | turb.f90:666-678 | VARDEN_DYNSMAG | 0.011 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 63 | L1385 | velo.f90:896-903 | DIRECT_FORCE | 0.011 | velocity | face loop; constant-subscript rank-1 array FVEC(n) as by-value scalar; ramp factor computed on the host | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 64 | L1386 | velo.f90:917-924 | DIRECT_FORCE | 0.011 | velocity | face loop; FVEC(2) as scalar | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 65 | L1387 | velo.f90:938-945 | DIRECT_FORCE | 0.011 | velocity | face loop; FVEC(3) as scalar | generator gap | in progress | claimed | AMR Pressure Solver Lead | - | R |
| 66 | L0890 | part.f90:605-607 | INSERT_VENT_PARTICLES | 0.008 | other | wall loop; wall gather; function callees; module scalar as argument | callee flatten (T3) | open | blocked | open — available | - | I |
| 67 | L1209 | pres.f90:65-228 | PRESSURE_SOLVER_COMPUTE_RHS | 0.008 | pressure | wall gather; neighbour-mesh (MESHES(NOM)%DX); function callees (EVALUATE_RAMP); vent table (VENTS(WC%VENT_INDEX)); wall-keyed 2-D outputs | neighbour mesh (T4) | not planned (host) | claimed | AMR Pressure Backend Implementer | - | R |
| 68 | L1224 | pres.f90:544-567 | TUNNEL_POISSON_SOLVER | 0.007 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | claimed | AMR Pressure Solver Lead | - | I |
| 69 | L1322 | turb.f90:570-581 | VARDEN_DYNSMAG | 0.005 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 70 | L1204 | pres.f90:1931-1994 | ULMAT_SOLVE_ZONE | 0.003 | pressure | wall loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 71 | L0596 | fire.f90:1902-1918 | COMBUSTION_BC | 0.002 | species/combustion | wall loop; wall gather; neighbour-mesh; derived-type designator table | neighbour mesh (T4) | not planned (host) | claimed | AMR Species & Combustion Lead | - | I |
| 72 | L1198 | pres.f90:1837-1846 | ULMAT_SOLVE_ZONE | 0.002 | pressure | cell/face loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 73 | L1211 | pres.f90:250-260 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop | none | in progress | claimed | AMR Pressure Backend Implementer | - | R |
| 74 | L1212 | pres.f90:267-277 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop; transposed output subscripts | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 75 | L1213 | pres.f90:282-292 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop; transposed output subscripts | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 76 | L1214 | pres.f90:297-307 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop; transposed output subscripts | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 77 | L1245 | radi.f90:5044-5059 | RADIATION_FVM | 0.002 | radiation | file output | host (I/O) | not planned (host) | host-side | AMR Radiation Lead | - | I |
| 78 | L1273 | soot.f90:412-472 | SOOT_SURFACE_OXIDATION | 0.002 | species/combustion | wall loop; wall gather; pointer alias | neighbour mesh (T4) | open | claimed | AMR Species & Combustion Lead | - | I |
| 79 | L0680 | func.f90:2845-2847 | COMPUTE_WIND_COMPONENTS | 0.001 | other | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | - | I |
| 80 | L1130 | pres.f90:5580-5601 | GET_H_REGFACES | 0.001 | pressure | wall gather; rank-4 LOGICAL table | table needed | open | claimed | AMR Pressure Solver Lead | - | R |
| 81 | L1154 | pres.f90:2031-2047 | ULMAT_GET_H_REGFACES | 0.001 | pressure | wall gather; rank-4 LOGICAL table | table needed | open | needs-feature | open — available | AMR Pressure Solver Lead | I |
| 82 | L1190 | pres.f90:1686-1694 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 83 | L1192 | pres.f90:1713-1720 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 84 | L1194 | pres.f90:1782-1790 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 85 | L1196 | pres.f90:1808-1815 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 86 | L1206 | pres.f90:729-742 | PRESSURE_SOLVER_CHECK_RESIDUALS | 0.001 | pressure | cell loop; bounds (WORK8 view, policy.arrays) | generator gap | open | claimed | AMR Pressure Solver Lead | - | R |
| 87 | L1208 | pres.f90:768-787 | PRESSURE_SOLVER_CHECK_RESIDUALS | 0.001 | pressure | cell loop; bounds (WORK8 view); reduction by the caller (MAXVAL, MAXLOC) | generator gap | open | claimed | AMR Pressure Solver Lead | - | R |
| 88 | L1248 | radi.f90:3651-3662 | INTERPOLATE_IL | 0.001 | radiation | wall loop; ragged per-wall table BR_ILW(NRA,bands,wall); ILW_OLD copy kept (ping-pong table or per-thread scratch); I_INTP accumulation ascending and serial per angle | table needed | translated on a branch (merge pending) | claimed | AMR Radiation Lead | - | I |
| 89 | L1318 | turb.f90:1054-1059 | TEST_FILTER | 0.001 | velocity | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | - | I |
| 90 | L1319 | turb.f90:1063-1068 | TEST_FILTER | 0.001 | velocity | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | - | I |
| 91 | L1320 | turb.f90:1072-1077 | TEST_FILTER | 0.001 | velocity | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | - | I |
| 92 | L0884 | part.f90:4606-4617 | CLIP_PARTICLE_DRAG | 0.000 | other | cell/face loop | none | open | translatable-now | open — available | - | I |
| 93 | L1203 | pres.f90:1923-1925 | ULMAT_SOLVE_ZONE | 0.000 | pressure | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | I |
| 94 | L1207 | pres.f90:758-764 | PRESSURE_SOLVER_CHECK_RESIDUALS | 0.000 | pressure | cell loop | none | in progress | claimed | AMR Pressure Backend Implementer | - | R |
| 95 | L1210 | pres.f90:238-245 | PRESSURE_SOLVER_COMPUTE_RHS | 0.000 | pressure | cell loop (K,I with J=1); non-perfect nest | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 96 | L1215 | pres.f90:394-400 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop | none | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 97 | L1216 | pres.f90:404-410 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop; transposed input subscripts | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 98 | L1217 | pres.f90:414-420 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop; transposed input subscripts | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 99 | L1218 | pres.f90:424-430 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop; transposed input subscripts | generator gap | not planned (retired) | claimed | AMR Pressure Backend Implementer | - | R |
| 100 | L1219 | pres.f90:438-440 | PRESSURE_SOLVER_FFT | 0.000 | pressure | non-perfect nest (section assignment); 1-D table H_BAR(I_OFFSET+I) | generator gap | open | claimed | AMR Pressure Solver Lead | - | R |
| 101 | L1220 | pres.f90:450-462 | PRESSURE_SOLVER_FFT | 0.000 | pressure | K,J-only nest; integer BC-code scalars as arguments; reads BXS/BXF | generator gap | not planned (host) | claimed | AMR Pressure Backend Implementer | - | R |
| 102 | L1221 | pres.f90:466-477 | PRESSURE_SOLVER_FFT | 0.000 | pressure | K,I-only nest; integer BC-code scalars as arguments | generator gap | not planned (host) | claimed | AMR Pressure Backend Implementer | - | R |
| 103 | L1222 | pres.f90:481-492 | PRESSURE_SOLVER_FFT | 0.000 | pressure | J,I-only nest; integer BC-code scalars as arguments | generator gap | not planned (host) | claimed | AMR Pressure Backend Implementer | - | R |
| 104 | L1223 | pres.f90:531-537 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather | neighbour mesh (T4) | open | claimed | AMR Pressure Solver Lead | - | I |
| 105 | L1225 | pres.f90:576-582 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | claimed | AMR Pressure Solver Lead | - | I |
| 106 | L1226 | pres.f90:583-589 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | claimed | AMR Pressure Solver Lead | - | I |
| 107 | L1227 | pres.f90:591-596 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | claimed | AMR Pressure Solver Lead | - | I |
| 108 | L1289 | turb.f90:882-884 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 109 | L1290 | turb.f90:887-889 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 110 | L1291 | turb.f90:892-894 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 111 | L1292 | turb.f90:897-899 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 112 | L1293 | turb.f90:904-906 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 113 | L1294 | turb.f90:909-911 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 114 | L1295 | turb.f90:914-916 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 115 | L1296 | turb.f90:919-921 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 116 | L1297 | turb.f90:926-928 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 117 | L1298 | turb.f90:931-933 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 118 | L1299 | turb.f90:936-938 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 119 | L1300 | turb.f90:941-943 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 120 | L1330 | vege.f90:618-623 | GET_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 121 | L1331 | vege.f90:624-629 | GET_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 122 | L1332 | vege.f90:630-639 | GET_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 123 | L1333 | vege.f90:678-699 | FILL_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; neighbour-mesh; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 124 | L1334 | vege.f90:885-894 | LEVEL_SET_ADVECT_FLUX | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | claimed | claimed | AMR Solid Phase Lead | - | I |
| 125 | L1335 | vege.f90:896-905 | LEVEL_SET_ADVECT_FLUX | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | claimed | claimed | AMR Solid Phase Lead | - | I |
| 126 | L1336 | vege.f90:907-914 | LEVEL_SET_ADVECT_FLUX | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | claimed | claimed | AMR Solid Phase Lead | - | I |
| 127 | L1340 | vege.f90:763-857 | LEVEL_SET_SPREAD_RATE | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | claimed | claimed | AMR Solid Phase Lead | - | I |
<!-- GENERATED-END:worklist -->

## 6. Species & Combustion

Scope: `species.f90` and the combustion modules (the survey records no loop there), `fire.f90` (two eligible loops, L0596 and L0608; five more are geometry-deferred), `soot.f90` (L1272 and L1273 eligible; one geometry-deferred), `mass.f90`, and the species parts of `divg.f90` (the S2, S3, S4 loops and the CHECK_MASS_DENSITY gather pieces). Proposed owner of every open loop in this section: **AMR Species & Combustion Lead**, working in own files and own tests only; they message the owner of a loop or a function before any shared commit (`shared:` rule above). The loops already claimed stay with the owners shown.

How the families map to this section (details and sign-offs in `blocked-loop-families.md`):
- **S1** (CHECK_MASS_DENSITY scatter, L0877 and L0876): D-083 CLOSED — hand-written K2 two-pass gather on s5-gen `08f6016e28` (`s5_clip_gather`; not markers-sidecar). Owner: **GPU Mesh Data Loops Engineer**. Species CLEARANCE docs `126662fc8f`; Integration `[[ported]]` restored docs `b46f6d0`; `ci_checks --strict` PASS. The gather adds in the serial source order `(I,J,K-1)`, `(I,J-1,K)`, `(I-1,J,K)`, self, `(I+1,J,K)`, `(I,J+1,K)`, `(I,J,K+1)`, interior targets only, with the same `SUM(MASS_N)` order; clip-gather conditions from `combustion/02` still apply. Earlier translated pieces (`delta_rho_zz_zero`, `rho_zz_clip_assign`) remain. This is **not** host-blocked S1 any more.
- **S2** (wall nests with pointer scratch and array constructors: L0880, L0401, L0882, L0403, L0398): the Legacy Mapper's package. The four finished mass-flux wall nests cover the L0880 and L0882 parts; the divergence-side nests L0401 and L0403 reuse the same nests and need the fourth `Z_TEMP` element padded with 0 (MP5 bit tests wait for upstream patches UP-0001 and UP-0002).
- **S3** (wall scatter into `U_DOT_DEL_RHO_*`, cell loops reading `CELL%WALL_INDEX`: L0398, L0399, L0405, L0406): L0398 stays with the Legacy Mapper; L0405 is translated by the GPU Wall Loops Engineer (`wall_spec_adv2`); L0399 and L0406 are open and need the six-flag mask table (`WALL_INDEX`), which exists as a `cell_int` table on the s5-wall branch and is not merged.
- **S4** (SETTLING_VELOCITY, L1272): translated on s5-gen (`fe8017c109`, settle_* kernels); host ACCEPT in `combustion/06` (docs `3e98931`). Companion host translations L0866 / L0875 / L0397 also ACCEPTed there (`tmp_pbar_*`, `rho_h_s_p`).
- Plain cell loops of `mass.f90` and `divg.f90` (L0865, L0874, L0872, L0879, L0878, L0881, L0402): the front end accepts them; what is missing is the sidecar entry and the bitwise test.
- Test-only: L0861 and L0870 (`IF (PERIODIC_TEST==7)`, guards at `mass.f90:483` and `mass.f90:665`) and L0384 (`divg.f90:693`).
- Not in this section: the DIVERGENCE_PART_1 nests (L0365, L0369, L0366, L0381) belong to the GPU Generator Engineer even though they handle species.

Classification of each loop (translatable-now / needs-feature / blocked / claimed) follows. Classes are the work classes of section 2.

<!-- GENERATED-BEGIN:species -->
8 not-yet-translated loops in the species/combustion and mass domains, 6.432 % modelled share.


#### claimed (6 loops, 6.090 %)

| Loop | file:lines | Routine | Share % | Features / blocker | Closest kernel | Status | Owner / proposed |
|---|---|---|---|---|---|---|---|
| L0401 | divg.f90:995-1088 | SPECIES_ADVECTION_PART_1_NEW | 3.260 | wall gather; pointer scratch/array constructors; function callees (GET_SCALAR_FACE_VALUE_PT); Z_TEMP pad of the fourth element / kernels rho_z_p_divg and spec_wall_zz are committed; the MP5 limiter reads a fourth wall element that the upstream wall nest leaves unset, so the bitwise claim for MP5 waits for upstream patc | gsfv_*, rho_z_p_divg | claimed | Legacy Mapper |
| L0877 | mass.f90:868-939 | CHECK_MASS_DENSITY | 1.967 | species loop N outermost; scatter into 7 cells (two-pass gather in source order); CYCLE; integer max-reduction CLIP_RHO_ZZ(N) / S1: D-083 hand-written K2 two-pass landed (s5_clip_gather on s5-gen 08f6016e28; not markers-sidecar) | rho_zz_clip_assign, delta_rho_zz_zero | translated on a branch (merge pending) | GPU Mesh Data Loops Engineer |
| L0876 | mass.f90:799-849 | CHECK_MASS_DENSITY | 0.599 | K,J nest; scatter into 7 cells (two-pass gather) / S1: D-083 hand-written K2 two-pass landed (s5_clip_gather; not markers-sidecar) | rho_zz_clip_assign | translated on a branch (merge pending) | GPU Mesh Data Loops Engineer |
| L0398 | divg.f90:835-939 | ENTHALPY_ADVECTION_NEW | 0.260 | wall gather; pointer scratch/array constructors; CSR cell-to-wall gather (ascending IW); pointer remap U_TEMP=>U_WORK / S2 and S3: accumulate through the wall gas-cell subscript | gsfv_*, CSR family (not built) | claimed | Legacy Mapper |
| L0596 | fire.f90:1902-1918 | COMBUSTION_BC | 0.002 | wall loop; wall gather; neighbour-mesh; derived-type designator table / needs neighbour-mesh or ragged per-wall data (inventory) | none (new); pvf_* shows the host launch rule | not planned (host) | AMR Species & Combustion Lead |
| L1273 | soot.f90:412-472 | SOOT_SURFACE_OXIDATION | 0.002 | wall loop; wall gather; pointer alias / needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open | AMR Species & Combustion Lead |

#### test-only (2 loops, 0.342 %)

| Loop | file:lines | Routine | Share % | Features / blocker | Closest kernel | Status | Owner / proposed |
|---|---|---|---|---|---|---|---|
| L0861 | mass.f90:484-495 | DENSITY | 0.171 | cell loop; source function / IF (PERIODIC_TEST==7) at mass.f90:483 | none | not planned (test-only) | open — available |
| L0870 | mass.f90:666-677 | DENSITY | 0.171 | cell loop; source function / IF (PERIODIC_TEST==7) at mass.f90:665 | none | not planned (test-only) | open — available |

Already translated in these domains: 35 loops, 8.805 %.
<!-- GENERATED-END:species -->

## 7. Regenerating

```
python3 tools/inventory/loop_work_list.py            # rewrite the generated blocks and loop_work_list.csv
python3 tools/inventory/loop_work_list.py --check    # fail if the committed files differ from a fresh run
python3 tools/inventory/loop_work_list.py --help
```

The script reads the survey CSV, `docs/amrex/loop_claims.csv` and the generator sidecar at the generator branch head (`--gen-repo`, `--gen-ref`). It never writes into the generator worktree. Options `--no-verify-lines` skips the check against the survey source.
