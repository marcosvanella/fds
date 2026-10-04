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
| Legacy Mapper | The partly translated nests (L0880, L0401, L0877), the S2 package (L0880, L0401, L0882, L0403, L0398), L0876 and L0877, the pointer, array-constructor and alias machinery, golden signatures and the manifest, `port_merge_check`, and the cell-loop batch L1355 (Vreman), L1315 and L1316 (TEST_FILTER), L1356 (WALE), L1321. |
| GPU Wall Loops Engineer | Wall hooks and tables, the CSR family, NEAR_SURFACE_GAS_VARIABLES per-surface tables, L1402, L1359, L0375, L1358. |
| GPU Mesh Data Loops Engineer | Neighbour-mesh (T4) loops: L1366, L1363, L1364, L1367, L1399; `s5_gsfv*` and `s5_pvf` builders (L1390 done). |
| AMR Solid Phase Lead | SP4 scratch-sum reference and tests (see wall split), then L0394 and L0405 as CSR consumers. |
| AMR Pressure Solver Lead, AMR Species & Combustion Lead, AMR Radiation Lead | Proposed owners of the open loops of their domains. They own nothing until they claim. |

Interpretation used for the registers: where a routine is owned as a whole (DIVERGENCE_PART_1), its nests belong to that owner even if they are also partly translated. The Legacy Mapper holds the partly translated nests outside that routine.

### Wall split (agreed between the AMR Solid Phase Lead, the GPU Wall Loops Engineer and the Legacy Mapper)

| Who | Takes |
|---|---|
| AMR Solid Phase Lead | The SP4 scratch-sum reference and its tests (`wall_sums.py`, `test_wall_sums.py`): OBSTRUCTION%MASS and the gas-cell D_SOURCE and M_DOT_PPP sums, added in ascending wall index with no atomics. The WALL_BC loops **L1486**, **L1487**, **L1488**, **L1489** (see the plan in `wall-bc-translation-plan.md`). **L0394** and **L0405** were released: the GPU Wall Loops Engineer translated them (`wall_bc_dp`, `wall_spec_adv2`). |
| GPU Wall Loops Engineer | The CSR family itself, **L1402**, **L1359**, **L0375**, **L1358**, **L0394**, **L0405**. |
| Legacy Mapper | The S2 package (**L0880**, **L0401**, **L0882**, **L0403**, **L0398**) and **L0876** and **L0877**. |

L1488, L1489, L1486, L1487 (WALL_BC) and L1485 (SURFACE_HEAT_TRANSFER) are not translated; they depend on the SP4 reference and on callee flattening (T3) or the neighbour-mesh data (T4). L1489 is the lateral thin-wall 1-D solve (wall.f90:192-194); `CALC_DEPOSITION` and `CALC_HVAC_BC` are called in L1488 (wall.f90:173, 178). L0823 and L0822 (`REASSIGN_WALL_CELLS`) are host-side.

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
| claimed | 45 | 22.243 |
| in progress | 5 | 15.268 |
| not planned (host) | 1 | 0.002 |
| not planned (test-only) | 9 | 1.836 |
| open | 109 | 5.169 |
| translated | 62 | 4.659 |

### Counts by work class

| Work class | Loops | Modelled share (% of total time) |
|---|---|---|
| blocked | 27 | 1.338 |
| claimed | 49 | 35.784 |
| host-side | 3 | 2.039 |
| needs-feature | 60 | 3.304 |
| test-only | 9 | 1.836 |
| translatable-now | 21 | 0.217 |
| translated | 62 | 4.659 |

### Counts by domain (not-yet-translated loops only)

| Domain | Open loops | Modelled share (%) | With a confirmed owner |
|---|---|---|---|
| other | 20 | 12.611 | 14 |
| mass | 14 | 8.471 | 10 |
| wall/BC | 8 | 7.517 | 7 |
| species/combustion | 9 | 5.181 | 5 |
| solid phase | 14 | 4.851 | 2 |
| velocity | 62 | 4.737 | 9 |
| pressure | 39 | 1.147 | 2 |
| radiation | 3 | 0.003 | 3 |

### Coverage figure
- Eligible denominator: 231 loops, 49.177 % of total modelled time (451 non-geometry loops, 69.444 %, minus 207 retired or replaced and 13 host-side). Geometry-deferred loops (378, 23.854 %) are outside every denominator.
- Inventory CSV status `translated+tested`: 46 loops, 2.553 %. Added by the sidecar at the generator ref (whole-loop kernels the CSV does not yet show: PATCH_VELOCITY_FLUX L1390 and the six GET_SCALAR_FACE_VALUE field kernels L0652-L0657): 16 loops, 2.106 %.
- **Translated at the generator ref: 62 loops, 4.659 % of total time (9.5 % of the eligible time).**
- The four finished mass-flux wall nests (spec_wall_zz, spec_wall_rmw, mass_wall_zz, mass_wall_rmw) complete L0880 (4.832 %, now partly translated through its cell nest) and L0882 (0.302 %, not translatable today). Counting them adds **5.134 points**, giving **9.793 % of total time (19.9 % of the eligible time)**, 64 loops. They are finished but not committed to s5-gen, so the committed figure above does not include them.
- The shares are inventory-model estimates; a loop's whole share is credited when a kernel covers its nests, so inner-nest time is an upper bound.

### Top 10 loops that are open and unclaimed (available)

| # | Loop | file:lines | Routine | Share % | Domain | Class | Blocker | Proposed owner |
|---|---|---|---|---|---|---|---|---|
| 12 | L1272 | soot.f90:54-175 | SETTLING_VELOCITY | 1.325 | species/combustion | needs-feature | table needed | AMR Species & Combustion Lead |
| 18 | L1305 | turb.f90:1796-1861 | SYNTHETIC_TURBULENCE | 0.572 | velocity | needs-feature | table needed | - |
| 19 | L1470 | wall.f90:3635-3677 | HT3D_TEMPERATURE_EXCHANGE | 0.500 | solid phase | blocked | neighbour mesh (T4) | AMR Solid Phase Lead |
| 20 | L1471 | wall.f90:3681-3717 | HT3D_TEMPERATURE_EXCHANGE | 0.465 | solid phase | blocked | neighbour mesh (T4) | AMR Solid Phase Lead |
| 24 | L1365 | velo.f90:1404-1459 | NO_FLUX | 0.347 | pressure | needs-feature | generator gap | AMR Pressure Solver Lead |
| 31 | L1317 | turb.f90:1027-1043 | TEST_FILTER | 0.284 | velocity | needs-feature | generator gap | - |
| 33 | L1452 | wall.f90:319-339 | ASSIGN_GHOST_VALUE | 0.199 | wall/BC | blocked | neighbour mesh (T4) | AMR Solid Phase Lead |
| 41 | L1360 | velo.f90:413-435 | COMPUTE_STRAIN_RATE | 0.089 | velocity | needs-feature | generator gap | - |
| 42 | L1121 | pres.f90:4012-4026 | CHECK_UNSUPPORTED_MESH | 0.086 | pressure | blocked | blocked | AMR Pressure Solver Lead |
| 43 | L0397 | divg.f90:814-822 | ENTHALPY_ADVECTION_NEW | 0.077 | species/combustion | needs-feature | generator gap | AMR Species & Combustion Lead |

### Top 10 translatable now (front end accepts; the missing piece is the bitwise test)

| # | Loop | file:lines | Routine | Share % | Domain | What | Proposed owner |
|---|---|---|---|---|---|---|---|
| 53 | L1314 | turb.f90:1458-1479 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | velocity | cell loop | - |
| 68 | L1391 | velo.f90:1243-1253 | VELOCITY_FLUX_CYLINDRICAL | 0.026 | velocity | cylindrical vorticity and stress (OMY, TXZ) | - |
| 73 | L1354 | velo.f90:187-196 | COMPUTE_VISCOSITY | 0.021 | velocity | cell loop | - |
| 77 | L1379 | velo.f90:970-978 | CORIOLIS_FORCE | 0.016 | velocity | CORIOLIS_FORCE: cell-centred velocities UP,VP,WP | - |
| 88 | L1276 | turb.f90:187-194 | COMPRESSION_WAVE | 0.011 | velocity | cell loop | - |
| 89 | L1277 | turb.f90:195-202 | COMPRESSION_WAVE | 0.011 | velocity | cell loop | - |
| 90 | L1278 | turb.f90:203-210 | COMPRESSION_WAVE | 0.011 | velocity | cell loop | - |
| 91 | L1279 | turb.f90:211-218 | COMPRESSION_WAVE | 0.011 | velocity | cell loop | - |
| 99 | L1280 | turb.f90:221-227 | COMPRESSION_WAVE | 0.005 | velocity | cell loop | - |
| 100 | L1281 | turb.f90:228-234 | COMPRESSION_WAVE | 0.005 | velocity | cell loop | - |

### Ranked list: 169 not-yet-translated loops

Columns: `Loop` is the survey id (sidecar and marker id) with `file:first-last` at FireX 36975d7; share is the inventory model, not a measurement; `Evid.` is `R` when the loop text was read for this list and `I` when the features come from the inventory CSV only.

| # | Loop | file:lines | Routine | Share % | Domain | Required generator features | Blocker | Status | Class | Owner | Proposed owner | Evid. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | L0365 | divg.f90:128-235 | DIVERGENCE_PART_1 | 7.679 | other | species loop N outermost; function callees (INTERPOLATE1D_UNIFORM, TENSOR_DIFFUSIVITY_MODEL); table column alias D_Z(:,N) | generator gap | in progress | claimed | GPU Generator Engineer | - | R |
| 2 | L0880 | mass.f90:65-192 | MASS_FINITE_DIFFERENCES | 4.832 | mass | wall gather; pointer scratch/array constructors (U_TEMP, F_TEMP, Z_TEMP); function callees (GET_SCALAR_FACE_VALUE_PT); off-wall face writes | none (finished, uncommitted) | in progress | claimed | Legacy Mapper | - | I |
| 3 | L0401 | divg.f90:995-1088 | SPECIES_ADVECTION_PART_1_NEW | 3.260 | species/combustion | wall gather; pointer scratch/array constructors; function callees (GET_SCALAR_FACE_VALUE_PT); Z_TEMP pad of the fourth element | pointer scratch (S2) | claimed | claimed | Legacy Mapper | - | I |
| 4 | L1488 | wall.f90:144-185 | WALL_BC | 2.688 | wall/BC | wall gather; function callees (SURFACE_HEAT_TRANSFER, CALCULATE_ZZ_F, ...); ragged per-wall table (BOUNDARY_ONE_D); scratch-sum reference (OBSTRUCTION%MASS, D_SOURCE) | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 5 | L0369 | divg.f90:287-421 | DIVERGENCE_PART_1 | 2.399 | other | species loop N outermost; cell loop nests | generator gap | in progress | claimed | GPU Generator Engineer | - | R |
| 6 | L0877 | mass.f90:868-939 | CHECK_MASS_DENSITY | 1.967 | mass | species loop N outermost; scatter into 7 cells (two-pass gather in source order); CYCLE; integer max-reduction CLIP_RHO_ZZ(N) | race-uniqueness | claimed | claimed | Legacy Mapper | - | I |
| 7 | L1485 | wall.f90:888-947 | SURFACE_HEAT_TRANSFER | 1.832 | solid phase | wall gather; neighbour-mesh (OMESH copies through EWC%IIO_MIN..KKO_MAX); ragged per-wall table | neighbour mesh (T4) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 8 | L0823 | init.f90:4940-5054 | REASSIGN_WALL_CELLS | 1.727 | solid phase | host-side event code (burn-away, control-driven creation or removal); allocation and derived-type copies | host-side event code | claimed | host-side | AMR Solid Phase Lead | - | I |
| 9 | L0381 | divg.f90:640-660 | DIVERGENCE_PART_1 | 1.491 | other | species loop N outermost; function callees (SPECIES_ADVECTION_PART_2, GET_SENSIBLE_ENTHALPY_Z); guarded call SET_EXIMRHOZZLIM_3D | generator gap | claimed | claimed | GPU Generator Engineer | - | R |
| 10 | L1355 | velo.f90:204-252 | COMPUTE_VISCOSITY | 1.467 | velocity | cell loop; private local arrays A_IJ(3,3), B_IJ(3,3) with inner 3x3 loops; live-out scalars (PRIVATE list); module scalar C_VREMAN as argument | generator gap | claimed | claimed | Legacy Mapper | - | R |
| 11 | L1489 | wall.f90:192-194 | WALL_BC | 1.370 | wall/BC | thin-wall loop (not a wall loop); function callees (SOLID_HEAT_TRANSFER); ragged per-wall table (BOUNDARY_ONE_D); one shared event (OBSTRUCTION%MASS=-1 for a burned-through thin wall) | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | R |
| 12 | L1272 | soot.f90:54-175 | SETTLING_VELOCITY | 1.325 | species/combustion | species loop N outermost; per-species table (SPECIES_MIXTURE(N)%...); function callees (GET_VISCOSITY, GET_CONDUCTIVITY, CUNNINGHAM); wall gather (constant stores); WORK7..WORK9 aliasing | table needed | open | needs-feature | open — available | AMR Species & Combustion Lead | I |
| 13 | L1363 | velo.f90:2891-3015 | MATCH_VELOCITY_FLUX | 0.994 | wall/BC | wall gather; neighbour-mesh (OMESH) | neighbour mesh (T4) | claimed | claimed | GPU Mesh Data Loops Engineer | - | I |
| 14 | L1367 | velo.f90:1858-1897 | VELOCITY_BC | 0.795 | wall/BC | wall gather; neighbour-mesh (OMESH US/VS/WS) | neighbour mesh (T4) | claimed | claimed | GPU Mesh Data Loops Engineer | - | I |
| 15 | L1399 | velo.f90:514-550 | VISCOSITY_BC | 0.710 | wall/BC | wall gather; neighbour-mesh (OMESH MU) | neighbour mesh (T4) | claimed | claimed | GPU Mesh Data Loops Engineer | - | I |
| 16 | L0384 | divg.f90:696-716 | DIVERGENCE_PART_1 | 0.658 | other | species inner loop; source function (VD2D_MMS_Z_SRC); pointer SM | test-only | not planned (test-only) | test-only | GPU Generator Engineer | - | R |
| 17 | L0876 | mass.f90:799-849 | CHECK_MASS_DENSITY | 0.599 | mass | K,J nest; scatter into 7 cells (two-pass gather) | race-uniqueness | claimed | claimed | Legacy Mapper | - | I |
| 18 | L1305 | turb.f90:1796-1861 | SYNTHETIC_TURBULENCE | 0.572 | velocity | cell/face loop; non-perfect nest (outer loop) | table needed | open | needs-feature | open — available | - | I |
| 19 | L1470 | wall.f90:3635-3677 | HT3D_TEMPERATURE_EXCHANGE | 0.500 | solid phase | wall gather; ragged per-wall table (BOUNDARY_THR_D); store into TMP | neighbour mesh (T4) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 20 | L1471 | wall.f90:3681-3717 | HT3D_TEMPERATURE_EXCHANGE | 0.465 | solid phase | wall gather; ragged per-wall table (BOUNDARY_THR_D) | neighbour mesh (T4) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 21 | L1364 | velo.f90:1376-1400 | NO_FLUX | 0.426 | pressure | wall gather; neighbour-mesh (OMESH(NOM)%H/HS); overlap sum in K,J,I order then divide | neighbour mesh (T4) | claimed | claimed | GPU Mesh Data Loops Engineer | - | R |
| 22 | L1486 | wall.f90:109-119 | WALL_BC | 0.407 | wall/BC | wall gather; function callees (ASSIGN_GHOST_VALUE, NEAR_SURFACE_GAS_VARIABLES) | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 23 | L1487 | wall.f90:124-137 | WALL_BC | 0.354 | wall/BC | wall gather; function callees | callee flatten (T3) | claimed | claimed | AMR Solid Phase Lead | - | I |
| 24 | L1365 | velo.f90:1404-1459 | NO_FLUX | 0.347 | pressure | obstruction box loop (outer N over OBSTRUCTION); CELL%SOLID mask on both cells; pointer re-aim HP=>H/HS (hoist); PREDICTOR switch as argument | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 25 | L1389 | velo.f90:1048-1054 | MMS_VELOCITY_FLUX | 0.325 | velocity | face loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 26 | L1388 | velo.f90:1040-1046 | MMS_VELOCITY_FLUX | 0.315 | velocity | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 27 | L0822 | init.f90:4917-4936 | REASSIGN_WALL_CELLS | 0.310 | solid phase | host-side event code; neighbour-mesh reads | host-side event code | open | host-side | open — available | - | I |
| 28 | L0882 | mass.f90:224-320 | MASS_FINITE_DIFFERENCES | 0.302 | mass | wall gather; pointer scratch/array constructors; whole-array assignment in the loop body; off-wall face writes | none (finished, uncommitted) | in progress | claimed | Legacy Mapper | - | I |
| 29 | L1315 | turb.f90:979-995 | TEST_FILTER | 0.284 | velocity | cell/face loop; live-out scalars (PRIVATE list) | generator gap | claimed | claimed | Legacy Mapper | - | I |
| 30 | L1316 | turb.f90:1003-1019 | TEST_FILTER | 0.284 | velocity | cell/face loop; live-out scalars (PRIVATE list) | generator gap | claimed | claimed | Legacy Mapper | - | I |
| 31 | L1317 | turb.f90:1027-1043 | TEST_FILTER | 0.284 | velocity | cell/face loop; rank/subscript contract | generator gap | open | needs-feature | open — available | - | I |
| 32 | L0398 | divg.f90:835-939 | ENTHALPY_ADVECTION_NEW | 0.260 | species/combustion | wall gather; pointer scratch/array constructors; CSR cell-to-wall gather (ascending IW); pointer remap U_TEMP=>U_WORK | pointer scratch (S2) | claimed | claimed | Legacy Mapper | - | I |
| 33 | L1452 | wall.f90:319-339 | ASSIGN_GHOST_VALUE | 0.199 | wall/BC | wall gather; neighbour-mesh | neighbour mesh (T4) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 34 | L0403 | divg.f90:1120-1183 | SPECIES_ADVECTION_PART_1_NEW | 0.193 | species/combustion | wall gather; pointer scratch/array constructors; whole-array assignment in the loop body | pointer scratch (S2) | claimed | claimed | Legacy Mapper | - | I |
| 35 | L1356 | velo.f90:256-279 | COMPUTE_VISCOSITY | 0.187 | velocity | cell loop; function callees (WALE_VISCOSITY, array dummy A_IJ); module scalar C_WALE as argument | generator gap | claimed | claimed | Legacy Mapper | - | R |
| 36 | L0861 | mass.f90:484-495 | DENSITY | 0.171 | mass | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 37 | L0870 | mass.f90:666-677 | DENSITY | 0.171 | mass | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 38 | L1321 | turb.f90:501-534 | VARDEN_DYNSMAG | 0.105 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | claimed | claimed | Legacy Mapper | - | I |
| 39 | L0878 | mass.f90:947-961 | CHECK_MASS_DENSITY | 0.100 | mass | cell loop; reductions (SUM, MAXLOC per cell); CYCLE on solid; array-section update | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 40 | L0364 | divg.f90:93-95 | DIVERGENCE_PART_1 | 0.099 | other | wall loop; wall gather; function callees; module scalar as argument | callee flatten (T3) | claimed | claimed | GPU Generator Engineer | - | I |
| 41 | L1360 | velo.f90:413-435 | COMPUTE_STRAIN_RATE | 0.089 | velocity | cell loop; live-out scalars (PRIVATE list) | generator gap | open | needs-feature | open — available | - | R |
| 42 | L1121 | pres.f90:4012-4026 | CHECK_UNSUPPORTED_MESH | 0.086 | pressure | mesh/rank loop; MY_RANK test | blocked | open | blocked | open — available | AMR Pressure Solver Lead | R |
| 43 | L0397 | divg.f90:814-822 | ENTHALPY_ADVECTION_NEW | 0.077 | species/combustion | cell loop; bounds -1:IBP1+1 (policy.arrays); function callees (GET_SENSIBLE_ENTHALPY); array-section copy ZZ_GET(1:N)=ZZP(I,J,K,1:N); DOT_PRODUCT over a table section | generator gap | open | needs-feature | open — available | AMR Species & Combustion Lead | R |
| 44 | L0865 | mass.f90:556-564 | DENSITY | 0.071 | mass | cell loop; CYCLE on solid; function callees (GET_SPECIFIC_GAS_CONSTANT); array-section copy | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 45 | L0874 | mass.f90:738-746 | DENSITY | 0.071 | mass | cell loop; CYCLE on solid; function callees (GET_SPECIFIC_GAS_CONSTANT); array-section copy ZZ_GET(1:N)=ZZ(I,J,K,1:N) | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 46 | L0881 | mass.f90:201-209 | MASS_FINITE_DIFFERENCES | 0.071 | mass | cell loop; bounds -1:IBP1+1; function callees (GET_MOLECULAR_WEIGHT); array-section copy | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 47 | L0371 | divg.f90:463-471 | DIVERGENCE_PART_1 | 0.068 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 48 | L1309 | turb.f90:1323-1345 | TENSOR_DIFFUSIVITY_MODEL | 0.068 | velocity | cell/face loop; live-out scalars (PRIVATE list); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 49 | L1310 | turb.f90:1347-1369 | TENSOR_DIFFUSIVITY_MODEL | 0.068 | velocity | cell/face loop; live-out scalars (PRIVATE list); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 50 | L1311 | turb.f90:1371-1393 | TENSOR_DIFFUSIVITY_MODEL | 0.068 | velocity | cell/face loop; live-out scalars (PRIVATE list); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 51 | L1312 | turb.f90:1414-1434 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | velocity | cell/face loop; live-out scalars (PRIVATE list) | generator gap | open | needs-feature | open — available | - | I |
| 52 | L1313 | turb.f90:1436-1456 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | velocity | cell/face loop; live-out scalars (PRIVATE list) | generator gap | open | needs-feature | open — available | - | I |
| 53 | L1314 | turb.f90:1458-1479 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 54 | L0402 | divg.f90:1097-1105 | SPECIES_ADVECTION_PART_1_NEW | 0.062 | species/combustion | cell loop; bounds -1:IBP1+1; function callees (GET_MOLECULAR_WEIGHT); array-section copy | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 55 | L1358 | velo.f90:306-351 | COMPUTE_VISCOSITY | 0.056 | velocity | wall loop; wall gather; function callees; module scalar as argument | table needed | claimed | claimed | GPU Wall Loops Engineer | - | I |
| 56 | L1366 | velo.f90:1463-1559 | NO_FLUX | 0.056 | pressure | wall gather; face write through wall subscripts; EXTERNAL_WALL(IW)%NOM designator (EW_NOM); B1 alias U_NORMAL/U_NORMAL_S; module scalar PRES_FLAG as argument; SELECT on IOR | race-uniqueness | in progress | claimed | GPU Mesh Data Loops Engineer | - | R |
| 57 | L1372 | velo.f90:1770-1778 | VELOCITY_CORRECTOR | 0.049 | pressure | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 58 | L1373 | velo.f90:1779-1787 | VELOCITY_CORRECTOR | 0.049 | pressure | face loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 59 | L1397 | velo.f90:1648-1656 | VELOCITY_PREDICTOR | 0.049 | pressure | cell loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 60 | L1398 | velo.f90:1657-1665 | VELOCITY_PREDICTOR | 0.049 | pressure | face loop; source function | test-only | not planned (test-only) | test-only | open — available | - | R |
| 61 | L0379 | divg.f90:608-616 | DIVERGENCE_PART_1 | 0.037 | other | cell/face loop; rank/subscript contract | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 62 | L1275 | turb.f90:819-831 | CALC_VARDEN_LEONARD_TERM | 0.037 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 63 | L1323 | turb.f90:633-645 | VARDEN_DYNSMAG | 0.037 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 64 | L0866 | mass.f90:570-577 | DENSITY | 0.029 | mass | cell loop; CYCLE on solid; rank-2 gather PBAR_S(K,PRESSURE_ZONE(I,J,K)) | generator gap | open | needs-feature | open — available | AMR Species & Combustion Lead | R |
| 65 | L0872 | mass.f90:709-716 | DENSITY | 0.029 | mass | cell loop; CYCLE on solid; array-section update ZZ(I,J,K,1:NS) | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 66 | L0875 | mass.f90:752-759 | DENSITY | 0.029 | mass | cell loop; CYCLE on solid; rank-2 gather PBAR(K,PRESSURE_ZONE(I,J,K)) | generator gap | open | needs-feature | open — available | AMR Species & Combustion Lead | R |
| 67 | L0879 | mass.f90:980-987 | CLIP_PASSIVE_SCALARS | 0.029 | mass | cell loop; CYCLE on solid; module integer ZETA_INDEX as argument | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 68 | L1391 | velo.f90:1243-1253 | VELOCITY_FLUX_CYLINDRICAL | 0.026 | velocity | cell loop; pointer aliases to WORK2, WORK5 | none | open | translatable-now | open — available | - | R |
| 69 | L0377 | divg.f90:573-582 | DIVERGENCE_PART_1 | 0.025 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 70 | L0380 | divg.f90:622-629 | DIVERGENCE_PART_1 | 0.025 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 71 | L0386 | divg.f90:757-767 | DIVERGENCE_PART_1 | 0.022 | other | wall loop; wall gather; live-out scalars (PRIVATE list); zone table | table needed | claimed | claimed | GPU Generator Engineer | - | I |
| 72 | L1325 | turb.f90:690-710 | VARDEN_DYNSMAG | 0.021 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 73 | L1354 | velo.f90:187-196 | COMPUTE_VISCOSITY | 0.021 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 74 | L0630 | func.f90:5509-5524 | BLOCK_CELL | 0.020 | other | cell/face loop; wall gather; derived-type designator table | neighbour mesh (T4) | open | blocked | open — available | - | I |
| 75 | L0817 | init.f90:4878-4889 | CREATE_OR_REMOVE_OBST | 0.017 | solid phase | cell/face loop; derived-type designator table | generator gap | open | needs-feature | open — available | AMR Solid Phase Lead | I |
| 76 | L0677 | func.f90:5183-5201 | PACK_CELL | 0.016 | other | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | - | I |
| 77 | L1379 | velo.f90:970-978 | CORIOLIS_FORCE | 0.016 | velocity | cell loop; output views WORK7..WORK9 (pointer aliases) | none | open | translatable-now | open — available | - | R |
| 78 | L1381 | velo.f90:992-1000 | CORIOLIS_FORCE | 0.016 | velocity | cell loop; constant-subscript rank-1 array OVEC(n) as by-value scalars | generator gap | open | needs-feature | open — available | - | R |
| 79 | L1382 | velo.f90:1006-1014 | CORIOLIS_FORCE | 0.016 | velocity | cell loop; constant-subscript rank-1 array OVEC(n) | generator gap | open | needs-feature | open — available | - | R |
| 80 | L1383 | velo.f90:1020-1028 | CORIOLIS_FORCE | 0.016 | velocity | face loop; constant-subscript rank-1 array OVEC(n) | generator gap | open | needs-feature | open — available | - | R |
| 81 | L0911 | part.f90:4643-4763 | PARTICLE_MOMENTUM_TRANSFER | 0.015 | other | wall loop; wall gather; neighbour-mesh; pointer alias | neighbour mesh (T4) | open | blocked | open — available | - | I |
| 82 | L1392 | velo.f90:1270-1302 | VELOCITY_FLUX_CYLINDRICAL | 0.014 | velocity | edge tables; K,I nest with J fixed | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 83 | L1393 | velo.f90:1306-1337 | VELOCITY_FLUX_CYLINDRICAL | 0.014 | velocity | edge tables; K,I nest with J fixed | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 84 | L0373 | divg.f90:499-505 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 85 | L0378 | divg.f90:591-597 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 86 | L0382 | divg.f90:668-674 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop | none | claimed | claimed | GPU Generator Engineer | - | I |
| 87 | L0383 | divg.f90:681-687 | DIVERGENCE_PART_1 | 0.012 | other | cell/face loop; rank/subscript contract | generator gap | claimed | claimed | GPU Generator Engineer | - | I |
| 88 | L1276 | turb.f90:187-194 | COMPRESSION_WAVE | 0.011 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 89 | L1277 | turb.f90:195-202 | COMPRESSION_WAVE | 0.011 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 90 | L1278 | turb.f90:203-210 | COMPRESSION_WAVE | 0.011 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 91 | L1279 | turb.f90:211-218 | COMPRESSION_WAVE | 0.011 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 92 | L1324 | turb.f90:666-678 | VARDEN_DYNSMAG | 0.011 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 93 | L1385 | velo.f90:896-903 | DIRECT_FORCE | 0.011 | velocity | face loop; constant-subscript rank-1 array FVEC(n) as by-value scalar; ramp factor computed on the host | generator gap | open | needs-feature | open — available | - | R |
| 94 | L1386 | velo.f90:917-924 | DIRECT_FORCE | 0.011 | velocity | face loop; FVEC(2) as scalar | generator gap | open | needs-feature | open — available | - | R |
| 95 | L1387 | velo.f90:938-945 | DIRECT_FORCE | 0.011 | velocity | face loop; FVEC(3) as scalar | generator gap | open | needs-feature | open — available | - | R |
| 96 | L0890 | part.f90:605-607 | INSERT_VENT_PARTICLES | 0.008 | other | wall loop; wall gather; function callees; module scalar as argument | callee flatten (T3) | open | blocked | open — available | - | I |
| 97 | L1209 | pres.f90:65-228 | PRESSURE_SOLVER_COMPUTE_RHS | 0.008 | pressure | wall gather; neighbour-mesh (MESHES(NOM)%DX); function callees (EVALUATE_RAMP); vent table (VENTS(WC%VENT_INDEX)); wall-keyed 2-D outputs | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | R |
| 98 | L1224 | pres.f90:544-567 | TUNNEL_POISSON_SOLVER | 0.007 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 99 | L1280 | turb.f90:221-227 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 100 | L1281 | turb.f90:228-234 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 101 | L1282 | turb.f90:235-241 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 102 | L1283 | turb.f90:242-248 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 103 | L1284 | turb.f90:250-256 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 104 | L1285 | turb.f90:257-263 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 105 | L1286 | turb.f90:264-270 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 106 | L1287 | turb.f90:271-277 | COMPRESSION_WAVE | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 107 | L1322 | turb.f90:570-581 | VARDEN_DYNSMAG | 0.005 | velocity | cell/face loop; bounds: policy.arrays line | generator gap | open | needs-feature | open — available | - | I |
| 108 | L1351 | velo.f90:115-121 | COMPUTE_VISCOSITY | 0.005 | velocity | cell/face loop | none | open | translatable-now | open — available | - | I |
| 109 | L1204 | pres.f90:1931-1994 | ULMAT_SOLVE_ZONE | 0.003 | pressure | wall loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 110 | L0596 | fire.f90:1902-1918 | COMBUSTION_BC | 0.002 | species/combustion | wall loop; wall gather; neighbour-mesh; derived-type designator table | neighbour mesh (T4) | open | blocked | open — available | AMR Species & Combustion Lead | I |
| 111 | L1198 | pres.f90:1837-1846 | ULMAT_SOLVE_ZONE | 0.002 | pressure | cell/face loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 112 | L1211 | pres.f90:250-260 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop | none | open | translatable-now | open — available | AMR Pressure Solver Lead | R |
| 113 | L1212 | pres.f90:267-277 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop; transposed output subscripts | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 114 | L1213 | pres.f90:282-292 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop; transposed output subscripts | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 115 | L1214 | pres.f90:297-307 | PRESSURE_SOLVER_COMPUTE_RHS | 0.002 | pressure | cell loop; transposed output subscripts | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 116 | L1245 | radi.f90:5044-5059 | RADIATION_FVM | 0.002 | radiation | file output | host (I/O) | not planned (host) | host-side | AMR Radiation Lead | - | I |
| 117 | L1273 | soot.f90:412-472 | SOOT_SURFACE_OXIDATION | 0.002 | species/combustion | wall loop; wall gather; pointer alias | neighbour mesh (T4) | open | blocked | open — available | AMR Species & Combustion Lead | I |
| 118 | L0680 | func.f90:2845-2847 | COMPUTE_WIND_COMPONENTS | 0.001 | other | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | - | I |
| 119 | L1130 | pres.f90:5580-5601 | GET_H_REGFACES | 0.001 | pressure | wall gather; rank-4 LOGICAL table | table needed | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 120 | L1154 | pres.f90:2031-2047 | ULMAT_GET_H_REGFACES | 0.001 | pressure | wall gather; rank-4 LOGICAL table | table needed | open | needs-feature | open — available | AMR Pressure Solver Lead | I |
| 121 | L1190 | pres.f90:1686-1694 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 122 | L1192 | pres.f90:1713-1720 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 123 | L1194 | pres.f90:1782-1790 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; derived-type designator table; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 124 | L1196 | pres.f90:1808-1815 | ULMAT_SOLVE_ZONE | 0.001 | pressure | cell/face loop; wall gather; zone table | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 125 | L1206 | pres.f90:729-742 | PRESSURE_SOLVER_CHECK_RESIDUALS | 0.001 | pressure | cell loop; bounds (WORK8 view, policy.arrays) | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 126 | L1208 | pres.f90:768-787 | PRESSURE_SOLVER_CHECK_RESIDUALS | 0.001 | pressure | cell loop; bounds (WORK8 view); reduction by the caller (MAXVAL, MAXLOC) | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 127 | L1248 | radi.f90:3651-3662 | INTERPOLATE_IL | 0.001 | radiation | wall loop; ragged per-wall table BR_ILW(NRA,bands,wall); ILW_OLD copy kept (ping-pong table or per-thread scratch); I_INTP accumulation ascending and serial per angle | table needed | claimed | claimed | AMR Radiation Lead | - | I |
| 128 | L1318 | turb.f90:1054-1059 | TEST_FILTER | 0.001 | velocity | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | - | I |
| 129 | L1319 | turb.f90:1063-1068 | TEST_FILTER | 0.001 | velocity | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | - | I |
| 130 | L1320 | turb.f90:1072-1077 | TEST_FILTER | 0.001 | velocity | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | - | I |
| 131 | L1375 | velo.f90:649-653 | VELOCITY_FLUX | 0.001 | velocity | 1-D loop 0:IBAR; function callees (EVALUATE_RAMP on a ramp table) | generator gap | claimed | claimed | GPU Generator Engineer | - | R |
| 132 | L0608 | fire.f90:1957-1961 | CONDENSATION_EVAPORATION | 0.000 | species/combustion | wall loop; constant store through B1 alias (idempotent) | none | claimed | claimed | AMR Species & Combustion Lead | - | R |
| 133 | L0884 | part.f90:4606-4617 | CLIP_PARTICLE_DRAG | 0.000 | other | cell/face loop | none | open | translatable-now | open — available | - | I |
| 134 | L1203 | pres.f90:1923-1925 | ULMAT_SOLVE_ZONE | 0.000 | pressure | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | I |
| 135 | L1207 | pres.f90:758-764 | PRESSURE_SOLVER_CHECK_RESIDUALS | 0.000 | pressure | cell loop | none | open | translatable-now | open — available | AMR Pressure Solver Lead | R |
| 136 | L1210 | pres.f90:238-245 | PRESSURE_SOLVER_COMPUTE_RHS | 0.000 | pressure | cell loop (K,I with J=1); non-perfect nest | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 137 | L1215 | pres.f90:394-400 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop | none | open | translatable-now | open — available | AMR Pressure Solver Lead | R |
| 138 | L1216 | pres.f90:404-410 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop; transposed input subscripts | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 139 | L1217 | pres.f90:414-420 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop; transposed input subscripts | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 140 | L1218 | pres.f90:424-430 | PRESSURE_SOLVER_FFT | 0.000 | pressure | cell loop; transposed input subscripts | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 141 | L1219 | pres.f90:438-440 | PRESSURE_SOLVER_FFT | 0.000 | pressure | non-perfect nest (section assignment); 1-D table H_BAR(I_OFFSET+I) | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 142 | L1220 | pres.f90:450-462 | PRESSURE_SOLVER_FFT | 0.000 | pressure | K,J-only nest; integer BC-code scalars as arguments; reads BXS/BXF | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 143 | L1221 | pres.f90:466-477 | PRESSURE_SOLVER_FFT | 0.000 | pressure | K,I-only nest; integer BC-code scalars as arguments | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 144 | L1222 | pres.f90:481-492 | PRESSURE_SOLVER_FFT | 0.000 | pressure | J,I-only nest; integer BC-code scalars as arguments | generator gap | open | needs-feature | open — available | AMR Pressure Solver Lead | R |
| 145 | L1223 | pres.f90:531-537 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 146 | L1225 | pres.f90:576-582 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 147 | L1226 | pres.f90:583-589 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 148 | L1227 | pres.f90:591-596 | TUNNEL_POISSON_SOLVER | 0.000 | pressure | cell/face loop; wall gather; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | AMR Pressure Solver Lead | I |
| 149 | L1243 | radi.f90:4961-4970 | RADIATION_FVM | 0.000 | radiation | wall loop; ragged per-wall table BR_ILW(NRA,bands,wall); per-band temporary T=0; T=T+ILW(N); Q=Q+T (flat chain is not bitwise); OPEN_BOUNDARY and IW<=N_EXTERNAL_WALL_CELLS conditions | table needed | claimed | claimed | AMR Radiation Lead | - | I |
| 150 | L1289 | turb.f90:882-884 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 151 | L1290 | turb.f90:887-889 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 152 | L1291 | turb.f90:892-894 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 153 | L1292 | turb.f90:897-899 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 154 | L1293 | turb.f90:904-906 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 155 | L1294 | turb.f90:909-911 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 156 | L1295 | turb.f90:914-916 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 157 | L1296 | turb.f90:919-921 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 158 | L1297 | turb.f90:926-928 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 159 | L1298 | turb.f90:931-933 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 160 | L1299 | turb.f90:936-938 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 161 | L1300 | turb.f90:941-943 | FILL_EDGES | 0.000 | velocity | cell/face loop; non-perfect nest (outer loop); reductions/CYCLE | generator gap | open | needs-feature | open — available | - | I |
| 162 | L1330 | vege.f90:618-623 | GET_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 163 | L1331 | vege.f90:624-629 | GET_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 164 | L1332 | vege.f90:630-639 | GET_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; function callees; non-perfect nest (outer loop) | callee flatten (T3) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 165 | L1333 | vege.f90:678-699 | FILL_BOUNDARY_VALUES | 0.000 | solid phase | cell/face loop; wall gather; neighbour-mesh; non-perfect nest (outer loop) | neighbour mesh (T4) | open | blocked | open — available | AMR Solid Phase Lead | I |
| 166 | L1334 | vege.f90:885-894 | LEVEL_SET_ADVECT_FLUX | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | AMR Solid Phase Lead | I |
| 167 | L1335 | vege.f90:896-905 | LEVEL_SET_ADVECT_FLUX | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | AMR Solid Phase Lead | I |
| 168 | L1336 | vege.f90:907-914 | LEVEL_SET_ADVECT_FLUX | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | AMR Solid Phase Lead | I |
| 169 | L1340 | vege.f90:763-857 | LEVEL_SET_SPREAD_RATE | 0.000 | solid phase | cell/face loop; non-perfect nest (outer loop) | generator gap | open | needs-feature | open — available | AMR Solid Phase Lead | I |
<!-- GENERATED-END:worklist -->

## 6. Species & Combustion

Scope: `species.f90` and the combustion modules (the survey records no loop there), `fire.f90` (two eligible loops, L0596 and L0608; five more are geometry-deferred), `soot.f90` (L1272 and L1273 eligible; one geometry-deferred), `mass.f90`, and the species parts of `divg.f90` (the S2, S3, S4 loops and the CHECK_MASS_DENSITY gather pieces). Proposed owner of every open loop in this section: **AMR Species & Combustion Lead**, working in own files and own tests only; they message the owner of a loop or a function before any shared commit (`shared:` rule above). The loops already claimed stay with the owners shown.

How the families map to this section (details and sign-offs in `blocked-loop-families.md`):
- **S1** (CHECK_MASS_DENSITY scatter, L0877 and L0876): the scatter writes one cell and its six neighbours, so it stays on the host until the two-pass gather exists. The pieces: pass 1, per source cell, stores the seven contributions; pass 2, per target cell, adds them in the serial source order `(I,J,K-1)`, `(I,J-1,K)`, `(I-1,J,K)`, self, `(I+1,J,K)`, `(I,J+1,K)`, `(I,J,K+1)`, interior targets only, with the same `SUM(MASS_N)` order; `CLIP_RHO_ZZ(N)` becomes an integer max-reduction; the test must compile host and gather with matching no-FMA flags. L0877 and L0876 are held by the Legacy Mapper (wall split above). The pieces that are already translated (`delta_rho_zz_zero`, `rho_zz_clip_assign`) count; the scatter credit is not claimed.
- **S2** (wall nests with pointer scratch and array constructors: L0880, L0401, L0882, L0403, L0398): the Legacy Mapper's package. The four finished mass-flux wall nests cover the L0880 and L0882 parts; the divergence-side nests L0401 and L0403 reuse the same nests and need the fourth `Z_TEMP` element padded with 0 (MP5 bit tests wait for upstream patches 0001 and 0002).
- **S3** (wall scatter into `U_DOT_DEL_RHO_*`, cell loops reading `CELL%WALL_INDEX`: L0398, L0399, L0405, L0406): L0398 stays with the Legacy Mapper; L0405 is translated by the GPU Wall Loops Engineer (`wall_spec_adv2`); L0399 and L0406 are open and need the six-flag mask table (`WALL_INDEX`), which exists as a `cell_int` table on the s5-wall branch and is not merged.
- **S4** (SETTLING_VELOCITY, L1272): open, needs a per-species table and a split into three kernels.
- Plain cell loops of `mass.f90` and `divg.f90` (L0865, L0874, L0872, L0879, L0878, L0881, L0402): the front end accepts them; what is missing is the sidecar entry and the bitwise test.
- Test-only: L0861 and L0870 (`IF (PERIODIC_TEST==7)`, guards at `mass.f90:483` and `mass.f90:665`) and L0384 (`divg.f90:693`).
- Not in this section: the DIVERGENCE_PART_1 nests (L0365, L0369, L0366, L0381) belong to the GPU Generator Engineer even though they handle species.

Classification of each loop (translatable-now / needs-feature / blocked / claimed) follows. Classes are the work classes of section 2.

<!-- GENERATED-BEGIN:species -->
23 not-yet-translated loops in the species/combustion and mass domains, 13.652 % modelled share.


#### claimed (15 loops, 11.846 %)

| Loop | file:lines | Routine | Share % | Features / blocker | Closest kernel | Status | Owner / proposed |
|---|---|---|---|---|---|---|---|
| L0880 | mass.f90:65-192 | MASS_FINITE_DIFFERENCES | 4.832 | wall gather; pointer scratch/array constructors (U_TEMP, F_TEMP, Z_TEMP); function callees (GET_SCALAR_FACE_VALUE_PT); off-wall face writes / S2: four wall nests (spec_wall_zz, spec_wall_rmw, mass_wall_zz, mass_wall_rmw) finished by the Legacy Mapper, not yet committed to s5-gen | gsfv_* and rho_z_p_mass (cell nest) | in progress | Legacy Mapper |
| L0401 | divg.f90:995-1088 | SPECIES_ADVECTION_PART_1_NEW | 3.260 | wall gather; pointer scratch/array constructors; function callees (GET_SCALAR_FACE_VALUE_PT); Z_TEMP pad of the fourth element / S2; MP5 bit tests wait for upstream patches 0001/0002 | gsfv_*, rho_z_p_divg | claimed | Legacy Mapper |
| L0877 | mass.f90:868-939 | CHECK_MASS_DENSITY | 1.967 | species loop N outermost; scatter into 7 cells (two-pass gather in source order); CYCLE; integer max-reduction CLIP_RHO_ZZ(N) / S1: the scatter stays on the host until the two-pass gather exists | rho_zz_clip_assign, delta_rho_zz_zero | claimed | Legacy Mapper |
| L0876 | mass.f90:799-849 | CHECK_MASS_DENSITY | 0.599 | K,J nest; scatter into 7 cells (two-pass gather) / S1 | rho_zz_clip_assign | claimed | Legacy Mapper |
| L0882 | mass.f90:224-320 | MASS_FINITE_DIFFERENCES | 0.302 | wall gather; pointer scratch/array constructors; whole-array assignment in the loop body; off-wall face writes / S2: covered by the finished mass_wall_* nests, not yet committed | gsfv_* | in progress | Legacy Mapper |
| L0398 | divg.f90:835-939 | ENTHALPY_ADVECTION_NEW | 0.260 | wall gather; pointer scratch/array constructors; CSR cell-to-wall gather (ascending IW); pointer remap U_TEMP=>U_WORK / S2 and S3: accumulate through the wall gas-cell subscript | gsfv_*, CSR family (not built) | claimed | Legacy Mapper |
| L0403 | divg.f90:1120-1183 | SPECIES_ADVECTION_PART_1_NEW | 0.193 | wall gather; pointer scratch/array constructors; whole-array assignment in the loop body / S2 | gsfv_* | claimed | Legacy Mapper |
| L0878 | mass.f90:947-961 | CHECK_MASS_DENSITY | 0.100 | cell loop; reductions (SUM, MAXLOC per cell); CYCLE on solid; array-section update / front end accepts; no bitwise test committed | flux_mw_fix, rho_sum | claimed | AMR Species & Combustion Lead |
| L0865 | mass.f90:556-564 | DENSITY | 0.071 | cell loop; CYCLE on solid; function callees (GET_SPECIFIC_GAS_CONSTANT); array-section copy / front end accepts; no bitwise test committed | rsum_pred (mass.f90:527-534), mu_dns | claimed | AMR Species & Combustion Lead |
| L0874 | mass.f90:738-746 | DENSITY | 0.071 | cell loop; CYCLE on solid; function callees (GET_SPECIFIC_GAS_CONSTANT); array-section copy ZZ_GET(1:N)=ZZ(I,J,K,1:N) / front end accepts; no bitwise test committed | rsum_pred (mass.f90:527-534), mu_dns | claimed | AMR Species & Combustion Lead |
| L0881 | mass.f90:201-209 | MASS_FINITE_DIFFERENCES | 0.071 | cell loop; bounds -1:IBP1+1; function callees (GET_MOLECULAR_WEIGHT); array-section copy / front end accepts; no bitwise test committed | mu_dns (private ZZ_GET, table callee) | claimed | AMR Species & Combustion Lead |
| L0402 | divg.f90:1097-1105 | SPECIES_ADVECTION_PART_1_NEW | 0.062 | cell loop; bounds -1:IBP1+1; function callees (GET_MOLECULAR_WEIGHT); array-section copy / front end accepts; no bitwise test committed | mu_dns, rho_z_p_divg (divg.f90:998-1004) | claimed | AMR Species & Combustion Lead |
| L0872 | mass.f90:709-716 | DENSITY | 0.029 | cell loop; CYCLE on solid; array-section update ZZ(I,J,K,1:NS) / front end accepts; no bitwise test committed | zz_corr (mass.f90:619-632), rho_sum | claimed | AMR Species & Combustion Lead |
| L0879 | mass.f90:980-987 | CLIP_PASSIVE_SCALARS | 0.029 | cell loop; CYCLE on solid; module integer ZETA_INDEX as argument / front end accepts; no bitwise test committed | rho_zz_clip_assign (mass.f90:931-937) | claimed | AMR Species & Combustion Lead |
| L0608 | fire.f90:1957-1961 | CONDENSATION_EVAPORATION | 0.000 | wall loop; constant store through B1 alias (idempotent) / front end accepts; no bitwise test committed | wall_b2_work1 (part.f90:3654-3659) | claimed | AMR Species & Combustion Lead |

#### needs-feature (4 loops, 1.460 %)

| Loop | file:lines | Routine | Share % | Features / blocker | Closest kernel | Status | Owner / proposed |
|---|---|---|---|---|---|---|---|
| L1272 | soot.f90:54-175 | SETTLING_VELOCITY | 1.325 | species loop N outermost; per-species table (SPECIES_MIXTURE(N)%...); function callees (GET_VISCOSITY, GET_CONDUCTIVITY, CUNNINGHAM); wall gather (constant stores); WORK7..WORK9 aliasing / S4: split into three kernels | mu_dns, wall_up_ghost | open | proposed: AMR Species & Combustion Lead |
| L0397 | divg.f90:814-822 | ENTHALPY_ADVECTION_NEW | 0.077 | cell loop; bounds -1:IBP1+1 (policy.arrays); function callees (GET_SENSIBLE_ENTHALPY); array-section copy ZZ_GET(1:N)=ZZP(I,J,K,1:N); DOT_PRODUCT over a table section / front end rejects the DOT_PRODUCT section H_SENS_Z(ITMP+1,1:NS) | mu_dns (ZZ_GET private, table callee) | open | proposed: AMR Species & Combustion Lead |
| L0866 | mass.f90:570-577 | DENSITY | 0.029 | cell loop; CYCLE on solid; rank-2 gather PBAR_S(K,PRESSURE_ZONE(I,J,K)) / rank-2 array PBAR_S indexed by a cell integer array | rho_sum | open | proposed: AMR Species & Combustion Lead |
| L0875 | mass.f90:752-759 | DENSITY | 0.029 | cell loop; CYCLE on solid; rank-2 gather PBAR(K,PRESSURE_ZONE(I,J,K)) / rank-2 array PBAR indexed by a cell integer array (inventory: rank contract) | rho_sum | open | proposed: AMR Species & Combustion Lead |

#### blocked (2 loops, 0.004 %)

| Loop | file:lines | Routine | Share % | Features / blocker | Closest kernel | Status | Owner / proposed |
|---|---|---|---|---|---|---|---|
| L0596 | fire.f90:1902-1918 | COMBUSTION_BC | 0.002 | wall loop; wall gather; neighbour-mesh; derived-type designator table / needs neighbour-mesh or ragged per-wall data (inventory) | none (new); pvf_* shows the host launch rule | open | proposed: AMR Species & Combustion Lead |
| L1273 | soot.f90:412-472 | SOOT_SURFACE_OXIDATION | 0.002 | wall loop; wall gather; pointer alias / needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open | proposed: AMR Species & Combustion Lead |

#### test-only (2 loops, 0.342 %)

| Loop | file:lines | Routine | Share % | Features / blocker | Closest kernel | Status | Owner / proposed |
|---|---|---|---|---|---|---|---|
| L0861 | mass.f90:484-495 | DENSITY | 0.171 | cell loop; source function / IF (PERIODIC_TEST==7) at mass.f90:483 | none | not planned (test-only) | open — available |
| L0870 | mass.f90:666-677 | DENSITY | 0.171 | cell loop; source function / IF (PERIODIC_TEST==7) at mass.f90:665 | none | not planned (test-only) | open — available |

Already translated in these domains: 20 loops, 1.585 %.
<!-- GENERATED-END:species -->

## 7. Regenerating

```
python3 tools/inventory/loop_work_list.py            # rewrite the generated blocks and loop_work_list.csv
python3 tools/inventory/loop_work_list.py --check    # fail if the committed files differ from a fresh run
python3 tools/inventory/loop_work_list.py --help
```

The script reads the survey CSV, `docs/amrex/loop_claims.csv` and the generator sidecar at the generator branch head (`--gen-repo`, `--gen-ref`). It never writes into the generator worktree. Options `--no-verify-lines` skips the check against the survey source.
