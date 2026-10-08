# Next work per role

Refreshed 2026-10-08 (after usage outage).

Owner rule: under the spend lock (CEO), only the active items below may proceed; everyone else is idle and must not take new claims. When the lock lifts, the standing rule returns: nobody sits idle; each role has an ordered list of ready items with a source reference (document path and item id). The Chief Architect owns these lists; the Spec & Program Lead refreshed an earlier version at spec v0.4.35; this is the Chief Architect's refresh (r7) after the usage outage. Where another list differs, this one wins. A role that finishes its active list tells the Chief Architect. "Dry risk" is the chance the role runs out of ready work within two days (low: five or more ready items; medium: three or four; high: fewer than three or all blocked). Items that wait for another role are counted as not ready when the level is set. Parked items below stay on the lists but are not active until the lock lifts.

Status as of spec v0.4.44 plus D-081 to D-089 and the local FireX merge into FDS-AMReX, next-work revision **r7** (spend lock; vortex owner-show abort + physics diagnosis; restart gate G1–G3; G1 PASS 2026-10-08 on Role 1 divfix, interim pe20 run at blocking 4/maxsize 8 launched, 8³ cube fix in progress). Prior r6 landed: GNU A-14 stretch `USE_AMREX=ON`+HDF5 closed; Species combustion/06 host ACCEPT; Radiation solid/08 R1–R4 and L1251 closed; Intel UP-0011 tip Rel+Dbg PASS; Mesh Data D-083 `[[ported]]` on s5-gen `08f6016e28`; Role 3 adaptive TAG_* green `370c4412cf` + Role 1 D-089 `8c8381b946`; Role 2 D-092+S6 accepted; Mesh Data L1305/L1317 claimed. Overnight owner-show MIRROR vortex aborted at step 1157/2000 (ENOSPC on plotfile s1158; physics already bad — see Active item 1).

Standing coding rule for implementers: `coding-rules.md` / D-084 (simplicity; minimum lines, FLOPs, conditionals, copies). Decisions since that note: D-085 (scheme-3 device `maxulp ≤ 8192`, named exception to D-070); D-089 (`FREEZE_VELOCITY` cold-flow exception for AMR adaptive demos).

Decisions since r5: D-083 (`CHECK_MASS_DENSITY`: GPU-only hand-written two-pass kernel in the K2 style, CPU loop unchanged, `UP-0003` tabled), D-082 (MP5 out of scope, `UP-0001` and `UP-0002` tabled; MP5-only kernels, markers, registry entries and tests are marked out of scope, not deleted, and do not count against coverage), D-081 (single-mesh FDS is the ground truth for step-sequence checks; multi-mesh first-step offset is a reference quirk). Earlier: D-080 (tag buffering, `MAXLOC` sign-offs, `WALL_LOOP_2` host path, test lock, `run_class_tolerance`, pressure conditions C1 to C3), D-079, D-078, D-077.

Landed since r5 (items removed from the lists, each checked in the repositories): the test-lock wrapper `test/with_s5lock.sh` with the 30-minute self-release and the lock rules in the generator README (s5-gen `ce4cb8704b`, `9f0b70924b`); the `rho_d_maxloc_fix` NaN limit in the kernel note and the tie test with three or more species that kills `maxloc_x_last_wins` and `maxloc_z_last_wins` (generator README, "Tie rule" and "Limit"; `test/s5_div1_bitwise.F90`); the HYPRE-off build fix (driver `HypreStub.cpp` and its CMake entry removed, src `fc49f96cc4`, the rewritten form of `c8e83f2b43`); Role 1's variable-density two-level case with discriminating overwrite on/off controls at 1 and 4 ranks (`617da72d85`), drafts of driver patches `0011` and `0012` with hooks and `run_mesh_dumps_check.sh` (`9442836334`, not yet run), the residual-check report of the solves (`91c75e9e49`) and the D-081 step-sequence check against single-mesh FDS (`e5164b00c6`, `9238d90e5f`); Role 2's conditions C1 to C3 (`d59803eb56`) and loop tests R2 and R3 (`8160c2cafc`); the `rad_wall_qin_zero` and `W_SURF_INDEX` fix and the wall-check hard stops, merged with `s5-four`, `s5-wall-next`, `s5-merge` and the `pvf_fvx/fvy/fvz` kernels on `s5-land` (`fd5b22dcae`, `fb990f42f6`), where `ci_checks --strict` passes (128 kernels, 111 in the golden); spec v0.4.41 with D-079 and D-080 and actions A-75 to A-77 (docs `ac04814`); the Pressure Solver Lead's sign-offs of L1211, L1207, L1220 to L1222 and L1209. A-77 closed: `WALL_LOOP_2` host path bitwise on `dec2_obst` with gather/scatter cost recorded (docs `d1e73c1`, stage1 §7.16). Also landed: GNU Build Chief items 2 and 6 (FireX tip `8ad685f435` HYPRE on/off with obst-removal and back-wall, refbin regenerated; `0011`/`0012` line check, mesh dumps, two-level `ns2d_16` np1/np4 PASS at src `c0d3573c7d`); Mesh Data COMPRESSION_WAVE batch on s5-gen `f72bd9b809` and D-083 two-pass `CHECK_MASS_DENSITY` on `s5-cmd-clip` (`8edd573d33`, Species CLEARANCE `126662fc8f`); GNU A-65 / patch 0010 host tip `047adca89bc6` (`vv-runs/a65-0010-host/RESULT.md`) and A-27 HYPRE pin plus `amrex-install-hypre232` (`vv-runs/a27-hypre-amrex/RESULT.md`); GNU UP-0011 stock Debug PASS (build+smoke; `vv-runs/up0011-stock-debug/RESULT.md`); GNU A-14 / D-019 HDF5/VTKHDF PASS on tip `047adca89bc6` Release `USE_AMREX=OFF` (scratch CMake patch `vv-runs/a14-hdf5-gnu/patches/cmake-with-hdf5.diff`, HDF5 1.14.5 at `firemodels-gnu/libs/hdf5/v1.14.5`, VTKHDF smoke green, hub `vv-runs/a14-hdf5-gnu/RESULT.md`, A-45 licence note at `firemodels-gnu/hdf5/COPYING`; show-and-tell viz remains AMReX plotfile per owner — A-14 is FR-074 background); optional 0010 Debug `-fcheck=all` PASS (`vv-runs/a65-0010-host/RESULT.md`); A-14 pointers folded into `vv/environment.md` (docs `02f7197`). Also: GNU A-14 stretch `USE_AMREX=ON`+HDF5 PASS (RESULT addendum `vv-runs/a14-hdf5-gnu/RESULT.md`; binary sha256 `e7646ad1f2af86f70cc2895e9e568a6633b905dfb4ce8d5697feb959ba2ea33c`); Species combustion/06 host-translation ACCEPT for L1272/L0866/L0875/L0397 (docs `3e98931`); Radiation solid/08 R1–R4 closed and L1251 ACCEPTed (docs `c588fbe`); Intel UP-0011 tip `8ad685f435` Rel+Dbg PASS (`vv-runs/upstream-patch-validation/intel/UP-0011-tip8ad685-rel-dbg.md`); Mesh Data D-083 `s5clip_*` `[[ported]]` restored on s5-gen `08f6016e28` with `ci_checks --strict` PASS (docs `b46f6d0`); Role 3 adaptive+plotfile TAG_* green `370c4412cf`; Role 1 D-089 refuse lift `8c8381b946`; Role 2 D-092 composite masked Built + S6 RHS fold accepted (selector flip still Role 2).

Removed by D-082 and D-083 (not landed, tabled): the revised `UP-0003` with persistent scratch and the density-loop patch for L0876 (Legacy Mapper r5 item 9), their review (Species r5 item 2), the `Z_TEMP` pad review for MP5 (Species r5 item 6), and every other MP5, `UP-0001`, `UP-0002` or `UP-0003` upstream item. The D-083 two-pass GPU kernel has landed (Mesh Data Loops Engineer; was item 1).

Numbering of patches (A-63): upstream-FDS patches are `UP-NNNN`; driver patches are the bare numbers `0001` to `0012`. The level-0 guard is `0010`; the output patches are `0011` and `0012`.

## Active under spend lock (priority order)

**SPEND LOCK in force (CEO).** Only the items in this section are active. Everyone else is idle; no new claims. Solid P2 retry is next after vortex and does not start until the CEO says.

1. **Vortex owner-show — L1 div / mass drift (Role 1 + Role 2; Role 3 diagnosis).** Case `vv-runs/owner-show/vortex_ring_re100_wall/` (TAG_VORT ring, floor INERT, sides PERIODIC, top MIRROR, `periodic=(110)` verified). Overnight run aborted at step 1157/2000: abort itself was ENOSPC creating plotfile s1158 (plot-every 2), but the run was already physically bad. Composite mass change grew ~linearly to +5.7%; L1 max|div u − D| grew from 0.125 at step 1 to 5 by step 143, peak ~40 near step 269, ~14–19 late; L0 stayed ~1e−4 to 0.13. Diagnosis in that folder's `ROLE3_PHYSICS_DIAG.md` (Role 3 authored): **primary cause B** = per-step composite projection not keeping L1 divergence-free (grows in `changed 0` stretches, so not regrid churn); **confirmed gap A** = fine levels have no wall cells (`fds_fine_level.f90` `N_EXTERNAL_WALL_CELLS=0`), only `mirror_domain_edges`, so INERT no-slip/impermeable is not enforced on L1 once L1 touches the floor (~s1000+); S6 fold not the cause; tag churn (blocking 4 / maxsize 8) rejected as primary. **Owners now: Role 1 (data layout) + Role 2 (pressure backend), woken by the CEO on L1 div/mass drift.**

2. **Vortex RESTART GATE (Architect owns readiness).** No full 2000-step restart until all gates pass, in order. Architect signs off G1/G2 before the CEO relays GO.
   - **G1** short rerun, same case, `--regrid-every 999999 --plot-every 0 --steps 50`: L1 max|div−D| must not grow step over step and must stay within 2× of L0 (and ≤ 0.2 absolute); the two relative criteria (no growth, ≤ 2× L0) apply only when L1 max|div−D| > 1e−10 (roundoff floor, Architect ruling 2026-10-08); composite mass relative change < 1e−6 over the 50 steps (closed domain: periodic x/y, INERT floor, MIRROR top, so mass should be conserved to solver tolerance).
   - **G2** 200 steps with regrid on (`--regrid-every 5`) and cubed patches `--blocking 8 --maxsize 8`, `--plot-every 20`: same div/mass criteria held across regrid steps, plus impermeability check max|W| on any L1 face at k=0 ~ 0 (floor), once L1 touches the floor. Fallback for cubes: `--blocking 16 --maxsize 16` if 8 gives too many boxes.
   - **G3** disk: at least 15 GB free before launch; plot-every 20 (pe20) only.
   - **Then full restart:** np=4, `--steps 2000`, `--plot-every 20`, `--tag-vort auto`, `--blocking 8 --maxsize 8`, `--regrid-every 5`.

3. **Solid P2 retry (Solid Phase Lead).** Next after vortex; not started until the CEO says.

## Parked until spend lock lifts

Keep existing role entries below; the following are **parked** (not active; no new claims; do not start):
- Generator: ORIG/K3D feature work and L1211 layout-contract feature.
- Legacy Mapper: WP9b 9b-10 follow-ups / FR-072 mesh-looping output routine list.
- Legacy Mapper: L0398 (waits S3 CSR).
- Mesh Data: L1305 O2 / L1317 (current claims stay claimed, work parked).
- Role 2: residuals (MLMG composite+mask; `face_gradient_composite`+mask) and optional H-dump — do not flip / do not start.
- GNU Build Chief item 15 (OPEN floor-INERT smoke after Role 1 TimeLoop die lift).
- Role 3: D-078(c) NODERIVE auto-fail assert (registry transfer check item).
- Integration Lead: movie refresh (held until a clean vortex run).

Summary of dry risk under spend lock: Role 1 and Role 2 **active** (vortex L1 div/mass); Solid Phase **queued** (P2 retry after vortex, CEO gate); all other roles **parked / idle** (no new claims). Prior ready-work counts remain in the per-role sections for when the lock lifts. FASTMATH GNU verify PASS (`vv-runs/d068-fastmath-pin/RESULT.md`). Adaptive owner-show path remains green at Role 3 `370c4412cf` + Role 1 D-089 `8c8381b946`, but the vortex_ring_re100_wall overnight is the spend-lock focus until the restart gate passes.

Shared references used below: the ranked loop list `amrex/loop-work-list.md` (and `amrex/loop_work_list.csv`; loop ids `Lnnnn`), the claim register `amrex/loop_claims.csv`, the stage-1 plan `amrex/stage1-gpu-spike-plan.md` (work packages `WPn`), the blocked-loop families `amrex/blocked-loop-families.md` (ids P1 to P4, S1 to S4, SP1 to SP4, R1, V1, V2, O1 to O4), the WP3 and WP9b note `amrex/wp3-generator-coverage.md`, the pressure loop table `amrex/pressure-velocity-h-loops.md`, the radiation candidate table `radiation/02-radiation-gpu-candidates.md`, the inventory `inventory/gpu_candidate_loops.csv`, the loop tracker `tracker/loop-tracker.md`, `tools/ported.toml`, `tools/kernel_registry.toml`, and the action log in `README.md` (ids `A-nn`). Line numbers are those of the loop work list unless an item says "after the merge".

## FireX merge: impact on loop translation

FireX `1fb5c300fa` is merged into FDS-AMReX locally (src `8ad685f435`; check report `inventory/upstream_merge_reports/2026-10-06.md`, not yet committed to the docs). It touches `func.f90`, `main.f90`, `init.f90`, `dump.f90`, `type.f90` and one line of `wall.f90`. No ported kernel body is affected.
- **Shifted line citations.** `func.f90` lines after 4532 move by +65 (new `PACK_WALL_LITE`): `PACK_CELL` L0677 is now 5248-5266, and `BLOCK_CELL` starts at 5570 with its nest L0630 now 5581-5627 (31 new lines inside). `main.f90` moves by +3 from line 1178, rising to +34 after line 3919 (`MESH_EXCHANGE`); `init.f90` by +11 at `CREATE_OR_REMOVE_OBSTRUCTIONS` and `REASSIGN_WALL_CELLS` (L0817 4889-4900, L0822 4928-4947, L0823 4951-5065); `dump.f90` by +1 to +5; `type.f90` by +3 after line 453 (`WALL_TYPE`). `divg.f90`, `mass.f90`, `velo.f90`, `wall.f90`, `pres.f90`, `radi.f90`, `turb.f90` and `soot.f90` keep their line numbers, so the generator markers (the `gsfv_*` markers in `func.f90` sit at 1348-1453, before the first hunk) do not move; comment citations into `type.f90` and `main.f90` in `s5_markers.toml`, kernel notes, the loop work list, and the driver patch anchors (`0010`, the `0011`/`0012` drafts, `upstream-patches/line-check-output-plan.md`, the `UP-0011` call blocks in `main.f90`) need a re-check.
- **New `BLOCK_CELL` branch.** A single-cell consumable-obstruction branch reads `CELL%WALL_INDEX(±1..±3)`, neighbour `TMP`, `PBAR_S`, `PRESSURE_ZONE`, `RSUM0` and `OBSTRUCTION%CONSUMABLE`, and writes `TMP`, `RHO`, `RHOS`. L0630 needs either a translation that takes these as kernel and shim arguments or a host path (it runs only when an obstruction is created or removed).
- **Exchange codes.** `MESH_EXCHANGE` code 6 no longer carries back-wall data; that moves to code 66 (full `PACK_WALL`, once at set-up) and 67 (`PACK_WALL_LITE`, every step). The driver's code 6 fills only `U`, `V`, `W`, `H` ghosts, so nothing breaks today; the driver's exchange-code mapping must keep 66 and 67 on the FDS host path and must follow the split if it ever takes over back-wall exchange.

Summary of dry risk (pre-lock counts preserved for when the lock lifts; under spend lock only Role 1 / Role 2 vortex work and queued Solid P2 are live): Role 1 low; Role 2 low; Role 3 low; Integration Lead low; Pressure Solver Lead low; Legacy Mapper low; Species and Combustion low; Solid Phase low; Radiation medium; V&V low; Spec Lead low; Generator Engineer low; Wall Loops Engineer low; Mesh Data Loops Engineer medium; GNU Build Chief medium; Intel Build Chief medium. Under spend lock: Role 1 + Role 2 **active** on vortex L1 div/mass drift; Solid P2 **queued**; all others **parked / idle** (see Active / Parked sections at top). GNU item 15 parked; Mesh Data L1305/L1317 parked; Mapper FR-072 / L0398 parked; Generator ORIG/K3D + L1211 parked; Role 3 D-078(c) assert parked; Integration movie refresh held until a clean vortex run. FASTMATH GNU verify PASS (`vv-runs/d068-fastmath-pin/RESULT.md`).

## Data layout, pressure and transport implementers

**Role 1, Data Layout. Dry risk: low (7 ready; under spend lock: vortex L1 wall/data-layout ownership ACTIVE).**
0. **ACTIVE (spend lock):** Vortex L1 div/mass drift with Role 2 — fine-level wall/data layout gap A (`fds_fine_level.f90` `N_EXTERNAL_WALL_CELLS=0`; INERT no-slip/impermeable not enforced on L1 once L1 touches the floor). Support restart gate G1/G2 readiness for the Architect. See Active section item 1. Other Role 1 items below are idle until the lock lifts.
1. Exchange-code mapping after the FireX merge: record in the driver that its code 6 fills `U`, `V`, `W`, `H` ghosts only, keep codes 66 and 67 (back-wall exchange) on the FDS host path, and add a test that a 66 or 67 call from `MAIN_LOOP` or `INITIALIZE_BACK_WALL_EXCHANGE` is never routed through the driver ghost fill (merge report, "Branch-patched Source files"; section above).
2. Driver patches `0011` and `0012` (drafts in `Source/driver/patches/`, `9442836334`): move the anchors to the merged `main.f90` and `dump.f90`, run `run_mesh_dumps_check.sh`, then ask for review; insertion points in `upstream-patches/line-check-output-plan.md`.
3. Three or more levels on the hierarchy (registry, bind per level, composite solve with Role 2, `EXACT_SUMS` on for decomposition comparisons; D-053).
4. Timing run for the Integration Lead on the two-level driver (stage-1 plan section 10).
5. Remaining WP1b dump points `CE_EDGE_INDEX`, `ED_OMEGA`, `ED_TAU` (`amrex/stage1-gpu-spike-plan.md` WP1b); record the routine that replaces `WALL_UVW_INTERP` with the Integration Lead.
6. FR-072 output plan with the Legacy Mapper's list of mesh-looping output routines (ADR-004, `amrex/fr072-output-amrex-notes.md`).
7. Check the E1 predictor `DIV1` gap against the two-level numbers (WP1b); close it or record what remains.
Blocked on: nothing for the idle list; spend-lock focus is item 0 with Role 2.
D-089: TimeLoop `FREEZE_VELOCITY` refuse lifted for cold-flow adaptive demos (src `8c8381b946`); reacting-flow production must not depend on freeze.

**Role 2, Pressure Backend. Dry risk: low (under spend lock: vortex L1 composite projection ACTIVE; residuals/H-dump PARKED).**
0. **ACTIVE (spend lock):** Vortex L1 div/mass drift with Role 1 — primary cause B from `ROLE3_PHYSICS_DIAG.md`: per-step composite projection not keeping L1 divergence-free (grows in `changed 0` stretches). Support restart gate G1/G2. See Active section item 1.
1. ~~D-092 composite masked Built (assembled HYPRE)~~ done — tip `774c804701` (Pressure Lead ACCEPT `pressure/08` §13.3; `pb_comp_masked` PASS). **In flight when lock lifts:** selector flip (Role 2 locks selector per D-091/D-092).
2. ~~S6 OPEN / inhomogeneous Poisson RHS fold~~ done — accepted; enabled the vortex overnight (S6 fold not the cause of L1 div per Role 3 diagnosis). Role 1 TimeLoop die lift for OPEN floor-INERT remains GNU item 15 (**parked**).
3. **PARKED (spend lock):** Residuals stay NotBuilt — MLMG composite+mask; `face_gradient_composite`+mask — do not flip those yet. Optional H-dump parked.
4. Tests R4 and R5 on real FDS files (`pressure/08`; `pressure/07` section 11.4); C1 to C3 and R2, R3 have landed and go to the Pressure Solver Lead for review. Idle under spend lock.
5. Upstream issue draft for the AMReX row-scaling problem, with the `habec_ijmat` symmetric-scaling candidate, in `upstream-issues/` (D-075 (4)). Draft only, not sent. Idle under spend lock.
6. Open `pres.f90` loops L1210 (`pres.f90:238-245`), L1212 to L1214 (`267-307`) and L1216 to L1218 (`404-430`): the Pressure Solver Lead's draft dispositions retire them (cylindrical deferred, transposed storage not used; draft `pressure/09-pressure-loop-dispositions.md`, not yet committed). Check that nothing in the backend depends on them and close them in `amrex/pressure-loops-r2-status.md` once the dispositions are committed. Idle under spend lock.
7. **PARKED (spend lock):** Layout contract for L1211 (`pres.f90:250-260`) with the Generator Engineer (`amrex/pressure-loops-generator-proposal.md`); L1212 to L1214 drop out if item 6 retires them. Generator ORIG/K3D + L1211 parked together.
8. Composite solve at 4 ranks, with a multi-box fine level, and with three or more levels (`pressure_backend` tests; needed by Role 1 item 3). Idle under spend lock.
9. Dry-run apply of `UP-0011` on FireX `1fb5c300fa` (its `main.f90` call blocks sit in `MAIN_LOOP`, which the merge changed) and update its target line in `upstream-patches/README.md`. Idle under spend lock.
Blocked on: nothing for the idle list; spend-lock focus is item 0 with Role 1. Residuals / H-dump parked (item 3).

**Role 3, Regrid and Multi-Level Transport. Dry risk: low (5 ready; under spend lock: IDLE except diagnosis already authored).**
0. Vortex physics diagnosis authored in `vv-runs/owner-show/vortex_ring_re100_wall/ROLE3_PHYSICS_DIAG.md` (primary B, gap A, S6 not cause, tag churn rejected). No further Role 3 spend-lock work unless the Architect asks; overnight restart is gated (see Active section).
1. Tag buffering done by the regrid code, with the zero-lost-cell assertion against a single-box reference on the 64^3 8-box case and the 3-level and ratio-4 cases (D-080 (1), A-75; `role3-regrid-transport-plan.md`, tagging section). Idle under spend lock.
2. Condition C2 in the converter: the backend already refuses `TUNNEL_PRECONDITIONER` (`d59803eb56`); add the input-converter refusal with a test so the input fails before set-up (D-080 (6), A-70). Idle under spend lock.
3. **PARKED (spend lock):** Registry transfer check against the D-078 (c) rule / NODERIVE auto-fail assert (`RSUM`, `MU`, `KRES`, `D`, `DS`, `H`, `HS` from the parent on a zero-initialised fine level), with a test that fails without the injection.
4. Two-level driver runs of the FR-016 gate cases on `--two-level-run`, then the moving blob (A-58, `vv/test-plan.md` 5.12.6; gate rows P3-B01 to B05; D-089: keep TAG_* emit into TwoLevelRun for adaptive owner-show / `blob*_amr`). Idle under spend lock.
5. t = 0 hierarchy, tagging and regrid on the real FDS stages (FR-012). Idle under spend lock.
6. Offload tagging test on the test machine (`RT_TAG_OFFLOAD`; Phase 4 item; waits for GPU time). Idle under spend lock.
Blocked on: spend lock for items 1–6; item 3 parked explicitly.
Adaptive / owner-show: TAG_* + set_tag_fields mid-run remesh green at src `370c4412cf` (with Role 1 D-089 `8c8381b946`); vortex_ring_re100_wall overnight aborted (see Active); V&V enabling adaptive as primary in `vv-runs/owner-show/`.

**AMReX Integration Lead. Dry risk: low (5 ready; under spend lock: IDLE; movie refresh PARKED).**
0. **PARKED (spend lock):** Owner-show / vortex movie refresh — held until a clean vortex run passes the restart gate.
1. ~~Encode D-085 scheme-3 device limit (`maxulp ≤ 8192` for `s5rad_plane3d_ild` under `-DS5_LIBM_CLASS`) in the device-tier judge / `kernel_categories.json` (V&V records in `vv/rad_scheme3_tolerance.md`).~~ done — gate `cb9da16` (`kernel_categories.json` + judge + selftest); follow-up rad set notes / device_logs on record.
2. Follow-on after A-77: K1 scatter for `WALL_LOOP_2`, driver wire-up and device run on the test machine (not a reopen of A-77; host path closed).
3. ~~FASTMATH CMake pin apply (D-068)~~ done — applied on s5-gen (`3546f746a7`); GNU verify PASS (`vv-runs/d068-fastmath-pin/RESULT.md`; `amrex/ported-candidates.md` §3).
4. Held-out kernels: real-field runs for `rho_d_interp`, `dp_div_heat` and a host-tier result for `d_z_max` and `del_rho_d_del_z`; `rho_d_maxloc_fix` stays out of `ported.toml` until A-76 (NS≥3 tie frequency) lands (`amrex/ported-candidates.md`, held out; Species ACCEPT-with-conditions).
5. WP9b pieces 9b-5 and 9b-7: device launch sequence part A, hook, part B (`amrex/wp3-generator-coverage.md` section 3.2).
6. Rule-7 macro form in `s4_omp.inc` and re-time (`amrex/s4-k2-single-source.md`; D-072 (d)).
7. Open actions A-28, A-31, A-44, A-48, A-49, A-51, A-54 (`README.md`).
8. ~~Sweep launch layer~~ done — re-pinned to s5-gen `d054f283f0` (ILD*/STORE tip), test-machine GPU re-run; scheme-1 costs unchanged within ~2% (RUN; stage-1 plan §7.14); 120/120 dumps sha256-equal to host. ILD*/STORE launch at Integration scale still open (Radiation §13.4 has the small-box subset).
9. WP8 flat tables on the device and WP11 wall kernels in the driver after generator change W1 (`amrex/stage1-gpu-spike-plan.md`). Waits for W1 and GPU time.
Blocked on: nothing for item 8; item 7 remaining design/write-up actions (A-48/A-49/A-51/A-54) and A-28 TinyProfiler install; A-44 full matrix and A-31 CUDA-aware check addressed on the test machine (see README).

**Pressure Solver Lead (loop duty). Dry risk: low (7 ready, 4 of them loop items).**
1. Review Role 2's conditions C1 to C3 (`d59803eb56`) and tests R2 and R3 (`8160c2cafc`); R4 and R5 when posted; rule on any written dispute with the Architect (D-080 (6); `pressure/08`).
2. A-56 backend comparison to the end: composite and masked cases with the HYPRE backend, report in `pressure/07` and `pressure/06`.
3. Non-matching Dirichlet wall data test on the hierarchy (D-075 (2): untested).
4. Loop L1365 (`velo.f90:1404-1459`, `NO_FLUX` face nests per obstruction): IDEMPOTENT marker form with `pvf_fvx` as the closest kernel, now generated (s5-gen `ce06e0c8dd`); claim it in the register (`amrex/pressure-velocity-h-loops.md`, NO_FLUX note and step 5).
5. Loops L1206 (`pres.f90:729-742`) and L1208 (`pres.f90:768-787`), residual diagnostics: the host code is written and bitwise on synthetic cases; finish the real-FDS check and commit the dispositions (draft `pressure/09-pressure-loop-dispositions.md`; `amrex/pressure-velocity-h-loops.md`, residual rows).
6. Loop L1130 (`pres.f90:5580-5601`, `GET_H_REGFACES`, P4 `LOG_INTWC`): host translation as the P4 sign-off allows, with a bitwise check (`amrex/blocked-loop-families.md` P4).
7. Loop L1121 (`pres.f90:4012-4026`, `CHECK_UNSUPPORTED_MESH`): host set-up check; state its AMR equivalent in the converter or retire it with a reason (`amrex/pressure-velocity-h-loops.md`).
8. Single-thread device add measurement for zone sums with the fixed-point mode (stage-1 plan, P1/V1 summation; `tools/README.md` zone-sum section). Waits for GPU time.
Blocked on: nothing for items 1 to 7.

## Domain leads

**Legacy Mapper (loop duty). Dry risk: low (after L1355+L1321; L1315/L1316 deferred on GE feature; A-34 deferred; L1314/L0884 reassigned; 9b-1/9b-10 done; under spend lock: IDLE; L0398 and FR-072/9b-10 PARKED).**
1. ~~Line citations after the FireX merge (section above): re-check the loop work list rows in `func.f90` (L0630, L0677), `init.f90` (L0817, L0822, L0823) and `main.f90`; the comment citations into `type.f90` and `main.f90` in `s5_markers.toml` when s5-gen re-pins; kernel notes; the anchors of driver patch `0010`, the `0011`/`0012` drafts and `upstream-patches/line-check-output-plan.md`. Write the result into the merge report's follow-up and the line-check file.~~ done — tip re-check on `315ced11efc5` (`upstream-patches/line-check-firex-1fb5c300.md` § Mapper tip re-check; `amrex/upstream-1fb5c300-loop-impact.md`); Intel UP-0011 Rel+Dbg PASS; GNU Debug CLOSED; L1356 tip 265-288 noted; no inventory pin rewrite.
2. ~~L0401 (`divg.f90:995-1088`, `SPECIES_ADVECTION_PART_1_NEW`) under D-082: run the bitwise tests of `rho_z_p_divg` and `spec_wall_zz` for the non-MP5 limiters, set the loop translated if they pass, and mark the MP5 wall element out of scope in the work list (`amrex/blocked-loop-families.md` S2).~~ done — non-MP5 PASS (`rho_z_p_divg` R2 O0; `spec_wall_zz` walltab/scratch quick); claim `translated`; S2 MP5 wall element out of scope (D-082); s5-gen tip `a0fa71a665`.
3. ~~L1391 (`velo.f90:1243-1253`, `VELOCITY_FLUX_CYLINDRICAL`): marker and bitwise test; translatable now (`amrex/wp3-generator-coverage.md` WP3 item 1).~~ done — marker `vflux_cyl_vort_tau`; R2 O0 dpd+bind PASS 1458/0/0 non-vacuous (OMY/TXZ scramble after `vflux_vort_tau`); s5-gen `a0fa71a665`; claim Legacy Mapper `translated` (AMR driver still aborts CYLINDRICAL per D-042; nest isolatable).
4. ~~Cell-loop batch L1355 / L1315 / L1316 / L1321; A-34.~~ partial — L1355 `visc_vreman` + L1321 `varden_uvw_s` translated on s5-gen `57965fcb99` (`test_mapper4_cell.py` 39/39). L1315/L1316 TEST_FILTER deferred (generator feature: variable stencil `ORIG(I+L,...)` + local PARAMETER `K3DT`/`K3DS`). **A-34 deferred** (still open/narrowed: FX/FY/FZ tangential J=0/K=0 / STORE_SPECIES_FLUX / CC_IBM; no new footprint evidence this turn). L1356 Wall `visc_wale` unchanged. Next: L0398 or remaining TEST_FILTER when GE lands the feature.
5. **PARKED (spend lock):** L0398 (`divg.f90:835-939`, `ENTHALPY_ADVECTION_NEW` wall loop): accumulate through the wall gas-cell subscript with the S3 cell table (`amrex/blocked-loop-families.md` S2, S3; `amrex/l0398-pre-analysis.md`). Waits S3 CSR; parked until lock lifts.
6. ~~L1314 (`turb.f90:1458-1479`) and L0884 (`part.f90:4606-4617`)~~ **reassigned (Architect/CEO 2026-10-06):** L1314 TENSOR_DIFFUSIVITY / private `tdm_*` stays Generator-owned — Mapper must not mark/claim; L0884 `CLIP_PARTICLE_DRAG` left for Mesh Data. Do not pick these up.
7. ~~Register and work list: move the L0876 and L0877 claims to the D-083 owner, L1356 to the Wall Loops Engineer, register the radiation sweep as L1242 with the FR-062 tag (D-073 (d)), paste the Species and Radiation R1 sign-offs, then regenerate the loop work list and the pressure-loop table (`tools/inventory/loop_work_list.py`; L1272, L0866, L0875, L0397 still read "claimed").~~ done — claims refresh + regen 2026-10-06 (L0876/L0877 `translated` hand-written D-083 on Mesh Data / s5-gen `08f6016e28`; L1356 Wall confirmed; L1242 FR-062 note refreshed; Species/Radiation R1 pastes; L1272/L0866/L0875/L0397 no longer read claimed).
8. ~~WP9b items 9b-1 (chain map of `DIVERGENCE_PART_1` against `run_divergence_part1`) and 9b-10 (patch 0009 line shift)~~ done — `amrex/wp9b-1-chain-map.md`; 9b-10 line-shift rule (+9 / +15) in that note §7 and `flux-hooks-kernel-review.md` line 6 (`amrex/wp3-generator-coverage.md` §3.2).
9. Patch table: add rows for the driver drafts `0011` and `0012` and keep `upstream-patches/README.md` in step (A-63).
10. **PARKED (spend lock):** List of every output routine that loops over meshes (FR-072, ADR-004, `amrex/fr072-output-amrex-notes.md`). WP9b 9b-10 follow-ups also parked.
11. Keep the input classification list current against the converter output (D-076; `amrex/level-jump-external-wall-check.md`).

**Species and Combustion Lead (loop duty). Dry risk: low (6 ready, 3 of them loop items; current: item 4 L0596).** CLIP_MASS/CLIP_CELLS landed `9b93ac753c` (C11/C12 with V&V next).
1. Sign the D-051 device sign-off for `rho_d_maxloc_fix` on the conditions of D-080 (2): the tie test and the NaN note have landed (generator README "Tie rule" and "Limit"); the real-data tie frequency is still to be reported (A-76; WP9b 9b-2).
2. Reviewer of the D-083 two-pass `CHECK_MASS_DENSITY` kernel (L0876, L0877): check it against the unchanged CPU loop under the D-031 rules, as the D-051 device sign-off (D-083; `combustion/02` checks 1 to 4 still describe the order and rounding to verify).
3. ~~Review the host translations of L1272 (`SETTLING_VELOCITY`), L0866, L0875 and L0397~~ done — ACCEPT per CEO (`combustion/06`; docs `3e98931`; tip check s5-gen `08f6016e28`).
4. **Current:** Loop L0596 (`fire.f90:1902-1918`, `COMBUSTION_BC`): same-mesh wall part with a per-wall table (D-055 excludes `NOM>0`), translate or host path with a bitwise test (`amrex/loop_work_list.csv`; S3).
5. Loop L1273 (`soot.f90:412-472`, `SOOT_SURFACE_OXIDATION`): per-wall data as for L0596 (`amrex/loop_work_list.csv`).
6. Loops L0866 (`mass.f90:570-577`) and L0875 (`mass.f90:752-759`) device form: `PBAR_S`/`PBAR` by `PRESSURE_ZONE` through an integer zone table (`[policy.tables]`, `type = "integer"`, generator README), bitwise against the host translation.
7. Species tagging cases for R4 and the clip-line budget test inputs (`combustion/03-combustion-under-refinement-and-tagging.md` section 4; with Role 3 and the V&V Lead). **Unblocked** by Role 3 adaptive TAG_* green at `370c4412cf`.

**Solid Phase Lead (loop duty). Dry risk: low (7 ready, 4 of them loop items; under spend lock: P2 retry QUEUED after vortex).**
0. **QUEUED (spend lock):** Solid P2 retry — next after vortex; not started until the CEO says. All other Solid items idle until then.
1. Forced-flip survey on the host to measure `E_class`, signed with the V&V Lead before the first device run (D-077 (1) V3; `solid/09-flip-budget-gate.md`). Idle under spend lock.
2. Flip-injection build and exposure tagger for the gate (`solid/09`; V4).
3. A1 predicate form in the cap design and the status-301 template in the driver abort path (D-077 A1, A2; `solid/06-solve-port-test-design.md`).
4. Loop L1488 (`wall.f90:144-185`, `WALL_BC`): SP4 hand port with generator help and the ordered-sum reference for the shared-state sums (`amrex/wall-bc-translation-plan.md`; SP4).
5. Loops L1487 (`wall.f90:124-137`) and L1486 (`wall.f90:109-119`): SP4 conditions; L1486 is generated by the Wall Loops Engineer with the call-site inlining of `NEAR_SURFACE_GAS_VARIABLES` (s5-wall-next `b4783bf0e1`); you hold the claims and sign.
6. Loops L1470 and L1471 (`wall.f90:3635-3717`, HT3D exchange): same-mesh part, `NOM>0` excluded (D-055); L1471 is a separate pass after the 1-D solve (`amrex/loop-work-list.md` section 3).
7. Loop L0817 (`init.f90:4889-4900` after the merge, `CREATE_OR_REMOVE_OBSTRUCTIONS`): flat `ZZ0` table as for `MW_SPEC` in `s5_markers.toml`; note the routine is now also called on restart (FireX merge).
8. Loops L0822 and L0823 (`init.f90:4928-5065` after the merge, `REASSIGN_WALL_CELLS`): host path for the event-driven re-assignment with the table refresh (`amrex/loop_work_list.csv`, host-side).
9. Follow-ups W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7; O3 (L1402) uniqueness review with the Wall Loops Engineer.
10. ~~L1485 retirement~~ done — retired in the AMR route (D-065 Q1 / D-077 (2) / A-65 green GNU+Intel); no device kernel; host abort guard patch 0010; evidence `vv-runs/a65-0010-host/RESULT.md` (tip `047adca89bc6`; P-3 / N-5 / N-5b / N-5c); claims `not planned (retired)`; FDS-only keeps wall.f90:888-947.

**Radiation Lead (loop duty). Dry risk: medium (setup L0825-27/L0676/L0664-65/L0006-07 confirmed host or GEOM-deferred in radiation/06 §18; L1240/L1244 host-translated on s5-gen tip `f9c9bed577` (L1240 host at `313c54ef43`), Wall BR-join confirmed; L1244 DLO mask hoist on `f9c9bed577`; duplicate `s5rad_lp_qin_zero` dropped in `d283813c68`; ready for ported.toml when the gate opens).**
1. Loop L1242 sweep kernels: landed on s5-land (`664917ad50`); no new s5-gen commits until Mapper posts FF. Remaining A1 pieces classified/feature-requested only (`radiation/06` §19; `radiation/02` A1).
2. Loops L1248, L1249, L1250, L1247, L1238, L1241: translated on s5-land (see radiation/06 §16-§17). L1240 / L1244 (both host): translated on s5-gen tip `f9c9bed577` (L1240 host `313c54ef43`; radiation/06 §16.3); L1244 DLO mask hoist on `f9c9bed577` (mask once per particle, reused across bands); duplicate bind(c) dropped in `d283813c68`; not in `ported.toml` until the gate opens.
3. Gate hookup for radiation kernels (`radiation/03-radiation-translation-notes.md`). Open; waiting Integration `rad_*` device driver (R1–R4 and L1251 closed).
4. ~~Follow-ups R1 to R4 from `solid/08-new-wall-kernel-review.md` section 7~~ done — closed on Radiation (`radiation/05` + `radiation/06` §14; docs `c588fbe`). L1251 ACCEPTed (Role 3 entry); waiting Integration `rad_*` device driver / gate hookup.
5. Ghost `UII`/`UIID`/`QR` readers: V&V closed D-073(b) with keep-ghost-bits for CASE(18)/(19) (`vv/d073b_ghost_readers.md`); A-50 box-split ghost clause closed. Remaining radiation verification variants as needed (FR-060/FR-062/FR-063).
6. Radiation under regrid (`radiation/04-radiation-phase4-notes.md`; FR-063). Waits on Role 3.
7. Setup/deferred closes confirmed: L0825-L0827 host (FR-062 exchange), L0676 host, L0664/L0665 host, L0006/L0007 GEOM deferred D-033 (`radiation/06` §18).

**V&V Lead. Dry risk: low.**
1. D-085: record ratification in `vv/rad_scheme3_tolerance.md` (Architect number `maxulp ≤ 8192` for scheme-3 `s5rad_plane3d_ild` under LIBM-CLASS); ask Integration to encode the limit in the device-tier judge / `kernel_categories.json`. Schemes 1–2 stay bitwise; host gate stays bitwise.
2. `run_class_tolerance` in `kernel_categories.json` and the `RUNCMP` lines of the harness (D-080 (5); `vv/dt_class_tolerance.md`).
3. Device-tier 2-ulp harness (D-070, D-072 (b)), including the run-level class-tolerance comparison for `wall_visc_les`.
4. Tests for D-080 and D-083: zero-lost-cell tag-buffer rows (P3 tagging), review of the landed `rho_d_maxloc_fix` tie test, the `WALL_LOOP_2` gate on `dec2_obst`, pressure tests R4 and R5, and the acceptance cases for the D-083 two-pass kernel (clipping exercised, bitwise against the CPU loop).
5. Drop MP5 from the case inventory and gate sets as out of scope (D-082); MP5 cases are not counted (`vv/case_inventory.md`, `vv/test-plan.md`).
6. P3-B09 distance d(n) = 4n in the test plan and `vv-runs/phase3/phase3_cases.csv` (A-61).
7. Flip-budget harness (D-077 (1); `solid/09`, test plan 5.10).
8. Gate rows of the Phase 3 moving blob; decomposition rows with `EXACT_SUMS` on (`vv/test-plan.md` 5.11, 5.12.6).
9. Gate hookup for the radiation kernels (`radiation/03-radiation-translation-notes.md`).
10. Debug `int_1to2` and the 13 skipped four-rank Intel cases when the test machine is idle (`vv/baseline_status.md`).
11. Held-out kernel reviews as evidence lands (`vv/ported-review.md`); A-46, A-47, A-26, A-39, A-42, A-49.
12. Owner-show adaptive primary: enable adaptive as primary under `vv-runs/owner-show/` (Role 3 `370c4412cf` + Role 1 D-089 `8c8381b946`; GNU rebuild smoke is GNU item 11).

## GPU and build roles

**GPU Generator Engineer (loop duty). Dry risk: low (8 ready, 4 of them loop items; under spend lock: IDLE; ORIG/K3D + L1211 PARKED).**
1. `DIVERGENCE_PART_1` small nests that the front end accepts but have no committed bitwise test: L0371 (`divg.f90:463-471`), L0373 (`499-505`), L0377 (`573-582`), L0378 (`591-597`), L0380 (`622-629`), L0382 (`668-674`) (`amrex/loop_work_list.csv`; `amrex/s4d-divergence-part1-map.md`).
2. Remaining nests of L0365 (`divg.f90:128-235`; interpolation nest 133-140 and the call at 187), L0369 (`287-421`; nests 390-402 and on, `WALL_LOOP_2` stays on the host) and L0381 (`640-660`; the outer loop) (`amrex/wp3-generator-coverage.md`).
3. ~~L1375 (`velo.f90:649-653`, `EVALUATE_RAMP` callee)~~ done — generated K2 `vflux_gravity` on s5-gen `eedf0c9a5f` (ramp inlined; `test_vflux_gravity.py` 46/46).
4. Feature-gated nests L0379 (`divg.f90:608-616`, rank-2 `R_PBAR`), L0383 (`681-687`, constant `GVEC` subscript), L0386 (`757-767`, `IPZ`), L0364 (`93-95`, callee flatten) (`amrex/loop_work_list.csv`).
5. Merge `s5-land` (`fb990f42f6`: `s5-four`, `s5-wall-next`, `s5-merge`, strict `ci_checks` pass) into `s5-gen` and regenerate once; bring in `s5-wall-sp`; then G2 sidecar justification and the G1 device-list structural test (`solid/08` G1 to G3).
6. D-082 in the generator: mark the `gsfv_mp5` marker, its golden entry and the MP5 cases and mutants of the gsfv tests out of scope (kept, not counted); drop MP5 from plans.
7. K2 conventions and harness support for the D-083 two-pass kernel (owner: Mesh Data Loops Engineer).
8. Generator change W1: `WLIST`/`NWL` indirection in the wall kernels (unblocks WP11).
9. **PARKED (spend lock):** Layout-contract feature for the pressure loops / L1211 (`amrex/pressure-loops-generator-proposal.md`). ORIG/K3D (variable stencil `ORIG(I+L,...)` + local PARAMETER `K3DT`/`K3DS` for Mapper L1315/L1316) also parked.
10. L1392, L1393 (K,I nest with J fixed, deferred O1); trace feature G1, status-301 predicate G2, ULP line format G3 (`solid/09`); fixed-point zone-sum how-to (`amrex/generator-howto.md`).

**GPU Wall Loops Engineer (loop duty). Dry risk: low (6 ready, 4 of them loop items).**
1. ~~L1356 (`velo.f90:256-279`, `COMPUTE_VISCOSITY`)~~ done — generated K2 `visc_wale` on s5-gen `8dcb39998c` (on tip `6dca8b2209`); D-086/D-088; claim translated (`amrex/loop_claims.csv`).
2. ~~L1354 / L1351 `COMPUTE_VISCOSITY`~~ done by Pressure Lead — generated K2 `mu_deardorff` / `mu_consmag` on s5-gen `6dca8b2209` (with L1360 `strain_rate`, L1379 `coriolis_up`); tests 36/36.
3. L1486 (`wall.f90:109-119`, `WALL_BC`): generate with the call-site inlining (`b4783bf0e1`) under the SP4 conditions; the Solid Phase Lead holds the claim and signs (`solid/04-blocked-loop-signoff.md`).
4. L1452 (`wall.f90:319-339`, `ASSIGN_GHOST_VALUE`): SP2 per-target gather (highest wall index wins), with the Solid Phase Lead (SP2; `solid/04`).
5. Run-level check that `dt_coupled` kernels are compared with `run_class_tolerance` (D-080 (5); `vv/dt_class_tolerance.md`).
6. Remaining W1 to W5 items from `solid/08` section 7 (the hard stops have landed in `fd5b22dcae`); SP4 scratch sums and the S2 guarded assertion.
7. O3 (L1402) uniqueness review with the Solid Phase Lead; WP8/WP11 wall-table gather with the Integration Lead (waits for W1).

**GPU Mesh Data Loops Engineer (loop duty). Dry risk: medium (L1305/L1317 claimed but PARKED under spend lock; device run B / `ci_checks` still waiting on GPU).**
1. ~~D-083 two-pass GPU kernel for `CHECK_MASS_DENSITY` (L0876/L0877)~~ done — fully closed: hand-written K2 two-pass on `s5-cmd-clip` (kernel `8edd573d33`; caller docs `f5cbeec5e8`; device evidence `b582e79027` / DEVICE.md); Species ACCEPT-with-conditions then CLEARANCE YES (docs `126662fc8f`); Integration restored `[[ported]]` on s5-gen `08f6016e28` (docs `b46f6d0`); `ci_checks --strict` PASS. Remaining wait: device run B / Wall merge when GPU free (nclip 0/1 and multi-box FillBoundary conditions still in force).
2. ~~L0630 `BLOCK_CELL`~~ closed — not planned (host); cite `amrex/l0630-block-cell-host.md` (no MESH_EXCHANGE 66/67 translation).
3. L1366 (`velo.f90:1463-1559`, `NO_FLUX` walls): translated on s5-gen `20fb4a1f03` (`amrex/loop_claims.csv`); keep thin-pair / gate evidence current.
4. ~~Neighbour-mesh L1363/L1364/L1367/L1399 (D-055)~~ closed — not planned (host); same-mesh empty; `amrex/o2-edge-exchange-scope.md`.
5. ~~`COMPRESSION_WAVE` batch L1276 to L1287 and L1379 `CORIOLIS_FORCE`:~~ done — CW markers `cw3_*`/`cw4p_*`/`cw4c_*` + r2 bitwise on s5-gen `f72bd9b809`; L1379 `coriolis_up` by Pressure Lead on s5-gen `6dca8b2209`.
6. **Current claim (PARKED under spend lock):** L1305 `SYNTHETIC_TURBULENCE` (`turb.f90` ~1802-1867, vent eddy contrib) O2 and L1317 `TEST_FILTER` (`turb.f90` ~1033-1049, midpoint 3×3×3); coordinate with Mapper L1315/L1316/L1321 when lock lifts. Claims stay claimed; work does not proceed. **Do not claim** `TENSOR_DIFFUSIVITY` L1309–L1314 (Generator private `tdm_*` under L0365). Optional warm-up (unclaimed until started): L0884 `CLIP_PARTICLE_DRAG`.
7. Keep `ported.toml` and the strict `ci_checks` pass current; `gsfv_mp5` stays listed but out of scope (D-082); `rho_d_maxloc_fix` only after the sign-off. Waits for GPU / Wall merge as needed.
8. Wire the fixed-point zone-sum mode into the `zone_sum_order` checks.
9. Device run B of the four registered wall nests (stage-1 plan 7.13b) and held-out kernel evidence with the Integration Lead. Waits for GPU time.

**GNU Build Chief. Dry risk: medium (~4 ready: items 11–14; under spend lock: IDLE; item 15 PARKED).**
1. ~~Validate patch 0010 on GNU and run P-3 and N-5 (N-5b on the bypass build) on the host build (A-65)~~ done — tip `047adca89bc6`; P-3 stock PASS (NIC_CHECK 256×3 / 512); N-5 and N-5b bypass PASS (ERROR(9001) hand tables); N-5c stock PASS (TRNX_ID converter refusal); RESULT `vv-runs/a65-0010-host/RESULT.md` (converter-supersession for stock N-5 noted there).
2. ~~Build and validate the merged tip `8ad685f435` (FireX `1fb5c300fa`) with HYPRE on and off~~ done — tip `8ad685f435` HYPRE on/off, obst-removal + back-wall PASS, Release/Debug refbin regenerated; UP-0011 unused-FN diagnostic filed (`upstream-patches/UP-0011-gnu-debug-validation.md`).
3. ~~Debug check of the `UP-0011` Debug build fix (`upstream-patches/README.md`); `UP-0008` Rel+Dbg already PASS.~~ done — stock Debug build+smoke PASS (`vv-runs/up0011-stock-debug/RESULT.md`; note updated); `UP-0008` Rel+Dbg already PASS.
4. ~~A-27 (GNU HYPRE pin and the second AMReX install)~~ done — HYPRE pin at `firemodels-gnu/libs/hypre/63331f19`; second AMReX install `amrex-install-hypre232`; RESULT `vv-runs/a27-hypre-amrex/RESULT.md`.
5. ~~FASTMATH pin on GNU (D-068)~~ done — GNU verify PASS on s5-gen tip `dcada34a2e` (pin `3546f746a7` ancestor): BF-05 NOTE + BF-07 pass; configure ON→FATAL_ERROR; RESULT `vv-runs/d068-fastmath-pin/RESULT.md`.
6. ~~Two-level `ns2d_16` at 1 and 4 ranks; line check of `0011`/`0012`~~ done — `0011`/`0012` line check, mesh dumps, two-level `ns2d_16` np1 and np4 all PASS (src tip `c0d3573c7d`).
7. ~~**A-14 / D-019 — HDF5-enabled GNU build**~~ done — PASS tip `047adca89bc6`, Release, `USE_AMREX=OFF`; scratch CMake patch `vv-runs/a14-hdf5-gnu/patches/cmake-with-hdf5.diff` (`USE_HDF5` + `HDF5_HOME`); HDF5 1.14.5 at `firemodels-gnu/libs/hdf5/v1.14.5`; VTKHDF smoke green; hub `vv-runs/a14-hdf5-gnu/RESULT.md`; A-45 licence note at `firemodels-gnu/hdf5/COPYING`. Show-and-tell viz remains AMReX plotfile per owner; A-14 is FR-074 background.
8. ~~Optional Debug `-fcheck=all` build of the tip with patch 0010 (host)~~ done — PASS (`vv-runs/a65-0010-host/RESULT.md`, Debug `-fcheck=all` row / §).
9. ~~Fold A-14 pointers into `vv/environment.md`~~ done — assigned fold landed (docs `02f7197`; §2 HDF5/VTKHDF note and §6 point 2).
10. ~~Stretch: `USE_AMREX=ON` + same HDF5 CMake patch / smoke~~ done — A-14 `USE_AMREX=ON` + HDF5 PASS (RESULT addendum `vv-runs/a14-hdf5-gnu/RESULT.md`; binary sha256 `e7646ad1f2af86f70cc2895e9e568a6633b905dfb4ce8d5697feb959ba2ea33c`).
11. **Next ready:** Rebuild Release `USE_AMREX=ON` at tip including Role 3 `370c4412cf` + Role 1 D-089 `8c8381b946` for adaptive owner-show; smoke `blob_16_l0 --two-level-run --steps 20 --regrid-every 5 --plotfile`; RESULT under `vv-runs/owner-show/`.
12. Overnight binary rebuild Release `USE_AMREX=ON` at tip ≥ `54f3441fb6` (R=0.30 a=0.20 MIRROR) for `vortex_ring_re100_wall`; record sha256.
13. Rebuild including Role 2 D-092 `ba730a0dfd` + S6 `bf0ba2d523`; link+smoke.
14. two-level `ns2d_16` np1/(np4) regression on that tip.
15. **PARKED (spend lock):** OPEN floor-INERT smoke after Role 1 TimeLoop die lift.

**Intel Build Chief. Dry risk: medium (3 ready; N-5b in progress).**
1. `-O2 -fp-model=precise` speed cost (A-73) and the cause of the `-prec-div` failures; build notes and scripts to the D-079 flags (`vv/environment.md`).
2. ~~UP-0011 tip Rel+Dbg on `8ad685f435`~~ done — PASS (`vv-runs/upstream-patch-validation/intel/UP-0011-tip8ad685-rel-dbg.md`). Still open: full Intel baseline subset on the merged tip under oneAPI.
3. **In progress / assigned:** N-5b control on Intel with the bypass build, keyed on messages (D-078 (d)).
4. A-36 optional ifx compile check for K2 on the generated kernels on `s5-land`.
5. The 13 skipped four-rank cases and Debug `int_1to2` with V&V; FASTMATH and radiation gate family; two-level `ns2d_16` under oneAPI. Wait on others.

## Program

**Spec and Program Lead. Dry risk: low (5 ready).**
1. Carry D-081 to D-083 into the spec: MP5 out of scope and not counted for coverage, `UP-0001` to `UP-0003` tabled in the roadmap and `owner-upstream-patch-table.md`, the D-083 GPU-only kernel, the single-mesh ground truth.
2. Record the FireX merge (`8ad685f435`) as the new base in the spec and roadmap; link the 2026-10-06 merge report once it is committed.
3. Roadmap dates from Role 1's two-level work and the remaining Role 1 items.
4. Regenerate `tracker/loop-tracker.md` after the Legacy Mapper regenerates the work list; the "ported" count treats `gsfv_mp5` as out of scope.
5. Fix the remaining `0010`/`0011` output-patch wording to `0011`/`0012` (A-63).
