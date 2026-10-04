
# Next work per role

Owner rule: nobody sits idle. Each role has an ordered list of ready items, each with a source reference (document path and item id). The Chief Architect owns these lists; the Spec & Program Lead refreshed an earlier version at spec v0.4.35, and this is the Chief Architect's refresh (r5) through D-080. Where another list differs, this one wins. A role that finishes its list tells the Chief Architect. "Dry risk" is the chance the role runs out of ready work within two days (low: five or more ready items; medium: three or four; high: fewer than three or all blocked). Items that wait for another role are counted as not ready when the level is set.

Status as of spec v0.4.40 plus D-080, next-work revision r5.

Decisions since r4: D-080 (complete tag buffering by the regrid code, `rho_d_maxloc_fix` sign-off conditions, `WALL_LOOP_2` host path, test-lock enforcement, `run_class_tolerance` for dt-coupled kernels, pressure sign-off conditions C1 to C3). Earlier: D-079 (Intel flags `-O2 -fp-model=precise`, converter TRN refusal), D-078 (libm entry for `wall_visc_les`, six `gsfv_*` kernels co-signed, patch 0010 on Intel), D-077 (flip-budget gate, NIC>1 guard ordering, reviewer of record, shared test lock), D-076 (level-jump wall check, abort guard, input converter), D-075 and D-074 (pressure results, 23 FDS-only inputs), D-073 and D-072 (radiation sweep design, libm list, "ported" definition).

Landed since r4 (items removed from the lists): S14.4; the L1209 OPEN-wall comparison closed (bitwise equal to FDS); the masked-MLMG comparison A-68 (131 of 131); the N-5b negative control and its test-plan text; the host translations of L1272, L0866, L0875 and L0397; the wall device driver; the converter refusal of `TRNX`/`TRNY`/`TRNZ` (D-079 (c)); `MASS_TOL` and `E1_ZZ_TOL` retightened to 1e-13; Role 1's fine-mesh velocity boundary-condition fix with its trigger test and the overwrite-off, multi-box and 4-rank controls on `ns2d_16`; the converter wiring in `main.cpp`; the `wall_visc_les` registry entry; the line check of the output plan (`upstream-patches/line-check-output-plan.md`); the `UP-0011` Debug build fix; the radiation sweep re-pin and 4M-cell level timing note (stage-1 plan 7.14).

Numbering of patches (A-63): upstream-FDS patches are `UP-NNNN`; driver patches are the bare numbers `0001` to `0012`. The level-0 guard is `0010`; the output patches are `0011` and `0012`.

Shared references used below: the ranked loop list `amrex/loop-work-list.md` (and `amrex/loop_work_list.csv`; loop ids `Lnnnn`), the stage-1 plan `amrex/stage1-gpu-spike-plan.md` (work packages `WPn`), the blocked-loop families `amrex/blocked-loop-families.md` (ids P1 to P4, S1 to S4, SP1 to SP4, R1, V1, V2, O1 to O4), the WP3 and WP9b note `amrex/wp3-generator-coverage.md`, the loop tracker `tracker/loop-tracker.md`, the ported candidates `amrex/ported-candidates.md`, and the action log in `README.md` (ids `A-nn`).

Summary of dry risk: Role 1 low (overloaded, not dry); Role 2 low; Role 3 low; Integration Lead low; Pressure Lead low; Legacy Mapper low; Species and Combustion low; Solid Phase low; Radiation medium; V&V low; Spec Lead low; Generator Engineer low; Wall Loops Engineer medium; Mesh Data Loops Engineer low; GNU Build Chief medium; Intel Build Chief medium. No role is high.

## Data layout, pressure and transport implementers

**Role 1, Data Layout. Dry risk: low (the list is long and the order matters).**
1. Non-trivial-density two-level case (variable density, with the interface overwrite on and off): the `ns2d_16` controls cannot discriminate the overwrite because density is constant (`vv/test-plan.md` 5.12.6; D-050).
2. Three or more levels on the hierarchy (registry, bind per level, composite solve with Role 2, `EXACT_SUMS` on for decomposition comparisons; D-053).
3. Timing run for the Integration Lead on the two-level driver (stage-1 plan section 10).
4. Output patches `0011` and `0012` (driver series): the line check is done; steps A and B of the output plan are driver patches `0011` and `0012`, with the insertion points in `upstream-patches/line-check-output-plan.md`. Use that numbering when renaming from the output plan's own labels.
5. Remaining WP1b dump points `CE_EDGE_INDEX`, `ED_OMEGA`, `ED_TAU` (not blocking; `amrex/stage1-gpu-spike-plan.md` WP1b); record the routine that replaces `WALL_UVW_INTERP` with the Integration Lead.
6. FR-072 output plan with the Legacy Mapper's list of mesh-looping output routines (ADR-004, `amrex/fr072-output-amrex-notes.md`).
7. Check the E1 predictor `DIV1` gap against the two-level numbers (`amrex/stage1-gpu-spike-plan.md` WP1b); close it or record what remains.
Blocked on: nothing.

**Role 2, Pressure Backend. Dry risk: low.**
1. Condition C1 of the Pressure Lead's sign-off: whole-extent assertion for the periodic wrap in `pres_h_bc_*` (D-080 (6); `pressure/08`).
2. Condition C3: the wind guard covers 0..KBP1 or validates KK per wall, with the KK=0 and KK=KBP1 synthetic throw case (D-080 (6); `pressure/08`).
3. Tests R4 and R5 on real FDS files (`pressure/08`; `pressure/07` section 11.4). Dispute any condition in writing with reasons if you disagree.
4. Upstream issue draft for the AMReX row-scaling problem, with the `habec_ijmat` symmetric-scaling candidate, in `upstream-issues/` (D-075 (4)). Draft only, not sent.
5. Open `pres.f90` loops: L1210 (cylindrical PRHS), L1212 to L1214 (transposed PRHS), L1216 to L1218 (transposed copies; retire like L1215) (`amrex/pressure-loops-r2-status.md`, "Left").
6. Layout contract for L1211 and L1212 to L1214 with the Generator Engineer: integer extents of an `exact` allocation as arguments to cell-loop kernels (`amrex/pressure-loops-generator-proposal.md`).
7. Composite solve at 4 ranks and with a multi-box fine level, and the composite solve with three or more levels (`pressure_backend` tests; needed by Role 1 items 1 and 2).
Blocked on: nothing.

**Role 3, Regrid and Multi-Level Transport. Dry risk: low.**
1. Tag buffering done by the regrid code (grow tags with a ghost exchange of the tag array to the buffer depth), with the zero-lost-cell assertion against a single-box reference on the 64^3 8-box case and the 3-level and ratio-4 cases (D-080 (1); `role3-regrid-transport-plan.md`, tagging section; V&V P3 tagging cases).
2. Condition C2: refuse `TUNNEL_PRECONDITIONER` in the input converter, with a test (D-080 (6); `pressure/08`).
3. Registry transfer check against the D-078 (c) rule: a zero-initialised fine level takes `RSUM`, `MU`, `KRES`, `D`, `DS`, `H` and `HS` from the parent, with a test that fails without the injection (Role 1's derive-injection rule is in the `LevelRegistry.H` header).
4. Two-level driver runs of the FR-016 gate cases on `--two-level-run`, then the moving blob (A-58 thresholds, `vv/test-plan.md` 5.12.6; gate rows P3-B01 to B05 in `test_blob_registry --gates`), with the multi-box fine level from Role 1 for the box-split rows.
5. t = 0 hierarchy, tagging and regrid on the real FDS stages (FR-012).
6. Offload tagging test on the test machine (`RT_TAG_OFFLOAD`, CUDA AMReX install; Phase 4 item).
Blocked on: nothing for items 1 to 3.

**AMReX Integration Lead. Dry risk: low.**
1. `WALL_LOOP_2` host path (`divg.f90:333-401`): one gather of `ZZP` and `TMP` at wall-adjacent cells per stage, loop body on the host, one scatter back; estimate 1.5 to 2.5 work-days; gates are bitwise equality to FDS on `dec2_obst`, gather/scatter cost measured per stage, and the D-064 host-side wall pass rule (D-080 (3); `amrex/wp3-generator-coverage.md` section 3.2).
2. FASTMATH CMake pin: apply `tools/patches/cmake-fastmath-pin.patch` to the generator worktree (D-068; `amrex/ported-candidates.md` section 3).
3. Held-out kernels: real-field runs for `rho_d_interp`, `dp_div_heat` and a host-tier result for `d_z_max` and `del_rho_d_del_z`; `rho_d_maxloc_fix` waits for the sign-off of D-080 (2) (`amrex/ported-candidates.md`, held out).
4. Sweep launch layer: one launch per wavefront plane across all boxes of a level (D-073 (f), D-075 (7); stage-1 plan 7.14).
5. WP8 flat tables on the device and WP11 wall kernels in the driver (after generator change W1, `WLIST`/`NWL`), `amrex/stage1-gpu-spike-plan.md`.
6. WP9b pieces 9b-5 and 9b-7: device launch sequence part A, hook, part B (`amrex/wp3-generator-coverage.md` section 3.2).
7. Rule-7 macro form in `s4_omp.inc` and re-time (`amrex/s4-k2-single-source.md`; D-072 (d)).
8. Open actions A-28, A-31, A-44, A-48, A-49, A-51, A-54 (`README.md`).
Blocked on: GPU time on the test machine for items 4 and 5.

**Pressure Solver Lead. Dry risk: low.**
1. Review of Role 2's response to the sign-off conditions C1 to C3 and tests R4 and R5; rule on any written dispute with the Architect (D-080 (6); `pressure/08`).
2. A-56 backend comparison to the end: composite and masked cases with the HYPRE backend, report in `pressure/07` and `pressure/06`; delete regenerable scratch afterwards.
3. Non-matching Dirichlet wall data test on the hierarchy (D-075 (2): untested).
4. Sign-off review of Role 2's host translations of L1211, L1207, L1220 to L1222 and L1209 against `pres.f90` (`amrex/pressure-loops-r2-status.md`; `amrex/pressure-velocity-h-loops.md`).
5. Single-thread device add measurement for zone sums, now with the fixed-point mode of the generator (stage-1 plan section on P1/V1 summation; `tools/README.md` zone-sum section).
6. Loop candidates for the pressure domain: L1365 and the other `pres.f90` loops of `amrex/pressure-velocity-h-loops.md`; claim them in the register.
Blocked on: nothing.

## Domain leads

**Legacy Mapper. Dry risk: low.**
1. Output-plan renumbering note for Role 1: steps A and B of the output plan are driver patches `0011` and `0012`, insertion points in `upstream-patches/line-check-output-plan.md`; write the note and keep the patch table in `upstream-patches/README.md` in step (A-63, D-080).
2. WP9b item 9b-1 (chain map of `DIVERGENCE_PART_1`, confirmed against `run_divergence_part1`) and 9b-10 (patch 0009 line shift, sidecar line check) (`amrex/wp3-generator-coverage.md` section 3.2).
3. Regenerate the loop work list and the pressure-loop table from `loop_claims.csv` (`tools/inventory/loop_work_list.py`) after the L1209, L1272, L0866, L0875 and L0397 closures.
4. L1391 (`VELOCITY_FLUX_CYLINDRICAL` cylindrical branch): add the marker and a bitwise test; translatable now, no owner (`amrex/wp3-generator-coverage.md` WP3 item 1).
5. List of every output routine that loops over meshes (FR-072, ADR-004, `amrex/fr072-output-amrex-notes.md`).
6. Cell-loop batch L1355, L1315, L1316, L1321 and S2 package L0398, L0401 (`amrex/loop-work-list.md` section 3); A-34.
7. Register the radiation sweep as claim L1242 with the hand-written FR-062 tag (D-073 (d)); paste the Species and Radiation R1 sign-offs into the claim register.
8. Keep the input classification list current against the converter output now that `main.cpp` calls it (D-076; `amrex/level-jump-external-wall-check.md`).
9. Revised UP-0003 (persistent scratch) and the density-loop patch for L0876 (`combustion/02-signoffs-patch0003-gather-mp5.md`, A-59).

**Species and Combustion Lead. Dry risk: low.**
1. Sign the D-051 device sign-off for `rho_d_maxloc_fix` as domain lead, on the conditions of D-080 (2): first-occurrence tie rule, tie test with three or more species and exact ties, both last-wins mutants killed, real-data tie frequency, NaN limit documented (WP9b item 9b-2 in `amrex/wp3-generator-coverage.md`).
2. Review the revised UP-0003 and the density-loop patch (L0876, `mass.f90:799-849`) when the Legacy Mapper posts them (A-59).
3. Review the translations of L1272 (SETTLING_VELOCITY), L0866, L0875 and L0397 that have landed on the host (`amrex/loop-work-list.md` section 6; `amrex/blocked-loop-families.md` S1 to S4).
4. Species tagging cases for R4 and the clip-line budget test inputs requested in `combustion/03-combustion-under-refinement-and-tagging.md` section 4 (with Role 3 and the V&V Lead).
5. S3 mask-table definition for L0399 and L0406 (six-flag `WALL_INDEX`); other open loops L0596 and L1273 (`amrex/blocked-loop-families.md` S3).
6. Review of the sweep and wall-kernel results that touch species fields: `Z_TEMP` pad for MP5 (`solid/08-new-wall-kernel-review.md` M1).
7. `CUNNINGHAM` libm entry: confirm the registry entry and the PM-07 callee resolution after `ci_checks` (D-077 (5); commit 0a24a1c).

**Solid Phase Lead. Dry risk: low.**
1. Forced-flip survey on the host to measure `E_class`, signed with the V&V Lead before the first device run (D-077 (1) V3; `solid/09-flip-budget-gate.md`).
2. Flip-injection build and exposure tagger for the gate (`solid/09`; V4 requires every exposure tag at least once in P1).
3. A1 predicate form in the cap design: finite state and `DT_BC` continue as upstream, otherwise status 301, same cap and text on the AMR CPU path; wire the 301 template into the driver abort path as for status 300 (D-077 A1, A2; `solid/06-solve-port-test-design.md`).
4. Follow-ups W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7 (W1 is ruled: abort with the report, D-075 (6)); answer the open section 6 items.
5. O3 (L1402) uniqueness-assertion review with the Wall Loops Engineer (`amrex/blocked-loop-families.md` O3; `WLIST_EXT`, `WLIST_INT` in `zone_sum_order.py`).
6. Claim the loops proposed to this lead: L1470, L1471 (HT3D exchange, T4), L1452, L0817 and the `vege.f90` loops L1330 to L1340; L1488, L1489, L1486 (callee flatten, T3) (`amrex/loop-work-list.md` section 3).
7. L1485 retirement: ready to merge once both Build Chiefs validate patch 0010 and P-3 and N-5 pass on the host build (D-077 (2)); prepare the change (D-065).
8. Review the generated SP2 and SP4 kernels and the wall device driver results as device run B lands.

**Radiation Lead. Dry risk: medium (1 is the long item in progress, 3 and 5 wait on design and V&V, 2, 4 and 6 are short).**
1. Sweep stages 1 to 3 (`radiation/06-fr062-sweep-kernel-design.md` section e): `BR_ILW` slot table and boundary kernels, 3D plane kernel on one box, full angle sequence; bitwise acceptance as listed. In progress.
2. Follow-ups R1 to R4 from `solid/08-new-wall-kernel-review.md` section 7 (tolerance class question for the open-boundary cases, row-ownership sentence in `radiation/05-br-ilw-table-spec.md`, owner rule OQ-S6).
3. Open sweep items: `ILD*`, `CC_IBM`, cylindrical (out of the first release), `STORE_RADIATION_TERMS`, launch latency (`radiation/06` section 8).
4. Radiation under regrid: wall `ILW` across levels, random rotation versus regrid, periodic self-coupling (`radiation/04-radiation-phase4-notes.md`; FR-063).
5. Confirm the ghost `UII`/`UIID` reader list from `test/sweep_ghost_reads.py` with the V&V Lead (D-073 (b), D-077 (6)).
6. Radiation verification variants (A-50) with the V&V Lead.

**V&V Lead. Dry risk: low.**
1. `run_class_tolerance` in `kernel_categories.json` and the `RUNCMP` lines of the harness: 1e-11 relative on the step sequence for runs up to 1,000 steps, provisional; add a device-versus-host run to the evidence when one exists (D-080 (5); `vv/dt_class_tolerance.md`).
2. Device-tier 2-ulp harness (D-070, D-072 (b); `vv/test-plan.md` device tier), including the run-level class-tolerance comparison for `wall_visc_les`.
3. Tests for D-080: zero-lost-cell tag-buffer assertion rows for Role 3 (P3 tagging), the tie test review for `rho_d_maxloc_fix`, the `WALL_LOOP_2` gate on `dec2_obst`, and the pressure tests R4 and R5 (`vv/test-plan.md`).
4. P3-B09 distance d(n) = 4n (frozen velocity) in the test plan and `vv-runs/phase3/phase3_cases.csv`, including the observed front (A-61 confirmed).
5. Flip-budget harness: PASS, FAIL and INCONCLUSIVE outcomes, 95% per run and 99% at phase exit, Level A and Level B, shortfall rule (D-077 (1); `solid/09`, test plan 5.10).
6. Review the gate rows of the Phase 3 moving blob against the two-level numbers; set the decomposition rows to `EXACT_SUMS` on (`vv/test-plan.md` 5.11, 5.12.6).
7. Gate hookup for the radiation kernels: apply the gate-rad patch (`radiation/03-radiation-translation-notes.md`).
8. Rerun the Debug `int_1to2` case and the 13 skipped four-rank Intel cases when the test machine is idle (`vv/baseline_status.md`).
9. Carry D-074 into the case inventory and the FR-003 list (`vv/case_inventory.md`).
10. Held-out kernel reviews as evidence lands (`rho_d_interp`, `dp_div_heat`, `d_z_max`, `del_rho_d_del_z`, and `rho_d_maxloc_fix` after the sign-off; `vv/ported-review.md`); A-46 CSV rerun, A-47 baselines, A-26, A-39, A-42, A-49.

## GPU and build roles

**GPU Generator Engineer. Dry risk: low.**
1. `rho_d_maxloc_fix` tie test with three or more species and exact ties, killing the mutants `maxloc_x_last_wins` and `maxloc_z_last_wins`; record the NaN limit in the kernel note (WP9b items 9b-3 and 9b-4; D-080 (2)).
2. Test-lock wrapper: run the job under `timeout 1800`, one full pass per engineer per hour, FIFO, scoped lock for shared build directories; maintain the lock and break a stale hold only after checking the holder process and logging it, then tell the holder (D-077 (4), D-080 (4)).
3. DIVERGENCE_PART_1 nests: L0365, L0369, L0381 and the small nests L0364, L0371, L0373, L0377 to L0380, L0382, L0383, L0386, with the blockers in `amrex/wp3-generator-coverage.md` and `amrex/s4d-divergence-part1-map.md`. `WALL_LOOP_2` stays on the host (D-080 (3)).
4. Generator change W1: `WLIST`/`NWL` indirection in the wall kernels (unblocks WP11, `amrex/stage1-gpu-spike-plan.md`).
5. Layout-contract feature for the pressure loops: integer extents of an `exact` allocation as arguments to cell-loop kernels (`amrex/pressure-loops-generator-proposal.md`).
6. L1392, L1393 (K,I nest with J fixed) and L1375 (`EVALUATE_RAMP` callee flatten) (`amrex/wp3-generator-coverage.md` WP3).
7. Merge `s5-four` and `s5-wall` onto one branch and regenerate once (`solid/08-new-wall-kernel-review.md` G3); sidecar justification for scatters (G2) and the "every array argument is in the device list" structural test (G1).
8. Trace feature G1, status-301 predicate form G2, and ULP line format and NaN/Inf semantics G3 (`solid/09`; D-077 (1)).
9. Document the fixed-point zone-sum mode in the generator how-to (`amrex/generator-howto.md`) and add the ordering statements for P1 and V1 (`amrex/blocked-loop-families.md`).

**GPU Wall Loops Engineer. Dry risk: medium (five items; 3 to 5 depend on other roles).**
1. Fix `rad_wall_qin_zero` and add `W_SURF_INDEX` to the device-address list so `ci_checks --strict` turns green; the ported entries count for the map only until then (D-077 (3), D-078 (b)).
2. Run-level check that `dt_coupled` kernels are compared at class tolerance, using `run_class_tolerance` of D-080 (5) (`vv/dt_class_tolerance.md`), and confirm `ci_checks` K2-10 passes with the `wall_visc_les` entry.
3. W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7 (face-write refusal as a hard stop in table build and refresh, test additions, `wall_checks.py` builders).
4. SP2 per-target gather, SP4 scratch sums and S2 guarded assertion (conditions from the Solid and Species leads, `solid/04-blocked-loop-signoff.md`).
5. O3 (L1402) uniqueness-assertion review with the Solid Phase Lead; WP8/WP11 wall-table gather on the device with the Integration Lead.

**GPU Mesh Data Loops Engineer. Dry risk: low.**
1. Device run B of the four registered wall nests (`amrex/stage1-gpu-spike-plan.md` 7.13b; `solid/08` section 2).
2. Neighbour-mesh (T4) loops: L1363, L1364, L1366, L1367, L1399 (`amrex/loop-work-list.md`, `amrex/o2-edge-exchange-scope.md`; D-055 excludes `NOM>0` branches).
3. Remaining O2 edge and exchange loops (`amrex/o2-edge-exchange-scope.md`).
4. Device and real-field evidence for the held-out kernels (`rho_d_interp`, `dp_div_heat`: dump `RHO_D` and the `D_Z` table and rerun `WP3_DIFF`; `amrex/ported-candidates.md`), with the Integration Lead.
5. Keep the `ported.toml` entries and the strict `ci_checks` pass current as kernels are added; `rho_d_maxloc_fix` is added only after the D-080 (2) sign-off (D-072 (c)).
6. Wire the fixed-point zone-sum mode into `zone_sum_order` checks (finish its tests and commit the wall-list allow-list change).

**GNU Build Chief. Dry risk: medium (items 4 and 6 wait for other roles).**
1. Validate patch 0010 on GNU and run V&V tests P-3 and N-5 (N-5b on the bypass build) on the host build (D-077 (2); Intel has passed).
2. Rerun validation on the tip after Role 1's changes to `fds_step.f90` and `fds_ghost_bc.f90` and the converter wiring in `main.cpp`; regenerate the Release and Debug `vv-runs/refbin` entries (`vv/environment.md`).
3. Debug check of upstream patch UP-0008, and of the `UP-0011` Debug build fix (`upstream-patches/README.md`).
4. Build and test the pinned `AMReX_CUDA_FASTMATH=OFF` change on GNU after the Integration Lead applies it (D-068).
5. A-27 (GNU HYPRE pin and the second AMReX install).
6. Run the two-level `ns2d_16` case at 1 and 4 ranks on GNU for the cross-compiler comparison (D-078 (c)); line check of driver patches `0011` and `0012` on GNU once Role 1 writes them.

**Intel Build Chief. Dry risk: medium (items 1, 5 and 6 wait for other roles).**
1. `-O2 -fp-model=precise` speed cost measurement and isolation of the cause of the `-prec-div` failures; update the build notes and scripts to the D-079 flags (`vv/environment.md`).
2. N-5b control on Intel with the bypass build, keyed on messages because aborting runs exit with status 1 under Intel MPI (D-078 (d)).
3. Run the 13 skipped four-rank Intel baseline cases and the Debug `int_1to2` rerun with V&V (`vv/baseline_status.md`).
4. A-36 optional ifx compile check for K2 (NFR-045) on the generated kernels, including `wall_visc_les` and the fixed-point zone-sum kernels.
5. Intel-side check of the FASTMATH pin and the radiation gate family once applied.
6. Run the two-level `ns2d_16` case under oneAPI at 1 and 4 ranks and compare with GNU (D-078 (c)).

## Program

**Spec and Program Lead. Dry risk: low.**
1. Carry D-079 and D-080 into the spec (Intel flag rule, converter refusals including `TUNNEL_PRECONDITIONER`, complete tag buffering and the withdrawn 5e-4 limit, `run_class_tolerance`, test-lock rule) and the roadmap (output patches `0011`/`0012`, A-63).
2. Update the roadmap dates from Role 1's two-level milestone and the remaining Role 1 items (non-trivial density, 3+ levels, timing).
3. Track the loop tracker and the "ported" listings (28 kernels set; D-072 (c)); regenerate `tracker/loop-tracker.md` after the Legacy Mapper regenerates the work list.
4. Close A-60, A-61 and the N-5 control repair in the action log; open actions for the D-080 items (tag buffering, `rho_d_maxloc_fix` sign-off, `WALL_LOOP_2` host path, C1 to C3).
5. Fix the remaining `0010`/`0011` output-patch wording in other documents to `0011`/`0012` (A-63).
