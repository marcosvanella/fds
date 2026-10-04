
# Next work per role

Owner rule: nobody sits idle. Each role has an ordered list of ready items, each with a source reference (document path and item id). The Chief Architect owns these lists; the Spec & Program Lead refreshed an earlier version at spec v0.4.35, and this is the Chief Architect's refresh (r4) through D-078. Where another list differs, this one wins. A role that finishes its list tells the Chief Architect. "Dry risk" is the chance the role runs out of ready work within two days (low: five or more ready items; medium: three or four; high: fewer than three or all blocked).

Status as of spec v0.4.38 plus D-078, next-work revision r4.

Decisions since r3: D-078 (libm entry for `wall_visc_les`, the six `gsfv_*` kernels co-signed as ported, Role 1 milestone rulings, patch 0010 result on Intel and the N-5 control repair). Earlier: D-077 (flip-budget gate, NIC>1 guard ordering, reviewer of record for ported status, shared test lock, `CUNNINGHAM` libm entry, interior-only ghost radiation compare), D-076 (level-jump wall check, abort guard, input converter), D-075 (pressure results), D-074 (the 23 pressure-code-0 inputs are FDS-only in every mode), D-073 (radiation sweep design), D-072 (libm list, "ported" definition, solid-phase rulings).

Landed since r3 (items removed from the lists): the fine-box domain-edge fix and the two-level run to round-off (Role 1, S14.1 to S14.3, S12, S13.1: mass and rho*Z drift 1.3e-14 over 40 steps); the D-076 input converter library and its `main.cpp` wiring; the exact fixed-point zone sums in the generator (`red_fixedsum`, serial FDS-order default, order-independent mode on a switch); the `wall_visc_les` kernel for L1358; the host translations of the pressure loops L1211, L1207, L1220 to L1222 and L1209 (L1215 retired); the FDS-derived frozen case (M4) and the pin-aware residual check; patch 0010 with its `NIC_CHECK level0:` line (Intel PASS); the six `gsfv_*` kernels ported on the V&V review (28 kernels set in `ported.toml`); `csmag_32_fishpak` in the kernel check; the Phase 3 moving-blob gate rows; the WP3 coverage note and the WP9b split (`amrex/wp3-generator-coverage.md`).

Numbering of patches (A-63): upstream-FDS patches are `UP-NNNN`; driver patches are the bare numbers `0001` to `0012`. The level-0 guard is `0010`; the output patches are `0011` and `0012`.

Shared references used below: the ranked loop list `amrex/loop-work-list.md` (and `amrex/loop_work_list.csv`; loop ids `Lnnnn`), the stage-1 plan `amrex/stage1-gpu-spike-plan.md` (work packages `WPn`), the blocked-loop families `amrex/blocked-loop-families.md` (ids P1 to P4, S1 to S4, SP1 to SP4, R1, V1, V2, O1 to O4), the WP3 and WP9b note `amrex/wp3-generator-coverage.md`, the loop tracker `tracker/loop-tracker.md`, the ported candidates `amrex/ported-candidates.md`, and the action log in `README.md` (ids `A-nn`).

Summary of dry risk: Role 1 low (overloaded, not dry); Role 2 low; Role 3 low; Integration Lead low; Pressure Lead low; Legacy Mapper low; Species and Combustion low; Solid Phase low; Radiation medium; V&V low; Spec Lead low; Generator Engineer low; Wall Loops Engineer medium; Mesh Data Loops Engineer low; GNU Build Chief medium; Intel Build Chief medium. No role is high. Items that wait for another role are counted as not ready when the level is set.

## Data layout, pressure and transport implementers

**Role 1, Data Layout. Dry risk: low (the list is long and the order matters).**
1. Explain the fine-mesh `VELOCITY_BC` effect: after `exchange(6)` the fine mesh overwrote level-0 periodic U/W ghosts, so `fill_omesh` and `VELOCITY_BC` are skipped on level>0 as an interim (`FDSTL_SKIPAFT=0` restores). Give a hypothesis and a reproducer before the Phase 3 gate (D-078 (c)).
2. Overwrite-off control for the two-level `ns2d_16` run: with the interface flux overwrite off the conservation drift must come back, so the control shows the overwrite is what holds it (D-050; `vv/test-plan.md` 5.12.6).
3. Multi-box fine level (more than one fine box) and the 4-rank two-level run, with `EXACT_SUMS` on for decomposition comparisons (D-053, `vv/test-plan.md` 5.11).
4. A test that triggers the `stage_boundary(1,3|6)` guard on a fine box (D-078 (c); A-60). The guard is unverified until a case reaches it.
5. Remaining WP1b dump points: `CE_EDGE_INDEX`, `ED_OMEGA`, `ED_TAU` (`amrex/stage1-gpu-spike-plan.md` WP1b). There is no `WALL_UVW_INTERP` routine in this tree: record the routine that replaces it with the Integration Lead.
6. Output patches `0011` and `0012` (driver series) once the Legacy Mapper has checked their line numbers; the earlier line check covered patch 0010 and `UP-0011` only (`upstream-patches/line-check-0010-0011.md`; `roadmap.md`, output patches row).
7. FR-072 output plan with the Legacy Mapper's list of mesh-looping output routines (ADR-004, `amrex/fr072-output-amrex-notes.md`).
8. Check the E1 predictor `DIV1` gap against the two-level numbers (`amrex/stage1-gpu-spike-plan.md` WP1b); close it or record what remains.
Blocked on: item 6 waits for the Legacy Mapper. Nothing else.

**Role 2, Pressure Backend. Dry risk: low.**
1. Upstream issue draft for the AMReX row-scaling problem (add the HYPRE row-scaling note with the `habec_ijmat` symmetric-scaling candidate to `upstream-issues/`, D-075 (4)). Draft only, not sent.
2. Open `pres.f90` loops of the work list: L1209 on OPEN walls against real FDS (the verbatim-Fortran test passes; find why the dump comparison differs), L1210 (cylindrical PRHS), L1212 to L1214 (transposed PRHS), L1216 to L1218 (transposed copies; retire like L1215) (`amrex/pressure-loops-r2-status.md`, "Left").
3. Layout contract for L1211 and L1212 to L1214 with the Generator Engineer: the extents of an `exact` allocation (ITRN, JTRN, KTRN) as integer arguments to cell-loop kernels (`amrex/pressure-loops-generator-proposal.md`).
4. A second FDS-derived frozen case beyond `csmag_32` periodic step 3: two meshes with wind and a pressure ramp, and an open-boundary case, so the fold is checked with nonzero wall data (`pressure/07-a56-backend-comparison-plan.md` section 11.4, A-09b).
5. "Hierarchy changed, rebuild" entry point and the composite solve with an arbitrary per-level right-hand side, now exercised by Role 1's two-level run; confirm with Role 1 that nothing else is missing (D-058, D-047).
6. Composite solve at 4 ranks and with a multi-box fine level (`pressure_backend` tests; needed by Role 1 item 3).
Blocked on: nothing.

**Role 3, Regrid and Multi-Level Transport. Dry risk: low (the domain-edge fix has landed).**
1. Tighten `MASS_TOL` and `E1_ZZ_TOL` to round-off against Role 1's first two-level numbers (mass and rho*Z drift 1.3e-14 over 40 steps; E3b drift 6e-16, E1 RHO/ZZ 3e-16; A-60, D-071 (g)).
2. Registry transfer check against the D-078 (c) rule: a new fine level takes `RSUM`, `MU`, `KRES`, `D`, `DS`, `H` and `HS` from the parent; add a test with a zero-initialised fine level that fails without the injection.
3. N-5 control input that passes the D-076 converter (an equal-level-0 layout with a NIC>1 face), with the V&V Lead, or a documented scratch build that bypasses the converter (D-078 (d); `upstream-patches/inputs/0010-neg-control-race_test_1.fds`).
4. Two-level driver runs of the FR-016 gate cases on `--two-level-run`, then the moving blob (A-58 thresholds, `vv/test-plan.md` 5.12.6; the gate rows P3-B01 to B05 are in `test_blob_registry --gates`).
5. Larger-case tag-buffer assertions (`role3-regrid-transport-plan.md`, tagging section; V&V P3 tagging cases).
6. t = 0 hierarchy, tagging and regrid on the real FDS stages (FR-012).
7. Offload tagging test on the test machine (`RT_TAG_OFFLOAD`, CUDA AMReX install; Phase 4 item).
Blocked on: item 4 needs the multi-box fine level from Role 1 for the box-split rows.

**AMReX Integration Lead. Dry risk: low.**
1. FASTMATH CMake pin: apply `tools/patches/cmake-fastmath-pin.patch` to the generator worktree (the scratch-copy check passes BF-05 and BF-07, `amrex/ported-candidates.md` section 3; D-068).
2. Held-out kernels: real-field runs for `rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat` and a host-tier result for `d_z_max` and `del_rho_d_del_z`, then add them to `ported.toml` (`amrex/ported-candidates.md`, held out).
3. Sweep launch layer: one launch per wavefront plane across all boxes of a level (D-073 (f), D-075 (7); stage-1 plan 7.14), with the 8-thread sweep run and the `-gpu=nofma` build on the test machine.
4. WP8 flat tables on the device and WP11 wall kernels in the driver (after generator change W1, `WLIST`/`NWL`), `amrex/stage1-gpu-spike-plan.md`.
5. WP9b pieces 9b-5 and 9b-7: device launch sequence part A, hook, part B, and the decision on the host pieces of part B (`STORE_SPECIES_FLUX` copies, `WALL_LOOP_2`) (`amrex/wp3-generator-coverage.md` section 3.2).
6. Rule-7 macro form in `s4_omp.inc` and re-time (`amrex/s4-k2-single-source.md`; D-072 (d)).
7. Open actions A-28, A-31, A-44, A-48, A-49, A-51, A-54 (`README.md`).
8. Refresh the kernel map after the `wall_visc_les` registration and the 28 ported entries (`inventory/port_kernel_status.md`).
Blocked on: GPU time on the test machine for items 3 and 4.

**Pressure Solver Lead. Dry risk: low (six ready items, all unblocked).**
1. A-56 backend comparison to the end: composite and masked cases with the HYPRE backend, report in `pressure/07` and `pressure/06`; delete regenerable scratch afterwards.
2. Second case family for P3-R02 so the bound can be reviewed for tightening (D-075 (3)).
3. Non-matching Dirichlet wall data test on the hierarchy (D-075 (2): untested).
4. Sign-off review of Role 2's host translations of L1211, L1207, L1220 to L1222 and L1209 against `pres.f90` (`amrex/pressure-loops-r2-status.md`; `amrex/pressure-velocity-h-loops.md`).
5. Single-thread device add measurement for zone sums, now with the fixed-point mode of the generator in place (stage-1 plan section on P1/V1 summation; `tools/README.md` zone-sum section).
6. Loop candidates for the pressure domain: L1365 and the other `pres.f90` loops of `amrex/pressure-velocity-h-loops.md`; claim in the register.
Blocked on: nothing.

## Domain leads

**Legacy Mapper. Dry risk: low.**
1. Line-number check for the output patches `0011` and `0012` (unblocks Role 1 item 6); the earlier check covered `0010` and `UP-0011` only.
2. WP9b item 9b-1 (chain map of `DIVERGENCE_PART_1`, confirmed against `run_divergence_part1`) and 9b-10 (patch 0009 line shift, sidecar line check) (`amrex/wp3-generator-coverage.md` section 3.2). Open question there: whether item 9b-1 is the intended split.
3. Regenerate the loop work list and the pressure-loop table from `loop_claims.csv` (`tools/inventory/loop_work_list.py`); the pressure statuses were written directly in the claims file (`amrex/pressure-loops-r2-status.md`, "Left").
4. L1391 (`VELOCITY_FLUX_CYLINDRICAL` cylindrical branch): add the marker and a bitwise test; it is translatable now and has no owner (`amrex/wp3-generator-coverage.md` WP3 item 1).
5. List of every output routine that loops over meshes (FR-072, ADR-004, `amrex/fr072-output-amrex-notes.md`).
6. Cell-loop batch L1355, L1315, L1316, L1321 and S2 package L0398, L0401 (`amrex/loop-work-list.md` section 3); A-34.
7. Register the radiation sweep as claim L1242 with the hand-written FR-062 tag (D-073 (d)); paste the Species and Radiation R1 sign-offs into the claim register.
8. Validation follow-up: check the validation inputs against the converter output now that `main.cpp` calls it; keep the input classification list current (D-076; `amrex/level-jump-external-wall-check.md`).
9. Revised UP-0003 (persistent scratch) and the density-loop patch for L0876 (`combustion/02-signoffs-patch0003-gather-mp5.md`, A-59).

**Species and Combustion Lead. Dry risk: low (six ready items; item 1 waits for the Legacy Mapper).**
1. Review the revised UP-0003 and the density-loop patch (L0876, mass.f90:799-849) when the Legacy Mapper posts them (A-59).
2. L1272 SETTLING_VELOCITY (1.325 %, S4): per-species table and the split into three kernels, using the three additive generator features (function callee, rank-2 table through a `PRESSURE_ZONE` cell array, `CYCLE` of the kernel's own innermost loop) (`amrex/loop-work-list.md` section 6).
3. Species tagging cases for R4 and the clip-line budget test inputs requested in `combustion/03-combustion-under-refinement-and-tagging.md` section 4 (with Role 3 and the V&V Lead).
4. S3 mask-table definition for L0399 and L0406 (six-flag `WALL_INDEX`; blocked-loop families S3); other open loops: L0397, L0866, L0875, L0596, L1273.
5. Sign-off for the `rho_d_maxloc_fix` kernel (`MAXLOC` first-maximum rule; no note names it, WP9b item 9b-2 in `amrex/wp3-generator-coverage.md`), including a tie case with at least three species.
6. `CUNNINGHAM` libm entry: append to `docs/tools/kernel_registry.toml` by path (the draft entry with the grep check is in the working copy), then run `ci_checks` (D-077 (5)).
7. Review of the sweep and wall-kernel results that touch species fields: `Z_TEMP` pad for MP5 (`solid/08-new-wall-kernel-review.md` M1).

**Solid Phase Lead. Dry risk: low.**
1. Forced-flip survey on the host to measure `E_class`, signed with the V&V Lead before the first device run (D-077 (1) V3; `solid/09-flip-budget-gate.md`).
2. Flip-injection build and exposure tagger for the gate (`solid/09`; V4 requires every exposure tag at least once in P1).
3. A1 predicate form in the cap design: finite state and `DT_BC` continue as upstream, otherwise status 301, same cap and text on the AMR CPU path; wire the 301 template into the driver abort path as for status 300 (D-077 A1, A2; `solid/06-solve-port-test-design.md`).
4. Follow-ups W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7 (W1 is ruled: abort with the report, D-075 (6)); answer the open section 6 items.
5. O3 (L1402) uniqueness-assertion review with the Wall Loops Engineer (`amrex/blocked-loop-families.md` O3). The wall-list allow-list in `zone_sum_order.py` (`WLIST_EXT`, `WLIST_INT`) is the mechanism under review.
6. Claim the loops proposed to this lead: L1470, L1471 (HT3D exchange, T4), L1452, L0817 and the `vege.f90` loops L1330 to L1340 (`amrex/loop-work-list.md` section 3, Wall split); L1488, L1489, L1486 (callee flatten, T3).
7. Review the generated SP2 and SP4 kernels as the device run B results land.
8. L1485 retirement: ready to merge once both Build Chiefs validate patch 0010 and P-3 and N-5 pass on the host build (D-077 (2)); prepare the change (D-065).

**Radiation Lead. Dry risk: medium (six items, but 1 is the long item in progress, 3 and 5 wait on design and V&V, 4 and 6 are short).**
1. Sweep stages 1 to 3 (`radiation/06-fr062-sweep-kernel-design.md` section e): `BR_ILW` slot table and boundary kernels, 3D plane kernel on one box, full angle sequence; bitwise acceptance as listed. Stages 1 to 3 are in progress.
2. Follow-ups R1 to R4 from `solid/08-new-wall-kernel-review.md` section 7 (tolerance class question for the open-boundary cases, row-ownership sentence in `radiation/05-br-ilw-table-spec.md`, owner rule OQ-S6).
3. Radiation under regrid: wall `ILW` across levels, random rotation versus regrid, periodic self-coupling (`radiation/04-radiation-phase4-notes.md`; FR-063).
4. Confirm the ghost `UII`/`UIID` reader list from `test/sweep_ghost_reads.py` with the V&V Lead (D-073 (b), D-077 (6)).
5. Radiation verification variants (A-50) with the V&V Lead.
6. `-gpu=nofma` build and the 8-thread sweep run on the test machine; restart checks compare ghost `UII`/`UIID`/`QR` interior-only (D-077 (6)).

**V&V Lead. Dry risk: low.**
1. Device-tier 2-ulp harness (D-070, D-072 (b); `vv/test-plan.md` device tier), including the run-level class-tolerance comparison for `dt_coupled` kernels such as `wall_visc_les` (D-078 (a)).
2. N-5 control repair with Role 3: a converter-passing variant with a NIC>1 face; tests key on messages, not exit status, and the Intel exit status of aborting runs is 1 (D-078 (d); test plan P-3, N-5).
3. P3-B09 distance d(n) = 4n (frozen velocity) into the test plan and `vv-runs/phase3/phase3_cases.csv`, including the observed front (A-61 confirmed).
4. Flip-budget harness: PASS, FAIL and INCONCLUSIVE outcomes, 95% per run and 99% at phase exit, Level A and Level B, shortfall rule (D-077 (1); `solid/09`, test plan 5.10).
5. Review the gate rows of the Phase 3 moving blob and the new two-level numbers; set the decomposition rows to `EXACT_SUMS` on (`vv/test-plan.md` 5.11, 5.12.6).
6. Gate hookup for the radiation kernels: apply the gate-rad patch to the gate script (adds the `rad` kernel family and its `KERNEL` result lines; `radiation/03-radiation-translation-notes.md`).
7. Rerun the Debug `int_1to2` case and the 13 skipped four-rank Intel cases when the test machine is idle (`vv/baseline_status.md`).
8. Carry D-074 into the case inventory and the FR-003 list (`vv/case_inventory.md`); the seven NIC>1 inputs per `amrex/level-jump-external-wall-check.md`.
9. Held-out kernel reviews as evidence lands (`rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat`, `d_z_max`, `del_rho_d_del_z`; `vv/ported-review.md`).
10. A-46 CSV rerun, A-47 baselines, A-26, A-39, A-42, A-49.

## GPU and build roles

**GPU Generator Engineer. Dry risk: low.**
1. DIVERGENCE_PART_1 nests: L0365, L0369, L0381 and the small nests L0364, L0371, L0373, L0377 to L0380, L0382, L0383, L0386, with the blockers in `amrex/wp3-generator-coverage.md` (`RHO_D = MAX(0,MU)*RSC_T` whole-array forms; `WALL_LOOP_2` array sections; `DEL_RHO_D_DEL_Z = 0`) and `amrex/s4d-divergence-part1-map.md`.
2. WP9b items 9b-3 and 9b-4: close or document the NaN limit of `rho_d_maxloc_fix`, and the `del_rho_d_del_z` kernel in the split chain with `DP` re-initialised in part A.
3. Generator change W1: `WLIST`/`NWL` indirection in the wall kernels (unblocks WP11, `amrex/stage1-gpu-spike-plan.md`).
4. Layout-contract feature for the pressure loops: integer extents of an `exact` allocation as arguments to cell-loop kernels (`amrex/pressure-loops-generator-proposal.md`, L1211 and L1212 to L1214).
5. L1392, L1393 (K,I nest with J fixed) and L1375 (`EVALUATE_RAMP` callee flatten); L1391 is the Legacy Mapper's marker item (`amrex/wp3-generator-coverage.md` WP3 item 5).
6. Merge `s5-four` and `s5-wall` onto one branch and regenerate once (`solid/08-new-wall-kernel-review.md` G3); sidecar justification for scatters (G2) and the "every array argument is in the device list" structural test (G1).
7. Trace feature G1, status-301 predicate form G2, and ULP line format and NaN/Inf semantics G3 (`solid/09`; D-077 (1)).
8. Three additive `s5gen.py` features for SETTLING_VELOCITY (function callee, rank-2 table gathered through a `PRESSURE_ZONE` cell array, `CYCLE` of the kernel's own innermost loop).
9. Shared test-lock implementation: one full pass per engineer per hour under `flock`, FIFO, 30-minute hold limit, scoped lock for shared build directories, stale-lock rule (D-077 (4)).
10. Document the fixed-point zone-sum mode in the generator how-to (`amrex/generator-howto.md`) and add the ordering statements for P1 and V1 (`amrex/blocked-loop-families.md`).

**GPU Wall Loops Engineer. Dry risk: medium (five items; 3 to 5 depend on other roles).**
1. `wall_visc_les` libm registry entry: append it to `docs/tools/kernel_registry.toml` by path with the status of an approved entry (not "proposed"), add the check that `dt_coupled` kernels are compared at class tolerance at run level, and run `ci_checks` so K2-10 passes (D-078 (a)).
2. Fix `rad_wall_qin_zero` and add `W_SURF_INDEX` to the device-address list so `ci_checks --strict` turns green; the ported entries count for the map only until then (D-077 (3), D-078 (b)).
3. W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7 (face-write refusal as a hard stop in table build and refresh, test additions, `wall_checks.py` builders).
4. SP2 per-target gather, SP4 scratch sums and S2 guarded assertion (conditions from the Solid and Species leads, `solid/04-blocked-loop-signoff.md`).
5. O3 (L1402) uniqueness-assertion review with the Solid Phase Lead; WP8/WP11 wall-table gather on the device with the Integration Lead.

**GPU Mesh Data Loops Engineer. Dry risk: low.**
1. Device run B of the four registered wall nests (`amrex/stage1-gpu-spike-plan.md` 7.13b; `solid/08` section 2).
2. Neighbour-mesh (T4) loops: L1363, L1364, L1366 (in progress), L1367, L1399 (`amrex/loop-work-list.md`, `amrex/o2-edge-exchange-scope.md`; D-055 excludes `NOM>0` branches).
3. Remaining O2 edge and exchange loops (`amrex/o2-edge-exchange-scope.md`).
4. Device and real-field evidence for the held-out kernels (`rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat`: dump `RHO_D` and the `D_Z` table and rerun `WP3_DIFF`; `amrex/ported-candidates.md`), with the Integration Lead.
5. Keep the `ported.toml` entries and the strict `ci_checks` pass current as kernels are added; kernel-map snapshot (D-072 (c)).
6. Wire the fixed-point zone-sum mode into `zone_sum_order` checks (the working copy has the wall-list allow-list change; finish its tests and commit).

**GNU Build Chief. Dry risk: medium (items 4 and 6 wait for other roles).**
1. Validate patch 0010 on GNU and run V&V tests P-3 and N-5 on the host build (D-077 (2)); Intel has passed (D-078 (d)).
2. Debug check of upstream patch UP-0008 (`upstream-patches/README.md`).
3. Keep the Release and Debug reference binaries current for the merged tip, including the converter wiring in `main.cpp` (`vv/environment.md`); regenerate `vv-runs/refbin` entries on a new tip.
4. Build and test the pinned `AMReX_CUDA_FASTMATH=OFF` change on GNU after the Integration Lead applies it (D-068).
5. A-27 (GNU HYPRE pin and the second AMReX install).
6. Run the two-level `ns2d_16` case at 1 and 4 ranks on GNU once Role 1's multi-box and 4-rank runs exist (cross-compiler comparison; D-078 (c)).

**Intel Build Chief. Dry risk: medium (items 1, 5 and 6 wait for other roles).**
1. N-5 control on Intel: run it with the repaired input or the converter-bypass scratch build, keyed on messages because aborting runs exit with status 1 under Intel MPI (D-078 (d)).
2. `-prec-div` speed cost measurement and build notes (D-069; `vv/environment.md`).
3. Run the 13 skipped four-rank Intel baseline cases and the Debug `int_1to2` rerun with V&V (`vv/baseline_status.md`).
4. A-36 optional ifx compile check for K2 (NFR-045) on the generated kernels, including the new `wall_visc_les` and the fixed-point zone-sum kernels.
5. Intel-side check of the FASTMATH pin and the radiation gate family once applied.
6. Run the two-level `ns2d_16` case under oneAPI at 1 and 4 ranks and compare with GNU (D-078 (c)).

## Program

**Spec and Program Lead. Dry risk: low.**
1. Carry D-076, D-077 and D-078 into the spec (FR-003 scope, flip-budget gate, NIC>1 guard ordering, reviewer of record, test lock, `dt_coupled` run-level tolerance, `DERIVE` rule for a new fine level) and the roadmap (output patches `0011`/`0012`).
2. Chase the Role 1 open items (overwrite-off control, multi-box and 4-rank two-level runs, the `VELOCITY_BC` cause before the Phase 3 gate) and update the roadmap dates.
3. Track the loop tracker and the "ported" listings (28 kernels set; D-072 (c)); regenerate `tracker/loop-tracker.md` after the Legacy Mapper regenerates the work list.
4. Close A-60 and A-61 in the action log; open actions for the N-5 control repair and the `wall_visc_les` entry.
5. Fix the remaining `0010`/`0011` output-patch wording in other documents to `0011`/`0012` (A-63).
