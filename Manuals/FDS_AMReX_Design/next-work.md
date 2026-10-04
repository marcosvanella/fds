
# Next work per role

Owner rule: nobody sits idle. Each role has an ordered list of ready items, each with a source reference (document path and item id). The Chief Architect owns these lists; the Spec & Program Lead refreshed the previous version at spec v0.4.35, and this is the Chief Architect's refresh (r3) through D-077. Where another list differs, this one wins. A role that finishes its list tells the Chief Architect. "Dry risk" is the chance the role runs out of ready work within two days (low: five or more ready items; medium: three or four; high: fewer than three or all blocked).

Status as of spec v0.4.37 plus D-076 and D-077, next-work revision r3.

Decisions since the previous list: D-076 (level-jump wall check, abort guard, input converter) and D-077 (flip-budget gate ratified, NIC>1 guard ordering, V&V Lead as reviewer of record for ported status, shared test lock, `CUNNINGHAM` libm entry, interior-only ghost radiation compare). Earlier: D-072 (libm list, "ported" definition, solid-phase rulings), D-073 (radiation sweep kernel design), D-074 (the 23 pressure-code-0 inputs are FDS-only in every mode; `fds_amr` aborts naming the FDS executable), D-075 (pressure results: fold sign confirmed, mixed Neumann/Dirichlet excess explained, P3-R02 bound kept, HYPRE Krylov defaults, face-write table abort, launch layer owner).

Done recently: patch 0010 (NIC>1 guard) drafted with its negative control; test plan v0.4.19 with P-3 and N-5; `ported.toml` aligned with the V&V review; Role 3's input converter library; the composite two-level MLMG solve and the HYPRE assembled backend; patches 0005 to 0009 validated on oneAPI and GNU Debug; the four mass-flux wall nests registered (86 to 90 kernels, 92 across the `s5-four` and `s5-wall` branches); Species and Radiation R1 sign-offs; Intel Release baseline; P3-R02 measured 7,000 to 62,000 times below its bound; shared-disk cleanup.

Shared references used below: the ranked loop list `amrex/loop-work-list.md` (and `amrex/loop_work_list.csv`; loop ids `Lnnnn`), the stage-1 plan `amrex/stage1-gpu-spike-plan.md` (work packages `WPn`), the blocked-loop families `amrex/blocked-loop-families.md` (ids P1 to P4, S1 to S4, SP1 to SP4, R1, V1, V2, O1 to O4), and the action log in `README.md` (ids `A-nn`).

Summary of dry risk: Role 1 low (overloaded, not dry); Role 2 low; Role 3 medium; Integration Lead low; Pressure Lead medium; Legacy Mapper low; Species and Combustion medium; Solid Phase low; Radiation medium; V&V low; Spec Lead low; Generator Engineer low; Wall Loops Engineer medium; Mesh Data Loops Engineer low; GNU Build Chief medium; Intel Build Chief medium.

## Data layout, pressure and transport implementers

**Role 1, Data Layout. Dry risk: low (the list is long and the order matters).**
1. Fine-box physical-domain ghost fill and wall cells on every domain-edge face: level-0 treatment, periodic and Neumann or mirror only, abort with a clear message otherwise (`README.md` A-60, D-071 (g), R-65). First priority; it gates Role 3's tolerance tightening.
2. `D_PBAR_DT` level-1 zone binding (`bind_level`, Phase 3 two-level transport test).
3. `stage_boundary(1,3|6)` abort for the unsupported boundary faces on fine boxes (A-60; same message style as item 1).
4. E1 predictor `DIV1` gap (`amrex/stage1-gpu-spike-plan.md` WP1b, Role 1's dump-point notes).
5. Call `install_post_regrid_projection` from the time loop (D-063; `POST_REGRID_PROJECTION = AUTO | ON | OFF`).
6. Apply the Integration Lead's WP1b dump-point patch (`amrex/stage1-gpu-spike-plan.md` WP1b).
7. Output patches 0010 and 0011, after the Legacy Mapper's line-number check (`roadmap.md`, output patches row).
8. FR-072 output plan with the Legacy Mapper's list of mesh-looping output routines (ADR-004, `amrex/fr072-output-amrex-notes.md`).
9. Apply Role 3's `main.cpp` input-converter patch (`Source/regrid_transport/notes/driver-patch-main-input-converter.patch`, Role 3 commit `6d8140a714`) after the domain-edge fix; `main.cpp` then calls `parse_amr_params` and the hierarchy builder and gives `assemble_level0` only level-0 meshes (D-076, D-077).
Blocked on: item 7 waits for the Legacy Mapper. Nothing else.

**Role 2, Pressure Backend. Dry risk: low.**
1. M4: FDS right-hand-side dump patch and an FDS-derived frozen case, so the fold and the solve are compared against FDS data and not only the dense matrix (`pressure/07-a56-backend-comparison-plan.md` section 11.4, A-09b, D-075 (1)). Delivered as a patch file (D-051).
2. Pin-row residual handling: exclude the pin row from the residual check and report it separately; `residual_tol` is not relaxed (D-075 (5)).
3. Loop translation: the pressure-side items proposed to the Pressure Lead in `amrex/pressure-velocity-h-loops.md` and the open `pres.f90` loops in `amrex/loop-work-list.md` (L1211 to L1222, L1206 to L1208 are the largest open ones); claim before starting (claim protocol, section 4 of the list).
4. Upstream issue draft for the AMReX row-scaling problem (PoissonHybrid singular mean is in `upstream-issues/`; add the HYPRE row-scaling note with the `habec_ijmat` symmetric-scaling candidate, D-075 (4)). Draft only, not sent.
5. "Hierarchy changed, rebuild" entry point (D-058; host rebuild allowed until Phase 11, D-047) and the composite solve with an arbitrary per-level right-hand side date (gates the full two-level `ns2d_16` run).
6. 2-D and singular-case mapping check (D-057).
Blocked on: nothing.

**Role 3, Regrid and Multi-Level Transport. Dry risk: medium (items 1 and 4 wait for Role 1).**
1. Confirm the P3-B09 (1) propagation distance d(n) = 4n against the driver ghost depth and extra reads (A-61); hand the value to the V&V Lead.
2. Tighten `MASS_TOL` and `E1_ZZ_TOL` to round-off once Role 1's fine-box domain-edge fix lands (A-60, D-071 (g)).
3. Larger-case tag-buffer assertions (`role3-regrid-transport-plan.md`, tagging section; V&V P3 tagging cases).
4. Align the plan wording to interface flux overwrite, no reflux, FR-016 (a) at step 1 (A-58 open part).
5. Two-level driver runs of the FR-016 gate cases, then the moving blob (A-58 thresholds, `vv/test-plan.md` section 5.12.6).
6. t = 0 hierarchy, tagging and regrid on the real FDS stages (FR-012).
Blocked on: Role 1 (items 2 and 5 need the edge fix and `bind_level`).

**AMReX Integration Lead. Dry risk: low.**
1. FASTMATH CMake pin: review and apply the patch (D-068, lint rules BF-05 and BF-07 fail until it is applied; `docs/tools/` registry).
2. Per-kernel tolerance class in the device harness (D-070, D-072 (b); `amrex/stage1-gpu-spike-plan.md` 7.10d), then set "ported" status in the kernel map per D-072 (c).
3. Sweep launch layer: one launch per wavefront plane across all boxes of a level (D-073 (f), D-075 (7); `radiation/06-fr062-sweep-kernel-design.md` option A).
4. WP8 flat tables on the device (owner proposal: Integration Lead) and WP11 wall kernels in the driver (after generator change W1, `WLIST`/`NWL`), `amrex/stage1-gpu-spike-plan.md`.
5. Rule-7 macro form in `s4_omp.inc` and re-time (`amrex/s4-k2-single-source.md`; D-072 (d)).
6. Open actions A-28, A-31, A-44, A-48, A-49, A-51, A-54 (`README.md`).
7. Update the `ported.toml` co-signature text to the D-077 (3) wording and keep the evidence scope in each note.
8. Add the six `gsfv_*` entries to `ported.toml` when the V&V Lead's review is accepted and listed (D-077 (3)).
9. Own the batched radiation sweep launch layer, with the 8-thread sweep run and the `-gpu=nofma` build on the GPU test machine (D-073, D-077 (6); stage-1 plan 7.14).
Blocked on: GPU time on the test machine for items 2 and 4.

**Pressure Solver Lead. Dry risk: medium.**
1. A-56 backend comparison to the end: composite and masked cases with the HYPRE backend, report in `pressure/07` and `pressure/06`; delete regenerable scratch afterwards.
2. Second case family for P3-R02 so the bound can be reviewed for tightening (D-075 (3)).
3. Non-matching Dirichlet wall data test on the hierarchy (D-075 (2): untested).
4. Single-thread device add measurement for zone sums (stage-1 plan section on P1/V1 summation).
5. Loop candidates for the pressure domain: L1365 and the `pres.f90` loops listed in `amrex/pressure-velocity-h-loops.md`; claim in the register.
Blocked on: nothing.

## Domain leads

**Legacy Mapper. Dry risk: low.**
1. Driver-code check that the AMR route never builds `EXTERNAL_WALL` at a coarse/fine level jump, and the AMR-mode status of the seven multi-mesh NIC>1 inputs; check the validation inputs (D-072 (f); `solid/05-fine-level-solid-plan.md`, `amrex/loop-work-list.md` L1485).
2. WP3 generator coverage in the order of the generator answers (`amrex/stage1-gpu-spike-plan.md` WP3, `amrex/stage1-generator-answers.md`) and WP9b (split the divergence chain at the DIF hook).
3. Register the radiation sweep as claim L1242 with the hand-written FR-062 tag (D-073 (d)); paste the Species and Radiation sign-offs into the claim register (A-59).
4. Line-number check for output patches 0010 and 0011 (unblocks Role 1).
5. List of every output routine that loops over meshes (FR-072, ADR-004, `amrex/fr072-output-amrex-notes.md`).
6. Cell-loop batch L1355, L1315, L1316, L1321 and S2 package L0398, L0401 (`amrex/loop-work-list.md` section 3); A-34.
7. Patch 0010 hook: the `NIC_CHECK level0:` assertion line ships in patch 0010; confirm the guard exits with status 0 and keys on `ERROR(9001)` (D-077 (2), `upstream-patches/0010-main-amr-level0-nic-guard.md`).
8. Validation follow-up: check the validation inputs against the converter output once Role 1 applies the patch; keep the input classification list current (D-076).

**Species and Combustion Lead. Dry risk: medium.**
1. Review the revised patch 0003 and the density-loop patch (L0876, mass.f90:799-849) when the Legacy Mapper posts them (`combustion/02-signoffs-patch0003-gather-mp5.md`, A-59).
2. Take the open loops proposed to this lead: L1272 SETTLING_VELOCITY (1.325 %, S4, needs a per-species table and a split into three kernels), L0397, L0866, L0875 (generator gaps), L0596 and L1273 (`amrex/loop-work-list.md` section 6).
3. Species tagging cases for R4 and the clip-line budget test inputs requested in `combustion/03-combustion-under-refinement-and-tagging.md` section 4 (with Role 3 and the V&V Lead).
4. S3 mask-table definition for L0399 and L0406 (six-flag `WALL_INDEX`; blocked-loop families S3).
5. Review of the sweep and wall-kernel results that touch species fields: `Z_TEMP` pad for MP5 (`solid/08-new-wall-kernel-review.md` M1).
6. Append the `CUNNINGHAM` libm entry to `docs/tools/kernel_registry.toml` by path (`func.f90:2129`, 2 ulp, `dt_coupled` false), with a recorded grep check that the settling velocity does not enter the time-step estimate, then run `ci_checks` (D-077 (5)).
7. SETTLING_VELOCITY landing (L1272): use the three additive generator features from the Generator Engineer (function callee, rank-2 table gathered through a `PRESSURE_ZONE` cell array, `CYCLE` of the kernel's own innermost loop).

**Solid Phase Lead. Dry risk: low.**
1. Fine-level solid plan follow-ups under D-064, D-065 and D-072 (e): flip-budget gate definition with the V&V Lead, 1-D solve cap check in the harness (`solid/05-fine-level-solid-plan.md`, `solid/06-solve-port-test-design.md`).
2. Follow-ups W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7 (W1 is now ruled: abort with the report, D-075 (6)); answer the open section 6 items.
3. O3 (L1402) uniqueness-assertion review, together with the Wall Loops Engineer (`amrex/blocked-loop-families.md` O3).
4. Claim the loops proposed to this lead: L1470, L1471 (HT3D exchange, T4), L1452, L0817 and the `vege.f90` loops L1330 to L1340 (`amrex/loop-work-list.md` section 3, Wall split).
5. Review the generated SP2 and SP4 kernels as the device run B results land.
6. Forced-flip survey on the host to measure `E_class`, signed with the V&V Lead before the first device run (D-077 (1) V3; `solid/09-flip-budget-gate.md`).
7. Flip-injection build and exposure tagger for the gate (`solid/09`; V4 requires every exposure tag at least once in P1).
8. Put the A1 predicate form in the cap design: finite state and `DT_BC` continue as upstream, otherwise status 301, with the same cap and text on the AMR CPU path (D-077 A1).
9. Wire the 301 message template into the driver abort path as for status 300 (D-077 A2; `solid/06-solve-port-test-design.md`).

**Radiation Lead. Dry risk: medium.**
1. Sweep stages 1 to 3 (`radiation/06-fr062-sweep-kernel-design.md` section e): `BR_ILW` slot table and boundary kernels, 3D plane kernel on one box, full angle sequence; bitwise acceptance as listed.
2. Follow-ups R1 to R4 from `solid/08-new-wall-kernel-review.md` section 7 (tolerance class question for the open-boundary cases, row-ownership sentence in `radiation/05-br-ilw-table-spec.md`, owner rule OQ-S6).
3. Radiation under regrid: wall `ILW` across levels, random rotation versus regrid, periodic self-coupling (`radiation/04-radiation-phase4-notes.md`, open items for the radiation ADR; FR-063).
4. Confirm with the V&V Lead that nothing reads ghost `UII` and `UIID` (D-073 (b)).
5. Radiation verification variants (A-50) with the V&V Lead.
6. Run the `nvfortran -gpu=nofma` build and the 8-thread sweep on the GPU test machine under the standing approval; ghost `UII`/`UIID`/`QR` are compared interior-only in restart checks (D-077 (6)).

**V&V Lead. Dry risk: low.**
1. Device-tier 2-ulp harness (D-070, D-072 (b); `vv/test-plan.md` device tier).
2. P3-B09 distance into the test plan and `vv-runs/phase3/phase3_cases.csv` after Role 3 confirms (A-61); include the observed front.
3. Gate hookup for the radiation kernels: apply the gate-rad patch to the gate script (adds the `rad` kernel family and its `KERNEL` result lines; `radiation/03-radiation-translation-notes.md`).
4. Rerun the Debug `int_1to2` case with a longer timeout and the 13 skipped four-rank Intel cases when the test machine is idle (`vv/baseline_status.md`).
5. Carry D-074 into the case inventory and the FR-003 list (`vv/case_inventory.md`, A-62 closed).
6. A-46 CSV rerun, A-47 baselines, A-26, A-39, A-42, A-49.
7. Flip-budget harness: PASS, FAIL and INCONCLUSIVE outcomes, 95% per run and 99% at phase exit, Level A and Level B, shortfall rule (D-077 (1); `solid/09`, test plan 5.10).
8. P-3 and N-5 on the host build, in the same series as patch 0010 (D-077 (2); test plan v0.4.18).
9. Review the six `gsfv_*` kernels on real-field device evidence and list them when accepted (D-077 (3)).
10. Ported review as reviewer of record: independent of the kernel authors, kernel owner notified (D-072 (c), D-077 (3)).

## GPU and build roles

**GPU Generator Engineer. Dry risk: low.**
1. DIVERGENCE_PART_1 nests: L0365, L0369, L0381, L0379 and the small nests L0371, L0377, L0380, L0386 (`amrex/loop-work-list.md` section 3; `amrex/s4d-divergence-part1-map.md`).
2. Generator change W1: `WLIST`/`NWL` indirection in the wall kernels (unblocks WP11, `amrex/stage1-gpu-spike-plan.md`).
3. Merge `s5-four` and `s5-wall` onto one branch and regenerate once (`solid/08-new-wall-kernel-review.md` G3); sidecar justification for scatters (G2) and the "every array argument is in the device list" structural test (G1).
4. L1375 EVALUATE_RAMP and the cell-list loop with `CYCLE` for the radiation plane body if wanted (`radiation/06-fr062-sweep-kernel-design.md` question 8).
5. P1 and V1 summation-order statements (`amrex/blocked-loop-families.md`).
6. Trace feature G1, status-301 predicate form G2, and ULP line format and NaN/Inf semantics G3 (`solid/09`; D-077 (1)).
7. Three additive `s5gen.py` features for SETTLING_VELOCITY: function callee, rank-2 table gathered through a `PRESSURE_ZONE` cell array, `CYCLE` of the kernel's own innermost loop.
8. Shared test-lock implementation: one full pass per engineer per hour under `flock`, FIFO, 30-minute hold limit, scoped lock for shared build directories, stale-lock rule; full pass before any commit that changes `s5gen.py` or a `kernel_lint`/`k2_ci_check` rule (D-077 (4)).

**GPU Wall Loops Engineer. Dry risk: medium.**
1. W1 to W5 from `solid/08-new-wall-kernel-review.md` section 7 (face-write refusal as a hard stop in table build and refresh, test additions, `wall_checks.py` builders).
2. O3 (L1402) uniqueness-assertion review with the Solid Phase Lead.
3. SP2 per-target gather, SP4 scratch sums and S2 guarded assertion (conditions from the Solid and Species leads, `solid/04-blocked-loop-signoff.md`).
4. L1358 COMPUTE_VISCOSITY wall part (table needed) and L1359, L0375 (`amrex/loop-work-list.md` section 3).
5. WP8/WP11 wall-table gather on the device with the Integration Lead.
6. Fix `rad_wall_qin_zero` and add `W_SURF_INDEX` to the device-address list so `ci_checks --strict` turns green; `ported` entries count for the map only until then (D-077 (3)).

**GPU Mesh Data Loops Engineer. Dry risk: low.**
1. Device run B of the four registered wall nests (`amrex/stage1-gpu-spike-plan.md` 7.13b; `solid/08` section 2).
2. Neighbour-mesh (T4) loops: L1363, L1364, L1366 (in progress), L1367, L1399 (`amrex/loop-work-list.md`, `amrex/o2-edge-exchange-scope.md`; D-055 excludes `NOM>0` branches).
3. "Ported" listing in `ported.toml` for the kernels with a device run on record, per D-072 (c), and the strict `ci_checks` pass.
4. Remaining O2 edge and exchange loops.

**GNU Build Chief. Dry risk: medium.**
1. Debug check of upstream patch UP-0008 (`upstream-patches/README.md`).
2. Add `csmag_32_fishpak` to `run_kernelcheck.sh` and fix the `FDS_K_MATCH` UBOUND bug (D-072 (g)).
3. Keep the Release and Debug reference binaries current for the merged tip (`vv/environment.md`); regenerate `vv-runs/refbin` entries on a new tip.
4. Build and test the pinned `AMReX_CUDA_FASTMATH=OFF` change on GNU after the Integration Lead applies it (D-068).
5. A-27 (GNU HYPRE pin and the second AMReX install).
6. Validate patch 0010 on GNU and run V&V tests N-5 and P-3 on the host build (D-077 (2)).

**Intel Build Chief. Dry risk: medium.**
1. `-prec-div` speed cost measurement and build notes (D-069; `vv/environment.md`).
2. Same-compiler baseline option for the output gate and `-prec-div` ifx kernelcheck references with `FDSREF_STEPS=2,3` (D-069, D-072 (g)).
3. Run the 13 skipped four-rank Intel baseline cases and the Debug `int_1to2` rerun with V&V (`vv/baseline_status.md`).
4. A-36 optional ifx compile check for K2 (NFR-045) on the generated kernels.
5. Intel-side check of the FASTMATH pin and the radiation gate family once applied.
6. Validate patch 0010 on oneAPI and run V&V tests N-5 and P-3 on the host build (D-077 (2)).

## Program

**Spec and Program Lead. Dry risk: low.**
1. Carry D-074 and D-075 into the spec (FR-003 list of the 23 FDS-only inputs; FR-006; FR-039 harness wording; HYPRE defaults) and close A-62.
2. Chase the Role 3 plan alignment (A-58), the A-61 confirmation, the Role 2 composite-solve date and the Role 1 domain-edge date.
3. Update the roadmap for dates given by Roles 1, 2 and 3.
4. Track the 90-kernel registration and the "ported" listings (D-072 (c)).

5. Carry D-076 and D-077 into the spec (flip-budget gate, NIC>1 guard ordering, reviewer of record, test lock) and the roadmap.