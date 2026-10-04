# FDS-AMR Verification & Validation Test Plan (DRAFT v0.4.21)

Owner: AMR V&V Lead · Status: **draft for review by the project owner**, 2026-09-25 · Nothing in this plan is approved yet. v0.3 brings the plan in line with `requirements.md` v0.3.1: v0.3, plus the bit-parity scope D-022/D-023 and FR-005, published while this revision was being written. It does not re-argue items the Spec Lead has decided (§9). The v0.2 text, including the full §9 rationale, is kept at `docs/vv/archive/test-plan-v0.2.md`.

**v0.4.21 (2026-10-04)** adds the second negative control N-5b next to N-5 in §5.4: a stretched-mesh variant of `race_test_1` that passes the input converter (which never reads grid stretching) and still has NIC 2 faces at level 0. Input from the Legacy Mapper (docs commit 8b1278e); copy in `vv-runs/negative/converter_pass_nic2/` with a README (sha256 recorded). N-5b is the control that exercises the guard on a stock `fds_amr` build; N-5 (ratio 5) is refused earlier by the converter on a stock build and reaches the guard only with a scratch bypass of the converter call. All counts are hand/script-derived and unrun. No other change. The text before this edit is archived at `docs/vv/archive/test-plan_v0.4.20.md`.

**v0.4.20 (2026-10-04)** aligns tests P-3 and N-5 (§5.4) with the assertion hook that ships in guard patch 0010 (D-076, D-077; docs commit 077a71e). The hook prints exactly one line, `NIC_CHECK level0: <n> walls checked`, to stderr; the longer proposed line format is replaced by it everywhere in the plan. P-3 now passes on that line present once with `n` equal to the hand count, no `ERROR(9001)` and the driver `level 0:` line present; N-5 additionally requires that no `NIC_CHECK` line is printed. The hook is owned by the Legacy Mapper. `vv-runs/negative/race_test_1_nic25/README.md` is updated to match. No other change. The text before this edit is archived at `docs/vv/archive/test-plan_v0.4.19.md`.

**v0.4.19 (2026-10-04)** records the V&V position on the device-tier flip-budget gate definition `docs/solid/09-flip-budget-gate.md` (commit 234163a) in §5.10, as a paragraph after the decision-flip rule: V1 one-sided exact (Clopper-Pearson) 95% bound gates with a three-way verdict, 99% reported as information; V2 denominator is libm-exposed distinct-input calls only, with the exposed predicate per decision site; V3 two levels (A gating, B gating once its class bound exists); V4 no waiver when the real population is short. It refines conditions 4 and 5 of the decision-flip rule (the point-estimate budget and the T1 1e-10 tolerance are replaced). **Pending Architect ratification.** The page has eleven sections, and V1 to V4 are its open questions in section 11. No other change. The text before this edit is archived at `docs/vv/archive/test-plan_v0.4.18.md`.

**v0.4.18 (2026-10-04)** adds the level-0 NIC=1 assertion test of D-076 (b) to §5.4 as P-3 (positive, `shunn3_4mesh_32` and `dec4_np4`) and N-5 (negative control, a copy of `race_test_1` with NIC 25). It rides on the guard patch `docs/upstream-patches/0010-main-amr-level0-nic-guard.patch` (D-065 Q1, drafted, not built), so it cannot run before that patch is in the AMR driver build, and all expected numbers are hand-derived. The negative control is judged on the message text and on the absence of the driver line `level 0: <n> box(es)`, not on the exit status (set-up errors end with a plain STOP, status 0). The input copy is `vv-runs/negative/race_test_1_nic25/`. The traceability row for FR-010 names P-3 and N-5. No other change. The text before this edit is archived at `docs/vv/archive/test-plan_v0.4.17.md`.

**v0.4.17 (2026-10-04)** marks the Intel baseline capture complete (V&V, `baseline_status.md` section "Intel / Intel MPI baseline"; `environment.md` v0.3.3). The Intel Release rebuild `impi_intel_firex-36975d7/impi_intel_rel` (rebuilt after an earlier environment reset) has G0, A-09 (Release and checked-debug), all Tier 1 cases and the A-24 pairs; the long four- and eight-rank cases and the `int_1to2` checked-debug run were run on a second machine (test machine). Compared with GNU: round-off or small differences for the non-chaotic cases, larger step-size-path differences for the chaotic and marginal cases, all FDS criteria that GNU meets are met. ULMAT references are recorded by pressure-solver library (HYPRE, MKL PARDISO), not by compiler. Still open on Intel: the restart-repeat calibration and the setup-only sweep. §3.1 and §3.2 (Intel status), §3.4 item 5 edited accordingly; nothing else changes. v0.4.16 is archived at `docs/vv/archive/test-plan_v0.4.16.md`.

**v0.4.16 (2026-10-04)** applies D-074 (requirements v0.4.36, A-62 closed): the 23 pressure-code-0 inputs (the 21 `soborot_*` and `Species/bound_test_1`, `_2`) are FDS-only in every mode, uniform and AMR. (1) §5.8: the "IN in uniform mode" wording and the open point about a no-pressure-solve driver mode are removed; the scope list moves the 23 from IN to OUT. FR-003 covers the verification set minus these 23 and the other FDS-only inputs. Denominators, the same in uniform and AMR mode: scope list 682 IN of 941 inputs (was 705 in uniform mode, 682 in AMR mode), OUT 49 (was 26); case inventory 152 IN of 187 rows (was 174 uniform, 152 AMR). G4: 847 of 861 firebot-listed cases with at most 8 ranks (was 870 of 884). (2) §5.12.5 wording follows. (3) `case_inventory.md`/`.csv` v0.6 and `scope_case_list.csv` v0.6 carry the change (update script `vv-runs/tools/update_inventory_v06.py`; `scope_filter.py` emits the same values). Not changed: the Tier 1/2 counts and run times in G2/G3/§6 (39 / 127 rows) still include the 22 inventory rows of these inputs as FDS-baseline reruns; for the AMReX code, Tier 1 is 37 runs and Tier 2 is 107 runs, to be refreshed at the next inventory render. The text before this edit is archived at `docs/vv/archive/test-plan_v0.4.15.md`.

**v0.4.15 (2026-10-04)** P3-B09 (1) is confirmed by Role 3: S = 2 stages, stencil half-width h = 2 cells (`GET_SCALAR_FACE_VALUE`, every limiter including MP5), frozen velocity, so d(n) = 4n fine cells inward from the patch outline after n steps. No other change. v0.4.14 is archived at `docs/vv/archive/test-plan_v0.4.14.md`.

**v0.4.14 (2026-10-04)** applies requirements v0.4.34 (A-61, A-62, D-071 (d)) and the GPU gate commit d6b421a. (1) P3-B09 (1): the propagation distance is d(n) = 4n fine cells (S = 2 stages, h = 2 cells), proposed in the spec and pending Role 3's confirmation against the driver (A-61); the report must print the observed front distance. (2) §5.8: the 23 FDS-only inputs are refused in AMR mode (D-071 (d)); they leave the AMR-mode denominators mechanically (`amr_scope`), and the answer on their uniform-mode scope is recorded (they stay IN: IR-001, FR-003, FR-006), with one open point on the driver. `case_inventory.md` and the CSVs go from v0.4 to v0.5. (3) §5.10: the gate side of the D-070 libm rule exists (`bitcmp.py --ulp`, `ULP` line format, `libm_check.py`, `kernel_categories.json`); the remaining gap is the driver side. v0.4.13 is archived at `docs/vv/archive/test-plan_v0.4.13.md`.

**v0.4.13 (2026-10-04)** applies the Architect rulings D-070 and D-071 and two Solid Phase Lead / Spec items. (1) D-071 (f): all PROVISIONAL Phase 3 thresholds are confirmed and the status is removed from §5.12: the P3-F02 corner limits are a gate (if exceeded, D-059 is reopened, the limit is not widened); P3-B09 (2) corner ≤ edge + 4 ulp; P3-B09 (1) uses the propagation distance from stages times stencil half-width per step (value supplied by the owner, TBD); the 0.5 discrimination rule on the tracer slice is a gate for P3-B02 and P3-B07 and report-only for P3-B08 (ratio 4); the discrimination part of P3-B05 is confirmed. (2) D-071 (f) and the Pressure Lead derivation for P3-R02: bound max\|div u − D − c\| ≤ 10·eps_rel·B + 20·eps_mach·U/dx_fine (derivation, not yet a measurement); working bound 1e-9·U/dx_fine accepted. (3) D-070: libm category in the GPU gate: kernels with libm transcendentals are compared at 2 ulp per value on the device, everything else bitwise; §5.10 and the G15 row; the host-side gate stays bitwise. (4) §5.10 records the V&V position on device decision flips in the 1-D wall solve (Solid Phase Lead question T4), pending Architect confirmation, and a note that V&V will be asked to sign the snapshot back-wall tolerances (D-065 Q2). v0.4.12 is archived at `docs/vv/archive/test-plan_v0.4.12.md`.

**v0.4.12 (2026-10-04)** applies the Spec Lead rulings of A-58 and requirements v0.4.33 to §5.12 and §5.11. (1) All PROPOSED Phase 3 thresholds of §5.12.6 become ACCEPTED-A-58 working thresholds, except the D-059 corner limits (P3-F02, P3-B09) and the use of the T3 discrimination rule on the tracer slice (P3-B02, P3-B07, P3-B08), which stay open until the Architect confirms. (2) FR-016(a) gate stays at step 1 with D and DS; Role 3 aligns. (3) New cases for D-063 (post-regrid composite projection; report-only in Phase 3, asserted once the composite MLMG is in the gate build; negative control with projection off) and D-062 (mass-weighted species transfer, negative control with linear Z): new §5.12.8, rows P3-R01 to R04, P3-S01, S02. (4) The exact-sum switch has its working name `EXACT_SUMS` in `&MISC` (logical, default F; T selects the exact fixed-point sum): §5.11, G11 and G14 rows, §5.3 and the §7 rows, and the Phase 3 rows now use it. (5) D-067 (pressure gauge: `sum(rho*V*(KRES-H)) = 0` per zone and component, applied to FFT cases too, `PRES` compared after removing a constant) replaces the open gauge question in §5.12.5. (6) The 23 FDS-only inputs stay IN with the sub-flag; no wording change. CSV: 41 rows, new input files `blob2d_mw_*`. v0.4.11 is archived at `docs/vv/archive/test-plan_v0.4.11.md`.

**v0.4.11 (2026-10-03)** adds one tooling bullet in §8 for `vv-runs/tools/patch_check.sh`, the standard bitwise behaviour-unchanged check for upstream patches (with `cmp_runs.py`, `pc_run_case.sh`). No other section changed. v0.4.10 is archived at `docs/vv/archive/test-plan_v0.4.10.md`.

**v0.4.10 (2026-10-03)** adds §5.12, the Phase 3 acceptance list for two-level runs (conservation, FR-016 mapped to concrete comparisons, moving-blob regrid tracking) and the D-057 `ns2d_16` check for the Pressure Backend Implementer, with the machine-readable `vv-runs/phase3/phase3_cases.csv` (35 rows) and the inputs, oracle and metrics under `vv-runs/phase3/`. G7 and G8 rows point to it. PROPOSED values are listed in §5.12.6 for the Spec Lead and Architect. No other section changed. v0.4.9 is archived at `docs/vv/archive/test-plan_v0.4.9.md`.

**v0.4.9 (2026-10-02)** (1) adds the GPU bitwise gate for the generated GPU kernels: gate row G15 in §5, new §5.10, tooling entry in §8 and files in §11. (2) applies the FR-005 (iii) test rule (requirements v0.4.31, D-053): new §5.11 gives every decomposition-invariance test (box split, rank count, thread count) its exact fixed-point sum switch setting, and the rows of G11, §5.3, G13, G14, §7 (FR-005, FR-062, FR-070, FR-076, NFR-010/011/012) carry it. MP5 `DIVG` GPU acceptance is tolerance-only until upstream patches UP-0001 and UP-0002 land (NFR-043). No other section changes. v0.4.8 is archived at `docs/vv/archive/test-plan_v0.4.8.md`.

**v0.4.9 addendum (2026-10-02).** §5.8 gets a one-line note on the 23 inputs that are FDS-only in AMR mode (no pressure solve, Chief Architect ruling); detail in `case_inventory.md` §9, column `amr_mode_status` in `case_inventory.csv` and `scope_case_list.csv`. No other change; the version number is unchanged. The text before this edit is archived at `docs/vv/archive/test-plan_v0.4.9_pre-amr-mode-status.md`.

**v0.4.8 (2026-09-26)** adds §5.9, the gas-gap face checks for the non-box level-0 ruling (N1-N3), with the vent list in `docs/vv/gap_face_vents.csv`. v0.4.7 is archived at `docs/vv/archive/test-plan_v0.4.7.md`.

**v0.4.7 (2026-09-26)** updates the tier sizes and runtimes to `case_inventory` v0.3: Tier 1 is 39 runs and about 17 min, measured, within the 30 min G2 budget; Tier 2 is 127 runs and at least 89 min. It also corrects the §5.8 wording on TUNNEL_PRECONDITIONER: the feature is OUT in AMR mode, but its inputs run with the keyword ignored. v0.4.6 is archived at `docs/vv/archive/test-plan_v0.4.6.md`.

**v0.4.6 (2026-09-26)** aligns the plan with requirements v0.4.23:
- D-040: the FR-014 order gate is "no worse than FDS on the same 2:1 case" (order ≥ FDS order − 0.1, L2 ≤ 1.1× FDS at each N); order ≥ 1.8 is reported, non-gating. Applied to G7 (c), G9, FR-039(c) and §7 rows FR-014/016/039. FR-039(c) also drops the stale maxorder-3 fallback (closed by D-032).
- D-039: the radiation stage is exempt from FR-005 (i) for box-split dependence only; rank independence at a fixed box layout and FR-005 (iv) run-to-run reproducibility still apply (G11, §7 FR-005, new FR-062 row).
- D-045 / ADR-004: `AMR_LEVEL` slice check (FR-078), `cmp`-identical Smokeview-format output across rank counts and redistributions (FR-076), AMR device values byte-identical across rank counts plus a MINLOC/MAXLOC tie case (FR-070), S-C Smokeview load (FR-072). Applied to G13 and §7.
- §5.8 note: Q4 answered (D-038 VARIABLE_THICKNESS IN, D-041 TUNNEL_PRECONDITIONER OUT in AMR mode with `tunnel_demo` kept as an MLMG regression/timing case, D-042 cylindrical DEFERRED in AMR mode). `scope_case_list.csv` counts are updated at the A-46 rerun.
- v0.4.5 is archived at `docs/vv/archive/test-plan_v0.4.5.md`.

**v0.4.5 (2026-09-25)** applies the Spec Lead rulings to §5.8 (requirements header read at v0.4.14). The section changed materially:
- GPU rule (D-035): CPU-inert GPU keywords no longer make a case OUT. The three `*_hypre` inputs without TUNNEL_PRECONDITIONER and the 9 CPU-runnable `GPU_Tests/` inputs are IN; `dancing_eddies_ulmat_hypre` is UNCLEAR (Q4).
- New counts: 705 IN / 175 DEFERRED / 26 OUT / 35 UNCLEAR. The `GPU_Tests/` part goes beyond D-035's four inputs, so the Spec Lead should acknowledge it.
- FR-040 R3 and the FR-034 leak-area summation are recorded as open points with the Chief Architect and Pressure Solver Lead.
- VARIABLE_THICKNESS, CYLINDRICAL and TUNNEL_PRECONDITIONER stay UNCLEAR pending Q4.
- New non-GEOM particle anchors proposed: FR-052 `VV/cascadempi_obst` (A-47), plus FR-050/051 replacements for `cascadempi` (`scope_alignment.md` §8).
- v0.4.4 is archived at `docs/vv/archive/test-plan_v0.4.4.md`.

**v0.4.4 (2026-09-25)** adds the FR-006 scope suite. FR-006 ("Verification-suite coverage of in-scope features") entered in requirements v0.4.10 under owner decision D-033 and V&V action A-41; `requirements.md` header read at v0.4.11. Changes:
- New §5.8 "Scope suite (requirement FR-006 from spec v0.4.10)". Its case list is `docs/vv/scope_case_list.csv` class IN. The section also sets pass criteria, how UNCLEAR cases get resolved, and that DEFERRED cases (GEOM/CC_IBM/HT3D) stay FDS-only baselines.
- New files: `docs/vv/scope_case_list.csv`, `docs/vv/scope_alignment.md`, `vv-runs/tools/scope_filter.py`.
- No other section changes. §6 tiers and `case_inventory.md/.csv` are unchanged; the status changes for inventory cases are listed in `scope_alignment.md` §6 for the next inventory refresh.
- v0.4.3 is archived at `docs/vv/archive/test-plan_v0.4.3_pre-scope.md`.

**v0.4.3 (2026-09-25)** aligns the plan with the requirements v0.4.3 items (the file header now reads v0.4.4, last read; v0.4.4 changes only the kernel-language items K2/IR-007/NFR-045, and none of them changes a test here). Changes:
- FR-016(a) covers the first ghost layer only (§5.2).
- FR-005 (ii): exact summation starts at the per-box sum. New FR-005 (v) test for setup areas and volumes (§5.3).
- FR-010 negative tests: ratio rejection, blocking-factor report, MLMG coarsening warning (§5.4).
- D-022 frozen-input kernel comparison excludes `STORE_SPECIES_FLUX` and CC_IBM cases until A-34 closes (§5.5).
- NFR-011 uses the A-35 derived race_test inputs (§5.6).
- S1/NFR-030 timing rules (§5.7).
- A-35 and the ratio-inventory cross-check (§6, §10).

v0.3 is archived at `docs/vv/archive/test-plan-v0.3.md`.

| Item | Value |
|---|---|
| Refactor base / acceptance baseline source | **FireX `36975d765f`** (`36975d765fcead401e14b094a04f910ac42eab8a`), local branch `AMReX`, this repository (read-only). Recorded with `git -C (repo root) rev-parse --short HEAD`. `git describe` = `FDS-6.11.1-1244-g36975d765f`. Working tree clean. |
| FireX reference builds | Build requests went to the **GNU and Intel Build Chiefs on 2026-09-25, GNU first** (spec §3.3). **GNU: delivered** (Release + Debug; §3.5). **Intel: built and baselined** (rebuild `impi_intel_firex-36975d7`, Release `-O2` + OpenMP, checked-debug and plain Debug; G0, A-09, Tier 1 and A-24 captured; `baseline_status.md`). The original request was blocked on the oneAPI reinstall (A-18, R-32). |
| Existing binaries | GNU/OpenMPI and Intel/Intel MPI builds of master **`ce1f659`** under `(local FDS master checkout)/Build/…`. **Smoke tests only** (toolchain, MPI, comparison tooling). GNU starts from a clean shell and passes the version print; Intel cannot start because oneAPI is missing (§3.1). Never used as references. |
| Requirements read | **Now: `docs/requirements.md` v0.4.4** (header; changelog in `docs/README.md`, last read). This plan implements the v0.4.3 items: D-028 (exact setup areas/volumes, FR-005 (v)), D-030 (legacy non-{2,4} inputs stay FDS-only), A-33 closed, A-34 narrowed, A-35, FR-010 blocking-factor and MLMG-coarsening checks, FR-016(a) first ghost layer, NFR-011 A-35 inputs, NFR-030 timing rules (1)-(6). requirements.md was **not** edited. *Previous read (v0.3):* `docs/requirements.md` **v0.3.1** (Spec & Program Lead, 2026-09-25; file still being edited, last read). v0.3.1 adds the bit-parity scope (D-022: T0/T1 against baseline only for explicit kernels on frozen input vs single-mesh FDS; whole runs T2), refluxing and tolerance-only multi-level comparisons (D-023), FR-005 (decomposition/rank/thread independence), A-24/A-25 and R-36. At the time of writing, `docs/README.md` (decision/action log, changelog) and `risks.md` still show v0.3 and do not yet define D-022, D-023, A-24, A-25 or R-36. The requirements text is authoritative. requirements.md was **not** edited. |
| Companion files | `environment.md`, `case_inventory.md/.csv`, `verification_case_survey.csv`: still v0.2 content (see §11); tools in `(local V&V run directory)/tools/`; A-19 inputs in `(local V&V run directory)/inputs/A-19/`; A-35 race_test remesh in `(local V&V run directory)/inputs/A-35/` (README); FM_Burner area check `docs/vv/fm_burner_area_check.md`. The Legacy Mapper inventory is at `docs/inventory/` (`mesh_ratio_cases*.csv`, `global_reductions.csv`, `upstream_issue_candidates.md`) |

Conventions: **T0/T1/T2/T3** are the tolerance classes of requirements.md §2.2. **eps_H** is the single-solve pressure tolerance of requirements.md §2.1a. **Tier 1/2/3** are case-set sizes (run frequency) and are unrelated to T-classes. Runtimes are model estimates (±3×) until measured. "Derived" inputs (`VV/…`) are V&V copies of committed inputs with a documented change. They are created under `vv-runs/`, never in the worktree. Item IDs "R-n"/"D-n" without a leading zero are this plan's §9 items. "D-0nn" are decision-log IDs in `docs/README.md`.

---

## 1. Purpose and scope

The refactor moves FDS `MESH` data onto AMReX grids, first in *uniform mode* (zero refinement levels) and then with static and dynamic refinement. V&V has to show, phase by phase (roadmap.md), that:
1. in uniform mode, results match baseline FDS to the class each requirement names (FR-001/002/003/037/038/061, IR-001/006), results do not depend on box split, rank or thread count beyond FR-005's limits, and the per-step projection-solver selection behaves as specified (FR-039);
2. with refinement, the code conserves (FR-020..024), reproduces FDS's existing coarse-to-fine coupling where that is the stated reference (FR-016), and is more accurate than the coarse grid and approaches the uniform-fine grid (FR-014, NFR-032);
3. restart, determinism, parallelism, output and performance requirements hold (FR-015/070-074/080, NFR-010..034).

This plan covers verification (code against analytic solutions, conservation, and baseline FDS). Validation against experiments is out of scope until compute beyond the 8-core development machine is available (charter Q6). The one validation input used, Heskestad `Qs=1_RI=10` for NFR-032, serves as a code-to-code AMR-vs-uniform-fine comparison, not as validation.

## 2. Tolerance classes as V&V will apply them

| Class | requirements.md v0.3 definition (short) | How V&V evaluates it | Tool |
|---|---|---|---|
| **T0 bitwise** | Byte-identical output, used only for (1) **kernel-level parity with single-mesh FDS**: output arrays of an explicit kernel (mass, velocity predictor, divergence) on frozen input (D-022), and (2) **self-comparisons** (restart vs continuous, IR-006, run-to-run reproducibility, FR-005, A-09b): byte-identical `CHID_devc/hrr/mass.csv`, timing/version columns excluded, same ranks, `OMP_NUM_THREADS=1`, on derived `&DUMP SIG_FIGS=17` copies (D-014) | (1) array-by-array comparison of kernel dumps; (2) CSV comparison, timing columns removed by name pattern plus a per-case list (e.g. DEVC `cpu` in `dancing_eddies_*`). **Never** applied to whole runs of the AMR code against FDS, or to multi-mesh FDS as a kernel reference. | kernel-dump comparator (planned); `compare_csv.py --class T0` |
| **T1 tight** | Same uses as T0: every kernel output field, or every CSV column at every output time: \|x−x_base\| ≤ 1e-10·max(‖x_base‖∞,1e-30) (proposed), on the same SIG_FIGS=17 copies | As defined; per column/field, no interpolation (output times must match) | `compare_csv.py --class T1` (`--t1-factor`) |
| **T2 verification** | FDS criterion e ≤ Tol_FDS always holds, and e ≤ (1+m)·e_base + 0.1·Tol_FDS with m = 5%, or max(5%, 3σ) for chaotic cases; plot-only/convergent series: observed order ≥ min(p_base, nominal) − 0.1 plus finest-grid L2 rule (provisional, D-013) | As defined. σ comes from the calibration ensemble (§3.4 step 5). | dataplot metrics re-implemented in numpy (no pandas/matplotlib on the development machine), `compare_csv.py` for series |
| **T3 physical** | Time-averaged plume ΔT and w 10%, HRR/MLR 2% (1% prescribed), layer height 5% or one fine cell, layer T 5%, probe mean 10%/RMS 25%; discrimination rule \|AMR − fine\| ≤ 0.5·\|coarse − fine\| always; calibration max(value, 2× GNU-vs-Intel spread) (provisional, D-013) | Window ≥ last 50% and ≥ 10 puffing periods; per-case scripts | T3 statistics script (planned) |
| **eps_H** (not a T-class) | `H` agreement for **same-discretisation, single solve on frozen input**: eps_H = max(1e-8, 2.4e-12·N²), N = largest cell count per direction, mean removed, relative L2 (D-012). Not valid for multi-step outputs (use T2/T3) or for coarse-fine comparisons against UGLMAT (FR-016(c)). | Used for FR-002 (composite single level vs baseline GLMAT), FR-037 (`amrex::FFT::Poisson` vs baseline FFT) and FR-039 (FFT vs MLMG). Needs a frozen-input solve harness in the code under test (§5.1). | frozen-input H comparator (planned) |

**Bit-parity scope (requirements §2.2, D-022/D-023), as V&V applies it:** whole runs of the AMReX code, single-level included, are judged at **T2** against FDS after the first pressure solve (Crayfishpak is replaced). Kernel parity is against **single-mesh** FDS only (committed, or a derived single-mesh copy per A-24). Multi-level runs are never compared bitwise with multi-mesh FDS. eps_H still applies to single solves on frozen input. Kernel-level parity needs frozen-input kernel dumps from both codes: from baseline FDS this is an instrumented build like A-09b (§3.6). **Until A-34 closes, the frozen-input kernel comparison excludes every case that sets `STORE_SPECIES_FLUX` or runs CC_IBM** (requirements §2.2 note, IR-007). The cases this removes are listed in §5.5. Whole-run T2 comparisons of those cases are unaffected.

Which requirement uses which class, and which class each case is judged at: §7 and the `tclass` column of `case_inventory.csv` (the inventory's `tclass` still reflects v0.2 and must be re-derived for D-022).

## 3. Baselines, reference results and the FireX build request

### 3.1 What exists today
- **Smoke-test binaries (ce1f659).** `fds_ompi_gnu_linux` (GNU Fortran 14.2.0, OpenMPI soname .40, OpenMP, HYPRE 3.0.0, SUNDIALS 7.5.0, no MKL) and `fds_impi_intel_linux` (ifx 2026.1.1, Intel MPI 2021.18, MKL static, no OpenMP). Details are in environment.md §2. They are used **only** to prove the toolchain, MPI launch and comparison tooling (G0 below).
- **Smoke status (corrected 2026-09-25):**
  - GNU: **PASS** (version print). From a clean shell, `ldd` resolves `libmpi_usempif08.so.40` to `/lib/x86_64-linux-gnu` (Debian `libopenmpi40`/`openmpi-bin` 5.0.7), and `/usr/bin/mpirun -np 1 fds_ompi_gnu_linux` prints the version banner (GCC 14.2.0, Open MPI 5.0.7, HYPRE 3.0.0, SUNDIALS 7.5.0). The earlier FAIL came from a non-clean shell environment and is withdrawn.
  - Intel: FAIL. `libimf.so` exists nowhere on disk and there is no `/opt/intel/oneapi/setvars.sh` or `mpiexec`. The oneAPI install is gone (the development machine rebooted Sep 23). Reinstalling oneAPI is a prerequisite for the Intel FireX build.

### 3.2 Where reference results must come from
**Reference (baseline) results for every acceptance test must come from FDS built from FireX `36975d765f`** (requirements.md §2.1: CMake, out of tree, HYPRE GPU options OFF, same compiler, MPI, flags, ranks and threads as the build under test; `ce1f659` results are not acceptance references). The `ce1f659` binaries can't serve even as a proxy: FireX changes radiation, init, pressure (GPU HYPRE paths, resource sets) and output code (`case_inventory.md` §1a), and it links a different HYPRE (§3.5). V&V does not build: A-08 records the build, it does not perform it.

**Build requests** went to the GNU Build Chief and the Intel Build Chief on **2026-09-25**, GNU first, with the spec in §3.3. Status:
- **GNU: delivered 2026-09-25** (§3.5). G0 and the Tier 1 capture are being run by another V&V worker, which records the results in `baseline_status.md`. This plan does not duplicate that record.
- **Intel: delivered as a rebuild and baselined (2026-10-04).** G0, A-09 on Release and checked-debug, Tier 1 and A-24 are captured (`baseline_status.md`, "Intel / Intel MPI baseline"). Acceptance items that name both toolchains now have both baselines; the remaining Intel gaps are the restart-repeat calibration for the T2 σ ensemble and the setup-only sweep.

### 3.3 Build request spec (to the GNU and Intel Build Chiefs; sent)

| # | Item | GNU build | Intel build |
|---|---|---|---|
| 1 | Source | this repository at `36975d765fcead401e14b094a04f910ac42eab8a` (branch `AMReX`). Check `git rev-parse HEAD` and an empty `git status --short` before and after. **The worktree is read-only**: nothing may be written into it. | same |
| 2 | Build system (primary, = baseline per requirements §2.1) | CMake ≥ 3.24, **out of tree**: `cmake -S (repo root) -B (local build directory)/firex-36975d7/ompi_gnu_rel -DCMAKE_BUILD_TYPE=Release`. Take the compiler wrappers from preset `ompi_gnu` (`CC=mpicc CXX=mpic++ FC=mpifort`). **Don't use the preset `binaryDir`** (`./Build/cmakeb/<preset>`): it is inside the worktree. | same, preset `impi_intel` wrappers (`mpiicx`/`mpiifx`), build dir `…/impi_intel_rel` |
| 3 | Compiler / MPI | Same as the existing binary: GNU Fortran 14.2.0 (Debian 14.2.0-19) + OpenMPI 5.0.x (soname .40), as NFR-020 requires | Same as the existing binary: ifx 2026.1.1 + Intel MPI 2021.18, with MKL visible to CMake (`find_package(MKL CONFIG)` → `WITH_MKL`) so the `*_pardiso` cases work as they do in the existing Intel binary |
| 4 | Options | `USE_HYPRE=ON`, `USE_HYPRE_NVIDIA/AMDGPU/INTELGPU=OFF`, `USE_SUNDIALS=ON`, `USE_OPENMP=ON` (CMake default). Record any deviation. | same. Note: the existing Intel makefile binary had **no** OpenMP; CMake enables it by default. Keep the CMake default and record it. |
| 5 | HYPRE | FireX CMake fetches hypre commit `63331f19` (requires ≥ 2.32.0). **Outcome (GNU):** fetched pin used (`USE_SYSTEM_HYPRE=OFF`). `63331f19c7` = `v2.32.0-24`, an unreleased commit between 2.32 and 2.33, **not 3.0.0**. The system v3.0.0 libs were rejected by `find_package(HYPRE 2.32.0)` (SameMajorVersion). See §3.5. | **Must use the same fetched commit `63331f19c7`**, so both toolchains share one HYPRE for FR-016/FR-038 references. Record its `git describe`. |
| 6 | SUNDIALS | v7.5.0 (fetched, or the existing libs); record. GNU used the system 7.5.0 libs. | same |
| 7 | HDF5 / VTK | Not available in CMake (no option). Leave out of the primary build; record "HDF5: no". **No HDF5 build exists yet**; it waits on D-019 (A-14, Q5). | same |
| 8 | Provenance (A-08) | `fds -V` must print a non-empty revision (CMake takes `git describe`; if building from a copy, pass `-DGIT_HASH=…`). Deliver: `CMakeCache.txt`, verbose compile lines (actual flags), `mpifort --version` / `mpifort -showme`, HYPRE/SUNDIALS versions, `ldd` output, sha256 of the binary, and the environment script that makes it run **on the development machine** | same, plus `mpiifx -show`, MKL version |
| 9 | Runnable on the development machine, reboot-proof | OpenMPI 5.0.7 runtime is present (Debian packages). Write down the environment setup (a script plus notes) so it survives the next reboot. | **Prerequisite:** reinstall oneAPI (ifx 2026.1.1, Intel MPI 2021.18, MKL) under `/opt/intel/oneapi` with `setvars.sh`; it is missing since the Sep 23 reboot. Write down the setup so it survives the next reboot. |
| 10 | Debug variant (requested, lower priority) | Same configuration, `CMAKE_BUILD_TYPE=Debug` (preset `ompi_gnu_db`), for A-09 NaN/bounds checks on the three AMR cases. **Delivered** (§3.5). | `impi_intel_db` |
| 11 | Optional: HDF5 variant for FR-074 / A-14 | Makefile target **`ompi_gnu_linux`** (`Build/ompi_gnu_linux/make_fds.sh`) with HDF5 1.14.5 (`build_thirdparty_libs.sh` builds it only if `$FIREMODELS/hdf5` is a clone at tag `hdf5_1.14.5`; otherwise it silently omits HDF5). Makefile builds run **inside** `Build/<target>`, so they must use an exported copy of the commit (e.g. `git archive 36975d765f`) outside the worktree. On hold pending D-019. | target **`impi_intel_linux`** (no OpenMP; `-O2 -ipo`) |
| 12 | Later (Phase 3): instrumented baseline for FR-016(a) | **A-09b**, placeholder spec in §3.6. | same |

Build names for reference: makefile targets `ompi_gnu_linux`, `impi_intel_linux` (+ `_db`, `_dv`, `impi_intel_linux_openmp`); CMake presets `ompi_gnu_rel`, `impi_intel_rel` (+ `_db`, `_dv`).

**V&V acceptance of a delivered build (G0):** `fds -V` shows `36975d765f`. `ns2d_16` runs on 1 rank and `obst_activation_default` on 4 ranks, each finishing normally in < 2 min. A setup-only sweep (`T_END=0` copies) over all 941 inputs records which ones stop at setup; that list becomes the IR-001 reference.

### 3.4 Reference capture procedure (after a build passes G0)
1. Order (requirements A-09, R-34): run `ns2d_16_int_1to2_refinement` and its derived `SOLVER='UGLMAT HYPRE'` copy on the FireX build **first**, on the Release and the Debug build of each toolchain, before FR-016 is treated as a gate. Then the two `emb` cases (A-09), Tier 1 and Tier 2. The int_1to2 pair is in Tier 1, so its GNU Release result comes with the GNU Tier 1 capture. The GNU Debug runs follow, and Intel follows A-18. **FR-016 does not gate anything until A-09 shows that baseline runs the case** (pass = runs to T_END with no NaN; analytic-error values recorded as baseline values). If the UGLMAT copy fails, FR-016 goes back to the Spec Lead and Pressure Lead for a replacement case (R-34).
2. Every run: `OMP_NUM_THREADS=1` unless the case is a threading case. Load < 2 and free RAM ≥ ranks×0.5 GB + 1 GB checked before start. `mpirun --bind-to none` (shared machine). Environment from the Build Chief's script (GNU: `(local build directory)/env/gnu_ompi_env.sh`). Watchdog per NFR-011 (R-10).
3. Every Tier 1 case runs twice: once as committed, and once as a SIG_FIGS=17 copy (D-014) for the T0/T1 self-comparisons that remain under D-022 (restart vs continuous, IR-006, reproducibility, A-09b acceptance). Chaotic cases (fire/LES) and the FR-002 cases also get the extra variants their class needs (GLMAT copy, alternate rank count).
4. Stored under `(local V&V run directory)/baseline/<toolchain>/<case>/` with a manifest: input sha256, binary sha256, `fds -V`, env, ranks/threads, launcher line, wall time, exit status, load/free-RAM at start. **The HYPRE version in the manifest comes from the build's `PROVENANCE.md`, never from the FDS run header** (the header prints a hard-coded "3.0.0"; §3.5). Outputs are read-only after capture.
5. T2 margin calibration ensemble (for chaotic cases): GNU vs Intel and firebot rank count vs one alternate rank count, giving σ for T2. The Intel baseline exists; the Intel restart-repeat calibration for chaotic cases is still to be run.
6. Measured wall times replace the model estimates, and cases are re-tiered.

### 3.5 FireX GNU build as delivered (GNU Build Chief, 2026-09-25)

| Item | Value |
|---|---|
| Release binary | `(local build directory)/firex-36975d7/ompi_gnu_rel/fds`, sha256 `9d3b5991447a70b8cd1d067e7db6fec8506479b32fd62961339f8d387e8050ec` (checked by V&V with `sha256sum`) |
| Revision | `FDS-6.11.1-1244-g36975d765f-AMReX` |
| Toolchain | GCC 14.2.0, Open MPI 5.0.7, OpenMP on, no MKL |
| Libraries | SUNDIALS 7.5.0 (system libs, `(local GNU third-party library tree)/libs/sundials/v7.5.0`, shared); HYPRE fetched commit `63331f19c7` (= `v2.32.0-24`, static); HDF5 not used |
| Fortran flags | `-O3 -cpp -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics -fopenmp` |
| Environment script | `(local build directory)/env/gnu_ompi_env.sh` (sets `HWLOC_LIBXML=0`) |
| Provenance | `(local build directory)/firex-36975d7/ompi_gnu_rel/provenance/PROVENANCE.md` (CMake cache, configure/verbose build logs, flags, ldd, third-party hashes) |
| Debug binary | `(local build directory)/firex-36975d7/ompi_gnu_db/fds`, built with `-O0 -ggdb -finit-real=snan -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` (for A-09 NaN/bounds checks) |
| HDF5 / VTK variant | none yet; waits on D-019 (A-14, Q5) |

**HYPRE note (risk to FR-016/FR-038 references; also NFR-021, R-30):**
- FireX's CMake pin `63331f19c7` is `v2.32.0-24`, an unreleased commit between 2.32 and 2.33. FireX requires HYPRE 2.x ≥ 2.32, so the development machine's system HYPRE 3.0.0 was rejected.
- FDS's run header still prints "Hypre library version: 3.0.0", because CMakeLists hard-codes that string. **Manifests take the HYPRE version from provenance, never from the header.**
- Consequences:
  - The FireX reference uses a different HYPRE from the old `ce1f659` smoke binary (a real 3.0.0). Nothing HYPRE-related may be compared across the two.
  - The **Intel build must use the same fetched commit**, or the GNU and Intel HYPRE references differ.
  - requirements.md §2.1/NFR-021 and risks R-30 describe `63331f19` as "labelled 3.0.0" and name the development machine's 3.0.0 as the HYPRE. The linked code is 2.32.0-24. Reported to the Spec Lead; not edited here.
  - If the AMR build later links a different HYPRE (e.g. one HYPRE for FDS and AMReX MLMG, R-30), FR-016 (UGLMAT-HYPRE reference) and FR-038 (T2 vs baseline HYPRE runs) are no longer like-for-like. Either keep `63331f19c7` in the build under test, or rebuild and re-capture every HYPRE baseline on the new version and record the difference (R-09). Trigger: any HYPRE change, or baseline `H` shifting beyond eps_H on a single frozen solve.

### 3.6 A-09b: instrumented baseline build (placeholder spec)
- **Purpose:** FR-016(a) needs fine-side ghost values of ρ, ρY_α, temperature and the ghost-averaged properties (`MU`, `KRES`, `D`/`DS`) on the first step of `ns2d_16_int_1to2_refinement`, on both sides of the coarse-fine interface. Baseline FDS doesn't write them. Coarse-fine (ratio-2) interface data also serves FR-016(b)/(c) diagnostics.
- **Build:** FireX `36975d765f` plus a V&V ghost-cell / coarse-fine dump patch that doesn't change arithmetic (writes only, after the ghost fill; no reordering of operations). Applied to a **scratch copy** exported outside the read-only worktree (e.g. `git archive 36975d765f` into `(local build directory)/scratch/…`), built out of tree with the same CMake configuration, compiler, flags and HYPRE commit as §3.5. **Never committed** to branch `AMReX` or anywhere else (D-005). The patch file and its sha256 are stored with the build provenance.
- **Owner:** AMR V&V Lead (with the GNU Build Chief), per requirements A-09b; due before the Phase 3 unit checks.
- **Acceptance:** the instrumented build is accepted only if it is **T0 against the uninstrumented build** of the same toolchain (§3.5) on **all Tier 1 CSVs** (SIG_FIGS=17 copies, same ranks, `OMP_NUM_THREADS=1`). Any T0 failure rejects the patch.
- **Open (detailed spec needed):**
  - which arrays and which code points (after which `MESH_EXCHANGE` code / ghost-averaging site in `docs/inventory/interface_averaging_sites.csv`);
  - the dump format and the matching dump in the AMR code;
  - whether H and the post-projection face velocities are dumped for FR-016(c)/FR-032 diagnostics.
- Possible extension (v0.3.1, D-022): the same scratch-copy mechanism could also dump the explicit-kernel outputs (mass, velocity predictor, divergence) on frozen input that the FR-001/FR-002 kernel-parity checks need from single-mesh FDS. To be decided with the spec.
- These need input from the **AMR Pressure Solver Lead** and the **AMReX Integration Lead**. Not yet requested.

## 4. Environment and execution rules
- Development machine: 8 cores (1 socket, AVX-512), 16 GB RAM with **no swap**, shared with other workloads (2.4–6.4 GB free observed). At most 8 ranks per case and at most 8 ranks in total across concurrent cases. No case above 12 GB (NFR-031).
- Python 3.13 + numpy only (no pandas, scipy or matplotlib). FDS's `Utilities/Python` plotting scripts can't run here as-is, so V&V re-implements the needed metrics in numpy (§8).
- Timing tests (NFR-030/033/034, and the NFR-032 wall-clock ratio) are valid only on an otherwise idle machine: load < 1, ≥ 8 GB free, no other V&V jobs. Otherwise the run is flagged invalid rather than failed (D-013).
- Smokeview: not installed. Reference version SMV 6.11.2 (FR-072).
- FireX runs use the Build Chief's environment script (GNU: `(local build directory)/env/gnu_ompi_env.sh`).

## 5. Test gates

| Gate | When | Content | Class / criterion |
|---|---|---|---|
| **G0** Toolchain smoke | Each new binary (incl. ce1f659 once runtimes exist) | `fds -V`; `ns2d_16` (1 rank); `obst_activation_default` (4 ranks); `compare_csv.py` self-test (identical, perturbed and re-timed synthetic CSVs) | Runs / exit 0. No T-class (ce1f659 output is never a reference). GNU FireX build: in progress by the capture worker (`baseline_status.md`). |
| **G1** Input compatibility (IR-001) | Phase 2, then each milestone | Setup-only (`T_END=0`) copies of all 941 inputs | Same set of setup stops and messages as baseline (message diff) |
| **G2** Tier 1 regression | Every merge to `AMReX` (≈ 17 min sequential, measured, budget ≤ 30 min) | 39 runs (§6) | Whole runs vs FDS: **T2** (single-mesh and same-resolution multi-mesh; FR-002 also T2 vs the single-mesh equivalent and vs multi-mesh FFT and GLMAT baselines) (D-022). Kernel parity (FR-001/FR-002): explicit kernels on frozen input at **T1, T0 per Q1**, vs single-mesh FDS (A-24 copies for multi-mesh anchors). Restart per FR-080 (G10). Determinism/reproducibility T0 (FR-005(iv)). |
| **G3** Tier 2 | Nightly while code changes; every phase exit | 127 runs (≈ 89 min sequential; a lower bound, because 118 rows are still model estimates and the 9 measured rows took about 4.6× their estimates) + optional 8 | T2 |
| **G4** Suite non-regression (FR-003) | Phase 2 exit, Phase 10 | All firebot-listed cases with ≤ 8 ranks, minus the 23 FDS-only inputs of D-074 (847 of 861) | T2. The 14 cases needing 9–64 ranks are exceptions recorded and approved by the project owner |
| **G5** Build-option equivalence (IR-006) | Phase 2 onward | FR-001 anchors with `USE_AMREX=OFF` | **T0** |
| **G6** FireX HYPRE options (FR-038) | Phase 2 | The 4 `Pressure_Solver/*_hypre.fds` at 1 and 4 ranks, `FDS_RANKS_PER_GPU` unset and =2 (the variable changes the matrix partition on CPU too: RS masters gather and solve, pres.f90:3411-3440) | **T2** vs baseline at the same setting (v0.3.1: whole runs, multi-mesh inputs, D-022); same HYPRE commit as baseline (§3.5) |
| **G7** Static two-level (FR-016) | Phase 3 (a–c), Phase 4 (d); **gating only after A-09 shows baseline runs the case (R-34)** | `ns2d_16_int_1to2_refinement` + `VV/…_uglmat` (`SOLVER='UGLMAT HYPRE'`) only (D-015). `ns2d_16_emb_1to2_refinement` is informational and non-gating ((a)-type comparison and FR-014 only). | (a) T0 on first-step fine-side ghosts, **first ghost layer only** (§5.2), vs the A-09b instrumented baseline (component-level only); multi-step and end-state comparisons vs `int_1to2` FDS are T2/T3, never bitwise (D-023); (b) coarse side not matched, FR-020/021 instead; (c) FR-032 normal-velocity mismatch at machine zero (target), `H` error against the exact `ns2d` solution no worse than UGLMAT-HYPRE's, observed order ≥ FDS's measured order − 0.1 and L2 ≤ 1.1× FDS's at each N on the same 2:1 `ns2d` case (u, w, `H`; N = 16/32/64; A-38 norm script, `H` mean removed, t ≈ 1; D-040), with order ≥ 1.8 reported but non-gating. **eps_H does not apply here** (not same-discretisation; replaces the v0.2 "10× solver tol" rule; A-16 closed); (d) T2 vs UGLMAT-HYPRE baseline. Phase 3 case list and thresholds: §5.12. |
| **G8** Conservation (FR-020/021/022/024) | Phase 3/4 | `species_conservation_1..4`, `energy_budget_*`, `simple_duct`, mass-balance cases with a refined patch | Requirement numbers (1e-12/step, 1e-10 cumulative; FR-022 with absolute floor) Phase 3 case list and thresholds: §5.12. |
| **G9** Accuracy (FR-014, NFR-032) | Phase 4, 10 | `ns2d_{8..64}` + static patch, `ns2d_16_*_refinement`; NFR-032 plume: Heskestad `Qs=1_RI=10` (D-008) on the trimmed 32×32×80 level 0, with the A-19 64×64×160 uniform-fine reference (§6) | FR-014: refined ≤ coarse, refined-region L2 ≤ 2.0× uniform-fine, observed order ≥ FDS's measured order − 0.1 and L2 ≤ 1.1× FDS's at each N on the same 2:1 `ns2d` case (u, w, `H`; N = 16/32/64; A-38 norm script, `H` mean removed, t ≈ 1; D-040), with order ≥ 1.8 reported but non-gating. NFR-032: T3 (Lf within one fine cell 0.057 m, HRR 1%, centreline ΔT and w 10%, discrimination rule) at ≤ 50% of uniform-fine wall time |
| **G10** Determinism / restart (FR-015, FR-080/081) | Phase 3, 9 | repeat runs, hierarchy dumps; `restart_test1a/b` + continuous copy, `device_restart_*` (+ `device_restart_base_case`), `restart_ulmat_*` | Hierarchy identical. FR-080 uniform mode: T0 proposed, **T1 minimum**. Only if baseline's own restart-vs-continuous difference misses T1 (measured at capture, A-08) does the criterion become "no worse than baseline's own restart-vs-continuous difference". AMR mode: T1 with an identical hierarchy after restart. FR-081: T2 |
| **G11** Parallel / threads / decomposition (NFR-010/011/012, FR-005) | Phase 2+ | Radiation stage per D-039: box-split dependence exempt, rank count at fixed box layout and run-to-run not exempt. Baseline FDS: {1,2,4,8} ranks ≤ mesh count (D-016), MULT-split copies for wider sweeps (T2 vs original). AMReX code: every anchor at 1, 2, 4 and 8 ranks, two `max_grid_size` values, two thread counts; a pressure-zone case (`zone_break_fast`) for FR-005(ii); FR-005(v) setup-area/volume cases (§5.3); 3 repeats per configuration; `race_test_1/4` (AMR mode: A-35 `_r4` copies; FDS baselines: committed originals; §5.6) | FR-005: reduction-free explicit stages byte-identical across box split/ranks/threads (i); zone integrals and gauge exact, with exact accumulation from the per-cell/per-box sum up, not only at the Allreduce (ii); `H` within eps_H across box split/ranks (iii); bitwise run-to-run at fixed ranks/threads/layout, `OMP_DYNAMIC=false` (iv); setup areas/volumes byte-identical across box split/ranks/threads (v). Whole runs T2 across ranks (uniform), T3 (AMR); NFR-011 watchdog **`EXACT_SUMS` (FR-005 (iii), §5.11):** every decomposition-invariance comparison runs with the `EXACT_SUMS=T` (exact fixed-point sum); the same configurations in default mode (`EXACT_SUMS=F`) are compared within eps_H only |
| **G12** Performance / memory (NFR-030/031/033/034) | Phase 2 (measure), 10 (meet) | `openmp_test64a`; one 8-rank multi-mesh anchor; NFR-033 `GPU_Tests/HYPRE_GPU_SCALING/test_8mesh_NOFRPG` 20-step copy | Timing ratios, median of 3, idle machine only; timing-run rules (1)-(6) of NFR-030 (§5.7) |
| **G13** Output (FR-070/071/072/074/076/078) | Phase 9 | header diff; Smokeview load; VTK (needs HDF5 build, D-019); AMR mode (D-045, ADR-004): S-A case at 1/2/4 ranks + one forced redistribution, `AMR_LEVEL` slice check, device values across rank counts, S-C Smokeview load | Exact header match (I); D for Smokeview/ParaView; AMR mode: Smokeview-format files `cmp`-identical across rank counts and redistributions (FR-076), `AMR_LEVEL` equals the level layout (FR-078), device values byte-identical (FR-070) Switch: ON for the rank-count and redistribution comparisons (§5.11) |
| **G14** Pressure single-solve checks (FR-002 `H`, FR-037, FR-039, FR-005(iii)) | Phase 2 (prototype P2 for FR-039), then Phase 4 | Frozen-input solves (§5.1) | **eps_H** (§2) on a single solve after gauge fixing; switch timing exact. Multi-step outputs of the same runs are judged T2, never eps_H `EXACT_SUMS`: not applicable to a single solve on frozen input; eps_H holds in both modes (§5.11) |
| **G15** GPU kernel bitwise gate (ADR-001 merge protocol item 4) | Quick tier: every commit that touches the generator, its markers or a marked loop. Full tier: each phase exit and each upstream FireX merge. Device tier: before a generated kernel is first accepted on the GPU and at each phase exit that changes a generated kernel | `vv-runs/gpu_gate/run_gpu_gate.sh --tier quick\|full` (§5.10); device tier is a documented command | **Bitwise** (T0 kernel parity: every output array and table, +0/−0 rule) on the host for every kernel; on the device bitwise too, except kernels in registry category `libm` (libm transcendentals), which are compared at **2 ulp per value** (D-070, §5.10); device decision flips in the 1-D wall solve follow the V&V position in §5.10 (pending Architect confirmation); golden signatures unchanged; every registered kernel has a bitwise test |

Negative tests (FR-004, FR-044, IR-002): one input per excluded feature, expecting `SETUP_STOP` and the named `ERROR(nnn)`. Written in Phase 3 as the supported feature set is defined. FR-010 refinement-ratio and blocking-factor negative tests: §5.4.

### 5.1 FR-039 (and FR-037/FR-002 `H`) proposed test

FR-039: the driver uses `amrex::FFT::Poisson` while the hierarchy is single-level and composite MLMG once a refined level exists. Both paths use the same pressure gauge, and MLMG uses `setMaxOrder(2)`. Owner: AMR Chief Architect / AMReX Integration Lead. Acceptance is part of the prototype P2 check (roadmap), then Phase 4. V&V supplies the cases, the frozen-input comparator and the switch check.

| Sub-test | Method | Cases (single mesh, one BC type per face, so FFT-able) | Criterion |
|---|---|---|---|
| FR-039(a) FFT vs MLMG, same input | Dump the Poisson RHS and BCs at steps 1 and ~N/2 of a uniform run. Solve once with `FFT::Poisson` and once with MLMG (maxorder 2, same gauge applied after the solve), then compare H (mean removed, relative L2) | `ns2d_16` (N=16), `ns2d_32` (32), `shunn3_32` (32; periodic), `csmag_32` (32; 3-D periodic), `dancing_eddies_1mesh` (N=300; inflow/open/wall faces + OBSTs, same all-cell operator in both). Optional: `shunn3_128` (128) | eps_H: 1e-8 for N ≤ 64; 3.9e-8 (N=128); 2.2e-7 (N=300) |
| FR-037 FFT fast path vs baseline | Same frozen input, `FFT::Poisson` vs baseline `pois.f90` FFT (FR-001 cases with one BC type per face) | as (a) | eps_H; plus D (FR-001 cases run with the FFT path selected) |
| FR-002 `H` vs baseline GLMAT | Same frozen input, composite single-level solve vs baseline GLMAT | `shunn3_4mesh_32` (M2a), then the other FR-002 anchors | eps_H |
| FR-039(b) switch | A derived `ns2d_32` copy with a static refined patch added at step n1 and removed at step n2 (needs the `&AMR` namelist, IR-003). Check the per-step solver log. | derived `VV/ns2d_32_switch` | Solver = FFT exactly on single-level steps and MLMG exactly on multi-level steps. Full-run outputs T2 vs baseline |
| FR-039(c) order at maxorder 2 (P2 item closed by D-032: MLMG stays at maxorder 2) | Refinement study at maxorder 2 on the 2:1 `ns2d` periodic patch | `ns2d_{16,32,64}` + patch | FR-014 gate (D-040): observed order ≥ FDS's measured order − 0.1 and L2 ≤ 1.1× FDS's at each N on the same 2:1 case (u, w, `H`; A-38 norm script, `H` mean removed, t ≈ 1), reference = A-38 FDS norms; order ≥ 1.8 reported, non-gating. A failure is a defect to fix at maxorder 2, not a trigger for maxorder 3 (D-032) |

Needs from the code under test: a frozen-input solve hook (RHS/BC dump + replay) and a per-step solver-selection log line. These go into the P2 interface. Not yet agreed with the Integration Lead.

### 5.2 FR-016(a) ghost-layer scope (v0.4.3)

The bitwise (T0) ghost-value check at a level jump covers **only the first fine-side ghost layer** of the transported and diffused fields (ρ, ρY_α, T, and the ghost-averaged `MU`, `KRES`, `D`/`DS`) on step 1, compared with the A-09b instrumented dump of `int_1to2`. The second ghost layer is not compared, either bitwise or against a separate tolerance.

**Rationale:**
- At a level jump FDS fills the second layer as a zero-gradient copy of the first (wall.f90:351-386; first order). The P1 fill-patch (`FillPatchTwoLevels` + `PCInterp`) is not required to reproduce that, so a bitwise check there would test an FDS artefact rather than a requirement.
- The second layer enters only through the wider advection stencils next to the interface, i.e. through the interface fluxes. Those are covered by two checks: the conservation/refluxing checks (G8, FR-020/021/024: round-off with refluxing, D-023) and the whole-run T2 comparison of `int_1to2` against its UGLMAT-HYPRE baseline (G7 (d)).
- A second-layer error large enough to matter shows up there. One too small to show up there isn't worth a dedicated check.
- Confirmed by the V&V Lead and recorded in requirements FR-016(a).

### 5.3 FR-005 (ii) and (v): exact sums from the per-box level (v0.4.3, D-028)

- **(ii) amended:** the exact (order-independent, fixed-point) accumulation of `DSUM`/`PSUM`/`USUM` and the pressure gauge starts at the per-cell/per-box sum, not only at the `MPI_Allreduce`. 68 of the 91 order-dependent rows in `inventory/global_reductions.csv` are per-process loops before the MPI call, and regridding reorders cells within a box.
  - Test (unchanged case, stronger check): on `zone_break_fast` (and on `zone_shape` for multiple zones), `D_PBAR_DT` and the gauge constant, dumped per step, are **byte-identical** across 1/2/4/8 ranks, two `max_grid_size` values and 1/2/4 threads. In AMR mode they also stay identical before and after a regrid that only reorders boxes. **`EXACT_SUMS=T`** (the byte-identity of zone sums holds only with the exact sum on, FR-005 (iii)); with `EXACT_SUMS=F` (default FDS order) the same configurations are compared within eps_H only (§5.11).
- **(v) new, setup-time areas and volumes:** the driver computes `FDS_AREA`/`AREA_ADJUST` and every setup-time area and volume sum in `global_reductions.csv` over the whole domain, with the exact sum.

| Sub-test | Input | Configurations | Criterion |
|---|---|---|---|
| (v-a) cylinder spanning boxes | Legacy Mapper MULT-voxel cylinder reproducer (R=H=0.1 m; variants a–e of upstream issue candidate #1: 1 mesh, 2 meshes/1 process, 2 meshes/2 processes, moved split, two cylinders sharing one MULT) as V&V copies | AMR code: level 0 split into 1, 2 and 4 boxes (two `max_grid_size`), 1/2/4/8 ranks, 1/2/4 threads | `AREA_ADJUST` per face class (top/side/bottom) and the resulting fuel MLR **byte-identical** across all configurations (T0), and equal to single-mesh FDS (see note) |
| (v-b) real case | `FM_15cm_Burner_CH4_2cm.fds` (MULT-voxel cylinder, top split over 4 meshes; §6 note and `fm_burner_area_check.md`) and its single-mesh copy `vv-runs/inputs/fm_burner_area/CH4_2cm_1mesh.fds` | as (v-a) | fuel MLR / (MASS_FLUX·πR²·tanh t) = 1.00000 in every configuration. FDS as committed gives 0.52083 on 1 process and 1.00000 with 1 mesh per process, so FDS multi-mesh is **not** a reference here. |
| (v-c) volumes | an anchor with pressure zones and OBSTs (`zone_shape`) plus `random_meshes` | as (v-a) | zone volumes and total gas volume byte-identical across configurations, equal to single-mesh FDS |

- Switch for (v-a) to (v-c): setup-time areas and volumes are exact in both modes (D-028); the tests run with **`EXACT_SUMS=T`** because the fuel MLR and zone integrals are observed after time steps that use per-step sums (§5.11).

- How the value is observed: FDS reports no adjusted area. Setup-only FDS never computes it (T_END=0 stops before `ADJUST_OBST_SHAPE_AREA`, init.f90:917; main.f90:238-248). The AMR code needs a setup diagnostic that prints `FDS_AREA`/`AREA_ADJUST` per shaped OBST/MULT and the zone volumes, **also in setup-only mode**. Until it exists, the fuel MLR of a few time steps is the observable, as in the FM_Burner check.
- **Note / clarification requested:** "reproduces single-mesh FDS" can hold bitwise only if FDS's sequential floating-point sum over one mesh happens to equal the exact sum rounded once. For equal voxel faces that is usually, but not always, the case. V&V applies **T0 across AMR configurations** and **T1 (1e-10 relative) against single-mesh FDS** unless the Spec Lead decides otherwise.

### 5.4 FR-010 negative and input-check tests (v0.4.3)

AMR mode accepts only refinement ratios 2 and 4, the same in every direction (D-030). Expected behaviour: a setup stop at input time with an error naming the mesh pair and the ratio. The blocking-factor check reports alongside it, and the MLMG coarsening check is a WARNING only.

| Test | Input | Expected |
|---|---|---|
| N-1 3:1 | `VV/ratio3_int.fds`: copy of `ns2d_16_int_1to2_refinement` with the fine patch `IJK=24,1,24` (3:1 against the 4×4 coarse meshes) | Rejected: error names the fine mesh and each coarse neighbour and "ratio 3" |
| N-2 direction-dependent | `VV/ratio_dirdep_int.fds`: same copy with the fine patch `IJK=16,1,32` (2:1 in x, 4:1 in z; each ratio is in {2,4} but they differ) | Rejected: error names the pair and both ratios |
| N-3 5:1 legacy | committed `Thread_Check/race_test_1.fds` | Rejected: meshes 1–3, 2–3, 3–4, 3–5, 3–6 at ratio 5 |
| N-4 legacy 3:1 and direction-dependent (secondary) | one each from the D-030 list, setup-only: `impinging_jet_*` (3:1) and `BST_FRS_6` (8:4) | Rejected with the named pairs |
| P-1 2-D positive control | committed `ns2d_16_int_1to2_refinement` | **Accepted**. Its single y cell (IJK(2)=1 on every mesh) must not count as a direction ratio of 1 (the inventory ignores y in 2-D inputs; inventory README §4 (f)). Otherwise the FR-016 gating case would be rejected. |
| P-2 A-35 remesh | `vv-runs/inputs/A-35/race_test_1_r4.fds` | Accepted for ratios (all five interfaces 4:1). With the default `blocking_factor` 8 it stops with the named blocking-factor error (level-0 domain 34×18×32; 34 and 18 not divisible by 8). With `amr.blocking_factor = 2 8` it starts. |
| W-1 MLMG coarsening warning | as P-2 | WARNING printed, coarsest MLMG grid **17×9×16** reported; the run continues |
| P-3 level-0 NIC=1 assertion (D-076 (b)) | (a) `Scalar_Analytical_Solution/shunn3_4mesh_32` (4 meshes, all `IJK=16,1,16`, cell size 0.0625×0.1×0.0625; Tier 2 in the inventory, its Tier 1 sibling `shunn3_4mesh_128` has the same layout) on 1 and 4 ranks; (b) `src/Source/driver/tests/cases/dec4_np4.fds` (16 meshes, all `IJK=8,1,8`, cell size 0.0625) on 4 ranks | `NIC_CHECK level0: <n> walls checked` once with `n` equal to the hand count, no `ERROR(9001)` (see below) |
| N-5 level-0 NIC 25 negative control (D-076 (a)) | `vv-runs/negative/race_test_1_nic25/amr_nic_guard_neg_race_test_1.fds` (copy of `Thread_Check/race_test_1`; only CHID and TITLE differ) on 1, 2 and 6 ranks | Set-up stops with `ERROR(9001)` before step 1 (see below) |
| N-5b level-0 NIC 2 negative control that passes the converter (D-076 (a)) | `vv-runs/negative/converter_pass_nic2/amr_nic_guard_neg_conv_pass.fds` (copy of the Legacy Mapper's `0010-neg-control-converter-pass.fds`: the `race_test_1` layout with mesh 3 `IJK=6,6,4`, mesh 5 with `TRNX_ID='XSTRETCH'` and three `&TRNX` lines, and `&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2 /`) on 1, 2 and 6 ranks | Set-up stops with `ERROR(9001)` before step 1 on a stock build with patch 0010 (see below) |

The derived N-1/N-2 copies are created under `vv-runs/inputs/FR-010/` in Phase 3, when the `&AMR` input exists (IR-003). They are Verification T (setup stop + message text).

**Level-0 NIC=1 assertion, P-3 and N-5 (v0.4.18, D-076; hook text v0.4.20, D-077).** In AMR mode level 0 has one cell size, so every mesh-to-mesh face of level 0 must have `EWC%NIC = 1`; a coarse/fine face is a refinement level, never a level-0 `EXTERNAL_WALL` with `NIC>1`. Evidence: `docs/amrex/level-jump-external-wall-check.md` (commit 410d370): level 0 refuses unequal cell sizes (`FdsAmr.cpp:30` at FDS-AMReX 991a759f79); FDS sets `NIC` at `main.f90:2184` (`main.f90:2111` at 36975d7); the argument that equal cell size and an aligned lattice give `NIC=1` is code reading of `init.f90:3147-3170` by the Legacy Mapper and has never been run. These two tests are what tests that reading. Placement: AMR set-up checks, run as setup-only copies (`T_END=0`, seconds, no time loop; the cheapest tier) at the G1/G2 cadence (every merge to `AMReX`); they need no baseline. **Prerequisite: the guard patch 0010 (D-065 Q1, `docs/upstream-patches/0010-main-amr-level0-nic-guard.patch`, extended with the `NIC_CHECK` hook, docs commit 077a71e; drafted, not built) must be applied to the AMR driver build. The guard and the hook are owned by the Legacy Mapper and ship in patch 0010. Until the patch is built in, neither test can pass, and neither has been run.**
- *P-3, what is asserted.* After set-up in AMR mode, the guard checks every `EXTERNAL_WALL` of every level-0 mesh: `NIC` is set only for walls that abut another mesh (`NOM>0`, `main.f90:2184`), so walls with `NOM>0` are the ones checked, and a wall with `NIC>1` aborts the set-up with `ERROR(9001)`. Domain-edge walls (`NOM=0`) carry no `NIC` and are not counted. The hook prints, once per run, to **stderr** (not stdout: unit 6 is reopened onto `CHID.out`), from rank 0: `NIC_CHECK level0: <n> walls checked`, where `n` is the number of `EXTERNAL_WALL` with `NOM>0` summed over ranks (`MPI_ALLREDUCE`). The line is not printed if the guard aborts. The count is what makes the test non-vacuous: a silent guard alone proves nothing about how many walls it looked at. A violation shows as the guard abort, so "0 violations" means no abort.
- *P-3, pass criterion.* In `stderr.txt`: the `NIC_CHECK level0: <n> walls checked` line is present **exactly once** with `n` equal to the hand count; no `ERROR(9001)`; and in stdout the driver line `level 0: <N> box(es)` is present (4 for `shunn3_4mesh_32`, 16 for `dec4_np4`). Hand count: `shunn3_4mesh_32` **256** (4 meshes × 64 walls); `dec4_np4` **512** (16 meshes × 32 walls). Periodic faces count, because the `PERIODIC` vents on the outer faces make each of those walls abut a mesh (`init.f90:3151-3160`); the y faces (one cell) are domain walls with `NOM=0` and are not counted. `n` must be the same at 1 and 4 ranks (the walls live on the owning rank). The counts are hand-derived and fixed at the first run; a different `n` means the derivation or the hook is wrong and the plan is corrected, not the check loosened. `race_test_1_r4.fds` is not a positive control (it still has a ratio-4 face; the guard aborts it too).
- *N-5, pass criterion (the guard works).* Set-up stops before step 1: stderr carries `ERROR(9001): AMR mode needs equal-resolution level-0 meshes. MESH <NM>: external wall cell <IW> (IOR=<IOR>) abuts MESH <NOM> with NIC=<NIC> (<count> such wall cells in this mesh). ...`, once per offending mesh per rank, then rank 0 prints `ERROR: FDS was improperly set-up - FDS stopped (CHID: amr_nic_guard_neg_race_test_1)` once; stderr has **no** `NIC_CHECK` line; stdout has **no** `level 0: <n> box(es)` line and no step 1. **The exit status is not a criterion:** set-up errors end with a plain `STOP` (`main.f90:2073`), status 0. The run must end inside 300 s on 1, 2 and 6 ranks (a timeout is a hang, FAIL). Hand-derived expected content (not run): the offending meshes are 1, 2, 5 and 6 with 24 wall cells each and mesh 4 with 36, all `NOM=3`, `NIC=25` (the ratio-5 faces of mesh 3, cell size 0.01 against 0.05; `NIC = 5×5`); mesh 3 is silent; with 1 rank only the first offending mesh prints. The text and numbers must match, otherwise the NIC reading or the guard is wrong. Before the guard exists the same input stops at `FdsAmr.cpp:30` with the cell-size message: that is a known abort, not a pass of N-5. On an `USE_AMREX=OFF` `fds` build the input runs normally (no guard), which checks the macro guard.
- *Which control runs on a stock build.* On a stock `fds_amr` the input converter refuses N-5 (ratio 5) before FDS set-up, so N-5 reaches the guard only on a scratch build in which the converter call is bypassed (`main.cpp:150`, the `prepare_amr_input` call replaced by `(void)conv`). **N-5b is the control that exercises the guard on a stock build with patch 0010.** Both are kept: N-5 has the ratio-5 NIC 25 face of the original legacy input, N-5b the NIC 2 face.
- *N-5b, why and what.* No plain `IJK`/`XB` layout can pass the converter and keep a level-0 NIC>1 face, so N-5b uses grid stretching, which the converter never reads: it sees six equal 0.05 meshes and exits 0 (6 level-0 meshes), while FDS set-up gives mesh 5 physical x cells 0.05, 0.025, 0.025, 0.10, 0.05, 0.05 and builds NIC=2 walls against it. *Pass criterion:* three `ERROR(9001)` lines then the `improperly set-up` stop; hand/script-derived content (unrun): mesh 3 first bad wall IW 50, IOR 2, abuts mesh 5, NIC 2, 4 wall cells; mesh 4 IW 338, IOR 2, mesh 5, NIC 2, 28 cells; mesh 5 IW 580, IOR -2, mesh 3, NIC 2, 32 cells (with 1 rank only mesh 3 prints); no `NIC_CHECK` line on stderr; no `level 0:` line and no step 1 on stdout; exit status 0 expected and not a criterion; 300 s limit. Without patch 0010 the driver stops at `FdsSetup.cpp:31` (TRN meshes unsupported) after set-up has passed: a known stop, not a pass. If the first run differs from these counts, the derivation is fixed, the check is not loosened.
- Input copy: `vv-runs/negative/race_test_1_nic25/` (input only and a README; the same file is used by the guard's own test, `docs/upstream-patches/inputs/0010-neg-control-race_test_1.fds`). N-5b: `vv-runs/negative/converter_pass_nic2/` (same arrangement, from `docs/upstream-patches/inputs/0010-neg-control-converter-pass.fds`).

### 5.5 D-022 frozen-input kernel comparison: cases excluded until A-34 closes

- **Rule:** exclude every case that sets `STORE_SPECIES_FLUX` or runs CC_IBM.
- **What triggers them:**
  - FDS sets `STORE_SPECIES_FLUX` internally when a DEVC/SLCF/BNDF asks for `ADVECTIVE|DIFFUSIVE|TOTAL MASS FLUX X|Y|Z` or `TOTAL MASS FLUX WALL` (read.f90:17011-17014).
  - Any `&GEOM` turns CC_IBM on for a time-stepping run (geom.f90:2698), as does `&MISC CC_IBM=T`.
- **Scan** (non-comment lines of the 180 committed inputs behind our 184 inventory rows; the 4 derived rows inherit from sources with none of these):

| Case | Tier / anchor | Trigger | Was it a kernel-comparison candidate? |
|---|---|---|---|
| `Sprinklers_and_Sprays/cascadempi.fds` | T2, anchor (FR-050/051/052) | `&GEOM` (lines 24-26) → CC_IBM | yes: FR-002 kernel vs a single-mesh/GLMAT copy (inventory `tclass`) |
| `Species/mass_flux_wall_yindex.fds` | T2 | `TOTAL MASS FLUX WALL` | yes (v0.2 `tclass` "T1 (T0 stretch)", single mesh) |
| `Species/mass_flux_wall_zindex.fds` | T2 | `TOTAL MASS FLUX WALL` | yes (as above) |
| `Species/mass_balance_gas_volume.fds` | T2 | `TOTAL MASS FLUX X/Y/Z` | yes (as above) |
| `Species/mass_balance_reac.fds` | T2 (8 meshes) | `TOTAL MASS FLUX WALL` | yes ("T1 first 100 steps"; would need an A-24 copy) |

- None of the FR-001 single-mesh anchors listed in §7 and none of the A-24 copies (`shunn3_4mesh_32`, `dancing_eddies_default`, `symmetry_test_mpi`) is affected. Their whole runs stay T2.
- Re-add these cases when A-34 closes: either `STORE_SPECIES_FLUX` output is shown not to read the J=0/K=0 rows, or those rows are restored.

### 5.6 NFR-011 race tests (A-35, D-030)

- **AMR mode** uses the A-35 derived inputs: `vv-runs/inputs/A-35/race_test_1_r4.fds` and `race_test_4_r4.fds`, with mesh 3 at `IJK=24,24,16` (dx 0.0125) so all five interfaces are 4:1.
  - Level 0: 34×18×32 at dx 0.05, `amr.blocking_factor = 2 8` (ratio-4 level 1) or `2`.
  - Pass: same metric as upstream (race_test_4 vs race_test_1 TMP/VEL at t=4.49–4.51 s, abs. error 0.01), no hang (watchdog).
  - The AMR reference values come from FDS runs of the `_r4` copies (1 and 4 threads), not from the originals.
- **FDS baselines** stay the committed `race_test_1.fds`/`race_test_4.fds`.
- The remesh raises the burner top and vent plane from 0.120 to 0.125 m, halves the burning-cell count (32 → 16; total fuel flow preserved by AREA_ADJUST) and moves the probe cell. So `_r4` results are never compared with the original baselines.
- Details: `vv-runs/inputs/A-35/README.md`.

### 5.7 S1 / NFR-030 timing-run rules (v0.4.3)

Every AMReX timing run (NFR-030, ADR-002 S1 / D-007) records and reports the following. V&V reads "acceptance requires all of" as: a run missing any item is invalid, not failed, like a run on a busy machine.
1. The pressure solver used (`amrex::FFT::Poisson` or MLMG), from the per-step solver log (G14).
2. Pressure-solve time separately from total time, from the TinyProfiler regions `FFT::Poisson::solve` and `MLMG::solve()`, inclusive and exclusive, per rank.
3. One level-0 (uniform) run forced onto MLMG, for a same-grid FFT-vs-MLMG comparison, so that AMR cost ratios don't mix the solver switch with the added cells.
4. Profiled builds link a separate AMReX install with `AMReX_TINY_PROFILE=ON`; the production-like install (`(local AMReX install)`) stays OFF. Headline NFR-030/D-007 ratios come from the clean install and the breakdown from the profiled one, unless the profiled build's overhead is measured to be under 2% on the same case (median of 3, idle machine).
5. The FR-010 MLMG coarsening WARNING, if any, and the coarsest MLMG grid size.
6. The rules apply to S1 once the AMReX driver exists, not to the A-26 FDS fine-reference runs.

### 5.8 Scope suite (requirement FR-006 from spec v0.4.10)

> **v0.4.6 note:** Q4 is answered. VARIABLE_THICKNESS cases are IN (D-038). The TUNNEL_PRECONDITIONER feature is OUT in AMR mode (D-041), but inputs that set it still run: the keyword is accepted, warned and ignored, so those cases stay in the suite and are compared T2 against baseline, and `tunnel_demo` stays as an MLMG regression and timing case. Cylindrical cases are DEFERRED in AMR mode and kept in uniform mode (D-042). The CSV counts below predate these rulings and are refreshed at the A-46 rerun.

> **FDS-only inputs in every mode (D-074; A-62 closed; case inventory v0.6):** 23 inputs (the 21 `soborot_*` and `Species/bound_test_1`, `_2`) have FDS pressure code 0 on non-periodic directions, so no pressure solve is set up. The Architect ruled that they are FDS-only in uniform mode and in AMR mode alike: `fds_amr` aborts on them with a message that names the FDS executable, and no no-pressure-solve driver mode is added. The driver sweep finds no verification input with a one-cell x or z direction, which the Architect confirms; the earlier one-cell Dirichlet description of this set no longer applies. FR-003 covers the verification set minus these 23 and the other FDS-only inputs (stretched TRN meshes, the four decided embedded/overlapping cases, VTK, missing files, DEFERRED GEOM/CC_IBM and HT3D). In the scope list their `cls` is OUT (reason column updated, `amr_scope` equals `cls`; `amr_mode_status` is kept), and in `case_inventory.csv` their 22 rows (21 inputs) carry `amr_scope` "FDS-only in every mode (D-074)". **Denominators (uniform mode and AMR mode are the same set):** scope list 682 IN of 941 inputs (175 DEFERRED, 49 OUT, 35 UNCLEAR); inventory 152 IN of 187 rows. They still run as FDS baselines only (T0 rerun against themselves) and are neither AMR nor uniform-mode failures. Of the 22 inventory rows, 2 are Tier 1 and 20 are Tier 2, so the AMReX-code Tier 1 is 37 runs and Tier 2 is 107 runs. Detail: `case_inventory.md` §9.

FR-006 (D-033, A-41) requires the AMReX code to run the FDS Verification inputs that cover in-scope features. The scope is:
- deferred, not dropped (FR-044): GEOM/CC_IBM and HT3D;
- must work with refinement: thin obstructions, HVAC, pressure ZONEs, level-set wildfire.

This section turns FR-006 into a run list and pass rules. The analysis behind it is `scope_alignment.md`.

**Case list.** The FR-006 case list is every row of `docs/vv/scope_case_list.csv` with `cls = IN`. That is 682 of the 941 `Verification/**/*.fds` inputs at FireX `36975d765f` (v0.4.5 rerun; 705 before D-074 moved the 23 FDS-only inputs from IN to OUT). The other classes:

| cls | Count | Meaning |
|---|---|---|
| IN | 682 | All features in scope and runnable on a CPU build here; same set in uniform and AMR mode |
| DEFERRED | 175 | Uses GEOM, CC_IBM (MISC or implied by GEOM) or HT3D |
| OUT | 49 | TRN stretching (IR-002/D-030); the four decided embedded/overlapping cases; VTK (D-019); missing external files (26 at v0.4.5); the 23 FDS-only pressure-code-0 inputs (D-074) |
| UNCLEAR | 35 | Needs a decision; reason in the `reason` column |

GPU rule (D-035):
- A case is OUT on GPU grounds only if it cannot run on a CPU build (a keyword in the script's `GPU_ONLY_KEYWORDS`; none at this pin).
- `&PRES HYPRE_DEVICE_RUN` is ignored on CPU builds (pres.f90:1177-1178).
- The ranks-per-GPU layout is the `FDS_RANKS_PER_GPU` job variable, not an input keyword.
- GPU runs of the AMReX code are verified under NFR-043/NFR-047, not here.
- The 9 `GPU_Tests/` inputs are IN at completes-only, run in G4 only (about 2,600 s each, inventory model). The four `_RSn` twins differ from their partners only in CHID.

Other rules for the list:
- `scope_case_list.csv` is generated. It is never edited by hand. Rerun `python3 vv-runs/tools/scope_filter.py` (static parse, about 10 s, no FDS run). The script re-checks its FireX source citations and reports any that no longer match.
- The list goes with the source commit. A new FireX base means a rerun and a diff of the `cls` column, recorded in §10.
- Derived inputs replace originals only where `needs_remesh` names a copy (at present `race_test_1`/`_4` → A-35 `_r4`, AMR mode only). The FDS baseline is always the committed original. The one exception is `stairwell`, whose baseline needs an MPI_PROCESS ≤ 7 derived copy (as `random_obstructions_fft` already has).

**Pass criteria for an IN case.**
1. *Completes:* normal exit at `T_END` with no NaN, trap, `SETUP_STOP` or watchdog stop (§3.4). Setup-only inputs (`T_END=0`) complete setup, and their G1 message diff is clean.
2. *Assigned gate:* the `tier`, `tolerance_class` and `gate` columns of the CSV.
   - Cases in the 184-row inventory keep their inventory tier and class (G2/G3 rules above).
   - Other IN cases are G4 (FR-003, T2 vs the FDS baseline), using the committed Verification-Guide metric where one exists. Ten must-work multi-mesh cases are promoted to T2 in G3.
   - Cases FDS itself does not run (absent from or commented out of `FDS_Cases.sh`) are completes-only.
   - G7 (FR-016), G11 (FR-005) and G13 (output) apply where §5 already assigns them.
   - `STORE_SPECIES_FLUX` cases take no part in the D-022 kernel comparison until A-34 closes (§5.5); their whole-run T2 is unaffected.
3. *Baseline failures:* an IN case that fails on the FireX baseline is recorded, not counted against the AMReX code, and reported to the project owner.
4. *Exceptions:* any other exception (a case dropped, a class loosened) needs the project owner's approval and is recorded in §10.

Milestones:
- Phase 2: every IN case in uniform mode.
- From Phase 3: the refined variants (below).
- Phase 10: complete.

**Refined variants (A-41).** None of the four must-work features has an IN Verification case at a resolution change, and level set has no multi-mesh case (`scope_alignment.md` §4). V&V therefore adds derived refined variants under `vv-runs/inputs/A-41/`, A-35-style: the original stays the FDS baseline, and each variant gets its own FDS baseline. First set:
- thin duct walls across a 2:1 interface (`duct_flow`);
- HVAC vents at a resolution change (`fan_test`);
- a zone breach across a level jump (`zone_shape`, `zone_break_fast`);
- a level-set front crossing a 2:1 interface (`LS4_ember_yield` or `level_set_fuel_model_1`).

These are T2 vs their own FDS baseline and gate from Phase 3.

**Resolving UNCLEAR cases.** Each UNCLEAR row names one reason group. Each group is closed by a single recorded decision (decision-log ID), after which the script rule is updated and rerun, and the cases move to IN, OUT or DEFERRED. Until then they run in the baseline only and do not gate.

| Group | Cases | Decided by |
|---|---|---|
| VARIABLE_THICKNESS vs deferred HT3D | 12 | Project owner / Spec Lead |
| CYLINDRICAL | 11 | Q4 |
| TUNNEL_PRECONDITIONER (includes `dancing_eddies_ulmat_hypre`) | 7 | Q4 |
| ZONE leakage over several meshes | 2 | FR-034 leak-area summation (open with the Chief Architect and Pressure Solver Lead) |
| Thin OBST at a static coarse/fine interface | 1 | FR-040 R3 ruling (open with the Chief Architect and Pressure Solver Lead) |
| Per-mesh CSVF restart files | 1 | Integration Lead |
| Disjoint multi-resolution meshes | 1 | A-35-style remesh |

The three Q4 groups (VARIABLE_THICKNESS, CYLINDRICAL, TUNNEL_PRECONDITIONER; 30 cases) stay UNCLEAR until the project owner answers Q4.

The 40 IN cases with a non-box mesh union (`nonbox_domain`) are IN on the assumption that inter-mesh gaps become solid level-0 cells. If the Integration Lead does not confirm this, they move to UNCLEAR.

**Particle anchors without GEOM (A-47, proposed).** `cascadempi` is DEFERRED, but FR-050/051/052 stay in scope. Proposed anchors (details, metrics and file:line evidence in `scope_alignment.md` §8):
- **FR-052:** derived `VV/cascadempi_obst` (the GEOM half of `cascadempi` removed). It has 3 meshes and 28,672 cells, and box tops are split by the x = 7 mesh face. Its 2:1 variant puts the level jump through the box tops. Metric: committed dataplot rows `OBST-S/OBST-M-water mass`, T2. Secondary: `VV/geom_sprk_mass_obst` (R-7).
- **FR-050:** `bucket_test_1` (4 meshes, spray from the 4-mesh corner), then `VV/cascadempi_obst`.
- **FR-051:** `energy_budget_particles` (unchanged), plus `VV/cascadempi_obst` (evaporation mass balance across meshes).

Each derived copy gets its own FDS baseline and gates from Phase 7.

**DEFERRED cases** stay FDS-only baselines (captured as usual, used for G1 and as FR-004/FR-044 negative-test sources) until the scope changes. With refinement they must stop at setup (FR-004). A-41 notes that `ht3d_energy_conservation_4` may stay in the uniform-mode FR-003 regression; that is outside FR-006. When a deferred feature returns to scope, its cases are reclassified by rerunning the script with the rule removed.

### 5.9 Non-box level 0: gas-gap faces (ruling `adr/drafts/ruling-nonbox-level0.md`, N1-N3)

In AMR mode, level 0 covers the bounding box of the level-0 meshes, and cells outside every mesh are static solid gap cells. A gas-gap face must behave exactly like the FDS exterior wall it replaces: default surface, plus any `&VENT` FDS snaps onto that face, including `OPEN` (N2). The list of affected vents is `docs/vv/gap_face_vents.csv` (per case: `gap_face_summary.csv`; method: `gap_face_vents.md`; tool: `vv-runs/tools/gap_face_vents.py`).
- Scope: 40 IN cases have a non-box level 0, all with one cell size, so every mesh is level 0. 21 of them have at least one vent on a gas-gap face: 20 `OPEN` vent lines (31 per-mesh vents) and 23 other lines (49 per-mesh vents; 5 MIRROR, the rest HOT, COLD, SLIP, GROUND, LS GRASS and velocity surfaces). The other 19 only have default-surface gap walls (15 of them set `&SURF DEFAULT=T`).
- Cases with `OPEN` on a gas-gap face (11): `Controls/bi_dir`, `Flowfields/velocity_bc_test`, `HVAC/qfan_multi`, `Heat_Transfer/back_wall_test`, `Pressure_Solver/hallways`, `Pyrolysis/shrink_swell`, `Radiation/hot_spheres`, `Restart/device_restart_a`, `device_restart_b`, `device_restart_base_case`, `WUI/LS4_ember_ignition`.
- Check G-1 (setup, T0): for every row of `gap_face_vents.csv`, AMR-mode setup produces a wall record on the same face with the same snapped extent and SURF_ID as FDS's `.smv` vent listing, and the same default surface on the remaining gap-face cells. This extends the FR-005 (v) setup comparison.
- Check G-2 (pressure single solve, eps_H): the N3 masked-MLMG solve against FDS `SOLVER='UGLMAT'` on derived copies of `hallways` (`OPEN` covering a whole gap face) and `device_restart_a` (`OPEN` given off-face and snapped by FDS's half-cell rule, extent also rounded), which covers the N3 [VERIFY] on where UGLMAT places the `OPEN` value. `simple_duct` stays in G-2 as the sealed-gap case only: its vents sit on interior planes of a mesh, not on gap faces, so it does not test `OPEN` placement.
- Check G-3 (whole runs, T2): the 21 cases with gap-face vents run in AMR uniform mode against baseline under §5.8. Partly obstructed `OPEN` gap faces (`velocity_bc_test`, `qfan_multi`, `back_wall_test`) act as `OPEN` only on their gas part.
- `stairwell` has no vent on a gas-gap face (its `OPEN` and extract vents lie on the bounding-box boundary), so the N6 fallback condition "no gap-face `OPEN` vent" holds for it.

### 5.10 GPU kernel bitwise gate (G15)

**Purpose.** Every GPU kernel generated by `s5gen.py` (and its builders) is checked against the verbatim upstream Fortran loop on the host, so that regression catches changed bits, an unreviewed generator change, and a kernel without a test. This implements the "bitwise gate" of the ADR-001 section "Upstream FireX merge protocol for ported kernels" (item 4): a merge that touches a ported routine and skips its bitwise test is not accepted. Entry point and details: `vv-runs/gpu_gate/run_gpu_gate.sh`, `vv-runs/gpu_gate/README.md`. The gate wraps the existing scripts of the generator worktree (`amrex/s4_mass/s5_gen`), works on a scratch copy and never edits the worktree.

**Tiers and where they run.**

| Tier | Content | Runs |
|---|---|---|
| quick | host only, no GPU, at most 3 cores under `nice`, minutes: golden-signature check and regeneration reproducibility; text goldens of the one-face and whole-field `GET_SCALAR_FACE_VALUE` kernels; bitwise drivers (round 1, rounds 2 to 7, one-face callee, whole-field kernels, scratch-pointer wall nests) for flag sets `O0` and `O2omp`, callee switch `dpd`; kernel-to-test coverage check | every commit that touches the generator, its sidecar markers, the generated files or a marked loop |
| full | quick plus all six flag sets (`O0 O2 O0omp O2omp O2omp_off O2omp_dpd`), both callee switches (`dpd`, `bind`), OpenMP sets at 4 and 8 threads (the one-face and whole-field drivers run 1, 4 and 8 threads inside the program), the optional wall tests of branch `s5-wall` once merged (an absent script is reported SKIPPED-absent, never PASS); mutation checks (a mutant that is not caught fails the gate) and the generator negative checks are options (`--mutants`, `--negative`), used at phase exit; the full tier queues on the generator lock and holds it for the whole run | each phase exit and each upstream FireX merge; for a merge, `run_gpu_gate.sh --merge OLD NEW` first runs the hunk classifier of ADR-001 item 3 (`port_merge_check.py`; exit 1 = a surface-changing hunk touches a marked kernel, reviewed by hand before the full tier) |
| device | documented command for the GPU machine (nvfortran, `-gpu=nofma`, gfortran reference, hash comparison of every array; logs checked by `check_device_logs.py`); never part of the host gate | before first acceptance of a kernel on the GPU; each phase exit that changes a generated kernel |

**Pass criteria.** (1) Every kernel's output arrays and tables are bit-for-bit equal to the reference in every configuration that ran; a sign-of-zero-only difference fails; a test case that changes nothing (vacuous) fails. (2) `s5gen.py` output equals `markers/golden_signatures.json` (otherwise exit 3, a review event) and the regenerated files equal the committed `generated/` files. (3) Every kernel registered in the sidecar files has at least one bitwise test and produced a result in every configuration that ran; a kernel without a test fails the gate loudly. (4) Mutation checks: every mutant is caught. (5) Tolerance is used only where a documented ruling names the kernel: stage-1 ruling (c) (sums are not bitwise against the CPU path; none of the generated reductions needs it, because `CHECK_STABILITY` and `CHECK_DIVERGENCE` extrema with locations are order independent and reproduce the serial tie rule) and D-070 (category `libm`, 2 ulp per value on the device tier only, see below; the host gate stays bitwise).

**Libm category (D-070, requirements v0.4.33 / README decision log).**
- Kernels whose arithmetic is only `+ - * /` and `sqrt` stay **bitwise** under the +0/−0 rule (criterion (1) above). Kernels that call libm transcendentals (`**` with a non-integer exponent, `exp`, `log`, trigonometric functions) are category `libm` in the kernel registry. On the **device tier** a `libm` kernel is compared at **2 ulp per value** against the host (the measured difference is 1 ulp); everything else on the device is compared bitwise. The measured value is reported with the verdict.
- The **host-side gate** (gfortran reference against the generated kernel, quick and full tiers) stays **bitwise for every kernel, `libm` ones included**, because reference and generated code share the same libm on the host.
- First registered case: `cfl_wall_max` / `UVWMAX` with `(ABS(Q)/RHO)**(1/3)`: the device differs from the host by 1 ulp for about 13% of the arguments (3 of 13 scenarios differ in the last bit). Dt-coupled quantities (`UVWMAX` feeds dt) use the run-level tolerance already used for cross-compiler runs, not the per-value rule.
- There is no device `pow` of our own: the host libraries (gfortran libm, Intel libm) already differ from each other.
- The K2 CI check flags every libm call, so none enters unlisted; a kernel that calls a transcendental and is not in category `libm` fails the gate.
- **Gate tooling (commit d6b421a in `vv-runs/gpu_gate`; gate side done, driver side open).** Done on the gate side: `bitcmp.py --ulp A B [--width 64|32] [--limit 2] [--emit KERNEL SIZE [--out NAME]]` compares two raw dumps value by value (ulp distance from the ordered-integer keys; +0 and −0 are distance 0; NaN must be bit-identical in position and payload; infinities must match; self-tested); `check_device_logs.py` judges a `libm` kernel only by its `ULP` lines, in the format `ULP <kernel> <NXxNYxNZ> maxulp N nvals M nspecial S [out NAME]` (PASS if maxulp ≤ 2, nspecial = 0, nvals > 0; a malformed line fails the log); `libm_check.py` (run by `gate.py` in the coverage, quick and full tiers) fails a kernel that calls a transcendental and is not in category `libm` (`SQRT` is excluded, it is correctly rounded; `**` with an integer literal or an INTEGER-declared exponent is not libm); `kernel_categories.json` registers `cfl_wall_max` as the only `libm` kernel (the sidecars have no `category` key; a sidecar `category = "libm"` is honoured if the generator owner adds it). Non-libm kernels stay bitwise and a `ULP` line of such a kernel is not evidence. **Remaining gap = driver side:** no device driver prints a `ULP` line yet (`s5_dev_hash.F90` and the round 4 to 7 driver print only 64-bit hashes), so `cfl_wall_max` is MISSING, which counts as a failure, on the device tier until the GPU engineers add the line (either the driver reads a host dump and applies the same ordered-integer rule, or both sides write raw dumps that `bitcmp.py --ulp` compares). Open point for the owner of D-070: all four outputs of `cfl_wall_max` (`UVWMAX` and the location) feed dt and use the run-level tolerance, so the per-value rule decides nothing for that kernel today; the `ULP` line must still exist so the deviation is on record. The generator is not edited by V&V.

**Device-tier rule for the 1-D wall solve: decision flips (V&V position, pending Architect confirmation).** Source: the Solid Phase Lead's test design (`docs/solid/06-solve-port-test-design.md`, section 10, question T4). On a device a 1-ulp libm difference can flip `NINT(TMP_RATIO)` (the sub-step count), a remesh threshold (the node count) or the Newton exit test (one more iteration) in the 1-D wall solve, so the value difference is then physical, not 2 ulp. V&V accepts this as a **device-tier rule** under these conditions:
1. each record-call is compared from identical frozen inputs; nothing accumulates across time steps;
2. the **decisions** (sub-step count, node count `NWP`, layer list, Newton iteration count) are compared **exactly first**;
3. where the decisions agree, the values are compared under the 2 ulp rule above;
4. at most **1 decision flip per 10**4 record-calls**; every flip is reported with the inputs that caused it (the statistic is refined by the position below);
5. a call with a flip is compared against a documented physical tolerance instead of 2 ulp: proposed **T1 class on temperature and heat flux** (|x − x_host| ≤ 1e-10·max(‖x_host‖∞, 1e-30)); the exact tolerance is to be confirmed with the Solid Phase Lead (replaced by the two levels of the position below);
6. the flip budget is a **gate** and is **measured, not assumed**; a flip rate above the budget goes back to the Solid Phase Lead and the Architect (it is not absorbed by widening a tolerance).
Status: V&V position, **pending Architect confirmation**. On the host there are no flips (shared libm), so the host gate is unchanged and bitwise.

**Flip-budget gate definition: V&V position on `docs/solid/09-flip-budget-gate.md` (commit 234163a; v0.4.19; V1 to V4 of its section 11; pending Architect ratification).** The page defines the flip classes D1 to D5, the populations (P1 captured real calls, P2 random calls, P3 adversarial calls) and the counterfactual host check. V&V takes these positions on its open questions; they refine conditions 4 and 5 above and are not in force until the Architect ratifies them.
- **V1, statistic and verdict.** The gate is the one-sided exact (Clopper-Pearson) bound at 95%, per population, with budget b = 1e-4. Three-way verdict: **PASS** if the upper bound ≤ 1e-4; **FAIL** if the lower bound > 1e-4; otherwise **INCONCLUSIVE**, which does not sign the device tier. The 99% bound is computed and reported as information, not gated. Minimum N to show 1e-4 with zero flips: **29,956 at 95%** (46,050 at 99%); with k flips the minimum grows (k = 1: 47,437 at 95%). A point-estimate reading (one flip per 10**4) is not used: N = 10**4 cannot demonstrate the budget even with zero flips.
- **V2, denominator.** N counts **libm-exposed calls with distinct inputs only** (distinct frozen images, by hash). The exposed predicate is **written per decision site** (for each of D1 to D5: which libm call or libm-fed value lies on the chain to that decision), and the call-level tag of the page (section 4.3: PYR, GEOM, STRETCH, HTC) is derived from the site predicates, not the other way round. The gate counts **exposed and non-exposed calls separately**: exposed calls enter the rate; non-exposed calls must have equal decisions and bitwise equal outputs, any flip or difference there is a failure, not budgeted. Non-exposed calls never enter the denominator.
- **V3, two levels, never one tolerance.** **Level A (gating, every flip):** the device result equals the host **forced to the same decisions** within 2 ulp per value (D-070). **Level B (gating once its class bound exists):** the physical size of each flip, measured on the host as the difference between the free host and the forced-flip host, must not exceed the class bound from the host forced-flip survey; **a flip larger than its class bound FAILs regardless of the count.** For classes without a solver error estimator (sub-step count, `REMESH_RATIO`, cell counts, delamination) an **empirical class bound is accepted if the survey states the number of forced flips and the largest effect seen**; a class without a signed bound has no PASS (its flips are recorded). Classes with a bound derived from a stopping tolerance (temperature extraction 1e-4 K; oxygen Newton 1e-6 kg/(m² s); layer-removal thresholds) use that bound. **Captured real calls (P1) and random calls (P2) are gated and reported separately and never pooled.** P3 flips are reported and never counted in a rate.
- **V4, shortfall of real calls.** If the real cases cannot reach the minimum N, the real-population verdict **stays INCONCLUSIVE; there is no waiver.** The response is to add real cases with distinct inputs (more supported inputs that have a 1-D wall). The random population is kept separate and does not stand in for it. V&V reports the **achieved upper bound** (and N and k) of the real population for the Architect to decide; V&V does not decide that the gate rests on P2.
- Not decided here: the factor applied to the largest surveyed effect when an empirical class bound is set (the page proposes twice the largest observed, signed by V&V and the Solid Phase Lead), and the survey size per class. The device tier cannot be signed on this rule before the harness, the forced-flip survey and the exposure predicates exist; nothing has been run.

**Recorded for later (D-065 Q2).** The back-wall heat transfer coefficient in the thick-wall and thin-wall passes uses a snapshot of the other side in the AMR route (FDS-only mode untouched). V&V will be asked to sign the tolerances for this snapshot mode on `back_wall_test` and `heat_conduction_a` once the thick and thin passes exist. No action now; the difference is measured as algorithm, not as port.

**Reporting a regression.** The gate writes `results/<run>/summary.json`, `kernels.csv` (kernel, PASS/FAIL/NO-TEST, comparison type and rule, threads, build configuration, run time, mutants caught) and the logs of every step. A failing run is reported to the generator owner and the Chief Architect with the kernel name, the failing configuration (flag set, callee switch, threads), the first failing case line and the run folder; V&V logs it on the sign-off list and does not edit the generator. A golden-signature change is accepted only after review and `--update-golden`, and the gate is then run again.

### 5.11 FR-005 (iii): switch settings of the decomposition-invariance tests (requirements v0.4.31, D-053)

Rule (FR-005 (iii)): decomposition-invariance tests (box split, rank count, thread count) run with the exact fixed-point sum selected: **`EXACT_SUMS=.TRUE.` in `&MISC`** (written "T" below). In default mode (`EXACT_SUMS=.FALSE.`, written "F") the per-step sums (`USUM`, `DSUM`, `PSUM` and the other zone sums, `RAD_Q_SUM`, `KFST4_SUM`, HVAC node sums) keep the FDS summation order, so agreement across layouts of results that depend on them is within eps_H only. Setup-time areas and volumes (D-028) are exact in both modes. Min and max reductions (`CHECK_STABILITY`) are order independent in both modes. `EXACT_SUMS` is the working name given by the Spec Lead in requirements v0.4.33 (IR-003; logical, default F: F keeps the FDS order for per-step sums, T selects the exact fixed-point sum; setup-time sums are always exact, D-028, D-053); it becomes final with the User Guide draft. Case copies for the rows marked T set `EXACT_SUMS=.TRUE.`; the rows marked F leave it at its default.

| Test (where) | Varies | `EXACT_SUMS` | Criterion |
|---|---|---|---|
| FR-005 (i) reduction-free explicit stages, frozen inputs (G11) | box split, ranks 1/2/4/8, threads | T | byte-identical (reduction-free stages do not depend on the switch; T is the recorded setting) |
| FR-005 (ii) zone sums and gauge, `zone_break_fast`, `zone_shape` (G11, §5.3) | box split (two `max_grid_size`), ranks, threads 1/2/4 | T | `D_PBAR_DT` and gauge constant byte-identical |
| same configurations, default mode (F) | as above | F (default) | zone sums and the results that depend on them within eps_H across layouts; not byte-identical; reported, gating on eps_H |
| FR-005 (iv) run-to-run, repeats at fixed ranks, threads, layout, `OMP_DYNAMIC=false` (G11) | repeat only | T and F | bitwise in both modes |
| FR-005 (v) setup areas and volumes, (v-a) to (v-c) (§5.3) | box split, ranks, threads | T (setup sums are exact in both modes) | T0 across configurations, T1 against single-mesh FDS |
| Radiation byte identity across ranks at fixed layout, FR-062 (G11, §7) | ranks | T | byte-identical; box-split dependence exempt (D-039) |
| AMR device values and Smokeview-format files across ranks and a redistribution, FR-070, FR-076 (G13) | ranks, redistribution | T | byte-identical / `cmp`-identical |
| Pressure single solves, FR-005 (iii), FR-037, FR-039 (G14) | box split, ranks | not applicable (frozen input, no per-step sum) | eps_H in both modes |
| Whole runs across ranks (G11, NFR-010/011/012) | ranks | F (default) | T2 uniform, T3 AMR |
| Explicit-stage bitwise comparison across ranks and threads (NFR-010/011/012) | ranks, threads | T | byte-identical; `H` within eps_H in either mode |

Rows that do not vary box split, rank count or thread count (G0 to G10, G12, §5.1, §5.2, §5.4 to §5.9) have no `EXACT_SUMS` setting. GPU acceptance: for the MP5 `DIVG` path it is tolerance-only, not bitwise against the baseline, until upstream patches UP-0001 and UP-0002 (the fourth `Z_TEMP` element) are committed upstream and merged (NFR-043); the host side of the G15 gate (§5.10) stays bitwise against the verbatim upstream loop text and is not affected.

### 5.12 Phase 3 acceptance: two-level runs and the D-057 `ns2d_16` check (v0.4.10)

Machine-readable twin: `vv-runs/phase3/phase3_cases.csv` (41 rows; columns id, phase, gate, title, purpose, input_source, tier_runtime, observables, comparator, tolerance_class, threshold, threshold_status, dependencies, owner, code_owner, spec_refs, notes). The CSV is generated by `vv-runs/phase3/make_phase3_cases_csv.py`; edit that script, not the CSV. Where this section and the CSV disagree, fix the script and regenerate.

#### 5.12.1 Scope, legend, and what can run when

- **Phase 3 content (roadmap):** regrid, interface transport and the FR-024 interface flux overwrite, with prescribed or simplified velocity; no walls, no OBST, no composite pressure solve. The composite pressure solve and the level-1 `TimeLoop` binding are Phase 4. The D-062 and D-063 regrid rulings have their own rows in §5.12.8. So the pressure halves of FR-016(c) and (d) and the named `species_conservation_1..4` / `Energy_Budget_*` refined cases are listed as later-phase rows, non-blocking now. Phase 3 stands on the transport, conservation and regrid checks plus the D-057 single-level pressure check.
- **Gate column:** *blocking* = must pass for the Phase 3 exit; *non-blocking (Phase N gate)* = written now, gates at Phase N.
- **Threshold status (v0.4.12):** **SPEC** = number or rule taken from the requirements (some SPEC numbers are themselves marked "proposed" in the spec: FR-020/021 1e-12 and 1e-10, FR-032 1e-8, the T1 constant). **ACCEPTED-A-58** = a V&V value that was PROPOSED in v0.4.10 and is now the working Phase 3 threshold by the Spec Lead ruling A-58. The two values that A-58 left open until the Architect confirmed them (D-059 corner limits; T3 rule on the tracer slice) are confirmed by the Architect (D-071 (f)) and carry **CONFIRMED-D-071**. Accepted does not mean frozen: a gate result that shows a threshold has no power or no margin is reported back.
- **Inputs:** `vv-runs/phase3/inputs/` (blob inputs, generated by `make_blob_inputs.py`; derived from `Verification/Scalar_Analytical_Solution/move_slug.fds`), `vv-runs/phase3/d057/` (D-057 inputs, oracle).
- **Not runnable yet in the AMR code (dependencies column in the CSV):** a prescribed-velocity (`FREEZE_VELOCITY`) mode in AMR mode (`TimeLoop::check_scope` rejects it); static patches from a finer `&MESH` (the driver aborts with the M2a "all meshes must have the same cell size" message); `TAG_*` keywords in `&AMR` (the parser accepts `MAX_LEVEL`, `REF_RATIO`, `REGRID_INTERVAL`, `BLOCKING_FACTOR`, `MAX_GRID_SIZE`, `N_ERROR_BUF`, `N_PROPER`, `GRID_EFF`, `OUTPUT_LEVEL_CAP`, `VELOCITY_TRANSFER`, nothing for tags); the level-1 `TimeLoop` binding; test hooks (flux overwrite off, composite-sum and face-flux dumps, hierarchy dump, RHS/H dump and replay). The input files carry no invented `TAG_*` keywords: working names `TAG_SPEC_ID='SLUG'` and `TAG_SPEC_DIFF` (undivided difference, threshold about 0.05) are noted in a comment in `blob2d_amr.fds`.
- **Reference runs done now (GNU Release `refbin/gnu_ompi_firex-36975d7/fds`, 1 rank, `nice`, `timeout`):** measured numbers below come from these. Wall times were taken with the machine load at 12 to 19 on 8 cores, so they overstate an idle-machine time.

#### 5.12.2 Conservation (FR-012, FR-020, FR-021, FR-024; D-050, D-061)

Quantities, all summed over **uncovered** cells of all levels with the D-028 exact sum: total mass M_tot = Σ ρV; species mass M_s = Σ ρ Y_s V for each species; energy / enthalpy E_h = Σ ρ h V. The interface check is the FR-024 flux identity: after the overwrite the coarse interface face flux equals the area-sum of the covered fine face fluxes, per stage, advective and diffusive, per species. There is no flux register in the design (D-050); the sum check replaces register bookkeeping. Within the periodic frozen-velocity test boxes the net boundary flux is zero, so the conserved totals must stay constant. Thresholds are the spec numbers (SPEC) except where marked.

| id | Gate | Check | Threshold | Status |
|---|---|---|---|---|
| P3-C01 | blocking | M_tot per step and cumulative, static patch (`blob2d_static_mr_amr`) | per step ≤ 1e-12; cumulative ≤ 1e-10 (relative) | SPEC (FR-020) |
| P3-C02 | blocking | M_s per species, same run | same numbers per species | SPEC (FR-021) |
| P3-C03 | blocking | E_h, isothermal frozen blob | same numbers | ACCEPTED-A-58 (FR-024 gives "as in FDS", no figure; the isothermal frozen case has no source terms, so the species-mass form applies) |
| P3-C04 | Phase 4 | energy closure, non-isothermal blob with pressure work | max(1.05·\|closure_base\|, 0.005·max\|Q_TOTAL\|) (FR-022 form) | SPEC form; no Phase 3 hot-blob number (accepted A-58: deferred to Phase 4) |
| P3-C05 | blocking | FR-024 face flux identity, ADV and DIF, both stages | (a) vs in-code ordered sum: bitwise; (b) vs test-side permuted sum: ≤ 1e-14 relative | ACCEPTED-A-58 (2 to 16 terms; a few ulp) |
| P3-C06 | blocking | composite sums identical across 1/2/4 ranks and two `max_grid_size` (`EXACT_SUMS=T`, §5.11) | bitwise | SPEC (FR-005, D-028) |
| P3-C07 | blocking | realizability: clip count; max \|ΣY−1\| | 0; ≤ 1e-14 | count SPEC (FR-025); ΣY ACCEPTED-A-58 ("at round-off" in the spec; a few ulp × species count) |
| P3-C08 | blocking | negative control: overwrite OFF must fail C01/C02 | cumulative species imbalance > 1e-6 (inverted verdict) | ACCEPTED-A-58 (4 orders above the pass limit) |
| P3-C09 | Phase 5 | `species_conservation_1..4` with a patch (FR-021 named) | ≤ 1e-2 (case Tol); baseline `species_conservation_1` 9.12e-3, margin 9% | SPEC; needs walls, VENT |
| P3-C10 | Phase 4–5 | `Energy_Budget_*`, `simple_duct`, `mass_balance_*` with patch | FR-022 form | SPEC |

eps_H is not used in this section: it covers a single solve on frozen input with the same discretisation, not conservation sums and not coarse-fine comparisons (§2). The tolerance classes enter as follows: T0 for decomposition independence and the identity (a); the REQ numbers of FR-020/021 for the sums; T2 later for whole-run budgets.

**Evidence that the check can fail (recorded, GNU Release, FDS own multi-mesh interface, `blob2d_static_mr`, 12 coarse 10×1×10 meshes + one 40×1×40 mesh, 678 steps):** Total mass drift 3.8e-14, BACKGROUND species 5.3e-3, SLUG species 5.1e-2 (max over the run). So the Total looks conserved while the species are not, which is why C02 is separate from C01 and why the species numbers are the sensitive observable. Role 3's unit test without the overwrite shows a tracer drift of 4.5e-4 over 128 steps.

#### 5.12.3 FR-016 (`int_1to2`) mapped to concrete comparisons

FR-016 as recorded in `requirements.md` (v0.4.32): the static two-level case is `Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement` (13 meshes: one 16×16 fine patch and 12 coarse 4×4 meshes, ratio 2, fully periodic, 2-D); the reference is its derived UGLMAT-HYPRE copy `_uglmat`. (a) first fine-side ghost layer, bitwise, first step, from the A-09b dump (§5.2); (b) the coarse side is **not** matched, conservation (FR-020/021) applies instead; (c) pressure: (i) interface normal-velocity mismatch at machine zero (FR-032, ≤ 1e-8 m/s), (ii) error against the exact solution no worse than UGLMAT, (iii) observed order ≥ FDS order − 0.1; compared at t≈1 and t≈2π only; (d) full-run T2, Phase 4. eps_H does not apply to (c) and (d).

| id | Gate | Comparison | Threshold | Status |
|---|---|---|---|---|
| P3-F01 | blocking | (a) first fine-side ghost layer, RHO, ZZ, TMP, RSUM, RHOS, ZZS, MU, KRES, D, DS on `int_1to2` vs A-09b dump, step 1 | bitwise; coarse side report-only | SPEC (step 1 with D, DS confirmed, A-58; Role 3 to align) |
| P3-F02 | blocking | D-059 shared corner cells: differing set ⊆ declared set | fine side: none; KRES coarse side: differing cells only in the declared D-059 set (expected 4 of 16 corner-zone cells); RHO, TMP coarse ≤ 3e-15 relative; a gate: if exceeded, D-059 is reopened and the limit is not widened | set rule SPEC (D-059); 3e-15 CONFIRMED-D-071 |
| P3-F03 | blocking | covers nothing: `ns2d_16_amr_notag` (AMR mode, no tags) vs same code without `&AMR`; vs recorded `ns2d_16__sf17` | bitwise; T2 (e ≤ 1.05·e_base + 0.1·Tol; e_base RMS u 0.1307); 0 level-1 boxes | SPEC (T0, T2) |
| P3-F04 | blocking | covers everything, transport: `blob2d` fully refined vs AMR-code uniform fine; level 0 vs average_down(level 1); AMR single-level vs FDS `blob2d_f80` | bitwise; ≤ 1e-14 relative; T1 (1e-10·max\|x_base\|) | ACCEPTED-A-58 (T0 inside the AMR code, T1 against FDS; FREEZE_VELOCITY, no pressure solve) |
| P3-F05 | Phase 4 | (c)(i) normal-velocity mismatch | target machine zero; limit ≤ 1e-8 m/s | SPEC (FR-032) |
| P3-F06 | Phase 4 | (c)(ii) L2 error vs exact ≤ 1.1 × UGLMAT | table below | SPEC (rule), recorded numbers |
| P3-F07 | Phase 4 | (c)(iii) observed order ≥ FDS order − 0.1 | table below | SPEC (rule), recorded numbers |
| P3-F08 | Phase 4 | (d) full-run T2, `int_1to2` vs UGLMAT baseline | RMS u (0 to 2π) ≤ 1.05 × 0.2981 = 0.3130 | SPEC |
| P3-F09 | Phase 4 | covers everything, pressure: `ns2d_16` fully refined vs baseline `ns2d_32` | RMS u ≤ 1.05 × 0.0319; H vs single-level fine solve within eps_H(32) = 1e-8 | SPEC |

"Covers nothing" and "covers everything" are therefore split by what the Phase 3 code can do: the empty hierarchy is a real `ns2d_16` run (about 0.85 s unloaded) and gates now; the full-coverage hierarchy gates now on transport only, with the pressure side in Phase 4 (P3-F09).

Recorded numbers for P3-F06 and P3-F07 (A-38, UGLMAT HYPRE, `vv-runs/A-38/norms.csv`; gate bound L2 ≤ 1.1 × these):

| Quantity | Time | L2 at N = 16 / 32 / 64 | Observed order 16→32, 32→64 |
|---|---|---|---|
| u | t≈1 | 4.152e-2 / 1.324e-2 / 4.929e-3 | 1.649, 1.426 |
| w | t≈1 | 3.774e-2 / 1.196e-2 / 3.937e-3 | 1.658, 1.603 |
| H (mean removed) | t≈1 | 1.107e-1 / 4.557e-2 / 2.108e-2 | 1.281, 1.112 |
| u | t≈2π | 2.524e-1 / 1.448e-1 / 1.0905e-1 | 0.802, 0.409 |
| w | t≈2π | 2.537e-1 / 1.316e-1 / 1.0236e-1 | 0.946, 0.363 |
| H | t≈2π | 5.225e-1 / 2.921e-1 / 2.240e-1 | 0.839, 0.383 |

FFT and UGLMAT agree to ≤ 3.9e-4 relative up to 2π. t = 30 is not usable (the solution is unresolved). Runtimes (recorded): `ns2d_16` 0.85 s, `ns2d_32` 4.0 s, `ns2d_64` 28.8 s; `int_1to2` 90.9 s (default solver), 29.7 s (UGLMAT); A-38 UGLMAT wall at N = 16/32/64: 30.5 / 141.3 / 917.2 s (1 rank).

**Two gaps to close for P3-F01** (from Role 3's ghost check, `src/Source/regrid_transport/notes/fr016-ghost-check.md`): its dumps use steps 2 and 3 (`FDSREF_STEPS=2,3`), while FR-016(a) and §5.2 say the first step; and it does not compare D and DS, which §5.2 lists. Ruling A-58: the gate stays at step 1 with D and DS (it matches the spec); Role 3 is to align the ghost check (`FDSREF_STEPS` must include step 1; D and DS added).

#### 5.12.4 Moving blob (regrid tracking)

**Input.** A tracer slab advected across a periodic box with frozen velocity, derived from `Verification/Scalar_Analytical_Solution/move_slug.fds` (`FREEZE_VELOCITY` plus `&WIND`, SUPERBEE, two slugs: SLUG mass fraction 1 on [0.125, 0.375]² and 0.5 on [0.5, 0.75]²; BACKGROUND elsewhere). The exact solution is the initial field shifted by u·t, so the field returns to the initial state at t = 1 for the 2-D case (u = w = 1). No existing FDS case serves unchanged: `move_slug` itself has the right physics but is a fixed-grid case without an `&AMR` line and without conservation output. Files, all in `vv-runs/phase3/inputs/`:

| File | Use | Reference-binary result (1 rank, loaded box) |
|---|---|---|
| `blob2d_c40`, `blob2d_f80` | 2-D coarse / fine single level | c40: 339 steps, 4.8 s stepping; f80: 678 steps, 124 s wall; mass drift ≤ 6.6e-15 / 3.3e-14 |
| `blob2d_c16`, `blob2d_f64` | ratio-4 pair | 0.8 s; 26.8 s wall; Total and SLUG drift ≤ 8.6e-15 / 1.7e-14 |
| `blob2d_amr` | 40×1×40 level 0 + `&AMR MAX_LEVEL=1, REF_RATIO=2, REGRID_INTERVAL=4, BLOCKING_FACTOR=4, N_ERROR_BUF=2` | AMR code only (the FDS binary does not read `&AMR`) |
| `blob2d_static_mr`, `blob2d_static_mr_amr` | FDS multi-mesh layout like `int_1to2` (12 coarse 10×1×10 + one 40×1×40 over [.25,.75]²); AMR variant with `REGRID_INTERVAL=0` | static_mr: 43.9 s wall; the C08 evidence above |
| `blob3d_c32`, `blob3d_f64`, `blob3d_amr` | 3-D, one slug, velocity (1, 0.5, 0.25) | c32: 19.3 s, drift ≤ 2.1e-15 (Total), 1.6e-14 (SLUG); f64: 483 s wall, drift ≤ 1.4e-14 |

All use `&DUMP SIG_FIGS=17, DT_SLCF=0.25, MASS_FILE=.TRUE., DT_MASS=0.05` and cell-centred SLUG slices. Tools: `run_blob_smoke.sh` (one rank, `nice -n 19`, `timeout`), `blob_metrics.py` (relative L2 and max error against the exact shifted field with box-overlap fractions, peak, mass, circular centroid, `_mass.csv` drift). Smoke runs show the FDS single-level runs conserve the mass file to round-off; relative L2 error against the exact field at t = 0.25 / 0.5 / 0.75 / 1.0 is c40 0.268 / 0.265 / 0.302 / 0.308 and f80 0.186 / 0.206 / 0.219 / 0.223; the 2-D ratio-4 pair c16 0.393 / 0.446 / 0.484 / 0.540, f64 0.210 / 0.239 / 0.222 / 0.248. The centroid error is a constant ~0.11 fine cell for both resolutions (a one-step time-stamp offset, not drift).

**Refinement criterion.** FR-011 / D-058: undivided difference of the SLUG species with `TAG_KEEP` hysteresis, threshold about 0.05 (working value, keyword TBD), buffer `N_ERROR_BUF = 2`. Buffer rule: `N_ERROR_BUF ≥ ceil(R · CFL per axis) + 1`; the blob moves about 0.5 fine cell per regrid interval here, so 2 is safe.

**Expected behavior.** Level 1 appears around the slabs at t = 0, follows the edges as they move, and tags never leave the finest level; the composite mass and species masses are constant to round-off including across each regrid.

| id | Gate | Observable and comparator | Threshold | Status |
|---|---|---|---|---|
| P3-B01 | blocking | per-regrid composite mass and per-species change; per-step and cumulative imbalance; clip count; containment count (cells the criterion would tag that are not on the finest level, sampled every step) | ≤ 1e-12; ≤ 1e-12 / ≤ 1e-10; 0; 0 | regrid and sums SPEC (FR-012, FR-020/021); containment 0 ACCEPTED-A-58 |
| P3-B02 | blocking | discrimination: relative L2 of SLUG slice, \|\|AMR−F\|\| ≤ 0.5·\|\|C−F\|\| at t = .25/.5/.75/1.0 | ≤ 0.084 / 0.089 / 0.089 / 0.090 as computed from the FDS C and F (‖C−F‖/‖F‖ = 0.168 / 0.177 / 0.178 / 0.180, coarse prolonged 2×2); recompute from AMR-code C and F at the gate | gate; CONFIRMED-D-071 (T3 rule on the tracer slice) |
| P3-B03 | blocking | circular centroid vs the fine run | ≤ 0.1 fine cell | ACCEPTED-A-58 |
| P3-B04 | blocking | determinism: hierarchy dumps and final fields across repeats and 1/2/4 ranks | identical / bitwise (`EXACT_SUMS=T`) | SPEC (FR-015, FR-005) |
| P3-B05 | blocking | negative control `REGRID_INTERVAL=0` | containment > 0 and discrimination bound exceeded (inverted verdict) | containment part ACCEPTED-A-58; discrimination part CONFIRMED-D-071 |
| P3-B06 | blocking | `RemakeLevel` on an unchanged grid | bitwise copy | SPEC (FR-012) |
| P3-B07 | blocking | 3-D blob (`blob3d_amr` vs `blob3d_c32`, `blob3d_f64`) | as B01 and B02 (discrimination is a gate) | conservation SPEC; discrimination CONFIRMED-D-071 (gate) |
| P3-B08 | blocking (conservation); discrimination report-only | ratio 4 (`blob2d_c16`, `blob2d_f64`) | as B01; discrimination (B02 rule) is **report-only** at ratio 4 | conservation SPEC; discrimination report-only (D-071 (f)) |
| P3-B09 | blocking | D-059 corner effect: (1) cells farther than d(n) from the patch outline bitwise equal to the uniform-fine AMR run for the first n steps, with d(n) = n × S × h in cells of the fine level, S = stages per step, h = dependence radius of one stage; **S = 2, h = 2, so d(n) = 4n fine cells, confirmed by Role 3 (A-61): S = 2 stages, h = 2 cells in `GET_SCALAR_FACE_VALUE` for every limiter including MP5, frozen velocity**, counted inward from the outline; the report must **print the observed front distance** so the bound can be checked for tightness; (2) corner-zone max difference ≤ edge-zone max difference at the same distance + 4 ulp | (1) bitwise outside d(n) = 4n fine cells; the observed front distance is printed in the report and must not exceed d(n); (2) corner max ≤ edge max + 4 ulp | (2) CONFIRMED-D-071; (1) rule CONFIRMED-D-071, value 4n confirmed by Role 3 (A-61) |

Notes on use: the comparisons of B02, B07, B08 are made between runs of the **AMR code** (level 0 only, fully refined, AMR); the FDS runs supply the exact-field errors and the check that the AMR-code uniform runs reproduce FDS single level (T1). A fixed FDS slice plane (3-D) is left by the blob after t = 0.25 (it moves in y and z), so the FDS-slice comparison in 3-D is valid at t = 0.25 only; the 3-D gate uses the AMR driver's volume output and composite sums. The metrics script normalises the 3-D plane error by the t = 0 plane norm for this reason. The FR-016 shared-corner effect (D-059) is covered by P3-F02 and P3-B09. FR-013: the same cases run once on a GNU Debug AMR build with zero assertion failures (P3-X01).

#### 5.12.5 D-057 check: 2-D `ns2d_16` (with Role 2, Pressure Backend Implementer)

**What D-057 decides** (README decision log and `requirements.md`): AMReX `PoissonHybrid` fails on fully singular problems, so all-periodic and all-Neumann cases use `FFT::Poisson`; a direction with a single cell is ignored by the FFT solver. Role 2 checks the mapping of FDS BC types (`FISHPAK_BC`, `&PRES`, surface and vent types) to the FFT and MLMG BCs for 2-D and singular cases, including mean removal; the Pressure Solver Lead reviews it; the V&V Lead supplies the 2-D test `ns2d_16`; Role 1 supplies the driver side.

**Inputs.** `ns2d_16` exists in the baseline set (`case_inventory.csv`, `vv-runs/baseline/gnu_ompi_firex-36975d7/ns2d_16` and `ns2d_16__sf17`): committed `Verification/NS_Analytical_Solution/ns2d_16.fds`, `PERIODIC_TEST=1`, 16×1×16, x and z periodic, y walls, dx = dz = 2π/16, dy = 0.1. The recorded `SIG_FIGS=17` outputs are reused (DEVC UVEL, PRES, VISC; 296 rows). New derived inputs, in `vv-runs/phase3/d057/` with derivation diffs: `ns2d_16_dy10.fds` (y extent × 10) and `ns2d_16_amr_notag.fds` (adds an `&AMR` line, no tags). 1 rank.

| id | Gate | Check | Comparator and threshold | Status |
|---|---|---|---|---|
| P3-D01 | blocking | BC mapping, setup only: `fds_amr ns2d_16.fds --pressure-bc` | prints `twod=1 n=16x1x16 codes=1,3,1 x=PP y=NN z=PP`, no ERROR; FFT chosen | SPEC; string verified once with a scratch driver binary (also for `ns2d_16_amr_notag`, `blob2d_amr`, `blob3d_amr`) |
| P3-D02 | blocking | single solve on frozen input: backend H (FFT and MLMG) vs `poisson_oracle.py`; vs the baseline FDS FFT H from the A-09b/refdump instrumented dump of step 1; FFT vs MLMG | relative L2 ≤ eps_H(16) = max(1e-8, 2.4e-12·256) = 1e-8, mean removed both sides | SPEC (eps_H) |
| P3-D03 | blocking | mean removal on the singular problem | removed_rel = \|removed mean\|/‖b‖ ≤ 1e-10 per solve; no true-residual warning | ACCEPTED-A-58 (1e-10 is the T1 scale; FR-039 asks for the diagnostic without a number) |
| P3-D04 | blocking | whole run vs the recorded baseline | T2 on UVEL, PRES, VISC DEVC (e ≤ 1.05·e_base + 0.1·Tol); RMS u error 0 to 2π ≤ 1.05 × 0.1307 = 0.1372; PRES compared after removing a constant (D-067) | SPEC (T2) |
| P3-D05 | blocking | thin-direction invariance: `ns2d_16_dy10` vs `ns2d_16` in the AMR code | T1, ≤ 1e-10·max\|x_base\| per series | SPEC (T1 constant is "proposed" in the spec) |

**Oracle and its validation.** `poisson_oracle.py` is an independent numpy solver for one cell-centred 7-point Poisson solve with periodic / Neumann / Dirichlet per direction (a one-cell direction contributes no term); it removes the mean on singular problems and prints `ORACLE rel_l2=... eps_H=... removed_rel=... verdict=PASS|FAIL`. Checked against Role 2's harness (`pb_harness`, 16×1×16): FFT periodic 4.2e-15, MLMG periodic 7.1e-14, FFT Neumann 7.6e-15, MLMG Neumann 1.3e-12 (the harness itself reports FFT vs MLMG 7.1e-14 and 1.3e-12 against eps_H 1e-8). Negative controls fail as required: wrong dx gives 0.20, BC mismatch gives 0.43.

**Thin-direction evidence (reference binary).** `ns2d_16_dy10` against the recorded baseline: UVEL differs by ≤ 3.1e-15 (max\|UVEL\| = 2.2), PRES by 2.2e-16, VISC by 2.6e-26; so FDS itself is invariant to dy at round-off, not bitwise: hence T1, not T0. FDS cannot run a one-cell x or z (ERROR(426), "Poisson initialization error") or a y-boundary VENT in a 2-D calculation (ERROR(809)); those variants are not FDS inputs. The one-cell, Dirichlet-refusal and y-open mapping cases stay in Role 1's unit test `pressure_bc_map` (21 checks). The 23 inputs with pressure code 0 on non-periodic directions are FDS-only in every mode (D-074), which also bounds what a blob input may contain: the blob inputs use code 0 only on periodic directions.

**Gauge (D-067, requirements v0.4.33).** The default gauge is `sum(rho*V*(KRES-H)) = 0` per pressure zone and per connected component, with the exact sum over uncovered cells (the ULMAT/UGLMAT convention); it is also applied to FFT-solved cases, where `PRES` is compared after removing a constant. The default mean removal is the composite volume-weighted mean (D-032); the FDS arithmetic removal remains a run-time parity switch. The recorded baseline FFT H slice has a non-zero mean (about 1.84 to 1.89 for `ns2d_16`) and UGLMAT's gauge is about 0 relative to the exact solution (FFT about +0.46 for `int_1to2`), so every H and `PRES` comparison here is made after removing a constant. The recorded H slices are float32 node slices (17×17): too coarse for eps_H = 1e-8, so P3-D02 needs the A-09b/refdump dump of `M%H` in float64.

**Questions for Role 2:**
1. Answered by D-067 (see Gauge above): the backend applies `sum(rho*V*(KRES-H)) = 0`; `PRES` is compared after removing a constant. Open for Role 2: is the composite `gauge_weight`/offset in place for the single-level FFT path, and does the harness print the applied offset?
2. How is a one-cell y handled in the MLMG backend: `setHiddenDirection` with `ref_ratio_vect = 2 1 2`, or plain Neumann?
3. Will the harness or backend expose an RHS and H dump / replay hook for `ns2d_16` at step 1?
4. Is `removed_rel` available per solve as a diagnostic?
5. Which instrumented-baseline H dump is used (the A-09b refdump)? Which step?
6. When does the composite path exist, and does a one-cell x/z with a Dirichlet face stay refused?

#### 5.12.6 Status of the former PROPOSED values (ruling A-58) and open decisions

A-58 accepted all values below as the working Phase 3 thresholds; the two rows it left open are confirmed by the Architect (D-071 (f)).

| Item | Row | Value | Status |
|---|---|---|---|
| Energy / enthalpy closure for an isothermal frozen blob | P3-C03 | same form as FR-020 (1e-12 per step, 1e-10 cumulative) | ACCEPTED-A-58 |
| Hot-blob energy number in Phase 3 | P3-C04 | none; FR-022 form in Phase 4 | ACCEPTED-A-58 |
| Interface flux identity | P3-C05 | bitwise vs in-code sum; ≤ 1e-14 vs permuted sum | ACCEPTED-A-58 |
| Species sum | P3-C07 | \|ΣY−1\| ≤ 1e-14 | ACCEPTED-A-58 |
| Negative-control margin | P3-C08 | cumulative imbalance without the overwrite > 1e-6 | ACCEPTED-A-58 |
| T0 vs T1 for full coverage and single-level blob comparisons | P3-F04 | T0 inside the AMR code, T1 against FDS | ACCEPTED-A-58 |
| D-059 corner limit | P3-F02, P3-B09 | coarse side 3e-15; corner max ≤ edge max + 4 ulp; P3-F02 is a gate (if exceeded, reopen D-059, never widen); B09 (1) distance d(n) = n × S × h = 4n fine cells (S = 2, h = 2; confirmed by Role 3, A-61); the observed front is printed | CONFIRMED-D-071 (rule); value 4n confirmed by Role 3 |
| Containment count and centroid bound | P3-B01, P3-B03 | 0; 0.1 fine cell | ACCEPTED-A-58 |
| Discrimination rule applied to the tracer slice | P3-B02, B07, B08 | ‖AMR−F‖ ≤ 0.5‖C−F‖: gate for B02 and B07, report-only for B08 (ratio 4) | CONFIRMED-D-071 |
| FR-016(a): step 1 with D and DS | P3-F01 | step 1 with D, DS (spec text) | SPEC, confirmed A-58; Role 3 to align |
| Mean-removal bound | P3-D03 | removed_rel ≤ 1e-10 | ACCEPTED-A-58 |
| Gate phase of F05 to F09 | P3-F05 … F09 | Phase 4 | ACCEPTED-A-58 |
| Numeric value of "pressure tolerance" for P3-R02 | P3-R02 | working bound 1e-9·U/dx_fine at solver tolerance 1e-12; derived bound in §5.12.8 | working bound accepted (D-071); derivation to be measured |

Spec-vs-plan wording to align: Role 3's plan (`docs/role3-regrid-transport-plan.md`) uses "YAFluxRegister / reflux"; the spec (D-050) uses the interface flux overwrite with no register. This plan follows the spec.

#### 5.12.7 Reference builds used and their status

All runs above used `vv-runs/refbin/gnu_ompi_firex-36975d7/fds` (Release). Verified: `sha256sum -c SHA256SUMS` passes for 36 entries, and a rerun of `ns2d_16` at `SIG_FIGS=17` gives devc and hrr CSVs byte-identical to the recorded baseline (P3-X02). The Debug binary (`fds_debug`) has landed and its checksums pass, but it has not been run by V&V. `refbin/impi_intel_firex-36975d7` has no `fds` binary yet (build in progress). `refbin/gnu_ompi_firex-bee11f0329` is not built. Open gate row P3-X01: a Debug AMR driver build is needed for FR-013 on the Phase 3 cases.

#### 5.12.8 Regrid rulings D-062 and D-063: species transfer and post-regrid projection (added v0.4.12; requirements v0.4.33, FR-012 (6) and (7))

Why the V&V blob inputs are not enough: `blob2d_*`/`blob3d_*` use a uniform divergence-free velocity, so a regrid creates no velocity seam (D-063 has nothing to correct), and their two species have the same molar mass, so density is uniform and mass-weighted and linear transfer of Z give the same numbers (D-062 cannot fail). New inputs for D-062: `blob2d_mw_c40`, `blob2d_mw_f80`, `blob2d_mw_amr` (SLUG molar mass 4, BACKGROUND air; density varies between 0.166 and 1.165 in the blob; same geometry and velocity as `blob2d_*`; generated by `make_blob_inputs.py`). Reference-binary smoke run of `blob2d_mw_c40` (GNU Release, 1 rank, loaded box): 339 steps, 18.6 s wall, Total mass drift 1.3e-14, SLUG drift 9.0e-15, divergence ≤ 5e-14. The relative L2 error of the SLUG mass-fraction slice against the shifted initial field is 0.48 / 0.49 / 0.51 / 0.52 at t = .25 / .5 / .75 / 1.0 (larger than the uniform-density case because the limiter acts on rho and rho·Y separately); these inputs are therefore used for the conservation and transfer checks only, not for the discrimination rule. For D-063 the seam needs a non-uniform divergence-free velocity on every level; the evidence for Phase 3 comes from Role 3's regrid tests (`src/Source/regrid_transport/tests/test_blob_registry.cpp`: 3-D 32³ two levels ratio 2, 2-D 48×1×48 three levels, 2-D 16×1×16 `ns2d_16` style, 2-D ratio 4), not from an FDS input.

| id | Gate | Check | Threshold | Status |
|---|---|---|---|---|
| P3-R01 | non-blocking (report-only in Phase 3) | After every regrid that retains old fine faces next to new ones: max\|div u − D\| over uncovered cells before and after the projection hook `project_after_regrid`, solver iterations, largest velocity change at the seam; also whether the hook ran or the AUTO setting skipped it | report only (the number is printed and recorded); values reported by Role 3, not rerun by V&V: 3-D 32³ ratio 2 1.0 → 7e-14 (≤ 9 iterations); 2-D 48×1×48 three levels 1.0 → 2.5e-13; 2-D 16×1×16 1.4 → 6e-13; ratio 4 3.5 → 1e-12 (report only: first-order coarse-fine gradient at ratio-4 corner patches) | SPEC (FR-012 (6) interim rule) |
| P3-R02 | blocking once the composite MLMG is in the gate build; report-only until then | same observable, asserted: max\|div u − D\| right after each regrid (with or without the hook running) ≤ the pressure tolerance on uncovered cells | ≤ the pressure tolerance, as the bound **max\|div u − D − c\| ≤ 10·eps_rel·B + 20·eps_mach·U/dx_fine** (Pressure Lead derivation; eps_rel = solver relative tolerance, B = max\|div u* − D\| measured before the projection, or U/dx_fine a priori, eps_mach = 2.2e-16, c = the constant removed by the D-032 mean removal). At eps_rel = 1e-12 and B ≈ U/dx_fine this is about 1e-11·U/dx_fine. The **working bound 1e-9·U/dx_fine is accepted (D-071)**; the measured-B form is recommended. The derivation is not a measurement yet: it is to be measured in the A-56 runs; it applies to the MLMG bottom solve, and a HYPRE bottom solver must have tolerance ≤ eps_rel | SPEC form; working bound ACCEPTED (D-071); derived bound to be confirmed by measurement |
| P3-R03 | same gate as P3-R02 | negative control: projection off (`POST_REGRID_PROJECTION='OFF'`, working name) and a solver that corrects nothing must each fail P3-R02 at every regrid with a seam | seam values 1.0 (3-D), 3.5 (ratio 4), 1.4 (16-cell case) exceed the bound; the control must be reported as FAIL | SPEC (FR-012 (6) negative control) |
| P3-R04 | blocking once P3-R02 asserts | the projection does not change composite mass or species mass | per regrid ≤ 1e-12 relative (recorded: 0 per regrid, 3e-16 over the run, with and without projection) | SPEC (FR-012 acceptance) |
| P3-S01 | blocking | D-062 mass-weighted species transfer: restriction, average-down, prolongation act on rho and rho·Z; Z derived. (a) unit level: Role 3's `test_species_avgdown` (bare level, registry, coarse-fine ghost hook); (b) end-to-end `blob2d_mw_amr` (regrid every 4 steps) when the driver can run it | per-regrid and cumulative change of M_tot and of each species mass ≤ 1e-12 per regrid, ≤ 1e-10 cumulative; rho, ZZ, ZZS on covered cells bitwise equal to the driver registry transfer (as reported by Role 3 for commit 7c85f23539) | SPEC (FR-012 (7), FR-020/021 numbers) |
| P3-S02 | blocking | negative control for P3-S01: linear averaging / interpolation of Z must fail P3-S01 | must fail the 1e-12 bound; expected size on Role 3's test data about 1e-1 relative (the earlier driver run showed about 1e-4); V&V records the value measured on `blob2d_mw_amr` at the gate | SPEC (FR-012 (7) negative control) |

Gate rules taken from the ruling: a regrid without projection support is not a gate case (masked cells, non-box level 0, domain walls other than periodic or Neumann, stretched or cylindrical meshes, and the 2-D three-level case with the mock solver are reported only); gate numbers come from ratio 2. In the driver, the hook is not wired into the time loop yet (dependency for the end-to-end rows). The ruling's default for `POST_REGRID_PROJECTION` is for the implementer to confirm.

## 6. Case selection (details and per-run estimates in `case_inventory.md`)

**Tier 1 (39 runs, ≈ 17 min sequential, all measured on the GNU FireX reference; `dancing_eddies_1mesh` moved to Tier 2 per spec R-0).**
- Interfaces/pressure: `obst_activation_{default,ulmat}`, `divergence_test_2/3`, `dancing_eddies_{1mesh,default,uglmat_refine}`, `obst_coarse_fine_interface`, `lapse_rate`.
- Refinement analogues: `random_meshes`, **`ns2d_16_int_1to2_refinement` + derived `VV/…_uglmat`** (FR-016 gating case, after A-09), **`ns2d_16_emb_1to2_refinement`** (informational, non-gating; FR-014).
- Interface invariance: `soborot_superbee_square_wave_128{,_1mesh}`, `shunn3_4mesh_128`.
- Convergence: `ns2d_{8..64}{,_nupt1}`, `saad_512_cfl_*`, `shunn3_{32,64,128}`.
- Conservation: `energy_budget_tmix`, `species_conservation_1/2`.
- Fire: `1_step_2_step_compare`.
- Restart: `restart_test1a/b` + derived continuous run.
- Determinism pair.

**Tier 2 (127 runs, ≈ 89 min sequential as a lower bound; now includes `dancing_eddies_1mesh` and the A-24 single-mesh copies of `dancing_eddies_default`, `shunn3_4mesh_32` and `symmetry_test_mpi`)** adds the pressure-solver variants, the remaining refinement analogues (`dancing_eddies_embed`, `duct_flow_uglmat_refine`, `ns2d_16_emb_1to1_refinement`, `porous_media`, `race_test_1` (AMR mode: A-35 `_r4` copy)), the full MMS / scalar / limiter convergence series, multi-mesh symmetry and volume-flow cases, pressure zones, mass-balance cases, fires, restarts, and every anchor case not in Tier 1 (`energy_budget_particles`, `zone_break_fast`, `zone_shape`, `cascadempi`, `cannon_ball`, the two radiation anchors, `openmp_test64a`, `device_restart_base_case`). **Optional (8):** `tunnel_demo`, `energy_budget_adiabatic_walls`, `race_test_4` (AMR mode: A-35 `_r4` copy), `shunn3_512`, `ht3d_energy_conservation_4`, `VV/test_8mesh_NOFRPG_short` (HYPRE cost, NFR-033), and the Intel-only `*_uglmat_pardiso` pair.

**Adaptive_Mesh_Refinement folder (scope item; A-09).** `ns2d_16_int_1to2_refinement` (13 meshes / 448 cells, fine 16×16 patch abutting 12 coarse 4×4 meshes, ratio 2) and `ns2d_16_emb_1to2_refinement` (2 meshes / 512 cells, fine mesh embedded in the coarse) are **Tier 1**. `ns2d_16_emb_1to1_refinement` (2 / 320, same-resolution control) is **Tier 2**. All are 2-D, periodic, T_END = 30 s, with `VELOCITY_TOLERANCE=1e-6` and up to 100 pressure iterations. Estimates: ~10–30 s each; the 13-mesh case is 0.5–5 min on 1 rank, depending on iteration counts. Memory is negligible.

*Why commented out:* all three lines (FDS_Cases.sh 4-6) come from commit `042aa2624d` (R. McDermott, 2024-11-08), *"updates to amr test cases, demonstrating errors in both embedded mesh strategy and interpolated refinement"*. That commit also created the `int` case. They were added deliberately disabled, as demonstrations of known error. They are not broken or slow. There is no dataplot row or Guide section, and nothing has changed since. A-09's "pass" is the V&V definition in §3.4 step 1 (`case_inventory.md` §4a).

**NFR-032 plume (Tier 3, Phase 1 S1 and Phase 10; D-008, A-19).**
- Level 0 / coarse (decided 2026-09-25 by Spec Lead, Integration Lead and Chief Architect): the committed `Validation/Heskestad_Flame_Height/FDS_Input_Files/Qs=1_RI=10.fds` (33×33×80) trimmed by one cell on the high x/y side to **32×32×80**, `XB=-1.8,1.690909,-1.8,1.690909,-0.45,8.59` (dx = dy = 3.6/33 and dz = 0.113 m unchanged; AMR level 0 `blocking_factor=8`, `max_grid_size = 16 16 8` per the D-024 amendment: 40 boxes, 5 per rank on 8 ranks). The trim keeps the burner at 9×9 cells; an evenly centred 32-cell grid would snap it to 10×10 (+23 % area, changes Q*). V&V input with the probes: `vv-runs/inputs/A-19/Qs1_RI10_coarse_32x32x80.fds` (81,920 cells), also the ADR-002 S1 input (`CFL_FILE`).
- Uniform-fine reference: `vv-runs/inputs/A-19/Qs1_RI10_fine_64x64x160.fds`, same XB, 655,360 cells in 8 meshes of 64×64×20 (81,920 each) stacked in z. Estimates (**not measured**): ≈ 2.3 h on 8 ranks (range 0.5–4.5 h), ≈ 1.6 GB total (range 1–2.5 GB); coarse ≈ 50 min on 1 rank.
- Burner check (setup-only FDS runs, GNU Release, 2026-09-25): committed 33×33×80 and new 32×32×80 both snap the burner to 9×9 cells (`.smv` OBST indices `12 21 12 21 3 4`), the fine input to 18×18×2 (`24 42 24 42 6 8`); footprint ±0.49091 m, 0.963967 m² on all three, so the area-adjusted HRR (1512.7 kW nominal) and Q* are unchanged. Details, changes, probe verification (`verify_A19_probes.py`) and assumptions: that folder's `README.md`. Superseded v1 (66×66×160 / 33×33×80) is in `vv-runs/inputs/A-19/superseded/`.
- Requirements v0.3.2 records this as D-024. Minor text points for the Spec Lead (A-19 README): `max_grid_size=16` gives 20 level-0 boxes, not 32 (now moot: the D-024 amendment sets `max_grid_size = 16 16 8`, 40 boxes); a 10×10 snap would change burner area and Q*, not total HRR (FDS area-adjusts HRRPUA); memory ≈ 1.4 GB there vs ≈ 1.6 GB here (both estimates).
- The V&V alternative `Fires/circular_burner.fds` (R-8) stays as a fallback only (D-020 open); requirements' own fallback is `McCaffrey_14_kW_5.fds`.

**FireX-only folders.**
- **`GPU_Tests`**: 9 inputs, all 1 M cells, propane fire, UGLMAT/HYPRE, plus HPC job scripts. The physics runs on CPU (`HYPRE_DEVICE_RUN` is compiled out without `WITH_HYPRE_DEVICE`). The job scripts don't apply. The 16/32-mesh variants need more than 8 ranks. Full runs take 20–45 min each, so Tier 3, except one 20-step copy (T2-opt) that serves NFR-033 and FR-038.
- **`VTK`**: 4 CPU inputs demonstrating the VTKHDF writer. Their data output needs an HDF5 build (D-019). They turn off `.smv` output and have no pass/fail criterion, so they are Demonstration only. Relevant to FR-072/074 and Q5.

**Refinement-ratio cross-check (A-33 inventory vs our lists; v0.4.3, D-030).** Sources: `docs/inventory/mesh_ratio_cases.csv` and `mesh_ratio_cases_rollup.csv`, crossed with the 184 rows of `case_inventory.csv` (read-only). 180 are in the rollup; the 4 derived rows are not, and they inherit from their sources.
- Worst class per case: single_mesh 99, uniform 70, supported_2_4 4, no_shared_faces 3, embedded_or_overlapping 2, other_nonpow2 2.
- **Only two of our cases use a ratio outside {1,2,4}:** `Thread_Check/race_test_1` (T2) and `race_test_4` (T2-opt), each with five 5:1 faces. They need the A-35 remesh for AMR mode (done, §5.6); the originals stay the FDS baselines.
- None of our cases has a 3:1, direction-dependent or TRN-stretched interface.
- The embedded/overlapping cases have no shared-face ratio and stay FDS-baseline-only for the multi-level comparison (D-015; FR-014/FR-043 use them as AMR-code checks, not as ratio inputs):
  - `ns2d_16_emb_1to2_refinement` (T1), `ns2d_16_emb_1to1_refinement` (T2);
  - `random_meshes` (T1; worst class uniform, but 8 overlapping pairs at 1:1 and 2:1);
  - `dancing_eddies_embed` (T2; 2:1 overlap).
- Supported 2:1 cases need no action: `ns2d_16_int_1to2_refinement` (2-D, y ignored), `dancing_eddies_uglmat_refine`, `obst_coarse_fine_interface`, `duct_flow_uglmat_refine`.
- Of the 43 requirements §2.3 anchors, only race_test_1/4 are outside {2,4}.

**Shaped-OBST area defect (FM_Burner; FR-005 (v)).** FDS's `AREA_ADJUST` for CYLINDER OBSTs split across meshes depends on the mesh-to-process map (upstream issue candidate #1).
- The 12 OBST-cylinder FM_Burner validation inputs deliver 0.52083 × the specified fuel when all meshes are on 1 process, and exactly 1.00000 with 1 mesh per process (the upstream layout) or on a single mesh.
- None of them is in our lists. They are FR-005 (v) test inputs (§5.3), not FDS baselines, unless run as a single-mesh copy or 1 mesh per process.
- Details: `fm_burner_area_check.md`.

Full assessment: `case_inventory.md` §6. **Excluded, with reasons:** `case_inventory.md` §7 (GEOM, >8 ranks, hour-plus runs, local physics sub-models covered by G4).

## 7. Requirement → test traceability (class used)

| Req | Test / cases | Class | Phase |
|---|---|---|---|
| FR-001 | single-mesh anchors (ns2d_16/32, shunn3_32, species_conservation_1, energy_budget_*, cannon_ball, radiation anchors, openmp_test64a, restart_test1*) | Explicit kernels on frozen input T1 (T0 per Q1); whole runs T2 (D-022) | 2 |
| FR-002 | multi-mesh same-res anchors (shunn3_4mesh_32 first at M2a, layer_4mesh, dancing_eddies_default, symmetry_test_mpi, duct_flow, …) | Kernels on frozen input T1 (T0 per Q1) vs an equivalent single-mesh input (committed `layer_1mesh`; derived copies for shunn3_4mesh_32, dancing_eddies_default, symmetry_test_mpi, A-24); whole runs T2 vs single-mesh and vs multi-mesh FFT and GLMAT baselines; `H` single frozen solve vs GLMAT within **eps_H** (G14) | 2 |
| FR-003 | G4 | T2 | 2, 10 |
| FR-005 | G11 (1/2/4/8 ranks, two `max_grid_size`, two thread counts, 3 repeats; `zone_break_fast`/`zone_shape` for (ii); §5.3 cylinder reproducer, FM_Burner CH4 2cm + single-mesh copy, `zone_shape`/`random_meshes` volumes for (v)) | (i) byte-identical explicit stages; (ii) exact zone sums/gauge from the per-box sum; (iii) eps_H; (iv) bitwise repeat; (v) setup areas/volumes T0 across configurations, T1 vs single-mesh FDS (§5.3). Radiation stage (D-039): exempt from (i) for level box-split dependence only; at a fixed box layout it stays rank-count independent, (iv) applies unchanged, and `RAD_Q_SUM`/`KFST4_SUM` use the (ii) exact sums; `EXACT_SUMS=T` for (i), (ii), (v) across layouts, default mode (`F`) within eps_H only, (iv) in both modes (§5.11) | 2 |
| FR-004, FR-044, IR-002 | negative tests | SETUP_STOP + message | 3 |
| FR-010/013/015 | hierarchy dumps, debug assertions, repeat runs; ratios {2,4} only (MLMG limit ≤ 4; D-030); §5.4 negative tests N-1..N-4 (3:1, direction-dependent, 5:1, legacy), positive controls P-1/P-2, blocking-factor report, MLMG coarsening warning W-1; level-0 NIC=1 assertion P-3 and NIC-25 negative control N-5 (D-076) | exact box lists (T0-equivalent); setup stop + named pair/ratio; warning text + coarsest grid | 3 |
| FR-011 | per-criterion tagging unit cases | exact cell sets | 3 |
| FR-012 | regrid mass/species deltas (P3-B01, P3-S01/S02); mass-weighted species transfer with linear-Z negative control (D-062); post-regrid composite projection, report in Phase 3, assert with the composite MLMG, projection-off negative control (D-063; P3-R01 to R04) (§5.12.8) | 1e-12 rel (req.); max\|div u − D\| ≤ pressure tolerance | 3 (report), 4 (assert) |
| FR-014 | ns2d series + patch, ns2d_16_*_refinement (AMR code only; baseline numbers are context) | refined ≤ coarse, factor 2.0; order ≥ FDS order − 0.1 and L2 ≤ 1.1× FDS at each N on the same 2:1 case (D-040; 1.8 reported, non-gating) | 4 |
| FR-016 | G7 (`int_1to2` + UGLMAT copy only, after A-09; `emb_1to2` informational) | T0 (a, component-level, first ghost layer only, via A-09b; §5.2); (c) FR-032 machine zero + exact-solution error ≤ UGLMAT + FR-014 order gate (D-040; 1.8 reported only), no eps_H; T2 (d); multi-step never bitwise (D-023); HYPRE commit as §3.5 | 3–4 |
| FR-020/021/024 | species_conservation_1..4, simple_duct, mass_balance_* with patch | req. numbers (round-off with refluxing, D-023) | 3 |
| FR-022 | energy_budget_tmix/particles/dns_100/adiabatic_walls with patch | max(1.05·\|closure_base\|, 0.005·max\|Q_TOTAL\|) | 4–5 |
| FR-023, FR-031, FR-032 | divergence_test_1..3, obst_activation, CHECK_POISSON output | req. numbers; FR-031 tolerance TBD(Pressure Lead). eps_H is not a residual criterion | 4 |
| FR-030/033/035 | duct_flow, tunnel_demo, dancing_eddies_default, ns2d periodic; pressure_iteration3d_default, random_obstructions_fft | T2 | 4–5 |
| FR-034 | zone_break_fast, zone_shape (+ zone_shape_2, zone_break_fast_uglmat_hypre) | T2 | 4/6 |
| FR-036 | shunn3_128 (PRESSURE_TOLERANCE 1e-6, analytic H) + helium_2d_isothermal | T2 + iterations ≤ baseline + 1 | 4 |
| FR-037 | FR-001 cases with the FFT path selected; single frozen solve vs baseline FFT (G14, §5.1) | D + eps_H | 2 |
| FR-038 | G6 | T2 (v0.3.1), same HYPRE commit (§3.5) | 2 |
| **FR-039** | G14 / §5.1: FFT vs MLMG single solve on ns2d_16/32, shunn3_32, csmag_32, dancing_eddies_1mesh; switch run `VV/ns2d_32_switch`; maxorder 2 order study (D-032) | eps_H (single solve, after gauge); exact switch steps; FR-014 order gate (D-040) | 2 (P2), 4 |
| FR-040/041/042/043 | mask checker; burning-OBST regrid case; obst_activation_default, box_burn_away1; random_meshes, duct_flow | req. / T2 | 5 |
| FR-050/051/052 | particle ledger; energy_budget_particles, cascadempi; FR-052 cascadempi + OBST-only geom_sprk_mass copy | exact counts / T2 | 7 |
| FR-060/061 | plate_view_factor_cart_30, radiating_polygon_square_20 | T2 (060); 061 via FR-001/002: kernels T1, whole runs T2 | 2, 8 |
| FR-062 | `radiation_gas_panel` split across boxes/ranks vs single-mesh FDS, `RADIATION_ITERATIONS` K = 1, 2, 3 (A-50); same case at 1/2/4 ranks with a fixed box layout, repeated | Box-split lag error: class set with FR-062's verification (D-039); fixed box layout: byte-identical across rank counts and run-to-run (FR-005 (iv)); `EXACT_SUMS=T` (§5.11) | 8 |
| FR-070 | all Tier 1 DEVC outputs; regrid-jump check; AMR mode (D-045, ADR-004 D6): one multi-level case at 1/2/4 ranks with a fixed box layout, plus a MINLOC/MAXLOC tie case | T2 uniform (D-022); T3 jump; AMR device values byte-identical across rank counts; ties resolved to the lowest global cell key; `EXACT_SUMS=T` for the rank-count comparison (§5.11) | 9 |
| FR-071 | header diff on all Tier 1/2 CSVs | exact (I) | 9 |
| FR-072 | Smokeview load of Tier 1 cases; AMR mode: ADR-004 spike S-C case (≥ 2 levels, ≥ 1 regrid in the output window) with slice, boundary, 3-D smoke and particle data in the refined region; format checker (every `.smv` file exists, header index bounds inside the named mesh's `GRID`) | D, SMV 6.11.2; zero load errors; checker passes | 9 |
| FR-076 | ADR-004 spike S-A (A-53): two-level case (static fine box + one regrid in the output window), fixed box layout, at 1, 2 and 4 ranks, plus one run with a forced redistribution between output times (D-045) | `.smv`, `.sf`, `.bf`, `.s3d`, `.q`, `.prt5` `cmp`-identical across all runs; `EXACT_SUMS=T` (§5.11) | 9 |
| FR-078 | `AMR_LEVEL` slice on the S-A case at every output time; Smokeview view (S-C) | Slice equals the level layout exactly (I); D for the Smokeview pattern | 9 |
| FR-074 | VTK/*.fds (HDF5 build; D-019 open) | identical VTKHDF | 9 |
| FR-080/081 | restart_test1a/b + continuous, device_restart_* (+ base case), restart_ulmat_*, csvf_restart_a | FR-080: T0 proposed, T1 minimum; baseline-difference fallback only if baseline misses T1; AMR T1 + identical hierarchy. T2 (081) | 9 |
| IR-001 | G1 | message diff | 2 |
| IR-004 | random_meshes, layer_4mesh (MPI_PROCESS inputs) | behaviour + warning text | 3 |
| IR-006 | G5 | T0 | 2 |
| NFR-010/011/012 | G11 (NFR-011 AMR mode on the A-35 race_test `_r4` copies, FDS baselines on the originals; §5.6) | explicit stages bitwise and `H` within eps_H across ranks (FR-005); whole runs T2 uniform / T3 AMR; watchdog; `EXACT_SUMS=T` for the bitwise explicit-stage comparison, `H` within eps_H in either mode (§5.11) | 2+ |
| NFR-020 | both toolchains build and pass FR-001 | blocked (A-18); GNU FireX baseline exists | 2 |
| NFR-021 | HYPRE/AMReX pins recorded; FireX HYPRE = `63331f19c7` (v2.32.0-24) per provenance | I | 1 |
| NFR-030/031 | openmp_test64a + 8-rank multi-mesh anchor; A-19 memory; timing rules (1)-(6) incl. solver, TinyProfiler pressure time, forced-MLMG level-0 run, separate profiled install (§5.7) | ratio, median of 3, idle machine; clean-install headline | 2, 10 |
| NFR-032 | Heskestad Qs=1_RI=10, level 0 32×32×80 trimmed + A-19 64×64×160 fine reference (D-008) | T3 + ≤ 50% wall | 10 |
| NFR-033 | test_8mesh_NOFRPG 20-step, MPI_PROCESS 1/2/4/8 ranks; secondary duct_flow_uglmat_refine | ≥ 60% efficiency (proposed) | 10 |
| NFR-040 | harness over the anchor set using T0–T3 | D | 2 |

## 8. Tooling (delivered and planned)
- **Delivered:**
  - `vv-runs/tools/patch_check.sh CONTROL_FDS PATCHED_FDS [--debug CTL_DEBUG PATCHED_DEBUG]`: standard bitwise behaviour-unchanged check for upstream patches (G0 pair, `shunn3_4mesh_32`, `csmag_32`, optional extra case; SIG_FIGS=17 copies; timing columns stripped); usage and pass criterion in the script header.
  - `vv-runs/gpu_gate/` (`run_gpu_gate.sh`, `gate.py`, `check_device_logs.py`, README, `results/`): the GPU kernel bitwise gate of §5.10.
  - `parse_fds_inputs.py`: static survey. Follows FDS `CHECKREAD` rules and expands MULT including SKIP ranges.
  - `make_inventory.py` and `render_inventory_md.py`: inventory, cost model, class mapping.
  - `compare_csv.py`: CSV comparison with `--class T0|T1` presets. Implements §2 T0/T1 exactly: timing-column exclusion patterns, T1 column-max scaling, `--exact`, per-column tolerances, interpolation for T2-style series, JSON reports. Tested on synthetic CSVs (`vv-runs/tools/testdata`); first real FDS output comes from the GNU capture.
  - `vv-runs/inputs/A-19/make_A19_inputs.py`: generator for the NFR-032 fine reference (64×64×160) and level-0/coarse input (32×32×80); `verify_A19_probes.py`: independent probe/device/burner-snap checker.
- **Planned (Phase 2, NFR-040):**
  - `run_case.sh` harness: pre-flight load/RAM check, watchdog, manifest (HYPRE version from provenance).
  - numpy re-implementations of the dataplot metrics (`end`, `max`, `mean`, `area`, `end_1_n`, `all`, `tolerance`, `slope`) and of the analytic-error scripts needed by Tier 1/2 (`ns2d.py` RMS error, `shunn_mms.py`, `saad_mms_temporal_error.py`, `soborot_mass_transport.py`, `mass_balance*.py`).
  - A T3 statistics script (time-window mean/RMS at probes; A-19 probe naming `T_z…`/`W_z…`).
  - A frozen-input `H` comparator (mean removal, relative L2, eps_H(N)) and a solver-selection log checker (G14).
  - A hierarchy-dump comparator once the dump format exists.

---

## 9. Status of the V&V items folded into requirements v0.3

The items below are this plan's v0.2 answers (R-n) and disagreements (D-n). Their rationale is in `archive/test-plan-v0.2.md` §9. Status is as recorded in requirements.md v0.3 and the `docs/README.md` decision log. v0.3.1 (D-022) changes the status of D-2 and D-3 (below). All accepted numbers are **provisional until calibrated** against the measured FireX baseline (A-08/A-10) and approved by the project owner.

### 9.1 Answers (R-n)

| # | Topic | v0.3 status | Where in v0.3 | Remaining V&V action |
|---|---|---|---|---|
| R-0 | Anchor runtimes, rank caveats | accepted | §2.3 runtime check; D-013 | Replace estimates with measured times |
| R-1 | T2 margin (floor, 3σ) | accepted | §2.2 T2; D-013 | σ ensemble (needs Intel) |
| R-2 | T3 values, discrimination rule | accepted | §2.2 T3, NFR-032; D-013 | Calibrate at Phase 4 |
| R-3 | FR-003 at two frequencies | accepted | FR-003; D-013 | — |
| R-4 | FR-014 factor 2.0 | accepted | FR-014; D-013 | — |
| R-5 | FR-022 floor | accepted | FR-022; D-013 | — |
| R-6 | FR-036 case shunn3_128 (+ helium_2d_isothermal) | accepted | FR-036; D-013 | — |
| R-7 | FR-052 case cascadempi (+ OBST-only geom_sprk_mass copy) | accepted | FR-052; D-013 | Create derived copy |
| R-8 | NFR-032 case circular_burner | **open, not adopted** | NFR-032 keeps Heskestad `Qs=1_RI=10` (D-008); **D-020** | A-19 inputs delivered (§6); alternative kept until A-19 / ADR-002 S1 report |
| R-9 | Smokeview 6.11.2 | accepted | FR-072; D-013 | Install needs authority over the development machine |
| R-10 | NFR-011 watchdog | accepted | NFR-011; D-013 | Implement in `run_case.sh` |
| R-11 | FR-080 restart rule | **adjusted** | FR-080: T0 proposed, **T1 minimum**; "no worse than baseline's own restart difference" only if baseline misses T1; D-013 | Measure baseline restart-vs-continuous at capture (G10) |
| R-12 | NFR-033 scaling case | accepted | NFR-033; D-013 | Create 20-step MPI_PROCESS copies |

### 9.2 Disagreements and clarifications (D-n)

| # | Topic | v0.3 status | Where in v0.3 | Remaining V&V action |
|---|---|---|---|---|
| D-1 | T0/T1 on `SIG_FIGS=17` derived copies | accepted | §2.2 T0/T1; **D-014** | Copies under `vv-runs/` |
| D-2 | Chaotic FR-001 anchors: T1 100 steps + T2 | **open in the decision log (D-017, Q1)**; **superseded in v0.3.1 text**: FR-001 whole runs are T2 for every anchor, and kernels are T1/T0 per Q1 (D-022) | FR-001; D-017 (log not yet updated) | Spec Lead to close D-017 or restate what Q1 still decides (kernel T0 vs T1) |
| D-3 | FR-002: T1 first 10 steps, then solver-tolerance bound | **open in the decision log (D-018)**; **resolved in v0.3.1 text**: full-run T1 vs GLMAT dropped ("resolves V&V's objection"), whole runs T2, kernels vs single-mesh T1/T0, `H` eps_H (FR-002, D-022) | FR-002; D-018 (log not yet updated) | Spec Lead to close D-018; FR-031 tolerance still needed for FR-023/FR-032 |
| D-4 | FR-016 = `int_1to2` only; `emb_1to2` informational | accepted | FR-016; **D-015** | — |
| D-5 | UGLMAT reference copy; instrumented baseline | accepted | FR-016, A-09 (UGLMAT copy run first), **A-09b** (instrumented build), R-34; D-015 | §3.4 step 1, §3.6 |
| D-6 | GPU_Tests inputs run on CPU | accepted | §2.3 note; D-013 | — |
| D-7 | VTKHDF baseline needs HDF5 build | **open** | **D-019**, depends on A-14/Q5; one option (HDF5 option in `AMReX` branch) conflicts with D-005 | No HDF5 build until decided |
| D-8 | Ranks = min(8, meshes) | accepted | §2.3, NFR-010/011; **D-016** | — |
| D-9 | T2 absolute floor / chaotic margin | accepted | §2.2 T2; D-013 | — |
| D-10 | FR-022 floor | accepted | FR-022; D-013 | — |
| D-11 | FR-080 T0 vs baseline | **adjusted** | as R-11 | as R-11 |
| D-12 | Timing valid only on idle machine | accepted | NFR-030/033; D-013 | — |
| D-13 | `device_restart_base_case` anchor; continuous restart companion | accepted | §2.3; D-013 | Create continuous run |

## 10. Open questions and actions for V&V
- **Blocking:**
  - (1) oneAPI runtime on the development machine (A-18, R-32). The OpenMPI runtime is fine (§3.1).
  - (2) Intel FireX build per §3.3 (Intel Build Chief), with the **same HYPRE commit `63331f19c7`** as GNU.
  - (3) Smokeview install (FR-072).
- **Needs a decision from others:**
  - Q1: now only kernel-level T0 vs T1 (FR-001/002). D-017/D-018 are superseded by D-022 in the v0.3.1 text; the decision log still lists them open.
  - Pressure Lead's FR-031 tolerance → FR-023, FR-032.
  - Kernel-parity dumps from baseline FDS (FR-001/FR-002 under D-022): extend the A-09b instrumented build or specify a separate one (Integration Lead, Pressure Lead, Legacy Mapper). A-24 derived single-mesh copies: owner not stated in the requirements; V&V can make them.
  - Q5 and A-14 → D-019 (VTK/HDF5 baseline).
  - D-020: Spec Lead, after the A-19/S1 report.
  - Q9 (MPI_PROCESS in AMR mode).
  - Spec Lead: HYPRE label in requirements §2.1/NFR-021 and R-30 ("3.0.0" vs the linked v2.32.0-24), and whether the AMR build keeps `63331f19c7` (§3.5).
  - Spec Lead/project owner: whether A-09 on the GNU Release + Debug builds is enough to start treating FR-016 as a gate while the Intel build is blocked.
  - ~~Integration Lead / Chief Architect: the D-008 level-0 grid 33×33×80 is odd in x/y~~ **Resolved 2026-09-25 (D-024, requirements v0.3.2)**: level 0 = 32×32×80 trimmed on the high side, `blocking_factor=8`, `max_grid_size = 16 16 8` (D-024 amendment, 40 boxes; §6).
  - Integration Lead: frozen-input solve hook and solver log for G14 (§5.1).
  - (v0.4.3) Spec Lead: FR-005 (v) says the exact setup sums "reproduce single-mesh FDS". V&V applies T0 across AMR configurations and T1 against single-mesh FDS, because FDS's sequential float sum need not equal the exact sum rounded once (§5.3). Confirm or tighten.
  - (v0.4.3) Integration Lead: a setup diagnostic that prints `FDS_AREA`/`AREA_ADJUST` per shaped OBST/MULT and the zone volumes, also in setup-only mode. FDS computes neither in setup-only runs (§5.3).
  - (v0.4.3) Integration Lead / input converter: the FR-010 ratio check must ignore the single y cell of 2-D inputs, or `ns2d_16_int_1to2_refinement` (FR-016) is rejected (§5.4 P-1).
  - (v0.4.3) Chief Architect, for the record: `amr.blocking_factor = 2 8` is valid only with a ratio-4 level 1 (AMReX needs bf0·ratio ≥ bf1: 2·2 < 8 fails). A ratio-2 hierarchy on the same level 0 needs `2 4` or `2` (§5.6; A-35 README).
- **V&V actions (from the requirements action log):**

| ID | Action | Status |
|---|---|---|
| A-08 | Record the FireX builds; capture baselines | GNU build recorded (§3.5); G0 + Tier 1 capture running (`baseline_status.md`); Intel blocked (A-18) |
| A-09 | Run `int_1to2` + UGLMAT copy first (Release + Debug, both toolchains), then the `emb` cases; FR-016 gates only after (R-34) | open; GNU Release/Debug binaries available |
| **A-09b** | Instrumented ghost-cell / coarse-fine dump build from a scratch copy, never committed; accepted only if T0 vs uninstrumented build on all Tier 1 CSVs (with GNU Build Chief) | open; placeholder spec §3.6; detailed spec needs Pressure Solver Lead + AMReX Integration Lead input |
| A-10 | Calibrate provisional numbers | open until measured baseline |
| A-14 | VTK/HDF5 decision support | open (D-019) |
| **A-19** | NFR-032 uniform-fine reference `Qs=1_RI=10` (with Pressure Solver Lead); ADR-002 S1 with `CFL_FILE`; uniform-fine runtime/memory | **v2 inputs generated and setup-checked (T_END=0), no time-stepping run**: `vv-runs/inputs/A-19/Qs1_RI10_coarse_32x32x80.fds` (level 0 / S1, 81,920 cells) and `Qs1_RI10_fine_64x64x160.fds` (655,360 cells, 8 × 64×64×20), trimmed XB `-1.8,1.690909,-1.8,1.690909,-0.45,8.59`; burner 9×9 / 18×18 cells, 0.963967 m² on committed, coarse and fine (§6). Estimates: fine ≈ 2.3 h on 8 ranks, ≈ 1.6 GB. v1 66×66×160 in `superseded/`. S1 and measurements wait for run authorisation and the A-26 exclusive window (D-024; path and grid now match requirements v0.3.2) |
| **A-35** | Remesh `race_test_1` (and hand-check/remesh `race_test_4`) for NFR-011 in AMR mode (D-030) | **inputs done 2026-09-25**, not yet run in AMR mode (needs M3). Ruling confirmed: mesh 3 → `IJK=24,24,16`, all 5 interfaces 4:1, faces tile. race_test_4 has identical meshes, so it gets the identical remesh (option A of 4). Level 0 34×18×32, `blocking_factor 2 8` or `2`. Setup-checked (28,656 cells vs 37,440). Side effects: burner top/vent 0.120 → 0.125 m, 32 → 16 burning cells with fuel flow preserved, probe cell moved. `vv-runs/inputs/A-35/README.md`. Next: FDS runs of the `_r4` copies at 1 and 4 threads as the AMR-mode reference (short runs only until run authorisation). |
| A-33 (input) | Legacy Mapper ratio scan, cross-checked against our lists | closed by the Legacy Mapper. V&V cross-check done (§6): only race_test_1/4 are outside {2,4}. |
| FM_Burner (R/FR-005 v) | Area-adjust check of the FM_Burner inputs | done: `fm_burner_area_check.md`. 0.52083 × fuel with all meshes on 1 process, 1.00000 with 1 mesh per process or single mesh. Not in our lists; used as an FR-005 (v) input. |

- **Also from v0.3.1:** A-24 (derived single-mesh copies of `shunn3_4mesh_32`, `dancing_eddies_default`, `symmetry_test_mpi`) and FR-005 verification (V&V is verifier for (i)-(iv)); A-25 (global-reduction list) is not V&V's.
- **V&V next steps:** GNU G0 → G1 sweep → A-09 (GNU Release + Debug) → Tier 1 capture (GNU) → Tier 2 → Intel once A-18 is done → calibration ensemble → replace model estimates with measurements and re-tier.

## 11. Files
- `docs/vv/test-plan.md` (this file, v0.4.3); `docs/vv/archive/test-plan-v0.3.md` (previous version); `docs/vv/archive/test-plan-v0.2.md` (full §9 rationale)
- `docs/vv/fm_burner_area_check.md` (FM_Burner / FR-005 (v) area check); work files in `vv-runs/inputs/fm_burner_area/`
- `vv-runs/inputs/A-35/{race_test_1_r4.fds, race_test_4_r4.fds, check_mesh_interfaces.py, rt*_check.txt, README.md, setup_check/}`
- `docs/vv/environment.md` (v0.2 content; still says no FireX binary exists and HYPRE 3.0.0. To be updated from §3.5)
- `docs/vv/case_inventory.md` and `case_inventory.csv` (v0.2 content; e.g. `emb_1to2` still marked "FR-016 disputed", no A-19 rows)
- `docs/vv/verification_case_survey.csv` (FireX) and `verification_case_survey_ce1f659.csv` (v0.1)
- `vv-runs/gpu_gate/{run_gpu_gate.sh, gate.py, check_device_logs.py, README.md, vendor/, results/}` (§5.10)
- `vv-runs/tools/{parse_fds_inputs.py, make_inventory.py, render_inventory_md.py, compare_csv.py, testdata/}`
- `vv-runs/smoke/{run_smoke.sh, inputs/, gnu/, intel/}`
- `vv-runs/inputs/A-19/{Qs1_RI10_fine_64x64x160.fds, Qs1_RI10_coarse_32x32x80.fds, make_A19_inputs.py, verify_A19_probes.py, README.md, setup_check/, superseded/}` (superseded v1 66×66×160 / 33×33×80 files in `superseded/`)
- Baseline status: `baseline_status.md` (maintained by the capture worker, not by this plan)
- `vv-runs/phase3/{phase3_cases.csv, make_phase3_cases_csv.py, inputs/ (blob inputs, make_blob_inputs.py), d057/ (ns2d_16_dy10, ns2d_16_amr_notag, poisson_oracle.py, make_d057_inputs.py), run_blob_smoke.sh, blob_metrics.py, runs/}` (§5.12)
