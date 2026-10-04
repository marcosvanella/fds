# Review of the ported-status candidates (D-072 (c)): device = host evidence

Owner: V&V Lead (reviewer of record for device results in `tools/ported.toml`; the Architect co-signs the list). Status: v0.1.
Subject: `amrex/ported-candidates.md` (docs commit 90e808a), section 2, and the evidence it cites (`amrex/stage1-gpu-spike-plan.md` 7.9c, 7.10a, 7.10c, 7.10d, 7.13, 7.13a, 7.13b; stage-1 repository `prototypes/s4_cuda_mass/stage1`, including its newer `device_gate/` logs; `vv-runs/gpu_gate`).
All source trees and the Integration Lead's files were read only. This file and the gate-side files named in section 6 are the only things written.

## 1. Result in one table

| # | Candidate group (section 2 of the candidates file) | Decision | One-line reason |
|---|---|---|---|
| 1 | `vpred_us/vs/ws`, `vcorr_u/v/w`, `div2_pred`, `div2_corr`, `baro_p_rrho`, `baro_fvx/fvy/fvz`, `vflux_vort_tau` (13) | **ACCEPT** (with the notes in 3.1) | Device = serial = 4 threads on every output array, real fields (`csmag_32`; the four baroclinic kernels also on the four `dec2_obst` boxes); the other nine kernels have no `dec2_obst` device run, which the candidates file and the draft map note overstate. |
| 2 | `vflux_fvx/fvy/fvz` | **ACCEPT** (notes in 3.2) | Real fields and real edge tables, device = serial = 4 threads = FDS on all faces; the device ran `csmag_32` and `dec2_obst` box 0 only. |
| 3a | `vn_max` | **ACCEPT** | Real fields, all outputs and locations equal; no untested branch. |
| 3b | `cfl_max` | **ACCEPT** (notes in 3.3) | Real-field runs cover norm 0 only; norms 1 to 3 are covered only by the synthetic rounds 4 to 7 hash run. |
| 3c | `div_extrema` | **ACCEPT** (notes in 3.3) | Independent numpy check agrees; cylindrical and stored-divergence branches are covered only by the synthetic hash run, and `CARTVELDIV` is all zero in every real-field case. |
| 4a | `d_z_max`, `del_rho_d_del_z` (SYNTHETIC `RHO_D`, `D_Z`) | **ACCEPT as device coverage, with the note in 3.4** | Bitwise on all outputs, non-vacuous on `dec2_obst`, host gate covers them; synthetic inputs only. |
| 4b | `rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat` (SYNTHETIC `RHO_D`, `D_Z`) | **CHANGES** | The host gate has no bitwise test for them (NO-TEST), so "device = host" is not tied to the FDS loop text; plus a real-field `RHO_D` / `D_Z` run (3.4). |
| 5 | `rho_d_dzd`, `h_rho_d_dzd` | **ACCEPT** (note in 3.4) | Same data and result as 4a (ZZ, TMP real; `RHO_D` synthetic); host gate PASS; also in gate set r1. |
| 6 | "WP1 non-wall and wall kernels" (266 and 76 arrays) | **CHANGES** | No kernel names; mixes real and substitute inputs; device outputs of the 7.9 run are no longer on the development machine, so I could not re-compare them (3.5). |
| 7 | `cfl_wall_max` (libm, 2 ulp) | **ACCEPT** for the ULP condition (section 4) | Gate ULP line: `maxulp 1` on the 32 real walls, `maxulp 2` on the 20,000 synthetic walls, `check_device_logs.py` says PASS. The other conditions are still open. |

Not candidates (agreed): the six `gsfv_*` kernels, loops without a kernel, `wall_rho_d_dzdn` `NIC>1` branch.

## 2. How I checked

For every device result I looked for: a device run on real fields; bitwise equality under the +0/-0 rule over all output arrays; host serial = host 4 threads = device; a tool that compares all bytes; a kernel that actually did work (non-vacuous, branches covered); inputs identified; build flags recorded.

1. **Independent re-comparison with the gate's own tool.** `vv-runs/gpu_gate/recheck_stage1_outputs.py` (new) walks every device output file of a stage-1 result tree and compares it, all bytes, with `bitcmp.py` (+0/-0 rule) against the host serial and the host 4-thread file; it also counts nonzero elements (vacuity). Output: `vv-runs/gpu_gate/device_logs/stage1_cflw/recheck_stage1_outputs.out`.

   | Result tree (`stage1/real/...`) | Device output files | Elements | Device vs serial failing elements | Device vs 4 threads | Serial vs 4 threads | Files that differ only in the sign of zero | All-zero device files |
   |---|---|---|---|---|---|---|---|
   | `results/wp3` (13 kernels) | 54 | 902,728 | 0 | 0 | 0 | 0 | 4 |
   | `res_wp3c` (8 + 2 kernels, WP3 item 3) | 155 | 490,764 | 0 | 0 | 0 | 22 (218,046 elements) | 12 |
   | `vfe_res` (`vflux_fv*`) | 12 | 241,656 | 0 | 0 | 0 | 0 | 0 |
   | `results/pf_all` (WP1 `dec2_obst`) | 44 | 74,768 | 0 | 0 | 0 | 0 | 8 |

   The gate's tool reproduces the harness's own counts (155 of 155 equal as values, 22 files / 218,046 elements sign-of-zero only). The harness comparisons themselves cover whole arrays (`compare_wp3.py`, `compare_wp3b.py`); the VFLUX summary quotes counts for the kernel region only, but the whole-file recheck above is also clean.
2. **Where the device ran.** From the result trees: `csmag_32` has device output for all six WP3 item 1 stages; `dec2_obst` boxes 0 to 3 have device output for `WP3_BARO` only. All other item 1 stages on `dec2_obst` are host only (serial and 4 threads).
3. **Gate state.** `check_device_logs.py` with the logs on record (section 5), and the latest host-tier runs under `vv-runs/gpu_gate/results` for the host leg.
4. **Build flags (read in `stage1/CMakeLists.txt`, `build_gpu.sh`, `README.md`).** Device: `nvfortran -O2 -mp=gpu -gpu=cc89,nofma -Minline`, `nvcc --fmad=false -Xcompiler=-ffp-contract=off`, cc 8.9; host: `gfortran -O2 -ffp-contract=off` (serial and `-fopenmp`, 4 threads). No fast-math option anywhere in these files. Not recorded per result file in the older results (it is in the build files and in 7.12); the newer `device_gate/*.log` headers record it per log. The FDS reference build is `-O3 -std=f2018 -frecursive -fopenmp` with no FMA instruction in the binary (7.9c).

## 3. Per group

### 3.1 Group 1: 13 velocity, divergence, baroclinic and vorticity kernels: ACCEPT
Evidence adequate for class bitwise:
- Real driver fields (`csmag_32`, 32^3, step 2, pass 1); device, serial and 4 threads equal on all 54 output files of `results/wp3` (recheck above); bit-equal to FDS where a truth exists (`vpred_us/ws` 31,744 of 31,744, `div2_pred` 32,768 of 32,768, `vcorr_u/w` 31,744 of 31,744 with corrector-time inputs, 7.9c).
- The generated source of these kernels has no branches (straight-line loops), so one non-vacuous real-field run exercises the arithmetic. Outputs on `csmag_32` are nonzero (`US/VS/WS` 33,792 elements, `OM*`, `T**` 35,937, `DIV` about 30,000).

Notes the Integration Lead must carry into the map and the `ported.toml` note:
1. **Wording is too broad.** The candidates file and the draft `ported.toml` note say "device = host serial = host 4 threads on every output array (`csmag_32`, four `dec2_obst` boxes)". On the device, `dec2_obst` ran only the baroclinic stage (4 boxes, 8 arrays each). `vpred_*`, `vcorr_*`, `div2_*`, `vflux_vort_tau` have device evidence on `csmag_32` only. Correct the note before it is entered.
2. **Baroclinic kernels are vacuous on `csmag_32`** (`FVX_B`, `FVY_B`, `FVZ_B` all zero, constant density). Their device evidence is the four `dec2_obst` boxes (`FVX_B/FVZ_B` 272 of 972, `FVY_B` 144 of 972 nonzero), inputs not physical in ghost and solid cells (equality test only, stated in 7.13).
3. **`vpred_vs`, `vcorr_v`, `div2_corr`, `baro_*`, `vflux_vort_tau` have no FDS truth** (no `FVY` dump; `FVY` input is zero, so the `FVY` term of `vcorr_v` is inactive). Their device = host result stands, but the host = FDS leg rests on the gate's host tier, which is PASS for all 13 in the latest quick run.
4. `vcorr_u/w` device inputs are the earlier FVX/FVZ snapshot; with the corrector-time values the host equals FDS and the kernel is unchanged (7.9c), so device = host on the old snapshot still holds.
5. The 7.13 compare names "arrays differing only in sign of zero" as none for this group; recheck agrees (0 files).

### 3.2 Group 2: `vflux_fvx/fvy/fvz`: ACCEPT
- Device = serial = 4 threads = FDS, all faces: `csmag_32` 33,792 of 33,792 per direction; `dec2_obst` box 0 272 / 512 / 272. Whole-file recheck: 12 of 12 files equal. Negative control with neutralised edge tables fails identically on device and host (10 of 272 faces), so the edge branch is exercised and decides the result at the obstruction edges.
- Notes: (1) **device ran `csmag_32` and `dec2_obst` box 0 only**; the 4-box statement in the candidates file is host. (2) `csmag_32` never takes the edge branch (0 of 65,345 table entries hold a value). (3) `dec2_obst` needs the MMS source switched off in a scratch tree; the run evolves differently; inputs and truth come from that run. (4) No case with open or interpolated-wall mesh-edge values (`ED_OMEGA` set by the wall BC); stated as not done in 7.13b. (5) `GX, GY, GZ, RHO_0` are constants (gravity zero in `dec2_obst`), so the gravity term is exercised only by `csmag_32`.

### 3.3 Group 3: `cfl_max`, `vn_max`, `div_extrema`: ACCEPT
- Real fields, `csmag_32` and four `dec2_obst` boxes, 10 device runs; all 155 outputs of WP3 item 3 equal (recheck: 0 failing elements), scalars and indices (`UVWMAX_TMP`, `ICFL_TMP`, `MU_MAX`, `I_VN`, `RESMAX`, `DIVMX/DIVMN` and locations) equal on all five cases. The tie rule (last or first index) is exercised in the sense that the locations agree; independent numpy check of `div_extrema` on `csmag_32` agrees.
- **Branches not run on real fields:** `cfl_max` only with `CFL_VELOCITY_NORM = 0` (norms 1, 2, 3 not run); `div_extrema` only `CYLINDRICAL = 0` and `STORE_CARTESIAN_DIVERGENCE = 0` (`CARTVELDIV` is all zero in all five cases, vacuous output). The rounds 4 to 7 hash run of the Mesh Data Loops work (`(local workspace)/scratch/devtier/r47_logs`, driver `gpu/mkdev_r47.py`: arguments `INORM`, `ICYL`, `IST`, scenarios with ties and signed zeros, sizes 96x40x56, 128^3, 192^3 and small scenario sizes, default and distribute-parallel-do builds) reports 0 DIFF for all 20 kernels; the gate reads those logs as PASS for `cfl_max`, `vn_max`, `div_extrema`. That evidence is SYNTHETIC and uses the hash form (`-0` and `+0` hash differently, so a pass is also a pass under the rule), and the build uses `mem:managed` (a bring-up aid, kernel_lint BF-06 WARN). I accept it as branch coverage, not as real-field evidence. If the owner of the map wants real-field branch coverage, a stage-1 rerun with the three scalar arguments set (device = host only) is cheap.
- Sign-of-zero rule applied: no file of this group differs only in the sign of zero.

### 3.4 Groups 4a, 4b and 5: synthetic `RHO_D` and `D_Z`
Facts: `RHO_D = 2e-5 x RHO` and the `D_Z` table are synthetic; ZZ, TMP, DS, RHO, H are real. All outputs bitwise (recheck: device = serial = 4 threads, 0 failing; 22 files with 218,046 sign-of-zero-only elements from `rho_d_maxloc_fix`, accepted by the rule). Non-vacuity on `dec2_obst` (2 species): `D_Z_MAX` 972 of 972, `DEL_RHO_D_DEL_Z` 512 of 1,944, `DP` 256 of 972, `RHO_D_DZD*` 1,090 / 196 / 1,090 of 1,944, `H_RHO_D_DZD*` 545 / 98 / 545 of 972. `csmag_32` is a weak case (uniform ZZ: all gradients zero, 7 of the 9 diffusion arrays all zero). Synthetic input counts only as device coverage. The Integration Lead's newer `device_gate/stage1_wp3_bitwise.log` (stage-1 commit 778034a) gives the gate-format line per kernel (arrays, elements, bit-equal, sign-only, diff, 330 input files compared by sha256 between the box and the test machine copy: 0 mismatches; the two lists differ only in order, I checked); the checker reads them as PASS.

Host leg (gate, latest runs under `vv-runs/gpu_gate/results`): `d_z_max` and `del_rho_d_del_z` have a test entry (`test/s5_r2_bitwise.F90`) but produced no result in the latest quick run (status FAIL, not NO-TEST: the run produced no result line for them; the coverage run lists them NOT-RUN); `rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat` are **NO-TEST** (no bitwise test against the verbatim FDS loop text exists; the checker flags them "no host test yet"). The stage-1 reference is the host run of the same generated source, with no FDS value.

Decisions:
- **`d_z_max`, `del_rho_d_del_z`: ACCEPT as device coverage**, with the note "inputs `RHO_D` and `D_Z` synthetic; device = serial = 4 threads bitwise on `csmag_32` and four `dec2_obst` boxes". Condition: the gate's host tier must show a result for both in the next full or quick run (a regenerate-and-run; nothing is missing in the test).
- **`rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat`: CHANGES.** (1) A bitwise host test in the gate (the generator owner adds the cases to `test/make_r2_tests.py`; the gate fails these kernels as NO-TEST until then). (2) A real-field run: dump `RHO_D` and the `D_Z` table of the driver (same scratch-tree dump route as the edge tables), rerun the `WP3_DIFF` stage on the test machine, and report device = serial = 4 threads. For `rho_d_interp` the point is the table lookup with the real table and real `TMP` including the end-of-table clamp (not shown for the synthetic table); for `rho_d_maxloc_fix` the branch with more than two species (`NS = 2` is the maximum run, so the species loop executes once) and a tie in `MAXV`; for `dp_div_heat` a non-synthetic `H_RHO_D_DZD*`.
- **`rho_d_dzd`, `h_rho_d_dzd`: ACCEPT** on the same evidence as 4a (the WP3 item 3 neighbours, nonzero on `dec2_obst`); host gate PASS (224 cases); also on record in gate set r1 (S4d harness, synthetic). The draft note "real fields; the device run is the one in the gate registry (set r1)" should say "ZZ and TMP real, `RHO_D` synthetic".
- The kernel map (`tools/device_runs.toml`) lists no device run for these five; PM-08 rejects entries without one. Either the overlay lists the stage-1 run, or the gate registry takes the log (`device_gate/stage1_wp3_bitwise.log` is in gate form).

### 3.5 Group 6: "WP1 non-wall and wall kernels": CHANGES
- The row names no kernels, so it cannot be mapped to `ported.toml`. List the kernels (the 23 non-wall and the 11 wall kernels), each with its truth status.
- 7.9 inputs are partly substitutes: `RHO_D` (synthetic), `DELTA_RHO_ZZ`, `FX_ZZ` etc. (copies of `FX`), `U_DOT_DEL_RHO_Z` (scratch value), the corrector stage uses the predictor flux, `Q` and `QR` are real zeros, `tg128` walls are reconstructed (SYNTH). `WP1_DIV1`, `WP1_MASS`, `WP1_CLIP` have no FDS truth. Kernels whose inputs are all real and that are bit-equal to FDS (`zzs_pred` ZZS, `kres`, `mu_dns` on fluid cells, `rhos` after `wall_uvw_interp`, 7.9a) are ready to be ACCEPTed once named.
- 7.10 wall kernels: real tables on `csmag_32` and `shunn3_32` with assumed `EW_NIC = 1`, `W_THIN = 0`, `W_SURF_INDEX = 0`; real on `dec2_obst` except `EW_NIC`; `wall_strain_rate` and `wall_us_pred` non-trivial only on `dec2_obst` (8 elements per box), with no FDS truth; the `NIC>1` branch is not run. After the shim, device = serial = 4 threads on 160 of 160 output arrays by sha256 (stricter than the rule). These are acceptable, but need to be named per kernel.
- **I could not re-compare the 7.9 device outputs:** under `results/test machine/res/gpu` only logs remain (raw outputs were removed for disk), and the wall results are `.bin` trees my script does not walk. Ask for the harness's own tables to be attached to the ported entry, or rerun on the test machine into a work folder if a bytes-level recheck is wanted.

## 4. `cfl_wall_max`: the gate's own ULP line

**Path used.** The stage-1 harness (7.10d) stores per-wall one-call results of both runners: `stage1/real/cflw_res/{cpu_serial,cpu_omp4,gpu}/<label>/uvw_wall.f64` (host `cflw_main.F90`, gfortran `-O2 -ffp-contract=off`; device `cflw_gpu.cpp`, nvfortran `-mp=gpu -gpu=cc89,nofma`, cc 8.9; one launch per box, exit 0). The device values were already on the development machine, so **no test machine run was needed**. I selected the `SOLID_BOUNDARY` walls (`BT == 1`, 8 per box; boxes 0 to 3 concatenated in order for the real case; all 20,000 for the synthetic one), wrote host and device raw dumps next to the log, and ran the gate's `bitcmp.py --ulp ... --limit 2 --emit`:

    ULP cfl_wall_max 32x1x1 maxulp 1 nvals 32 nspecial 0
    ULP cfl_wall_max 20000x1x1 maxulp 2 nvals 20000 nspecial 0

Files (gate repository, new): `vv-runs/gpu_gate/device_logs/stage1_cflw/device_cflw_ulp.log` (these two lines plus provenance comments), `host_real32.f64`, `dev_real32.f64`, `host_synth20000.f64`, `dev_synth20000.f64`. Reproduce: `python3 bitcmp.py --ulp device_logs/stage1_cflw/host_real32.f64 device_logs/stage1_cflw/dev_real32.f64 --limit 2 --emit cfl_wall_max 32x1x1`.

**`check_device_logs.py` on it:** log verdict PASS (2 passing lines, 0 failing); per kernel `cfl_wall_max  r5  libm  PASS  [max ulp 2]` (full output in `check_device_logs.out`; the other 113 registered kernels read MISSING there because only this log was given, as expected). No `out` name is attached, so both lines are judged against the per-value limit (`UVWMAX/ICFL/JCFL/KCFL` are dt-coupled and not judged by the gate).

**Cross-check.** The Integration Lead's `device_gate/stage1_cfl_wall_max_ulp.log` (stage-1 commit 778034a) carries the same values (`32x1x1 maxulp 1`, `20000x1x1 maxulp 2`, with `out UVW_WALL`), produced by their own script; mine come from the gate's tool on the stored dumps, so the two routes agree. Their log also has the full-call `UVWMAX` lines (584 walls: 0, 1, 1, 0 ulp for boxes 0 to 3; synthetic 1 ulp), which the gate shows but does not judge.

**Reading, and what it does not show.**
- 4 of 32 real values differ (box 1: 2, box 2: 1, box 3: 1), each by 1 ulp (7.1e-15 to 1.4e-14 absolute at values near 43 to 71); 28 of 32 are bit-equal. Matches 7.10d.
- The 32 real values contain only 19 distinct values (boxes 1 and 2 are mirror images) in the range 42.7 to 71.0, a narrow argument range. The synthetic sweep (SYNTHETIC, eleven decades) reaches exactly 2 ulp (635 of 20,000 values) and none at 3. The bound is a measured one for this toolchain (nvfortran 26.9 and its device libm against the host gfortran libm; the host compiler version is not recorded in the result files), not a guarantee; rerun the ULP lines after a compiler or libm change.
- The real case is a scratch hot-obstruction variant of `dec2_obst`, not a verification case; the value is the FDS-build expression, not the value inside `CHECK_STABILITY`.
- Run-level tolerance for `UVWMAX` (differs by 1 ulp in 2 of 4 boxes) is not judged by the gate and still needs an owner for the test (open point already in the gate README).

Verdict: the **ULP-gate condition is met** for the stage-1 evidence (`ulp_gate_passed = true` is justified from the V&V side). It does not make `cfl_wall_max` ported by itself (section 7).

## 5. Gate state with all logs on record

Run: `check_device_logs.py` on `device_cflw_ulp.log`, `stage1_cfl_wall_max_ulp.log`, `stage1_wp3_bitwise.log` (Integration Lead) and `bitwise_{dflt,dpd}_s96.txt` of the rounds 4 to 7 run: 25 of 114 registered kernels PASS, 0 FAIL, 89 MISSING (the other kernels have no log in this set; the tier stays red until they do). The 13 + 3 + 1 kernels of groups 1, 2, 3 and `cfl_wall_max`, and the five synthetic-input kernels, read PASS. `rho_d_dzd` and `h_rho_d_dzd` are in set r1 and read PASS (6 passing lines each) when the S4d logs `s5perf_{BIND,DPD}_gen1.log` and the `bitwise_*.txt` files of `prototypes/s4_cuda_mass/s5_gen_34/logs` are added (S4d logs alone: 34 of 114 registered kernels PASS). Note that these checker PASS lines say nothing on real versus synthetic input; that judgement is the one in this file.

## 6. Files written (gate side, `vv-runs/gpu_gate`, committed there)
`device_logs/stage1_cflw/` (ULP log, the four dumps, `check_device_logs.out`, `recheck_stage1_outputs.out`) and `recheck_stage1_outputs.py`. No tracked file of the gate was edited (the gate's `README.md`, `check_device_logs.py`, `gate.py` carry uncommitted edits of someone else and were left alone). The README line "the `ULP` line is not printed by any driver yet" is out of date; it is the gate owner's text to update.

## 7. Independent check of the D-072 (c) wording

Sources read: `README.md` D-072 (a) to (g), `inventory/port_kernel_status.md`, `tools/README.md` (PM-08, "Ported"), `amrex/ported-candidates.md`.

| Condition | D-072 (c) as ruled | Where else it is stated | Covered by this review? |
|---|---|---|---|
| 1. Device run on record on real fields, within the class tolerance (bitwise; 2 ulp for libm) | yes | `port_kernel_status.md`, `tools/README.md`, PM-08 ("a kernel with no device run on record") | **Yes**: adequacy of the device = host evidence, section 3; with the notes about synthetic and partial coverage. |
| 2. Reviewed by the kernel owner | yes ("the kernel owner has reviewed it") | PM-08 (`reviewed_by`, `set_by` present) | **No.** My review is the independent review of the device results, not the kernel owner's review. `ported.toml` draft uses `reviewed_by = "V&V Lead (reviewer of record ...)"`, which is a substitution of the owner by the reviewer of record; the Architect's co-signature must say that it accepts that substitution, or the owner's review is added. |
| 3. V&V ulp gate (libm kernel only) | **not in (c)**; (b) says the 2 ulp gate belongs to the V&V device-tier harness and the registry only records the tolerance | `port_kernel_status.md` ("a libm kernel can be ported once the V&V Lead's ulp gate ... exists and passes"), PM-08 (`ulp_gate_passed = true`) | **Yes** for `cfl_wall_max`: section 4. |
| 4. `ci_checks.sh --strict` passes | yes | `tools/README.md` | **No**, but checked: see below. |

So the "four conditions" are a correct combination of (c), (b) and PM-08, but D-072 (c) itself has three conditions; the ulp gate is a fourth only through (b) and PM-08. The candidates file's text is not wrong, but the ruling text should be quoted as three plus the libm rule.

`ci_checks.sh` run in this review (read only, output not kept in the trees): **FAIL**. `k2_ci_check` 2 FAIL and `kernel_lint` 1 FAIL (CL-07), all on the radiation golden `s5gen_rad_wall_qin_zero` (`W_SURF_INDEX` missing from the device-address list; `rad_wall_qin_zero` is not a candidate); `zone_sum_order` PASS; `port_kernel_map` PASS (114 kernels, self-check PASS). Warnings: BF-06 managed memory in two device build scripts (bring-up aid; an acceptance run must not use it: the rounds 4 to 7 hash runs are built that way, so they stay coverage evidence only). The same finding as in the candidates file. As D-072 (c) says "passes `ci_checks --strict`" without a per-kernel scope, **no kernel can be set ported until the script is green**, or the Architect rules that the condition applies to the findings that concern the kernel itself.

Two more points on the draft `tools/ported.toml` (uncommitted change of the docs working tree): (1) it already names the V&V Lead as `reviewed_by` for four groups, before this review; it must be revised to this review (drop or hold 4b and group 6, correct the group 1 note, "ZZ and TMP real, `RHO_D` synthetic" for group 5) before it is committed. (2) It does not contain `cfl_wall_max`, consistent with section 4.

## 8. What is needed to close the CHANGES

1. Group 4b: bitwise host tests for `rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat` in the gate (generator owner); real `RHO_D` and `D_Z` dump and a test machine rerun of `WP3_DIFF` (Integration Lead). `d_z_max` and `del_rho_d_del_z`: a host-tier result in the next gate run.
2. Group 6: kernel names with truth status; ideally the 7.9 device outputs kept or rerun for a bytes-level recheck.
3. Group 1: correct the note (device on `csmag_32`; baroclinic also on `dec2_obst`).
4. Group 3 (optional): real-field rerun with `CFL_VELOCITY_NORM` 1 to 3, `CYLINDRICAL = 1`, `STORE_CARTESIAN_DIVERGENCE = 1` (device = host only).
5. Condition 2 (kernel-owner review) and condition 4 (`ci_checks --strict` green) before any entry is committed; Architect's ruling on the reviewer substitution.

## 9. The six whole-field `gsfv_*` kernels (follow-up; the Architect's ruling D-077 names the V&V Lead reviewer of record in place of the kernel owner)

Kernels: `gsfv_central`, `gsfv_godunov`, `gsfv_superbee`, `gsfv_minmod`, `gsfv_charm`, `gsfv_mp5` (whole-field `GET_SCALAR_FACE_VALUE`, one kernel per limiter).

### 9.1 Correction of the candidates file
`amrex/ported-candidates.md` (v0.1) lists them as "not candidates (no device run)". That is wrong, and was stale when written: a device run of all six is on record in two forms, and both pass.
- **Form A, the gate's device driver** (`gpu/s5_dev_hash.F90`, hash form): test machine runs on 2026-10-03, logs in `(local workspace)/scratch/devtier/results/run1..run3/logs` (not under `vv-runs/gpu_gate/device_logs`, which holds only the `cfl_wall_max` material of section 4; the gate does not look for logs in a fixed place, it is given paths). Grids (`NX x NY x NZ`, one cell block, scalar `UU(-1:NX+2,..)`, face velocity `A(0:NX+1,..)`): 24x16x12, 96x40x56, 128x128x128 (two runs), 192x192x192; two directive builds (`dflt`, `dpd`), plus a host serial and a host all-threads run per grid. Inputs: SYNTHETIC, generated by the integer random generator of the gate tests (exact operations only; hash of the inputs is equal on host and device in every run). Three input classes, each with all six limiters and `IOR` = 1, 2, 3: class 0 (A > 0 only, scalar random in [-1, 1]), class 1 (A of mixed sign including +0, -0, +-1e-300 and +-(2 x random), scalar scaled by 2^-12 .. 2^12), class 2 (A random sign, scalar integer-valued in -2..2, so many exact ties). What is compared: a 128-bit running hash over the 64-bit pattern of EVERY element of the whole output array `F` (including the cells outside the box, preset to a sentinel), device generated kernel (nvfortran `-O2 -mp=gpu -gpu=cc89,nofma,mem:managed -Minline`) against the VERBATIM upstream `GET_SCALAR_FACE_VALUE` text compiled with gfortran (reference compiled without FMA). Hashing the raw pattern is stricter than the +0/-0 rule (it would flag a sign-of-zero difference; none occurred). `cmp_hash.py`: all 6 kernels BITWISE EQUAL, host gen == ref True, inputs identical True, in 8 + 2 + 2 files (`bitwise_{dflt,dpd}_*`). Run 2 also contains the launch-order mutant build (only meaningful for the `pvf_*` kernels).
- **Form B, stage-1 harness on real FDS calls** (plan section 7.13c; stage-1 repository `gsfv/`, commit 9715a9e; results `gsfv/results/`): the inputs `A`, `U`, `F` (before and after) of every whole-field `GET_SCALAR_FACE_VALUE` call, taken from a scratch copy of the driver tree. 1,296 dumps: `shunn3_32` (32x1x32, 2 species, 42 dumps per limiter, A > 0 only), `dec2_obst` (four 16x1x16 boxes with an obstruction, 168 per limiter, A of both signs: 2,072 negative cells per limiter run), `csmag_32` (32^3, 6 per limiter, **uniform scalar**), each run once per limiter (the six limiters), so every kernel has an FDS result on its own limiter. Every kernel ran on every dump: 7,776 kernel runs, device (one job, exit 0, cc 8.9 test-machine GPU, nvfortran `-mp=gpu -gpu=cc89,nofma`), host serial and host 4 threads (gfortran `-O2 -ffp-contract=off`). Result: the device `F` is byte-identical to host serial and host 4 threads in 7,776 of 7,776, and for the kernel whose limiter matches the FDS limiter of the dump, 1,920,768 of 1,920,768 box cells are bit-equal to FDS for each of the three runners (cells outside the box unchanged from `F` before, as in FDS). Pinned source: generator branch `s5-gen` commit 3546f746a7 (committed text; the six kernels in the working tree are text-identical).

### 9.2 What I checked myself
1. **Gate checker** (`check_device_logs.py` on the Form A logs): all six `gsfv_*` kernels read PASS (24 passing lines each, 0 failing), set `gsfv_field`. The gate's own README row for this set already said "yes"; only the `SETS` label inside the checker still says "no device run on record" (stale, same type of fix as the README line in section 6).
2. **Form B logs re-read** (`gsfv_full.log.gz`, `dev_sha.txt.gz`, `dev_run.log.gz`): 7,776 kernel lines, 7,776 with `device==serial True serial==4thr True device==4thr True`, 0 DIFF; the device sha256 list has 7,776 entries (6 outputs x 1,296 dumps); device job log ends `gsfv_gpu done`. Consistency checks on the device sha list that I computed: on all 36 `csmag_32` dumps the first five limiters give the same output hash and `gsfv_mp5` differs from them (as the stage-1 text states for a uniform scalar), and on every `shunn3_32` and `dec2_obst` dump the six limiters give 3 to 6 distinct outputs, so the data discriminates the limiters and the comparison is not vacuous. Negative controls are in the stage-1 results (another limiter against the FDS result agrees on 0 to about 37,000 of 58,240 / 59,136 cells).
3. **Not possible:** a re-comparison of the raw device `F` arrays with `bitcmp.py`. The dumps and the host outputs were deleted after the results were recorded (disk), the device files stay on the test machine, and only the sha256 lists came back. sha256 over whole files compares all bytes and is stricter than the +0/-0 rule, and the stage-1 compare reports no sign-of-zero case. A rerun would need the scratch driver dumps again (about 2 GB for `tg128`; the three cases used here are smaller) and an idle machine; not done. This is an evidence limit, not a defect: the gate-form Form A run covers the same comparison with the gate's own driver.
4. **Gate host steps** (`vv-runs/gpu_gate/results`): `bitwise_gsfv_point` and `bitwise_gsfv_field` PASS in every run on record, including the latest quick run (2 flag sets) and the last full run (6 flag sets, threads 1/4/8 inside the program); `drift_gsfv_field` / `drift_gsfv_point` PASS in the last full run and in the earlier quick runs, and TIMED OUT (exit 124 after the 900 s limit) in the latest quick run only. **I did not rerun them:** at the time of this review the machine load average is about 57 on 8 cores and another gate run is in progress, so a rerun would time out again and prove nothing. The timeout is a load effect (the same steps take 20 to 80 s on an idle machine, golden `test/gsfv_field_kernels.golden` last changed on 2026-10-02, before the PASS runs), not a content failure. Rerun when the load is low: `run_gpu_gate.sh --tier quick --only drift_gsfv_field,drift_gsfv_point` (expected PASS); I would not enter the kernels into the map as "final" without that line, but it does not change the device-result judgement.
5. **Build flags and arch:** device `-gpu=cc89,nofma` (cc 8.9), `nvcc --fmad=false` for the harness; host `gfortran -O2 -ffp-contract=off`; no fast-math option in either build file (read). FDS reference build `-O3 -std=f2018 -frecursive -fopenmp` with no FMA instruction in its binary (plan 7.9c). Managed memory (`mem:managed`) is used by the Form A build only (bring-up aid; Form B does not use it). `gsfv_mp5` is the most arithmetic-heavy kernel and the first to expose a contraction difference; it is bit-equal in both forms with `nofma`, so the no-FMA setting is part of the evidence scope.

### 9.3 Decision per kernel

| Kernel | Decision | Reason (one line) |
|---|---|---|
| `gsfv_central` | **ACCEPT** | Real calls (3 cases) bit-equal to FDS, device = serial = 4 threads; gate-form run bitwise on 4 grid sizes (5 runs) x 2 builds. |
| `gsfv_godunov` | **ACCEPT** | Same; upwind choice exercised for both signs of A by `dec2_obst` and by the synthetic classes 1 and 2. |
| `gsfv_superbee` | **ACCEPT** | Same; limiter ratio branches exercised by real varying scalars and by the tie class. |
| `gsfv_minmod` | **ACCEPT** | Same. |
| `gsfv_charm` | **ACCEPT** | Same (CHARM is the default limiter of `shunn3_32` and `dec2_obst`, so its evidence is the largest). |
| `gsfv_mp5` | **ACCEPT** | Same; on the uniform-scalar `csmag_32` it differs from the other limiters by rounding and still equals FDS when FDS ran MP5; on the varying cases it equals FDS bitwise. |

All six: ACCEPT as device evidence. The `drift_gsfv_*` rerun (9.2 item 4) is the one open gate-side line, and the two project-wide conditions of section 7 (`ci_checks --strict` green; D-077 now covers the reviewer substitution) are unchanged. The `ci_checks.sh` result at the time of this review is the same as in section 7 (k2_ci_check and kernel_lint FAIL on `rad_wall_qin_zero`; zone_sum_order and port_kernel_map PASS).

### 9.4 Coverage limits (carry these in the note)
- **Negative `A` with a varying scalar** on real data only in `dec2_obst` (2,072 cells per limiter run); `shunn3_32` has A >= 0 only and `csmag_32` has a uniform scalar (no gradient: a test of plumbing and of MP5 rounding only). The synthetic Form A classes 1 and 2 supply negative and zero `A`, +0 and -0, tiny values and ties on 4 grid sizes.
- **`A == 0` cells** take the `A <= 0` branch (all y-faces of the 2-D cases).
- **Signed zeros and near-ties** are covered by the synthetic class only (real data does not control them).
- **One rank, one level**; no `tg128` run (uniform scalar; no new information on the limiters).
- **Interface condition, not a result defect:** FDS passes the face velocity with lower bound -1 in its own direction, the kernel declares `A(0:IBAR+1,..)`. The harness copies index by index (kernel `A(0,0,0)` = FDS `A(0,0,0)`) and passes. A caller that passes the START of the FDS array, or of a MultiFab fab beginning at -1, is shifted by one cell silently. The device driver must pass the address of element (0,0,0), or the generator must declare the lower bound -1; this has to be fixed or tested at the first driver call site (request recorded in plan 7.13c).
- The device evidence is on the stage-1 and gate-driver call interface, not yet inside the driver; the scalar-face hooks (flux-hook review) and the wall-seam order are separate work.

### 9.5 Exact `ported.toml` note wording (evidence scope kept in the note)
One `[[ported]]` entry for the six kernels (same evidence), `reviewed_by` as for the other entries, with a note of this form (the Integration Lead owns the file; text below is the scope, not an edit):

    kernels = ["gsfv_central", "gsfv_godunov", "gsfv_superbee", "gsfv_minmod", "gsfv_charm", "gsfv_mp5"]
    note = "Bitwise class. Device (cc 8.9, nvfortran nofma) = host serial = host 4 threads, byte-identical F, on 7,776 kernel runs over 1,296 REAL GET_SCALAR_FACE_VALUE calls (shunn3_32, dec2_obst four boxes, csmag_32; stage-1 gsfv run, stage1 repo 9715a9e), and bit-equal to FDS on 1,920,768 of 1,920,768 box cells for the kernel of the matching limiter. Gate device driver (SYNTHETIC random inputs incl. signed zeros and ties, grids 24x16x12 to 192^3, two directive builds): hash of the whole output array equal to the verbatim upstream text compiled with gfortran, 0 DIFF. Limits: negative A with a varying scalar only in dec2_obst (2,072 cells); csmag_32 has a uniform scalar (plumbing and MP5 rounding only); one rank, one level; the caller must pass the address of element (0,0,0) of the FDS face velocity (lower bound -1 in the flow direction). Raw device arrays are not kept (sha256 and hash comparison only); gate drift_gsfv_* steps timed out once under machine load (rerun pending)."

If the Integration Lead prefers one entry per kernel, repeat the same note; the evidence is identical for the six. Remove the line "gate drift_gsfv_* steps timed out ... (rerun pending)" once the rerun passes.

### 9.6 Corrections to carry back
1. `amrex/ported-candidates.md` section 2, "Not candidates": the six `gsfv_*` kernels have device runs (Forms A and B above); they are now candidates, ACCEPT.
2. `check_device_logs.py` `SETS` label for `gsfv_field` still says "NO device run on record" and the gate README "Device-tier gaps" paragraph does not mention the `ULP` line now on record (section 6): gate-owner text to refresh.
3. `tools/device_runs.toml` already lists the stage-1 gsfv run (stage1 repo 9715a9e); the map shows a device run for all six.
