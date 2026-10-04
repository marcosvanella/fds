# Run-level class tolerance for `dt_coupled` kernels (D-070, D-078 a): measured proposal

Owner: V&V Lead. Status: v0.1 (proposal; `run_class_tolerance` in `vv-runs/gpu_gate/kernel_categories.json` stays `null` until the owner decides). Date of the measurement: 2026-10-04.
Subject: the number that `check_device_logs.py` compares with `X` of the `RUNCMP` line (gate commit 38e94ae), `X = max over steps of |dt_dev - dt_host| / |dt_host|`, for kernels of class `dt_coupled` (today `cfl_wall_max`; `wall_visc_les` when registered). `c75997b` and test plan 5.10 say only "the run-level tolerance already used for cross-compiler runs", with no number.

## 1. Proposal

**`run_class_tolerance = 1e-11`** (relative, on every step of the compared run, last step included), for compared runs of up to 1,000 steps. For longer runs, the same model is `1e-14 x N` (N = number of steps compared); the gate has no N-dependent rule today, so a run longer than 1,000 steps needs the owner's decision first.

Basis (sections 3 to 5):
- The measured Intel-against-GNU `X` over the cases that agree is at most **4.3e-13** (all steps, including the step shortened to reach `T_END`) and **8.4e-14** when that last step is left out. 1e-11 is 23 times the first and 120 times the second.
- No baseline case lies between 4.3e-13 and 2.6e-4. The cases whose step path diverges start at 2.6e-4 (smallest first deviation seen). Any value between 1e-12 and 1e-4 gives the same verdict on every case on record; 1e-11 is 1.4 orders of magnitude above the largest agreeing case and 7.4 orders below the smallest deviation of a diverging case (2.6e-4).
- The tolerance already used for cross-compiler runs is the T1 rule, 1e-10 of the column inf-norm (`requirements.md` T1; `intel_vs_gnu.py` class ROUNDOFF). That is a field-and-column rule, not a time-step rule, but it is the nearest existing number. Every agreeing case also satisfies 1e-10 (by a factor of 230 or more), so adopting 1e-10 instead would change no verdict on the data; it is 10 times looser than the measurement supports. If the owner prefers to keep one number across both rules, 1e-10 is acceptable; I propose 1e-11 because it is the number the step-sequence data themselves support.
- Linear model: 1e-14 per step is 3 times the largest per-step rate seen in any agreeing case (3.5e-15 per step, set by the last step) and about 10 times the typical rate (1e-15 per step, section 4). At 1,000 steps it gives 1e-11.

## 2. What was measured

Data: `vv-runs/baseline/gnu_ompi_firex-36975d7` (GNU Release, Open MPI) and `vv-runs/baseline/impi_intel_firex-36975d7` (Intel Release, Intel MPI, including the 17 runs made on the second machine), both FireX 36975d7, 110 runs present in both; the class of each comes from `vv-runs/analysis/impi_intel_firex-36975d7/intel_vs_gnu.json`. Read only; no simulation was run.

**Where a per-step time sequence exists.** FDS does not write DT per step at full precision:
- `_steps.csv` and the `.out` step lines give the step size with 3 digits (`0.446E-01`) and the total time with 7: too coarse for round-off (resolution about 1e-3 for the step size).
- The `_hrr.csv` time column is single precision (its values are float32 numbers printed with 17 digits), so it cannot resolve double-precision differences; it is also written at fixed intervals once a run has more than about 1,000 steps.
- The `_devc.csv` time column is double precision and has one row per time step when `DT_DEVC` is smaller than the step (the `SIG_FIGS=17` copies set `DT_DEVC=1e-5`). Rows equal steps + 1 for runs up to about 1,000 steps. This is the record used: `dt_n = T_n - T_(n-1)` from consecutive rows, with `X = max_n |dt_Intel - dt_GNU| / dt_GNU`.
- Alignment check: row n must equal the total time printed on the `.out` step line n (to the printed digits) for every printed step, on both builds; otherwise the file is rejected as not per step. Mass files were accepted only on the same test.

Which of the 61 `SIG_FIGS=17` runs qualified (the as-committed copies print 8 digits and cannot show round-off; they are used only in section 5):

| Result of the selection | Runs |
|---|---|
| Per-step double-precision times, same step count, class IDENTICAL / ROUNDOFF / SMALL ("agree") | **11** |
| Per-step double-precision times, same step count, class LARGER (step path diverges; context only) | 6 (4 distinct cases; the `_np1` and `_glmat` copies repeat them) |
| Time column is float32 (`saad_512_cfl_*` x4, `shunn3_*` x4): `X = 0` here only means equal at float32 resolution (about 6e-8 of T) | 8 (not used in the statistics) |
| Time column not one row per step (more than about 1,000 steps, or an output interval) | 31 |
| No time column (`lapse_rate`, 4 steps) | 3 |
| Step count differs between Intel and GNU | 2 |

The 11 agreeing cases: `csmag_32`, `csmag_32` with `FISHPAK_BC` periodic, `divergence_test_2`, `ns2d_8`, `ns2d_8_nupt1`, `ns2d_16`, `ns2d_16_nupt1`, `ns2d_32`, `ns2d_32_nupt1`, `symmetry_test_mpi` and its single-mesh copy (the 8-rank `symmetry_test_mpi` ran on the second machine; 1,000 steps). Steps per case 28 to 1,000. Step counts are equal on both builds in every case used.

## 3. Result: the cases that agree

`X` per case (`SIG_FIGS=17` copy). "X excl. last" leaves out the final step, which FDS shortens to land on `T_END`. "q" is the quantisation of the record itself: `ulp(T_n) / dt_n`, the smallest relative dt difference that two time records can show at that step; `max X/q` says how far above the record's resolution the result lies.

| Case | Steps | X | X excl. last | max over steps of X/q | X in steps 1 to 10 | X in steps 1 to 100 |
|---|---|---|---|---|---|---|
| `csmag_32` | 28 | 7.9e-15 | 2.2e-15 | 2 | 8.0e-16 | 7.9e-15 |
| `csmag_32` (`FISHPAK_BC` periodic) | 28 | 0 | 0 | 0 | 0 | 0 |
| `divergence_test_2` | 57 | 0 | 0 | 0 | 0 | 0 |
| `ns2d_8` | 156 | 2.3e-14 | 2.3e-14 | 1 | 1.3e-15 | 1.2e-14 |
| `ns2d_8_nupt1` | 120 | 4.3e-13 | 1.3e-14 | 4 | 0 | 9.1e-15 |
| `ns2d_16` | 295 | 1.4e-13 | 4.4e-14 | 3 | 6.8e-16 | 1.0e-14 |
| `ns2d_16_nupt1` | 227 | 2.3e-13 | 2.8e-14 | 6 | 6.2e-16 | 1.6e-14 |
| `ns2d_32` | 568 | 4.3e-13 | 8.4e-14 | 4 | 0 | 1.1e-14 |
| `ns2d_32_nupt1` | 440 | 2.7e-13 | 5.3e-14 | 5 | 1.3e-15 | 4.6e-15 |
| `symmetry_test_mpi` | 1000 | 0 | 0 | 0 | 0 | 0 |
| `symmetry_test_mpi` (1 mesh) | 1000 | 0 | 0 | 0 | 0 | 0 |

**Distribution of X over the 11 cases:** maximum **4.3e-13** (`ns2d_8_nupt1` and `ns2d_32`), median **2.3e-14**; 4 cases are exactly 0 (the whole time sequence is bitwise equal); the 7 non-zero cases have median 2.3e-13. Leaving out the last step: maximum **8.4e-14**, median 1.3e-14.

Readings:
1. **The maximum of every non-zero case falls on the last step** (6 of 7; the exception is `ns2d_8`, step 119 of 156, where it equals the record's resolution). The last step is the one shortened to reach `T_END`; its `dt` is a remainder `T_END - T`, so one ulp of accumulated T difference becomes `ulp(T)/dt_last` relative. It is an artifact of the clamp, not a dt-law difference. A device driver that compares the final clamped step inherits the same amplification; this is why the proposal keeps a margin of 23 and not 2.
2. **The result is at the resolution of the record.** `max X/q` is 6 or less in every case; the cumulative time difference `|T_Intel - T_GNU|` stays within 6 ulp(T) (relative 9e-16 at most) over runs of up to 568 steps. The measured `X` is therefore an upper bound limited by how finely two time records can show a dt difference, not a measured dt difference. It does not reveal whether the true per-step difference is 1e-16 or 1e-14.
3. **Indicative size of the true difference.** In the non-zero cases the time sums differ in the last bit at 2.5 to 8 % of the steps (7 to 17 of 120 to 568; `csmag_32`: 41 %). If the true relative dt difference is `delta`, the rounding of `T + dt` flips with probability about `delta * dt / ulp(T)`, so `delta` is about (fraction of steps that flip) x (mean q). That gives **about 1e-15 relative per step** (5e-16 to 1.2e-15 over the 6 non-zero double-precision cases other than `csmag_32`, and 7e-16 for it), constant in step index. This estimate assumes uniform rounding and is only an order of magnitude. The fields of these cases differ by 1e-14 of the column norm in the same runs (class ROUNDOFF), which is larger than a 2-ulp effect.

**Growth against step index** (largest relative dt difference in any agreeing case, last step left out):

| Steps | Largest X | Cases |
|---|---|---|
| 1 to 10 | 1.3e-15 | 11 |
| 11 to 100 | 1.6e-14 | 11 |
| 101 to 300 | 4.4e-14 | 8 |
| 301 to 1000 | 8.4e-14 | 4 (2 of them exactly 0) |

The growth is that of the quantisation floor, which rises with T because `ulp(T)` does and `dt` stays about the same (the floor `q` runs from 4e-15 in `csmag_32` to 1.1e-13 for 1,000 steps in `symmetry_test_mpi`). Dividing by q, the ratio stays between 0 and 6 at all step indices. No exponential growth is visible, and there is no sign of chaotic amplification in these cases. They are laminar or low-Re 2-D cases and a short 3-D LES case; a chaotic fire case with equal step counts and agreement is not available (the chaotic cases are in the diverging set).

## 4. Result: the cases whose step path diverges (context, not used for the tolerance)

Same record, class LARGER, equal step count, `SIG_FIGS=17` copy:

| Case | Steps | First step with X above 1e-12 | X at that step | X in steps 1 to 10 | X in steps 1 to 100 | X over the run |
|---|---|---|---|---|---|---|
| `obst_activation_default` (and `_glmat`) | 144 | 9 | 1.9e-3 | 1.9e-3 | 1.0e-1 | 0.64 |
| `obst_coarse_fine_interface` (also 1 rank) | 43 | 3 | 5.6e-4 | 1.5e-3 | | 0.15 |
| `species_conservation_2` | 187 | 1 | 2.6e-4 | 2.6e-4 | 1.0e-1 | 0.99 |

Pairs with a different step count (`obst_activation_ulmat`, the restart set, `1_step_2_step_compare`, `layer_1mesh`, `ns2d_16_emb_1to2`, `species_conservation_1`) are in the classification of `baseline_status.md` and have no common step sequence.

Reading: these runs do **not** show gradual amplification of a round-off difference. `obst_activation_default` is bitwise equal for 8 steps (X = 0), then jumps to 1.9e-3 in one step (the step 9 size differs, 0.03848 against 0.03855); `obst_coarse_fine_interface` jumps at step 3; `species_conservation_2` differs by 2.6e-4 from step 1. A deviation of 1e-4 or more appearing in one step is a difference in a discrete decision or in a property evaluation, not round-off growth. I did not diagnose their cause here (`baseline_status.md` records them as LARGER, with the interface and species cases as known differences). This matters for the tolerance only as the upper side of the gap: the smallest deviation seen in a diverging run is 2.6e-4, 7 orders above the tolerance.

## 5. Coarse check on the whole agree set (as-committed copies, `.out` step lines)

For every as-committed case of class IDENTICAL, ROUNDOFF or SMALL with equal step counts (32 cases, including `dancing_eddies_default` and its variants, SMALL, 4,737 to 4,751 steps, 4 ranks, second machine), the total time and the 3-digit step size printed on every printed step line (12 to 66 lines per case, up to step 4,751) are **equal on Intel and GNU in every line** (relative difference 0 at 7 and 3 digits). This only bounds the cumulative difference by 1e-7 of T; but it shows that the step sequence of the SMALL class cases does not drift apart over thousands of steps. In most of the diverging cases from section 4 these lines differ (relative 8e-6 to 5e-2 in the total time).

## 6. What the tolerance can and cannot prove for a device kernel that differs by 2 ulp per value

A kernel of class `libm` is allowed to differ by up to 2 ulp per value (about 4.4e-16 relative). When such a kernel feeds `dt` (`cfl_wall_max` feeds `UVWMAX`), the first-order change of dt in one step is of that size or smaller (a maximum over walls changes by at most the ulp error of the controlling value): **about 1e-16 to 4e-16 relative at the first step**, as noted in the request. The measured Intel-GNU differences (about 1e-15 per step, resolution-limited bound 4e-13) are in the same range or above, on fields that differ by 1e-14, so the CPU-CPU data are a reasonable but not equal analogue: they cover compiler and libm differences in the whole code, not one kernel on a device.

**It can prove:**
- That in the compared run the device step sequence follows the host sequence to within 1e-11 at every step: the same step count, no time step set by a different limiter, no shifted event time. A gross error in the dt-feeding quantity (relative 1e-11 or more, including any wrong branch or wrong constant) is caught; every diverging case on record is caught at its first deviation (2.6e-4 or more, 7 orders above the tolerance).
- That a 2-ulp per-value difference does not, over the compared steps, grow into a step-path change in these cases (the cases on record show none up to 1,000 steps).

**It cannot prove:**
- **That the kernel is within 2 ulp.** That is the `ULP` line, judged value by value. `RUNCMP` is a second, independent line; passing it does not substitute for the per-value rule, and 1e-11 is 4 to 5 orders of magnitude looser than the 2-ulp effect it accompanies. A kernel with an error of 1e-12 relative (2,000 times 2 ulp) in a dt-feeding value passes `RUNCMP` and fails `ULP`.
- **That other outputs agree.** dt is a minimum or maximum over cells (a norm), so it is blind to errors in cells that do not control it. Field-level agreement is judged by the bitwise or ULP lines of the other kernels and by the case tolerances (T1, T2).
- **That the resolution of the record is enough to show 1e-16.** In the baseline data the smallest difference that can be seen is about 1e-15 to 1e-13 (section 3). A device driver that prints `dt` directly, not time differences, does not have this floor; for such a driver the values of `X` will be smaller (expected about 1e-16 to 1e-14 from the estimate above) and the tolerance will have more headroom than the baseline data show. The baseline gives an upper limit for what a round-off-level difference produces, not the value the device will give.
- **Anything beyond about 1,000 steps** (the longest agreeing record) or about 570 steps with a non-zero difference; and anything in a chaotic case. Chaotic runs (fire, LES with strong feedback) amplify any round-off difference until the step sequences separate, with Intel and GNU as the example (section 4: `layer_1mesh`, the restart set); a device run of such a case would not give `X` below the tolerance for reasons unrelated to the kernel. `RUNCMP` cases should be chosen from the laminar or short runs of section 3 and not from the chaotic set.
- **That the clamped last step is compared fairly**: see section 3 reading 1.

## 7. Notes for the gate and the driver

1. `RUNCMP ... steps N maxrel X`: use the time steps the solver takes (the `dt` values as the solver computes them), not differences of an accumulated T; if a differenced T is the only record, the floor of section 3 applies (up to 1e-13 for T about 6).
2. Choose the compared run so that the kernel is exercised in the controlling position of the dt minimum for a good fraction of the steps (for `cfl_wall_max`: a wall-dominated `UVWMAX`), otherwise `X = 0` is trivially met. The baseline cases with `X = 0` (4 of 11) show what an inert comparison looks like.
3. The 4 cases with float32 time columns and `X = 0` do not enter the statistics; they are listed in section 2 so that nobody reads them as a bitwise result.
4. Files: the analysis scripts, the per-case table and the summary are in `vv-runs/analysis/dt_class_tolerance/` (`load.py`, `ana2.py`, `stats.py`, `stats2.py`, `per_case_sf17.txt`, `summary.txt`). `kernel_categories.json` and the gate were not edited.

## 8. Limits of this measurement

- Eleven cases, 7 of them non-zero, all CPU-CPU, all GNU against Intel (Intel Release is `-O2`, GNU `-O3`; across-compiler results are not bitwise comparable). The 2-ulp libm effect on a device is not measured by any of them.
- The record resolves dt only to 1e-15 to 1e-13 (section 3), so the numbers are upper bounds.
- Only runs up to about 1,000 steps have a per-step record; the longer SMALL cases (4,700 steps) are bounded only to 1e-7 of T at printed steps (section 5).
- No chaotic case with equal step counts and agreement exists on record, so the tolerance says nothing about growth in chaotic runs.
- The explanation of the three diverging cases (section 4) is not established here.
