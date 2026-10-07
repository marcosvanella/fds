# 04 — D-051 device sign-off: `rho_d_maxloc_fix` (D-080 (2))

Owner: AMR Species & Combustion Lead · Status: ACCEPT-with-conditions · Kernel: `rho_d_maxloc_fix` (divg.f90 DIVERGENCE_PART_1 MAXLOC species-sum fix).

## Verdict

**ACCEPT-with-conditions** for the device tier of `rho_d_maxloc_fix` under D-080 (2).

Conditions:
1. Real-data exact-tie frequency for `N_TRACKED_SPECIES >= 3` remains **PENDING** (A-76). Available dumps only cover NS = 1 or 2 (see below). Unblock: a `WP1_DIV1/ZZP.f64` (or equivalent stage dump with `meta.json`) from a case with NS >= 3, without writing under `vv-runs/refbin`.
2. NaN in `ZZP` stays outside the bitwise contract (README Limit).

## Evidence

| Item | Result |
|------|--------|
| Tie rule / NaN Limit text | s5-gen commit `9f0b70924b` (README "Tie rule" and "Limit" for `rho_d_maxloc_fix`) |
| Head under test | s5-gen `ce06e0c8dd` (private clone `work-species/ci_head`) |
| Lock | `test/with_s5lock.sh` (flock + timeout 1800) |
| `run_div1` O0 dpd/bind | PASS, 7668 cases each (rho_d_maxloc_fix 2160, failing 0, vacuous 0) |
| `run_div1` O2omp dpd/bind @ 4 and 8 threads | PASS, 7668 cases each |
| `mutation_div1` | 32/32 caught |
| Mutants named in D-080 (2) | `maxloc_x_last_wins` caught (2160 cases, 1365 failing); `maxloc_z_last_wins` caught (2160 cases, 1365 failing) |
| Tie test | `test/s5_div1_bitwise.F90` routine `t_maxloc` (NS = 1, 2, 3, 4, 8; all-equal, one-hot, last-species tie with earlier) |

Logs: `work-species/logs_s4/r6/run_div1.log`, `work-species/logs_s4/r6/mutation_div1.log`.

## First-occurrence tie rule

Upstream `MAXLOC(..., 1)` keeps the first species on an exact tie of the two-cell species sum. The generated kernel keeps the first with a strict `>` (species loop from 1, then 2..NS). The two last-wins mutants (`> MAXV` -> `>= MAXV` on the x and z pairs) are caught by the tie scenarios in `t_maxloc`.

## NaN outside contract

A NaN in `ZZP` is outside the bitwise claim: serial `MAXLOC` and the kernel's strict `>` disagree on NaN, and MAX/MIN behaviour differs across optimisation levels and compilers. Not tested. Stated in the generator README Limit (`9f0b70924b`).

## Real-data tie frequency (A-76)

Measured on existing `WP1_DIV1/ZZP.f64` dumps under `(local prototypes directory)/s4_cuda_mass/stage1/` (read-only; `vv-runs/refbin` not touched), using `work-species/maxloc_ties/count_ties.py` (IEEE float64 sum of the two face-neighbour cells, count of pairs where two or more species share the exact maximum):

| Case | NS | Face pairs | Exact ties | All-equal | NaN |
|------|----|------------|------------|-----------|-----|
| shunn3_32 (data_wp1 + run) | 2 | 2 x 6534 | 0 | 0 | 0 |
| csmag_32 (data_wp1 + run) | 1 | 2 x 107811 | 0 | 0 | 0 |
| dec2_obst_b0..b3 (real/data_pf) | 2 | 4 x 1734 | 0 | 0 | 0 |
| **TOTAL** | | **235626** | **0 (0 %)** | **0** | **0** |

NS = 1 cannot form a multi-species tie. NS = 2 cases show no exact ties in these dumps. **PENDING:** the same count on a dump with NS >= 3 (where first-occurrence among three or more species is exercised on real fields). Log: `work-species/logs_s4/tie_freq_wp1_div1.txt`.

## Notify

Integration Lead and Architect: this file is the record (do not duplicate).
