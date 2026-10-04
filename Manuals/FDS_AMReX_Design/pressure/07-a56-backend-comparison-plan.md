# 07. A-56 backend comparison plan (MLMG vs assembled HYPRE vs FFT) and first results

**Status: plan complete; first accuracy set measured; the R=64 cause established (section 8); repeat-timing tables in sections 7.4 and 7.5 are filled as the runs finish (marked PENDING until then).**
Measured = run. Read = source only. Unverified = neither.
Software: sections 1 to 10 (except where stated) use the common-layer tip of the pressure backend (tree this repository, commit `1f3f4de4c4`; latest `Source/pressure_backend` commit `29c2f9aa4a`, which contains the D-067 defaults from `a4a1073952`), read only and built unchanged. Sections 11 and 12 use the later tip `626c5f4b7f` (adds mixed open/closed faces, `fold_boundary_data`, the D-057 mapping check), built unchanged in a scratch directory. AMReX 26.09 and HYPRE v2.32.0-24 for that tip; the assembled-HYPRE and stretched/masked harnesses (docs 05, 06) use AMReX 99ddfda and the same HYPRE.
Machines: CPU timing and all GPU runs on an owner-provided NVIDIA test machine (pinned performance cores, 1 thread per rank, one GPU, managed arena). CPU accuracy checks of the tip ran on the shared build box.

## 1. Answers in short

1. **Plan:** A-56 is checked per zone by one number, `eps_H = max(1e-8, 2.4e-12·N²)` (requirements §2.1a), on a single frozen solve on the same discretisation, after removing the volume-weighted mean per zone, relative L2 (section 2). Cases, sizes, rank counts and pass/fail rules are in sections 3 to 6.
2. **What is measurable today:** the product backend (common layer + MLMG + FFT) can be compared MLMG-vs-FFT on the uniform box only. Assembled HYPRE is not in the product backend; it exists only in the study harnesses. Masked, stretched and mixed-face cases are `NotBuilt` in the product backend (section 9), so those rows use the harnesses and are labelled as such.
3. **First set, accuracy (all pass):** MLMG vs FFT on the product tip, 32³ to 128³, Neumann, periodic and Dirichlet, 1 and 2 ranks: worst 7.0e-14 against eps_H 1e-8 to 3.9e-8 (section 7.1). D-067 mode checks: 65 of 65 checks pass (section 7.2). Harness cases (uniform 100³ Dirichlet, hallways open, hallways sealed, stairwell with three zones, stretched 64³ at R=1, 8, 64): every backend pair on CPU and GPU agrees with the assembled-HYPRE reference to at most 5.0e-11, six orders below eps_H (section 7.3).
4. **R=64 slowness (task 2): cause established.** AMReX's HYPRE interface hands HYPRE the row-scaled matrix `D⁻¹A` (unit diagonal), which is not symmetric when the diagonal varies. On a stretched grid the diagonal varies (7.3× at R=64), and the FDS option set uses PCG, which needs a symmetric matrix. The first bottom solve then diverges until its 200-iteration cap; later solves recover. It is not an anisotropy/coarsening problem: BoomerAMG on the same unscaled matrix needs 27 iterations. Recommended setting: `hypre.hypre_solver=BiCGSTAB` with the bottom tolerance 1e-11 (about 8× faster at 64³, and the only converging setting at 1M and R=64 on CPU), until AMReX can pass the unscaled or symmetrically scaled matrix (section 8).
5. **Stairwell repeats (task 3):** section 7.5.
6. **Fold sign check (task 4):** the three signs of `fold_boundary_data` (low-side Neumann `+g/h`, high-side Neumann `-g/h`, Dirichlet `-2 H_b/h²`) are **right**: derived from `pres.f90` (section 11.1) and confirmed by a nonzero-data run against a dense matrix built the FDS way (24 of 24 cases, worst 4.4e-14; each sign flip gives differences of 0.31 to 3.6 relative to max |H|) and by a manufactured-solution convergence test (section 11.2 and 11.3). No FDS run with nonzero wall data and an H dump exists in the V&V area (section 11.4).
7. **Mixed N/D error excess (task 5):** not the solver tolerance and not the wall stencil. It is a smooth, domain-wide discretisation effect: the coarse truncation error outside the refined patch is no longer partly cancelled by the patch's own truncation error, and the C/F interface adds a smaller second-order remainder (section 12).
8. **Not measurable yet:** items in section 9 (HYPRE inside the product backend, composite masked/stretched/mixed-face cases, GPU build of the product tip, FDS H references).

## 2. Requirement and metric

- **A-56 (requirements, assumption list):** the backends agree within eps_H on frozen single solves for the rectangular-domain, masked and composite branches, and this check runs in CI.
- **eps_H (requirements §2.1a):** `eps_H = max(1e-8, 2.4e-12·N²)`, N = largest cell count per direction (floor governs up to N ≈ 64). Volume-weighted mean removed, relative L2. Valid only for one solve on frozen input on the same discretisation; not for multi-step outputs.
- **Implementation of the metric used here (all rows):** for each zone (connected pressure component), remove from each solution its volume-weighted mean over that zone's cells, then `rel_L2 = ‖a − b‖₂ / ‖b‖₂` over the zone's cells, with `b` the reference. Where cells have unequal volume the volume-weighted norm is also printed; the two differ by at most 1.6× in the stretched rows and the larger is quoted. The pinned cell of a sealed zone is excluded from the harness dumps (it carries the gauge value, not an unknown); this changes the mean by one cell in 10⁶. Scripts: `a56_cmp.py` (harness dumps), `pb_harness` CMP line (product tip).
- **eps_H values for the sizes used:** N = 32, 64: 1.0e-8; 100: 2.4e-8; 128: 3.9e-8; 320 (hallways refined 2.5×): 2.5e-7; 552 (stairwell, padded): 7.3e-7.
- **Same discretisation:** MLMG with `setMaxOrder(2)`, HYPRE assembled with the identical face conductances, FFT with the cell-centred stencil. For stretched cells the harness gives MLMG the scaled operator `K = Σ (AF/DX1)(x_i − x_j)` (unit-spacing geometry, FDS convention), so MLMG and HYPRE solve the same matrix.

## 3. Backends and what each may legally take

| Backend | Product backend (tip)? | Legal cases | Notes |
|---|---|---|---|
| FFT (`FFT::Poisson`) | yes | full uniform box, homogeneous Dirichlet/Neumann/periodic per face pair (plus lifted boundary data) | reference for the uniform box; 1-cell direction with Dirichlet is wrong (doc 06 §9) |
| MLMG, HYPRE bottom (`Mf`: FDS BoomerAMG options) | yes (uniform, composite) | box and masked (overset mask), stretched only through the harness operator | one MG level on pinned singular components |
| Assembled HYPRE (`H`: PCG + BoomerAMG as FDS sets it) | **no** (harness only) | any gas-cell set, stretched, masked, mixed faces | the FDS-equivalent reference |
| MLMG, HYPRE bottom with BiCGSTAB (`Mb`, new, section 8) | harness only | as Mf | recommended bottom for stretched cells |

## 4. Cases

| ID | Geometry | Size | Zones | Boundary | Backends legal | Rank counts | Status |
|---|---|---|---|---|---|---|---|
| U1 | uniform box | 32³, 64³, 128³ | 1 | Neumann, periodic, Dirichlet | FFT, MLMG (product tip); H, Mf in harness | 1, 2 (tip); 1, 8 (harness) | measured (7.1) |
| U2 | uniform box, manufactured solution | 100³ | 1 | Dirichlet | FFT, Mf, Mb, H | 1, 8, GPU | measured (7.3, 7.4) |
| H1 | hallways, refined 2.5× (320×160×160, 1,152,000 gas unknowns, 84 boxes) | 1.15M | 1 | walls + open vent slab | Mf, Mb, H | 1, 8, GPU | measured |
| H2 | hallways sealed, same mesh | 1.15M | 1 sealed, pinned | walls | Mf, Mb, H | 1, 8, GPU | measured |
| S1 | stairwell union (192×176×552 padded, `ba=drop`, two z-planes cut) | 1.63M | 3 (1 open, 2 sealed) | walls + vent | Mf, Mb, H | 1, 4, 8, GPU | measured (7.3, 7.5) |
| R1 | stretched z (geometric, last/first = R), pure Neumann (pinned) and Dirichlet | 64³ and 100³ | 1 | as stated | Mf, Mb, H; FFT only at R=1 | 4 (64³), 8 and GPU (100³) | measured |
| C1 | composite (2 and 3 levels), masked-composite | 32³ to 64³ | per level | Neumann, periodic, Dirichlet | MLMG (tip) vs uniform-fine FFT reference | 1, 2 | tip self-checks pass (7.2); comparison with H not measurable |
| C2 | mixed open/closed faces in the product backend | | | | | | not built (9) |

## 5. Metrics and pass/fail criteria

1. **Accuracy (hard criterion):** every pair of backends, per zone, `rel_L2 ≤ eps_H` (section 2). Also each solution's own true relative residual ≤ 1e-10 (the solve tolerance used for all backends) and ≥ 1e-14 reported for transparency.
2. **Rank and decomposition independence:** the same case on 1 and N ranks and two `max_grid_size` values agrees within eps_H (requirements FR iii); CPU and GPU runs of the same backend agree within eps_H.
3. **Gauge:** with D-067 defaults the removed mean, the pinned value and the zone gauge `Σρ·V·(KRES−H)=0` are checked by the tip's `meankind` and `comp_gauge` modes. The parity switch (`ScaledArithmetic`) must reproduce FDS's arithmetic mean removal (negative control in the same mode).
4. **Time (report, with regression triggers; not an A-56 pass/fail):** re-solve median over at least 10 timed solves per run, at least 3 runs per row for rows quoted as medians (2 for the wide first-set rows, marked), spread = (max − min)/median. A row with spread above 10% is rerun. Also first-solve time, setup time, iterations (MLMG V-cycles, bottom iterations where known) and peak memory. Trigger for investigation: Mf or Mb more than 3× slower than H on the same case, or any non-converged run.
5. **CI subset (A-56 "runs in CI"):** `pb_harness solve` at 32³ and 64³ (Neumann, periodic, Dirichlet, 1 and 2 ranks, both mean kinds) plus the `comp_gauge`, `meankind`, `exactsum` and `selector` modes: 22 solves and 65 checks, a few minutes on the shared machine. The harness cases (H1, H2, S1, R1) are too slow and need a GPU/large-memory runner, so they are nightly or manual, until HYPRE is in the product backend.

## 6. Protocol (timing rows)

- Quiet machine before every run: the pinned performance cores at least 85% idle over a 3 s sample, at most 1.5 busy CPUs overall, package temperature under 70 °C. Package temperature, pinned-core clocks and GPU clock are sampled once per second during every run.
- Abort and rerun a row if the package stays at or above 95 °C for 30 s, if the median pinned-core clock falls below 2.8 GHz (8 ranks), 3.3 GHz (4 ranks) or 4.2 GHz (1 rank), or if the GPU SM clock median falls below 2.0 GHz. Up to three attempts per run; every rejected attempt is logged.
- Short package-temperature peaks (setup of the assembled-HYPRE rows reaches 92 to 97 °C for a few seconds) are reported, not rejected, as the rule requires 30 s.
- Ranks pinned to performance cores 0, 2, 4, …; one rank is pinned to one core (the earlier doc 05/06 one-rank runs left the rank free to move between performance cores).
- Solve-only timing: device synchronised and a barrier around each solve; MLMG with `recompute_preconditioner=0`; H resets x to 0 and calls the PCG solve; stop at true relative L2 residual ≈ 1e-10 for all.
- Reuse: the doc 05 tables (mms 100³, hallD, sealD; CPU and GPU, power profile "high performance") and doc 06 §7 (stretched 1M) remain valid for solve time because the solve path of the harnesses did not change; they are quoted as "doc 05/06" where used and re-measured in 7.4.

## 7. Results

### 7.1 Uniform box, product tip: MLMG vs FFT (measured)

Harness `pb_harness solve`, D-067 defaults (`MeanKind::Volume`). Rows with the parity switch (`mean_kind=scaled`) are identical to the printed digits on a uniform grid (volume weight constant), which is itself a check. Product tip built unchanged from a copy of the committed tree.

| N | BC | Ranks | max_grid_size | FFT true rel. residual | MLMG true rel. residual | MLMG V-cycles | MLMG vs FFT rel. L2 | eps_H | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| 32^3 | neumann | 1 | 16 | 1.95e-14 | 1.09e-13 | 12 | 3.98e-14 | 1.00e-08 | PASS |
| 32^3 | periodic | 2 | 16 | 1.27e-14 | 1.64e-13 | 11 | 3.82e-14 | 1.00e-08 | PASS |
| 32^3 | dirichlet | 1 | 16 | 1.03e-14 | 4.02e-13 | 10 | 7.02e-14 | 1.00e-08 | PASS |
| 64^3 | neumann | 1 | 32 | 8.44e-14 | 2.95e-13 | 12 | 4.94e-14 | 1.00e-08 | PASS |
| 64^3 | neumann | 2 | 32 | 8.44e-14 | 2.95e-13 | 12 | 4.94e-14 | 1.00e-08 | PASS |
| 64^3 | periodic | 2 | 32 | 5.46e-14 | 4.15e-13 | 11 | 4.00e-14 | 1.00e-08 | PASS |
| 64^3 | dirichlet | 1 | 32 | 4.39e-14 | 3.02e-14 | 12 | 8.60e-15 | 1.00e-08 | PASS |
| 64^3 | dirichlet | 2 | 32 | 4.39e-14 | 3.02e-14 | 12 | 8.60e-15 | 1.00e-08 | PASS |
| 128^3 | neumann | 2 | 64 | 3.58e-13 | 1.22e-13 | 13 | 6.96e-14 | 3.93e-08 | PASS |
| 128^3 | periodic | 2 | 64 | 2.31e-13 | 9.95e-14 | 12 | 1.22e-14 | 3.93e-08 | PASS |
| 128^3 | dirichlet | 2 | 64 | 1.84e-13 | 4.95e-14 | 13 | 6.50e-14 | 3.93e-08 | PASS |

Worst MLMG-vs-FFT difference: 7.0e-14 (128³ Neumann), against eps_H 1.0e-8 (N ≤ 64) and 3.9e-8 (N = 128). The same rank counts (64³ Neumann on 1 and 2 ranks; 64³ Dirichlet on 1 and 2 ranks) give the same difference to three digits. This confirms and extends the earlier finding (requirements, P2 evidence) with the D-067 tip: the periodic 128³ and Dirichlet 128³ rows are new.

### 7.2 D-067 mode checks on the tip (measured)

`pb_harness` modes `meankind`, `comp_gauge`, `comp`, `comp_sel`, `exactsum`, `selector`, 2 ranks: **65 checks pass, none fail.** Examples: the removed mean is bitwise independent of the gauge fields; the weighted gauge constant differs from the plain one (the test is sensitive); the `ScaledArithmetic` solution satisfies the independent divergence formula to 0 and the negative control (Volume zero mode) does not (0.0145 against 2e-14); the composite solve's true residual is below eps_H; covered cells equal the average-down of the next finer level; the selector returns `NotBuilt` for covered cells, cylindrical geometry, stretched cells (also with MLMG explicitly requested), variable coefficients and mixed open/closed faces. Raw output: `scratch/pressure-signoff/a56/res/d067_modes.txt`.

### 7.3 Accuracy across cases, backends, ranks and devices (measured, harnesses)

Reference: assembled HYPRE (H), one CPU rank, one solve, tolerance 1e-10. Entries are the worst zone `rel_L2` against that reference after volume-weighted per-zone mean removal (section 2). `c1`, `c4`, `c8` = CPU ranks; `g` = GPU. Mf = MLMG with the FDS BoomerAMG options; Mb = the same with a BiCGSTAB bottom and bottom tolerance 1e-11 (section 8). The stretched rows are 64³ at R = 1, 8, 64 (offset 0, `mean=vol`). The 64³ stretched H reference is also one rank on the CPU.

| Case | eps_H | Mf c1 | Mf c4/c8 | Mb c4/c8 | Mf g | Mb g | H c4/c8 | H g | FFT c1/c4 | FFT g |
|---|---|---|---|---|---|---|---|---|---|---|
| U2 mms 100³ Dirichlet | 2.4e-8 | 5.8e-13 | 5.8e-13 | 5.8e-13 | 5.8e-13 | 5.8e-13 | 6.6e-13 | 6.5e-13 | 3.7e-13 | 3.7e-13 |
| H1 hallways open | 2.5e-7 | 5.0e-11 | 5.0e-11 | 5.0e-11 | 5.0e-11 | 5.0e-11 | 9.5e-13 | 7.6e-13 | - | - |
| H2 hallways sealed | 2.5e-7 | 2.2e-11 | 2.2e-11 | 3.1e-11 | 2.2e-11 | 2.2e-11 | 6.5e-13 | 4.9e-13 | - | - |
| S1 stairwell, 3 zones (worst) | 7.3e-7 | 1.2e-11 | 1.2e-11 | 1.7e-11 | 1.2e-11 | 1.6e-11 | 1.3e-12 | 2.8e-12 | - | - |
| R1 neu R=1 | 1.0e-8 | - | 1.5e-12 | 1.5e-12 | 1.5e-12 | - | 1.3e-12 | 4.4e-12 | 1.5e-12 | - |
| R1 neu R=8 | 1.0e-8 | - | 6.8e-12 | 8.2e-12 | 6.7e-12 | - | 5.9e-12 | 5.5e-12 | - | - |
| R1 neu R=64 | 1.0e-8 | - | 2.9e-12 | 4.6e-12 | 2.9e-12 | - | 5.6e-12 | 3.6e-12 | - | - |
| R1 dir R=1 | 1.0e-8 | - | 3.2e-12 | 3.2e-12 | 3.2e-12 | - | 9.9e-13 | 1.4e-12 | 9.9e-13 | - |
| R1 dir R=8 | 1.0e-8 | - | 1.5e-12 | 1.5e-12 | 1.5e-12 | - | 1.0e-12 | 7.7e-13 | - | - |
| R1 dir R=64 | 1.0e-8 | - | 1.2e-12 | 1.2e-12 | 1.2e-12 | - | 9.9e-13 | 9.1e-13 | - | - |

- **Every entry passes eps_H by at least four orders** (largest 5.0e-11 against 2.5e-7). The size of the differences is the iterative stopping error (Mf and Mb stop at max-norm 1e-10) rather than a discretisation difference; H on 8 ranks or on the GPU differs from H on one rank by 5e-12 or less.
- **Stairwell, per zone** (one OPEN zone 545,550 cells, two sealed zones 204,593 and 884,511 unknowns): Mf vs H 8.4e-12, 5.8e-13, 1.2e-11 on one rank; Mb 1.2e-11, 7.8e-13, 1.7e-11; GPU Mf 9.3e-12, 8.5e-13, 1.2e-11.
- **Reuse of earlier data:** the stairwell one-rank Mf vs H agreement was first reported in doc 06 §8 (same order of magnitude); doc 05 reported only residuals and the manufactured-solution error (identical for all backends to four digits, max 8.22e-5), so this table is the first backend-vs-backend metric on hallways. The doc 06 64³ stretched error table (Mf and H identical error to the printed digits) remains valid.
- Raw dumps and the table generator: `scratch/pressure-signoff/a56/test machine/logs/a56/acc_table.txt`, `a56_cmp.py`.

### 7.4 Solve time, first set (re-measured with the quiet-machine protocol)

PENDING: filled from the repeat runs (uniform Dirichlet 100³, hallways open and sealed, stretched 100³ R=8 Neumann and Dirichlet; Mf, H and, where it matters, Mb; 8 CPU ranks and 1 GPU; 2 repeats). Until then the doc 05 tables (mms, hallD, sealD) and doc 06 §7 (stretched 1M) stand; the single-run R=64 results at 1M are in section 8.

### 7.5 Masked stairwell repeats (task 3)

PENDING: Mf, H and Mb, 1/4/8 CPU ranks and GPU, 3 repeats each, median and spread; replaces the single-repeat table in doc 06 §8.

## 8. Why Mf is slow at stretching ratio R = 64 (task 2)

Case: `h2h_str case=neu` (pure Neumann, pinned last cell, `mean=vol`), 64³ unless stated, 4 CPU ranks, solve tolerance 1e-10. "Mf" = MLMG, HYPRE bottom with the FDS options (PCG, BoomerAMG, PMIS 8, l1-Jacobi 18, strong threshold 0.25, interpolation 6, one cycle per PCG iteration). Each experiment is one first solve plus one re-solve; times are re-solve medians of single runs (the same experiment varies by about 20%), so differences smaller than 1.3× are not significant. Scripts and raw logs: `r64_exp.sh`, `r64_batch2.sh`, `r64_batch3.sh`, `logs/pr06/run/r64_*.out` on the test machine (copies under `scratch/pressure-signoff/a56/test machine/`).

### 8.1 Instrumentation: where the time goes

- **MLMG:** `mg_levels = 1`. The pinned (masked) cell prevents coarsening, so the bottom solver is the whole solve. MLMG timers: bottom 2.62 s of 2.63 s. R=1 and R=8 need 3 V-cycles, R=64 needs 5.
- **Bottom (HYPRE PCG) iterations per MLMG iteration:**

| R | iteration 1 | 2 | 3 | 4 | 5 | total | re-solve (s) |
|---|---|---|---|---|---|---|---|
| 1 | 7 | 8 | 8 | – | – | 23 | 0.17 |
| 8 | 10 | 9 | 8 | – | – | 27 | 0.21 |
| 64 | **200 (cap)** | 9 | 29 | 21 | 18 | 277 | 2.0-2.5 |

  The first bottom solve at R=64 never converges: its preconditioned residual `‖r‖_C/‖b‖_C` starts at 0.42, falls to a minimum of 0.285, then grows steadily to 26 by iteration 200, where the 200-iteration cap (MLMG's bottom maximum) stops it. The correction it returns is wrong by orders of magnitude (fine-grid residual/bnorm 2.7e3 after MLMG iteration 1, from an initial 1.7e-3). The next outer iteration starts from a residual with different content and converges (0.08 to 6e-5 in 9 iterations), so the penalty is the first solve: 72% of all PCG iterations, and the bottom solve is 2.62 s of the 2.63 s MLMG time. At R=8 every bottom solve converges in 8 to 10 iterations.
- **Larger R:** R=512 does not converge at all (19 MLMG iterations, 32 s, true residual 7.5, `converged=0`). The Dirichlet case (no pin, 6 MG levels) does not show this: R=64 needs 117 V-cycles (point-smoothed geometric multigrid struggles with the strongly different cell widths) but only 0.45 s, and R=512 40+ V-cycles without converging within 0.18 s of 40 cycles. That is a smoother/anisotropy effect, separate from the Neumann slowness (see 8.3).
- **At 1M cells (100³, 8 ranks, single runs):** R=8 Mf 1.51 s (4 V-cycles); R=64 Mf **does not converge: 200 V-cycles, 54 s per solve, residual 2.8e20 (diverged)**. On the GPU R=64 Mf takes 5.2 s (6 V-cycles) against 0.145 s for H. The 2.3 to 7.2 s range in doc 06 is reproduced: 2.0 to 2.7 s per 64³ solve at 4 ranks, and a similar 7.2 s when the bottom tolerance is tightened with PCG (section 8.2).

### 8.2 Cause: the matrix AMReX gives HYPRE is not symmetric

Read in the AMReX HYPRE interface (`AMReX_HypreABecLap3.cpp`, `AMReX_Habec_3D_K.H`, `habec_ijmat`): every row is scaled by the inverse of its own diagonal (`sten[m] *= diaginv`, centre entry set to 1) and the right-hand side is multiplied by the same `diaginv`. For the symmetric operator `A` this hands HYPRE `D⁻¹A`, which is symmetric only when all diagonal entries are equal. On a stretched grid they are not (diagonal max/min: 2.0 at R=1, 2.6 at R=8, **7.3 at R=64**, 38 at R=512; the 2.0 at R=1 is the boundary rows). PCG assumes a symmetric positive definite operator and preconditioner; with `D⁻¹A` it breaks down.

Controls that isolate this (assembled-HYPRE harness, same 64³ matrix, same BoomerAMG options, PCG; measured):

| Matrix handed to PCG + BoomerAMG | R=8 | R=64 | R=512 |
|---|---|---|---|
| `A` (what H does) | 22 it (doc 06) | 27 it, 0.128 s | 39 it, 0.215 s |
| `D⁻¹A` (what AMReX does) | 107 it, 0.52 s | **1000 it, not converged, 5.2 s** | **1000 it, not converged** |
| `D⁻½ A D⁻½` (symmetric scaling) | 28 it, 0.13 s | 38 it, 0.20 s | 60 it, 0.32 s |

So: the same matrix, the same AMG options and the same pinning converge in 27 iterations when it is symmetric and fail when only the row scaling is added. A symmetric scaling (the variable diagonal handled both sides) costs 40% more iterations at R=64 but converges. In addition, a local AMReX build in which `diaginv` is forced to 1 (the one change, under a compile flag, in `habec_ijmat`) makes the **Mf route itself** converge in 3 V-cycles at every R: 0.26, 0.28, 0.38, 0.55 s at R = 1, 8, 64, 512 (64³, 4 ranks).

**Ruled out as the cause:**
- Anisotropy and the BoomerAMG strength/coarsening settings: BoomerAMG handles the unscaled stretched matrix in 27 to 39 iterations up to R=512 with the same threshold 0.25 and PMIS coarsening.
- Strength threshold 0.5 and 0.8, Falgout (6) and HMIS (10) coarsening, relaxation 6 and 8, direct interpolation, two AMG cycles per PCG iteration, all with PCG on the scaled matrix: none converges within 12 V-cycles (22 to 38 s), except direct interpolation (10 V-cycles, 12.8 s) and two cycles (6 V-cycles, 7.8 s), both slower than the baseline. Each is worse than or equal to the baseline, so the cause is not in the AMG hierarchy.
- `max_coarsening_level = 3`, agglomeration and consolidation: no effect (`mg_levels` stays 1 because the pinned cell blocks level-0 coarsening); 2.56 s against 2.47 s.
- MLMG's own bottom solvers: its BiCGSTAB without AMG does not converge in 12 V-cycles (true residual 1e-4); the smoother alone stops at residual 0.72.
- Tightening the PCG bottom tolerance to 1e-11 (scaled matrix, PCG): worse (4 V-cycles, 7.2 s).

### 8.3 Remedies tested (64³, R=64, 4 ranks)

| Setting | V-cycles | Re-solve (s) | Result |
|---|---|---|---|
| Mf baseline (PCG, bottom tolerance 1e-4) | 5 | 2.0-2.5 | converges, slow |
| `hypre.hypre_solver=GMRES` | 3 | 0.29 | converges |
| `hypre.hypre_solver=FlexGMRES` | 3 | 0.27 | converges |
| `hypre.hypre_solver=BiCGSTAB` | 3 | 0.29-0.32 | converges |
| BiCGSTAB, bottom tolerance 1e-8 / 1e-11 | 2 / **1** | 0.46 / **0.255-0.275** | best single setting |
| GMRES / FlexGMRES, bottom tolerance 1e-11 | 2 / 2 | 0.47-0.51 / 0.59 | converges |
| BoomerAMG as the solver (no Krylov) | 3 | 0.51 | converges |
| BiCGSTAB 1e-11 plus threshold 0.5 / 0.6 | 2 / 2 | 0.64 / 0.68 | worse |
| BiCGSTAB 1e-11 plus relaxation 6 / HMIS | 1 / 1 | 0.32 / 0.31 | no gain |
| assembled HYPRE (H), reference | 27 PCG it | 0.128 | |
| local AMReX build without the row scaling, PCG | 3 | 0.38 (0.277 with bottom tolerance 1e-11) | converges |

Bottom tolerance 1e-11 with BiCGSTAB makes the MLMG converge in one V-cycle at every R tested (R=1, 8, 64, 512: 1 V-cycle; 0.19, 0.19, 0.28, 0.35 s). Without the Krylov swap (still PCG) the same tolerance is slower.

**At 1M cells (100³, single runs, same machine, pure Neumann R=8 and R=64):**

| Setting | CPU 8 ranks R=8 | CPU R=64 | GPU R=8 | GPU R=64 |
|---|---|---|---|---|
| Mf (PCG) | 1.51 s (4) | not converged (200, 54 s/solve) | 0.365 s (3) | 5.19 s (6) |
| Mb: BiCGSTAB, bottom tolerance 1e-4 | 1.05 s (3) | 1.65 s (3) | 0.252 s (3) | 0.347 s (3) |
| Mb: BiCGSTAB, bottom tolerance 1e-11 | 0.94 s (1) | 1.44 s (1) | 0.218 s (1) | 0.298 s (1) |
| local AMReX build without row scaling, PCG | 1.00 s (3) | 1.35 s (3) | - | - |
| H | 0.43 s (24 PCG) | 0.62 s (32 PCG) | 0.111 s | 0.145 s |

(V-cycles or PCG iterations in parentheses.) With the fix the Mf route no longer depends on R, but it stays about 2× slower than H at 1M: BiCGSTAB spends two preconditioner applications per iteration, and with a PCG bottom MLMG restarts the Krylov solve three times at the loose bottom tolerance. This is the same pinned-singular-component weakness doc 05 found for sealed hallways (consistent with 4 restarts at 1e-4 against one H solve at 1e-10; Mb on sealed hallways at 8 ranks: 1 V-cycle, 1.90 s, against Mf 2.30 s and H 0.65 s).

### 8.4 Recommended setting

1. **Now (no AMReX change):** for any level whose diagonal varies (stretched cells, and in general any non-uniform or masked operator), set `hypre.hypre_solver=BiCGSTAB` and MLMG's bottom tolerance to 1e-11 (`MLMG::setBottomTolerance`). Keep the FDS BoomerAMG options. Expect R-independent iteration counts and about 2× the cost of assembled HYPRE at 1M.
2. **Product fix (candidate, measured with a local build only):** make AMReX's HYPRE interface optionally skip the row scaling (or scale symmetrically, `D⁻½AD⁻½`, with the matching vector scaling) so that PCG can be kept. The change is one block in `habec_ijmat` (about 5 lines; the right-hand-side multiplication by `diaginv` then becomes 1). It is not yet a runtime option and has not been run through the AMReX tests; propose upstream only after that.
3. For the uniform box nothing changes: R=1 Mf with PCG converges in 3 V-cycles.
4. The Dirichlet stretched case (no pin) is not slow in absolute terms but needs 117 V-cycles at R=64 because point-smoothed geometric multigrid does not cope with the cell-width ratio. Setting `max_coarsening_level = 0` with the BiCGSTAB bottom turns it into one HYPRE solve (1 V-cycle, 0.20 s at 64³ against 0.45 s) and keeps R=512 at 0.22 s; semi-coarsening is not needed because BoomerAMG's strength-based coarsening already handles the anisotropy on the whole level.

Not established: why the first solve diverges rather than merely stalls (the preconditioner-induced indefiniteness of `D⁻¹A` is the mechanism the controls support; the PCG internals were not instrumented beyond the residual history).

## 9. What is and is not measurable yet

| Item | State |
|---|---|
| MLMG vs FFT, uniform box, product backend | measurable now (7.1, CI subset in section 5) |
| MLMG vs assembled HYPRE | harness only: HYPRE is not a backend in the product tip |
| Masked domains (hallways, stairwell) in the product backend | `NotBuilt` (covered cells, variable coefficient `a`). At the tested commit mixed open/closed faces were also `NotBuilt`; at the later tip they are built (read from the commit log; used in section 11 and 12) |
| Stretched cells in the product backend | `NotBuilt` (selector refuses non-uniform widths, also for explicit MLMG) |
| Composite branch against H | not possible: no composite H assembly exists; the tip's own composite checks (7.2) and the uniform-fine FFT reference are the available evidence |
| GPU runs of the product tip | not measured: the tip was built for CPU on the shared machine only; GPU rows use the harnesses (same HYPRE, AMReX 99ddfda with CUDA) |
| FDS H as the reference | not available to me as an executed run: see section 11.4 for what exists in the V&V area; the assembled-HYPRE harness reproduces FDS's matrix and options but is not FDS output |
| Full-size CI | not before HYPRE is in the product backend (section 5) |

## 10. Failures, caveats and open items

- One diagnostic batch hung on a non-converging PCG variant (22 s per run times 12 V-cycles); I stopped it and repeated those variants with an iteration cap of 12. The test machine connection dropped twice during waits; no data were lost.
- Section 8 experiments are single runs on a shared machine (not the quiet-machine protocol); only the factor differences of 2× or more are significant. The 1M R=64 numbers are single runs.
- The GPU rows of the R=64 batch ran with the host package at 77 to 92 °C (GPU plateau) and low host clocks; GPU solves are device-bound, but treat them as indicative.
- The local AMReX build used for the unscaled experiment is a copy of the tree at 99ddfda with one compile flag; it is not committed anywhere.
- Open: runtime option or symmetric scaling in the AMReX HYPRE interface and its upstream test; a product-side guard that selects the Krylov method by the diagonal spread; HYPRE inside the product backend (needed for A-56 in CI); H dump from a real FDS case for the final A-56 sign-off.

## 11. Fold sign check (task 4): `fold_boundary_data` against FDS's own H

Tip `626c5f4b7f` (read only, built unchanged). Raw outputs: `fold_matrix.txt`, `fold_mms.txt` in the scratch results directory; driver and reference are `a56_fold.cpp` and `fold_ref.py` there.

### 11.1 Derivation from the FDS source (read)

Source: `Source/pres.f90`. FDS assembles `K x = F` for the pressure unknown `x`, with `H = -x` (`HP = -X_H`).

| Item | FDS | Consequence for the unit-operator form `laplacian(H) = b` |
|---|---|---|
| Matrix | `K = sum over faces of (AF/DX)(x_i - x_j)`, positive semi-definite; a Dirichlet face adds `2·AF/DX` to the diagonal | `V·laplacian_hom(H) = F`, so `b = F/V` (the sign flip of `x = -H` is absorbed) |
| Interior RHS | `F = V·PRHS` | `b = PRHS` |
| Neumann data | `BXS` and `BXF` are `dH/dx` at the low and high x face, built from `-FVX +/- DUNDT`; always the derivative along the increasing coordinate | `g = dH/dx_d` at both faces |
| Neumann at the low face | `F += +BXS·AF` | `b += g/h` |
| Neumann at the high face | `F += -BXF·AF` | `b -= g/h` |
| Dirichlet (either side) | `F += -2·(1/dx)·AF·BCV`, `BCV` = wall value of H | `b -= 2 H_b/h²` |
| Ghost cells in the FFT path | low N: `HP(0) = HP(1) - dx·BXS`; high N: `HP(IBP1) = HP(IBAR) + dx·BXF`; D: `HP(0) = -HP(1) + 2·BXS` | identical to the backend's documented ghost values |

AMReX `MLPoisson` solves `L phi = rhs` with `L = div grad` (the positive-sign Laplacian), and the FFT Poisson solver the same equation, so the backend's `phi` is H directly, not `-H`. The three folded signs therefore follow from the FDS assembly with no extra sign change. The one condition: the caller must pass `b = PRHS` with `phi = H`. If a driver passes the FDS unknown `x = -H`, every sign flips.

### 11.2 Nonzero-data run against a dense FDS-style reference (measured)

- Backend side: `a56_fold.cpp` builds the problem on a uniform cubic grid (n = 10), fills a smooth right-hand side, calls `fold_boundary_data`, then `solve_pressure` (MLMG and FFT, relative tolerance 1e-13).
- Data: every face carries nonzero data; constants of 0.5, -0.8, 1.1, 0.7, -0.4, 0.9 on the six faces ("const"), or the same constants plus a smooth tangential variation of amplitude 0.3 ("var", passed as a cell-centred slab `MultiFab`). Neumann faces get `g`, Dirichlet faces get `H_b`.
- Reference: `fold_ref.py` builds the dense matrix `K` and vector `F` exactly as the table above (each interior face adds `+h` on both diagonals and `-h` off-diagonal; Dirichlet adds `2h` on the diagonal; `F = V·PRHS` plus the three BC terms with `V = h³`, `AF = h²`), solves `K x = F` and takes `H = -x`. Singular (all-Neumann) systems: mean removed from `F` and from both fields, minimum-norm solution.
- Face-type sets: `NN,NN,NN` (singular), `DD,DD,DD`, `ND,DN,NN`, `DN,ND,DD`, `NN,DD,DN`, `ND,NN,NN`; each with const and var data, MLMG and FFT: **24 of 24 cases pass**, relative max difference between 1.2e-15 and 4.4e-14 (largest: `ND,NN,NN`, MLMG).
- Negative controls (the same dense reference with one sign flipped; difference to the backend result, relative to max |H|), on the cases that contain the flipped face type:

| Flipped | Difference to the backend, smallest to largest over the cases that have that face type |
|---|---|
| low-side Neumann `g/h` | 0.31 to 3.0 (20 cases) |
| high-side Neumann `g/h` | 0.38 to 3.6 (20 cases) |
| Dirichlet `2 H_b/h²` | 0.78 to 2.1 (20 cases) |
| all three | 1.7 to 2.0 (24 cases) |

(On cases without the flipped face type the control is trivially unchanged; those are not counted.) The folded backend therefore matches the FDS-convention reference to round-off, and the test is sensitive to each sign.

### 11.3 Independent check by manufactured solution (measured)

`pb_harness mode=bcdata part=mms` (in the tip's harness) uses `u = cos(1.3x+0.2) cosh(0.7y) sin(0.9z+0.4)` with exact wall values and exact normal derivatives, which does not use the backend's ghost convention. L2 error at n = 8, 16, 32 (identical for MLMG and FFT):

| Face types | n = 8 | n = 16 | n = 32 | Order |
|---|---|---|---|---|
| `NN,NN,NN` | 3.04e-4 | 7.62e-5 | 1.91e-5 | 2.0 |
| `DD,DD,DD` | 8.64e-4 | 2.23e-4 | 5.63e-5 | 1.95 to 1.98 |
| `ND,DN,NN` | 3.94e-4 | 9.89e-5 | 2.48e-5 | 2.0 |

A wrong sign on any data class would not converge to `u`. The tip's own `part=exact` test is circular (it builds the right-hand side with the same ghost convention that the fold uses); it checks the algebra of the fold, not the sign against FDS. The new test in 11.2 closes that gap because its reference comes from the FDS matrix assembly.

### 11.4 What exists in the V&V area for an FDS H comparison (read)

- A reference FDS binary exists in the V&V area.
- No run there dumps H with nonzero wall data. The H-dump request for A-09b is recorded as not yet requested. The A-38 runs are periodic (no wall data), and the duct-flow case is multi-mesh with HVAC.
- I did not invent a run. **Open:** an FDS case with a nonzero normal velocity or `DUNDT` on a wall and an H dump, to repeat 11.2 with the true FDS H.

### 11.5 Result

All three signs are right: low-side Neumann `rhs += g/h`, high-side Neumann `rhs -= g/h`, Dirichlet `rhs -= 2 H_b/h²`, with `g = dH/dx_d` along the increasing coordinate at both faces, `H_b` the wall value, and `phi = H`.

## 12. Mixed Neumann/Dirichlet: fine-level error larger than the coarse one (task 5)

Source of the question: the harness comment in `composite_modes.cpp` (mode `comp`, mixed faces) and the note on composite faces in the tip's notes. The harness relaxes the strict check "composite error on the fine level is below the uniform coarse error" for mixed faces, saying that the manufactured solution need not vanish at the patch edges, and checks the order instead. Tip `626c5f4b7f`, MLMG composite, 1 rank. Work files (scripts, numbers): `nd_excess.py`, `nd_superpos.py`, `nd_resid.py` and their outputs in the scratch results directory.

Setup (harness, unchanged): unit cube, two levels, the refined patch is the middle half of the domain in each direction (coarse cells n/4 to 3n/4), ratio 2, manufactured modes with half-integer wave numbers for N/D directions (so the exact solution matches the wall data), right-hand side = analytic Laplacian at the cell centres. n = 16, 32, 64 coarse cells. I dumped the solutions and rebuilt the errors in numpy (checked: my uniform-coarse error reproduces the harness dump to 1e-14).

### 12.1 Is it the solver tolerance? No (measured)

`ND,DN,NN`, n = 32, relative tolerance 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14 (4 to 10 MLMG iterations; the 1e-14 run stopped at the 200-iteration cap with a true residual of 4e-15): the error to the exact solution is 2.01846e-3 at 1e-4 and 2.01853e-3 for every tolerance from 1e-6 down to 1e-14. The error is fixed by the discretisation from the tolerance 1e-6 on; the default 1e-12 is far below it.

### 12.2 Order of convergence per region (measured)

L2 errors to the exact solution at n = 16, 32, 64, and the order between consecutive grids. "Fine patch" = composite fine level; "outside" = uncovered coarse cells; "uniform" = a single-level solve on the same cells restricted to the same region.

| Face types | Region | n = 16 | n = 32 | n = 64 | Orders |
|---|---|---|---|---|---|
| `ND,DN,NN` | composite, fine patch | 6.82e-3 | 1.73e-3 | 4.38e-4 | 1.98, 1.98 |
| | composite, coarse outside | 8.32e-3 | 2.06e-3 | 5.12e-4 | 2.02, 2.01 |
| | uniform fine, same patch | 1.37e-3 | 3.42e-4 | 8.54e-5 | 2.00, 2.00 |
| | uniform coarse, same outside | 5.51e-3 | 1.37e-3 | 3.42e-4 | 2.01, 2.00 |
| `ND,ND,ND` | composite, fine patch | 4.20e-3 | 1.08e-3 | 2.76e-4 | 1.95, 1.97 |
| | uniform fine, same patch | 1.61e-3 | 4.01e-4 | 1.00e-4 | 2.00, 2.00 |
| `DN,ND,DD` | composite, fine patch | 9.43e-4 | 2.55e-4 | 6.62e-5 | 1.89, 1.94 |
| | uniform fine, same patch | 1.37e-3 | 3.42e-4 | 8.54e-5 | 2.00, 2.00 |
| `NN,NN,NN` | composite, fine patch | 1.07e-3 | 2.80e-4 | 7.12e-5 | 1.94, 1.97 |
| `DD,DD,DD` | composite, fine patch | 3.05e-4 | 7.41e-5 | 1.86e-5 | 2.04, 2.00 |

Every region converges at second order. The excess is therefore a constant factor, not a loss of order. The ratio composite-fine / uniform-fine on the patch is constant in n: `ND,DN,NN` 4.98, 5.07, 5.13; `ND,ND,ND` 2.61, 2.70, 2.75; `NN,NN,NN` 1.71 to 1.81; `DN,ND,DD` 0.69 to 0.78; `DD,DD,DD` 0.22 to 0.23. For `ND,DN,NN` the composite fine error is 1.24 to 1.28 times the uniform coarse error on the same cells, which is what the harness's strict check flags. The sign of the effect depends on the face types: for `DD` the composite fine error is below even the uniform fine error.

Controls on the geometry: a fine level over the whole domain (`full=1`) gives exactly the uniform fine error (difference 1e-15), so the excess needs a C/F interface; a patch in the low corner (`layout=1`) gives a fine-level error of 3.03e-4, below the uniform one; ratio 4 raises the fine-level error to 2.31e-3 (ratio 2: 1.73e-3), which is the behaviour of a coarse-driven error (the coarse region is the same, the fine region gets no better).

### 12.3 Where the excess sits (measured)

- It is not near the walls. The refined patch does not touch the domain faces. The composite error in the uncovered coarse cells outside the patch is 1.50 times the uniform coarse error there (constant for n = 16, 32, 64), and it is flat from the wall layer to the interior of the domain (n = 16: 8.19e-3 within two cells of a wall, 8.57e-3 deeper). The excess in the fine patch is also flat: at n = 64 the rms error is 4.7e-4 in the two fine cells next to the interface, 4.65e-4 at 2 to 3 cells, 4.5e-4 at 4 to 7 cells and 4.0e-4 in the core. About 90 to 91% of the squared error of the whole composite lies in the coarse region outside the patch and 9 to 10% in the fine patch.
- The truncation error at the first cell layer at the walls is only 1.2 times that of the interior (rms 0.75 against 0.62 at n = 16), because the half-integer modes satisfy the wall data exactly. The cell-centred ghost-cell Dirichlet closure (ghost = 2 H_b - phi_1) does not produce a first-order local term here.
- Mechanism, tested by superposition. On the uniform coarse grid the error is `e = -L⁻¹ tau`, with `tau` the truncation error of the exact solution (my numpy operator reproduces the solver's coarse error exactly). Split `tau` into the patch footprint and the rest. For `ND,DN,NN` the error caused by the truncation outside the patch alone has rms 8.41e-3 at n = 16 (the measured composite error outside the patch is 8.32e-3), while the part caused by the truncation inside the patch is 6.0e-3 in the outside region and anti-correlated with it (correlation -0.56 at every n). In the uniform solve these two parts partly cancel (5.5e-3). Replacing the coarse cells of the patch with fine cells removes most of the patch part, so the cancellation is lost and the error outside rises to 1.5 times. Removing the patch's truncation entirely reproduces the measured composite error to a residual of 15% of its size (13% for the best-fit patch share of 0.08 instead of 0); the residual is the C/F interface and fine-grid part. It is largest in the cells touching the interface and falls away from it (n = 64, rms: outside, first cell 1.5e-4, second 1.2e-4, deeper 4.6e-5; inside, first cell 1.6e-4, second 1.4e-4, deeper 1.1e-4; wall layer 2.2e-5), and it also converges at second order (relative size 13 to 15% at all three n).
- The same mechanism explains the other signs: where the patch truncation reinforces rather than cancels the outside one (`DD`, `DN,ND,DD`), the composite fine error is lower than the uniform fine error.

### 12.4 Answer

- **Not the solver tolerance** (section 12.1).
- **Not primarily the boundary stencil order:** the wall layers carry a small share of the residual and the truncation error there is close to the interior's (12.3); the patch is far from the walls.
- **A discretisation effect that is smooth and domain-wide.** The fine-level error of a composite solve is not bounded by the fine grid's own accuracy; it carries the elliptic influence of the coarse truncation error outside the patch. With the half-integer N/D modes the coarse-outside and patch parts partly cancel in the uniform solve; the composite removes the cancellation (this reproduces about 85% of the composite error), and the C/F interface and the fine-grid part contribute the remaining 13 to 15%, concentrated in the cells next to it. All of it converges at second order (orders 1.89 to 2.04).
- **Consequence for the harness:** the relaxed check (order instead of "fine below uniform coarse") is right, but the stated reason ("the solution need not vanish at the patch edges") is not what the data show. A sharper acceptance test is the order per region plus the ratio to the uniform fine error constant in n (the values in 12.2), or comparing with a uniform solve on the same cells restricted to the same region.
- **Open:** I did not test non-matching wall data (a solution with nonzero second derivative at a Dirichlet wall), where the ghost-cell closure has an O(1) local truncation term; the global order there is the known second order, but I have not measured it.
