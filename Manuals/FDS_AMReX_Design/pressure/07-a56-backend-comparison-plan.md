# 07. A-56 backend comparison plan (MLMG vs assembled HYPRE vs FFT) and first results

**Status: plan complete (sections 1 to 6); first accuracy set measured (7.1 to 7.3, 7.6); first-set timing rows complete with the corrected driver (7.4; two strD8 CPU rows unstable after 5 repeats); R=64 cause established (section 8); stairwell repeats done (7.5); fold sign check (section 11) and mixed N/D excess (section 12) done; P3-R02 measured on the composite MLMG, family 1 (section 13) and family 2 with a moved patch (section 16); product HYPRE against MLMG and FFT, single level and composite measured, masked not built (section 14); residual_tol size scaling measured with a recommendation (section 15). Not committed.**
Measured = run. Read = source only. Unverified = neither.
Software: sections 1 to 10 (except where stated) use the common-layer tip of the pressure backend (the read-only source tree, commit `1f3f4de4c4`; latest `Source/pressure_backend` commit `29c2f9aa4a`, which contains the D-067 defaults from `a4a1073952`), read only and built unchanged. Sections 11 and 12 use the later tip `626c5f4b7f` (adds mixed open/closed faces, `fold_boundary_data`, the D-057 mapping check), built unchanged in a scratch directory. AMReX 26.09 and HYPRE v2.32.0-24 for that tip; the assembled-HYPRE and stretched/masked harnesses (docs 05, 06) use AMReX 99ddfda and the same HYPRE.
Machines: CPU timing and all GPU runs on an owner-provided NVIDIA test machine (pinned performance cores, 1 thread per rank, one GPU, managed arena). CPU accuracy checks of the tip ran on the shared build box.

## 1. Answers in short

1. **Plan:** A-56 is checked per zone by one number, `eps_H = max(1e-8, 2.4e-12·N²)` (requirements §2.1a), on a single frozen solve on the same discretisation, after removing the volume-weighted mean per zone, relative L2 (section 2). Cases, sizes, rank counts and pass/fail rules are in sections 3 to 6.
2. **What is measurable today:** the product backend (common layer + MLMG + FFT + assembled HYPRE) can be compared FFT/MLMG/HYPRE on the uniform box and HYPRE against MLMG on composite hierarchies, including mixed face types (section 14). Masked cells, covered cells on the single-level API, stretched cells and cylindrical geometry are `NotBuilt` for every product backend (section 9, 14.3), so masked hallways, the stairwell and the stretched cases use the study harnesses and are labelled as such.
3. **First set, accuracy (all pass):** MLMG vs FFT on the product tip, 32³ to 128³, Neumann, periodic and Dirichlet, 1 and 2 ranks: worst 7.0e-14 against eps_H 1e-8 to 3.9e-8 (section 7.1). D-067 mode checks: 65 of 65 checks pass (section 7.2). Harness cases (uniform 100³ Dirichlet, hallways open, hallways sealed, stairwell with three zones, stretched 64³ at R=1, 8, 64): every backend pair on CPU and GPU agrees with the assembled-HYPRE reference to at most 5.0e-11, six orders below eps_H (section 7.3).
4. **R=64 slowness (task 2): cause established.** AMReX's HYPRE interface hands HYPRE the row-scaled matrix `D⁻¹A` (unit diagonal), which is not symmetric when the diagonal varies. On a stretched grid the diagonal varies (7.3× at R=64), and the FDS option set uses PCG, which needs a symmetric matrix. The first bottom solve then diverges until its 200-iteration cap; later solves recover. It is not an anisotropy/coarsening problem: BoomerAMG on the same unscaled matrix needs 27 iterations. Recommended setting: `hypre.hypre_solver=BiCGSTAB` with the bottom tolerance 1e-11 (about 8× faster at 64³, and the only converging setting at 1M and R=64 on CPU), until AMReX can pass the unscaled or symmetrically scaled matrix (section 8).
5. **Stairwell repeats (task 3):** section 7.5. 3 repeats per row, 36 runs, none rejected; H is 2.6× faster than Mf on CPU and 6.2× on GPU (8 CPU ranks: Mf 2.70 s, H 1.04 s; GPU: Mf 1.19 s, H 0.193 s); BiCGSTAB (Mb) does not help on this uniform-spacing case.
6. **Fold sign check (task 4):** the three signs of `fold_boundary_data` (low-side Neumann `+g/h`, high-side Neumann `-g/h`, Dirichlet `-2 H_b/h²`) are **right**: derived from `pres.f90` (section 11.1) and confirmed by a nonzero-data run against a dense matrix built the FDS way (24 of 24 cases, worst 4.4e-14; each sign flip gives differences of 0.31 to 3.6 relative to max |H|) and by a manufactured-solution convergence test (section 11.2 and 11.3). No FDS run with nonzero wall data and an H dump exists in the V&V area (section 11.4).
7. **Mixed N/D error excess (task 5):** not the solver tolerance and not the wall stencil. It is a smooth, domain-wide discretisation effect: the coarse truncation error outside the refined patch is no longer partly cancelled by the patch's own truncation error, and the C/F interface adds a smaller second-order remainder (section 12).
8. **P3-R02 pressure tolerance (section 13):** on the composite MLMG (two levels, D-067 defaults) the measured `max|div u - D - c|` after projection is **2.8e-12 at eps_rel 1e-12 and 6.9e-9 at 1e-9** for n = 32 (U/dx_fine = 64), i.e. 4.4e-14 and 1.1e-10 times U/dx_fine, against the working bound 1e-9·U/dx_fine (met with margins of 22,000 and 9); every projected case satisfies the derived bound `10·eps_rel·B + 20·eps_mach·U/dx_fine` (tightest margin 11 at 1e-12), and the unprojected control fails it. The working bound is broken only at eps_rel 1e-6.
9. **Timing, first set (section 7.4):** re-measured with the corrected driver, 2 repeats (5 for the unstable rows), quiet-machine protocol, no run rejected. H is faster than Mf on the GPU in every case (1.8 to 6.2×: sealed hallways 6.2×, stretched Dirichlet 4.3×) and on the CPU for sealed hallways (3.4×), stretched Neumann (3.3×) and Dirichlet (1.5×); slower on the CPU for the uniform and open-hallways cases (0.6×). The `sealD` GPU Mb value (7.5 s, 2 V-cycles) is **confirmed** by two more repeats (0.1% spread), so the row is real. The `strD8` CPU Mb and H rows are unstable (spread 363% and 27% over 5 repeats; medians 0.96 s and 0.39 s, orientation only).
10. **Product HYPRE (section 14):** against FFT on the uniform box (15 cases, 32³ to 160³, Neumann, periodic, Dirichlet) the worst relative L2 difference is 2.0e-12 (MLMG against FFT 1.5e-13), 3e4 or more below eps_H. Against MLMG on 46 composite solves (two and three levels, ratio 2 and 4, Neumann, periodic, Dirichlet, four mixed-face sets, corner and two-patch layouts, up to 19.1 million unknowns) the worst difference is 7.1e-12, 12,500 or more below eps_H; MLMG hit its 200-cycle cap in the five three-level ratio-4 cases (true residual still 6e-14 to 2.5e-13), HYPRE converged in all 46 (23 to 43 GMRES iterations). Masked cases are `NotBuilt` for every product backend, so no masked comparison exists at product level.
11. **`residual_tol` (section 15):** the relative residual of a backend at its round-off floor grows like N² for smooth data (FFT 1.5e-14 at 32³ to 3.9e-13 at 160³; HYPRE Neumann 9.8e-14 to 1.5e-12, crossing 1e-12 between 128³ and 160³), because `||A||·||H||/||b||` grows like N²; the normwise backward error does not grow (FFT 1.1 u, HYPRE 4e-16 to 6e-16). The solution error stays at most 5.2e-4 of eps_H in all 96 solves, also where the residual is above 1e-12. Recommendation: keep 1e-12 as base and use `max(1e-12, c·2^-53·||A||·||H||_2/||b||_2)` with the measured constant at the floor at most 5.6 (HYPRE Neumann; at most 1.4 otherwise), so c = 10 (the value already coded in the informational `residual_floor`). Three MLMG overshoots (1.2e-12 to 1.3e-12, 3 of 24 default solves) are V-cycle granularity, not size, and vanish at `tol_rel` 5e-13.
12. **P3-R02 family 2 (section 16):** noisy velocity (seam 14 to 111 before reprojection), three levels, patch moved by n/8 level-0 cells, nonuniform D with nonzero mean: at eps_rel 1e-12 the error is 1.0e-13 to 1.6e-13 times U/dx_fine in the main cases (2.1e-14 to 3.9e-13 over all variations), margins 13 to 102 against the derived bound and 2,500 to 47,000 against the working bound 1e-9·U/dx_fine; at 1e-9 margins 10 to 160 (derived) and 5.7 to 72 (working, tightest 5.7 for a patch moved by 0.19 of the domain); at 1e-6 the working bound fails (0.01 to 0.09) while the derived bound holds.
13. **Not measurable yet:** items in section 9 (masked and stretched cases in the product backend, GPU build of the product tip, FDS H references).

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
| Assembled HYPRE (`H`: PCG + BoomerAMG as FDS sets it) | **yes, explicit request only** (single level PCG, composite GMRES(30)); masked, stretched, covered cells on the single-level API are `NotBuilt` there, so those rows stay on the harness | product: uniform box and composite hierarchies with mixed faces; harness: any gas-cell set, stretched, masked | the FDS-equivalent reference (harness); section 14 (product) |
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
| C1 | composite (2 and 3 levels, ratio 2 and 4), masked-composite | 32³ to 96³ coarse | per level | Neumann, periodic, Dirichlet | product HYPRE vs MLMG (FFT not legal on composite); masked composite not built | 4 | measured (14.2); masked not built (14.3); tip self-checks pass (7.2) |
| C2 | mixed open/closed faces in the product backend | 32³, 64³ coarse | per level | `ND,DN,NN`, `DD,NN,NN`, `NN,DD,DN`, `PP,NN,DN` | product HYPRE vs MLMG | 4 | built at the later tip; measured (14.2, section 12) |

## 5. Metrics and pass/fail criteria

1. **Accuracy (hard criterion):** every pair of backends, per zone, `rel_L2 ≤ eps_H` (section 2). Also each solution's own true relative residual ≤ 1e-10 (the solve tolerance used for all backends) and ≥ 1e-14 reported for transparency. The product's own residual warning (`residual_tol`) is judged with the size-aware limit of section 15.5, `max(1e-12, 10·2^-53·||A||·||H||_2/||b||_2)`, on the pin-excluded residual; it is a warning criterion, never a replacement of the eps_H solution criterion.
2. **Rank and decomposition independence:** the same case on 1 and N ranks and two `max_grid_size` values agrees within eps_H (requirements FR iii); CPU and GPU runs of the same backend agree within eps_H.
3. **Gauge:** with D-067 defaults the removed mean, the pinned value and the zone gauge `Σρ·V·(KRES−H)=0` are checked by the tip's `meankind` and `comp_gauge` modes. The parity switch (`ScaledArithmetic`) must reproduce FDS's arithmetic mean removal (negative control in the same mode).
4. **Time (report, with regression triggers; not an A-56 pass/fail):** re-solve median over at least 10 timed solves per run, at least 3 runs per row for rows quoted as medians (2 for the wide first-set rows, marked), spread = (max − min)/median. A row with spread above 10% is rerun. Also first-solve time, setup time, iterations (MLMG V-cycles, bottom iterations where known) and peak memory. Trigger for investigation: Mf or Mb more than 3× slower than H on the same case, or any non-converged run.
5. **CI subset (A-56 "runs in CI"):** `pb_harness solve` at 32³ and 64³ (Neumann, periodic, Dirichlet, 1 and 2 ranks, both mean kinds) plus the `comp_gauge`, `meankind`, `exactsum` and `selector` modes: 22 solves and 65 checks, a few minutes on the shared machine. The harness cases (H1, H2, S1, R1) are too slow and need a GPU/large-memory runner, so they are nightly or manual (masked and stretched cases cannot go through the product API at all). The product HYPRE backend can now join the subset: the uniform solves of section 14.1 (`backends="fft mlmg hypre"`) at 32³ and 64³ and the composite `hypre_cmp` cases at 32³ coarse (two and three levels, ratio 2, Neumann, periodic, Dirichlet and the mixed-face sets) each take under 1 s on 4 ranks, with the eps_H criterion of section 2.

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
- Raw dumps and the table generator: `acc_table.txt` in the scratch results directory of the test machine logs, `a56_cmp.py`.

### 7.4 Solve time, first set (re-measured with the quiet-machine protocol)

Cases: uniform Dirichlet 100³ (`mms`), hallways open (`hallD`) and sealed (`sealD`), stretched 1M at R = 8 Neumann (`strN8`) and Dirichlet (`strD8`). Backends Mf, H and Mb (Mb: Mf with a BiCGSTAB solver). 8 CPU ranks (pinned performance cores) and 1 GPU; each repeat is a fresh process, interleaved by row, and the quiet-machine gate and rejection rules are as in 7.5 (at most 1.5 busy CPUs from other processes, all performance cores at least 85% idle, package below 70 °C before each run; a run is rejected for 30 s at or above 95 °C, a collapsed core clock or a GPU clock below 2000 MHz). No run was rejected or skipped for load in the second pass (package at most 96 °C, never 30 s at or above 95 °C; held clocks 3.4 to 3.6 GHz on 8 ranks in the second pass, GPU SM 2430 MHz).

**History of the rows.** The first pass had a driver fault: a shell variable holding the case name was overwritten by the package-temperature variable of the quiet-machine check, so most H and Mb output files were named after a temperature and overwritten. The second pass uses a corrected driver (case name in its own variable, one output file per run, tag `<device>_<case>_<backend>_r<repeat>`), 2 repeats of every affected H and Mb row, plus the `mms` H, `sealD` GPU Mb and `strN8` GPU rows for confirmation. The Mf rows and the `sealD` CPU Mb and `strN8` CPU rows come from the first pass, whose files for those rows were unique and valid. Two rows (strD8 on 8 CPU ranks, H and Mb) are marked unstable after 5 repeats each (below); every other row has its repeats.

| Device | Case | Backend | Repeats | Re-solve per repeat (s) | Median (s) | Spread | Iterations | Pkg max (C) | Clock (MHz) |
|---|---|---|---|---|---|---|---|---|---|
| 8 CPU ranks | mms | Mf | 2 | 0.2415, 0.2373 | 0.2394 | 1.7% | 23 | 93 | 4343 |
| 8 CPU ranks | mms | H | 2 | 0.3691, 0.4239 | 0.3965 | 13.8% | 21 | 92 | 3600 |
| 1 GPU | mms | Mf | 2 | 0.1678, 0.1678 | 0.1678 | 0.0% | 23 | 92 | 2430 |
| 1 GPU | mms | H | 2 | 0.09429, 0.09426 | 0.09427 | 0.0% | 21 | 86 | 2430 |
| 8 CPU ranks | hallD | Mf | 2 | 0.3744, 0.3721 | 0.3733 | 0.6% | 15 | 94 | 4382 |
| 8 CPU ranks | hallD | H | 2 | 0.6511, 0.6301 | 0.6406 | 3.3% | 26 | 95 | 3589 |
| 8 CPU ranks | hallD | Mb | 2 | 0.4934, 0.4967 | 0.4950 | 0.7% | 15 | 93 | 3478 |
| 1 GPU | hallD | Mf | 2 | 0.2413, 0.2411 | 0.2412 | 0.1% | 15 | 95 | 2430 |
| 1 GPU | hallD | H | 2 | 0.1194, 0.1193 | 0.1193 | 0.0% | 23 | 89 | 2430 |
| 1 GPU | hallD | Mb | 2 | 0.3012, 0.3013 | 0.3013 | 0.0% | 15 | 93 | 2430 |
| 8 CPU ranks | sealD | Mf | 2 | 2.295, 2.274 | 2.284 | 0.9% | 4 | 94 | 3588 |
| 8 CPU ranks | sealD | H | 2 | 0.6729, 0.6630 | 0.6679 | 1.5% | 27 | 96 | 3594 |
| 8 CPU ranks | sealD | Mb | 2 | 1.917, 1.913 | 1.915 | 0.2% | 1 | 96 | 3589 |
| 1 GPU | sealD | Mf | 2 | 0.8349, 0.8350 | 0.8350 | 0.0% | 4 | 92 | 2430 |
| 1 GPU | sealD | H | 2 | 0.1346, 0.1346 | 0.1346 | 0.0% | 26 | 90 | 2430 |
| 1 GPU | sealD | Mb | 2 (+1 earlier) | 7.501, 7.509 (earlier 7.502) | 7.505 | 0.1% | 2 | 95 | 2430 |
| 8 CPU ranks | strN8 | Mf | 2 | 1.479, 1.479 | 1.479 | 0.0% | 4 | 93 | 3552 |
| 8 CPU ranks | strN8 | H | 2 | 0.4484, 0.4481 | 0.4482 | 0.1% | 24 | 93 | 3575 |
| 8 CPU ranks | strN8 | Mb | 2 | 0.9216, 0.9237 | 0.9227 | 0.2% | 1 | 92 | 3544 |
| 1 GPU | strN8 | Mf | 2 | 0.3615, 0.3615 | 0.3615 | 0.0% | 3 | 90 | 2430 |
| 1 GPU | strN8 | H | 2 | 0.1103, 0.1104 | 0.1104 | 0.1% | 24 | 86 | 2430 |
| 1 GPU | strN8 | Mb | 2 | 0.2166, 0.2167 | 0.2166 | 0.0% | 1 | 86 | 2430 |
| 8 CPU ranks | strD8 | Mf | 2 | 0.5938, 0.5992 | 0.5965 | 0.9% | 54 | 93 | 3331 |
| 8 CPU ranks | strD8 | H | 5 | 0.3761, 0.3966, 0.3876, 0.4811, 0.3880 | 0.3880 | 27% (**above 10%**) | 21 | 93 | 3595 |
| 8 CPU ranks | strD8 | Mb | 5, **unstable** | 0.9434, 4.147, 0.6660, 0.9594, 1.441 | 0.9594 | 363% | 54 | 94 | 3389 to 4274 |
| 1 GPU | strD8 | Mf | 2 | 0.3979, 0.3974 | 0.3977 | 0.1% | 54 | 89 | 2430 |
| 1 GPU | strD8 | H | 2 | 0.09209, 0.09211 | 0.09210 | 0.0% | 20 | 85 | 2430 |
| 1 GPU | strD8 | Mb | 2 | 0.6541, 0.6554 | 0.6547 | 0.2% | 54 | 91 | 2430 |

Reading (re-solve medians, speed-up = Mf time / other time; values below 1 mean slower than Mf):

| Case | CPU H | CPU Mb | GPU H | GPU Mb |
|---|---|---|---|---|
| mms | 0.60 | not run | 1.78 | not run |
| hallD | 0.58 | 0.75 | 2.02 | 0.80 |
| sealD | 3.42 | 1.19 | 6.20 | 0.11 |
| strN8 | 3.30 | 1.60 | 3.27 | 1.67 |
| strD8 | 1.54 | 0.62 (unstable) | 4.32 | 0.61 |

- H is faster than Mf on the GPU in every case (1.8 to 6.2×) and on the CPU where Mf is slow (sealed hallways 3.4×, stretched Neumann 3.3×, stretched Dirichlet 1.5×); it is slower on the CPU for the uniform and the open hallways case (0.60 and 0.58 of Mf's speed), where Mf needs 15 to 23 V-cycles of cheap smoothing.
- The 13.8% spread of the `mms` H CPU row (two repeats, 0.369 and 0.424 s) is the largest among the valid rows; its median is a coarse number. The earlier single value 0.392 s lies between the two.
- The `sealD` GPU Mb row is **confirmed**, not a driver artefact: three runs (7.501, 7.509, earlier 7.502 s) agree within 0.1%, 2 V-cycles, so each V-cycle costs about 3.7 s against about 0.1 to 0.2 s per V-cycle on the other GPU rows. The cause was not investigated (the sealed case is singular and the GPU BiCGSTAB path with the sealed-component handling is the obvious suspect, but that was not tested); the row is a real property of that configuration, and Mb should not be recommended on the GPU for sealed domains.
- Mb on the GPU is slower than Mf for open hallways (0.80) and Dirichlet stretched (0.61), faster for stretched Neumann (1.67); on the CPU it is faster only for sealed hallways (1.19) and stretched Neumann (1.60).
- The `strD8` CPU rows (H and Mb) are **unstable and unresolved**, not missing. After the first two repeats of Mb disagreed (0.943 and 4.147 s) three more repeats of Mb and of H were run in a later quiet window (gate passed before each run, no run rejected; `time_rep4.log`): Mb re-solve medians 0.666, 0.959 and 1.441 s, H 0.388, 0.481 and 0.388 s. The first solve of Mb is stable (0.643 to 0.668 s in all five runs) while its re-solve time drifts upward within and between runs (p90 of the timed re-solves up to 5.6 s in the second run, 2.9 s in the fifth), with the same 54 iterations and held clocks. The medians (Mb 0.96 s, H 0.39 s) are given for orientation only; the spread (363% and 27%) exceeds the 10% rerun limit of section 5 and the cause of the drift on the CPU Dirichlet stretched case was not found (the GPU rows of the same case have 0.0 to 0.2% spread).
- **Cause analysis of the `strD8` CPU rows (log-based, no new run).** Ruled out with the per-run logs: (a) *iteration counts and arithmetic*: every Mb run has 54 iterations and every H run 21, and the true residual, error norms and pin-row residual are identical to all printed digits in all five repeats of each backend, so the numerical work is the same; (b) *first-solve effects*: the first solve is stable (Mb 0.643 to 0.668 s, H 0.305 to 0.322 s, spread under 5%) and the instability sits only in the timed re-solves (Mb median 0.666 to 4.147 s, with at least half of the 20 re-solves of the second run at 4.1 s or more: median 4.147 s, p90 5.554 s); (c) *memory*: high-water mark per rank 87.7 to 93.4 MB in all Mb runs; (d) *clock and temperature*: the slowest Mb run held 3600 MHz on all eight pinned cores (per-second samples, no core below about 3200 MHz after the first two seconds) at 71 C, while the fastest ran at 4274 MHz and 92 C; the time multiplied by the median pinned-core clock is 3200, 14929, 2846, 3250 and 4323 (Mb) and 1352, 1418, 1393, 1676, 1393 (H), so clock explains the fast and the ordinary runs (three of five Mb runs within 14%, four of five H runs within 5%) but not the slow ones (Mb second run 4.6 times, fifth run 1.4 times; H fourth run 1.2 times the ordinary value); (e) *GPU activity*: the slowest Mb run had an idle GPU, so it is not GPU-job interference; (f) *hardware threads*: the odd-numbered logical CPUs (the sibling threads of the pinned performance cores) are offline on the test machine, so no hyper-thread sharing; (g) *the quiet gate*: it passed before every run (0.5 to 0.8 busy CPUs, performance-core idle 0.90 to 0.99).
  **Not ruled out, and the best-supported hypothesis:** foreign CPU load on a pinned core *during* the run. The gate samples three seconds before the start and nothing logs per-CPU load during the run; eight synchronously communicating ranks run at the speed of the most contended core, and an 80 s slowdown that starts after a normal first solve and varies within the run (p90 2.9 s against a 1.4 s median in the fifth run) fits that pattern better than any property of the solver. At the time of this analysis a foreign process held about 99% of one of the pinned performance cores (1.6 busy CPUs, above the gate), which shows the situation occurs on this machine; it was not caused by these runs and cannot be tied to the earlier repeats. The check to run when the machine is quiet: repeat Mb `strD8` eight times with `/proc/stat` sampled per CPU at 1 Hz during the run, and reject a run when a pinned core shows more than 2% non-rank use or a foreign process is on the pinned list; record the per-solve times, not only the median and p90. This was not run: the gate failed when it was checked once for this analysis, and no waiting was done. The rows stay marked unstable.
- Doc 05 tables and doc 06 §7 are superseded by this table for the rows above; the single-run R = 64 results at 1M are in section 8.

### 7.5 Masked stairwell repeats (task 3) (measured, quiet-machine protocol)

Case as in doc 06 §8: stairwell union, 1,634,654 unknowns, three components (one open, two sealed), synthetic right-hand side, tolerance 1e-10, no warm-up. Rows: Mf (MLMG + HYPRE with the FDS option set, PCG), H (assembled HYPRE, PCG), Mb (Mf with `hypre.hypre_solver=BiCGSTAB` and bottom tolerance 1e-11, the section 8.4 setting). **3 repeats per row, interleaved by repeat** (all rows once, then again), each repeat a fresh process. Before every run the machine had to be quiet (pinned performance cores at least 85% idle, at most 1.5 busy CPUs overall, package below 70 °C); a run was to be rejected and repeated if the package stayed at or above 95 °C for 30 s, if the median pinned-core clock fell below 4200 (1 rank), 3300 (4 ranks) or 2800 (8 ranks) MHz, or the GPU SM clock median below 2000 MHz. **No run was rejected: all 36 were accepted at the first attempt.** The package touched 97 °C once (a short setup spike on one H run, never 30 s), otherwise at most 94 °C. The median is of the three per-repeat re-solve medians; "spread" is (max - min)/median over the three repeats.

| Backend | Ranks | Repeats | Re-solve median (s) | Min-max (s) | Spread (max-min)/median | Iterations | First solve median (s) | Setup median (s) | Timed solves per run | Package max (°C) | Core clock median (MHz) / GPU SM (MHz) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Mf | 1 CPU | 3 | 5.785 | 5.704-5.864 | 2.8% | 3 | 9.907 | 4.2 | 10 | 94 | 5200 / - |
| H | 1 CPU | 3 | 2.216 | 2.214-2.258 | 2.0% | 27 | 2.225 | 3.25 | 20 | 92 | 5200 / - |
| Mb | 1 CPU | 3 | 7.413 | 7.396-7.524 | 1.7% | 1 | 11.63 | 4.21 | 10 | 88 | 5200 / - |
| Mf | 4 CPU | 3 | 2.781 | 2.779-2.992 | 7.6% | 3 | 4.458 | 1.58 | 20 | 92 | 4192 / - |
| H | 4 CPU | 3 | 1.09 | 1.084-1.166 | 7.5% | 25 | 1.009 | 1.46 | 40 | 94 | 4200 / - |
| Mb | 4 CPU | 3 | 3.772 | 3.552-3.831 | 7.4% | 1 | 5.257 | 1.71 | 20 | 94 | 4162 / - |
| Mf | 8 CPU | 3 | 2.7 | 2.538-2.768 | 8.5% | 3 | 3.336 | 0.632 | 20 | 94 | 3574 / - |
| H | 8 CPU | 3 | 1.036 | 0.9443-1.042 | 9.4% | 26 | 0.795 | 0.989 | 40 | 97 | 3540 / - |
| Mb | 8 CPU | 3 | 3.421 | 3.123-3.62 | 14.5% | 1 | 3.929 | 0.508 | 20 | 93 | 3565 / - |
| Mf | 1 GPU | 3 | 1.191 | 1.191-1.191 | 0.1% | 3 | 1.858 | 0.667 | 30 | 92 | - / 2430 |
| H | 1 GPU | 3 | 0.1934 | 0.1934-0.1936 | 0.1% | 26 | 0.1952 | 0.223 | 50 | 93 | - / 2430 |
| Mb | 1 GPU | 3 | 1.263 | 1.262-1.263 | 0.1% | 1 | 1.929 | 0.666 | 30 | 92 | - / 2430 |

Reading: (1) The CPU repeats agree within 1.7 to 14.5% (largest at 8 ranks, where 8 pinned cores share a package at 3.5 GHz); the single-rank and GPU rows within 3%. (2) H is faster than Mf by 2.6× (8 ranks), 2.6× (4), 2.6× (1 rank) and 6.2× (GPU). (3) On this case Mb (BiCGSTAB) is **not** faster than Mf: 1 V-cycle instead of 3, but 1.3× slower per solve on CPU and 1.06× on GPU. The BiCGSTAB setting is a remedy for stretched cells (section 8); on this uniform-spacing masked case the diagonal does not vary, so Mf's PCG is not affected and BiCGSTAB only adds cost. (4) Compared with the single-repeat table of the first version of doc 06 §8, the GPU medians agree to 1%; the CPU medians differ by 5 to 18% (the first version ran at lower, uncontrolled core clocks: the old single-rank Mf row, 6.52 s, is now 5.79 s at a held 5.2 GHz, and the old single-rank H row 2.70 s is now 2.22 s).

Raw log: `stair_rep.log` and the table `stair_rep_table.md` in the scratch results directory.

### 7.6 Rerun of the uniform CI subset and the D-067 modes on the later tip (measured)

Tip `626c5f4b7f` (mixed faces, `fold_boundary_data`, D-057 mapping check), built unchanged. The same 11 configurations as in 7.1 with both mean kinds (22 runs, MLMG vs FFT, 1 and 2 ranks): **22 of 22 pass, worst relative L2 difference 7.0e-14** (eps_H 1e-8 to 3.9e-8), identical to the earlier tip within round-off. D-067 mode checks (`meankind`, `comp_gauge`, `comp`, `comp_sel`, `exactsum`, `selector`, 2 ranks): **65 of 65 checks pass, none fail.** One batch of the 128³ cases was killed by memory pressure on the shared machine and was repeated alone; the repeated runs pass.

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
| MLMG vs FFT, uniform box, product backend | measurable now (7.1, CI subset in section 5; 15 cases to 160³ in 14.1) |
| Product HYPRE vs FFT / MLMG, uniform box | measured (14.1): worst difference to FFT 2.0e-12 |
| Product HYPRE vs MLMG, composite 2 and 3 levels, ratio 2 and 4, mixed faces, corner and two-patch layouts | measured (14.2): 46 solves, worst 7.1e-12; FFT is not legal on a composite hierarchy, MLMG is the reference |
| Assembled HYPRE on masked hallways and stairwell, and on stretched cells | harness only (7.3 to 7.5, 8): the product backends return `NotBuilt` for masked cells, covered cells on the single-level API, non-uniform cell widths and variable coefficients (14.3) |
| Masked or stretched composite cases | not possible (`NotBuilt`); the missing piece is the masked branch in the common layer, not a HYPRE feature |
| Mixed open/closed faces in the product backend | built at the later tip (read from the commit log); measured with HYPRE and MLMG (14.2, sections 11 and 12) |
| Composite branch against H | not possible: no composite H assembly exists; the tip's own composite checks (7.2), the uniform-fine reference (14.2, `full=1`) and HYPRE against MLMG are the available evidence |
| GPU runs of the product tip | not measured: the tip was built for CPU on the shared machine only; GPU rows use the harnesses (same HYPRE, AMReX 99ddfda with CUDA) |
| FDS H as the reference | not available to me as an executed run: see section 11.4 for what exists in the V&V area; the assembled-HYPRE harness reproduces FDS's matrix and options but is not FDS output |
| Residual criterion | measured (15); the recommended size-aware limit is not yet in the code |
| Full-size CI | the product HYPRE can join the CI subset (section 5); the harness cases (hallways, stairwell, stretched 1M) stay nightly or manual |

## 10. Failures, caveats and open items

- One diagnostic batch hung on a non-converging PCG variant (22 s per run times 12 V-cycles); I stopped it and repeated those variants with an iteration cap of 12. The test-machine connection dropped twice during waits; no data were lost.
- Section 8 experiments are single runs on a shared machine (not the quiet-machine protocol); only the factor differences of 2× or more are significant. The 1M R=64 numbers are single runs.
- The GPU rows of the R=64 batch ran with the host package at 77 to 92 °C (GPU plateau) and low host clocks; GPU solves are device-bound, but treat them as indicative.
- The local AMReX build used for the unscaled experiment is a copy of the tree at 99ddfda with one compile flag; it is not committed anywhere.
- Scratch cleanup: the large regenerable files (rebuilt FDS binaries, objects, restart and run dumps of the mean-removal study; raw field dumps of the D-067 mode checks and of the HYPRE composite runs; build trees of the study drivers; the transferred R=64 archive) were deleted once their results were in docs 06 and 07; scripts, logs and result text are kept.
- Timing rows, second pass (7.4): the first pass had the driver fault described there; in the second pass the quiet gate was passed before every run and none was rejected. The accuracy matrices of sections 14 to 16 were also gated (at most 1.5 busy CPUs from others) with a 20-minute waiting budget; part of the HYPRE matrix was skipped for load in the first pass (19 rows) and rerun in a second pass without skips (`hy_retry.log`), so all rows of section 14 are present. The two `strD8` CPU rows remain unstable after 5 repeats (cause not proven: the log analysis in 7.4 rules out iterations, first-solve, memory, clock, temperature, GPU activity and hyper-threading; foreign load on a pinned core during the run is the open hypothesis).
- Section 14 timings are single runs, 4 ranks, indicative only; section 14.2 iteration counts and differences are exact. MLMG reports `NotConverged` in five three-level ratio-4 cases with a true residual below 3e-13; I did not investigate why its own residual stagnates there.
- Section 15: one right-hand-side family per spectrum type, uniform single-level only; the composite and masked residual floors are not measured. Section 16: single rank, hand-written face interpolation for the regrid, synthetic velocity.
- Open: runtime option or symmetric scaling in the AMReX HYPRE interface and its upstream test; a product-side guard that selects the Krylov method by the diagonal spread; the masked branch of the common layer (needed for any masked backend comparison at product level); the size-aware residual limit in `evaluate_residual` (section 15.6); the cause of the `sealD` GPU Mb cost per V-cycle (7.4) and of the `strD8` CPU drift; H dump from a real FDS case for the final A-56 sign-off.

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

## 13. P3-R02: divergence error after projection on the composite MLMG (measured)

Question (test plan P3-R02): the bound `max|div u - D - c| <= 10·eps_rel·B + 20·eps_mach·U/dx_fine`, with eps_rel the solver relative tolerance, B = max|div u* - D| before the projection, eps_mach = 2.2e-16, c the constant removed by the mean removal, compared with the working bound 1e-9·U/dx_fine (accepted at eps_rel = 1e-12).

Setup (tip `626c5f4b7f`, D-067 defaults, volume-weighted mean removal; driver `a56_proj.cpp` in the scratch directory; results `p3r02_results.txt`):
- Two-level composite hierarchy on the unit cube, refined patch = the middle half of the coarse domain in each direction, uncovered coarse cells plus fine cells.
- Face velocity `u*` sampled analytically on every level (smooth field with zero normal component at closed walls, amplitude U); coarse faces under the fine level replaced by the average of the fine faces (`average_down_faces`), the same rule the backend uses for the gradient. Cell source `D = U·(0.3 + 0.5 sin(3πx) cos(2πy) sin(πz))`.
- Projection: `b = div u* - D` on uncovered cells (composite divergence), `solve_pressure` for `lap H = b`, `face_gradient_composite` for the face gradients of H, `u = u* - grad H`. The observable is `max|div u - D - c|` over the uncovered cells of both levels, where c is the removed mean (closed box: the volume mean of b, which equals the value the backend reports; open or mixed faces: no singularity, c = 0).
- B (a priori maximum of b) is 9.7 for U = 1; U/dx_fine = 64 (n = 32, ratio 2), 128 (n = 64 or n = 32 with ratio 4).

| Case (face types) | eps_rel | V-cycles | max abs(div u - D - c) | divided by U/dx_fine | Derived bound divided by U/dx_fine | Working bound 1e-9 met | Derived bound met |
|---|---|---|---|---|---|---|---|
| n=32, ratio 2, `NN,NN,NN` | 1e-12 | 11 | 2.84e-12 | 4.4e-14 | 1.5e-12 | yes | yes |
| | 1e-9 | 8 | 6.85e-9 | 1.1e-10 | 1.5e-9 | yes | yes |
| | 1e-6 | 6 | 1.47e-6 | 2.3e-8 | 1.5e-6 | no | yes |
| | none (no projection) | - | 9.39 | 0.15 | 1.5e-12 | no | **no** |
| n=64, ratio 2, `NN,NN,NN` | 1e-12 | 12 | 2.01e-12 | 1.6e-14 | 7.6e-13 | yes | yes |
| | 1e-9 | 9 | 2.03e-9 | 1.6e-11 | 7.6e-10 | yes | yes |
| | none | - | 9.42 | 0.074 | 7.6e-13 | no | **no** |
| n=32, ratio 4, `NN,NN,NN` | 1e-12 | 11 | 5.48e-12 | 4.3e-14 | 7.6e-13 | yes | yes |
| | 1e-9 | 9 | 7.70e-10 | 6.0e-12 | 7.6e-10 | yes | yes |
| n=32, ratio 2, `DD,DD,DD` | 1e-12 | 10 | 8.87e-12 | 1.4e-13 | 1.5e-12 | yes | yes |
| | 1e-9 | 8 | 2.02e-9 | 3.2e-11 | 1.5e-9 | yes | yes |
| | none | - | 9.69 | 0.15 | 1.5e-12 | no | **no** |
| n=32, ratio 2, `ND,DN,NN` | 1e-12 | 12 | 1.82e-12 | 2.9e-14 | 1.5e-12 | yes | yes |
| | 1e-9 | 9 | 3.91e-9 | 6.1e-11 | 1.5e-9 | yes | yes |
| n=32, ratio 2, `NN,NN,NN`, U = 100 | 1e-12 | 11 | 2.78e-10 | 4.4e-14 | 1.5e-12 | yes | yes |
| | 1e-9 | 8 | 6.85e-7 | 1.1e-10 | 1.5e-9 | yes | yes |
| n=64, ratio 2, `NN,NN,NN`, 2 ranks | 1e-12 | 12 | 2.07e-12 | 1.6e-14 | 7.6e-13 | yes | yes |

(The derived-bound column shows the bound itself divided by U/dx_fine, to be compared with the measured column. At eps_rel = 1e-12 the bound is dominated by the 10·eps_rel·B term (9.7e-11); the machine-precision term 20·eps_mach·U/dx_fine is 2.8e-13 for n = 32, ratio 2, i.e. 0.3% of it.)

Findings:
1. **At eps_rel = 1e-12 the measured error is 1.6e-14 to 1.4e-13 times U/dx_fine** (2.0e-12 to 8.9e-12 absolute for U = 1), which is 7,000 to 62,000 times below the working bound 1e-9·U/dx_fine. It is also inside the derived bound in every case, with margins of 11 (`DD,DD,DD`, the tightest), 18 (ratio 4), 34 (n = 32), 49 (n = 64) and 53 (`ND,DN,NN`). The measured value is 0.2 to 0.9 of `eps_rel·B`, so the factor 10 in the derived bound is about right in size (the Dirichlet case uses 0.9 of it).
2. **At eps_rel = 1e-9 the measured error is 6.0e-12 to 1.1e-10 times U/dx_fine**, still inside the 1e-9·U/dx_fine working bound (margin 9 to 170) and inside the derived bound (measured 0.07 to 0.4 of it). At 1e-6 the working bound is exceeded (2.3e-8 times U/dx_fine) while the derived bound still holds, as expected: the working bound 1e-9·U/dx_fine is valid for solver tolerances of about 1e-9 or tighter, not looser.
3. **Scaling:** the error is proportional to U (U = 100 gives 100 times the absolute error and the same ratio to U/dx_fine) and follows the solver tolerance times B, not the grid: from eps_rel 1e-12 to 1e-9 the error grows by 140 to 2,400 times (a factor 1,000 in the tolerance; the V-cycle count steps in whole cycles, so the measured value jumps with it). The same holds at n = 32 and 64, ratio 2 and 4, and with 2 ranks.
4. **Negative control:** with no projection (phi = 0, the solver correcting nothing) the observable is B itself (9.4 to 9.7, i.e. 0.07 to 0.15 times U/dx_fine) and fails the bound by about 10 orders at eps_rel 1e-12; the test is sensitive.
5. **Mean removal:** the constant c in the closed case is -0.29995 (n = 32), the volume mean of b, equal to the backend's reported `removed_mean` to all printed digits; subtracting it is what makes the observable small. With open faces c = 0 as expected.

Limits: the field `u*` is smooth and synthetic and the hierarchy is static (a centred patch), not a regrid seam of a real run; the observable is therefore the pure solver-plus-composite-gradient part, which is what the bound is derived for. A regrid-seam test with a moved patch, noisy velocity and three levels is in section 16 (a hand-written regrid, not the code's own; a seam from a real run is still open). The bound's factor of 10 and the U/dx_fine floor were not stressed (eps_mach term 0.3% of the bound at 1e-12); a case with B much smaller than U/dx_fine would test the floor.

## 14. Product HYPRE backend against MLMG and FFT: single level, composite, masked (measured)

Scope. The assembled HYPRE backend is now part of the product backend (explicit request `BackendKind::HYPRE`; `Auto` never selects it). Per its notes it solves single-level problems with PCG and composite problems with GMRES(30) on the same operator as MLMG, with mixed face types, inhomogeneous face data through the fold, and the pin handling for singular problems. Tip `99c55244e1`, CPU, 4 ranks, on the owner-provided NVIDIA test machine (gate: at most 1.5 busy CPUs from other processes before each run). Harness modes `solve` (backends `fft mlmg hypre`), `hypre_cmp` (HYPRE against MLMG on the same composite problem, with the true residual evaluated by MLMG's operator) and `comp` (composite against a uniform solve); case matrix and raw output in the scratch results directory (`hy_matrix.log`, `hy_retry.log`, `hy_table.md`). The criterion is the relative L2 difference after mean removal against `eps_H = max(1e-8, 2.4e-12·N²)` with N the cells per direction of the finest level.

### 14.1 Single level, uniform box: FFT, MLMG and HYPRE (measured)

N = 32, 64, 96, 128, 160, Neumann, periodic, Dirichlet (15 cases, three backends each, `max_grid_size` 32 to 80, 4 ranks):

| Comparison | Cases | Worst relative L2 difference | Worst relative max difference | eps_H at the worst case | Result |
|---|---|---|---|---|---|
| MLMG against FFT | 15 | 1.5e-13 (160³ Neumann) | 1.3e-13 | 6.1e-8 | all pass, margin at least 4e5 |
| HYPRE against FFT | 15 | 2.0e-12 (160³ periodic) | 2.8e-10 | 6.1e-8 | all pass, margin at least 3e4 |

HYPRE iterations (PCG) 21 to 36 against MLMG 10 to 13 V-cycles. The HYPRE true residual (relative 2-norm, pin row included) is 3e-14 to 1.8e-13 up to 64³, 7.7e-12 (Neumann) and 7.4e-11 (periodic) at 96³, and 4.5e-11 (Neumann) and 4.3e-10 (periodic) at 160³; the pin-aware check of section 15 judges the non-pin value. The largest difference to FFT (2.0e-12, 160³ periodic) belongs to a singular, pinned case; Neumann is at most 1.3e-12, Dirichlet at most 4.6e-13.

### 14.2 Composite two- and three-level hierarchies: HYPRE against MLMG (measured)

FFT is not legal on a composite hierarchy (it needs one box with uniform spacing), so MLMG is the reference; the product HYPRE backend was run on the same `PressureProblem` (same hierarchy, same right-hand side). 46 solves:

| Group | Solves | Levels | Finest N | HYPRE GMRES iterations | MLMG V-cycles | Worst relative L2 difference | Worst relative max difference | Worst HYPRE true residual |
|---|---|---|---|---|---|---|---|---|
| Neumann, ratio 2 | 5 | 2, 3 | 64 to 256 | 24 to 41 | 10 to 11 | 7.1e-12 | 4.6e-10 | 2.0e-10 |
| Neumann, ratio 4 | 4 | 2, 3 | 128 to 1024 | 27 to 43 | 11 to 200 | 6.8e-13 | 1.0e-11 | 4.4e-12 |
| Periodic, ratio 2 | 4 | 2, 3 | 64 to 256 | 23 to 29 | 9 | 2.0e-13 | 2.1e-13 | 1.3e-13 |
| Periodic, ratio 4 | 4 | 2, 3 | 128 to 1024 | 27 to 39 | 10 to 200 | 2.2e-12 | 2.9e-12 | 2.2e-12 |
| Dirichlet, ratio 2 | 5 | 2, 3 | 64 to 256 | 23 to 31 | 10 to 12 | 9.5e-14 | 9.1e-14 | 1.6e-13 |
| Dirichlet, ratio 4 | 4 | 2, 3 | 128 to 1024 | 27 to 41 | 10 to 200 | 6.3e-13 | 7.4e-13 | 2.0e-12 |
| Mixed faces (`ND,DN,NN`, `DD,NN,NN`, `NN,DD,DN`, `PP,NN,DN`), ratio 2 | 16 | 2, 3 | 64 to 256 | 23 to 28 | 9 to 12 | 3.9e-13 | 1.8e-13 | 2.2e-13 |
| Corner patch and two-patch layouts, Neumann, ratio 2 | 4 | 2 | 64 to 128 | 23 to 35 | 10 to 17 | 5.3e-14 | 1.1e-12 | 4.2e-13 |

- **Agreement:** 46 of 46 differences are between 9.9e-15 and 7.1e-12, at least 12,500 times below eps_H (smallest margin: the 96³ coarse, 192³ fine Neumann case, 7.1e-12 against 8.8e-8; the finest-N 1024 cases have eps_H 2.5e-6 against at most 2.2e-12). HYPRE converged and passed its own true-residual check (relative 2-norm by MLMG's operator, at most 2.0e-10 with the pin row; its non-pin value is lower) in all 46, and every HYPRE composite solve was bitwise repeatable (fresh workspace and reused set-up) with the set-up reuse flag set.
- **MLMG did not converge in 5 of the 46** (`NotConverged` after the 200-V-cycle cap): all three-level, ratio-4 cases (n = 32 Neumann and Dirichlet, n = 64 Neumann, periodic, Dirichlet). Its true residual in those cases is 6e-14 to 2.5e-13 and the HYPRE solution agrees to 3.5e-13 to 2.2e-12, so the MLMG solution is right; the status is its own residual measure stagnating. HYPRE converged in 31 to 43 GMRES iterations in the same cases. (The five harness `CHECK FAIL` lines of the run are exactly these MLMG status lines, nothing else failed.)
- **Cost (CPU, 4 ranks, single runs, indicative; first solve with set-up / solve reusing the set-up):** where MLMG converges, HYPRE takes 1.4 to 3.8 times as many iterations and a median 7.8 times the time of MLMG for the first solve (3.1 to 14.5 times; 14.5 at 192³ Neumann, 3.3 s against 0.23 s). Where MLMG hits its cap HYPRE is faster: 32³ three-level ratio 4 (2.4 million unknowns) MLMG 13.4 s against HYPRE 4.5 s first, 2.1 s reused; 64³ three-level ratio 4 (19.1 million unknowns) MLMG 88 to 115 s against HYPRE 52 to 56 s first, 24 to 29 s reused. Per the HYPRE notes, GMRES and the assembled matrix are not meant to beat MLMG on a well-behaved hierarchy.
- **Fine level over the whole domain (`full=1`), two levels ratio 2:** the composite solution on the fine level against a uniform single-level solve of the same backend at the fine resolution (n = 32 and 64, Neumann and Dirichlet, both backends): 7e-16 to 2.5e-14 relative L2 (HYPRE 7.6e-15, 8.2e-15, 1.4e-15, 2.5e-14; MLMG 7.3e-16 to 6.3e-15). For HYPRE this is a composite-against-single-level check inside the same backend; the single-level HYPRE against FFT is 14.1. Without `full` the difference between composite and uniform fine solution is the discretisation difference (1.8e-4 to 7.5e-4 relative), identical for both backends to six digits, which shows that HYPRE and MLMG solve the same composite problem.

### 14.3 Masked cases: not legal for the product backend, so not compared there

The selector checks of the harness (`mode=selector`, all pass) show that an explicit product HYPRE request returns `NotBuilt` for masked cells (`cell_class != 0`: "masked cells ... are not built (masked branch)"), covered cells on the single-level API, cylindrical geometry, non-uniform cell widths (stretched meshes) and variable coefficients; the `Auto` request (which would pick FFT or MLMG) returns `NotBuilt` for masked and covered cells as well. The HYPRE notes list anisotropic refinement ratios and ratios other than 2 and 4 as `NotBuilt` too (read, not run here). Masked hallways and the stairwell therefore cannot be run through the product API with any backend, and the product HYPRE backend cannot take the stretched 1M cases. The masked and stretched comparisons of sections 7.3 to 7.5 remain on the study harnesses (assembled HYPRE and AMReX MLMG on the FDS-style operator); they are labelled as such there. What is missing for a product-level masked comparison is the masked branch in the common layer, not the HYPRE backend.

### 14.4 Not run

GPU (the product tip has no GPU build here); composite with masks or stretching (NotBuilt); ratios other than 2 and 4 (NotBuilt); more than 4 ranks on composite cases; composite n = 96 Dirichlet and larger HYPRE composite runs than the 19.1 million unknown case; FFT as a reference for composite problems (not legal).

## 15. Does `residual_tol` need a size-scaled form? (measured)

Question: the pin-aware residual check of the HYPRE backend reports that the non-pin residual of every backend grows slowly with N and that the fixed `residual_tol = 1e-12` would be crossed above about 128³. Measure the residual against N for the three backends, find out what it is relative to, compare it with the real solution error, and recommend a form with a measured constant.

### 15.1 Method

Product tip `99c55244e1` (contains the pin-aware check), built on the owner-provided NVIDIA test machine, 4 CPU ranks, unit cube, N cells per direction N = 32, 48, 64, 96, 128, 160 (160³ = 4.1 million cells), backends FFT, MLMG and HYPRE, all-Neumann (singular, with the common-layer mean removal and, for HYPRE, the identity pin) and all-Dirichlet. Driver `a56_resid.cpp` (scratch directory): pick a known solution H*, set `b = L H*` with the same 7-point operator and ghost rules as the backends, solve with the product API at `tol_rel` 1e-12 (the default), and read back the product's own residual quantities plus the error against H*.
- H* "smooth": two low-wavenumber modes plus a constant (so `||b||_2` per cell is 31 to 36 for every N, `max|b|` 129 to 435); H* "noisy": the same plus grid-scale noise of amplitude 0.3 (`||b||_2` per cell 560 to 14,000, growing like N²).
- Reported per solve: `rel2 = ||r||_2/||b||_2` (the quantity `residual_tol` is compared with), the same with the pin cell excluded (`rel2_nopin`, the quantity actually checked for singular problems), `relmax_nopin = max|r|/max|b|`, the normwise backward error `||r||_2/(||b||_2 + ||A||·||H||_2)`, the solution error `||H - H*||_2/||H*||_2` (mean removed for the closed box) and `eps_H = max(1e-8, 2.4e-12·N²)`.
- 96 solves: 72 at `tol_rel` 1e-12 (FFT, MLMG, HYPRE, smooth and noisy, Neumann and Dirichlet) and 24 at `tol_rel` 1e-15 for smooth data (MLMG and HYPRE iterate to their own floor), plus three MLMG reruns at 5e-13 and 2e-13 and two negative controls (`tol_rel` 1e-8, 1e-10). Raw outputs: `resid_matrix.log`, `resid_table.md`, `resid_const.md`, `resid_mlmg_tol.log`, `resid_negctl.log` in the scratch results directory. Gate: at most 1.5 busy CPUs from other processes before each run; none was skipped.

### 15.2 Residual against N

Neumann, smooth right-hand side, `tol_rel` 1e-12: relative residual `||r||2/||b||2` (bold: above 1e-12)

| N | FFT | MLMG | HYPRE (non-pin) | HYPRE (full, pin row included) |
|---|---|---|---|---|
| 32 | 1.5e-14 | 4.2e-13 | 9.8e-14 | 1.1e-13 |
| 48 | 3.4e-14 | 3.8e-13 | 1.5e-13 | 1e-11 |
| 64 | 6.6e-14 | **1.2e-12** | 2.4e-13 | 3e-13 |
| 96 | 1.4e-13 | 1.1e-13 | 5.3e-13 | 1.4e-10 |
| 128 | 2.5e-13 | 2.6e-13 | 9.6e-13 | 1e-12 |
| 160 | 3.9e-13 | 2.5e-13 | **1.5e-12** | 7.3e-10 |

Dirichlet, smooth right-hand side, `tol_rel` 1e-12

| N | FFT | MLMG | HYPRE (non-pin) |
|---|---|---|---|
| 32 | 1.2e-15 | 3.3e-13 | 8.2e-14 |
| 48 | 1.5e-15 | **1.2e-12** | 8.7e-14 |
| 64 | 1.9e-15 | 1.9e-13 | 5.6e-14 |
| 96 | 2.4e-15 | 4.4e-13 | 8.1e-14 |
| 128 | 2.9e-15 | 7.2e-13 | 5e-14 |
| 160 | 3.3e-15 | 9.7e-13 | 4.4e-14 |

Neumann, noisy right-hand side (grid-scale noise, 0.3), `tol_rel` 1e-12

| N | FFT | MLMG | HYPRE (non-pin) |
|---|---|---|---|
| 32 | 9.4e-16 | 2.9e-13 | 3.5e-14 |
| 48 | 9.5e-16 | 1.3e-13 | 5.5e-14 |
| 64 | 1e-15 | 1.9e-13 | 6.2e-14 |
| 96 | 1.1e-15 | 1.1e-13 | 5.8e-14 |
| 128 | 1.2e-15 | 1.4e-13 | 5.2e-14 |
| 160 | 1.1e-15 | 9.1e-14 | 7.7e-14 |

Dirichlet, noisy right-hand side, `tol_rel` 1e-12

| N | FFT | MLMG | HYPRE (non-pin) |
|---|---|---|---|
| 32 | 9.2e-16 | 2.3e-13 | 7.1e-14 |
| 48 | 1e-15 | 7.3e-13 | 7e-14 |
| 64 | 1.2e-15 | **1.3e-12** | 4.1e-14 |
| 96 | 1.2e-15 | 2.1e-13 | 5.5e-14 |
| 128 | 1.3e-15 | 3e-13 | 3.1e-14 |
| 160 | 1.4e-15 | 3.7e-13 | 9.3e-14 |

Neumann, smooth, `tol_rel` 1e-15 (solver iterates to its own floor; MLMG stops at its 200-cycle cap)

| N | MLMG | HYPRE (non-pin) |
|---|---|---|
| 32 | 2.9e-15 | 6e-14 |
| 48 | 6.2e-15 | 1.4e-13 |
| 64 | 9.8e-15 | 2.4e-13 |
| 96 | 2e-14 | 5.4e-13 |
| 128 | 3.5e-14 | 9.7e-13 |
| 160 | 5.2e-14 | **1.5e-12** |

Dirichlet, smooth, `tol_rel` 1e-15

| N | MLMG | HYPRE (non-pin) |
|---|---|---|
| 32 | 1.5e-15 | 1.5e-15 |
| 48 | 6.5e-16 | 1.9e-15 |
| 64 | 1.3e-15 | 2.2e-15 |
| 96 | 5.8e-16 | 3e-15 |
| 128 | 8e-16 | 3.4e-15 |
| 160 | 1e-15 | 4e-15 |


Fitted exponents, `d ln rel2_nopin / d ln N` (6 sizes): FFT Neumann smooth 2.02 and noisy 0.13; HYPRE Neumann smooth 1.75 (2.01 at `tol_rel` 1e-15) and noisy 0.33; HYPRE Dirichlet -0.4 to 0.6 (flat); MLMG: no trend (-0.6 to 0.4 at 1e-12; 1.77 for Neumann at 1e-15, where the solver reaches the floor and the 200-cycle cap).

What the numbers say:
1. **The growth is real but belongs to smooth data.** For a smooth right-hand side the relative residual of a backend that reaches round-off (FFT directly; HYPRE, which stops at 4e-16 backward error) grows like N²: FFT Neumann 1.5e-14 at 32 to 3.9e-13 at 160, HYPRE Neumann 9.8e-14 to **1.5e-12**. HYPRE crosses 1e-12 between 128³ (9.6e-13) and 160³ (1.5e-12), in agreement with the backend developer's estimate; FFT would cross at about 260³. For noisy data and for Dirichlet data the relative residual does not grow (FFT 1e-15 to 3e-15, HYPRE 3e-14 to 9e-14).
2. **The reason is the scale, not the solver.** `||b||_2` is N-independent for smooth data while the matrix norm `||A|| = 12 N²` grows, so `||A||·||H||_2/||b||_2` and with it the round-off size of the residual `~ u·||A||·||H||/||b||` (u = 2⁻⁵³) grows like N². The normwise backward error `||r||/(||b|| + ||A||·||H||)` does not grow with N for any backend (log-log slope over N: FFT 0.09 to 0.14, HYPRE -0.9 to 0.34, MLMG -2.5 to -0.1; FFT 0.9e-16 to 1.3e-16, i.e. 1.1 u; HYPRE Neumann smooth 4.0e-16 to 6.2e-16, flat). So a relative-to-`||b||` threshold fixed in N is reached by the data, not by a loss of accuracy.
3. **MLMG does not follow the round-off floor at the default tolerance; it follows its stopping rule.** Its residual sits between 9e-14 and 1.3e-12 with no trend in N, and **3 of its 24 default-tolerance solves exceed 1e-12** (N = 48 Dirichlet smooth 1.18e-12, N = 64 Neumann smooth 1.22e-12, N = 64 Dirichlet noisy 1.32e-12). All three are one V-cycle short of the next step: at `tol_rel` 5e-13 each takes one more cycle (11, 12, 11 instead of 10, 11, 10) and gives 8.9e-14, 1.0e-13 and 1.1e-13, no warning (`resid_mlmg_tol.log`); 2e-13 gives the same cycles. This is a granularity effect of a max-norm stopping rule against a 2-norm check (a factor of up to 1.3 overshoot), unrelated to N.
4. **The pin row is large and already excluded.** The full HYPRE Neumann residual (pin row included) is 1.1e-13, 1.0e-11, 3.0e-13, 1.4e-10, 1.0e-12, 7.3e-10 at N = 32, 48, 64, 96, 128, 160 (largest absolute max-norm entry 5.4e-5 at 160³, always the pin row), irregular in N, against the non-pin 9.8e-14 to 1.5e-12 smooth. The fixed 1e-12 against the full residual would warn at 48, 96, 128 (1.03e-12) and 160; the pin-aware check (already in the tip) removes that.

### 15.3 What the residual should be relative to: the constant

Three candidate scalings, each fitted as a constant over the six sizes (`resid_const.md`): `c_A = rel2_nopin / (u·||A||_2·||H||_2/||b||_2)` (round-off of a floating-point residual, relative to `||b||_2`), `c_N = rel2_nopin / (u·N²)` (size only), `c_S = relmax_nopin / (eps_mach·sqrt(N_cells))` (max-norm, relative to `max|b|`, the form `c·eps_mach·sqrt(N_cells)·|b|_inf`). Range over N, per series; "floor" = a run limited by round-off (FFT, HYPRE, MLMG at `tol_rel` 1e-15); "stop" = limited by the solver tolerance:

| Series | Type | c_A | c_N | c_S |
|---|---|---|---|---|
| FFT, Neumann, smooth | floor | 0.87 to 1.0 | 0.13 to 0.14 | 0.52 to 0.73 |
| FFT, Neumann, noisy | floor | 0.95 to 1.2 | 0.0004 to 0.008 | 0.006 to 0.047 |
| FFT, Dirichlet, smooth / noisy | floor | 1.0 to 1.2 | 0.0005 to 0.011 | 0.001 to 0.013 |
| HYPRE, Neumann, smooth (`tol_rel` 1e-12 and 1e-15) | floor | **3.4 to 5.6** | 0.52 to 0.86 | 1.9 to 3.8 |
| HYPRE, Dirichlet, smooth (`tol_rel` 1e-15) | floor | 1.3 to 1.4 | 0.001 to 0.013 | 0.003 to 0.024 |
| MLMG, Neumann, smooth (`tol_rel` 1e-15) | floor | 0.14 to 0.17 | 0.018 to 0.026 | 0.18 to 0.31 |
| MLMG, Dirichlet, smooth (`tol_rel` 1e-15) | floor | 0.27 to 1.3 | 0.0004 to 0.013 | 0.001 to 0.017 |
| MLMG, all, `tol_rel` 1e-12 | stop | 0.66 to 1,300 | 0.03 to 4.6 | 0.4 to 9.2 |
| HYPRE, noisy or Dirichlet, `tol_rel` 1e-12 | stop | 16 to 85 | 0.02 to 0.72 | 0.1 to 1.7 |

- `c_A` is the scaling that is flat in N (within 1.7 times in every floor series except MLMG Dirichlet, 0.27 to 1.3) and nearly flat across right-hand sides and boundary conditions: **at the round-off floor c_A is at most 1.4 for FFT, MLMG and HYPRE-Dirichlet and at most 5.6 for HYPRE-Neumann**. `c_N` varies by a factor of 2,000 between smooth and noisy data (it is only right for smooth data), `c_S` by a factor of 3,000.
- For runs limited by the solver tolerance the residual is whatever the stopping rule leaves (HYPRE noisy or Dirichlet at `tol_rel` 1e-12: 3e-14 to 9e-14, i.e. `c_A` 16 to 85, far below 1e-12); those are not floor-limited and never approach the check, except for the MLMG granularity cases in item 3.

### 15.4 Solution error against the residual

Neumann smooth: solution error `||H - H*||2/||H*||2` (mean removed)

| N | FFT | MLMG | HYPRE |
|---|---|---|---|
| 32 | 2.3e-15 | 1.7e-13 | 1.6e-14 |
| 48 | 6.1e-15 | 5.2e-14 | 9.9e-13 |
| 64 | 8.8e-15 | 2.3e-13 | 1.3e-14 |
| 96 | 1.6e-14 | 6e-15 | 5.2e-12 |
| 128 | 3.6e-14 | 2.3e-14 | 1.2e-14 |
| 160 | 8.5e-14 | 6e-15 | 1.4e-11 |

Dirichlet smooth: solution error

| N | FFT | MLMG | HYPRE |
|---|---|---|---|
| 32 | 4e-15 | 6.5e-13 | 3.6e-14 |
| 48 | 6.2e-15 | 2.5e-12 | 6.2e-14 |
| 64 | 6.4e-15 | 4.3e-13 | 5.9e-14 |
| 96 | 3.4e-14 | 1.1e-12 | 1.5e-13 |
| 128 | 4.5e-14 | 1.9e-12 | 1.3e-13 |
| 160 | 9e-14 | 2.5e-12 | 2.4e-13 |

Neumann noisy: solution error

| N | FFT | MLMG | HYPRE |
|---|---|---|---|
| 32 | 2.2e-15 | 2e-12 | 7.9e-14 |
| 48 | 5.7e-15 | 6.6e-13 | 1.2e-12 |
| 64 | 8.8e-15 | 2.3e-12 | 4.6e-13 |
| 96 | 1.5e-14 | 9.1e-13 | 5.4e-12 |
| 128 | 3.5e-14 | 3e-12 | 6.5e-13 |
| 160 | 8.3e-14 | 9.1e-13 | 1.6e-11 |

Dirichlet noisy: solution error

| N | FFT | MLMG | HYPRE |
|---|---|---|---|
| 32 | 3.8e-15 | 6.5e-13 | 4.6e-14 |
| 48 | 6.1e-15 | 2.4e-12 | 7.6e-14 |
| 64 | 6.3e-15 | 5.2e-12 | 8.1e-14 |
| 96 | 3.3e-14 | 1.1e-12 | 2.1e-13 |
| 128 | 4.4e-14 | 1.8e-12 | 1.8e-13 |
| 160 | 8.9e-14 | 2.5e-12 | 5e-13 |

- In all 96 solves the error is at most **5.2e-4 of eps_H** (smallest margin 1,930 times: MLMG Dirichlet noisy at N = 64, error 5.2e-12 against eps_H 1e-8). The error is not tied to the residual size: HYPRE Neumann has an error of 1.4e-11 at 160³ with residual 1.5e-12 (ratio 9), FFT has error 8.5e-14 with residual 3.9e-13 (ratio 0.2); over all runs error / residual is 0.004 to 200, which is conditioning and norm scatter, not growth with N.
- The residual **does** cross 1e-12 in 5 of the 96 solves (the three MLMG cases above and HYPRE Neumann 160³ at both tolerances, 1.52e-12 and 1.53e-12), and in all five the solution is far inside eps_H: errors 2.5e-12, 2.3e-13, 5.2e-12 and 1.4e-11 (twice) against eps_H of 1e-8 to 6.1e-8. A residual of 1.5e-12 relative to `||b||` therefore does not mean a solution error near eps_H.
- Negative controls (the check must still warn on a really under-converged solve), N = 96 Neumann smooth: MLMG at `tol_rel` 1e-8 gives 2.7e-9 and at 1e-10 gives 1.6e-11; HYPRE at 1e-8 gives 8.9e-10 and at 1e-10 gives 9.1e-12 (non-pin); errors 1.5e-10 and 9.4e-13 (MLMG), 5.5e-11 and 5.5e-12 (HYPRE). All four are 6 to 1,900 times the scaled limit of 1.4e-12 at that size (section 15.5) and still warn with it. The solution errors are 1.5e-10 or lower, still inside eps_H (2.2e-8), which shows that `residual_tol` is a much stricter requirement than A-56 needs.

### 15.5 Recommendation

Keep `residual_tol = 1e-12` as the base and make the limit size-aware through the quantity the common layer already computes:

    residual_limit = max( residual_tol , c · 2^-53 · ||A|| · ||H||_2 / ||b||_2 ),   c = 10,   ||A|| = 4 · sum_d (1/dx_d^2)

with `||H||_2` and `||b||_2` the weighted 2-norms of the returned solution and of the (mean-removed) right-hand side, and the same pin-excluded `rel2_nopin` as the compared quantity. This is today's informational `residual_floor` (it is already coded with `kResidualRoundoff = 10`), promoted from "reported for singular problems only" to "applied as a lower bound of the limit for every component, singular or not".
- **Constant:** measured `c_A` at the floor is at most 5.6 (HYPRE Neumann; at most 1.4 for every other backend and boundary condition), so c = 10 leaves a factor 1.8 above the largest measured value and 7 above the typical one. The measured flatness over N = 32 to 160 (within 1.7 times for FFT, HYPRE and MLMG-Neumann) is what justifies a constant that does not depend on N.
- **What it does in numbers:** for the smooth data used here `||A||·||H||/||b|| = 157·(N/32)²`, so the limit is max(1e-12, 1.7e-16·N²): 1.7e-13 at 32, 1.4e-12 at 96 (the base 1e-12 applies up to about N = 77), 2.5e-12 at 128 and 3.7e-12 at 160 (the measured values of the floor; 1.7e-16·N² is a slightly high fit). With it none of the 96 solves exceeds the limit except the three MLMG granularity cases of item 3 (the two HYPRE 160³ cases pass: 1.52e-12 against 3.7e-12). It is not a size-only formula; if a closed form in N alone is wanted for documents, `2·2^-53·N²` (N = cells per direction of the finest level; `N_cells^(2/3)` for non-cubic) bounds the floor-limited residual of all backends on smooth data (largest `c_N` 0.86), but it overstates the limit for noisy data by up to 1,000 times and should not be used as the code's check.
- **Not recommended:** the max-norm form `c·eps_mach·sqrt(N_cells)·|b|_inf` has no usable constant (0.001 to 3.8 at the floor depending on the data; 9.2 for MLMG at its stopping point), and a plain relative-to-`||b||_2` threshold fixed in N is the present behaviour that reaches 1e-12 at 160³.
- **MLMG granularity:** independently of the size scaling, either run MLMG with `tol_rel` 5e-13 (one more V-cycle in the three affected cases, none lost elsewhere in this set; not measured across all cases) or accept the base 1.5e-12 (1.14 times the largest overshoot seen, 1.32e-12 in 24 default-tolerance solves). The first removes the spurious warnings at their source and is the better choice if the solve time allows.

### 15.6 What changes in the pass/fail rules

1. **A-56 solution criterion unchanged:** `eps_H = max(1e-8, 2.4e-12·N²)` on the relative L2 difference. All 96 solves above pass it by at least 1,930 times, including those whose residual exceeds 1e-12, so the residual criterion must not be turned into a solution failure.
2. **Residual criterion (section 5 and the residual rows of every case table):** replace "no residual warning at `residual_tol` = 1e-12" by "`residual_rel2_nopin <= residual_limit`" with the limit of 15.5; a case with `residual_rel2_nopin` between 1e-12 and `residual_limit` is a pass, a case above it is a warning that must be explained. The negative controls (loose `tol_rel`) must still produce the warning; they do (15.4).
3. **Reporting:** every case table adds `residual_floor` and `residual_backward` next to the residual (they are in `PressureResult` already). A backward error that grows with N for FFT or HYPRE would be a real loss of accuracy independent of the data scale; none occurred in this set (log-log slopes at most 0.34; HYPRE Neumann smooth 4.0e-16 to 6.2e-16, HYPRE Dirichlet 1.7e-15 to 8.5e-15, MLMG 7e-17 to 1.3e-13 and tolerance-limited).
4. **Code change (not made here, for the backend owner):** `evaluate_residual` in the common layer compares `residual_check` with `o.residual_tol`; the proposed rule compares with `max(o.residual_tol, residual_floor)` and computes the floor for non-singular components as well (today it is 0 there). `PressureResult::residual_limit` should then report the effective limit.

### 15.7 Limits

One right-hand-side family per spectrum type (two modes plus a constant; the same plus noise), unit cube, uniform grid, single level, 4 ranks; no composite or masked case. HYPRE was run with its default setup; the identity pin sits at a fixed cell. The conclusion for composite problems (where `||A||` is taken at the finest level) is not measured. The floor constant was measured on six sizes up to 160³; it is flat there, but extrapolation beyond 160³ is not tested. The Neumann series at `tol_rel` 1e-15 for MLMG stops at the 200-cycle cap (status `NotConverged`) and is used only as a floor measurement, not as a recommended setting.

## 16. P3-R02, second case family: noisy velocity, three levels, patch moved between projections (measured)

Question: does the bound of section 13 (`max|div u - D - c| <= 10·eps_rel·B + 20·eps_mach·U/dx_fine`, working bound 1e-9·U/dx_fine at eps_rel = 1e-12) still hold when the velocity is not smooth, the hierarchy has three levels, the refined patches move between projections so that a regrid seam is exercised, the cell source D is nonuniform with a nonzero mean, and the constant c is nonzero?

Setup (product tip `99c55244e1`, D-067 defaults, one rank, one box per level; driver `a56_regrid.cpp` in the scratch directory; raw results and the table generator in the scratch results directory, `rg_matrix.log`, `rg_table.md`):
- **Hierarchy A.** Level 1 is the middle half of the unit cube (in level-0 cells), level 2 is the middle half of level 1; refinement ratio 2 (also ratio 4 and two levels in the variations). Face velocity on every level = smooth field (zero normal component at closed walls) plus **grid-scale noise**: a hash-based value in [-0.5, 0.5] times `noise·U` on every face (noise 0.3 in the main runs, 1.0 and 0 in variations; the noise is not put on closed-wall faces). The noise gives `max|div u - D|` of 102 on A for n = 32, three levels, i.e. 0.8 times U/dx_fine. D = `U·(dmean + 0.5 sin(3πx) cos(2πy) sin(πz) + 0.2 cos(5πx) sin(4πz))` evaluated analytically on each level (nonuniform; mean dmean·U plus the oscillating part), dmean = 0.3 in the main runs. Project on A (`solve_pressure` + `face_gradient_composite`, `u = u* - grad H`), then average the fine faces down to the coarse faces.
- **Regrid to hierarchy B.** Level 1 patch shifted by `shift` level-0 cells in every direction (n/8 in the main runs, i.e. a quarter of the patch width; variations 0, 1 and 6 at n = 32), deeper levels follow as the middle half of the shifted level. Faces that exist in A are copied; faces of cells that are new on a fine level are interpolated trilinearly from the coarser level of B (face-normal linear in the normal direction, bilinear in the tangential ones); coarse faces under fine cells are the average of the fine faces. D is re-evaluated analytically on B.
- **Reprojection on B** with the same solve, then the observable `max|div u - D - c|` over the uncovered cells of B. c is the volume mean of `b = div u - D` for the closed box (the removed constant); c = 0 when any face is Dirichlet. B in the bound is `max|b|` on B before the reprojection. "Seam before reprojection" is `max|div u - D - c|` on B right after the regrid and before reprojecting; it doubles as the negative control (no projection).
- eps_rel = 1e-12, 1e-9 and 1e-6; U/dx_fine is 128 for n = 32 with three levels of ratio 2 (finest dx = 1/128) and for two levels of ratio 4, 256 for n = 64, 384 for n = 96, 512 for n = 32 with three levels of ratio 4.

Results (error divided by U/dx_fine, then margins = bound / measured error, for the derived bound and for the working bound 1e-9·U/dx_fine; a margin below 1 is a failure):

| n | levels | ratio | faces | shift | noise | D mean | U | seam before reprojection | B | c | err/(U/dx) at 1e-12 | margin derived / working | err/(U/dx) at 1e-9 | margin derived / working | err/(U/dx) at 1e-6 | margin derived / working |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | 3 | 2 | `NN,NN,NN` | 4 | 0.3 | 0.3 | 1 | 25 | 25 | -0.3 | 1.3e-13 | 16 / 8e+03 | 4.6e-11 | 43 / 22 | 2.3e-08 | 83 / 0.043 |
| 32 | 3 | 2 | `DD,DD,DD` | 4 | 0.3 | 0.3 | 1 | 25 | 25 | 0 | 1.1e-13 | 18 / 9.1e+03 | 4.4e-11 | 44 / 23 | 2.3e-08 | 84 / 0.044 |
| 32 | 3 | 2 | `ND,DN,NN` | 4 | 0.3 | 0.3 | 1 | 25 | 25 | 0 | 1.1e-13 | 17 / 9e+03 | 4.4e-11 | 43 / 23 | 2.3e-08 | 82 / 0.043 |
| 64 | 3 | 2 | `NN,NN,NN` | 8 | 0.3 | 0.3 | 1 | 59 | 59 | -0.3 | 9.8e-14 | 24 / 1e+04 | 4.8e-11 | 48 / 21 | 2.7e-08 | 85 / 0.037 |
| 64 | 3 | 2 | `DD,DD,DD` | 8 | 0.3 | 0.3 | 1 | 59 | 59 | 0 | 1e-13 | 23 / 9.9e+03 | 4.8e-11 | 48 / 21 | 2.7e-08 | 84 / 0.036 |
| 64 | 3 | 2 | `ND,DN,NN` | 8 | 0.3 | 0.3 | 1 | 59 | 59 | 0 | 1e-13 | 23 / 9.9e+03 | 4.8e-11 | 48 / 21 | 2.7e-08 | 84 / 0.036 |
| 32 | 3 | 2 | `NN,NN,NN` | 4 | 0 | 0.3 | 1 | 0.05 | 0.35 | -0.3 | 3.3e-16 | 96 / 3e+06 | 2.7e-13 | 1e+02 / 3.7e+03 | 5e-11 | 5.4e+02 / 20 |
| 32 | 3 | 2 | `NN,NN,NN` | 4 | 1 | 0.3 | 1 | 82 | 82 | -0.3 | 3.9e-13 | 16 / 2.5e+03 | 1.5e-10 | 42 / 6.6 | 7.8e-08 | 82 / 0.013 |
| 32 | 3 | 2 | `NN,NN,NN` | 4 | 0.3 | 0 | 1 | 25 | 25 | 1e-05 | 1.3e-13 | 15 / 8e+03 | 4.6e-11 | 42 / 22 | 2.3e-08 | 82 / 0.043 |
| 32 | 3 | 2 | `NN,NN,NN` | 4 | 0.3 | 1 | 1 | 25 | 26 | -1 | 1.2e-13 | 16 / 8.1e+03 | 4.6e-11 | 44 / 22 | 2.3e-08 | 85 / 0.043 |
| 32 | 3 | 2 | `NN,NN,NN` | 4 | 0.3 | 0.3 | 100 | 2.5e+03 | 2.5e+03 | -30 | 1.3e-13 | 16 / 8e+03 | 4.6e-11 | 43 / 22 | 2.3e-08 | 83 / 0.043 |
| 32 | 3 | 2 | `NN,NN,NN` | 0 | 0.3 | 0.3 | 1 | 3.6e-11 | 0.3 | -0.3 | 1.9e-15 | 15 / 5.4e+05 | 1.1e-16 | 2.1e+05 / 8.8e+06 | 2.2e-14 | 1.1e+06 / 4.6e+04 |
| 32 | 3 | 2 | `NN,NN,NN` | 1 | 0.3 | 0.3 | 1 | 25 | 25 | -0.3 | 1.5e-13 | 13 / 6.7e+03 | 4.4e-11 | 45 / 23 | 2.4e-08 | 80 / 0.041 |
| 32 | 3 | 2 | `NN,NN,NN` | 6 | 0.3 | 0.3 | 1 | 23 | 23 | -0.3 | 2.1e-14 | 85 / 4.7e+04 | 1.8e-10 | 10 / 5.7 | 8.5e-08 | 21 / 0.012 |
| 32 | 2 | 2 | `NN,NN,NN` | 4 | 0.3 | 0.3 | 1 | 14 | 14 | -0.3 | 2.2e-14 | 1e+02 / 4.6e+04 | 5.4e-11 | 41 / 19 | 2.6e-08 | 84 / 0.039 |
| 32 | 2 | 4 | `NN,NN,NN` | 4 | 0.3 | 0.3 | 1 | 29 | 29 | -0.3 | 6.4e-14 | 36 / 1.6e+04 | 1.4e-11 | 1.6e+02 / 72 | 9.9e-08 | 23 / 0.01 |
| 32 | 3 | 4 | `NN,NN,NN` | 4 | 0.3 | 0.3 | 1 | 1.1e+02 | 1.1e+02 | -0.3 | 7.4e-14 | 29 / 1.3e+04 | 2.8e-11 | 77 / 36 | 1.2e-08 | 1.9e+02 / 0.086 |
| 32 | 2 | 2 | `DD,DD,DD` | 4 | 0.3 | 0.3 | 1 | 14 | 14 | 0 | 4.6e-14 | 47 / 2.2e+04 | 2.8e-11 | 76 / 36 | 2e-08 | 1e+02 / 0.049 |
| 96 | 3 | 2 | `NN,NN,NN` | 12 | 0.3 | 0.3 | 1 | 85 | 85 | -0.3 | 1.6e-13 | 14 / 6.4e+03 | 6.3e-11 | 35 / 16 | 3.2e-08 | 69 / 0.031 |

Findings:
1. **At eps_rel = 1e-12 every case passes both bounds.** The error is 1.0e-13 to 1.6e-13 times U/dx_fine in the main cases (2.1e-14 to 3.9e-13 over all moved-patch variations; 3.9e-13 is the noise = 1.0 case), which is 2,500 to 47,000 times below the working bound, and the margin to the derived bound is 13 to 102 for the moved-patch runs (tightest 13 for shift 1, 14 for n = 96, 15 to 16 for the main closed-box cases). The tightest margins match family 1 (11 to 53).
2. **At eps_rel = 1e-9 every case passes the derived bound (margins 10 to 160) and the working bound (margins 5.7 to 72; the two smallest are 5.7 for shift 6 and 6.6 for noise = 1.0).** The working bound is thus not violated at 1e-9 here but is tighter than in family 1; the most stressed case (shift 6, a patch moved by 0.19 of the domain) has an error of 1.8e-10 times U/dx_fine, which is 0.97 of `eps_rel·B`.
3. **At eps_rel = 1e-6 the working bound fails in every case with a real seam** (margins 0.01 to 0.09; the noise-free case, with a seam of only 0.05, still passes at 20) while the derived bound still holds (margins 21 to 540), as in family 1: the working bound is valid for solver tolerances of about 1e-9 and tighter.
4. **The seam does not break the bound.** The seam before reprojection is 14 to 111 (U = 1; 0.1 to 0.9 times U/dx_fine) in every moved-patch case, about 1e11 times the derived bound at 1e-12 (negative control), and the projection removes it down to the level of the solver tolerance. The error is not concentrated at the coarse-fine interface: the maximum over the interface cells (within one cell of a coarse-fine boundary) and the maximum over all other cells differ by less than a factor of 2.4 in every moved-patch case (interior maximum / interface maximum 0.59 to 2.3; for the closed box the interface cells are the larger by up to 1.1 times at 1e-12, for sets with Dirichlet faces the interior is larger). The static patch (shift 0) starts with a seam of only 3.6e-11 because the fields are already a projected field, and gives 1.9e-15 times U/dx_fine at 1e-12; it is not a seam test.
5. **Nonzero mean c.** c = -0.3 (D mean 0.3 U), -1.0 (D mean U) and -30 (U = 100) in the closed box, 1.0e-5 (D mean 0, the oscillating part only); the error does not depend on c (1.2e-13 to 1.3e-13 times U/dx_fine for c = -1, -0.3, 0 at 1e-12) and scales with U (U = 100: error 100 times larger, same ratio to U/dx_fine). With open faces (`DD,DD,DD`, `ND,DN,NN`) c = 0 and the results match the closed box within 15%.
6. **Floor term.** The noise = 0 case has B = 0.35 and a seam of only 0.05, so the machine-precision term of the bound (20·eps_mach·U/dx_fine = 5.6e-13) is 14% of the bound at 1e-12 (4.1e-12); the measured error is 4.2e-14, margin 96. The floor term exceeds the first term only when B < 2·(eps_mach/eps_rel)·U/dx_fine, i.e. B below 0.056 for U/dx_fine = 128 at 1e-12; no case here is below that, so the floor term was not stressed.
7. **Levels and ratios.** Two levels, ratio 2: 2.2e-14 times U/dx_fine at 1e-12 (margin 102); two levels, ratio 4: 6.4e-14 (36); three levels, ratio 4: 7.4e-14 (29); no degradation with ratio or with depth (three levels, ratio 4: margin 29; three levels, ratio 2: 16). Errors grow only with the tolerance times B.
8. **Size.** n = 32, 64 and 96 with three levels, ratio 2, closed box, noise 0.3: error 1.3e-13, 9.8e-14, 1.6e-13 times U/dx_fine, margins 16, 24, 14 to the derived bound, 6,400 to 10,000 to the working bound: no trend with size in this range.

Limits: single rank only; the regrid is a hand-written interpolation of the face velocity (trilinear from the coarser level), not the AMR code's own interpolator, and the velocity is a synthetic smooth-plus-noise field, not a flow. The bound is a statement about the projection, so the pressure-solve error is what is tested; whether a real regrid in the code produces seams larger than the 0.9 U/dx_fine of this family is not covered. The noise amplitude (up to 1.0) and the patch shift (up to 0.19 of the domain) were pushed until the margin shrank (5.7 at 1e-9) but not until the derived bound failed.
