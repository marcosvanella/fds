# 06. Mean removal, solution gauge, stretched cells, masked components, and the D-057 2-D review

**Status: complete for the items below; open items are listed in section 10.** Everything marked "measured" was run; everything marked "read" is from source only; "unverified" is neither.
FDS `file:line` cites are against FireX `36975d765f` (`src/Source/pres.f90`); the same blocks are in master `ce1f659cd4` at the lines given in patch 0005. AMReX is 99ddfda (26.09), HYPRE v2.32.0-24.
Timings: CPU runs on an owner-provided NVIDIA test machine (pinned P-cores, 1 thread per rank) and its GPU (1 rank, managed arena, HYPRE on the device for H); stretched accuracy runs at 64³ are not timing-grade.

## 1. Answers in short

1. **Study verdict (backend study, commit `6f414db72f`): confirmed**, with one nuance (section 2). FDS removes the arithmetic mean of the volume-scaled RHS `F = V·b`; on stretched cells that is *not* the volume-weighted mean of `b`; both give a compatible system, they pick different compatible RHS when the input is not compatible, and they agree to round-off when it is.
2. **A real FDS defect was found and patched**: `GLMAT_SOLVER` does not remove the mean of `F_H` (or of `X_H` in the periodic-test-7 fallback) on one MPI rank. Patch `docs/upstream-patches/0005-glmat-singlerank-mean-removal.patch` (section 4).
3. **Gauge recommendation** for the Architect: section 5.
4. **ADR wording** (draft): section 6.
5. **Stretched and masked runs**: sections 7 and 8. **D-057 review**: section 9.

## 2. Review of the mean-removal study (FDS source, read)

| Item | FDS fact | Cite |
|---|---|---|
| ULMAT RHS | `F_H = PRHS·DX·DY·DZ` (volume-scaled) | `pres.f90` ≈1497-1503 |
| ULMAT mean | arithmetic, count-weighted mean of the scaled `F_H`, no volume weight | 1683-1733 (`H_INDEFINITE_IF_1`) |
| ULMAT pin | `F_H(NUNKH)=0` | 1748 |
| ULMAT solution mean / gauge | solution-mean removal, then `SHIFT = Σ(V·ρ·(KRES+X_H)) / (Σ(V·ρ)+20ε)` | 1779-1828, 1830-1878 |
| Matrix and sign | entries `AF/DX1`; `K x = F`, `H = −x` | ≈2585-2610 |
| GLMAT whole-domain mean | `SUM_FH(2)` is set only inside `IF (N_MPI_PROCESSES>1)`; on one rank `MEAN_FH = 0` | 3370, 3401-3406 |
| GLMAT X-mean fallback | same pattern with `SUM_XH(2)` | 3550-3556 |
| GLMAT zone branch | fine (element 2 holds the cell count) | 3371-3400 |

**Operator and compatibility per backend (read, and checked by the dense script below):**

| Backend | Operator it solves | Stretched cells | Native singular handling |
|---|---|---|---|
| FDS ULMAT/GLMAT/UGLMAT | scaled `K` (symmetric, `1ᵀK=0`) | yes | arithmetic mean of `F`, pin |
| MLMG (`MLPoisson`, `MLABecLaplacian`) | uniform `dx` from `Geometry` | only via `MLABecLaplacian` on unit spacing with `beta = AF/DX1` (what the harness does, "Mf") | `makeSolvable` (`AMReX_MLMG.H:2355`, called at 1555 when `isSingular(0)` and `enforceSingularSolvable`, default true, `AMReX_MLLinOp.H:888`) subtracts the cell-count mean of the RHS (`AMReX_MLCellLinOp.H:1899-1932`). With a pin cell (`a>0`) `isSingular(0)=0` and it does not run |
| Assembled HYPRE | `K` assembled directly | yes | none (caller removes the mean and pins) |
| `FFT::Poisson` | unscaled physical Laplacian | **no**, uniform cells only | zero mode dropped |

So on a unit-spacing scaled operator MLMG's native removal equals FDS's; on the unscaled `b` with different `dx` per level it is not volume-weighted, and the common layer must remove the mean first (D-032, `docs/README.md` line 78, deliberately uses the composite volume-weighted mean).

**Dense numerical check** (measured; `scratch/pressure-signoff/meanremoval/check_meanremoval.py`, `.out`; three small stretched cases, 270 and 40 unknowns, V max/min 4.95 to 5.78): `K` is symmetric, `1ᵀK = 0` (4e-16), so `K x = F` is solvable exactly when `ΣF = 0`; the unscaled operator `L = −V⁻¹K` has left null vector `V`, i.e. the same condition `Σ V b = 0`.
- Removing the arithmetic mean of the *unscaled* `b`: residual 5.6 (3-D cases) and 2.2 (1-D-like): **not compatible**.
- FDS (arithmetic mean of `F`) and volume-weighted (`F − cV`, `c = ΣF/ΣV`): residuals 1e-14 to 1e-13, both compatible.
- On an incompatible RHS the two solutions differ by 0.70 (3-D) and 0.94 (1-D-like) (volume-weighted relative L2). On a compatible RHS (`|mean F|/rms F` 6e-17 to 2e-16) they agree to 3e-17 to 3e-16.
- The final `ρV` gauge removes the pin value (max |g(x)−g(x+3.7)| 5e-16) and leaves `Σ ρV (KRES−H)` at 1e-16 or less.

**Nuance:** "volume-weighted mean is wrong versus FDS" is true for FDS *parity*, not for compatibility. Both are valid; they differ in how an incompatible part of the RHS is attributed (a constant per cell versus a constant per unit volume).

## 3. Where the incompatibility comes from, and which removal is right

If the imbalance is spread per unit volume (a net dilatation error, an inconsistent zone sum), the volume-weighted removal removes exactly it; the arithmetic removal moves a different, cell-count-uniform amount and the solution is wrong by an amount that grows with stretching (section 7). If the imbalance is truly per cell, the opposite holds. FDS's own RHS is compatible up to its zone-sum correction, so in normal runs both are within round-off of each other. The recommendation (section 5) is therefore: default to the volume-weighted removal, keep arithmetic as an FDS-parity option.

## 4. GLMAT single-rank defect and patch 0005

**Defect (read, confirmed by runs):** on one rank `MEAN_FH` and `MEAN_XH` are 0 in the whole-domain branch (section 2 table). With 2 or more ranks the `ALLREDUCE` fills element 2 and the mean is removed.

**Patch:** `docs/upstream-patches/0005-glmat-singlerank-mean-removal.patch`, two added lines (`SUM_FH(2) = SUM_FH(1)`, `SUM_XH(2) = SUM_XH(1)`), target both lines. `git apply --check` passes on FireX `36975d765f` and master `ce1f659cd4` (offsets -63 and -107). Index row added in `docs/upstream-patches/README.md`. The X-mean fix the Chief asked for (pres.f90 ≈3553-3555) is in the same patch.

**Behavior-unchanged check (what was run; full text in the patch file):**

| Test | Result |
|---|---|
| Standard `patch_check.sh` release set plus an extra GLMAT case on 2 ranks (`gl_closed2`, two meshes) | **PASS, all five bitwise identical** (23, 51, 91, 41, 19 files). The multi-rank path is untouched: the added statement is overwritten by the `ALLREDUCE` |
| 1 rank, mean already exactly zero (scratch hook sets the last `F_H` entry so the sequential sum is exactly 0.0; checked 0.0 in all four solves) | **byte-identical**: dumps, restart, s3d, hrr, steps, `.out`/`.smv` (11 identical, 4 identical after timing-strip; only captured stdout differs, it has wall-clock lines) |
| 1 rank, real compatible RHS (`gl_closed1`) | **not bitwise**: the real sum is rounding-level, never exactly 0. Unpatched removes exactly 0, patched removes −1.3e-16 of the F rms. First solve: `F_H` differs ≤ 2.1e-16 rms, final H ≤ 7.4e-15 of max\|H\|; restart file differs at round-off (heuristic max relative word difference 5.8e-10). Declare as a round-off change |
| 1 rank, incompatible RHS (hook adds 0.3·V; mean = 7.3e-2 of rms; N=1920) | **corrected**: unpatched leaves Σ F at 7.3e-2 (minimum possible residual 3.2 rms; H off by 0.73 relative L2 from the dense reference); patched leaves 4e-16 and matches the dense arithmetic-removal reference to 4.2e-14 (`patchwork/ref_check.py`, `.out`) |
| X-mean branch (`gl_shunn3`, 32×1×32, periodic test 7, 1 rank) | **declared behaviour change**: patched − unpatched H is a pure constant (spatial variation 1.3e-14). Unpatched mean(H) = +0.0304 (the pin gauge), patched −3e-17. Over the 0.05 s run the H slice differs by that constant and the mass-fraction slice by up to 1.5e-5 (float32 view), temperature 9.6e-8, hrr `Q_PRES` 9.3e-6 of the column norm |

Not run: master build and run (only `git apply --check`), a `-fcheck=all` build, GLMAT with cut cells or several zones on one rank.

## 5. Gauge recommendation (for the Architect)

1. **Default gauge: `Σ ρ·V·(KRES − H) = 0` per pressure zone and per connected component**, computed with the exact (decomposition-independent) sum over uncovered cells. This is the UGLMAT/ULMAT convention (`SHIFT`, pres.f90 1830-1878).
2. **Always apply it, even when the solve is "exact up to a constant"**: `p = ρ·(H − KRES)` carries the constant into the baroclinic pass-2 `PRHS` (variable ρ makes `∇(ρc) = c∇ρ ≠ 0`), so the gauge is not cosmetic.
3. **Evidence that the constant matters:** in a one-rank periodic-test-7 run the only non-rounding difference between unpatched and patched GLMAT is the H constant (+0.0304 versus 0), and it changed the mass-fraction slice by up to 1.5e-5 in 0.05 s (section 4). The same case on 2 ranks already had zero mean, so FDS currently gives rank-count-dependent gauges.
4. **Mean removal default: the composite volume-weighted mean of `b`** (`F − cV`, D-032). It equals FDS's arithmetic removal on uniform cells, agrees with it to round-off on any compatible RHS (dense check 3e-17 to 3e-16), and corrects a per-volume imbalance (offset 0.3 at 64³, R=8: error 9.85e-4 versus 3.42e-2 with the arithmetic removal; at 1M cells 4.02e-4 versus 3.40e-2).
5. **FDS-parity option:** arithmetic removal of the scaled `F` (runtime switch). Do not rely on native `makeSolvable` (count-weighted, skipped when a pin is present).
6. **Open for the Architect:** FFT-solved FDS cases have an arbitrary zero-mode constant (baseline `ns2d_16` H mean ≈ 1.88). Whether the AMR path should reproduce that or apply item 1 also to FFT cases is a decision; item 1 is recommended, and PRES comparisons to the FFT baseline are then made after removing a constant (as the test plan already does).

## 6. Draft ADR wording

> Pressure mean removal and gauge. The backend-independent common layer removes the incompatible part of the RHS per connected component by subtracting a constant per unit volume (`F − cV`, the composite volume-weighted mean, D-032) with decomposition-independent exact sums over uncovered cells; on uniform cells this equals FDS's arithmetic removal, and an FDS-parity arithmetic option is kept. After the solve the gauge is fixed by `Σ ρ V (KRES − H) = 0` per zone and component and is always applied, because `p = ρ(H − KRES)` feeds the baroclinic second pass. Native MLMG `makeSolvable` is not relied on.

## 7. Stretched-grid runs (measured)

Harness `scratch/pressure-signoff/pressure_1M/src/h2h_str.cpp`; raw results `results_pr06/`. System as FDS builds it: `K x = F`, `K = Σ (AF/DX1)(x_i − x_j)`, `F = V·b`, unit-spacing geometry. z widths geometric with last/first ratio `R` (`R=1` uniform). `case=dir`: `u = sin·sin·sin`, all faces Dirichlet (no mean). `case=neu`: `u = s(x)s(y)s(z)`, `s(t) = cos(π t²)`, zero flux, pin on the last cell; the solution is not symmetric, so the midpoint sum of `F` is O(h²) from zero (a naturally slightly incompatible RHS). `offset=a` adds `a·V` to `F` (a per-unit-volume imbalance). `mean=arith` removes the mean of `F` (FDS), `vol` removes `cV`, `none` removes nothing. Backends: Mf = MLMG on `MLABecLaplacian` with HYPRE bottom (FDS BoomerAMG options), H = assembled HYPRE PCG+BoomerAMG, FFT only for `R=1`. Error = V-weighted relative L2 against `u_ex` up to the volume-weighted constant. Tolerance 1e-10.

**64³, 4 ranks, CPU** (Mf and H give identical error to the printed digits in every row):

| case | R | offset | mean | error | Mf iters | H iters |
|---|---|---|---|---|---|---|
| neu | 1 | 0 | arith (= vol) | 3.67e-4 (FFT, Mf, H) | 3 | 22 |
| neu | 1 | 0.3 | arith / vol | 3.67e-4 / 3.67e-4 | 3 | |
| neu | 8 | 0 | arith | 9.77e-4 | 3 | 22 |
| neu | 8 | 0 | vol | 9.85e-4 | 3 | 22 |
| neu | 8 | 0 | none | 1.26e-3 | 3 | 22 |
| neu | 8 | 0.3 | arith | **3.42e-2** | 3 | 22 |
| neu | 8 | 0.3 | vol | **9.85e-4** | 3 | 22 |
| neu | 8 | 0.3 | none | 0.247 | 3 | 22 |
| neu | 64 | 0 | arith / vol / none | 2.84e-3 / 2.85e-3 / 3.24e-3 | 5 | 27 |
| neu | 64 | 0.3 | arith / vol | **5.58e-2 / 2.85e-3** | 6 / 5 | 27 |
| dir | 1 | | | 2.01e-4 (FFT, Mf 29 it, H 19 it) | | |
| dir | 8 | | | 2.74e-4 (Mf 44 it, H 19 it) | | |
| dir | 64 | | | 5.92e-4 (Mf 117 it, H 19 it) | | |

**About 1M cells (100³, R=8), error and solve-only median per solve** (CPU 8 ranks, 10 timed solves; GPU 1 rank, 20 timed solves):

| case | mean | error | CPU Mf s (it) | CPU H s (it) | GPU Mf s (it) | GPU H s (it) |
|---|---|---|---|---|---|---|
| neu, offset 0 | arith | 3.99e-4 | 1.47 (4) | 0.472 (27) | 0.274 (3) | 0.111 (24) |
| neu, offset 0 | vol | 4.02e-4 | 1.48 (4) | 0.465 (27) | 0.274 (3) | 0.111 (24) |
| neu, offset 0.3 | arith | **3.40e-2** | 1.11 (3) | 0.472 (27) | 0.274 (3) | 0.111 (24) |
| neu, offset 0.3 | vol | **4.02e-4** | 1.49 (4) | 0.473 (27) | 0.274 (3) | 0.111 (24) |
| dir | n/a | 1.12e-4 | 0.845 (32) | 0.372 (21) | 0.396 (54) | 0.097 (21) |
| neu, R=1 | arith | 1.50e-4 | 0.881 (3) | 0.439 (26) | 0.243 (3) | 0.109 (24) |
| neu, R=1, FFT | arith | 1.50e-4 | 0.0044 (1) | | 0.0087 (1) | |

Findings: (a) with no imbalance the two removals differ by under 1% of the error; (b) with a per-volume imbalance the volume-weighted removal restores the discretisation-level error and the arithmetic one is 20 to 85 times worse (R=64 at 64³: 20; R=8 at 64³: 35; R=8 at 1M: 85; the gap widens as the grid is refined because the discretisation error falls while the imbalance error does not); (c) CPU and GPU, Mf and H agree on the error; (d) `none` is 1.3× worse with no imbalance and 250× worse at offset 0.3; the post-removal sums are checked per run (`ΣF/Σ|F|` ≤ 2e-14 for arith and vol); (e) at `R=64` the MLMG route is much slower than H (Mf 2.3 to 7.2 s per 64³ solve against H 0.14 s; at `R=8` 0.2 to 0.3 s against 0.11 to 0.14 s); cause not investigated; (f) FFT is two orders faster but legal only at `R=1`.

## 8. Masked stairwell, three components (measured)

Stairwell union (11 meshes, 192×176×552 padded, `dx` = 0.1, `ba=drop`, 1,637,828 gas cells, one connected component as built). `cutk=60 190` turns two z-planes into one-cell gaps, giving **three components: one OPEN (545,550 cells, vent on the low-y face) and two sealed (204,594 and 884,512 cells)**, 1,634,654 unknowns. The harness runs a union-find check (`COMPCHECK`) that agrees with the flood fill, shows one pin per sealed component and none in the open one. RHS synthetic, zero mean per sealed component, uniform `dx` (so arithmetic and volume-weighted removals coincide here). Tolerance 1e-10; no warm-up; solve-only median over the stated number of solves; max_grid_size 32 (CPU) / 64 (GPU).

| Backend | Ranks | Iters | Setup s | Re-solve median s (3 repeats) | Min-max s | Spread | Timed solves |
|---|---|---|---|---|---|---|---|
| Mf | 1 CPU | 3 | 4.2 | 5.785 | 5.704-5.864 | 2.8% | 10 |
| H | 1 CPU | 27 | 3.25 | 2.216 | 2.214-2.258 | 2.0% | 20 |
| Mf | 4 CPU | 3 | 1.58 | 2.781 | 2.779-2.992 | 7.6% | 20 |
| H | 4 CPU | 25 | 1.46 | 1.09 | 1.084-1.166 | 7.5% | 40 |
| Mf | 8 CPU | 3 | 0.632 | 2.7 | 2.538-2.768 | 8.5% | 20 |
| H | 8 CPU | 26 | 0.989 | 1.036 | 0.9443-1.042 | 9.4% | 40 |
| Mf | 1 GPU | 3 | 0.667 | 1.191 | 1.191-1.191 | 0.1% | 30 |
| H | 1 GPU | 26 | 0.223 | 0.1934 | 0.1934-0.1936 | 0.1% | 50 |

3 repeats per row, interleaved, on a quiet machine (protocol and the Mb rows in doc 07 §7.5); no run was rejected. Replaces the single-repeat table of the first version (the GPU medians agree within 1%; the CPU medians differ by 5 to 18% because the first version ran at lower, uncontrolled core clocks, e.g. one-rank Mf 6.52 s, now 5.79 s). Setup is the median set-up time of the repeats.

Single open component (no cut, same geometry): Mf 8 CPU 2.63 s, H 8 CPU 0.988 s, Mf GPU 1.29 s, H GPU 0.209 s. Splitting into three components costs nothing measurable (−2% and −4% on CPU, within the p90 spread).

Per-component true residuals (relative L2, one rank run; all three components, both backends): Mf 6.6e-12 (c1) to 2.8e-11 (c2); H 1.3e-11 to 6.9e-11; the 8-rank H run reaches 1.5e-10 on the large sealed component (requested 1e-10 is global). The post-removal `ΣF/Σ|F|` per sealed component is 1e-14 or less in every run. **Mf and H agree per component** (one-rank dump, constant removed per component): relative L2 8.4e-12 (open), 5.8e-13 and 1.2e-11 (sealed).

Caveats: Mf on a pinned singular component builds one MG level (the pin blocks coarsening), so an MLMG iteration is a full HYPRE solve, as in doc 05. In the repeat runs the package reached at most 97 °C (one short setup spike), never 30 s at or above 95 °C, and the held core clocks were 5.2 GHz (1 rank), 4.2 GHz (4) and 3.5 GHz (8). The Mf-versus-H gap is not a clock artefact.

## 9. Task 4: D-057 2-D (single-y-cell) review from the pressure side

D-057 (`docs/README.md` line 103, `requirements.md` line 183): all-periodic and all-Neumann use `FFT::Poisson`; a one-cell direction is ignored by the FFT solver; the backend owner checks the BC mapping including mean removal. Sources: AMReX `AMReX_FFT_Poisson.H`, `AMReX_MLLinOp.H`, `AMReX_MLPoisson.H`, `AMReX_MLABecLaplacian.H`; the driver note `src/Source/driver/notes/fft-thin-direction-check.md`; the backend `pressure_backend/` and its `frozen/composite-notes.md`; my harness (`h2h_str ny=1`, n=64, 2 ranks).

| Question | Finding | Status |
|---|---|---|
| FFT ignores a one-cell direction | `dxfac[idim]=0` for every direction with `length==1`, whatever the BC (`AMReX_FFT_Poisson.H:317-324`). Matches FDS: `JBAR==1` sets `TWO_D` (`read.f90:703`) and the 2-D FFT call `H2CZSS` has no y term (`pres.f90:345-348`). Run: 64×1×64 Neumann, FFT true residual 1.5e-13 against an assembled operator with no y term; MLMG and HYPRE agree with it (error against the exact field identical to 9 digits) | **confirmed** (read + run) |
| FFT with a Dirichlet face in the one-cell direction | the term is dropped, so the operator is wrong: FFT true residual 0.10 against the operator that includes the y wall term, MLMG/HYPRE 1e-12. FDS ignores y in `TWO_D`, so this is right for y of a TWO_D case and wrong for a one-cell x or z; the driver refuses the latter | **confirmed**; relies on the driver refusal |
| Assembled HYPRE on n×1×n | no y couplings are generated for Neumann y; PCG+BoomerAMG converges in 18 to 21 iterations (R=1, R=8) with pin-row residual ≤ 7e-13 and the same solution as MLMG. Periodic one-cell y (self-neighbour entries) was not run | Neumann **confirmed (harness)**; periodic one-cell y **unverified**; no production HYPRE backend exists to review yet |
| Mean removal in 2-D | the removal is the same code as in 3-D (components by flood fill over `beta ≠ 0` faces, pin at the lowest index, exact sums); the one-cell direction adds no face. Backend `comp_ns2d` reproduced here: `status=Ok`, 9 iterations, true residual 1.5e-13, `removed_rel` 7e-17, y gradient exactly 0 | **confirmed** (single level and the two-level `ns2d_16` pattern) |
| Zone sums with ratio 1 in y | per-level cell volume is `dx·dy·dz` of that level and uncovered counts follow the coarsened fine BoxArray (`ExactSum.H` takes a per-cell weight; backend test: 192 coarse + 256 fine cells, i.e. 4 fine per coarse in 2-D, not 8). The driver's own zone sums (`USUM`, `DSUM`, `PSUM`) with ratio (2,1,2) were not run | pressure-side sums **confirmed by backend test**; driver zone sums **unverified** |
| MLMG hidden-direction route (`setHiddenDirection(1)`, `ref_ratio_vect = 2 1 2`) | the feature exists (`AMReX_MLLinOp.H:89,951,1177-1179,1219`) but `m_amr_ref_ratio` is a single int per level (`rr[0]` when no hidden dimension, `:1223`), so an anisotropic ratio **requires** the hidden direction; and its kernels are implemented in `MLPoisson`/`MLALaplacian` only, not in `MLABecLaplacian` (one reference, ratio handling at `AMReX_MLABecLaplacian.H:718`), so the unit-spacing stretched/masked operator cannot use it. The backend therefore does not use it: it solves an extruded periodic copy (4 isotropic cells × product of ratios in the thin direction, ratio in that direction must be 1) after finding hidden-direction slow (73 to 200 iterations against 9) or divergent with more than one level (reported in `frozen/composite-notes.md`, not reproduced by me). Extruded `ns2d_16` composite equals a real 4-cell-y 3-D problem to 8.6e-17 (reproduced) | feature gap **confirmed (read)**; slow/divergent behaviour **reported, not reproduced**; extruded route **confirmed (run)** at 4× cells in y |
| Single-level hidden direction, periodic x,z (my `hid_test`, 64×1×64, 2 ranks) | solution equals the plain solve to 1.2e-15 (relative 1e-13); 8 iterations against 3 | **measured** (periodic only) |
| Single-level hidden direction with Neumann x,z | in my test it returned a zero solution after one iteration (residual reported 0), but the same test's *plain* Neumann solve also disagrees with the periodic-y solve by 75× while the backend harness's plain Neumann MLMG matches FFT to 1e-12, so my test setup is suspect | **unverified; do not rely on it** |
| Rest of D-057: BC-type mapping | the driver mapping table and its 21-check unit test are in the driver note; one-cell x/z with a Dirichlet face is refused; mixed open/closed face sets are refused by the selector (31 verification inputs); 23 code-0 inputs stay refused (Architect ruling) | read only (driver work); not re-run here |

Pressure-side conclusion: the single-cell-y handling works for FFT (Neumann/periodic y), MLMG single level, assembled HYPRE (Neumann y) and mean removal; the multilevel case works through extrusion, **not** through the hidden direction; the zone sums are right at the backend level and unverified at the driver level.

## 10. Failures, caveats and open items

- First stairwell `COMPCHECK` attempt failed (rc=6) because `cutk` was given comma-separated; AMReX arrays are space-separated. Re-run, no effect on results.
- Stretched 64³ runs use 2 timed solves at 4 ranks: accuracy, not timing. The 1M stretched rows have 10 (CPU) and 20 (GPU) timed solves, single repeats; the stairwell rows (§8) have 3 repeats.
- The cause of the Mf slowdown at R=64 is established in doc 07 §8 (AMReX hands HYPRE the row-scaled, nonsymmetric matrix; PCG then diverges at the first bottom solve); remedies and the recommended setting are there.
- The synthetic RHS means iteration counts on real FDS flow fields are still not measured (as in doc 05).
- Hidden-direction + Neumann behaviour not resolved (section 9). Periodic one-cell y in the assembled matrix not run. Driver zone sums in a refined 2-D case not run. Patch 0005 not built against master, no debug build.
- Decision for the Architect: gauge for FFT-solved cases (section 5, item 6).

Reproduce: `scratch/pressure-signoff/meanremoval/` (dense check, `patchwork/` binaries, cases, `ref_check.py`), `scratch/pressure-signoff/pressure_1M/` (`src/h2h_str.cpp`, `src/h2h.cpp` with `cutk`/`COMPCHECK`, `src/hid_test.cpp`, `pr06_*.sh`, `results_pr06/`).
