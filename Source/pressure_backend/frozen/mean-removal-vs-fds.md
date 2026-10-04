# Mean removal and gauge: FDS ULMAT/UGLMAT/GLMAT versus the pressure backend

Plain-language note. FDS line numbers are `Source/pres.f90` of this worktree (identical to the FireX reference
commit used for the runs); our line numbers are in `Source/pressure_backend`. Numbers come from
`stretched_study_results.txt` (all cases) and from the CTests named below.

## 1. What FDS does (ULMAT/UGLMAT and GLMAT)

| item | ULMAT / UGLMAT (`ULMAT_SOLVE_ZONE`, 1463-2001) | GLMAT whole domain (`GLMAT_SOLVER`, 3299-3658) |
|---|---|---|
| RHS | `F_H = PRHS * (DY*DX*DZ)`, times `RC` when cylindrical (1502); Neumann data added as `+-BXS/BXF*AF`, Dirichlet data as `-2*IDX*AF*BCV` (1603-1679). Rows are volume-scaled; matrix entries are `AF/DXN` (area over centre distance, 2596-2601). | same scaling, built by `GET_FH_FROM_PRHS_AND_BCS` (`ccib.f90` 2096) |
| singular? | `MTYPE = SYMM_INDEFINITE` unless some wall cell is OPEN or INTERPOLATED (2723-2785). A periodic face on a single mesh is INTERPOLATED, so it is a Dirichlet face for ULMAT: that zone is positive definite and has no zero mode. | `SYMM_INDEFINITE` for closed or periodic zones (matrix built globally, periodic faces are real couplings) |
| RHS mean | only for `SYMM_INDEFINITE` (1683-1733): **arithmetic, count-weighted mean of the volume-scaled `F_H`**, over gas Cartesian cells of the zone (`MUNKH>0`, same pressure zone) plus cut cells (1689-1708). Solid cells are not unknowns and are not counted. No volume weight. | zone mode (UGLMAT): same arithmetic mean (3371-3400). Whole-domain mode (GLMAT default, `func.f90` 7165): `SUM(F_H)/NUNKH_TOTAL`, **but `SUM_FH(2)` is only set inside `IF (N_MPI_PROCESSES>1)` (3402-3404), so on one MPI rank the mean is 0 and nothing is removed** (see section 6). |
| compatibility | `K x = F` with `K` symmetric, row sums zero. Solvable iff the **plain sum of the scaled `F_H` (boundary terms included) is zero**. | same |
| pin (HYPRE) | `F_H(N) = 0` (1748), last row identity and last column zeroed (2965-2983): the reduced system for `x_1..x_{N-1}` with `x_N = 0`. PARDISO uses its indefinite solver with no pin. | `F_H(last) = 0` (3447), same matrix treatment |
| solution mean | arithmetic mean of `X_H` subtracted (1779-1828). It is redundant: the gauge below overrides any constant. | none in the main path; arithmetic mean only for periodic test 7 or multi-zone whole domain (3549-3558) |
| final gauge | `SHIFT = sum(VOL*RHOP*(KRES + X_H)) / (sum(VOL*RHOP) + 20 eps)`, `X_H -= SHIFT` (1830-1878); `RHOP = RHO` (predictor) or `RHOS` (corrector); `VOL` includes `RC`. Result: `H = -X_H`, so **`sum(rho*V*(KRES - H)) = 0`**. | same formula (3495-3548) when `.NOT.PRES_ON_WHOLE_DOMAIN .OR. N_ZONE<=1` and `PERIODIC_TEST/=7` |
| sign | `K x = F` with `K = -V L`, `H = -X`, so `H` solves `L H = PRHS`, the same sign as our `phi` | same |

## 2. Our side, uniform cells

* Our right-hand side is the **unscaled** `PRHS` (`rhs`, valid cells). The operator is the unscaled 7-point `L` on uniform
  cells (`apply_operator`, `CommonLayer.cpp`). For uniform cells the scaled system `K x = V b` is the same system times `V`.
* `remove_mean` (`CommonLayer.cpp` 112-186): on uniform cells subtracts the exact volume-weighted mean of `b`,
  which equals `mean(V b)/V`, so `sum(V b') = 0`: the same system as ULMAT's arithmetic mean removal of `F_H`.
  It runs on the RHS before the backend (`PressureIface.cpp` ~119), is idempotent, per singular component, over
  uncovered cells.
* `apply_gauge` (`CommonLayer.cpp` 188-218): after the solve subtracts `sum(W (phi - g)) / sum(W)` per singular
  component with `W = V * gauge_weight`, `g = gauge_offset`. With no optional input it is the exact volume-weighted
  mean of `phi` (zero volume-weighted mean per component). So **yes, the common layer shifts the solution after the
  solve**; that single shift plays the roles of FDS's X-mean removal and gauge. Non-singular components are untouched.
* Tests (independent numpy: dense `K`, identity-row pin, arithmetic mean of `F` and of `x`): `pb_ulmat_equiv_neumann`,
  `pb_ulmat_equiv_periodic` on 8x6x5 cells with a compatible RHS and with an RHS offset (mean(F)/rms(F) = 0.98).
  FFT agrees to 1e-15 and MLMG to 2e-12 (relative L2) with the ULMAT-like result; our removed-mean diagnostic equals
  numpy's `|mean F|/rms F` to 1e-9.

## 3. Stretched z (non-uniform cell volumes)

Can the backends solve it? **No.** `FFT::Poisson` assumes uniform spacing, and `MLPoisson` takes `dx` from the
`Geometry` (no face coefficients). **Decision (a): the selector returns `NotBuilt`** for non-uniform
`PressureProblem::cell_width` (`PressureIface.cpp` 71-78, test `pb_selector`); equal widths must match the geometry
(else `InvalidInput`). The volume-scaled form (b) is the natural form of the planned assembled-matrix HYPRE backend
(and of `MLABecLaplacian` on a unit grid with face coefficients `AF/DXN`, `a = 0`, RHS `V b`).

For that path the common layer already takes per-cell volumes (`MeanKind`, `CommonLayer.H`):

* `ScaledArithmetic` (default when a volume field is given): `b_k -= mean(F)/V_k`, `F = V b`. This is FDS.
* `Volume`: `b -= sum(V b)/sum(V)`.

Both make `sum(F) = 0`; they differ by a multiple of `(V_k - mean V)`. Test `pb_meankind_stretched` (12x10x9, volume
ratio 4.7): both kinds equal the numpy formulas to 1e-14, the sum of the scaled RHS is zero, both are idempotent
(bitwise), results are bitwise identical across rank counts, and the two kinds differ by 1.4 relative for the offset RHS used (so the test
is sensitive).

### Study against FDS (12x10x16, z stretched with cell-volume ratio 3, 60 C upper half so rho varies, floor burner)

Numpy rebuilds the FDS matrix from the dumped widths, solves with the identity-row pin, removes the zero mode of `F`
either way, applies the FDS gauge, and compares with the FDS `H` of the same solve (scratch dump hook, below). The
RHS includes FDS's own boundary terms. "offset" is a scratch constant added to `PRHS` before the mean removal, to make
the RHS incompatible (as FDS computes it the mean is round-off, ~1e-17 of rms, and every weighting agrees).

All values are relative L2 versus FDS (H after FDS's gauge; "grad" is the difference of face gradients, i.e. the
velocity correction, relative to FDS's). Residual = `||K x - F1|| / ||F1||` of the scaled system.

| case | zero mode removed from F | mean(F)/rms(F) removed | residual | H mean-removed | H raw | grad |
|---|---|---|---|---|---|---|
| A closed, ULMAT, offset 0 | arithmetic | 2e-17 | 2.5e-14 | 3.6e-14 | 3.6e-14 | 1.9e-13 |
| A closed, ULMAT, offset 0 | volume | 3e-17 (c/rms b) | 2.5e-14 | 3.6e-14 | 3.6e-14 | 1.9e-13 |
| B closed, ULMAT, offset 0.3 | **arithmetic** | 7.3e-2 | 2.4e-14 | **4.2e-14** | **4.2e-14** | **2.0e-13** |
| B closed, ULMAT, offset 0.3 | volume | 1.09e-1 (c/rms b) | 2.5e-14 | 1.25e-1 | 1.30e-1 | 9.2e-2 |
| I same as B, second solve | arithmetic / volume | 7.2e-2 | 2.5e-14 | 5.1e-14 / 1.2e-1 | 5.2e-14 / 1.2e-1 | 2.1e-13 / 9.1e-2 |
| E closed, GLMAT with the single-rank fix, offset 0.3 | **arithmetic** | 7.3e-2 | 2.4e-14 | **4.2e-14** | 4.1e-14 | 2.0e-13 |
| E closed, GLMAT with the fix, offset 0.3 | volume | 1.09e-1 | 2.5e-14 | 1.25e-1 | 1.30e-1 | 9.2e-2 |
| H periodic x,y + solid z, GLMAT with the fix, offset 0.3 | **arithmetic** | 7.3e-2 | 3.4e-14 | **3.9e-14** | 4.0e-14 | 2.3e-13 |
| H periodic x,y + solid z, GLMAT with the fix, offset 0.3 | volume | 1.09e-1 | 3.4e-14 | 1.24e-1 | 1.29e-1 | 9.2e-2 |
| D closed, GLMAT as shipped (1 rank), offset 0.3 | none (FDS removes nothing) | 0 | 3.2 (row N dropped) | 3.2e-14 | 7.2e-14 | 1.5e-13 |
| G periodic x,y, GLMAT as shipped, offset 0.3 | none | 0 | 3.2 | 1.1e-14 | 2.4e-14 | 5.7e-14 |

Further observations:

* FDS's dumped `F_H` after its mean removal equals numpy's arithmetic result to 4e-16; the volume-weighted result is
  3e-2 away.
* ULMAT's X after its own arithmetic mean removal matches numpy's arithmetic removal to 4e-14 and the volume-weighted
  one only to 0.25 (before the gauge). After the rho-volume gauge the X-mean choice is irrelevant.
* Mean of the pinned solution relative to its rms is 0.93 to 0.99 in every case (the pin value is arbitrary); this is what the
  gauge removes.
* ULMAT on the periodic input is not comparable: the periodic faces are Dirichlet-coupled interpolated faces
  (`MTYPE = 2`, no zero mode), so the periodic study uses GLMAT, which assembles the true periodic coupling.
* To get ULMAT to build its matrix on a regular single mesh at all (otherwise FFT is used, 1221-1267), the scratch
  copy has `FDS_DBG_FORCE_ULMAT`.

**Which matches FDS, and why.** The arithmetic (count) mean of the volume-scaled RHS, because the compatibility
condition of the system FDS solves is the plain sum of the scaled `F_H`; FDS subtracts the same number from every
row. Removing a volume-weighted mean of `b` instead subtracts `c V_k` from `F_k`, which also gives `sum F = 0` but
adds a different, non-constant source `(c V_k - mean(F))`, so the solution differs by about 12% of its norm in our
case. When the RHS is already compatible (round-off mean) they are indistinguishable, which is why this never shows
in ordinary runs. After the solve the weight does not matter as long as the final rho*volume gauge follows, because
it fixes the free constant.

## 4. Cylindrical and masked cells

* Cylindrical (`&MESH CYLINDRICAL=T`, rows and RHS scaled by `RC`, `VOL` with `RC`): not built.
  `PressureProblem::cylindrical` (or a non-Cartesian `Geometry`) returns `NotBuilt` with an explicit message
  (`PressureIface.cpp` 69-71; selector test).
* Masked/obstructed cells: FDS ULMAT drops solid cells from the unknowns and from the mean (1689, 1716); GLMAT in zone
  mode does the same, whole-domain GLMAT keeps them as unknowns. The backend returns `NotBuilt` for `cell_class != 0`
  ("obstructed/solid, known-value or pinned cells ... masked branch"), covered cells, driver component ids, composite
  levels and variable coefficients. The mean must use the same cell set as the solved system; the common layer already
  takes `uncovered` and the component label for that (D4).

## 5. Final gauge with rho and KRES (implemented)

`PressureProblem::gauge_weight` (rho) and `gauge_offset` (KRES): after the solve `phi -= sum(V rho (phi - KRES)) /
sum(V rho)` per singular component, computed with the exact sum (decomposition independent). FDS passes `RHO`
(predictor) or `RHOS` (corrector) and `KRES`. Null means weight 1 and offset 0, the plain volume mean. FDS adds
`20 eps` to the denominator; we do not. Tests `pb_ulmat_gauge_neumann` and `pb_ulmat_gauge_periodic` (variable rho
and KRES, FFT and MLMG against the independent numpy FDS formula, rel L2 <= 1e-9, achieved 1e-16 (FFT) and 2e-13 (MLMG); also
`sum(rho (H - KRES)) = 0`). The weighted gauge constant differs from the plain gauge by 3.8e-2 in the test (sensitive).

## 6. Finding in FDS (not changed here)

`pres.f90` 3402-3404: in whole-domain GLMAT (`PRES_ON_WHOLE_DOMAIN`, the GLMAT default) `SUM_FH(2)` is set only by
`MPI_ALLREDUCE` when `N_MPI_PROCESSES>1`; on one rank `SUM_FH(2)` stays 0 (set at 3370), so `MEAN_FH = 0` and the RHS
mean is **not** removed. The fallback X-mean at 3553-3555 has the same pattern. Observed: GLMAT on one rank equals the
"no removal" variant to 1e-14 (cases D, G) and the "arithmetic" variant only once the sum is fixed (cases E, H).
With a compatible RHS (round-off mean) the effect is invisible; with an incompatible RHS the HYPRE pin makes the
solution depend on which row is dropped. The fix is one line (`SUM_FH(2) = SUM_FH(1)` for one rank). To be flagged to the
Pressure Solver Lead and the Chief Architect; the reference tree is untouched.

## 7. Reproduce

Backend tests: `ctest -R "ulmat|meankind|selector"` in the build directory (see `README.md`).

FDS study (scratch copy only; the hook edits a copy of `pres.f90`):
```
tar xzf <reference tarball> -C <scratch> ; python3 frozen/fds_dump_hook.py <scratch>/Source/pres.f90
# configure as the reference build (HYPRE ON) and build the copy, then, in a run directory:
FDS_DBG_FORCE_ULMAT=1 FDS_DBG_RHS_OFFSET=0.3 <scratch fds> in.fds      # in.fds from frozen/fds_cases/ (set &PRES SOLVER and T_END)
python3 frozen/stretched_study.py dbg_ULMAT_1.txt                        # periodic: --periodic-xy dbg_GLMAT_1.txt
FDS_DBG_FIXMEAN=1 ...                                                    # single-rank GLMAT with the mean fix
```
`dbg_*_1.txt` is the first predictor solve, `_2` its corrector, `_3`, `_4` the next step.
