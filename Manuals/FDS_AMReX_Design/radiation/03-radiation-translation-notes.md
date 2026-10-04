# Radiation loops through the generator: translation notes

Owner: AMR Radiation Lead. Worktree: `(local project directory)/src-s5gen` (branch `s5-gen`), directory `amrex/s4_mass/s5_gen/`. Reference tree: `Source/radi.f90` at `bee11f0329` (upstream `afb5e31a48` merged). Line numbers below are that file.

Companion docs: `02-radiation-gpu-candidates.md` (which loops, with modelled shares), `00-r1-signoff.md` (wall loops), `04-radiation-phase4-notes.md`.

## 1. Result in one paragraph

17 kernels were generated from `radi.f90`. They cover the per-band source `RTE_SOURCE`, the extinction and `UIIOLD` copies, the `UIID` accumulation, `UII`, `QR`, `QR_W`, the emission and absorption fields, three `KFST4_GAS` branches (optically thin, corrected, gray) and the wall loop L1239. All are element-wise cell or wall kernels. Each is compared bit for bit with the verbatim Fortran, compiled with gfortran, in six flag sets and at 1, 4 and 8 threads. The angular sweep (L1242 proper), the open-boundary wall loop L1243, the interpolation wall loop L1248, the wavelength-band and WSGG source branches, and the droplet branch were **not** translated; the generator refuses each, for the reason in §5. L1243 has a hand-written specification kernel and a test with a mutant, and is waiting for a wall-table feature.

The translated kernels are element-wise pieces of L1242's body and of its neighbours; the sweeps, which are the bulk of the body, are not among them. The modelled share of L1242 is 4.576% of the run (rank 7 of 829); all other radiation loops are 0.024% or less. See `02-…`.

## 2. What was translated

Each kernel is an entry in `markers/rad_kernels.toml`. Lines are `radi.f90` at `bee11f0329`.

| Kernel | Lines | Statement | Form |
|---|---|---|---|
| `rad_extcoe` | 4206 | `EXTCOE = KAPPA_GAS + KAPPA_PART + SCAEFF + SCAEFF_G` (whole array) | element-wise over `0:IBAR+1` etc. |
| `rad_uiiold_wb` | 4213 | `UIIOLD = UIID(:,:,:,IBND)` | element-wise copy |
| `rad_uiiold_gray` | 4215 | `UIIOLD = UII` | element-wise copy |
| `rad_rte_source` | 4221-4222 | `RTE_SOURCE = KFST4_GAS + KFST4_PART + RSA_RAT*(SCAEFF+SCAEFF_G)*UIIOLD` | element-wise; the new precompute |
| `rad_uiid_wb` | 4818 | `UIID(:,:,:,IBND) = UIID(:,:,:,IBND) + WEIGH_CYL*RSA(N)*IL` | element-wise |
| `rad_uiid_gray` | 4820 | `UIID(:,:,:,ANGLE_INC_COUNTER) = ...` | element-wise |
| `rad_qr_wb` | 4949 | `QR = QR + KAPPA_GAS*UIID(:,:,:,IBND) - KFST4_GAS` (left to right: `(QR + K*U) - F`) | element-wise |
| `rad_qrw_wb` | 4954 | `QR_W = QR_W + KAPPA_PART*UIID(:,:,:,IBND) - KFST4_PART`, particles | element-wise, host guard |
| `rad_qr_gray` | 4982 | `QR = KAPPA_GAS*UII - KFST4_GAS`, gray | element-wise |
| `rad_qrw_gray` | 4983 | `QR_W = QR_W + KAPPA_PART*UII - KFST4_PART`, gray | element-wise, host guard |
| `rad_emis_wb` | 4951 | `RADIATION_EMISSION = RADIATION_EMISSION + KFST4_GAS` | element-wise |
| `rad_abs_wb` | 4952 | `RADIATION_ABSORPTION = RADIATION_ABSORPTION + KAPPA_GAS*UIID(:,:,:,IBND)` | element-wise |
| `rad_uii_sum` | 4963 | `UII = SUM(UIID,DIM=4)` | ordered inner sum over `1..UIIDIM` |
| `rad_thin_kfst4` | 4131-4141 | optically thin `KFST4_GAS` | cell loop with solid test |
| `rad_corr_kfst4` | 4151-4161 | `RTE_SOURCE_CORRECTION` `KFST4_GAS` | cell loop with clip |
| `rad_gray_kfst4` | 4099-4117 | gray `KFST4_GAS` | cell loop, reductions on host |
| `rad_wall_qin_zero` | 3887-3893 (L1239) | reset wall `Q_RAD_IN` | wall loop, needs `B1_PRESENT` |

Notes.

- **Operation order is reproduced exactly.** The text rewrite keeps the original expression and only adds the indices, so the association of `A+B+C*D*E` is the original. `-ffp-contract=off` is used in the tests so that FMA does not hide a difference. On the device, FMA contraction must be turned off for these kernels if the bitwise claim is to hold against gfortran.
- **`SUM(UIID,DIM=4)`**: gfortran evaluates it as a sequential chain over `n=1..UIIDIM` starting from +0. The kernel is the same chain; the mutants for a reversed order and for a start at -0 are both detected. ifx may vectorise `SUM`, so the claim is against gfortran only.
- **`rad_rte_source`** is the new precompute. The old inline expression had the same grouping (`02-…` Part B), so the in-loop value and the precompute are bit-identical; the test compares the precompute with the verbatim statement.
- **UIID accumulation on the device.** It runs once per angle; fuse it into the sweep kernel rather than launching `rad_uiid_*` separately per angle.

## 3. What each test covers

Test program `test/s5_rad_bitwise.F90`; cases built by `test/make_rad_tests.py`; reference module is verbatim `radi.f90` text.

- **Sizes.** From 1x1x1 to 24x20x16, including degenerate directions of size 1 and non-cubic boxes.
- **Layers.** All kernels run over the full ghost-inclusive range `0:IBAR+1` etc. wherever the original does. Ghost values are filled with data (not zero) and compared.
- **Data classes** for every real input: random, wide magnitude (many decades), exact zero, +0 and -0, integer-valued ties, denormal-range values, constants, and cancellation cases (terms that nearly cancel, so an association change shows).
- **Solid cells.** `CELL_INDEX`/`CELL` with solid and mixed cells for the three `KFST4_GAS` kernels.
- **Tied and odd cases.** `QR_CLIP` equal to `CHI_R*Q` (clip tie), `WEIGH_CYL` equal to 1 and 2 and random, `NS` in 1, 2, 5, 6, 12 (one, few and many angle-set bins), `N` over the angle range, `RSA_RAT` positive and negative.
- **Guarded kernels.** The cases for `rad_qrw_*` are skipped (and counted as guard-skipped, 148 in the run) when the original `IF` would not run them; the test does not claim anything about them in that case.
- **Wall loop L1239.** The generated `rad_wall_qin_zero` against the original loop with: walls with and without `B1`, `NULL_BOUNDARY` walls, and `TMP_GAS_FRONT` of a positive value, zero, -0 and -1e-300.
- **Threads and flags.** Serial and OpenMP at 1, 4 and 8 threads; flag sets O0, O2, O0omp, O2omp, O2omp_off (`-DS4_OFFLOAD`), O2omp_dpd (`-DS4_OFFLOAD -DS5_FORCE_DPD`). All with `-ffp-contract=off`. Exact comparison of the bit patterns of every element of every output array (`same3`, `same4`, `same1`), ghost layers included, so an extra or missing write shows.
- **Pass status.** See §7 for the final run.

Python-side tests (`test/test_mesh_rad.py`):

- generation: every kernel is produced and equals the golden text (`test/rad_kernels.golden`);
- write sets: each kernel writes exactly the arrays the original statement writes;
- rewrite checks: the rewrite keeps the number of lines; it touches only the lines of the entries and the one `ALLOCATE(MOLD=)` line; an entry that overlaps a rewritten range is refused; an edited anchor is refused (drift); a `[[file_rewrite]]` whose text is not on its line, or on the wrong line, is refused;
- nine loops that the generator refuses today are pinned with their refusal message (§5), so a generator change that starts accepting them shows up as a test failure to review, not a silent change;
- `--mutants`: 17 mutants of the generated kernels (regrouped sums such as `A + (B + C)`, a neighbour index, a sign, a reversed or -0-started `SUM`, an inverted solid test, `>=` for `>` at the clip tie, overwritten instead of accumulated, three wall-condition mutants) and one mutant of the specification kernel. The bitwise test must fail on every one.

## 4. Preconditions of the kernels

A kernel is bitwise the original only under these conditions. The driver (or the sidecar) must hold them.

1. `CC_IBM=.FALSE.` for `rad_thin_kfst4`, `rad_corr_kfst4`, `rad_gray_kfst4`. The `CCVAR`/cut-cell branches are blanked in the sidecar. Cut cells are out of scope until the cut-cell design.
2. `rad_gray_kfst4`: the reductions `RAD_Q_SUM` and `KFST4_SUM` are excluded from the kernel (the generator cannot do an ordered scalar reduction in a kernel today). The host keeps them as ordered sums (D-053).
3. `rad_qrw_wb`, `rad_qrw_gray`: launched only if the original `IF (N_LP_ARRAY_INDICES>0 ...)` holds. The `IF` is stripped from the kernel.
4. `RSA(N)` is passed as a scalar; `UIIDIM` equals `NS` (`rad_uii_sum`).
5. `rad_wall_qin_zero`: the wall table carries `B1_PRESENT`.
6. `RTE_SOURCE` is a driver-owned array of the same bounds as `KFST4_GAS`. The kernel does not allocate it. The original `ALLOCATE(RTE_SOURCE, MOLD=KFST4_GAS)` is rewritten to an explicit-bounds `ALLOCATE` only so that the file parses; it is host code.

## 5. What the generator refuses, and why it matters

These are pinned in `test_mesh_rad.py` with the exact refusal text.

| Loop | Lines | Refusal | Needed |
|---|---|---|---|
| 2D sweep | 4472-4492 | loop with a step | variable-step loop; upwind recurrence `IL(I-ISTEP,...)`; FR-062 design for boxes |
| cylindrical sweep | 4433-4468 | loop with a step | same |
| 3D slice sweep | 4496-~4830 | not attempted (cell-list wavefront `IJK_SLICE`, `N_SLICE`, `M_IJK`) | wavefront scheduling plus the per-box design |
| `UIID` IF construct | 4817-4821 | neither a K,J,I nest nor a whole-array fill: `If_Construct` | the two arms are generated separately (`rad_uiid_wb`, `rad_uiid_gray`); the host selects |
| L1243 open boundary | 4965-4974 | pointer assignment `BR => BOUNDARY_RADIA(WC%BR_INDEX)` is not a recognised wall alias | derived-type record per wall; ragged `BR_ILW(NRA,NSB,wall)` |
| L1248 `INTERPOLATE_IL` wall loop | 3651-3662 | same alias refusal | same, plus per-thread scratch for `ILW_OLD` |
| wide-band `KFST4_GAS` | 4027-4046 | reference to `BLACKBODY_FRACTION`: a source function or unknown name | callee support |
| RADCAL `KAPPA` | 4005-4016 | reference to `GET_KAPPA`: unknown name | callee support |
| WSGG | 4051-4082 | callee `GET_MASS_FRACTION` uses module variable `Z2Y` | module-variable support for callees; `A_WSGG` is recursive |
| droplet | 3959-3972 | callee `INTERPOLATE1D`: dummy `X` has an assumed-shape form | assumed-shape rank-1 dummy |
| rank-1 `BR_ILW` | n/a | "rank-1 array BR_ILW in a wall loop is not declared in `[policy.wall.arrays]`"; a rank-3 `BR_ILW` is silently mapped to a mesh-field shape | a ragged per-wall table in the wall policy |
| reduction scalars | various | `IBND` etc: "scalar written in the loop and may be read afterwards" | worked around with `private = [...]` in the entry |

**L1243 (open-boundary `Q_RAD_IN`).** The test has a hand-written specification kernel `spec_open_qin` that implements the agreed form (per-band temporary starting at +0, `T=T+ILW(N)` over angles, then `Q = Q + T`) and compares it with the original loop in a multi-band case. A mutant (`-DRAD_MUT_FLATSUM`, the flat chain `Q = Q + ILW(N)` over all bands) fails on wide-magnitude data. This is the evidence for the "accept with change" sign-off in `00-r1-signoff.md`. It stays a specification until the generator can express the wall table.

## 6. Missing generator features

Concrete descriptions, in order of use to radiation.

1. **F2008 parse.** `ALLOCATE(RTE_SOURCE, MOLD=KFST4_GAS)` (`radi.f90:4220`) is Fortran 2008, and the shared front end uses `std="f2003"`, so `radi.f90` does not parse. Switching to `std="f2008"` is not drop-in: fparser2 then returns `Fortran2008.If_Stmt` and `Block_Nonlabel_Do_Construct`, so the `type(n) is F.If_Stmt` tests in `s5gen.py` stop matching. Needed: accept the F2008 subclasses (use `isinstance`), or parse `MOLD=` as a statement. Workaround: a `[[file_rewrite]]` text replacement.
2. **Whole-array and section statement expansion with ghost bounds.** `X = A + B*C` where all are `(0:IBP1,0:JBP1,0:KBP1)`, and `UIID(:,:,:,IBND) = ...`. Needed: expand to a K,J,I nest over the declared bounds of the arrays (not `1:IBAR`). Workaround: in-place text rewrite (`rad = "expand"`).
3. **`SUM(A, DIM=n)`** as an ordered inner loop. Workaround: `rad = "sum4"`.
4. **Local allocatable used as driver scratch.** `RTE_SOURCE` is `ALLOCATABLE` in `INTENSITY_UPDATE`; the kernel needs it as an argument of fixed shape. Today it is registered as an OUT dummy. Needed: a way to declare "local allocatable becomes a persistent per-box work array" (the scratch-pointer support doc covers per-thread scratch, not this).
5. **Derived-type component that exists only in the flat wall table** (`B1_PRESENT`): the generator should emit a presence flag for `ASSOCIATED(WC%B1)`-style tests. Today: `type_comps` in the sidecar.
6. **Ragged per-wall tables and record aliases.** `BR => BOUNDARY_RADIA(WC%BR_INDEX)` and `BR%BAND(N)%ILW(A)`. Needed: a wall array policy for a table `(NRA, NBANDS, wall)` addressed by `BR_INDEX`, and recognition of the alias. The wall policy today accepts only per-wall tables indexed by `IW`.
7. **`CC_IBM`/`CCVAR` branches.** Today blanked under the precondition `CC_IBM=.FALSE.`. Needed: a compile-time or launch-time switch that keeps both arms.
8. **Ordered scalar reductions** (`RAD_Q_SUM`, `KFST4_SUM`, D-053): deterministic chain or fixed-point sum inside a kernel.
9. **Callees that are source functions**: `GET_KAPPA`, `BLACKBODY_FRACTION`, `A_WSGG` (recursive), `KAPPA_WSGG`, `GET_VOLUME_FRACTION`, `GET_MASS_FRACTION` (uses module `Z2Y`). Each needs either inlining support or a `!$omp declare target` path with the module data passed in.
10. **Assumed-shape rank-1 dummy** (`INTERPOLATE1D`).
11. **Loop with a variable step, loop-carried upwind recurrence with `CYCLE`, and the cell-list wavefront** (`IJK_SLICE`, `N_SLICE`, `M_IJK`). This is the sweep; the generator cannot express the data dependence. It needs the FR-062 per-box design first (the dependence is then local to a box). Not a generator task to do blind.
12. **Pitfall in my own front-end hook (fixed).** The text rewrite is global to the file text; an entry that overlaps a rewritten range would silently see blanked lines. `configure()` now refuses it.

## 7. Test run and commits

Run on the worktree at `743dbce961` plus the uncommitted radiation files, with `(local workspace)/s5venv/bin/python` and gfortran 14.

| Check | Result |
|---|---|
| Fast checks (`test_mesh_rad.py`: generation, golden, write sets, rewrite/overlap/drift refusals, nine pinned refusals) | all PASS, no FAIL line |
| Bitwise program, six flag sets (O0, O2, O0omp, O2omp, O2omp_off, O2omp_dpd), 1/4/8 threads, 8632 cases each | **PASS in all six**: 0 mismatches, 0 kernels with too few non-vacuous cases, 148 guard-skipped cases; run by `test_mesh_rad.py --all` under the shared lock `.s5gen.lock` ("radiation bitwise test passes in all six flag sets") |
| Spec mutant (flat chain `Q += ILW(n)` for L1243) | **detected**: clean build 0 mismatches, `-DRAD_MUT_FLATSUM` build 51 mismatches (private output directory, O2omp) |
| 17 generated-kernel mutants | the sweep runs under the same lock after the bitwise stage; see the line below |

Test bug found and fixed on the way: the serial flag sets (O0, O2) failed to link because the generated test called `omp_set_num_threads` without the `!$` sentinel. Fixed in `make_rad_tests.py` (two lines). The five sets other than the first O2omp run were first run in private output directories while the lock was held by other jobs for over an hour; the final claim above is from the run under the lock.

Mutant sweep status and commit hashes: see the report; commits are local only.

## 8. Shared-file requests (not made by me; I own none of these)

For the GPU Generator Engineer:

- `new_builder` dispatch for `rad` entries (or merge `markers/rad_kernels.toml` into the shared sidecar). Until then use `test/rad_gen.py`.
- Fix the F2008 parse (feature 1).
- Extend `[policy.wall.arrays]` for a ragged `BR_ILW` (feature 6).
- Add the `B1_PRESENT` flag to the wall table.

For the Wall Loops Engineer: `B1_PRESENT` and the `BR_ILW` table `(NRA, NBANDS, wall)` addressed by `BR_INDEX` with the host-side assertion `0 < BR_INDEX <= N_BOUNDARY_RADIA_DIM`.

For the Legacy Mapper: the list in `02-…` Table A1 and A2 is the proposed split. Radiation covers the elementwise kernels above and the R1 wall loops; the Mapper keeps the callee-dependent branches and the sweep until the FR-062 design is implemented.

## 9. Messages not sent

This session had no tool to message the other roles. The drafts below are ready to forward.

**To the Legacy Mapper (claim per `loop-work-list.md` §4).** My uncommitted radiation sidecar (`markers/rad_kernels.toml`) and builder (`s5_rad.py`) cover these loop ids of `radi.f90`:

- **L1239** (RADIATION_FVM wall `Q_RAD_IN` reset): kernel `rad_wall_qin_zero`, bitwise test written; status in progress until committed. Needs `B1_PRESENT` in the wall table.
- **L1243** (open-boundary `Q_RAD_IN`) and **L1248** (INTERPOLATE_IL wall loop): claimed, blocked on the ragged `BR_ILW(NRA,NBANDS,wall)` table and the `BR=>BOUNDARY_RADIA(WC%BR_INDEX)` alias. L1243 has a specification kernel and a mutant test. L1248 is low priority (runs only at a new angle cycle).
- **L1245** (RADF file output): host side, not planned.
- **L1242** (RADIATION_FVM, 4.576%): not claimed as a whole loop. Sixteen element-wise sub-nests inside it (and next to it) are translated; they are not whole-loop kernels, so the register should not count them as L1242. Their ids are the kernel names in §2. The sweep inside L1242 waits for the FR-062 per-box design.

The survey numbers (36975d7) for L1239, L1243, L1245, L1248 are 3886-3892, 4961-4970, 5044-5059, 3651-3662; the same loops in the reference tree at `bee11f0329` are 3887-3893, 4965-4974, 5048-5063, 3651-3662. The loop-work-list row for L1243 says "neighbour mesh (T4)"; the actual blocker is the ragged per-wall `BR_ILW` table, not neighbour-mesh data. The row for L1239 says the `B1_INDEX` gather is the issue; it is, plus the `B1_PRESENT` flag.

**To the GPU Generator Engineer.** I own `s5_rad.py`, `README_rad.md`, `markers/rad_kernels.toml`, and `test/rad_gen.py`, `make_rad_tests.py`, `s5_rad_tutil.F90`, `s5_rad_bitwise.F90`, `run_rad.sh`, `test_mesh_rad.py`, `rad_kernels.golden`. I did not edit `s5gen.py`, `s5_markers.toml` or any other shared file. Requests: the four in §8.

**To the Wall Loops Engineer.** Radiation needs `B1_PRESENT` in the flat wall table and a ragged `BR_ILW(NRA,NBANDS,wall)` table. Details in `00-r1-signoff.md`. The L1243 specification test is ready (`spec_open_qin`) and is the acceptance test for the generated kernel.

**To the Chief Architect.** The R1 sign-off text is in `00-r1-signoff.md` (accept L1239; accept with change L1243 and L1248; L1245 host; L1242 not covered). It is for you and the Legacy Mapper to apply to `docs/amrex/blocked-loop-families.md`.
