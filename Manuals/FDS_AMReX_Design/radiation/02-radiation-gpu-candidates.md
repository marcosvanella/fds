# Radiation loops: GPU candidates and review of the RTE_SOURCE change

Owner: AMR Radiation Lead. Status: draft, not reviewed.

**Pins.** The inventory CSV (`docs/inventory/gpu_candidate_loops.csv`) and the Legacy Mapper's survey (`docs/inventory/gpu_callgraph_survey.md`) are on FireX `36975d765f`. The code and the line numbers quoted here ("merged file") are `Source/radi.f90` at the local merge of upstream `afb5e31a48` (s5-gen reference commit `bee11f0329`). Between the two, `radi.f90` grew by the `RTE_SOURCE` change and other upstream edits; the shift is +1 for lines 3798-4217, about +8 for 4218-4456, then falling to **+4 from 4943 to the end of `RADIATION_FVM`**. A loop id (`L12xx`) is the CSV id; the lines are merged-file lines.

## Part A. Candidate table

**Where the time is.** The survey gives the RADI timer as 4.6% / 5.0% / 3.6% of run time in the three pooled runs and models loop shares inside it. One loop carries almost all of it: **L1242, the `BAND_LOOP` of `RADIATION_FVM` (radi.f90:3913-4957), 4.576% of total time, rank 7 of the 829 device-eligible loops**, a tree loop (13 callees reachable, recursion in `A_WSGG`), blocked and geometry-deferred (O4: the `CCVAR`/`CC_SOLID` test). Everything else in the CSV is 0.024% or less. The CSV has no rows for the pieces inside L1242, so every row below marked "in L1242" shares that 4.576% and has no number of its own. The radiation stage has 29 loops in the CSV: 13 candidates (10 leaf, 3 with callees, 3 particle loops), all 13 with a body blocker, 16 non-candidates (survey section 3).

**What is plain and what is not.** The per-cell source terms, the `RTE_SOURCE` precompute, the `UIID` accumulation, the `UII` sum and the `QR` updates are elementwise kernels with no dependence between cells. The angular sweep itself is not: `IL(I,J,K)` is computed from the upwind cells `IL(I-ISTEP,J,K)`, `IL(I,J-JSTEP,K)`, `IL(I,J,K-KSTEP)` (radi.f90:4435-4437, 4474-4475, and the 3D forms from 4494), so cells must be visited in upwind order. The 3D sweep does this with hyperplanes: a list `IJK_SLICE` of the cells of each diagonal plane (`N_SLICE`, `M_IJK`), one `!$OMP PARALLEL DO` per plane, `CYCLE` on solid cells, and `CELL_ILW` overrides from wall-adjacent faces. That part needs the FR-062 per-box lagged-sweep design (one kernel per box and angle with the wavefront inside it, box-face intensities double-buffered) and is **not** a generator task at its present capability.

**Verdict key.** *Translated* = done with the generator and bitwise-tested (see `03-radiation-translation-notes.md`). *Generator feature* = translatable once a named feature exists. *FR-062* = needs the per-box lagged sweep design. *Host* = stays on the host.

### A1. Inside L1242 (`RADIATION_FVM` band loop, radi.f90:3913-4957, 4.576% as a whole)

| Piece | Lines (merged) | Blockers | Verdict | Feature needed |
|---|---|---|---|---|
| Band initialisation, whole-array zeroing of `KFST4_GAS`, `KFST4_PART`, `SCAEFF`, `SCAEFF_G`, `KAPPA_PART` | 3914-3920 | whole-array statements | translatable now (fills) | none; not done (trivial) |
| Liquid-droplet absorption and scattering per particle class | 3959-3972 | `CALL INTERPOLATE1D(LPC%R50, LPC%WQABS(:,IBND), ...)` has an assumed-shape dummy; per-class host loop | generator feature | callee with an assumed-shape rank-1 table argument (generator refuses: "assumed-shape form not rank-one `LOWER:`") |
| Solid-particle deposit scatter (`GET_IJK`, add into cell arrays) | 3981-3991 | scatter into cells from particles: a race on the device | FR-062 / FR-005 (4a) design | deterministic deposition (sorted segmented sum or fixed-order per-box accumulation); no atomics today, none allowed |
| RADCAL `KAPPA_GAS` loop | 4002-4019 | `GET_KAPPA` is a source function (loop over `N_RADCAL_ARRAY_SIZE`, `DOT_PRODUCT`, `GET_MOLECULAR_WEIGHT`, `GET_MASS_FRACTION_ALL`); private `ZZ_GET` allocated per thread; `CC_IBM` branch | generator feature | source-function callees with module tables; per-thread automatic arrays (the scratch-pointer mechanism of `amrex/scratch-pointer-support.md` is for fixed-size private arrays and may carry `ZZ_GET`) |
| Wide-band `KFST4_GAS` and `RAD_Q_SUM`/`KFST4_SUM` | 4027-4046 | `BLACKBODY_FRACTION` (source function, `INTERPOLATE1D_UNIFORM` on the module table `BBFRAC`); module-scalar sums; `CC_IBM` branch | generator feature | source-function callees; module table `BBFRAC`; ordered scalar reduction (D-053) |
| WSGG `KAPPA_GAS`, `KFST4_GAS`, sums | 4051-4082 | `A_WSGG` is `RECURSIVE`; `GET_VOLUME_FRACTION`, `GET_MASS_FRACTION` (module array `Z2Y`); `CALL`s in the body; module-scalar sums | generator feature | recursion removal (the CSV's own blocker `recursion`), module arrays, ordered reduction |
| Gray gas `KFST4_GAS` fill with RTE-correction sums | 4099-4117 | `RAD_Q_SUM_PARTIAL` per thread then `!$OMP CRITICAL` (thread-count dependent, FR-005 (ii)); `CC_IBM` branch | **translated** for the fill (kernel `rad_gray_kfst4`, reductions removed, CC_IBM off); sums need the feature | ordered scalar reduction (D-053) for `RAD_Q_SUM`/`KFST4_SUM` |
| Optically thin `KFST4_GAS = CHI_R*Q + KAPPA_GAS*UII` | 4131-4141 | `CC_IBM` branch | **translated** (`rad_thin_kfst4`, CC_IBM off) | `CC_IBM`/`CCVAR` handling for the cut-cell case |
| RTE-source correction scaling of `KFST4_GAS` | 4149-4163 | `CC_IBM` branch | **translated** (`rad_corr_kfst4`, CC_IBM off) | as above |
| `ADD_VOLUMETRIC_HEAT_SOURCE` (CSV L1238, 0.024%) | call 4167, routine 5098-5146 | pointer-heavy tree loop over `INIT` regions (`alloc_ptr_local` 26, `reduction_accum` 4); only with `INIT_HRRPUV` | host | none; 0.024% does not justify it |
| Condensed-species absorption | 4171-4202 | not examined in detail | open | to be classified when L1242 is taken up |
| `EXTCOE = KAPPA_GAS + KAPPA_PART + (SCAEFF+SCAEFF_G)*RSA_RAT` | 4206 | whole-array statement | **translated** (`rad_extcoe`) | whole-array expansion (done by the rewrite in `s5_rad.py`) |
| `UIIOLD` copies | 4213, 4215 | whole-array / array section | **translated** (`rad_uiiold_wb`, `rad_uiiold_gray`) | as above |
| **`RTE_SOURCE` precompute** | 4220-4222 | `ALLOCATE(MOLD=)` (F2008) | **translated** (`rad_rte_source`) | F2008 parse (done by a one-line text rewrite); local allocatable as driver scratch |
| Wall `OUTRAD` loops and `INRAD_W` | 4229-4258, 4262-4280 | `BOUNDARY_RADIA` records, `SF` gathers | generator feature | same wall-table features as R1 (ragged `BR_ILW`, `BR` alias) |
| Wall ghost intensities (`WALL_LOOP1`, `CFACE_LOOP1`), `CELL_ILW` | 4328-4371, 4375-4384 | per-wall `ILW` records; write to `CELL_ILW(IC,ABS(IOR))` | generator feature + FR-062 | ragged wall table; `CELL_ILW` is the sweep's boundary input |
| **Angular sweep**, cylindrical | 4431-4468 | variable step `K=KSTART,KEND,KSTEP`; reads own output upwind; named `CYCLE`; `CELL_ILW`; `RFPI`, `DPHI0` | **FR-062** | per-box sweep kernel; loop-carried upwind recurrence with `CYCLE` |
| **Angular sweep**, 2D | 4470-4492 | same | **FR-062** | same |
| **Angular sweep**, 3D (STEP, DIAMOND, EXPONENTIAL, alternate schemes) | 4494-4718 | hyperplane wavefront (`IJK_SLICE`, `N_SLICE`, `M_IJK`), `!$OMP PARALLEL DO` per plane, `ILDX/Y/Z` | **FR-062** | same; plane lists built on the host or by a device scan |
| Cylindrical copy, `WALL_LOOP2/3`, `CFACE_LOOP2/3` | 4724-4743, 4749-4788, 4791-4812 | `BR` records, `ILW(ANGLE_INC_COUNTER)` updates | generator feature | ragged wall table |
| **`UIID` accumulation** (two forms) | 4817-4821 | array section | **translated** (`rad_uiid_wb`, `rad_uiid_gray`) | section expansion; `RSA(N)` passed as a scalar |
| `IL_S` interpolation to other meshes | 4826-4843 | gather into `OMESH%IL_S`; skips `NM==NOM` | replaced by the FR-062 face exchange | not a kernel in the AMR driver |
| Oriented-particle intensity | 4848-4881 | particle loop | host | none |
| `RADF` save | 4886-4895 | file output buffer | host | none |
| Wall incoming flux, cface incoming flux | 4902-4920, 4922-4940 | wall/cface records, `CC_IBM` | generator feature | wall-table features |
| **`QR` band update** and emission/absorption stores | 4949-4954 (wide band / WSGG), 4981-4988 (gray) | array sections, `IF` guard | **translated** (`rad_qr_wb`, `rad_qrw_wb`, `rad_emis_wb`, `rad_abs_wb`, `rad_qr_gray`, `rad_qrw_gray`) | section expansion |
| **`UII = SUM(UIID,DIM=4)`** | 4963 | intrinsic reduction over a dimension | **translated** (`rad_uii_sum`, ordered inner sum from +0) | `SUM(...,DIM=)` expansion |

### A2. Separate CSV loops

| CSV id | Routine, lines (merged) | Modelled share | Blockers (CSV) | Verdict | Feature needed |
|---|---|---|---|---|---|
| L1242 | `RADIATION_FVM` band loop, 3913-4957 | **4.576%** (rank 7) | alloc_ptr_local 462, dt_alloc_component 93, pointer_reassign 63, reduction_accum 23, module_global_write 6; callees: recursion | **not covered by the R1 sign-off**; geometry-deferred (O4); pieces above | FR-062 for the sweep; the rest per A1 |
| L1238 | `ADD_VOLUMETRIC_HEAT_SOURCE`, 5098-5146 | 0.024% | alloc_ptr_local 26, dt_alloc_component 6, reduction_accum 4 | host | none |
| L1239 | `RADIATION_FVM`, 3887-3893 (zero `Q_RAD_IN`) | 0.000% | alloc_ptr_local 8, pointer_reassign 3 (`B1_INDEX` as a value) | **translated** (`rad_wall_qin_zero`) | wall-table flag `B1_PRESENT` |
| L1240 | particle zeroing, 3894-3901 | 0.000% | alloc_ptr_local 10 | translatable with the particle table | particle table with a presence flag |
| L1241 | cface zeroing, 3902-3908 | 0.000% | alloc_ptr_local 8 | translatable (geometry) | cface table |
| L1243 | open-boundary `Q_RAD_IN`, 4965-4974 | 0.000% | pointer assignment `BR => BOUNDARY_RADIA(WC%BR_INDEX)`, `SUM` | **blocked**; hand-written specification kernel tested (`open_qin_spec`) | ragged per-wall table `BR_ILW(NRA,NSB,wall)` and a `BR` alias |
| L1244 | solid-particle flux, 4992-5037 | 0.001% | particle tree, alloc_ptr_local 38 | host | none |
| L1245 | `RADF` write, 5048-5063 | 0.002% | io 7, string_op 3 | host | none |
| L1246 | `GET_KAPPA` loop, 5174-5212 | 0.000% | reduction_accum, reduction_intrinsic | with the RADCAL `KAPPA_GAS` row | source-function callee |
| L1247 | `INTERPOLATE_IL` angle weights, 3623-3645 | 0.000% | write_indirect (`IDX`), reduction_intrinsic | host (run at new angle cycles only) | none |
| L1248 | `INTERPOLATE_IL` wall loop, 3651-3662 | 0.001% | `BR` alias, dt_alloc_component 12 | host, or table + per-thread scratch | ragged wall table, scratch of NRA values |
| L1249, L1250 | `INTERPOLATE_IL` cfaces 3665-3676, particles 3679-3695 | 0.000% | same | host | as L1248 |
| L1251 | `INTERPOLATE_IL` `OMESH%IL_R`, 3698-3713 | 0.000% | `OMESH` records | replaced by the FR-062 exchange | none |
| L1230-L1233 | `A_WSGG`, 5261-5312 | 0.000% | recursion, reduction_accum | with the WSGG row | recursion removal |
| L1252 | `KAPPA_WSGG`, 5217-5256 | 0.000% | reduction_accum | with the WSGG row | source-function callee |
| L1234-L1237 | `CALCULATE_DIRECTION_COEFFICIENTS`, 3431-3596 (loops 3510-3512, 3534-3542, 3547-3562, 3581-3592) | 0.000% | write_derived_index, reduction_intrinsic | host (startup, mesh-independent: nothing to port and nothing to rebuild at regrid) | none |
| L0658 | `INTERPOLATE1D` | 0.000% | exit_out_of_nest | with the droplet row | assumed-shape table |
| L0825-L0827 | `ALLOCATE_RADIATION_RECV_PKG` / `SEND_PKG` (main.f90) | 0.000% | MPI calls, allocation | host; replaced by the FR-062 exchange | none |
| L0006, L0007 | `CCCOMPUTE_RADIATION` (ccib.f90), cut faces | 0.000% (no GEOM timer in the pooled runs; use rank score) | dt_alloc_component | geometry-deferred | cut-face tables |

**Summary.** Translated now: 17 kernels (16 element-wise cell kernels and one wall kernel; Table A1 marks them) plus one hand-written specification for L1243. Together they are far below 1% of run time, and none of them is the cost centre: **the cost is the sweep inside L1242**, so the GPU benefit of radiation depends on the FR-062 kernel, not on these. The translated kernels still matter for a different reason: they are the non-sweep half of the stage, they all sit in the same routine, and after they run on the device the whole stage stays resident between sweeps. The `UIID` accumulation runs once per angle, so its launch count equals the sweep's; on the device it should be fused into the sweep kernel rather than launched alone.

## Part B. Review of the RTE_SOURCE change (upstream `afb5e31a48`)

**What changed.** Upstream added a per-band precompute of the angle-independent source term and replaced the inline expression in the sweeps by a read of it.

- Declaration `REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: RTE_SOURCE` (radi.f90:3798).
- Per band, inside `INTENSITY_UPDATE` (4210): `ALLOCATE(RTE_SOURCE, MOLD=KFST4_GAS)` (4220), then `RTE_SOURCE = KFST4_GAS + KFST4_PART + RSA_RAT*(SCAEFF+SCAEFF_G)*UIIOLD` (4221-4222). `DEALLOCATE` at 4942, so one allocate/free per band per call.
- Used in the cylindrical sweep (4465), the 2D sweep (4489) and the 3D sweeps (4628 STEP; 4672 and 4708 for the other schemes), each as `... VC*RSA(N)*RFPI*RTE_SOURCE(I,J,K)`. The upstream commit changed three sites (cylindrical, 2D, 3D STEP); the local merge converted the two FireX-only alternate-scheme sites in the same way.
- `ADD_VOLUMETRIC_HEAT_SOURCE` (4167) and the condensed-species block (4171-4202) run before it, so their contributions are included.

**Bitwise equality with the expression it replaced: yes, same grouping.** The inline text was `( KFST4_GAS(I,J,K) + KFST4_PART(I,J,K) + RSA_RAT*(SCAEFF(I,J,K)+SCAEFF_G(I,J,K))*UIIOLD(I,J,K) )` inside `RAP*(AIU_SUM + VC*RSA(N)*RFPI*( ... ))`. Fortran evaluates left to right at equal precedence, so both the old inline term and the new whole-array statement are `(KFST4_GAS + KFST4_PART) + ((RSA_RAT*(SCAEFF+SCAEFF_G))*UIIOLD)`, and in the sweep the factor `VC*RSA(N)*RFPI` multiplies the same value either way, `((VC*RSA(N))*RFPI)*RTE_SOURCE`. The only difference is that the sum is now stored to memory (a rounded double) and read back, instead of being held in a register or contracted: with `-ffp-contract=off` (the test flags) the value is identical; with default gfortran flags an FMA contraction of `RSA_RAT*(...)*UIIOLD + (A+B)` could occur in either form and may differ between the old inline expression and the precompute. The generated kernel is tested bitwise against the verbatim statement (`rad_rte_source`, all six flag sets). Not tested: the old sweeps against the new sweeps (the old text no longer exists in the tree to build); the argument above is by operation order.

**Array bounds and ghost cells.** `KFST4_GAS`, `UIIOLD`, `KFST4_PART`, `SCAEFF`, `SCAEFF_G` are `WORK1`, `WORK3`, `WORK7`, `WORK6`, `WORK9`, all allocated `(0:IBP1,0:JBP1,0:KBP1)` (init.f90:691-699; the preamble aliases are radi.f90:3832-3840), and not initialised at allocation. `MOLD=KFST4_GAS` therefore gives `RTE_SOURCE` the same bounds, ghost layers included. Every operand is fully defined before the statement: per band `KFST4_GAS`, `KFST4_PART`, `SCAEFF`, `SCAEFF_G` and `KAPPA_PART` are set to 0 as whole arrays (3914-3920) and then filled in the interior only, `UIIOLD` is assigned whole from `UII` or `UIID(:,:,:,IBND)` (4213-4216; `UII`/`UIID` are `(0:IBP1,...)` and initialised, init.f90:850-855, 901-902). So the ghost values of `RTE_SOURCE` are `0 + 0 + RSA_RAT*(0+0)*UIIOLD(ghost) = 0` (or `-0`/`NaN` only if `UIIOLD` held `Inf`/`NaN` in a ghost cell, which the initialisation does not allow). The old inline expression read only the interior cell; the new statement reads and writes the ghost layers too but never uses them: **the sweeps read only `RTE_SOURCE(I,J,K)` of the interior cell being solved** (cell ranges 1..IBAR etc., checked at the five use sites). The generated kernel covers `0:IBAR+1` etc. to match the whole-array statement; a launch over the interior only would also be correct for the sweeps and is cheaper.

**AMR implications.**

1. *Allocation.* In GPU code do not allocate and free per band. `KFST4_GAS` and the other terms are already driver work arrays (`WORK1`..`WORK9`); `RTE_SOURCE` should be one more driver-owned work array per box, allocated with the box (same shape and ghost width as `WORK1`) and reused for all bands and passes. It is not a per-thread scratch (the scratch-pointer mechanism in `amrex/scratch-pointer-support.md` is for fixed-size private arrays), it is a box-sized field.
2. *Halo.* No ghost value is used, so no fill, no exchange and no coarse-fine treatment is needed for `RTE_SOURCE`. At a box face or a coarse-fine face the sweep takes upwind intensity from `IL` ghost cells (the FR-062 face buffers), never from the source.
3. *FR-062 per-box lagged sweep.* The source is angle-independent and box-local: compute it once per band per box before the K passes' angle loops. It depends on `UIIOLD`, which is the previous pass's `UII` or `UIID`: that is the lagged quantity of the scheme, so `RTE_SOURCE` is rebuilt at the start of each pass (each call of `RADIATION_FVM`) from the box's own previous-pass values, with no data from a neighbour box. This keeps the pass independent of box ownership and order (FR-062 rule (1)).
4. *Memory.* One 3-D array per box: for a 32³ box with ghost layers, 34³·8 B ≈ 0.31 MB; it replaces nothing (the inline form needed none). Negligible next to `IL`, `UIID` and the exchange buffers.
5. *Regrid.* Nothing persistent: it is rebuilt each call. No transfer, no interpolation.
6. *Determinism.* Elementwise, no reduction, no atomics: byte-identical across ranks and threads.

**Conclusion.** The change helps. It turns the source into a clean elementwise kernel that is separate from the sweep (translated and tested here), removes the repeated recomputation inside each of the NRA angles (the old form recomputed the sum per cell per angle), and adds no new obstacle for the AMR design. The one cost for the port is a host-side mechanical one: `ALLOCATE(MOLD=)` is Fortran 2008 and the shared generator front end parses with `std="f2003"`, so `radi.f90` at this commit does not parse until the front end handles it (see the notes).

## Updates to the spec

`docs/radiation/01-radiation-amr-spec.md` cites the merged-file line numbers for all `radi.f90` lines above 3797 (re-pin note in its header) and §1 lists `RTE_SOURCE`.
