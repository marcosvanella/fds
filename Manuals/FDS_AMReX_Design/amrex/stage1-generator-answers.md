# Stage 1 GPU spike: answers to the generator questions (Q1 to Q4)

Answers the "Legacy Mapper (generator)" questions in `docs/amrex/stage1-gpu-spike-plan.md` section 9 (Q1 loops beyond the 34, Q2 box lo/hi offsets, Q3 reductions, Q4 flat tables). Source line numbers are at the inventory baseline `36975d7` (Fortran `Source/`) unless marked otherwise; the current FDS-AMReX line has `#ifdef WITH_AMREX` guards that shift some lines in `velo.f90` and `main.f90`. Loop ids and time shares come from `docs/inventory/gpu_generator_loop_classes.csv` (modelled share of total run time, percent). **All ETAs are estimates in working days for one engineer, not measurements.** Nothing here is implemented.

## Q1. Loops a periodic step needs beyond the 34 kernels

The 34 kernels cover about 1.2% of modelled time (`docs/amrex/stage1-gpu-spike-plan.md` section 1). The loops the plan names add up to about 7.2% (all rows below, 7.153). Status is from the loop-class CSV; "front end ok" means the front end accepts the loop but no kernel or bitwise test exists yet.

| Routine | Loops | Lines | Time % | Status and blocker |
|---|---|---|---|---|
| `VELOCITY_FLUX` | L1374 | velo.f90:619-639 | 0.079 | front end ok, untested |
| `VELOCITY_FLUX` | L1375 | velo.f90:649-653 | 0.001 | callee reads derived-type components |
| `VELOCITY_FLUX` | L1376, L1377, L1378 | velo.f90:662-714, 720-772, 778-830 | 0.226 each | loop body reads `CELL(IC)%EDGE_INDEX(n)`, `EDGE(IE)%OMEGA(+-1,+-2)` and `%TAU` (velo.f90:678-697): no edge tables exist |
| `VELOCITY_FLUX_CYLINDRICAL` | L1391 ok; L1392, L1393 blocked | velo.f90:1243-1253; 1270-1302, 1306-1337 | 0.026; 0.014 each | same `CELL%EDGE_INDEX` access |
| `MASS_FINITE_DIFFERENCES` | L0880 | mass.f90:65-192 | 4.832 | the flux nest, the largest single loop here. The species `N` loop remaps pointers (`FX_P(LBOUND(FX,1):,...) => FX(:,:,:,N)`, mass.f90:81-83) and calls `GET_SCALAR_FACE_VALUE`; the wall nest (mass.f90:93-188) has `Z_TEMP` array constructors and a pointer callee. The inner `RHO_Z_P` product is already translated (`rho_z_p_mass`) |
| `MASS_FINITE_DIFFERENCES` | L0882 | mass.f90:224-320 | 0.302 | whole-array assignment inside the loop body |
| `MASS_FINITE_DIFFERENCES` | L0881, L0883 | mass.f90:201-209, 328-350 | 0.071; 0.128 | L0883 translated and tested; L0881 front end ok, untested |
| `GET_SCALAR_FACE_VALUE` | L0652 to L0656 | func.f90:1348-1435 | 0.005, 0.016, 0.053, 0.053, 0.074 | classified translatable; front end stops at the `POINTER` dummy arguments ("pointer U associated outside the routine") |
| `GET_SCALAR_FACE_VALUE` | L0657 | func.f90:1439-1453 | 0.074 | whole-array assignment (MP5 branch) |
| `VELOCITY_PREDICTOR` | L1394 to L1396 | velo.f90:1603-1629 | 0.005 each | front end ok, untested |
| `VELOCITY_PREDICTOR` | L1397, L1398 | velo.f90:1648-1665 | 0.049 each | function reference `VD2D_MMS_U/V`. **Dead code**: guarded by `PERIODIC_TEST==7 .AND. .FALSE.` (velo.f90:1647) |
| `VELOCITY_CORRECTOR` | L1369 to L1371 | velo.f90:1724-1750 | 0.005 each | front end ok, untested |
| `VELOCITY_CORRECTOR` | L1372, L1373 | velo.f90:1770-1787 | 0.049 each | same MMS function reference, same dead guard (velo.f90:1769) |
| `DIVERGENCE_PART_2` | L0395, L0396 | divg.f90:1611-1630 | 0.012 each | front end ok, untested |
| `DIVERGENCE_PART_2` | L0394 | divg.f90:1574-1604 | 0.042 | alias `B1=>BOUNDARY_PROP1(WC%BC_INDEX)` has no table; the loop writes `DP` at internal-solid cells and copies `DP` at the gas cell (read/write overlap to analyse) |
| `CHECK_STABILITY` | L1347 | velo.f90:3059-3080 | 0.074 | **translated and bitwise-tested** (`cfl_max`, argmax with the serial `>=` tie rule, `MAXVAL` expanded; see Q3 and `blocked-loop-families.md` V1) |
| `CHECK_STABILITY` | L1348 | velo.f90:3093-3108 | 0.014 | **translated and bitwise-tested** (`cfl_wall_max`; continues the accumulator of the cell pass) |
| `CHECK_STABILITY` | L1349 | velo.f90:3121-3134 | 0.037 | **translated and bitwise-tested** (`vn_max`; `CYCLE I_LOOP` is the kernel's own innermost loop) |
| `BAROCLINIC_CORRECTION` | L1343 to L1346 | velo.f90:3254-3305 | 0.011 each | front end ok, untested |
| `CHECK_DIVERGENCE` | L0363 | divg.f90:1675-1711 | 0.272 | `RESMAX` self-update and `CYCLE LOOP1` out of a nest (listed for completeness; not named in the plan) |

The plan's "VELOCITY_FLUX: 4 of 5 loops not translatable" is L1375 to L1378 above (L1374 is the one that is accepted).

**Proposed order** (by cost of unlocking, then share):
1. *Front end ok, add markers and bitwise tests only*: L1374, L1394-L1396, L1369-L1371, L0395, L0396, L0881, L1391, L1343-L1346. About 0.27% of time in 17 loops. Estimate 3 working days.
2. *`GET_SCALAR_FACE_VALUE` (0.275%)*: make the `POINTER` dummies explicit-shape arguments through the callee-directive route already used for 4 callees, handle the MP5 whole-array line (L0657) by an element loop. Estimate 3 to 4 days. This unlocks the species `N` loop of L0880.
3. *`MASS_FINITE_DIFFERENCES` L0880 (4.832%)*: split the cell flux nest from the wall nest as for `DIVERGENCE_PART_1`; replace the `FX_P` remap with a species-index argument and the `Z_TEMP` constructors with scalar locals. Estimate 4 to 6 days, the biggest item. L0882 (0.302%, whole-array assignment to element loops) 1 to 2 days.
4. *`CHECK_STABILITY` and `CHECK_DIVERGENCE`* (0.397% together): `CHECK_STABILITY` (0.125%, three loops) is done: about one working day including tests. `CHECK_DIVERGENCE` L0363 (0.272%) is done in round 7 (`red_extrema.py`, kernel `div_extrema`). It has a `RESMAX` self-update, a `CYCLE LOOP1` out of the nest and a minimum with `<`; the argmax builder handles one maximum with `>=`, so a second builder handles three accumulators and a `<` rule. Took about one working day with tests.
5. *Edge-table loops L1376 to L1378, L1392, L1393* (0.706% together): need a decision on edge flat tables (`EDGE_INDEX`, `OMEGA`, `TAU`) the way the wall tables were decided. Estimate 2 days for the design note plus 4 to 5 days to implement and test.
6. *L0394 (0.042%)*: needs the alias-to-table decision for `B1`; do it with the wall-table owner. 1 to 2 days.
7. *L1397, L1398, L1372, L1373 (0.196%)*: dead code. Leave on the host (no work); drop them from the denominator of "still to translate".

Owners: the GPU Generator Engineer takes items 1 to 4 and 6; items 5 and 6 need the Integration Lead for the table layout, as for the wall tables. The unblock list for the periodic step is therefore about 25 to 30 working days of estimated effort, which is a planning figure only.

## Q2. Ruling (b): lo/hi offsets as kernel arguments

**Status (owner decision, see `generator-decisions.md`, decision A): the lower-bound argument change is deferred. Nobody starts it; it goes back to the CEO only when a loop truly needs it.** The proposal below is kept as the record of what it would cost. The working assumption is now "bounds change not done".
**Consequence, plainly:** with the kernels as they are, a mesh split into several boxes works only if the driver hands each box to the kernel as its own array with local origin 1 (box-local arrays with the FDS ghost widths, box-local wall tables, sliced 1-D metric arrays, box-local coordinates for position-dependent callees). A box whose lower corner is not 1 cannot be run in place. Whole-mesh (single-box) runs are unchanged. This is the "What therefore works with no change" case below.

**What the generator does today** (`amrex/s4_mass/s5_gen/s5gen.py`, generated kernels in `generated/`):
- Mesh extents are integer arguments: `IBAR`, `JBAR`, `KBAR` are always added to the argument list (s5gen.py:1275-1289). They are declared `integer(c_int), intent(in)` in every kernel (for example `s5gen_zzs_pred(IBAR,JBAR,KBAR,NS,NL,DT,...)`, `generated/s5gen_k2.F90:322-324`), so they are run-time values, not constants.
- Arrays are explicit-shape dummies whose bounds are the FDS bounds, taken from the `ALLOCATE(M%X(...))` text (layout "exact", s5gen.py:1669-1690; IBP1 is renamed to `IBAR+1`, around s5gen.py:471-490): for example `UU(-1:IBAR+1,0:JBAR+1,0:KBAR+1)`, `RHO(-1:IBAR+2,-1:JBAR+2,-1:KBAR+2)`, 1-D metrics `RDX(0:IBAR+1)` (`generated/s5gen_k2.F90:327-330`). The round-1 default ("amrex0") is `0:IBAR+1,0:JBAR+1,0:KBAR+1` (s5gen.py:1696-1702).
- Loop bounds are copied from the upstream text (`DO K=1,KBAR`). In the 34 kernels `IBAR/JBAR/KBAR` appear only in dummy declarations and in 66 `DO` bound expressions; no kernel body uses them otherwise (checked by search of `generated/s5gen_k2.F90`).

**What therefore works with no change.** A box of any size `nx,ny,nz` can be run if the driver (a) passes `IBAR=nx, JBAR=ny, KBAR=nz`, (b) passes FABs whose memory layout equals the dummy shape, that is the same ghost widths as the FDS allocation (low side 1 or 2 cells as in the declarations above), and the data pointer at the lowest declared element, and (c) treats the box as having its own local index space with first cell index 1. The 1-D metric arrays must then be passed as a pointer offset into the mesh-wide array. Tables keyed by cell index (wall tables `BC_IIG`, `BC_JJG`, `BC_KKG`) must hold box-local indices. This is the "single box with lower corner 1" case, which is what the WP0 and WP1 harnesses use.

**What does not work.** (1) A kernel cannot be told the global index of its box, so the driver must build box-local wall tables and slice every 1-D array (`R`, `RDX..RDZN`, `X`, `XC`, `Z`, `ZC`); any host or device error in that slicing is silent. (2) Functions of position called inside a kernel (for instance the MMS functions or position-dependent callees) get box-local coordinates. (3) Ghost layers of an interior box face must be filled by `FillBoundary` before the kernel; a dummy declared with FDS bounds has no check that the FAB ghost width matches (a mismatch gives a shifted read with no diagnostic).

**Smallest change that satisfies ruling (b)** (proposal, option A; deferred, not started): keep FDS (mesh-global) indexing, and replace the mesh extents by box bounds. Add six integer arguments `ILO, IHI, JLO, JHI, KLO, KHI` (global FDS indices of the box's first and last interior cell; `1` and `IBAR` for a whole mesh) and rewrite only the literal forms the generator already knows:
- dummy lower bound literal `c` becomes `ILO+(c-1)`; upper bound `IBAR+c` becomes `IHI+c` (same for J, K);
- `DO` bounds `c` and `IBAR+c` follow the same rule;
- wall tables stay mesh-global, so the driver only has to select the walls whose gas cell is in the box (already the design in `docs/adr/drafts/gpu-generator-design.md`, "Interaction with the AMReX box layout").
Because `ILO=1`, `IHI=IBAR` reproduces the current text exactly, single-box numerics are unchanged, and the existing 408-case driver and the WP0/WP1 outputs are the regression set. Extra tests needed: the same kernel on a sub-box of a larger mesh, compared with the full-mesh run on that sub-box (bitwise), plus a ghost-width mismatch case.

*Effort estimate*: generator change about 2 days (bounds substitution in `exact_dims`/`dims_of` and the loop emitter, 6 arguments in `s5gen_args.json` and the header), regenerate and re-check goldens 1 day, new sub-box tests 1 to 2 days: **about 4 to 5 working days** to a build the harness can take. Option B (index-offset argument, box-local indexing kept) is smaller in the generator (about 2 days) but moves the offset bookkeeping into the driver and every table, which is the failure mode above; not recommended.

A residual condition: loops whose range is `0..IBAR+1` (ghost cells) now run on box ghost cells; this is correct only when the ghost data are valid there, as `MESH_EXCHANGE` made them in FDS. Domain-edge special cases (an `IF (I==IBAR)` style test) do not occur in the current 34 kernels but will in the `VELOCITY_*` and `CHECK_*` loops, which then need the mesh extent as well as the box bounds.

## Q3. Order-independent reductions (proposal; owned by the GPU Generator Engineer)

What exists: ruling (c) in the plan (min and max as they are; sums by a fixed-order per-box tree, partials combined in global box-index order; not bitwise against the FDS CPU path). The marker attribute `REDUCTION(op:var)` is parsed and emitted as an OpenMP reduction clause, no loop uses it (`docs/inventory/gpu_generator_coverage.md` section 5). `ExactSum.H` (`Source/pressure_backend/ExactSum.H:1-8`) already implements a decomposition-independent sum with 128-bit integers on a fixed power-of-two grid; plan section 2 notes it is unchecked on the device.

Options:
1. *Max or min with a tie rule* (`CHECK_STABILITY` L1347, L1349, `CHECK_DIVERGENCE` L0363). The value is order independent. The recorded index is not: `IF (UVW>=UVWMAX_TMP)` (velo.f90:3072-3077) keeps the last tied cell in loop order within a thread, and the merge in the OMP CRITICAL (velo.f90:3082-3089) uses `>`, so across threads the first merged thread to arrive wins; the FDS index on a tie therefore depends on the thread count. Rule adopted and implemented (`s5_gen/red_argmax.py`, tested): the device result equals the SERIAL loop, which is a lexicographic maximum of (value, linear iteration index in K,J,I order) over the iterations that pass `V >= ACC`, with the corrections found by the tests: the incoming accumulator takes part with `>=` (a cell equal to it still writes its location, an incoming value above every cell leaves everything untouched), NaN never wins, the stored value is the winner's own value (a tie between -0 and +0 keeps the winner's sign), and the wall pass continues the accumulator of the cell pass (a wall beats an equal cell). It is not the first-merged-thread rule of the threaded FDS merge, and the strict `>` merge with the initial 0 (velo.f90:3085) is left to the caller: an all-zero field leaves `ICFL` unchanged there, but writes the last fluid cell in the `>=` loops L1347 (cell pass) and L1349. One kernel call covers one box: several boxes must be merged by (value, GLOBAL linear index), not by box order. Implementation: three passes (max of values, max of the linear index among the equal ones, one iteration re-runs the body at the winner), no atomics. No floating-point sum is involved.
2. *Sums over zones and mass* (`DSUM/PSUM/USUM(IPZ)` divg.f90:739, 750, 765-766; `DELTA_RHO_ZZ`, `CELL_COUNTER` families). **Owner decision (see `generator-decisions.md`, decision B): the default keeps the FDS order. The GPU computes the per-cell and per-wall terms (elementwise, bitwise equal to FDS) and the additions are serial in FDS order, on the host or in a single-thread device pass; this is bitwise equal to FDS for one process. In addition a compile-time or macro switch selects exact fixed-point sums (option 2a below) on the GPU; it is off by default, tested separately, independent of box layout and thread count, and documented as not bitwise equal to FDS. The GPU Generator Engineer implements the switch.** Option 2a as originally recommended: fixed-point exact accumulation because the result is independent of box layout and thread count, which is what the plan's acceptance needs across layouts; use the per-box fixed-order tree of ruling (c) (2b) only where exactness cannot be applied. Both change the summation order against FDS.
3. *Compensated summation*: reduces error but is not order independent; not recommended.
4. *OpenMP or atomic reduction*: order depends on the schedule; not acceptable for the acceptance rule, use only for debug.

Recommendation (amended by decision B): option 1 for the maxima; for sums the default is the serial FDS-order add of GPU-computed terms, and 2a is the optional switched path. For the switched path, use a host reference that sums in the same fixed point so the device result is compared bit for bit with the host fixed-point result (the FDS floating-point value is compared with tolerance only). The device `__int128` support must be checked before 2a is relied on; fallback is two 64-bit limbs. Per-family detail is in `docs/amrex/blocked-loop-families.md`.

## Q4. Who fills and refreshes the flat tables on the device

What is written down (`docs/adr/drafts/gpu-generator-design.md`, section "WALL flat tables", lines 43-48 and 80):
- *Owner*: the driver shim (AMReX side) allocates and owns every flat table in device memory, sized by the box's wall count, and passes them as arguments. The generator emits only the kernel and the argument manifest and never copies.
- *Static between geometry changes*: `BC_II/JJ/KK/IIG/JJG/KKG/IOR` (init.f90:3330-3339, 3437-3443), `W_THIN` (init.f90:3315), `EW_NIC` (main.f90:2111), `MW_SPEC`. Re-gather after `REASSIGN_WALL_CELLS` (init.f90:4898, called from main.f90:1793) and when `WALL` is reallocated (func.f90:4003-4014, steps of 1000).
- *Geometry-event tables*: `W_BOUNDARY_TYPE` (init.f90:5076-5107, 5192; main.f90:2772), same re-gather.
- *Time-varying, one writer per step*: `EW_BOUNDARY_TYPE_PREVIOUS` (velo.f90:2667), `B1_RHO_D_DZDN_F` (divg.f90:223-228 and wall.f90:946, 952), `UVW_SAVE`/`U_GHOST`/`V_GHOST`/`W_GHOST` (velo.f90:2732-2832 and the deferred ccib.f90 writers), `UN_WALLS`. The shim scatters host-written fields in before a kernel that reads them and scatters back after a kernel that writes them only if host code reads them.
- *Known risk*: a stale `W_BOUNDARY_TYPE` after an obstruction event gives a wrong result with no diagnostic; the design note advises a step-level checksum of the tables against the host records.

What is **not decided**:
- The species and property tables `MU_RSQMW_Z`, `K_RSQMW_Z`, `CP_Z`, `H_SENS_Z`, `RSQ_MW_Z`, `MW` (kernel arguments such as `MU_RSQMW_Z:table(0:I_MAX_TEMP,NS)`, `markers/golden_signatures.json:597`; declared at cons.f90:489 for `MU_RSQMW_Z`). The design note lists only `MW_SPEC` (type.f90:510). The bitwise harness fills them synthetically (`test/make_r2_tests.py:350-371`). Proposal: they are set at initialisation and constant (no reaction or species change at run time), so the driver uploads them once per run and again only if `I_MAX_TEMP` or the species set changes; owner is the Integration Lead's driver shim with the generator listing each table in `s5gen_args.json`.
- A refresh trigger signal: today it is implicit in the host call to `REASSIGN_WALL_CELLS`. Proposal: a counter incremented by that routine (patch file in `docs/upstream-patches/`, not an edit of upstream) that the shim compares each step.
- The staleness check (checksum of static tables per step) is advised but not specified or implemented.
- Edge tables (Q1 item 5) have no entry yet.
