# L0398 pre-analysis: the wall loop of `ENTHALPY_ADVECTION_NEW` (divg.f90:835-939)

Status: analysis only. No kernel, marker entry or test exists for this loop. The loop is 0.260 % of the modelled time (loop work list, row L0398). It belongs to the
packages S2 (pointer scratch and array constructors) and S3 (wall scatter through the gas-cell subscript); see `blocked-loop-families.md`. Line numbers are those
of the pinned upstream sources used throughout these notes (`36975d7`).

## 1. What the loop does

`ENTHALPY_ADVECTION_NEW(U_DOT_DEL_RHO_H_S)` computes the advection term of the enthalpy equation, `u . grad(rho h_s)`, in five steps:

1. (divg.f90:814-822, loop L0397, a cell loop, separate work item) `RHO_H_S_P = RHOP * h_s(ZZP, TMP)` for every cell including one ghost layer.
2. (divg.f90:829-831) three whole-field calls `GET_SCALAR_FACE_VALUE(UU,RHO_H_S_P,FX_H_S,...)` (and VV/FY, WW/FZ) fill the face values `FX_H_S`, `FY_H_S`, `FZ_H_S`.
3. **The loop of this note** (divg.f90:835-939, `WALL_LOOP`), once per wall cell, in ascending wall index `IW`, serial in the source (there is no OpenMP directive on it):
   * skips `NULL_BOUNDARY` walls;
   * **(a) off-wall face overwrite** (859-916), skipped for `INTERPOLATED_BOUNDARY` and `OPEN_BOUNDARY`: when the flow next to the wall points away from it and the
     face is not itself a wall face, it recomputes the first face off the wall with the one-face form of the limiter and overwrites `FX_H_S`, `FY_H_S` or `FZ_H_S` there.
     Six cases by `BC%IOR` (+-1, +-2, +-3). Each case: a guard on the sign of the face velocity and on `CELL(...)%WALL_INDEX(+-n) > 0`, a strip of `RHO_H_S_P` copied
     into the scratch array `Z_TEMP` by an array constructor, the velocity into `U_TEMP(1,1,1)`, `CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,IOR,LIMITER)`,
     and the store of `F_TEMP(1,1,1)`;
   * **(b) wall term** (840-857 and 919-937): the sensible enthalpy `H_S` of the wall-face state (`GET_SENSIBLE_ENTHALPY(B1%ZZ_F, ..., TMP_F_GAS)`, with the gas
     temperature used for a solid wall with outflow), the normal velocity `UN` (read from `UU/VV/WW`, from `B1%U_NORMAL_S`/`B1%U_NORMAL` for a solid wall, from
     `UVW_SAVE(IW)` for an interpolated wall), `DU = (B1%RHO_F*H_S - RHO_H_S_P(gas cell))*UN`, and
     `U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) = U_DOT_DEL_RHO_H_S(IIG,JJG,KKG) - SIGN(1,IOR)*DU*B1%RDN`.
4. (946-966) a cell loop adds the interior flux differences on top: `U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + (DU_P-DU_M)*RDX + ...`, using
   `WALL_INDEX(n)==0` to skip faces that are walls. The wall terms are therefore the first terms of each sum and the cell term is added last.

The two halves of the loop body are independent: (a) writes only the face arrays, which nothing in the wall loop reads; (b) writes only `U_DOT_DEL_RHO_H_S`, which
nothing in (a) reads. They may be run as two loops.

## 2. Data accessed

Read: `WALL(IW)%BOUNDARY_TYPE`, `BC%II,JJ,KK,IIG,JJG,KKG,IOR`; `B1%U_NORMAL_S`, `B1%U_NORMAL`, `B1%TMP_F`, `B1%ZZ_F(1:N_TRACKED_SPECIES)`, `B1%RHO_F`, `B1%RDN`; the
module flag `PREDICTOR` (and `CORRECTOR`); `UVW_SAVE(IW)`; `TMP(IIG,JJG,KKG)`; `UU`, `VV`, `WW` (module pointers set in `DIVERGENCE_PART_1`); `RHO_H_S_P` (`WORK_PAD`,
bounds -1:IBP1+1); `CELL(CELL_INDEX(i,j,k))%WALL_INDEX(+-1..3)`; the table `H_SENS_Z(0:I_MAX_TEMP, 1:N_TRACKED_SPECIES)` and the integer `I_MAX_TEMP` (inside
`GET_SENSIBLE_ENTHALPY`); the integer `I_FLUX_LIMITER`.
Written: `FX_H_S`, `FY_H_S`, `FZ_H_S` (`WORK2..WORK4`), at the first face off the wall; `U_DOT_DEL_RHO_H_S` (`WORK6`), at the gas cell of the wall. Private: the pointers `WC, BC, B1`, the
scratch arrays `U_WORK, Z_WORK (-1:3)^3`, `F_WORK (0:3)^3` and the pointers aimed at them, the scalars `UN, UN_P, TMP_F_GAS, DU, H_S`, the allocatable `ZZ_GET(1:N_TRACKED_SPECIES)`.

## 3. What blocks the generator today

| # | blocker | lines | state |
|---|---|---|---|
| 1 | pointer remap `U_TEMP=>U_WORK` etc. and uses of the scratch pointers | 863-865 | `legacy_ptr.py` rewrites them to private scalars `US`, `FS` and `ZS(0:3)` (`scratch = true` in the sidecar); written for the nests of L0401/L0403, not run on this loop |
| 2 | array constructors fill three of the four elements of `Z_TEMP` | 872, 882, 889, 896 and the y, z cases | pad to 0 in the rewrite (`Emitter.assign_constructor`); the fourth element is only read by MP5, and upstream patch UP-0002 would pad it in the source (not applied upstream) |
| 3 | pointer callee `GET_SCALAR_FACE_VALUE` called for one face | 874 etc. | the one-face callee `GET_SCALAR_FACE_VALUE_PT(A,Z,F,LIMITER)` from `s5_gsfv.py` exists; its registration (and the `VALUE` dummy intent in `leaf_events`) is open |
| 4 | `CELL(CELL_INDEX(a,b,c))%WALL_INDEX(n)` | 871, 881, ... | table `WALL_INDEX(0:IBAR+1,0:JBAR+1,0:KBAR+1,-3:3)` exists (policy.arrays) |
| 5 | accumulate through the wall's gas-cell subscript | 936-937 | needs the cell-to-wall table in ascending wall index and a per-gas-cell gather (design of the Wall Loops Engineer; the same family as L0375, L0394, L0405) |
| 6 | derived-type pointers `WC`, `BC`, `B1` | 837-840 | flat per-wall tables (`W_BOUNDARY_TYPE`, `BC_II`.., `BC_IIG`.., `BC_IOR`, `B1_RHO_F`, `B1_TMP_F`, `B1_RDN`, `B1_ZZ_F`) exist in the driver of the round-2 tests and in the wall family; `B1%U_NORMAL_S` and `B1%U_NORMAL` need two more flat tables |
| 7 | `GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP)`: an allocatable private array, an array-section copy `ZZ_GET(1:N)=B1%ZZ_F(1:N)` and `DOT_PRODUCT` over table sections | 849-857 | the same gap as L0397 (open, needs a feature); a hand-written device callee would do (section 5) |
| 8 | `UN` is a scalar assigned in a `SELECT CASE` with a default branch, then used | 919-934 | fine as long as every path assigns it (see risks) |
| 9 | the loop is serial upstream | 835 | not a blocker, but it fixes the reference order of the accumulate |

Nothing here is specific to a missing language feature except 5 and 7; 1 to 4 and 6 are covered by the S2 machinery once it is run on this loop.

## 4. What it needs from the callee and the constructor machinery

* One mechanism pads the unassigned fourth element of the stencil, and it sits in the callee (`GET_SCALAR_FACE_VALUE_PT`, MP5 branch). The constructor rewrite must not add a
  second pad. A pad written as `0._EB` in both places would be harmless, but the test below checks that only one exists.
* The constructor rewrite must fill the axis `|IOR|` of the scratch strip, `ZS(0:2)` for the positive and `ZS(1:3)` for the negative directions, exactly as for L0401.
* The one-face callee gets the limiter as an argument (`I_FLUX_LIMITER`); no per-limiter kernels are needed for step (a). The six whole-field limiter kernels serve only
  step 2 above (the `field` call site of `callsite/s5_gsfv_dispatch.F90`, routine `s5_gsfv_faces3`, which takes the rank-3 fields of this routine).
* The `WALL_INDEX` table and the per-wall tables from the wall family; the cell-to-wall list in ascending wall index, which must not depend on the decomposition into boxes.

## 5. Proposed kernel shape

Two kernels, in this order, after the cell loop of L0397 and the three whole-field calls, and before the final cell loop:

1. **`enth_wall_offface`** (step a). One thread per wall, `S5_LOOP_WALL` region. For each wall that is not `NULL`, `INTERPOLATED` or `OPEN`: the six-case guard, the strip into
   `ZS(0:3)` with the pad, `CALL GET_SCALAR_FACE_VALUE_PT(US,ZS,FS,LIM)`, one store to the face arrays. Stores are unique per face by the sign guards; the face-write check
   of `massflux-wall-table-spec.md` section 4 (the O-writer rules, no W-writers here because this loop writes no wall-face value) is run once per table rebuild.
2. **`enth_wall_gather`** (step b). One thread per gas cell that has walls, over the cell-to-wall list in ascending wall index. The thread starts from the value already in the cell
   (zero after the initial clear), and for each of its walls in list order computes `TMP_F_GAS`, `H_S`, `UN`, `DU` and subtracts `SIGN(1,IOR)*DU*RDN` from a register, writing the
   result back once. This reproduces the serial order of the source exactly, with no atomics. The expression is kept verbatim (`SIGN(...)*DU` first, then `*RDN`, then the subtraction).
   `H_S` comes from a small device callee `s5_sens_enthalpy(htab, imax, n, z, tmp)` that keeps the source's order: `ITMP = MIN(I_MAX_TEMP-1, INT(TMP))`, the two dot products as
   sequential sums over the species index (the second over the differences `H(ITMP+1,:)-H(ITMP,:)`), compiled without fused multiply-add. Until the generator takes L0397's callee, this
   is a hand-written file in the same style as the other callee files, not generated text.

Generator side: the body of the original loop has to be emitted twice with different statement subsets (a sub-nest, as the round-8 kernels of `DIVERGENCE_PART_1`). The split is legal
because of the independence stated in section 1. If the generator cannot split a loop, the alternative is one gather kernel that also performs step (a) per wall inside the
same list pass; this changes nothing numerically but mixes a scatter store into the gather, so the split is preferred.

Arguments (sketch): `IBAR,JBAR,KBAR, NWALL, PREDICTOR, LIMITER, N_TRACKED_SPECIES, I_MAX_TEMP`, the wall tables listed in section 2, `UU,VV,WW,TMP,RHO_H_S_P,UVW_SAVE`, `WALL_INDEX`,
`H_SENS_Z`, the list arrays, and `FX_H_S,FY_H_S,FZ_H_S,U_DOT_DEL_RHO_H_S` as explicit-shape or assumed-shape dummies with the bounds of `WORK_PAD` (-1:IBP1+1).

Scope under the neighbour-mesh ruling (D-055, no kernels for the `NOM>0` branches): this loop has no `NOM>0` branch of its own. Its `INTERPOLATED_BOUNDARY` walls (the walls
shared with another mesh) are excluded from step (a) and take `UN = UVW_SAVE(IW)` in step (b). Step (b) for such a wall is translatable if the host fills `UVW_SAVE` first
(it is written by the host velocity matching and is a plain input array here). A single-mesh run has no such walls. Tests use `UVW_SAVE` as an input; no neighbour-mesh data is touched.

## 6. Risks

1. **MP5 reference.** The verbatim loop reads a stale fourth element with MP5, so the bitwise reference for MP5 has to be the loop with the pad of patch 0002. Limiters 0 to 4
   are not affected. If the patch is not applied upstream, MP5 results of the device code differ from the unpadded baseline by construction.
2. **Order of the accumulate.** Cells with more than one wall (corner cells, thin obstructions seen from both sides) are summed in ascending wall index; a different list order changes
   the last bits. The list must be built in wall order, must exclude `NULL` walls (the source skips them, and including them with a zero term would not change a sum but would change
   the count the checks use), and must not depend on the decomposition. Start the sum from the cell's current value, not from zero, so a later change that adds earlier terms stays exact.
3. **Floating point.** `DOT_PRODUCT` order and fused multiply-add in `H_S`; the build must use no-FMA flags on host and device. `INT(TMPG)` truncation and `MIN` are kept as written.
4. **Stale scalar `UN`.** `UN` is assigned on every path only if `BOUNDARY_TYPE` is one of the covered cases and (for solid walls) exactly one of `PREDICTOR`, `CORRECTOR` is true. The kernel
   must assign `UN` from a defined value for the uncovered combinations or refuse them in the table check; the host reference leaves the previous wall's value, which a per-wall kernel
   cannot reproduce.
5. **Two walls naming one off-wall face.** A gap two cells wide gives two writers (left wall and right wall) whose sign guards exclude each other; a one-cell gap can put an off-wall target on
   another wall's face, which the `WALL_INDEX` guard is meant to remove. The face-write check must run on the built tables, and the test must contain both gap sizes.
6. **Bounds.** The strip reads `RHO_H_S_P(II-2)` and `RHO_H_S_P(II+2)`: the arrays must carry the bounds of `WORK_PAD`, and a mesh with `IBAR = 1` is a boundary case (`II+2 = IBP1+1`).
7. **`GET_SENSIBLE_ENTHALPY` is shared with L0397.** The two should use one callee, and its change of status (generated or hand-written) changes both loops.
8. **Size of the gain.** 0.260 % of the modelled time; most of it is the wall loop at about one walls-worth of work per step. The effort is mostly shared with L0401/L0403 (step a)
   and L0375/L0394/L0405 (step b), so the new work is small once those exist.

Effort estimate, assuming the cell-to-wall list, the scratch machinery run on L0401 and the L0397 callee exist: about 0.5 day for step (a), 1 day for step (b), 1 to 1.5 days for the test
driver and mutants, 0.5 day for the device run. Without the list or the callee, the loop waits for them.

## 7. Bitwise test plan

Reference: the verbatim text of steps 2 to 4 of the routine (the three whole-field calls, `WALL_LOOP`, the final cell loop), compiled and run on the host with the upstream
pointer-based callee; for MP5 the same text with the pad of patch 0002.

* Data: random `UU, VV, WW, TMP, RHO_H_S_P` (with values of very different magnitude, so that a changed summation order shows in the last bits), random wall sets with all
  boundary types (`NULL`, `SOLID`, `OPEN`, `INTERPOLATED` with a host-supplied `UVW_SAVE`), external walls on all six faces, internal obstructions with thin walls, gaps of one and
  of two cells, corner cells with up to three walls, gas cells with no wall, `PREDICTOR` true and false, solid walls with `U_NORMAL` of both signs and zero, `IBAR` = 1, 2, 3 and larger.
* Compared: `FX_H_S`, `FY_H_S`, `FZ_H_S` and `U_DOT_DEL_RHO_H_S` over the whole array including the ghost layer (a sentinel in untouched elements), after step 3 and after step 4.
* Limiters 0 to 5, 6 flag sets (O0, O2, O0 with OpenMP, O2 with OpenMP, O2 with OpenMP and offload disabled, O2 with OpenMP and DPD directives), 1, 4 and 8 threads.
* Mutants that must be caught: wall list in descending order or sorted by cell; the sum started from zero instead of the cell's value; `SIGN` factor dropped or applied after `*RDN`;
  a guard sign flipped; `WALL_INDEX` test removed; one face written twice; pad non-zero; `ITMP` rounded instead of truncated; `H_S` computed with the wall temperature instead of the gas
  temperature for outflow; `UN` taken from the wrong side (`II-1` vs `II`).
* Refusals that must be reported and write nothing: a list not ascending, a wall outside the box, a gas cell outside 1..IBAR, a face-write conflict with the same sign (check of section 4 of
  the wall-table spec), a species count larger than the callee's table.
* Device: the same random state hashed on host and on the GPU in the style of `gpu/s5_dev_hash.F90` (input, reference, generated hashes; bitwise equal, timings), at three sizes.
