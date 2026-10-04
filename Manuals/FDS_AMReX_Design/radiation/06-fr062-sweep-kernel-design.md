# FR-062 per-box sweep kernel design (L1242, RADIATION_FVM)

Owner: AMR Radiation Lead. Status: design for review; every number marked (est.) is an estimate, not a measurement.

Inputs: `docs/adr/drafts/rulings-IR008-FR062.md` §2.1, `01-radiation-amr-spec.md` §2 and §4, `04-radiation-phase4-notes.md`, `05-br-ilw-table-spec.md`, `loop-work-list.md`, `generator-howto.md`, and the sweep text in `Source/radi.f90` at `bee11f0329`. All line numbers below are that file.

## Summary

1. FDS already works the way FR-062 asks: each mesh sweeps from `IL_R` (filled at the last exchange) and writes `IL_S`; `MESH_EXCHANGE(2)` moves the buffers between passes. The per-box design formalises that and uses it for boxes of one level. A box-split run is therefore, by construction, the same arithmetic as a multi-mesh FDS run with the same meshes. This gives a stronger test than a tolerance check (§e).
2. Inside a pass, boxes are independent and, for the 3D sweep, all cells of one diagonal plane inside a box are independent. Angles are **not** independent in their bookkeeping: several quantities are accumulated in processing order over the angles (§b). Floating-point sums over angles must keep that order.
3. Recommendation: stage 1 is **plane-parallel across all boxes, serial over angles**, as a hand-written K2 kernel whose arithmetic is a verbatim copy of the cell body, not a generated one. Stage 2 (angle batching with per-angle slabs and an ordered reduction) only if measurement shows launch latency dominates. The generator can take the cell body later when it supports loops over a cell list with `CYCLE` (§e).
4. No atomics, no floating-point reductions across threads inside a pass, so thread count independence (FR-005 (i)) and run-to-run reproducibility (FR-005 (iv)) hold. The exemption is for box-split dependence only (D-039).

---

## (a) Data dependence of the sweeps

Common frame. Per band, per update (`UPDATE_LOOP`, 4287-), the angle loop is `DO N = NSTART,1,-ANGLE_INCREMENT` with `NSTART = NRA - ANGLE_INC_COUNTER + 1` (4315-4317): each update sweeps the angles of one residue class, in **descending** `N`. The default `ANGLE_INCREMENT` is `MIN(5, NRA/15)` (read.f90:10301), so about 20 angles of 100 per update; the initial update and `UPDATE_ALL_ANGLES` run all `ANGLE_INCREMENT` classes.

Per angle the direction signs `ISTEP,JSTEP,KSTEP = sign(DLX,DLY,DLZ)` pick the **octant**: eight sign cases. The sweep starts at the corner `(ISTART,JSTART,KSTART)` of the box and moves downwind (4388-4424).

**Ghost cells.** Before the sweep, `WALL_LOOP1` (4328-4371) writes `IL` in the ghost cell `(II,JJ,KK)` of every wall whose `DLN(IOR,N) >= 0` (the radiation enters the box through it): open wall, mirror, interpolated (from `IL_R`) or solid. Those are the only ghost values the sweep reads. `IL(:,:,:)` is set to the constant `BBFA*RPI_SIGMA*TMPA4` once per update (4319), not per angle.

### 3D Cartesian (4494-4718)

- Recurrence: `IL(I,J,K)` depends on `IL(I-ISTEP,J,K)`, `IL(I,J-JSTEP,K)`, `IL(I,J,K-KSTEP)` (the three upwind neighbours), the cell data (`EXTCOE`, `RTE_SOURCE`, `DX,DY,DZ`), `RSA(N)`, and the `CELL_ILW` overrides. Update: `IL = MAX(0, RAP*(AIU_SUM + VC*RSA(N)*RFPI*RTE_SOURCE))` with `AIU_SUM = AXU*ILXU + AYU*ILYU + AZU*ILZU (+AILFU)`, `A_SUM = AXD+AYD+AZD (+AFD)`, `RAP = 1/(A_SUM + EXTCOE*VC*RSA(N))`.
- Wavefront: `N_SLICE = ISTEP*I + JSTEP*J + KSTEP*K` runs from the value at the start corner to the value at the far corner (4496-4497). For a box of `n^3` that is `3n-2` planes. For each plane, `IJK_SLICE(:, 1:M_IJK)` is filled **serially** (4498-4515) with the `(I,J,K)` of the plane, then `SLICE_LOOP` runs over it in an OpenMP parallel loop (4517-4715). Cells of one plane do not depend on each other: each reads only cells of plane `N_SLICE-1`. This is why the CPU result does not depend on thread count.
- Solid cells: `IF (CELL(IC)%SOLID) CYCLE` before any write (4544): `IL` of a solid cell is never written in the update. Its value stays at the constant of 4319. The cell data of a gas cell next to a solid wall carry `CELL_ILW(IC,1:3)` overrides: if `CELL_ILW(IC,d) > -1e6` the upwind value in direction `d` is replaced by the wall intensity (4545-4547). `CELL_ILW` is reset to `-HUGE` at the start of **every angle** (4323) and filled by `WALL_LOOP1` for solid walls (4365).
- `RAD_DIFF_SCHEME > 1` (diamond or exponential): the upwind values come from `ILDX/ILDY/ILDZ` (4536-4540), which are three more fields holding the downwind face values (set at 4675-4709; initialised as copies of `IL` at the start of each angle, 4425-4429). The recurrence is then on three arrays; a negative-intensity correction loop (`NEG_ITER=1..3`) runs per cell. The three arrays are per-angle state; the dependence structure is the same plane-to-plane.
- `IL_UP(I,J,K)` (`SOLID_PARTICLES`) is written in the same cell update, `MAX(0, AIU_SUM/A_SUM)`; it is read after the sweep by the particle loops of the same angle.
- Cut cells (`CC_IBM`): extra face terms (4560-4617). Geometry-deferred; not designed here.

### 2D Cartesian (4470-4492)

`J=1`. Loops `K=KSTART,KEND,KSTEP`, `I=ISTART,IEND,ISTEP`. `IL(I,K)` depends on `IL(I-ISTEP,K)` and `IL(I,K-KSTEP)`. The natural wavefront is the anti-diagonal `ISTEP*I+KSTEP*K = const`: `nx+nz-1` planes. The same solid and `CELL_ILW` rules; no `ILD*` arrays in this branch. The two y-direction walls are skipped (`.NOT.TWO_D .OR. ABS(IOR)/=2` at 4339, `TWO_D .AND. .NOT.CYLINDRICAL .AND. ABS(IOR)==2` at 4756).

### Cylindrical (4431-4468)

`J=1`; three terms: `ILXU = IL(I-ISTEP,J,K)`, `ILYU = IL(I,J-JSTEP,K)`, `ILZU = IL(I,J,K-KSTEP)`. **`ILYU` is a ghost value**: `J-JSTEP` is `0` or `2`, outside the one-cell-wide box. It is set by `WALL_LOOP1`'s cylindrical branch (`IL(II,JJ,KK) = ILW(N)`, 4367-4369) from `ILW(N)`, which the previous angle's post-sweep copy wrote (`BR_UP%ILW(N-1) = BR_DOWN%ILW(N)` or `ILW(N)`, 4722-4742). The comment at 4312 says it: "in cylindrical case the Nth angle boundary condition comes from the (N+1)th angle". So in the cylindrical case **angle N depends on the result of angle N+1** (processed just before it): a chain of angles within each azimuthal row, `NRP(1)` long [VERIFY the exact row structure in `INIT_RADIATION`]. Angle-level parallelism is limited to different rows. `WEIGH_CYL=2` multiplies `DLN` and the `UIID` weight.

### Couplings between angles (what stops independent angle launches)

| Quantity | Where | Coupling |
|---|---|---|
| `UIID(:,:,:,IBND)` or `UIID(:,:,:,ANGLE_INC_COUNTER)` | 4818, 4820 | `UIID = UIID + WEIGH_CYL*RSA(N)*IL` for each `N` of the update, in descending `N`, starting from the zero set at the update start (4298-4303). A left-to-right chain over angles per cell. Whole-array: ghost cells and solid cells included. |
| `INRAD_W(IW)`, `INRAD_F` | 4761, 4774 (walls), 4799, 4801 (CFACEs) | Per wall: `+DLN*ILW_old` then `-DLN*ILW_new` for each incoming angle, in `N` order. A chain over angles per wall. Read at the next angle's `WALL_LOOP1` (solid walls) |
| `ILW(ANGLE_INC_COUNTER)` of open walls and CFACEs | zeroed 4305-4310, accumulated 4779-4787, 4803-4812 | `ILW(c) -= DLN*IL(IIG,JJG,KKG)` over all `N` of the update. A chain over angles per wall. |
| `ILW(N)` mirror | 4351 | `ILW(N) = ILW(DLM(N,\|IOR\|))`: reads **another angle** of the same wall. `DLM` is the nearest reflected direction; it need not be in the same residue class, and when it is, it is read before or after its own update depending on the order. |
| Cylindrical Y-ghost and `BR_UP/BR_DOWN` copy | 4367-4369, 4722-4742 | angle `N` reads what angle `N+1` wrote (above). |
| `IL` working array | 4319 | One array for all angles, reset once per update. Gas cells are overwritten every angle, solid cells never. **Ghost cells that no wall sets for the current angle keep the previous angle's value**, and `UIID` sums over ghost cells too. [VERIFY] whether anything reads `UII` in ghost cells; if not, the ghost values of `UII` need not be bit-reproduced (§e question 2). |
| `CELL_ILW` | reset per angle (4323) | per-angle state (small: only cells next to solid walls) |
| `IL_UP` | set in the cell update (e.g. 4628) | per-angle state, consumed by the particle loops of that angle |
| `IL_S(LL,N,IBND)` | 4831-4840 | written at index `N`: **disjoint** between angles. No coupling. |
| `ILW(N)` of solid walls and particles | 4363, 4873, 4878 (`ILW(N)`) | written at index `N`: disjoint between angles (the read of `INRAD_W` is the coupling). |

Between **bands**: all of the above are indexed by band or reset per band, so bands are independent, except that `QR` and `UII` are formed from all bands after the loop.

Conclusion: angles are independent as arithmetic on `IL`, but not as bookkeeping. Any scheme that runs angles concurrently has to keep the per-angle `IL` (and `CELL_ILW`, `IL_UP`, `ILD*`) separate and then apply the chains above in the original order.

---

## (b) Per-box kernel decomposition and summation order

Options, with what each costs.

| Option | Parallelism per launch | Memory | Chains over angles | Verdict |
|---|---|---|---|---|
| **A. Plane-parallel, serial over angles, all boxes in one launch** | cells of plane `d` in every box: `nbox * ~89` (16^3) or `nbox * ~349` (32^3) threads (est., §d) | one `IL` (plus `ILD*`) per box | exactly as the original: after each angle the boundary updates and the `UIID` update run in order | **Stage 1.** Bitwise equal to the CPU order by construction. |
| B. Plane-parallel, angles batched (`B` angles, per-angle slabs of `IL`, `CELL_ILW`, `IL_UP`, `ILD*`), ordered reduction after the batch | A times `B` | `B * cells * 8` bytes per array: 170 MB per 1M cells at `B=21` (est.; matches §4 of the spec) | one reduction kernel per chain, per cell or wall running over the `B` angles in descending `N`: the same left-to-right chain as the original | Stage 2, only if A is launch-bound. Bitwise equal if every chain keeps its order and mirror and cylindrical dependences are kept (see below). |
| C. One thread-block per box, the block runs the whole plane loop | one block per box | tiny | as A | Not recommended: a block has 1 to 1024 threads, one plane of a 16^3 box has up to 192 cells, 32^3 up to 768; fine, but the occupancy is then `nbox` blocks only, and the 46 or 94 planes need `__syncthreads` between planes. Acceptable only when `nbox` is in the thousands. Not portable to K2 `omp target teams` without a team-level barrier. |
| D. Order-free (exact fixed-point) `UIID` and `INRAD` sums | any | as B | order-independent | Changes the bits against the CPU. Allowed only under the D-053 switch (default off), as a performance option; not needed. |

**Rules for B (if built):**
1. Chains: `UIID`, `INRAD_W/F`, `ILW(ANGLE_INC_COUNTER)` are each reduced over the batch in the original order (descending `N`), one thread per cell or wall, serial over the batch. No atomics.
2. Mirror boundary: batches must be split so that an angle and its `DLM` partner are not in the same batch unless the order of their update is reproduced; simplest: run mirror cases with `B=1`.
3. Cylindrical: `B=1` within a row; angle batching only across independent rows [VERIFY].
4. Particle loops after the sweep (4848-4881, which read `IL`/`IL_UP`) run once per angle in the original order, from the slab of that angle.

**Atomics policy.** None. The only places FDS uses a parallel reduction in this routine are `RAD_Q_SUM` and `KFST4_SUM` (`!$OMP CRITICAL`, 4120-4123), which are per-step sums under FR-005 (ii) (FDS order by default, fixed-point sum as the switch) and are outside the sweep. The sweep and its boundary updates write disjoint elements (unique slots, unique cells), and `UIID` etc. are chains. Under (i) the radiation stage is exempt only for the box split (D-039); thread-count independence and (iv) run-to-run reproducibility still apply, and A satisfies them.

---

## (c) Face double-buffering and exchange (FR-062 rules (1) to (5))

Mapping to the FDS structures:

| FR-062 element | FDS today | AMR box |
|---|---|---|
| incoming face buffer (read by the sweep) | `OMESH(NOM)%IL_R(NIC_R, NRA, NSB)`, read in `WALL_LOOP1` (4354-4359) | `FACE_IN(box)`: one entry per ghost cell of the box that borders another box, per angle and band |
| outgoing face buffer (written by the sweep) | `OMESH(NOM)%IL_S(NIC_S, NRA, NSB)`, written after each angle (4818-4833) | `FACE_OUT(box)` |
| swap | `MESH_EXCHANGE(2)`: same rank: array copy; other rank: `MPI_ISEND` of packed `IL_S` for the angle subset (`ALLOCATE_RADIATION_SEND/RECV_PKG`, main.f90:4053-4132) | one mechanism for both |
| K passes | `DO ITER=1,RADIATION_ITERATIONS`, exchange after each pass when K>1 (main.f90:1095-1117), at the end of the step when K=1 (main.f90:1183-1189) | same |
| value in a ghost cell from a neighbour of a different size | `IL(II,JJ,KK) = SUM(IL_R(NIC_MIN:NIC_MAX))/(NIC_MAX-NIC_MIN+1)` (4355-4359) | the same average gives fine-to-coarse restriction; a fine ghost reading a coarse cell is the `NIC=1` case (constant prolongation) |

Design points:

1. **Double buffering is not new code.** The sweep writes only `FACE_OUT`; it reads only `FACE_IN`. A box can never see a value from the pass in progress. Same-rank neighbours go through the buffers too (no direct read of the neighbour's `IL`). This gives rule (1): results do not depend on rank, ownership or processing order.
2. **Exchange only the angles that were swept.** The send buffer covers the angle subset of the update (`ANGLE_INCREMENT` classes), as `ALLOCATE_RADIATION_SEND_PKG` does; face values for other angles stay as left by the earlier exchange.
3. **K passes.** Pass 1 reads the values left by the last update; pass `k` reads those of pass `k-1`. `INITIAL_RADIATION_ITERATIONS=3` for the initial fill. The exchange count per pass is one (the `IF (ICYC>1) EXIT` at main.f90:1114 means one exchange per pass after the first cycle).
4. **After regrid.** New faces are filled from the upwind interior or from the coarse level (rule 3, `04` §1.4), never zero. `UPDATE_ALL_ANGLES` for the next update.
5. **Determinism.** Packing and unpacking are index copies; no sums. The average over `NIC` cells is a fixed-order sum (ascending `LL`), part of the ghost fill, kept in the same order on all ranks.
6. **Escalation.** The larger-box radiation BoxArray (rule 5, first step) is the same design with different boxes; the sweep kernel does not change. An ordered sweep in global box order would change it (a box must wait for its upwind boxes), so it is a separate scheduling layer and is not designed here.

Memory for the face buffers (est., 8-byte reals): a 16^3 box has up to `6*256 = 1536` boundary ghost cells; at `NRA=100`, `NSB=1`: 1.2 MB per buffer, three buffers (`IL_S`, `IL_R`, `IL_R_OLD`, main.f90:2392-2400) 3.7 MB per box per band. Storing only the outgoing half of the angles would halve it; FDS does not.

---

## (d) Memory and occupancy (est.)

Box of `n^3` cells, `3n-2` planes, from a direct count:

| | 16^3 | 32^3 |
|---|---|---|
| planes per angle | 46 | 94 |
| average cells per plane | 89 | 349 |
| largest plane | 192 | 768 |
| `IL` with one ghost layer, 8 bytes | 47 kB | 314 kB |
| `IL` + `ILDX/Y/Z` (RAD_DIFF_SCHEME>1) | 190 kB | 1.3 MB |
| angle-batched, `B=21` | 1.0 MB | 6.6 MB |

Occupancy. A GPU needs of the order of 10^4 to 10^5 resident threads to hide latency (rule of thumb, est.). One box in one plane gives 89 or 349 threads, so **one-box launches are too small by two to three orders of magnitude**. What fills the device, in order of preference:
1. **All boxes of the level (and of the rank) in one launch per plane** (option A): `nbox * 89`. 1000 boxes of 16^3 (4M cells) give 89k threads (est.); 100 boxes of 32^3 give 35k. This needs plane `d` of every box in one launch, so the plane index must be relative to each box's start corner: `d = |i-i0|+|j-j0|+|k-k0|`, which makes the plane count the same for all boxes of the same size and octant (the cell list of the original is replaced by arithmetic on `d`: thread to `(i,j,k)` with the octant flips).
2. **Angles of different octants together** (option B): the planes `d` of all angles in flight run together, so `nbox * B * 89` threads.
3. **Bands** (independent): multiplies by `NSB`; memory scales with it, so only for wide-band runs and with care.
4. Boxes of different sizes in one launch: pad to the maximum plane length with a bounds test.

Launch count per band per update, option A: `46` planes (16^3) per angle, about 21 angles, so about 970 plane launches plus about 6 small boundary and accumulation launches per angle (est.: 1100 launches). At 5 to 10 microseconds each this is 5 to 11 ms (est.). The arithmetic is small: 4M cells * 21 angles * ~50 flops = 4 GFlop, a few ms (est.). **Stage 1 is probably launch-bound by a factor of two to four on a large device; not memory- or flop-bound** (est., to be measured). Options: B (fewer launches by a factor `B`), CUDA graphs (the launch sequence is identical every update), or fusing the boundary updates into the plane kernel. These are measurement-driven decisions.

When boxes are few (a refined level of a few 16^3 boxes on one GPU), no scheme fills the GPU; the radiation stage stays latency-bound and the answer is to put several levels or several bands in flight, or accept it, since radiation updates happen once per `TIME_STEP_INCREMENT` steps and only for one residue class (about 1/5 of the angles).

---

## (e) What the generator can do, and the plan

**Now (verified in the worktree, `03` §5).** The generator refuses the 2D and cylindrical sweeps ("loop with a step") and has no support for: a variable step; the loop-carried recurrence `IL(I-ISTEP,...)` with `CYCLE`; the cell-list wavefront (`IJK_SLICE`, filled serially); `PRIVATE` locals such as `NEGATIVE_MASK(3)` and `PLX(3)` with `SELECT CASE`; the `CC_IBM` branches (blanked in my entries only).

**Recommendation.** For the sweep cell update: **hand-written K2 kernel** in a radiation-owned file, with the arithmetic copied **verbatim** from the original text by line range (the arithmetic, not a rewrite), because:
- the data dependence (plane-to-plane) is scheduling, which is not what the generator models;
- the cell body is long and contains the exact operation order that makes the result bitwise (`AIU_SUM`, `RAP`, `MAX`), so a text copy with a drift check against `radi.f90` is safer than a re-expression;
- adding the generator features (loop over a list with `CYCLE`, variable step, local arrays with `SELECT CASE`) is worthwhile only if the Generator Engineer wants them for other loops. If they do, the plane body (`SLICE_LOOP` body) becomes a generated per-cell kernel over a cell list, and the plane loop and the octant arithmetic stay hand-written. The test below is the same either way.

**Staged plan with acceptance tests.**

| Stage | Content | Bitwise acceptance |
|---|---|---|
| 0 (done) | Source and boundary pieces generated: `RTE_SOURCE`, `EXTCOE`, `UIIOLD`, `UIID`, `UII`, `QR`, `KFST4` branches, L1239 | `test_mesh_rad.py` (six flag sets) |
| 1 | `BR_ILW` slot table and the sweep boundary kernels (`05`) | `open_qin_spec` generated; WALL_LOOP1/2/3 against the verbatim loops for one angle and for a full angle sequence |
| 2 | 3D plane kernel, one box, one angle, `RAD_DIFF_SCHEME=1`, `CC_IBM` off | `IL` array byte-identical to the verbatim original (`SLICE_LOOP` run on the CPU) for all eight octants, solid cells, `CELL_ILW` overrides, ghost values from `WALL_LOOP1`, random data, at 1/4/8 threads, six flag sets |
| 3 | One box, full angle sequence with boundary kernels and `UIID` chain | `IL`, `UIID`, `ILW`, `INRAD_W`, `IL_S` byte-identical to the original for the update; cases: gray, wide band (6 bands), `UPDATE_ALL_ANGLES` |
| 4 | Box-batched launches (option A), several boxes of different sizes | same data as stage 3 per box, boxes interleaved; result independent of box order and of the batch composition |
| 5 | Face buffers and exchange, K passes | **Box split against multi-mesh FDS**: the same domain run as 2 (and 8) FDS meshes and as 2 (and 8) boxes of the AMR radiation module, same `K`: `IL`, `UIID`, `QR` byte-identical cell for cell. This is expected because the design reproduces FDS's exchange (§c.1). The FR-062 (4) lag check against single-mesh (`radiation_gas_panel`, 16^3 and 32^3, K=1,2,3) then measures the physical lag and sets the tolerance. Also byte-identical across 1, 2, 4 ranks |
| 6 | 2D (anti-diagonal planes) | as stage 2 for `TWO_D` |
| 7 | `RAD_DIFF_SCHEME` 2 and 3 | as stage 2 with `ILD*` compared as well |
| 8 | Angle batching (option B) | identical to stage 3 for every batch size; mirror and cylindrical cases forced to `B=1` and checked |
| 9 | Cylindrical | as stage 3; angle chain kept |
| later | `CC_IBM`, particles on the device | with the geometry refactor and the deposition kernel |

Mutants for every stage: a different upwind index, a swapped octant sign, `MAX` dropped, `UIID` accumulation in ascending `N`, solid cells written, `CELL_ILW` ignored, face buffer read in place. Each must fail.

Reference build: the original Fortran compiled with gfortran, `-ffp-contract=off`, as in the existing radiation tests. FMA contraction must also be off in the device kernel for the bitwise claim; if the device compiler cannot guarantee that, the claim becomes "within rounding" and the acceptance test becomes a tolerance on `IL` (that choice belongs to V&V).

---

## (f) Open questions

For the Chief Architect:
1. Is plane-parallel across boxes, serial over angles (option A) accepted as the baseline, with angle batching deferred until measured? It keeps all chains in their original order.
2. Ghost-cell `UII`/`UIID`: the original sums `IL` over ghost cells, where `IL` can be stale from the previous angle. If nothing reads those values, may the AMR port define them as "not bit-reproduced" (or zero)? V&V needs to confirm that nothing reads ghost `UII` (output or wall code).
3. The comparison with multi-mesh FDS in stage 5: is that an accepted T0 test, given that FR-005 exempts the box-split dependence? It tests our kernels, not the exemption.
4. Is a hand-written K2 kernel acceptable for this loop (the generator rules say generated where possible), given the reasons in §e?
5. Cylindrical geometry: in scope for the first release, or deferred? It is the only branch with an angle-to-angle chain.

For the Legacy Mapper:
6. Should the sweep be registered as a loop claim (L1242) with a "hand-written, FR-062" tag, or stay outside the register as the sweep and its pieces? My claim covers the sub-nests only.
7. Do the box-level launch batching and the plane arithmetic belong to the driver layer owned by the AMReX Integration Lead? I assume yes.

For the Generator Engineer:
8. Is a loop over a cell list with `CYCLE` and a bounded `PRIVATE` set (the `SLICE_LOOP` body) on your list? If yes, the cell body can be generated and the hand-written part shrinks to the plane scheduling.
