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
| swap | `MESH_EXCHANGE(2)`: same rank: array copy; other rank: `MPI_ISEND` of packed `IL_S` for the angle subset (`ALLOCATE_RADIATION_SEND/RECV_PKG`, main.f90:4087-4166 at FireX `8ad685f435`; host; AMR -> FR-062 face exchange, §18) | one mechanism for both |
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

---

## Addendum 1 (D-073): rulings, kernel descriptor, angle loop

Status: stages 1 to 3 in progress (kernels written, drift check and bitwise test in `amrex/s4_mass/s5_gen/s5_rad_sweep/`, branch `s5-gen`). The kernel file is `s5rad_sweep_k2.F90` (module `s5rad_sweep_k2`); the host helpers are in `s5rad_sweep_host.F90`. Every signature below is the one in that file; where this text and the file differ, the file is right and this text is to be corrected.

### 1. Rulings recorded (answers to §f)

| § f item | D-073 ruling | Effect on this design |
|---|---|---|
| 1 option A baseline | Accepted: plane-parallel across boxes, serial over angles. Angle batching only after measurement shows launch latency dominates, with the ordered-chain rules for option B; `B=1` for mirror (and cylindrical). | Stage 1 is option A. Stage 8 (batching) stays closed until measured. |
| 2 ghost-cell `UII`/`UIID` | May be "not bit-reproduced" provided the V&V Lead confirms that nothing reads them (outputs, wall code, devices). The check is listed in the stage-3 test. | `test/sweep_ghost_reads.py` lists every reader (§6). The single-box kernel path does reproduce the ghost cells; the box-split path may not. |
| 3 T0 test | Accepted: box split against multi-mesh FDS, byte-identical `IL`/`UIID`/`QR` on 2 and 8 boxes and 1/2/4 ranks. It does not test the FR-005 exemption. | Stage 5 as written in §e. |
| 4 hand-written K2 kernel | Accepted with: cell-body arithmetic copied verbatim by line range; drift check against `radi.f90` in CI; K2 CI check and kernel lint passing; registry entry tagged hand-written/FR-062; no FMA contraction (`-gpu=nofma`). Regenerate the cell body later if the generator supports cell-list loops with `CYCLE`. | Done in the file (marks `!S5RAD_COPY radi.f90:a-b`), `check_sweep_drift.py`, `s5rad_sweep.toml`. |
| 5 cylindrical | Out of the first release. Cylindrical meshes stay FDS-only in AMR mode. | Stage 9 and the cylindrical rows of §a and §b are closed. The host refuses `CYLINDRICAL=T`. The kernels contain no cylindrical term. |
| 6 register | The sweep is registered as claim L1242, tagged hand-written (the Legacy Mapper owns the register). | Registry data in `s5rad_sweep.toml` (`[meta]`: claim, tag, FR-062, D-073). |
| 7 box-level batching | Belongs to the Integration Lead's driver layer. Per-box plane index ranges (box id, plane index, lo/hi in the two in-plane directions) as a small device array uploaded once per regrid; K2 launches are blocking and host-synchronous, no stream. | §2 and §3 are written to that. |

### 2. Kernel descriptor

**Launch unit.** One blocking launch per (band, update, angle `N`, plane `D`). The plane kernels' thread space is `(box, k', j')` in 3D and `(box, k')` in 2D, where `(i',j',k')` are corner-relative cell offsets (see below). The host passes the widths `KW`, `JW` of the thread space for plane `D` (the maximum of the in-plane extents over all boxes, from `s5rad_plane_widths3d` / `s5rad_plane_widths2d`); a thread outside its box's rectangle, or with `D >= NPL(box)`, returns without writing.

**What the plane kernel indexes by**

| Index | Where it comes from | Meaning |
|---|---|---|
| box id `IB`, 1..`NB` | thread index; last (trailing) dimension of every field array | the box. `BOX_IBAR(IB)`, `BOX_JBAR(IB)`, `BOX_KBAR(IB)` give its size. Boxes smaller than `NXM,NYM,NZM` are padded. Neither the plane kernels nor the cell kernels (`s5rad_il_fill`, `s5rad_uiid_zero`, `s5rad_uiid_acc`, `s5rad_uii_sum`, `s5rad_qr_*`) touch the padding: each cell kernel loops over the padded nest and skips a cell with `I > BOX_IBAR(IB)+1 .OR. J > BOX_JBAR(IB)+1 .OR. K > BOX_KBAR(IB)+1`. The three arrays are int32 `(NB)` (in 2D `BOX_JBAR = 1`, `NYM = 1`). |
| plane index `D`, 0..`NPL(IB)-1` | scalar argument | `D = i'+j'+k'`. `NPL(IB) = IBAR+JBAR+KBAR-2` (2D, `J=1`: `IBAR+KBAR-1`). The host launches `D = 0,1,..,max(NPL)-1` in ascending order, one blocking launch each. |
| in-plane lo/hi per box | `PLN(1:4,D,IB) = (k' lo, k' hi, j' lo, j' hi)` int32, `PLN(4,0:NPLMAX-1,NB)`; 2D: `PLN2(1:2,D,IB) = (k' lo, k' hi)` | bounding rectangle of plane `D` in the two in-plane directions (octant independent). `k' lo = max(0, D-(IBAR-1)-(JBAR-1))`, `k' hi = min(KBAR-1, D)`, `j' lo = max(0, D-(IBAR-1)-(KBAR-1))`, `j' hi = min(JBAR-1, D)`. The third offset follows: `i' = D-j'-k'`; the thread skips the cell if `i'` is outside `0..IBAR-1`. The table is built by `s5rad_plane_table3d/2d` and uploaded once per regrid (together with `NPL`, `BOX_*BAR`). |
| octant signs `ISTEP,JSTEP,KSTEP` | scalar arguments, each +1 or -1 | as upstream 4386-4429: `-1` where `DLX(N)<0`, `DLY(N)<0`, `DLZ(N)<0`. They map the offsets to the cell: `I = 1+i'` if `ISTEP>0`, else `IBAR-i'` (same for `J`, `K`). The upwind cells are `(I-ISTEP,J,K)`, `(I,J-JSTEP,K)`, `(I,J,K-KSTEP)`. `s5rad_octant` computes them on the host. |
| angle index `N`, 1..`NRA` | scalar argument | selects `DLX(N)`, `DLY(N)`, `DLZ(N)`, `RSA(N)` (device arrays of length `NRA`, read as in the original: `ABS(DLX(N))`, `RSA(N)`); `DLN(-3:3,N)` and `DLM(N,1:3)` in the wall kernels. |
| band `IBND` | scalar to the wall kernels, `s5rad_ils`, `s5rad_open_ilw_zero`, `s5rad_uiid_acc` (as the `ISL` slice) | the plane kernels do not take the band: the host passes the band's `EXTCOE` and `RTE_SOURCE` (arrays `(0:NXM+1,0:NYM+1,0:NZM+1,NB)` of that band). The band indexes `BR_ILW(:,:,IBND)`, `ILR(:,:,IBND)`, `ILS(:,:,IBND)`. |

**Signatures** (all `bind(c)`; `integer(c_int)` and `real(c_double)`; scalars by value; arrays explicit-shape):

```
s5rad_plane3d(NB,NXM,NYM,NZM,NPLMAX,KW,JW,D,ISTEP,JSTEP,KSTEP,SOLID_PARTICLES,NCELL,NRA,N,RFPI,
              BOX_IBAR,BOX_JBAR,BOX_KBAR,NPL,PLN,DX,DY,DZ,EXTCOE,RTE_SOURCE,CELL_INDEX,CELL_SOLID,CELL_ILW,
              DLX,DLY,DLZ,RSA,IL,IL_UP)
s5rad_plane2d(NB,NXM,NZM,NPLMAX,KW,D,ISTEP,KSTEP,SOLID_PARTICLES,NCELL,NRA,N,RFPI,
              BOX_IBAR,BOX_KBAR,NPL2,PLN2,DX,DZ,EXTCOE,RTE_SOURCE,CELL_INDEX,CELL_SOLID,CELL_ILW,
              DLX,DLZ,RSA,IL,IL_UP)                     ! J = 1; the y extent of the arrays is 0:2
s5rad_il_fill(NB,NXM,NYM,NZM,BBFA,RPI_SIGMA,TMPA4,BOX_IBAR,BOX_JBAR,BOX_KBAR,IL)
s5rad_cell_ilw_reset(NCELL,CELL_ILW)
s5rad_uiid_zero(NB,NXM,NYM,NZM,NS,ISL,BOX_IBAR,BOX_JBAR,BOX_KBAR,UIID)
s5rad_open_ilw_zero(NWE,NBR,NRA,NSB,IBND,ANGLE_INC_COUNTER,W_BTYPE,W_BR_INDEX,BR_ILW)
s5rad_wall1(NB,NXM,NYM,NZM,NW,NCELL,NBR,NRA,NSB,NICR,N,IBND,TWO_D,BBFA,RPI_SIGMA,TMPA4,RPI,
            W_BOX,W_BTYPE,W_BR_INDEX,W_IOR,W_II,W_JJ,W_KK,W_IIG,W_JJG,W_KKG,W_EMISSIVITY,W_TMP_EXTERIOR,W_RAMP_VAL,
            W_NIC_MIN,W_NIC_MAX,DLN,DLM,TMP_0,OUTRAD_W,INRAD_W,ILR,CELL_INDEX,IL,CELL_ILW,BR_ILW)
s5rad_wall2(NB,NXM,NYM,NZM,NW,NBR,NRA,NSB,N,IBND,TWO_D,W_BOX,W_BTYPE,W_BR_INDEX,W_IOR,W_IIG,W_JJG,W_KKG,DLN,IL,INRAD_W,BR_ILW)
s5rad_wall3(NB,NXM,NYM,NZM,NW,NBR,NRA,NSB,N,IBND,ANGLE_INC_COUNTER,W_BOX,W_BTYPE,W_BR_INDEX,W_IOR,W_IIG,W_JJG,W_KKG,DLN,IL,BR_ILW)
s5rad_uiid_acc(NB,NXM,NYM,NZM,NS,NRA,N,ISL,WEIGH_CYL,BOX_IBAR,BOX_JBAR,BOX_KBAR,RSA,IL,UIID)
s5rad_ils(NB,NXM,NYM,NZM,NICS,NRA,NSB,N,IBND,ILS_BOX,IIO_S,JJO_S,KKO_S,IL,ILS)
s5rad_uii_sum(NB,NXM,NYM,NZM,NS,BOX_IBAR,BOX_JBAR,BOX_KBAR,UIID,UII)                                       ! radi.f90:4963
s5rad_qr_wb(NB,NXM,NYM,NZM,NS,ISL,BOX_IBAR,BOX_JBAR,BOX_KBAR,KAPPA_GAS,KFST4_GAS,UIID,QR)                  ! radi.f90:4949
s5rad_qr_w_wb(NB,NXM,NYM,NZM,NS,ISL,BOX_IBAR,BOX_JBAR,BOX_KBAR,KAPPA_PART,KFST4_PART,UIID,QR_W)            ! radi.f90:4954
s5rad_qr_gray(NB,NXM,NYM,NZM,BOX_IBAR,BOX_JBAR,BOX_KBAR,KAPPA_GAS,KFST4_GAS,UII,QR)                        ! radi.f90:4982
s5rad_qr_w_gray(NB,NXM,NYM,NZM,BOX_IBAR,BOX_JBAR,BOX_KBAR,KAPPA_PART,KFST4_PART,UII,QR_W)                  ! radi.f90:4983
```

`NXM,NYM,NZM` are the maximum box sizes (`NYM=1` in 2D). `ISL` is the `UIID` slice: `IBND` for wide band and WSGG, `ANGLE_INC_COUNTER` for gray. `TWO_D` is 0 or 1. The wall kernels run one thread per wall over all boxes (`W_BOX(IW)` is the box). `BOX_IBAR/JBAR/KBAR` are int32 `(NB)`. The `s5rad_qr_*` and `s5rad_uii_sum` kernels are the box-stacked forms of the statements 4949, 4954, 4963, 4982, 4983; the host keeps the upstream conditions: the band forms after the update of every band, `s5rad_uii_sum` when `UPDATE_INTENSITY`, `s5rad_qr_gray` only when not wide band and not WSGG, `s5rad_qr_w_wb` only when `N_LP_ARRAY_INDICES>0 .AND. PARTICLES_EXISTED`, `s5rad_qr_w_gray` only when `N_LP_ARRAY_INDICES>0`. `STORE_RADIATION_TERMS` (4950-4953, 4984-4987) is not covered.

**Array layouts** (Fortran column-major; box is always the last index of a field):

| Array | Shape | Notes |
|---|---|---|
| `IL`, `IL_UP`, `EXTCOE`, `RTE_SOURCE` | `(0:NXM+1, 0:NYM+1, 0:NZM+1, NB)` | one ghost layer per box at `0` and `BOX_*BAR+1`. `IL` is a working array, not an output. |
| `UIID` | `(0:NXM+1, 0:NYM+1, 0:NZM+1, NS, NB)` | `NS` = `UIIDIM` (number of bands for wide band and WSGG, `ANGLE_INCREMENT` for gray). Slice `ISL` is zeroed and accumulated; the other slices are not touched. |
| `UII`, `QR`, `QR_W`, `KAPPA_GAS`, `KFST4_GAS`, `KAPPA_PART`, `KFST4_PART` | `(0:NXM+1, 0:NYM+1, 0:NZM+1, NB)` | box last, like `IL`. `KAPPA_GAS`, `KFST4_GAS`, `KAPPA_PART`, `KFST4_PART` are the band's gas and particle absorption and emission arrays as in the upstream (`KAPPA_GAS` is not `EXTCOE`, which also holds the particle and scattering terms; the test sets them equal only for the gas array it feeds to the sweep). The host chooses the band (a section `(...,:,IBND)` of a `(...,NB,NSB)` array is contiguous). Made by `s5rad_uii_sum`, `s5rad_qr_wb`, `s5rad_qr_gray`, `s5rad_qr_w_wb`, `s5rad_qr_w_gray`; the older generated single-box kernels `rad_uii_sum`, `rad_qr_wb`, `rad_qr_gray` remain valid for one box. |
| `DX`, `DY`, `DZ`, `TMP_0` | `(0:NXM+1,NB)`, `(0:NYM+1,NB)`, `(0:NZM+1,NB)` | box-wise grid data |
| `CELL_INDEX` | `(0:NXM+1,0:NYM+1,0:NZM+1,NB)` | global cell number `IC` (0 = no record); `CELL_SOLID(IC)` is 0/1, `CELL_ILW(IC,1:3)` |
| wall tables | `W_*(NW)`, `NW = sum over boxes of NWE+NWI` | flat over all boxes; `W_BOX` gives the box; `W_II..W_KKG` are the ghost cell `(II,JJ,KK)` and the gas cell `(IIG,JJG,KKG)` of the wall; `W_RAMP_VAL` is `EVALUATE_RAMP(T-T_BEGIN,...)` done on the host once per radiation update (T and T_BEGIN do not change during the call, so the value is the same for every angle and band), for open walls with `TMP_EXTERIOR>0`; it is not evaluated per angle |
| `BR_ILW` | `(NBR, NRA, NSB)` | slot-fastest table of `05`; `W_BR_INDEX(IW)` is the slot (global over boxes) |
| `DLN`, `DLM` | `(-3:3,NRA)`, `(NRA,3)` | as `DLN(IOR,N)`, `DLM(N,\|IOR\|)` |
| `ILR`, `ILS` | `(NICR,NRA,NSB)`, `(NICS,NRA,NSB)` | the `IL_R` of all other meshes concatenated per box (wall positions `W_NIC_MIN/MAX`); `IL_S` of all boxes concatenated, `ILS_BOX(LL)` is the box of entry `LL`. `IIO_S, JJO_S, KKO_S(LL)` are indices in the SENDING box's own index space (`MESHES(NOM)%OMESH(NM)%IIO_S = MESHES(NM)%OMESH(NOM)%IIO_R`, main.f90:1584): the sender's own cells that coincide with the receiver's ghost cells (interior face cells at equal resolution), not cells outside the sender. The read is `IL(IIO_S(LL),JJO_S(LL),KKO_S(LL))` (radi.f90:4831). |

**Driver duty for `s5rad_open_ilw_zero`.** The upstream zeroing loop (radi.f90:4305-4310) runs over `N_EXTERNAL_WALL_CELLS` only. The kernel runs over the flat table, so it is equivalent only if no internal wall row is `OPEN_BOUNDARY`. That is not guaranteed (init.f90:5085 copies `OPEN_BOUNDARY` from the original surface of an external wall cell; nothing excludes an internal row). The driver must pass to this kernel a `W_BTYPE` copy in which the internal rows are not `OPEN_BOUNDARY` (for example `NULL_BOUNDARY`), and the true `W_BTYPE` to the other wall kernels. The bitwise test passes such a copy and contains internal OPEN rows; the mutant `open_zero_internal` (true `W_BTYPE` passed) must fail. Second consequence: an internal row that is `OPEN_BOUNDARY` writes the ghost cell `IL(II,JJ,KK)` in `s5rad_wall1` (the solid cell it lies in), and several internal rows can share one solid cell. Two active rows on one cell are written in table order by the original (the last wins) and race in the kernel; the host assertion below (a unique ghost cell for every active non-NULL row) must therefore cover the internal rows, or the host must refuse internal OPEN rows. Solid-type rows do not write `IL(II,JJ,KK)`. The test gives a second OPEN row on a solid cell the type NULL; with the duplicates left in, the 2, 3 and 4 thread runs of the test failed (found while adding the internal OPEN rows).

**Host duties and refusals.** The host refuses `CYLINDRICAL`, `CC_IBM`, `RAD_DIFF_SCHEME>1`, CFACE walls and the particle ORIENTATION loops in this release. It asserts, per table build: every non-NULL wall owns a unique ghost cell `(II,JJ,KK)` and a unique `BR_ILW` slot (no two threads write one element); `CELL_INDEX /= 0` for the gas cell of every solid-type wall; `D < NPL(IB)` ranges are those of `PLN`. Zero direction cosines together with solid walls on opposite faces of a one-cell-wide pocket would make two walls write one `CELL_ILW` element in the original as well; the host does not need to handle it but the test avoids it.

### 3. Angle-loop structure (per band, per update)

```
DO IBND = 1, NSB                       ! 1 only for gray
  DO I_UIID = 1, N_UPDATES             ! ANGLE_INCREMENT for UPDATE_ALL_ANGLES, else 1
    ANGLE_INC_COUNTER = MOD(ANGLE_INC_COUNTER,ANGLE_INCREMENT)+1            ! host
    s5rad_uiid_zero(ISL), s5rad_open_ilw_zero, s5rad_il_fill                ! 3 launches
    DO N = NRA-ANGLE_INC_COUNTER+1, 1, -ANGLE_INCREMENT                     ! descending, host loop
      s5rad_cell_ilw_reset
      s5rad_wall1                      ! ghost cells, CELL_ILW, mirror/open/interpolated/solid
      octant signs from DLX(N),DLY(N),DLZ(N)
      DO D = 0, max(NPL)-1             ! ascending, one blocking launch per plane
        s5rad_plane3d (or plane2d)
      s5rad_wall2, s5rad_wall3         ! INRAD_W, ILW chains in N order
      s5rad_uiid_acc                   ! UIID chain in descending N
      s5rad_ils                        ! IL_S(:,N,IBND), disjoint per angle
    ENDDO
  ENDDO
  (wide band / WSGG) s5rad_qr_wb, and s5rad_qr_w_wb if N_LP_ARRAY_INDICES>0 .AND. PARTICLES_EXISTED   ! after the band's update
ENDDO
IF (UPDATE_INTENSITY) s5rad_uii_sum
(gray) s5rad_qr_gray, and s5rad_qr_w_gray if N_LP_ARRAY_INDICES>0
```

Launches per angle: `max(NPL) + 6`; per update 3 more. Everything is serial over angles and planes and host-synchronous: no stream, no event. The chains over angles (`UIID`, `INRAD_W`, `ILW(ANGLE_INC_COUNTER)`, the mirror read of `ILW(DLM)`) therefore keep the original order by construction. Box-level launch batching (several boxes per launch) is already in the kernels (box id in the thread space); which boxes go into one launch, and when the tables are uploaded, is the driver's decision. The angle loop and the octant arithmetic stay on the host; they are not in the kernel file.

Arithmetic rules for the kernel file: FMA contraction off (gfortran `-ffp-contract=off`, nvfortran `-gpu=nofma`); no regrouping of the sums in the copied lines (`AIU_SUM`, `A_SUM`, the `RAP*(...)` product); `MAX` as written. The test mutants include FMA contraction and regrouped sums and must fail.

### 4. Stage 2/3 tests (what is checked)

`test/run_sweep.sh` builds the verbatim upstream sweep (radi.f90 lines 4287-4899 minus the ranges listed as omitted in `s5rad_sweep.toml`, cut by `make_sweep_tests.py`; nothing is retyped) and the kernel path, runs the same random state through both and compares bit for bit: `IL`, `IL_UP`, `UIID` (interior and full array), `UII`, `QR`, `QR_W`, `CELL_ILW`, `INRAD_W`, every `ILW` slot (all walls, angles, bands, including the open-wall accumulations), `IL_S`, the angle counter and the completed flag. Cases: 3D and 2D; gray, WSGG and wide band; all eight octants; solid cells and cells without record; open, mirror, interpolated, solid and NULL walls; zero direction cosines; `UPDATE_ALL_ANGLES`; stacked boxes of different sizes (including a group with one-cell and large boxes in one stack); internal OPEN wall rows; particle source terms present, absent and `PARTICLES_EXISTED` false; plane tables partition every box. The padding of every box smaller than `NXM,NYM,NZM` (cells beyond `0:IBAR+1` etc. of `IL`, `IL_UP`, `UIID`, `UII`, `QR`, `QR_W`) must still hold the value it was allocated with, bit for bit (the run prints `PADDING cells checked=... touched=0`). `UII`, `QR`, `QR_W` are made for all boxes in one launch by the stacked kernels and compared with the verbatim statements 4949, 4954, 4963, 4982, 4983. Flag sets as the other radiation tests, 1 and 4 threads. Mutants (`run_sweep_mutants.sh`): wrong octant sign, wrong plane order, ascending angle order, FMA contraction, regrouped sums, solid cell written, `CELL_ILW` ignored, wrong upwind cell, wrong mirror partner, a cell kernel that writes the padding (one mutant per kernel), a wrong or off-by-one extent guard, a descending `UII` sum, regrouped, re-signed or non-accumulating `QR` and `QR_W` forms, the `QR_W` host condition ignored, the true `W_BTYPE` passed to `s5rad_open_ilw_zero`; each must fail. The drift check has its own negative test (a changed post-sweep body, a removed guard, a wrong declaration of the extent arrays).

### 5. CI

`check_sweep_drift.py` (marked blocks equal to the upstream lines at the pinned commit through the documented renames; omitted ranges pinned by hash; every statement of each sweep region is copied, omitted or glue), `docs/tools/k2_ci_check.py` and `docs/tools/kernel_lint.py --clauses-only` on `s5rad_sweep_k2.F90`. The default K2 globs do not include the new directory; CI must pass the file explicitly (or the glob must be extended; an edit of the shared tools is not mine to make).

### 6. Ghost cells of UII/UIID/QR (D-073 (2) / D-077 (6)): reader list for V&V confirmation

`test/sweep_ghost_reads.py` runs `git grep` for `UII`, `UIID`, `QR` over `Source/*.f90` at the pinned commit and classifies every line; an unclassified line fails the check. Re-run at FireX pin `bee11f0329` (2026-10-06): **77** mentioning lines; classes `N=32`, `W=14`, `I=11`, `G=1`, `V=9`, comments `-=10`; **0 UNCLASSIFIED**. Exit `GHOSTREAD OK`.

**Already decided (not open):**

- Interior physics readers are interior or gas-cell: `divg.f90` / `radi.f90` `RAD_Q_SUM` / `KFST4` loops (`I`); wall `UII(BC%IIG,BC%JJG,BC%KKG)` (`G`, radi.f90:5035); whole-array `W` statements only feed other ghost elements inside `RADIATION_FVM`.
- **Restart-file comparison is interior-only** for `QR`/`UII`/`UIID` ghost cells (D-077 (6); Architect-accepted; `vv-runs/gpu_gate/README.md` § "Rule D-073 (2)"). A ghost-only difference is a NOTE, not a FAIL. **Implementation** of the mask waits on a restart-file comparison path (none exists yet).

**Nine `V` lines still need an explicit V&V Lead confirmation** so D-073 (b) / D-077 (6) / A-50 can close the "nothing reads ghost UII/UIID/QR" clause for the box-split path (which may leave those ghosts not bit-reproduced):

| # | Pin `bee11f0329` | Class of read | Proposed confirmation |
|---|---|---|---|
| 1–6 | `dump.f90:3918,3920,3921` WRITE and `:4107,4109,4110` READ of full `QR`/`UII`/`UIID` | restart I/O (ghost bytes in the file) | Covered by D-077 (6) interior-only restart compare. Confirm: no other restart consumer requires ghost bit-identity of these three arrays. |
| 7–8 | `dump.f90:8767` INTEGRATED INTENSITY `UII(II,JJ,KK)`; `:8769` RADIATION LOSS `QR(II,JJ,KK)` | device/slice output; unlike `CASE(16)` these do **not** clamp to `1..IBAR` | Confirm: `II,JJ,KK` for these quantities never address the ghost layer (device init / `GET_IJK` / slice index stay in `1..*BAR`, or out-of-mesh devices are refused). If a ghost address is possible, name the rule (clamp, refuse, or keep ghost bit-reproduced). |
| 9 | `dump.f90:9184` thermocouple `UII(II,JJ,KK)` | particle/thermocouple device cell | Same as 7–8 for the device cell of `CASE(110)`. |

Single-box kernel path: ghost cells of `UIID`/`UII`/`QR` stay bit-identical (sweep test `GHOST` line). Box-split path may differ there only after V&V signs the table above.

**Reader-list confirmation (Radiation Lead draft, for V&V Lead to accept or amend):**

1. Fields that need ghost fill for radiation **physics** (must stay bit-reproduced): none of `UII`/`UIID`/`QR` ghosts — only `IL` ghosts from `WALL_LOOP1` / face exchange (already in the sweep design).
2. Readers of `UII`/`UIID`/`QR` ghosts: restart I/O only (handled by interior-only compare); no interior physics or wall reader; device/slice/thermocouple readers assumed interior pending V&V sign-off of rows 7–9.
3. Who reads them: `dump.f90` restart and (if unconfirmed) gas-phase output / thermocouple; Radiation Lead owns the sweep write; V&V Lead owns the restart compare rule; Integration Lead owns any future restart-dump path that applies the mask.

**V&V co-sign (2026-10-06, FireX `bee11f0329`; `docs/vv/d073b_ghost_readers.md` @ `575eb02`):** D-073 (b) / D-077 (6) verdicts (1)–(4) **confirmed**. Radiation Lead elects **R-D073b-out (A)**: keep `UII`/`UIID`/`QR` ghosts bit-reproduced when any SLCF/ISOF/PLOT3D/user SM3D requests INTEGRATED INTENSITY or RADIATION LOSS; otherwise those ghosts need not be bit-reproduced. The A-50 / box-split "nothing reads ghost UII/UIID/QR" clause is closable on this basis.

### 7. Status of the stages and what is open

| Stage | State |
|---|---|
| 1 kernels | 3D and 2D Cartesian; K2 CI and clause lint pass; drift check passes (see the commit of the kernel files) |
| 2/3 tests | bitwise test: all six flag sets (O0, O2, O0omp, O2omp, O2omp_off, O2omp_dpd) PASS at 1 and 4 threads (gfortran 14), **1246 cases in 19 groups**, 0 mismatches, 0 vacuous, UIID ghost cells differ in 0 of 1122 full-array compares, padding cells checked 295,906, touched 0; O2omp and O2omp_off also PASS at 8 threads. Mutation run: **69 of 69 mutants killed** (the 43 earlier ones plus 26 for the ILD and STORE kernels, §9.3 and §10). The reorder `WEIGH_CYL*RSA*IL` is bit-equivalent when `WEIGH_CYL` = 1 and is not counted as a mutant (cylindrical is out of the first release). Device: nvfortran 26.9 `-mp=gpu -gpu=cc89,nofma` (RTX 4070), `OMP_TARGET_OFFLOAD=MANDATORY`, 1 and 8 host threads: 1228 of 1246 cases bitwise equal; the other 18 all have `RAD_DIFF_SCHEME` 3 and are in the libm class (§9.4: 212 differing fields, at most 5632 ulp); no other difference, padding untouched, ghost compares 0 differ; the same with the host side built with and without `-Mnoflushz -Mnodaz` (§8). Not rerun for the 1246-case set: the nvfortran host-only build. Build rule for the host side of nvfortran builds: `-Mnofma -Mnovect`, and `-Mnoflushz -Mnodaz` where the data reach the denormal range (§8) |
| D-079 (b) items | ILD* kernels and `STORE_RADIATION_TERMS` kernels: done (§9, §10). `CC_IBM`: design only (§11). Cylindrical: design only, out of the first release (§12). Launch latency: status from Integration Lead SYNTHETIC measurements (§13, cite `docs` `85b4785` stage-1 plan §7.14); ILD/STORE subset re-checked here. |
| open | V&V confirmation of §6; registration of claims L1242/L1248/L1249/L1250/L1238 and of the `EXP` libm entry by the Legacy Mapper; V&V tolerance for scheme-3 libm class (§9.4); L1240/L1244 translated on s5-gen (§16.3; L1240 host, ramp `313c54ef43`, DLO mask hoist f9c9bed577); Wall BR-join review confirmed (`b175a7386b`, no layout change); `CC_IBM` kernels (§11); real-geometry / multi-rank / ILD*/STORE at scale (§13.5) |

## 8. nvfortran build rule found in the D-077 run

Host code compiled by nvfortran at `-O2` vectorizes reductions with reassociated partial sums (`SUM()` and loops of the form `X = X + A(L)`), which gfortran does not. Without `-Mnovect` the host reference differs from the device (and from gfortran) by 1 ulp in the INTERPOLATED wall average (`IL(II,JJ,KK) = sum of ILR / n`) and in the `SUM()` of the open-boundary `Q_RAD_IN` loop (L1243, radi.f90:4965-4974): 34 mismatches of 151 sweep cases and 45 of 160 `open_qin_spec` cases. The device code keeps the source order, so the device result is the ascending-order result. Rule: the host reference of an nvfortran build is compiled with `-Mnofma -Mnovect`; a bitwise claim for any summation loop is made against the ascending source order, and a build that vectorizes the reference is not a valid reference. With `-Mnovect` all suites pass at 1 and 8 threads on host and on the device.

**Denormals (found in the L1248 device run, §15).** nvfortran 26.9 at `-O2` on x86 runs the host code with flush-to-zero and denormals-are-zero (`-Mflushz`, `-Mdaz` behaviour by default: a denormal operand is read as 0 and a denormal result is written as 0); `-Mnoflushz` alone keeps denormal results but still reads denormal operands as 0, only `-Mnoflushz -Mnodaz` gives IEEE gradual underflow (checked with a three-line program). The device keeps denormals in double precision (`nofma` does not change that). With denormal table entries the L1248 kernel on the device differed from the nvfortran host reference in up to 4677 of 11232 elements per case; with `-Mnoflushz -Mnodaz` on every host file all 138 cases are bitwise equal, as they are with gfortran (which keeps IEEE underflow). Rule: an nvfortran host reference is built with `-Mnofma -Mnovect -Mnoflushz -Mnodaz`. The sweep suite gives the same result with and without the two flags (its data do not reach the denormal range), so no earlier sweep result changes. Denormal test data are made from bit patterns (`transfer`), because nvfortran also flushes denormal literals.

## 9. RAD_DIFF_SCHEME > 1: the ILD* kernels (D-079 (b), hand-written, FR-062)

Upstream (`radi.f90` at `bee11f0329`): `RAD_DIFF_SCHEME` is 1 (STEP), 2 (DIAMOND) or 3 (EXPONENTIAL) (`read.f90:10267-10271`); `ILDX`, `ILDY`, `ILDZ` (the directional derivatives of the intensity, allocated only for `RAD_DIFF_SCHEME>1`, 3826-3829) are set to `IL` per angle before the sweep (4425-4429), used by the 3D cell body (upwind read 4533-4541, STEP branch 4621-4633, ALT branch 4635-4714 with the three weight formulas and the negative-correction loop `NEG_ITER = 1..3`), by the wall update WALL_LOOP2 (4764-4773) and by the interpolation into `IL_S` (4832-4841).

### 9.1 Kernels and launch order

All kernels are `bind(c)` in `s5rad_sweep_k2.F90`, with the per-box extent guard and the arrays `(0:NXM+1,0:NYM+1,0:NZM+1,NB)` of the other kernels. Exact signatures:

| Kernel | Upstream lines | Signature |
|---|---|---|
| `s5rad_ild_init` | 4426-4428 | `(NB,NXM,NYM,NZM,BOX_IBAR,BOX_JBAR,BOX_KBAR,IL,ILDX,ILDY,ILDZ)` |
| `s5rad_plane3d_ild` | 4529-4556, 4616, 4621-4714 | `(NB,NXM,NYM,NZM,NPLMAX,KW,JW,D,ISTEP,JSTEP,KSTEP,SOLID_PARTICLES,RAD_DIFF_SCHEME,NCELL,NRA,N,RFPI,BOX_IBAR,BOX_JBAR,BOX_KBAR,NPL,PLN,DX,DY,DZ,EXTCOE,RTE_SOURCE,CELL_INDEX,CELL_SOLID,CELL_ILW,DLX,DLY,DLZ,DLANG,RSA,IL,IL_UP,ILDX,ILDY,ILDZ)` |
| `s5rad_wall2_ild` | 4751-4774 | `(NB,NXM,NYM,NZM,NW,NBR,NRA,NSB,N,IBND,TWO_D,RAD_DIFF_SCHEME,W_BOX,W_BTYPE,W_BR_INDEX,W_IOR,W_IIG,W_JJG,W_KKG,DLN,IL,ILDX,ILDY,ILDZ,INRAD_W,BR_ILW)` |
| `s5rad_ils_ild` | 4831-4841 | `(NB,NXM,NYM,NZM,NICS,NRA,NSB,N,IBND,RAD_DIFF_SCHEME,ILS_BOX,IIO_S,JJO_S,KKO_S,IOR_S,IL,ILDX,ILDY,ILDZ,ILS)` |

`IOR_S(NICS)` is the orientation `IOR_S` of the other-mesh record used in 4832-4841 (the table the `ILS` kernel did not need). `DLANG(NRA)` is the table `DLANG` of the upstream angle setup. `RAD_DIFF_SCHEME` is passed by value; `plane3d_ild` is valid for 1, 2 and 3 (scheme 1 takes the STEP branch 4621-4633 of the same kernel), so the host may use either `plane3d` (scheme 1, smaller argument list) or `plane3d_ild`.

Dispatch (host) for `RAD_DIFF_SCHEME>1`, per angle: `cell_ilw_reset`, `wall1`, **`ild_init`**, the planes with **`plane3d_ild`** (one launch per plane as before), **`wall2_ild`** instead of `wall2`, `wall3`, `uiid_acc`, **`ils_ild`** instead of `ils`. That is one more launch per angle than scheme 1 (53 instead of 52 for 16^3 boxes). The 2D plane kernel is unchanged (below). `ild_init` is a whole-array copy and carries the per-box extents like the other whole-array kernels, so the padding rule of Addendum 1 §2 holds (the tests check the padding of `ILDX/Y/Z`).

`plane3d_ild` has local arrays (`PLX(3)`, `NEGATIVE_MASK`), so it uses the loop macro `S5R_LOOP3A` (`target teams distribute parallel do collapse(3)`, host: `parallel do collapse(2)`) instead of `target teams loop`: a `loop` construct may not privatise arrays that way (K2 rule 7, `generator-howto.md`). nvfortran 26.9 accepts it. The flag `APPLY_ALT_SCHEME = RAD_DIFF_SCHEME>1` (4616; the CC_IBM clause 4617-4619 is omitted, §11) is set before the loop directive (lint rule CL-03).

### 9.2 Upstream quirk reproduced, not changed (for the FDS owners)

The 2D Cartesian sweep (4469-4492) and the cylindrical sweep (4431-4468) never write `ILDX/ILDY/ILDZ`. With `TWO_D` (or cylindrical) and `RAD_DIFF_SCHEME>1`, WALL_LOOP2 (4764-4773) and the `IL_S` interpolation (4832-4841) therefore read the copy of `IL` taken before the sweep (ambient or wall values), not the swept intensity. The kernels reproduce this bit for bit (the group `ild_2d` and the 2D cases of `ild_3d_band`, `ild_stacked`, `ild_zero_cos` test it, with a coverage line `ild_2d_boxes`). It looks unintended; it is the owners' call and nothing in the AMR code changes it.

### 9.3 Tests and mutants

The reference is the verbatim upstream sweep (now 410 lines, scheme selected at run time). New groups: `ild_3d_gray`, `ild_3d_band`, `ild_2d`, `ild_stacked` (three boxes of different extents, mixed gray/band), `ild_zero_cos` (direction cosines of exactly zero: the EXPONENTIAL optical depth is infinite for `EXTCOE>0` and zero for `EXTCOE=0`), `ild_contrast` (199 small boxes with cell optical depths up to 10 and small sources, scheme 2 and every tenth case scheme 3: cascades of negative downwind intensities in which the third iteration of the negative-correction loop changes the result; seed 19120, a one-cell box whose intensities all end clipped to 0, is skipped as vacuous and the test now names any vacuous case), and the scheme-1 angles of the existing groups that use `plane3d_ild` on every second angle. Coverage lines: DIAMOND boxes, EXPONENTIAL boxes, 2D ILD boxes, cells where the negative-correction loop binds (`ILD=0` with `IL>0`, 1505 cells at the final angle; DIAMOND 37 boxes, EXPONENTIAL 55, 2D 39, STORE 62), data that hit the `TWO_EPSILON_EB`/`TWENTY_EPSILON_EB` thresholds of the EXPONENTIAL weights (the `pos()` generator alternates 1e-300 and 1e-15 optical depths). **Surviving mutants of the interrupted round.** The earlier round ended with 3 survivors out of 69; their names were not recorded. With the present tests all 69 are killed. `neg_third_iteration_missing` (`DO NEG_ITER = 1, 2`; the old name `neg_iterations_one` was wrong) is killed only by `ild_contrast`: without that group the third iteration never changes a result and the mutant survives. That makes it a real defect the test can see, not an equivalent mutant, so it was kept and not replaced. No mutant had to be proved equivalent in this round. 26 mutants cover: missing `ild_init` (also the Z component alone), wrong weight (DIAMOND Z, EXPONENTIAL axis, small threshold, missing `HUGE_EB` guard), two instead of three negative-correction iterations (`neg_third_iteration_missing`), a mask that does not persist across iterations, a correction not applied in Z, `ILDY` taken from X, upwind read from `IL` instead of `ILD*`, regrouping of the STEP sums and of the ALT upwind sums, `wall2_ild` axis swapped or ignoring `ILD*`, `ils_ild` with a signed `IOR_S` or swapped axis, padding writes, and for the STORE kernels (§10) overwrite instead of accumulate, wrong slice, swapped factors, only one of the two arrays zeroed.

### 9.4 EXP: libm class

`plane3d_ild` with scheme 3 is the first sweep kernel with a libm call (`EXP`, 4660-4690). Per D-072 a libm function is bitwise only against the same libm; the device `EXP` is not the host `EXP`. Measured on the test machine's GPU in the earlier round (nvfortran 26.9, 2,000,000 arguments in [-20,0]; not repeated): 9,970 differ, maximum 1 ulp. In the sweep tests the host builds (gfortran at 1, 4 and 8 threads on the 1246-case set; the nvfortran host-only build was last run on the earlier 1047-case set) are bitwise equal to the reference, and 18 of the 1246 cases differ on the device: all 18 have `RAD_DIFF_SCHEME` 3 (3 of the 10 scheme-3 cases of `ild_3d_gray`, 5 of 8 of `ild_3d_band`, 2 of 3 of `ild_stacked`, 6 of 18 of `ild_zero_cos`, 2 of the 20 scheme-3 cases of `ild_contrast`; the 2D cases never call `EXP` and are bitwise equal), every difference is in a floating-point field, padding and counters are exact. With `-DS5_LIBM_CLASS` (device runs only) the test reports these cases as `LIBM-CLASS` with the largest ulp distance over `IL`, `UIID`, `UII`, `QR`, `QR_W`, `ILW`, `INRAD_W`, `IL_S`: **5632 ulp (about 1e-12 relative)**, 212 differing fields in 18 cases; the 1-ulp `EXP` difference is amplified by the recurrence along the sweep and by the negative-correction clipping. Rule proposed: scheme 1 and 2 are bitwise on the device; scheme 3 is bitwise on the host and in the libm class on the device, tolerance to be set by the V&V Lead (question in the report). Registry: the draft entry `[[libm]] kernel = "s5rad_plane3d_ild"`, ruling D-079, status proposed, is in `s5rad_sweep.toml` (`[meta.libm_registry]`); the shared `kernel_registry.toml` is the Legacy Mapper's file and does not have it yet (the K2 CI check passes against a copy that has it, and prints the missing entry as INFO against the shared one).

SASS check of the device object (`cuobjdump -sass`, earlier round, not repeated): the `DFMA` instructions are inside the IEEE double division expansions (`__cuda_sm20_div_rn_f64_full`, reciprocal slow path, a few in `plane3d`, `plane2d`, `wall1`) and, in `plane3d_ild`, inside the inlined `EXP`; the kernel arithmetic itself is contracted nowhere (build flag `nofma`, bitwise equal to the host for every scheme-1 and scheme-2 case).

## 10. STORE_RADIATION_TERMS kernels (D-079 (b), hand-written, FR-062)

Upstream: `RADIATION_EMISSION` and `RADIATION_ABSORPTION` (`mesh.f90:49-50`, `0:IBP1,0:JBP1,0:KBP1`) are the diagnostic fields filled only when `STORE_RADIATION_TERMS`: zeroed at the start of the update for wide band or WSGG (3877-3880, statements 3878-3879), accumulated per band after each band (4950-4953: `RADIATION_EMISSION += KFST4_GAS`, `RADIATION_ABSORPTION += KAPPA_GAS*UIID(:,:,:,IBND)`), and assigned in the gray form after `QR` (4984-4987).

| Kernel | Upstream | Host condition | Signature |
|---|---|---|---|
| `s5rad_emis_zero` | 3878-3879 | `STORE_RADIATION_TERMS .AND. (WIDE_BAND_MODEL .OR. WSGG_MODEL)`, once per update before the first band | `(NB,NXM,NYM,NZM,BOX_IBAR,BOX_JBAR,BOX_KBAR,RADIATION_EMISSION,RADIATION_ABSORPTION)` |
| `s5rad_emis_wb` | 4951-4952 | `STORE_RADIATION_TERMS`, after the band is finished (after `qr_wb`) | `(NB,NXM,NYM,NZM,NS,ISL,BOX_IBAR,BOX_JBAR,BOX_KBAR,KAPPA_GAS,KFST4_GAS,UIID,RADIATION_EMISSION,RADIATION_ABSORPTION)` |
| `s5rad_emis_gray` | 4985-4986 | `STORE_RADIATION_TERMS .AND. .NOT.(WIDE_BAND_MODEL .OR. WSGG_MODEL)`, after `qr_gray` | `(NB,NXM,NYM,NZM,BOX_IBAR,BOX_JBAR,BOX_KBAR,KAPPA_GAS,KFST4_GAS,UII,RADIATION_EMISSION,RADIATION_ABSORPTION)` |

All are one-thread-per-cell, no reduction, per-box extents, padding untouched. `KAPPA_GAS` and `KFST4_GAS` in the kernel calls are the same per-band slices that `qr_wb` / `qr_gray` receive. The group `store_terms` and the `store` flag of the other groups run them with scheme 1, gray, wide band and stacked boxes; the mutants `emis_zero_missing`, `emis_zero_only_emission`, `emis_wb_overwrite`, `emis_wb_wrong_slice`, `emis_gray_swapped`, `store_flag_ignored` (driver) and the two padding mutants `pad_emis_wb`, `pad_emis_gray` kill the host-condition, zero-order, slice and padding errors. Bitwise on host and device in every flag set. Note for the driver: the band form accumulates in band order; the order of the bands is the host loop and must not be changed.

## 11. CC_IBM: design (not implemented; design only until the geometry tables exist)

Scope decision: host-only is not possible (the cut-cell terms sit inside the cell body), a kernel is needed, and it needs data the AMR code does not yet carry (`CCVAR`, `FCVAR`, `CUT_FACE`, `CUT_CELL`, `CFACE` with `BR_INDEX`). It is therefore a design item. No ADR is needed for the kernels themselves (same rules as D-073), but the data layout needs the Integration Lead and the geometry owner to agree to the table below; the host refuses `CC_IBM=T` until the kernels exist.

### 11.1 What CC_IBM adds to the sweep (upstream lines)

1. Cell body, 4561-4613: a solid cut-region cell is skipped (`CCVAR(I,J,K,CC_CGSC)==CC_SOLID`); for a cell with a cut face (`INDCF=CCVAR(I,J,K,CC_IDCF)>0`) the faces of `CUT_FACE(INDCF)` are walked in `IFACE` order: with `DLF = DOT_PRODUCT(GET_CFACE_RAD_NVEC(CFA),DLA)`, an upwind face (`DLF>0`) adds `DLF*AREA*ILW(N)` to `AILFU`, a downwind face (`DLF<0`) adds `ABS(DLF)*AREA` to `AFD`. The six Cartesian areas are scaled by the aperture `ALPHA_CF` of the face (or set to zero when the face is solid) and `VC` by `ALPHA_CC`; the cell is skipped when both `EXTCOE` and the downwind area sum are below `TWO_EPSILON_EB`.
2. The weights: gas cut cells use the STEP form with `AILFU` and `AFD` (4616-4619 sets `APPLY_ALT_SCHEME = .FALSE.` for `CCVAR(I,J,K,CC_IDCC)>0`, so ILD* never applies there; 4623-4624 add `AFD` to `A_SUM` and `AILFU` to `AIU_SUM`).
3. CFACE walls: `CFACE_LOOP1` (4375-4384) sets `ILW(N)` of every upwind internal CFACE from `OUTRAD_F` and `INRAD_F` before the sweep; `CFACE_LOOP2` (4791-4802) and `CFACE_LOOP3` (4803-4812) update `INRAD_F` and `ILW(N)`/`ILW(ANGLE_INC_COUNTER)` after the sweep from the intensity of the gas cell `(IIG,JJG,KKG)`. The setup of `OUTRAD_F` (4251-4274) and the extinction setup in the translated kernels (4009-4186 `IF (CC_IBM)` skips and `ALPHA_CC` scaling) are in the generated `rad_*` kernels and have their own CC_IBM clauses.

### 11.2 Tables (device, built per regrid)

| Table | Content | Used by |
|---|---|---|
| `CC_CGSC(I,J,K,IB)` int8 | solid/gas/cut flag of the cell | cell body skip |
| `CC_CUT_START(I,J,K,IB)`, `CC_CUT_LIST(*)` | compressed list: for every cell with `CCVAR(...,CC_IDCF)>0` the offset into a per-face list, in the **upstream `IFACE` order** | `AILFU`/`AFD` sums |
| per-list entry `CF_ID(*)` | index of the CFACE (slot in `BR_ILW` through `CF_BR_INDEX`), `CF_AREA`, `CF_NX/NY/NZ` (from `GET_CFACE_RAD_NVEC`, precomputed in the same expression order) | face sums |
| `AP_FACTOR(face,axis)` | the six aperture factors per cell, already combined (`ALPHA_CF`, or 0 for a solid face), stored per cell as `AXU..AZD` multipliers | areas |
| `ALPHA_CC(I,J,K,IB)` | volume factor, 1 for non-cut cells | `VC` |
| `CELL_IS_GASCUT(I,J,K,IB)` | `CCVAR(I,J,K,CC_IDCC)>0` | STEP-only flag, 4617-4619 |

The cell loop then becomes the existing `plane3d`/`plane3d_ild` cell body with an extra block: a loop over the cell's face list in list order (a sequential per-thread loop; the order of the additions is the order of the upstream loop, so the sums are bitwise, the same argument as for the wall tables), the multipliers, and the skip test. `GET_CFACE_RAD_NVEC` and `DOT_PRODUCT(...,DLA)` are evaluated per angle per face: the dot product is `NX*DLX+NY*DLY+NZ*DLZ` in the order of the intrinsic (to be pinned by the drift tool with an omit entry and a test that fixes the order for gfortran and nvfortran).

### 11.3 Kernels (planned)

1. `s5rad_cfaces_pre` (CFACE_LOOP1): one thread per CFACE, race-free (each CFACE has its own `ILW` slot).
2. `s5rad_plane3d_cc` / `_cc_ild` (cell body above), one launch per plane as now; the ILD variant keeps the STEP form for gas cut cells.
3. `s5rad_cfaces_post` (CFACE_LOOP2 and 3): one thread per CFACE: it updates its own `INRAD_F` and its own `ILW` slots; the read is `IL` of the gas cell, no two CFACES write the same element, so no atomics. The three loops of upstream are three passes; the kernel keeps the two-step order (add, then replace, then subtract) of 4799-4801.

Dependencies: (a) the CFACE and cut-cell arrays must exist in the AMR geometry (not in the current fork plan), (b) CFACE `BR_INDEX` slots must join the `BR_ILW` table (doc 05 §3 rows for CFACE), (c) the generated extinction and source kernels must have their CC_IBM clauses (the translated `rad_*` kernels blank them today: a CC_IBM run is refused), (d) the omit pin for 4617-4619 in `s5rad_sweep.toml` is replaced by copy blocks. Test: the same bitwise scheme against the verbatim sweep with random cut-cell data (solid, gas, cut cells with 1 to 4 faces, apertures between 0 and 1, `DLF` of both signs and exactly zero), mutants for the face order, the sign of `DLF`, the aperture on the upwind versus downwind side, `ALPHA_CC` and the STEP-only flag.

## 12. Cylindrical: design only (status unchanged)

D-042 defers cylindrical in AMR mode and D-073 (5) puts it out of the first release (cylindrical meshes stay FDS-only in AMR mode); neither ruling has been changed, so nothing here is implemented, and the host refuses `CYLINDRICAL=T`. For the record, the pieces a later implementation needs: the axisymmetric sweep 4431-4468 (one plane `J=1`, `VC`, `AXU`, `AXD` from `R(I-1)`, `R(I)` and `DPHI0`, the `MODULO(N,NRP(1))` symmetric-overhang zeroing of `AYU/AYD`, no ALT scheme and no ILD* write, §9.2), the Y-ghost copy 4722-4743 for the angles that wrap around the axis (a data dependence between angles of one `IPHI` family: the angle chain of option B, closed by D-073 (1) `B=1` for cylindrical), the CYLINDRICAL branch of WALL_LOOP1 (`DLB`, the `RPI` weights) and `WEIGH_CYL` in `uiid_acc` (already an argument, fixed 1 for Cartesian). A `s5rad_plane_cyl` kernel would be 2D-like (thread space `(box,k')`, plane index along `I+K`), but the angle-to-angle Y-ghost copy forces `B=1` and a host loop in `IPHI` order; no kernel work is planned until a ruling reopens the item.

## 13. Launch latency: status (D-079 (b); cite Integration Lead)

### 13.1 Source of the numbers

The authoritative SYNTHETIC launch-cost measurements are the Integration Lead's stage-1 harness (`stage1/rad_sweep/`, `docs/amrex/stage1-gpu-spike-plan.md` §7.14), re-pinned to s5-gen `d054f283f0` (ILD*/STORE tip; scheme-1 launch costs re-run, unchanged within ~2% vs prior pin `3814e1f893`). Inputs are synthetic gray 3D boxes (`RAD_DIFF_SCHEME = 1`), internal wall rows included, blocking launches, managed memory, test-machine GPU (RTX 4070, cc 8.9, nvfortran 26.9). Byte identity of the batched path against box-by-box and against a natural-order sweep holds on host and device (120 of 120 dump files). Those figures supersede any earlier private bench table for the same configurations.

### 13.2 Per-angle cost (Integration Lead, SYNTHETIC)

| boxes | cells | total per angle | notes |
|---|---|---|---|
| 1 of 16^3 | 4,096 | 0.361 ms | plane floor ~6.7 us |
| 8 of 16^3 | 32,768 | 0.356 ms | launch-bound |
| 32 of 16^3 | 131,072 | 0.423 ms | |
| 64 of 16^3 | 262,144 | 0.574 ms | |
| 128 of 16^3 | 524,288 | 1.018 ms | |
| 256 of 16^3 | 1,048,576 | 2.593 ms | 2.359 ms without internal rows |
| 512 of 16^3 | 2,097,152 | 6.316 ms | |
| **1024 of 16^3** | **4,194,304** | **13.13 ms** | ~2 GB managed |
| 128 of 32^3 | 4,194,304 | 13.39 ms | same cell count, different box shape |

At 4M cells the launch floor is 2.7% of the time (52 × 6.7 us ≈ 0.35 ms of 13.1 ms); the cost is in the kernels. Pre/post wall kernels are about 25% of the pass at that size. Scaling from 1.05M to 4.2M cells is slightly super-linear (cause not investigated by the harness). Memory estimate: ~490 bytes/cell on this synthetic level; 4.19M cells fit in one pass on an 8 GB card.

### 13.3 What the harness did not cover

Real geometry (wall counts and `NPL` shapes of a production level), multi-rank (`IL_S`/`IL_R` exchange, shared `NICR`), a full 24-angle pass at 4M cells, FMA and device negative controls, cylindrical variants, and the cause of the super-linear scaling above 1M cells. Also not in that harness pin: `RAD_DIFF_SCHEME` 2/3 (ILD*), `STORE_RADIATION_TERMS`, CFACE, and `CC_IBM`.

### 13.4 ILD* / STORE subset check (this work)

`test/s5_sweep_bench.F90` reuses the bitwise case builder and times the full launch sequence of Addendum 1 §2–§3 (including `ild_init`, `plane3d_ild`, `wall2_ild`, `ils_ild`, `emis_*`). On the same test-machine GPU, three quiet runs, medians of the per-angle sum of blocking launch times (synthetic walls denser than the Integration Lead's generator; not a replacement for §13.2):

| config | launches / angle | sum per angle (median of 3) |
|---|---|---|
| 8 × 16^3, scheme 1 | 52 | 0.405 ms |
| 8 × 16^3, scheme 2 | 53 | 0.709 ms |
| 8 × 16^3, scheme 3 | 53 | 0.944 ms |
| 64 × 16^3, scheme 1 | 52 | 0.821 ms |
| 64 × 16^3, scheme 2 | 53 | 1.642 ms |
| 64 × 16^3, scheme 3 | 53 | 2.210 ms |

Readings: (1) scheme-1 at 8 boxes agrees with the Integration Lead order of magnitude (0.36 ms there; the excess here tracks the denser synthetic wall set). (2) ILD* (scheme 2/3) raises the per-angle sum by about 1.7× / 2.3× at 8 boxes and 2.0× / 2.7× at 64 boxes against scheme 1: the cell body is heavier and the argument list longer. (3) At 64 boxes of 16^3 the ILD plane kernels are already compute-bound relative to the launch floor.

Two measurement pitfalls that apply to any driver using managed memory: (a) first-touch warm-up of at least three angles on the same arrays; (b) no host reads of managed direction cosines or plane tables inside the angle loop (they migrate pages; keep launch parameters in host-only memory).

### 13.5 What remains to measure

(a) real geometry and real wall / `NPL` shapes; (b) multi-rank `IL_S`/`IL_R` exchange and shared `NICR`; (c) ILD* / `STORE` / CFACE / `CC_IBM` at 4M cells (the subset above stops at 64 boxes); (d) a full 24-angle pass at 4M cells; (e) device-resident (`omp_target_alloc`) vs managed memory; (f) the angle-batching decision (stage 8 stays closed): at the Integration Lead numbers the launch-bound region extends to about 64 boxes of 16^3 for scheme 1; above about 128 boxes the plane and wall kernels dominate; at 4M cells the launch count is irrelevant. That is a fact for the Architect to rule on, not a recommendation made here.

## 14. Closing the open points of the Integration Lead's descriptor questions and of the review follow-ups

The seven descriptor questions (stacked `UII`/`QR`, meaning of the `ILS` indices, `W_NIC_MIN`/`W_NIC_MAX`, internal OPEN walls, `W_RAMP_VAL`, the 2D path, the padding cost of the whole-array kernels) were answered in the stage 1 to 3 work (`UII`/`QR`/`QR_W` stacked kernels; the `ILS` indices as in §2; `W_NIC_MIN/MAX` in the `wall1` signature; the internal OPEN rows in `open_ilw_zero` via `W_BTZ`; `W_RAMP_VAL` once per angle in `wall1`; the 2D kernels; per-box extents `BOX_*BAR` so that padding is not computed). Nothing is open on them. The padding-cost question is answered by the per-box extents: the whole-array kernels loop over `BOX_IBAR/JBAR/KBAR(IB)`, and the tests check that the padding cells are untouched.

Follow-ups R1 to R4 of `solid/08-new-wall-kernel-review.md` §7: **all closed on Radiation's side** (docs `05` Additional rules + this paragraph; no `s5_rad_*` code change; Wall owns wall tables). **R2** row ownership, **R3** OQ-S6 owner+fan-out (emissivity per surface), **R4** pure `W_B1_PRESENT` are in `05-br-ilw-table-spec.md`. **R1** (compiler scope of the L1243 `SUM`-order claim): the claim "the ascending-order sum is the result" holds for gfortran at every flag set tested and for ifx `-O2 -fp-model=precise` (D-079 (a): 0 mismatches); it holds for nvfortran only with `-Mnovect` on the host (§8: without it 34 of 151 sweep cases and 45 of 160 `open_qin_spec` cases differ by 1 ulp from the device and from gfortran) and always on the device. The claim therefore has this scope: gfortran 13 (all flag sets), ifx `-O2 -fp-model=precise`, nvfortran 26.9 with `-Mnofma -Mnovect` (host) and on the device. A compiler or flag set outside that list needs the `open_qin_spec` test before it is used as a reference; with D-079 flags no oneAPI/ifx tolerance class is required for open-boundary cases.

## 15. INTERPOLATE_IL wall part (claim L1248) and L1245 note

### 15.1 Scope

Upstream `INTERPOLATE_IL` (`radi.f90` at `bee11f0329`): angle tables 3623-3645 (host), wall loop 3651-3662 (L1248), CFACE 3665-3676 (L1249), particles 3679-3695 (L1250), neighbour-mesh `IL_R` 3698-3713 (unclaimed). This section covers L1248 only. Code: `amrex/s4_mass/s5_gen/s5_rad_interp/`.

### 15.2 Why a ping-pong

Angle `N` reads the old values of other angles through `IDX(N,:)`. The upstream whole-array form keeps `ILW_OLD = BR%BAND(IBND)%ILW(:)` then writes `BR%BAND(IBND)%ILW`. On the slot table that is `BR_ILW_OLD(NBR,NRA,NSB)` allocated by the host for the call: copy then interpolate, in that order. The sum is ascending from `+0`, product grouped `(W*ILW_OLD)*RSA`, then `/RSA(N)` (same as radi.f90:3656-3660).

### 15.3 Signatures

```
subroutine s5rad_interp_weights(NRA,DLANG,DLANG_OLD,IDX,W)                                 ! host (3623-3645)
subroutine s5rad_interp_ilw_host(NW,NBR,NRA,NSB,W_BR_INDEX,IDX,W,RSA,BR_ILW)                ! host (3651-3662)
subroutine s5rad_ilw_old_copy(NW,NBR,NRA,NSB,W_BR_INDEX,BR_ILW,BR_ILW_OLD) bind(c)          ! K2
subroutine s5rad_ilw_interp(NW,NBR,NRA,NSB,W_BR_INDEX,IDX,W,RSA,BR_ILW_OLD,BR_ILW) bind(c)  ! K2
```

`N_INTP = 4` (line 3602). `IDX(NRA,N_INTP)`, `W(NRA,N_INTP)`. Drift: rename sets `ident`, `intp_host`, `intp_elem`; omit hashes pin CFACE, particles, `IL_R`.

### 15.4 BR_ILW interface

s5-wall `[policy.wall.slot_tables]` (`45389a05ec`): extent `NBR,NRA,NSB`, slot fastest, gather `W_BR_INDEX`. The generator treats the table as read-only; L1248 writes it (hand-written kernels cover the write). Driver: allocate `BR_ILW_OLD`, launch copy then interp (or call the host routine alone). Slot uniqueness: `wall_checks.check_br_index`.

### 15.5 Tests (verified)

138 cases in 7 groups (random, tiny_nra NRA 2–4, identity, permuted, ties, special, big); six flag sets PASS at 1 and 4 threads; 22 of 22 mutants killed (including driver `MUT_NOCOPY` / `MUT_COPY_AFTER` and FMA with `-mtune=haswell`). Device (RTX 4070, nvfortran 26.9, `OMP_TARGET_OFFLOAD=MANDATORY`, 1 and 8 threads): PASS once the host reference uses `-Mnofma -Mnovect -Mnoflushz -Mnodaz` (§8: denormals from bit patterns). Drift, K2 CI, lint, registry: PASS.

### 15.6 Blocked / not owned here

- CFACE loop (L1249) and particle loop (L1250): same arithmetic on `CF_BR_INDEX` / `LP_BR_INDEX`; next.
- `IL_R` interpolation (3698-3713): device-resident under a device sweep; no owner yet (ask Legacy Mapper / Integration Lead).
- Generator write access to `BR_ILW` and the `BR_ILW_OLD` allocation: Wall Loops Engineer / driver.

### 15.7 L1245 (RADF file write, radi.f90:5044-5059)

Excluded from the device plan (`loop_claims`: not planned (host)). It is **not** purely host-only in AMR mode: it needs `RF%IL_SAVE`, filled per angle by 4884-4895 inside the sweep (omitted from the sweep kernels). Options for the driver: (a) a per-angle gather kernel into a host buffer when `N_RADF>0`, or (b) refuse `N_RADF>0` in AMR mode. Also needs `TMP` on the host and a mapping of RADF `(I,J,K)` to boxes. Upstream quirks kept as-is if (a) is chosen: `IL_SAVE` keeps the last band's values and stale angles. No kernel work here.

## 16. INTERPOLATE_IL CFACE/particles (L1249, L1250), ADD_VOLUMETRIC_HEAT_SOURCE (L1238), and particle-table blocks (L1240, L1244)

### 16.1 L1249 / L1250 (host first)

Same arithmetic as L1248 on the slot table through `CF_BR_INDEX` / `LP_BR_INDEX`. Host routines in `s5_rad_interp/`:

```
subroutine s5rad_interp_ilw_cf_host(NCF,NBR,NRA,NSB,CF_BR_INDEX,IDX,W,RSA,BR_ILW)
subroutine s5rad_interp_ilw_lp_host(NLP,NBR,NRA,NSB,LP_BR_INDEX,LP_ORIENTATION,LP_SOLID_OR_TARGET,IDX,W,RSA,BR_ILW)
```

Particle filters are flat: `LP_SOLID_OR_TARGET(IP)/=0` and `LP_ORIENTATION(IP)>0`. Drift covers regions 3665-3676 and 3679-3695; glue for pointer aliases and the two IFs. Tests: the L1248 suite extended (138 cases; CFACE and particle slots exercised; host_cf_lp path PASS). Geometry SoA / particle SoA production tables are not required for the arithmetic test (gathers are synthesised); production `CF_BR_INDEX` joins when CFACE slots exist (doc 05).

### 16.2 L1238 ADD_VOLUMETRIC_HEAT_SOURCE (host, gather/scatter)

Code: `s5_rad_hrr/`. For each `INITIALIZATION` with `HRRPUV>0`: host evaluates `TIME_RAMP_FACTOR`; gather touched `(I,J,K)` by DEFAULT box or CONE; scatter MODE update onto `Q` / `KFST4_GAS`. CC_IBM cut-cell branch (5132-5140) omitted (geometry deferred). Tests: 40 cases, modes 0/1/2, box and cone, PASS. **Superseded by §17.7**: this first version is not bitwise and its test reference was not verbatim. The corrected version is on s5-land (`6ccf137752`).

```
subroutine s5rad_hrr_gather_box(...)
subroutine s5rad_hrr_gather_cone(...)
subroutine s5rad_hrr_scatter(IBAR,JBAR,KBAR,MODE,NTOUCH,IT,JT,KT,DELTA,CHI_R_DELTA,Q,KFST4_GAS)
```

### 16.3 L1240 / L1244 landed (particle tables from Mesh Data)

Mesh Data delivered production flat gathers in `s5_part_gather/` (`6eaad9ea70` / `dd02eef1c5` / `0185448a9c`). Radiation owns the loop landing under `s5_rad_lp/` (does not edit `s5_part_gather`).

- **L1240** (radi.f90:3894-3901 at pin; inventory 3893-3900): **host** `s5rad_lp_qin_zero` (bind(c), name and signature unchanged), verbatim copy in `s5rad_lp_host.F90` next to L1244. Moved off the device by the simplicity rule: L1244 writes the same `LP_Q_RAD_IN` on the host, so a K2 would copy the field over and back only to zero it. `s5rad_lp_k2.F90`, its toml row and the K2 test path are removed (no `k2_ci_check` / `DEFAULT_K2_GLOBS` entry covered `s5_rad_lp`). Tables: `LP_B1_PRESENT`, `LPC_SURF_INDEX`, `SF_TMP_GAS_FRONT`. Ours is the only definition of the C symbol; Mesh Data dropped the twin in `d283813c68`.
- **L1244** (radi.f90:4993-5038 at pin; README 4989-5034): host `s5rad_lp_qin_in` (doc 05 Q6). Filter `LP_SOLID_OR_TARGET`, `SF_TMP_GAS_FRONT` gate, oriented (`BR_ILW(LP_BR_INDEX,:,IBND)` masked `SUM`) vs isotropic (`0.25*LP_UII + LP_EFLUX`). `BR_ILW` layout unchanged (Wall `45389a05ec`); Wall reviews the particle join only. **Orientation ramp (5014-5030) implemented**: verbatim in the per-particle loop (`IS_PARTICLE_ORIENTATION_RAMP`, local `COS_DLO_ARR(N) = DOT_PRODUCT(IN_TEMP_ORIENTATION(1:3,N),DLANG(1:3,N))`); nothing shared is overwritten, so particles of one orientation in different ramp states are independent (the earlier "host may overwrite a `DLO_ORIENTATION` row" idea is wrong in exactly that case and is dropped). New tables `LP_INIT_INDEX`, `IN_ORIENTATION_RAMP_INDEX(3,NIN)`, `IN_TEMP_ORIENTATION(3,NIN)`, `DLANG(3,NRA)`. Host glue: `EVALUATE_RAMP` → `LP_EFLUX` (5000-5009, pinned omit); `UII` → `LP_UII`; 5019 (pinned omit) `CALL CALC_PARTICLE_TEMP_ORIENTATION(T,IN%ORIENTATION_RAMP_INDEX(1:3),IN_TEMP_ORIENTATION(:,N))` once per ramped `INITIALIZATION` record (depends only on `T` and the record, so equal to upstream's call per particle). `LP_EFLUX` / `LP_UII` are gathered for filtered-out particles too (fine on the host; a device version must filter first). **DLO mask hoist (Wall-requested, owner's simplicity rule; supersedes prior decline):** radi.f90:5025-5032 omitted as a mechanical refactor — orientation mask computed once per oriented particle into local `DLO_MASK(NRA)` and reused as `SUM(BR_ILW(S,:,IBND), MASK=DLO_MASK)` for every band (same FLOPs per band SUM, no repeated comparisons; `BR_ILW` join and ramp logic unchanged). Device K2 deferred until deposition design settles.
- **Evidence (s5-gen `313c54ef43`; DLO mask hoist f9c9bed577).** Drift OK (12 copy blocks, 3 pinned omits including radi.f90:5025-5032 hoist). `test/run_lp.sh all nv`: L1240 120 cases and L1244 120 cases × 5 flag sets (gfortran O0 / O2 / O0omp / O2omp with `-ffp-contract=off`; nvfortran `-O2 -mp -Mnofma -Mnovect -Mnoflushz -Mnodaz`) × 1 / 4 threads, bitwise PASS (hoist bitwise-equivalent to per-band mask). L1244 cases carry `ORIENTATION_RAMP_INDEX /= 0` (records (0,0,0), (0,0,k), two different all-positive ramps, rotating over record numbers) and particles of one orientation in different ramp states in both orders (coverage counters for each). Mutants 18/18 killed at O2omp and nv (6 L1240, 6 L1244, 6 ramp: dropped ramp, `ANY`→`ALL`, mask sign, INIT gate, wrong record, flag not reset per particle). A one-off check with the row-overwrite design (DLO row replaced by the last ramped particle's cosines) fails all 120 cases.

**Register (Legacy Mapper, 2026-10-06, updated for s5-gen f9c9bed577).** loop_claims L1240/L1244: `translated`, owner AMR Radiation Lead; the flat gather tables (`s5_part_gather/`) stay owned by the GPU Mesh Data Loops Engineer and are consumed read-only. Code on s5-gen f9c9bed577 (ancestor of s5-gen; ramp landing `313c54ef43`; first landing `b175a7386b`; gathers `6eaad9ea70` / `dd02eef1c5` / `0185448a9c` and Mesh Data pack asserts `a71c662376` also ancestors). Pins rechecked at `bee11f0329`: 3894-3901 is the `LAGRANGIAN_PARTICLE` / `LPC` `Q_RAD_IN` zero loop; 4993-5038 is `PARTICLE_LOOP` (solid / massless-target incoming flux). Routes: L1240 host (K2 dropped; `s5rad_lp_k2.F90` gone at the s5-gen tip), L1244 host-first (doc 05 Q6). Evidence as registered: the bullet above. Review items from the first entry (`7538bb4`): (1) `LP_BR_INDEX` → `BR_ILW` join: Wall review done, confirmed; `BR_ILW` stride kept; asserts `a71c662376`. (2) Orientation ramp: implemented and tested; **DLO mask hoist done** (Wall follow-up under owner's simplicity rule; supersedes prior decline). (3) Duplicate bind(c) `s5rad_lp_qin_zero`: closed. Radiation keeps it in `s5rad_lp_host.F90`; Mesh Data dropped the twin in s5-gen `d283813c68` (single definition at `s5rad_lp_host.F90:21`). Not in `ported.toml` (port gate open); L1240 and L1244 are ready for `ported.toml` when the gate opens.

L1245 remains excluded (not planned host; see §15.7).

## 17. Next round (loop_claims 670de64): L1251, L1247, L1241, L1246 and WSGG, L1234-L1237; L1238 correction

All loops in this round take about 0.000% of run time in the inventory profile, so this section records dispositions and tests, not hot kernels. Lines are at `bee11f0329`.

### 17.1 Inventory line offsets

The inventory and loop_claims lines match the pin for 3500-3713 (L1234-L1237, L1247-L1251) and for L1238 (5098-5146). Three groups are off:

| Claim | Inventory | Pin |
|---|---|---|
| L1240 / L1241 | 3893-3900 / 3901-3907 | particle loop 3894-3901, CFACE loop 3902-3908 (one line lower) |
| L1246 `GET_KAPPA` | 5190-5204 | `DO N` loop 5194-5208 (function 5174-5212) |
| L1252 `KAPPA_WSGG` | 5244-5246 | `DO NN` loop 5248-5250 |
| L1230-L1233 `A_WSGG` | 5271-5302 | `DO MM` loops 5275-5277, 5286-5288, 5293-5295, 5300-5306 |

The code in this round uses the pin lines. Drift checks assert the first statement of every copied range.

### 17.2 L1251: neighbour-mesh `IL_R` interpolation (3698-3713)

**Disposition: host now (FDS-mode path unchanged). AMR replacement: Role 3 lagged `ILR` ghost table (`fdsrt::build_il_ghost_table` + `fdsrt::exchange_il_ghost`), no Radiation kernel.** Sweep contract unchanged: `s5rad_ils` packs; AMR lag via `exchange_il_ghost` (FDS-mode / single-box harness: whole-array `ILS→ILR` copy); `s5rad_wall1` reads the mean of `ILR(W_NIC_MIN:W_NIC_MAX,N,IBND)`.

What upstream does:

- After an angle rotation, every `OMESH(NOM)%IL_R(LL,N,IBND)` is re-interpolated to the new angles with the same weights and order as L1248: copy to `IL_R_OLD`, sum ascending from `+0` of `(W*IL_R_OLD(LL,IDX))*RSA(IDX)`, then `/RSA(N)`.
- `IL_R` holds the intensities that the neighbour mesh sent (`IL_S` → `IL_R` in `MESH_EXCHANGE`), lagged by one radiation update.
- The sweep reads them at an `INTERPOLATED_BOUNDARY` wall: the ghost value is the mean of `IL_R(NIC_MIN:NIC_MAX,N,IBND)`, at 4357-4361. This is already copied in `s5rad_sweep_k2.F90` on a flat `ILR(NICR,NRA,NSB)`.

Why no kernel now:

- In AMR mode there are no `OMESH` arrays. Patch-to-patch intensity at box and level boundaries is Role 3's lagged `ILR` table (not a live FillPatch of `IL` into the sweep).
- The per-element arithmetic of 3707-3709 equals that of L1248. Angle rotation re-interpolates `ILR` with the existing `s5rad_ilw_old_copy` / `s5rad_ilw_interp` kernels (§15.3) and an identity gather (`W_BR_INDEX(L)=L`, `NBR=NICR`).
- The loop order (`LL`, then band, then `I_INTP`) does not matter: every element is independent.

#### Role 3 entry point — Radiation review (2026-10-06)

Delivered on FireX `cdb0e8d45d` (`IlGhost.H`/`.cpp`, `notes/radiation-il-ghost.md`, `test_il_ghost`): `fdsrt::build_il_ghost_table` builds index tables only (no intensity, no average).

| Contract item | Evidence | Verdict |
|---|---|---|
| Layout `ILR(NICR,NRA,NSB)` ghost-fastest | notes §Contract; matches `s5rad_wall1` / `s5rad_ils` | **Accept** |
| `W_NIC_MIN`/`W_NIC_MAX` + sender lists in FDS `IIO_R` order | positions grouped by (recv,send) then K,J,I with I fastest; test compares to coordinate `IIO_R` reimplementation | **Accept** |
| Raw lagged fill; mean stays in `s5rad_wall1` | table fills no value; test packs as `s5rad_ils`, copies `ILS→ILR`, checks wall1 mean == FDS mean bitwise | **Accept** |
| Lag = whole-array `ILS→ILR` (not pointer swap) | notes: unswept angles under `ANGLE_INCREMENT>1` keep old values | **Accept** |
| Equal-level bitwise; c/f ratios 2 and 4; three-level | `test_il_ghost`: 7370 checks, 0 failed; wrong-layer mutant 660 failures | **Accept** |
| Multi-rank send/recv | `fdsrt::exchange_il_ghost` @ `b3e1a08275`; `regrid_transport_il_ghost` 7521/0 + `regrid_transport_il_ghost_np2` 7522/0 bitwise | **Green (Role 3)** |

**Review verdict: accept** the entry point and host test. Multi-rank send/recv delivered by Role 3 at `b3e1a08275` (below). No Radiation kernel change; host call contract + stub in `s5rad_sweep_host.F90`.

#### Multi-rank ILS/ILR distribution choice (Radiation → Role 3) — **implemented green**

**Choice: global NICR indices, pair-contiguous one-way slabs** (published earlier; Role 3 landed it).

1. Keep the global position numbering from `build_il_ghost_table` (FDS `IIO_R` order). `W_NIC_MIN`/`W_NIC_MAX` passed to `s5rad_wall1` / `s5rad_ils` are those global 1-based indices on every rank — no remapping, so equal-level bitwise is unchanged.
2. Each radiation rank that owns any box appearing as sender or receiver allocates `ILS(NICR,NRA,NSB)` and `ILR(NICR,NRA,NSB)` (ghost index fastest). Rows this rank never touches may stay undefined; kernels only read `W_NIC_*` ranges of walls on local boxes and write `ILS` rows for local sender packs.
3. Message unit = one `pair_lo[p] .. pair_lo[p+1]-1` range (already one contiguous (receiver box, sender box) block). Payload is the Fortran slab `ILS(lo:hi, 1:NRA, 1:NSB)` (first dimension contiguous ⇒ no gather). Direction: sender-box rank → receiver-box rank into the matching `ILR(lo:hi, :, :)` slab. Same-rank pairs: copy those slabs `ILS→ILR` (or one whole-array copy covering all same-rank pairs).
4. Timing: after each radiation update's pack of swept `(N,IBND)` into `ILS`, exchange/copy the full `(NRA,NSB)` for every pair range (unswept angles keep prior `ILS` values — the FDS lag). Do not pointer-swap.
5. Rejected alternatives: rank-local renumbering of `LL` (breaks bitwise / remaps `W_NIC_*`); per-angle messages (extra sync, worse than one slab per neighbour pair per update).

**Role 3 delivery (`b3e1a08275`):** `fdsrt::exchange_il_ghost(t, box_rank, nra, nsb, ils, ilr)` in `IlGhost.H`/`.cpp`; `pair_recv` on the table; AMReX `ParallelDescriptor` Asend/Arecv for cross-rank pair slabs; same-rank `memcpy` of those slabs. Host checks: `regrid_transport_il_ghost` 7521/0; `regrid_transport_il_ghost_np2` 7522/0 bitwise. Notes: `Source/regrid_transport/notes/radiation-il-ghost.md`.

#### Host call contract (Radiation → Integration / AMReX driver)

Production AMReX radiation host is not yet linked for a live call (still blocked on Integration `rad_*` device driver / gate). Call order is fixed here and mirrored as a stub comment in `s5_rad_sweep/s5rad_sweep_host.F90`. Single-box / FDS-mode sweep bitwise tests are unchanged (they keep the whole-array `ILS→ILR` lag and do not link IlGhost).

**Entry points (Role 3):**
- Once per regrid: `fdsrt::build_il_ghost_table(boxes, ratio, walls)` → `IlGhostTable t` (`cdb0e8d45d`).
- Per radiation update / K-pass, **after** all `s5rad_ils` packs of that pass and **before** the next pass's `s5rad_wall1` reads `ILR`: `fdsrt::exchange_il_ghost(t, box_rank, NRA, NSB, ILS, ILR)` (`b3e1a08275`).

**Call order:**
```
regrid:     t = build_il_ghost_table(...)
radiation update / K-pass:
  DO IBND / DO N (angle subset):
    ... s5rad_wall1 (reads ILR from prior exchange) ...
    ... plane sweep ...
    s5rad_ils   ! pack ILS(:,N,IBND) only
  ENDDO
  exchange_il_ghost(t, box_rank, NRA, NSB, ILS, ILR)
  ! replaces whole-array ILS→ILR lag; same-rank slabs memcpy; cross-rank pair slabs MPI
  ! do NOT pointer-swap; do NOT re-mean (mean stays in s5rad_wall1)
```

Radiation does not duplicate Role 3's mean or remapping. FDS-mode multi-mesh path keeps host `MESH_EXCHANGE` unchanged.

loop_claims: L1251 `not planned (host)`, AMR table+exchange owned by Role 3; Radiation accepted `cdb0e8d45d` + `b3e1a08275`; host call order documented §17.2 + `s5rad_sweep_host` stub.

### 17.3 L1247: angle weights (3623-3645)

Scope confirmed: 3623-3645 is exactly `s5rad_interp_weights(NRA,DLANG,DLANG_OLD,IDX,W)` (host, committed with L1248).

- Every interp test case compares `IDX` and `W` bit for bit with the verbatim reference.
- Seven weights mutants are killed: `no_floor`, `maxloc_back`, `not_normalised`, `dist_linear`, `no_mask`, `cosang_swapped`, `wsum_start`.
- `wsum_start` (the sum starting from the first weight) survived until the `antipodal` group was added. That group has old angles near `-DLANG` with `NRA` 2-8, so the weights are small and the rounding of the `WSUM` start shows.

**Closed: translated (host).** No device kernel: the routine is `O(NRA^2)` once per rotation, and `IDX`/`W` are uploaded with the other angle tables (§17.6).

### 17.4 L1241: CFACE `Q_RAD_IN` zero (3902-3908)

**Disposition: hand-written with a drift check, not the generator.** The generator (`s5gen.py`) has wall policies only (`[policy.wall.*]`) and no CFACE tables or sidecar policy, so a sidecar entry would not generate.

The loop is the CFACE twin of `rad_wall_qin_zero` (L1239).

```
subroutine s5rad_cface_qin_zero(NCF,NSF,CF_B1_PRESENT,CF_SURF_INDEX,SF_TMP_GAS_FRONT,CF_Q_RAD_IN) bind(c)   ! K2, s5_rad_cface/
```

- Copy blocks 3902, 3904 and 3907-3908 (rename `cfqin`); glue 3903 and 3905-3906 (pointer aliases); region 3902-3908.
- In AMR mode the CFACE tables do not exist and CC_IBM is refused, so `NCF=0`. The kernel is ready for when CFACE slots exist (doc 05).
- Tests:
  - 120 cases, 0 vacuous;
  - six flag sets PASS at 1 and 4 threads;
  - 6 of 6 mutants killed: `no_presence_test`, `strict_less`, `negative_zero_store`, `surface_by_icf`, `loop_one_short`, `ge_instead`;
  - drift, K2 CI and lint clean.
- Code is on s5-land (`664917ad50`).

### 17.5 L1246 `GET_KAPPA`, L1252 `KAPPA_WSGG`, L1230-L1233 `A_WSGG`: host

Call sites (all in `radi.f90`):

- `GET_KAPPA` at 4013 (the RADCAL `KAPPA_GAS` cell loop 4005-4016);
- `KAPPA_WSGG` at 4071;
- `A_WSGG` at 3931 and 3949 (band setup), 4068 (WSGG `KAPPA_GAS`/`KFST4_GAS` cell loop 4051-4082) and 4243 (wall `OUTRAD_W` setup).

None of these loops is a device kernel. The translated set (`rad_kernels.toml`: 3887-3893, 4099-4161, 4206-4222, 4818/4820, 4949-4983) and the sweep kernels never call them. So a device callee has no caller.

**Disposition: host**, with these reasons:

1. The callers run on the host. Their outputs (`KAPPA_GAS`, plus `KFST4_GAS` for WSGG) are uploaded per band per radiation update, and the K2 extinction and source kernels read them.
2. Blockers for a future device callee:
   - `GET_KAPPA` calls `GET_MOLECULAR_WEIGHT` / `GET_MASS_FRACTION_ALL`, uses an automatic `Y_OUT(N_SPECIES)`, an `INTENT(INOUT) Z_IN`, and `LOG`;
   - `A_WSGG` is `RECURSIVE`, but only one level deep (`JWSGG==5` sums `MM=1..4`), so it could be unrolled;
   - both use `**REAL(n,EB)` powers, which are libm-class `pow` calls (K2-10, a registry entry and the D-079 tolerance class);
   - the WSGG cell loop also calls `GET_VOLUME_FRACTION` / `GET_MASS_FRACTION` and an ordered `RAD_Q_SUM` reduction.
3. Revisit trigger: a profile of a RADCAL wide-band or WSGG case where the 4005-4016 or 4051-4082 loop, or the per-band upload, shows up. The kernel would then be the whole cell loop with these functions as device callees (generator `[[libm]]` entries), not the functions alone.

loop_claims: L1246, L1252, L1230-L1233 `not planned (host)`.

### 17.6 L1234-L1237: `CALCULATE_DIRECTION_COEFFICIENTS` (3510-3592)

**Disposition: host setup, not planned.**

Call sites:

- `main.f90:383` at initialisation;
- `main.f90:573` and `1114` on angle rotation (`ITER==1`, immediately before `INTERPOLATE_IL`).

Each call is `O(NRA^2)` or less. After each call the driver uploads the angle tables:

- `DLX/DLY/DLZ`, `DLN`, `DLANG`, `DLM`, `DLO_ORIENTATION`, `VIEW_ANGLE_FACTOR`, `NEAREST_RADIATION_ANGLE`;
- and, after L1247, `IDX`/`W`.

It keeps the host copies. `DLANG_OLD` stays on the host as the input of L1247.

### 17.7 Correction of L1238 and strengthening of L1249/L1250 (landed on s5-land)

Two defects in the version first committed with §16 (s5-gen `59ea7dbffa`, then merged to s5-land):

1. `s5rad_hrr_scatter` received a precomputed `DELTA` and `CHI_R_DELTA=CHI_R*DELTA`. Upstream evaluates `CHI_R*TIME_RAMP_FACTOR*HRRPUV*VOLUME_ADJUST` left to right, so the regrouped product is not bitwise.
2. The test reference was a restatement, not the upstream lines, so it could not catch defect 1.

Corrected version (`s5_rad_hrr/`, s5-land `6ccf137752`):

```
subroutine s5rad_hrr_gather(IBAR,JBAR,KBAR,XC,YC,ZC,SHAPE,X1,X2,Y1,Y2,Z1,Z2,X0,Y0,Z0,RADIUS,HEIGHT,NTOUCH,IT,JT,KT)
subroutine s5rad_hrr_scatter(IBAR,JBAR,KBAR,MODE,NTOUCH,IT,JT,KT,TIME_RAMP_FACTOR,HRRPUV,VOLUME_ADJUST,CHI_R,Q,KFST4_GAS)
```

- Copy blocks 5110-5121, 5141-5143 and 5122-5129, under a drift check `check_hrr_drift.py`. The CC_IBM block 5132-5140 is omitted and pinned by hash. Glue is 5105-5109 and 5144 (record loop, `HRRPUV` filter, host ramp).
- The test builds its reference from the verbatim upstream lines (`make_hrr_tests.py`):
  - 150 cases, with `Q` and `KFST4_GAS` compared bit for bit (ghost layer included);
  - coverage: overlapping records, `HRRPUV<=0`, ramps, MODE 0/1/2, cone levels clipped to `RR=0`, bounds exactly on cell centres, and centres exactly on the cone surface;
  - four host flag sets PASS at 1 and 4 threads (host only: no device code, so the offload flag sets do not apply). Re-verified after the land.
- 12 of 12 mutants killed, including `chi_regroup`, `q_regroup`, both strict bounds, `cone_ge` and `reverse_record_order`.

L1249/L1250 were committed with the L1248 suite only. They now have:

- a 150-case suite (new `antipodal` group; particles skipped only by `ORIENTATION_INDEX=0`);
- six flag sets PASS at 1 and 4 threads;
- 31 of 31 mutants killed, 9 of them new: `cf_gather_by_icf`, `cf_inplace`, `cf_band_fixed`, `lp_no_orientation_filter`, `lp_no_class_filter`, `lp_inplace`, and the three weights mutants of §17.3.

### 17.8 K2-09 in discovery mode

`k2_ci_check.py` in discovery mode (`DEFAULT_K2_GLOBS`) does not list `s5_rad_sweep/`. The `[[libm]]` entry for `s5rad_plane3d_ild` (kernel_registry 374a64e) is therefore reported as stale (K2-09 WARN). The sweep test now runs discovery with `s5rad_sweep_k2.F90` added to the globs. It asserts no K2-09 for `s5rad_plane3d_ild` and no FAIL in the file; the K2-10 `EXP` finding is WAIVED by the registry entry. The shared default globs still need `s5_rad_sweep/` (tool owner). Discovery append is on s5-land (`664917ad50`); shared DEFAULT_K2_GLOBS still need the tool owner.

## 18. Setup and GEOM-deferred loop closes (confirmed)

Legacy Mapper claimed these in `amrex/loop_claims.csv` (docs `c4f7c7f` / pin fix `ea8a62d`). Radiation Lead confirms the dispositions below. Pins rechecked at FireX `8ad685f435`. Inventory pins for L0825-L0827 had drifted into particle unpack inside `MESH_EXCHANGE` (`PARTICLE_RECV_BUFFER` / CODE 11 near main.f90:3985-4055); use the FireX lines here, not the inventory CSV.

### 18.1 L0825 / L0826 / L0827: `ALLOCATE_RADIATION_RECV_PKG` / `SEND_PKG` (host)

| Claim | Routine | FireX lines | What |
|---|---|---|---|
| L0825 | `ALLOCATE_RADIATION_RECV_PKG` | main.f90:4092 | `DO I=1,N_REQ5` / `MPI_REQUEST_FREE` |
| L0826 | same | main.f90:4094-4123 (AIC/LL/N nests 4104-4112) | recv buffer size + `MPI_RECV_INIT` |
| L0827 | `ALLOCATE_RADIATION_SEND_PKG` | main.f90:4136-4163 (AIC/LL/N nests 4146-4154) | send buffer size + `MPI_SEND_INIT` |

Coarse span covering both routines: main.f90:4087-4166 (task handout used 4101-4180).

**Disposition: not planned (host).** MPI persistent-request setup and OMESH pack buffers. In AMR mode there is no `OMESH` radiation exchange; box-face intensity exchange is the FR-062 mechanism (`radiation/02` A2; §b face buffers). FDS-mode multi-mesh path keeps these host routines unchanged. No device kernel.

### 18.2 L0676: `PACK_BOUNDARY_RADIA` (host)

FireX: func.f90:5193-5218; band loop `DO NN=1,NUMBER_SPECTRAL_BANDS` at 5208-5211 (then `BR%IL` equate). Inventory pin 5143-5146 is stale.

**Disposition: not planned (host).** Restart / wall `ILW` pack-unpack via `EQUATE` into `STORAGE_TYPE`. Setup and I/O path; not a radiation update kernel. Slot tables for the sweep already use flat `BR_ILW` (doc 05); this pack stays on the host for FDS-mode restart compatibility.

### 18.3 L0664 / L0665: `ALLOCATE_BOUNDARY_RADIA_ARRAYS` (host)

| Claim | FireX lines | What |
|---|---|---|
| L0664 | func.f90:4355-4358 | open-slot search `DO I=NAS,M%N_BOUNDARY_RADIA_DIM` (grow path continues 4362-4371) |
| L0665 | func.f90:4378-4381 | band `ILW` alloc/zero; `BR%IL` init follows at 4383-4384 |

Routine span: func.f90:4347-4386 (handout 4355-4385). Inventory 4351-4354 / 4374-4377 are stale by a few lines.

**Disposition: not planned (host).** Storage grow and per-band init for `BOUNDARY_RADIA`. Host setup only.

### 18.4 L0006 / L0007: `CCCOMPUTE_RADIATION` (GEOM deferred)

| Claim | FireX lines | What |
|---|---|---|
| L0006 | ccib.f90:7139-7143 | no-radiation branch: cut-cell `QR = -CHI_R*Q` |
| L0007 | ccib.f90:7147-7155 | radiation on: copy Cartesian `QR(I,J,K)` onto cut cells |

**Disposition: GEOM / CC_IBM deferred (D-033).** Host until cut-face / cut-cell tables exist for AMR. Same deferral as sweep §11 (`CC_IBM` refused). No kernel now. loop_claims status remains `not planned (host)` with the D-033 note.

### 18.5 Register

loop_claims for L0825-L0827, L0676, L0664, L0665, L0006, L0007: `not planned (host)` confirmed; pins as above; see this section.

## 19. L1242 remaining open pieces (classify / feature request only)

No new `s5-gen` commits. Land pause until Legacy Mapper posts FF (s5-land at `664917ad50`). Classification of the still-open A1 rows in `radiation/02` (merged line pins; refine against FireX when a piece is taken):

| Piece (radiation/02 A1) | Lines (merged) | Classification | Feature / owner request |
|---|---|---|---|
| Band-init whole-array fills (`KFST4_*`, `SCAEFF*`, `KAPPA_PART`) | 3914-3920 | **translatable now (fills)**; not started | none; trivial whole-array / fill rewrite already used elsewhere; do after FF open, not before |
| Liquid-droplet `INTERPOLATE1D` per particle class | 3959-3972 | **generator feature** | assumed-shape rank-1 table callee (same blocker as L0658); request generator support or keep host per-class loop |
| Solid-particle deposit scatter | 3981-3991 | **design (FR-062 / FR-005 (4a))** | deterministic deposition; tables landed with L1240/L1244 (§16.3); deposit scatter still design |
| RADCAL `KAPPA_GAS` / `GET_KAPPA` | 4002-4019 | **generator feature** (aligned §17.5 host for now) | source-function callees + module tables; revisit with wide-band profile |
| Wide-band `BLACKBODY_FRACTION` / `INTERPOLATE1D_UNIFORM` sums | 4027-4046 | **generator feature** | source-function + ordered scalar reduction (D-053) |
| WSGG cell loop / `A_WSGG` / `KAPPA_WSGG` | 4051-4082 | **host now** (§17.5); device later needs features | recursion removal / unroll; `[[libm]]` pow; ordered reduction |
| Condensed-species absorption | 4171-4202 | **generator feature** (same family as droplet row) | `INTERPOLATE1D` on `SPECIES_MIXTURE%R50` / `WQABS`/`WQSCA`; `CC_IBM` skip; per-species host loop acceptable until assumed-shape callee exists; was "open" in radiation/02 |
| Wall `OUTRAD` / `INRAD_W` | 4229-4258, 4262-4280 | **generator feature** (wall tables) | ragged `BR_ILW` / `BR` alias (R1 / doc 05); partial coverage already in sweep wall kernels |
| Wall ghost / `CELL_ILW` (`WALL_LOOP1`, CFACE twin) | 4328-4384 | **FR-062 + wall tables** | largely covered by sweep `s5rad_wall1` / `cell_ilw_*`; CFACE twin waits on geometry tables |
| Post-sweep wall/CFACE `INRAD` updates | 4749-4812 | **generator feature** / FR-062 | wall-table features; cylindrical copy stays out (D-073 (5)) |

**Not in this table (already closed or owned elsewhere):** angular sweep and `UIID`/`QR`/`UII` (FR-062 / translated); `IL_S` mesh pack (FR-062 exchange, §18.1); L1238/L1247-L1250/L1240/L1244 translated; L1245 host; L1251 Role 3 ghost fill.

**Next-work status (2026-10-06):** (1) gate hookup — partial / blocked on Integration device `rad_*` driver + `ported.toml`; host steps in `radiation/03` §10. (2) ghost `UII`/`UIID` reader list — **closed** (V&V co-sign; R-D073b-out (A)). (3) L1242 A1 remainders: feature requests only until taken. (4) solid/08 R1–R4 — **closed** on Radiation (`05` + §14). (5) L1251 — **multi-rank green**: `build_il_ghost_table` (`cdb0e8d45d`) + `exchange_il_ghost` (`b3e1a08275`); host call order in §17.2 + `s5rad_sweep_host` stub; live AMReX link waits Integration driver. (6) FR-063 regrid — unblocked in principle by L1251 green; short draft in `radiation/04`; live work still gated on device `rad_*` driver.
