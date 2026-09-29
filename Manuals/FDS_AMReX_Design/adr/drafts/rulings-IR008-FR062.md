# Rulings: IR-008 refinement outside the refinable region, FR-062 (d) radiation sweep order

Owner: Chief Architect. Status: IR-008 ruling final. FR-062 ruling final on the architecture side; requirement wording to be confirmed by the AMR Radiation Lead before the Spec & Program Lead applies it.

## 1. IR-008: may tagging refine outside the declared refinable region?

**Ruling: no.** In AMR mode every grid at level 1 and above lies inside the refinable region (the union of the finer `&MESH` static boxes and the declared region boxes, IR-008). Tags outside the region are discarded.

Mechanism (AMReX Integration Lead confirms at implementation):
- R1. Region boxes are aligned to the blocking factor at every level they permit (IR-008 already requires this for output tiles).
- R2. Tags outside the region at the level being refined are cleared in `ManualTagsPlacement`, which AMReX calls after tag buffering (`AMReX_AmrMesh.cpp:697,727`), so buffer cells cannot leak out.
- R3. Clustering can still produce boxes that include untagged cells outside the region. New grids are therefore intersected with the region before they are accepted; proper nesting is checked after that step (a setup/regrid assertion in debug builds).
- R4. Diagnostics: when tags are discarded, the run reports once per level the count of discarded tags at the first regrid where it happens, and a total at the end, so the user knows to enlarge the region. No per-regrid log spam.
- R5. An empty region means no refinement beyond level 0; AMR mode then warns once (IR-008 already warns for output).

Rationale:
- ADR-004 D1 writes output at `L_out` only inside the region; refinement outside it would be computed and then invisible in Smokeview output.
- FR-041b ruling condition A2 bounds finest-ever wall records to the region; G2a's memory bound assumes it.
- The output size report (FR-079) and memory planning on small machines stay predictable.
- The owner's §8 answer rules out any output that follows refinement, so the region cannot grow at run time either.

Rejected alternatives:
- Unrestricted tagging, with level-0 output outside the region. Refined results would be hidden from output, A2 and the G2a memory bound would no longer hold, and memory use would depend on the flow.
- Growing the region at run time to follow tags. That changes the output mesh count within a run, which Smokeview cannot load (owner, ADR-004 §8).

## 2. FR-062 (d): upwind sweep order for boxes on the same rank

**Ruling: (d) is not adopted.** The default and only required mode is the fully parallel per-box lagged sweep (S-B, D-039), with `RADIATION_ITERATIONS`. This matches the owner's concern that ordered sweeps where boxes wait on each other will be slow.

Proposed requirement wording for FR-062 (replaces additions (a)-(d) and the open point on (d)):

> (a) Each radiation pass sweeps every box of a level independently, and in any order or concurrently, using box-face intensities from the previous exchange only. Box faces are exchanged after each pass by one uniform mechanism, whether the neighbouring box is on the same rank or another; a box never reads a face value written in the current pass (double-buffered faces). `RADIATION_ITERATIONS=K` gives K passes per radiation step, with an exchange after each pass; defaults `RADIATION_ITERATIONS=1` and `INITIAL_RADIATION_ITERATIONS=3`, as in FDS.
> (b) After a regrid or rebalance, new box faces are filled from upwind or coarse-level intensities, never left at zero.
> (c) Radiation results at a fixed box layout are independent of rank count, thread count, box-to-rank assignment and load balancing (FR-005 (iv); D-039 note in FR-005). They depend on the box split only, within the tolerance class set by the verification below.
> (d) An ordered (upwind, dependency-ordered) sweep across boxes is an optional mode, added only if the lag-error check (V&V item (e)) shows the lagged sweep with K ≤ 3 misses the agreed tolerance on the anchor cases. If added, its order is defined on the level's global box graph for each angle, never per rank, so it keeps (c); per-rank ordering is not allowed.
> (e) V&V: lag-error tolerance against the single-mesh run on `radiation_gas_panel` split into several boxes, with K = 1, 2 and 3 (A-50), plus byte identity across 1, 2 and 4 ranks at a fixed box layout.

Rationale:
- Per-rank upwind order makes the lag depend on which boxes share a rank, so results would change with rank count and with load balancing at a fixed layout. That breaks D-039's rank-independence note and FR-076's rank-independent output.
- The independent sweep has no waiting between boxes, so it maps directly onto GPU execution (one kernel launch over all boxes) and onto AMReX's `FillBoundary`-style exchange.
- Double-buffered faces are what make (c) hold on a single rank: without them, the result would depend on the order in which boxes on one rank are processed.

Rejected alternatives:
- (d) as proposed, per-rank upwind order. Rank-count dependent; conflicts with D-039 and FR-076.
- Global ordered sweep (S-A) as default. Boxes wait on upstream boxes for every angle, which serialises work across ranks and on the GPU; the owner rejected this as the default on speed grounds.

### 2.1 Agreed wording (merged with the Radiation Lead's draft; supersedes the block above)

> (1) Default: a fully parallel per-box sweep. Every box sweeps its angle subset independently using lagged box-face intensities, as FDS meshes do. There is no ordering between boxes, on the same rank or across ranks. A box reads only face intensities from the previous exchange, never values written in the current pass (double-buffered faces), and same-rank and cross-rank faces are exchanged by one mechanism. At a fixed box layout, results do not depend on rank count, ranks per GPU, box ownership or processing order. D-039 exempts only dependence on the box split.
> (2) `RADIATION_ITERATIONS` keeps FDS semantics: K passes per step, with a box-face exchange after each pass. Defaults are `RADIATION_ITERATIONS=1` and `INITIAL_RADIATION_ITERATIONS=3`.
> (3) After a regrid or rebalance, new faces are filled from upwind interior values or from the coarse level, never left at zero.
> (4) V&V: the lag-error check runs `radiation_gas_panel` split into 16^3 and 32^3 boxes at K = 1, 2, 3 against the single-mesh result and sets the recommended K (A-50). At a fixed box layout, radiation outputs are byte-identical across 1, 2 and 4 ranks and run to run.
> (4a) Byte identity across ranks and threads also needs deterministic radiation source inputs. Particle absorption deposition (droplets, vegetation) is deterministic in FireX today (no OpenMP atomics in `part.f90`, `radi.f90` or `vege.f90`) and shall stay deterministic in the GPU port: no order-varying floating-point atomics; per-box accumulation in a fixed order or the D-028 fixed-point sum. The one thread-order-dependent sum in `radi.f90` (the `!$OMP CRITICAL` accumulation of `RAD_Q_SUM`/`KFST4_SUM`, `radi.f90:4119-4122`) is already covered by FR-005 (ii). The byte-identity check in (4) therefore applies to particle cases too.
> (5) If K = 1..3 fails the tolerance, escalate in this order: first, an optional radiation BoxArray with larger boxes, built from the level's geometry and a fixed box size only (never from rank count or ownership) and exchanged as in (1); second, an optional dependency-ordered sweep, which must use the global box order per angle and never a rank-local one.

Cost basis (Radiation Lead estimate): total work is the same for any ordering; the dependency chain per angle is 46 planes for a 16^3 box at K=1 regardless of rank count, against about 190 for an ordered sweep over a 4^3 arrangement of such boxes.
