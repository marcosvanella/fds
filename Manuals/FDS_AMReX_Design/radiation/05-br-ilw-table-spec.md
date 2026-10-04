# BR_ILW table: spec for the generator and wall-table work

Owner: AMR Radiation Lead. Audience: GPU Generator Engineer (policy and builder), GPU Wall Loops Engineer (wall tables and driver fill).

Line numbers are `Source/radi.f90`, `func.f90`, `type.f90` at `bee11f0329` (the s5-gen reference). Companion docs: `00-r1-signoff.md` (why L1243 and L1248 need this table), `03-radiation-translation-notes.md` §5 (the refusals), `06-fr062-sweep-kernel-design.md` (the sweep that updates it).

## 1. What the data is today

- `BOUNDARY_RADIA(BR_INDEX)` is one record per wall, CFACE or Lagrangian particle (type.f90:381-384, comment "Angular radiation intensities associated with a WALL, CFACE, or LAGRANGIAN_PARTICLE"). It holds `IL(1:NSB)` (output only) and `BAND(1:NSB)%ILW(1:NRA)` (type.f90:183-185).
- **All records have the same shape.** `NSB = NUMBER_SPECTRAL_BANDS` and `NRA = NUMBER_RADIATION_ANGLES` are global run constants, set during input and `INIT_RADIATION` (radi.f90:2788-); every record is allocated with them (func.f90:4373-4380; `ILW=0`, `IL=SIGMA*TMPA4/PI`). So the table is **dense** in (NRA, NSB). "Ragged" in the earlier notes meant "array of allocatable components, not a rectangular array", not variable lengths. An offset/count scheme is not needed; see §4 for what to do if that ever changes.
- **One slot space per mesh.** Walls, CFACEs and particles take slots from the same list: `BOUNDARY_RADIA_OCCUPANCY`, `NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT`, `N_BOUNDARY_RADIA_DIM` (func.f90:4343-4382). Slots are unique per owner and are reused when a particle is deleted. The list grows by `N_NEW_STORAGE_SLOTS` by reallocation (func.f90:4364-4365). The index is stored on the owner: `WALL%BR_INDEX` (type.f90:395), the particle (type.f90:441), `CFACE%BR_INDEX` (type.f90:1367).
- Every SURFACE has `INCLUDE_BOUNDARY_RADIA_TYPE=.TRUE.` (type.f90:1001), so every wall, including `NULL_BOUNDARY` walls, has a slot. Particle classes default to `.FALSE.` (type.f90:164).
- Sizes (estimates, 8-byte reals): `NBR*NRA*NSB*8` bytes. A 32^3 box has about 6.1e3 boundary walls with no obstructions, so about 4.9 MB at NRA=100, NSB=1, and about 29 MB at NSB=6. Particles with a record add 800 B each per band: 1e5 radiating particles are about 80 MB per band.

## 2. One table family or three

**Decision: one table, indexed by the slot `BR_INDEX`, plus one gather array per owner kind.**

- `BR_ILW(1:NBR, 1:NRA, 1:NSB)` real, with `NBR = N_BOUNDARY_RADIA_DIM`.
- `BR_IL(1:NBR, 1:NSB)` real (the output-only `IL`; only particles with a massless target write it, radi.f90:4875 in the particle loop).
- Gather arrays: `W_BR_INDEX(NWE+NWI)` (int), `CF_BR_INDEX(NCF)` (int, geometry-deferred), `LP_BR_INDEX(NLP)` (int, only when the particle kernels exist).

Reasons:
1. The slot space is shared, so three tables would need three occupancy lists and the host would have to translate slots on every growth and every particle deletion. One table keeps the host logic of func.f90:4343-4382 as it is.
2. Slots are unique, so a launch over walls (or over particles) never has two threads writing the same `ILW(:, N, IBND)` element: no race, no atomics.
3. A table indexed by wall (`IW`) would need a gather and scatter against the slot table every time step, or a second master copy. The slot table is the master on the device, and wall, CFACE and particle kernels address it through their own gather array.

The existing wall policy gathers `WALL`, `BOUNDARY_PROP1` etc. into tables indexed by `IW` (`[policy.wall]`, s5_markers.toml:44-50). `BR_ILW` is the one record type that does not follow it, because it is shared between owners and persistent. That is the feature the generator needs (§6).

## 3. Layout: slot-fastest

**Decision: `BR_ILW(slot, angle, band)`, slot index first (stride 1).**

Access patterns, with thread = wall (or particle):

| Loop | Pattern | Slot-first | Angle-first `(angle, band, slot)` |
|---|---|---|---|
| Sweep boundary updates WALL_LOOP1/2/3, CFACE loops, particle loop (radi.f90:4328-4371, 4749-4812, 4848-4881): one launch per angle `N`, band `IBND`; thread per wall reads/writes `ILW(N)` | threads in a warp read consecutive slots at fixed `(N, IBND)` | **coalesced** (walls take slots in wall order at set-up, so adjacent walls have adjacent slots; to be checked by the host assertion in §5) | stride `NRA*NSB*8` bytes between threads: uncoalesced |
| L1243 (4965-4974): thread per wall, `DO IBND; DO N` sums | at each `(N, IBND)` step a warp reads adjacent slots | **coalesced** | each thread walks its own contiguous `NRA` run: good for the thread, poor across the warp |
| L1248 (3651-3662), once per angle cycle | thread per (slot, band, angle) | coalesced on slot | coalesced on angle |
| MPI pack of particle `ILW` (func.f90:5128-5153), restart | whole record | strided gather (rare) | contiguous |

The hot loops are the per-angle boundary updates: they run `NRA/ANGLE_INCREMENT` times per update per band. L1243 and L1248 run once per time step or per angle cycle and are tiny (0.000% and 0.001% in the survey). So the layout is chosen for the sweep boundary updates. L1243's sum over angles is a serial loop per wall in either layout, and its **order** is the same: ascending `N` inside a per-band temporary (rule in §7).

If a later measurement shows that a kernel does better with angle-first, the generator interface does not change: it is only the order of the extents in the declaration.

## 4. Offset/count scheme

Not needed while NRA and NSB are global. State it in the builder as an assertion: every record has `SIZE(BAND)==NSB` and `SIZE(BAND(n)%ILW)==NRA`.

If a future change makes either per-record (for example per-particle-class band sets), use a CSR layout: `BR_OFF(1:NBR+1)` int, `BR_NA(1:NBR)`, `BR_NB(1:NBR)`, and data `BR_ILW_DATA(BR_OFF(slot) + (IBND-1)*BR_NA(slot) + N)`. The sum rule and the ascending order stay the same. Not recommended now: it costs an extra read per access and loses coalescing.

## 5. Host builder, ownership, write-back

**Ownership.** On the device, `BR_ILW` is the master for the whole radiation update. The host `BOUNDARY_RADIA` records are filled at set-up and after a regrid, and read back only when the host needs them (restart dump, particle MPI transfer, output). There is no per-step gather and scatter.

**Fill (host to table).**
1. `NBR = N_BOUNDARY_RADIA_DIM`. Allocate `BR_ILW(NBR,NRA,NSB)`, `BR_IL(NBR,NSB)`.
2. For every occupied slot (`BOUNDARY_RADIA_OCCUPANCY(slot)==1`) copy `BAND(n)%ILW(:)` and `IL(n)`. Unoccupied slots: zero (never read).
3. Build `W_BR_INDEX(IW) = WALL(IW)%BR_INDEX` for `IW=1..NWE+NWI`.
4. **Assertion (host, once per build and after regrid or growth):** for every wall with `BOUNDARY_TYPE/=NULL_BOUNDARY` (and in fact for all walls), `0 < BR_INDEX <= N_BOUNDARY_RADIA_DIM`, the slot is occupied, and the slot is unique among walls. No device guard is needed (00-r1-signoff.md; unguarded uses at radi.f90:4266, 4308, 4332 rely on the same fact).
5. Growth of `N_BOUNDARY_RADIA_DIM` (a particle appears) reallocates the table: copy old to new on the device, zero the new slots.

**Who writes `ILW` in a time step** (all disjoint by slot, so no race within one launch):

| Writer | Where | What |
|---|---|---|
| Zeroing | radi.f90:4305-4310 | open walls: `ILW(ANGLE_INC_COUNTER)=0` per update |
| WALL_LOOP1 | 4328-4371 | mirror: `ILW(N)=ILW(DLM(N,\|IOR\|))` (reads **another angle** of the same wall); solid: `ILW(N)=OUTRAD_W+RPI*(1-EMIS)*INRAD_W`; also sets `CELL_ILW` |
| CFACE_LOOP1 | 4375-4385 | `ILW(N)` from `OUTRAD_F`, `INRAD_F` (geometry-deferred) |
| cylindrical copy | 4722-4742 | `BR_UP%ILW(N-1)` or `(N)` from `BR_DOWN%ILW(N)` (reads another record; two records of different walls) |
| WALL_LOOP2 | 4749-4775 | `INRAD_W` step 1; `ILW(N)=IL(IIG,JJG,KKG)` (or `ILD*` for `RAD_DIFF_SCHEME>1`); `INRAD_W` step 2 |
| WALL_LOOP3 | 4779-4787 | open walls: `ILW(ANGLE_INC_COUNTER) -= DLN*IL` (**accumulation over the angles N**, ascending in the processing order) |
| CFACE_LOOP2/3, particle loop | 4791-4812, 4848-4881 | the same for CFACEs and particles with an ORIENTATION |
| INTERPOLATE_IL | 3651-3695 | rewrites every angle of every record at a new angle cycle |
| Allocation, regrid, particle insert | func.f90:4343-4382 | zero fill (replace by the fill rule of 04 §1.4 after a regrid) |

`INRAD_W(IW)` and `OUTRAD_W(IW)` (pointers to `WALL_WORK2`, `WALL_WORK1`, radi.f90:3841-3842) are plain per-wall mesh arrays indexed by `IW` with extent `NWE+NWI`. They fit the existing `[policy.wall.arrays]` as they are (the entries `UVW_SAVE`, `U_GHOST` ... are of the same kind). `INRAD_W` is updated per angle in angle order: that is a fixed-order accumulation (§7, and 06 for the sweep).

**Write-back (table to host).** At restart dump and before a particle leaves the rank or the mesh: scatter `BR_ILW(slot,:,:)` back into `BOUNDARY_RADIA(slot)%BAND(n)%ILW`. This is a one-record gather per event, not per step. It is also the path the AMR code uses to move wall records at regrid (04 §2).

## 6. Interaction with the generator

Current state (verified in the worktree): `[policy.wall.arrays]` accepts rank-1 mesh arrays indexed by `IW`; a rank-1 `BR_ILW` in a wall loop is refused ("not declared in `[policy.wall.arrays]`"); a rank-3 `BR_ILW` is silently mapped to a mesh-field shape, which is wrong here (03 §5). The alias `BR => BOUNDARY_RADIA(WC%BR_INDEX)` is refused as "not one of the recognised wall aliases".

Needed, in order:

1. **Alias recognition.** `BR => BOUNDARY_RADIA(WC%BR_INDEX)`, where `WC => WALL(IW)`, is a gather alias like `B1 => BOUNDARY_PROP1(WC%B1_INDEX)`. The generator replaces `BR%BAND(IBND)%ILW(N)` by `BR_ILW(W_BR_INDEX(IW), N, IBND)` and `BR%IL(IBND)` by `BR_IL(W_BR_INDEX(IW), IBND)`. `W_BR_INDEX` is a normal flat wall table (prefix `W_`, component `BR_INDEX`, int, extent `NWE+NWI`): no new feature beyond the table name being known.
2. **A slot-table declaration** in the sidecar, for example
   `[policy.wall.slot_tables]` `BR_ILW = { shape = ["NBR","NRA","NSB"], gather = "W_BR_INDEX", component = "BAND(*)%ILW" }`.
   The kernel gets `BR_ILW` as a rank-3 dummy with explicit extents; the golden signature records `BR_ILW:slot(NBR,NRA,NSB):real:inout`.
3. **Several records in one loop** (mirror boundary, cylindrical copy): two gather aliases in a loop. Do this only after (1) and (2) work for single-record loops.
4. **Ranges that reach another angle** (`ILW(DLM(N,ABS(IOR)))`): index expressions on the angle extent are fine, but the generator must not claim the loop is independent across angles. Mark such kernels "sequential in angle" (06 §b).

Not needed: derived-type pointer records on the device; `ALLOCATE(MOLD=)`; CSR offsets.

## 7. Rule for L1243 (per-band temporary)

```
Q = 0
DO IBND = 1, NSB
   T = 0                                  ! +0, not the first element
   DO N = 1, NRA
      T = T + BR_ILW(W_BR_INDEX(IW), N, IBND)
   ENDDO
   Q = Q + T
ENDDO
B1_Q_RAD_IN(IW) = Q
```

Conditions kept from the original (radi.f90:4965-4974): `IW=1,N_EXTERNAL_WALL_CELLS` only, `BOUNDARY_TYPE==OPEN_BOUNDARY` only. `B1_INDEX==0` cannot occur for an open wall; the original has no guard, so the kernel has none, but the flag `B1_PRESENT` (below) is read for L1239. A flat chain over bands and angles is not bitwise equal for NSB>1 (tested: 51 mismatches). Compile with FMA contraction off for these kernels if bitwise equality against gfortran is claimed.

### Additional rules for L1243 and B1_PRESENT

- **Ownership of Q_RAD_IN.** Covered coarse OPEN walls are not in the owner list (row ownership): a wall row covered by a finer box has no owner row of its own and gets no `Q_RAD_IN` from the kernel.
- **Summation order.** The L1243 bitwise claim holds for gfortran's `SUM` order only (ascending, one accumulator per band). ifx may vectorise `SUM` and change the order; for ifx the claim becomes "within rounding" unless the kernel's explicit chain is the reference (the kernel's chain is the one written in the rule above).
- **B1_PRESENT is the pure flag.** `W_B1_PRESENT(IW)` is exactly `B1_INDEX /= 0`; it does not include the `NULL_BOUNDARY` test. The L1243 kernel keeps its own `NULL_BOUNDARY` test (it does not take it from `B1_PRESENT`), and L1239 combines both as the original does (`B1_INDEX==0 .OR. BOUNDARY_TYPE==NULL_BOUNDARY`).
- **Q_RAD_IN belongs to the owner row, conditional on EMISSIVITY staying per surface** (func.f90:4920). If `EMISSIVITY` ever becomes per cell or per band-dependent state outside the surface record, the owner-row rule for `Q_RAD_IN` has to be re-examined.

## 8. B1_PRESENT

`W_B1_PRESENT(IW) = 1` if `WALL(IW)%B1_INDEX /= 0`, else 0 (int, extent `NWE+NWI`). Used by L1239 (radi.f90:3887-3893), where the original test is `B1_INDEX==0 .OR. BOUNDARY_TYPE==NULL_BOUNDARY` then `SF%TMP_GAS_FRONT<=0 -> B1%Q_RAD_IN = 0`. Semantics: a gather index of 0 means "no `BOUNDARY_PROP1` record", and the flat tables `B1_*` are then undefined for that wall, so the kernel must not read them. For all walls with `B1_PRESENT=0` the kernel does nothing. The generator now reads `SURFACE(WALL%SURF_INDEX)%TMP_GAS_FRONT` as the per-surface table `SF_TMP_GAS_FRONT(0:N_SURF+N_SURF_RESERVED)` indexed through the wall table `W_SURF_INDEX(IW)` (generated interface: `IBAR,JBAR,KBAR,NWI,NWE,N_SURF,N_SURF_RESERVED,B1_Q_RAD_IN,SF_TMP_GAS_FRONT,W_B1_PRESENT,W_BOUNDARY_TYPE,W_SURF_INDEX`); an earlier version of this note and of the test assumed a per-wall gather `SF_TMP_GAS_FRONT(IW)`. The generated kernel `rad_wall_qin_zero` already uses `W_B1_PRESENT`; the sidecar injects the component through `type_comps` until the wall table has it.

## 9. L1248 (INTERPOLATE_IL) and the ping-pong

Original (3651-3662): per wall, per band: `ILW_OLD = ILW` (automatic copy of the `NRA` values, radi.f90:3613, 3655); `ILW = 0`; for `I_INTP=1..N_INTP`: `ILW(:) = ILW(:) + W(:,I_INTP)*ILW_OLD(IDX(:,I_INTP))*RSA(IDX(:,I_INTP))` (product grouped `(W*ILW_OLD)*RSA`); then `ILW(:) = ILW(:)/RSA(:)`.

- The copy **must stay**: angle `N` reads the old values of other angles. In place, thread order would change the result and race.
- Device form: before the call, copy `BR_ILW` to `BR_ILW_OLD` (a device-to-device copy of the whole table; same shape; allocated only for the call), then a kernel with one thread per (slot, band, angle) reads `BR_ILW_OLD` and writes `BR_ILW`. The `I_INTP` sum is a serial ascending loop per (slot, band, angle), starting at +0, so the bits equal the original. `IDX`, `W` are the host-built `NRA x NINTP` tables (loop at 3623-3645, 0.3% of a tiny share).
- Alternative: per-thread scratch of `NRA` reals with thread = (slot, band). 100 reals per thread is large for registers; prefer the ping-pong.
- It runs once per angle cycle (main.f90:560, 1101) under `ALLOW_RANDOM_RADIATION_ROTATION`, so host-first is acceptable and it is low priority.
- The CFACE (3665-3676) and particle (3679-3695) loops are the same pattern on other gathers; the `IL_R` loop (3698-3713) is on neighbour-mesh buffers, not on `BR_ILW`.

## 10. Acceptance tests

1. **`open_qin_spec`** (exists in `test/make_rad_tests.py`, in the committed worktree): the verbatim L1243 loop against the specification kernel, 1 to 6 bands, 3 to 104 angles, `-DRAD_MUT_FLATSUM` must fail. The generated kernel replaces the specification kernel in the same test, with the same data, and must pass in the six flag sets at 1/4/8 threads.
2. **Layout test:** the same data stored in a slot-first table and read through the gather must give the same bits as the record form. Mutants: slot and angle swapped, gather by `IW` instead of `W_BR_INDEX`.
3. **Sweep boundary updates:** WALL_LOOP1/2/3 against the verbatim loops for one angle, with solid, open, mirror and interpolated walls; accumulations compared bit for bit over a full angle sequence (the order of `N`).
4. **L1248:** verbatim loop against the ping-pong kernel; an in-place mutant must fail.
5. **Assertion test:** a record with `BR_INDEX=0`, one beyond the dimension, and a duplicated slot must be refused by the host builder.
6. **Growth test:** table reallocation after a slot is added preserves all old values and zeros the new ones.

## 11. Open questions

1. Generator Engineer: is a new `[policy.wall.slot_tables]` section acceptable, or should `BR_ILW` be modelled as an IW-indexed table with a gather and scatter per step (simpler for the generator, costs a copy per step and loses sharing with particles)?
2. Generator Engineer: can the alias `BR => BOUNDARY_RADIA(WC%BR_INDEX)` be handled as a second gather alias in the same loop (needed for the mirror boundary and the cylindrical copy), or only one alias per loop?
3. Wall Loops Engineer: do walls get consecutive slots in wall order after a regrid or an obstruction change? If not, the coalescing argument in §3 weakens (the correctness does not). Can the driver renumber slots at table build?
4. Wall Loops Engineer: who owns `W_BR_INDEX` and `W_B1_PRESENT` in the driver, and is the `NULL_BOUNDARY` wall slot kept (it is allocated today)?
5. Both: should `INRAD_W`/`OUTRAD_W` be declared in `[policy.wall.arrays]` now (extent `NWE+NWI`)? They are needed by every sweep boundary kernel.
6. Both: particles. Does the first version carry particle records in the table at all? Proposal: no. Allocate walls (and CFACEs later) only; keep the particle `BR` on the host until the deposition and particle design are settled.
7. Generator Engineer: kernels whose loop reads `ILW` of another angle (mirror) need a "sequential in angle" tag. Is there an existing tag for loop-carried dependence across launches, or should the host driver enforce the order?
