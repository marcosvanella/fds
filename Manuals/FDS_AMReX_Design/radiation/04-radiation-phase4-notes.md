# Radiation implementation notes: FR-062 sweep and regrid

Owner: AMR Radiation Lead. Companion to `01-radiation-amr-spec.md` (design), `02-radiation-gpu-candidates.md` (loop table) and `03-radiation-translation-notes.md` (generator work).

Sources: `docs/adr/drafts/rulings-IR008-FR062.md` §2.1 (agreed text), requirements FR-062 rules (1) to (5) and (4a), `docs/next-work.md` (Radiation Lead row, items 2 and 3).

Naming. The next-work row calls these "Phase 4 notes". In the roadmap, radiation on levels is Phase 8 (M8). The per-box sweep and its exchange are needed earlier, at M2, for box-split uniform runs (FR-061, FR-002), and are then carried unchanged onto levels. The notes below cover both.

Pins. Citations are to the local reference tree at `bee11f0329` (FireX, upstream `afb5e31a48` merged), the same tree s5-gen uses. `radi.f90` lines are the post-RTE_SOURCE-merge numbers. `main.f90` is 73 lines longer than at the older pin `36975d765f` from `~1000` on, so older notes that cite `main.f90:1022` mean `:1095` here; `01-radiation-amr-spec.md` has been re-pinned to match. `mesh.f90`, `dump.f90`, `func.f90`, `type.f90`, `init.f90` are unchanged between the two pins.

---

## Part 1. FR-062 implementation notes

### 1.1 What is built (rule 1)

One sweep kernel per (box, band, angle). A box is one FAB of the level's radiation BoxArray. By default that is the level's normal BoxArray; see 1.6 for the escalation.

- **No ordering between boxes.** Every box sweeps from its own interior and from its face buffers. No box waits for another box on the same pass. Results do not depend on rank count, box ownership or processing order.
- **Face buffers are double-buffered.** Each box face has an "incoming" buffer (read by the sweep) and an "outgoing" buffer (written by the sweep). They swap only at the exchange. A box can never see a neighbour's value from the current pass. That is what removes the dependence on processing order.
- **One exchange mechanism for every face.** Same-rank faces, cross-rank faces and coarse/fine faces all go through the same pack, copy or MPI move, and unpack. This replaces the FDS split between same-rank copy and `MPI_ISEND` of `IL_S`. Do not add a same-rank shortcut that reads the neighbour's live array; it would break the rule above.
- **Inside a box the sweep is exact.** The upwind recurrence within a box (the 3D slice or wavefront order, `radi.f90` ~4290-4830) is unchanged and is bitwise the single-mesh sweep for that box. The only approximation is across box faces, and it is the lag.

### 1.2 K passes per step (rule 2)

- Defaults follow FDS: `RADIATION_ITERATIONS=1` (`cons.f90:509`), `INITIAL_RADIATION_ITERATIONS=3` (`read.f90:10224`).
- Pass structure per radiation update: `DO ITER=1,K` { sweep all boxes; exchange faces }. FDS has the same shape: `main.f90:555` (initial) and `:1095` (steady), with `MESH_EXCHANGE(2)` after each pass only when K>1 (`main.f90:1110-1116`). With K=1 the exchange runs once at the end of the step (`main.f90:1183-1185`).
- The exchange is after the sweep, so the face values a box reads on pass 1 are the ones left by the previous update. That is the lag. With K=1 the lag is one radiation update; with K=3 the faces have been refreshed twice before the step ends.
- FDS exchanges only once per pass after the first cycle (`IF (ICYC>1) EXIT`, `main.f90:1112-1114`, similarly `:575`). Keep that: one exchange per pass, not one per angle increment.
- `RAD_CALL_COUNTER` advances once per update, on the last pass (`radi.f90:3866`). The angle subset and the random rotation depend on it, so it must be one global value, not per box or per rank.

### 1.3 Source terms and the lagged quantity

- The per-band source `RTE_SOURCE` (`radi.f90` ~4220-4222) is angle-independent and box-local. Compute it once per band per box before the K passes' angle loops, from `KFST4_GAS`, `KFST4_PART`, `SCAEFF`, `SCAEFF_G` and `UIIOLD`.
- `UIIOLD` comes from the previous update. It is the lagged quantity, and it is already lagged in FDS. The box split adds nothing to it. `UIID` accumulates within the update (`radi.f90` ~4818-4820) and is never read across a face.
- On the device the `UIID` accumulation should be fused into the sweep kernel (it runs once per angle). See `02-radiation-gpu-candidates.md`.

### 1.4 Regrid and new faces (rule 3)

- After a regrid, any face with no valid history (a new box face, or a face that was previously interior or on another level) is filled before the first sweep. Never from zero.
- Fill source, in order:
  1. The upwind interior value of the same box and band: for each angle, the cell on the upwind side of that face. This uses data already present.
  2. If the face is a coarse/fine face, or the upwind interior is itself new, the value prolonged from the coarse level (constant or linear: a choice for the radiation ADR; constant is the proposed default because it cannot make an intensity negative).
- Why not zero: a zero incoming intensity is a cold black wall. For the first pass it would remove most of the incoming radiation at every new face and cause a large spurious drop in `QR` near the face. FDS itself starts from ambient blackbody (`IL_S = IL_R = SIGMA*TMPA4/PI`, `main.f90:2395, 2399-2400`), which is acceptable only at initialisation.
- Wall `ILW` of a newly allocated wall slot starts at zero in FDS (`func.f90:4375-4376`). Do not let a regrid expose this. A new wall slot gets its `ILW` and `IL` from the nearest old wall record of the same surface, or from the coarse level, before the first sweep.
- FDS precedent for "geometry changed": when an obstruction is created or removed, `UPDATE_ALL_ANGLES` is set (`main.f90:1861-1868`) and the next radiation update sweeps all angles, then returns to the normal subset (`main.f90:1121`). Use the same switch after a regrid: one full-angle update. Cost: `NUMBER_RADIATION_ANGLES` is 100 by default (`read.f90:10225`) and `ANGLE_INCREMENT` is `MIN(5, NRA/15)` = 5 (`read.f90:10301`), so an ordinary update sweeps one fifth of the angles and a full-angle update is about five ordinary updates, paid once per regrid. A regrid interval that is short compared with an angle cycle makes this significant; check it against NFR at the Phase 8 gate.

### 1.5 Determinism (rule 4 and 4a)

All must hold, and are checked by V&V:

- Output byte-identical across 1, 2 and 4 ranks, and run to run.
- No dependence on box ownership or processing order inside a pass. This follows from 1.1.
- **(4a) Particle deposition is deterministic.** Current deposition into `QR_W` and `UIID` has no atomics in `part.f90` or `radi.f90` (sequential per mesh). The rule is that it must stay deterministic in the GPU port. Options, in order: per-box fixed-order accumulation, then a fixed-point sum (D-028). No floating-point atomics. The earlier gas-only scoping of this rule is dropped, or kept only until a GPU deposition kernel exists.
- **The CRITICAL sums.** `RAD_Q_SUM` and `KFST4_SUM` (CRITICAL section, `radi.f90:4120-4123` merged; FR-062 text cites the baseline 4119-4122) are scalar reductions, used by `CALCULATE_RTE_SOURCE_CORRECTION_FACTOR` (`main.f90:1823-1848`). They are covered by FR-005 (ii), not by rule 4. On the development machine set they must be summed in a fixed order (per-box partial sums, then reduced in global box order). The factor is state (`dump.f90:3950`, `cons.f90:513-518`), so restart must restore it.
- **No order-dependent fills.** The face fills of 1.4 must use only data that exist on every rank count (upwind interior or coarse values), not "whatever a neighbour had in its buffer".

### 1.6 V&V check and escalation (rules 4 and 5)

- Lag-error check: `radiation_gas_panel` split into 16^3 and 32^3 boxes, at K=1, 2, 3, against the single-mesh result. Tolerance is set by V&V under FR-005 (i) (exemption covers box-split dependence only; it does not cover the rest of radiation).
- If K<=3 misses tolerance, escalate in this order, and only these two:
  1. An optional larger-box radiation BoxArray, built from the geometry and a fixed box size only (never from rank count or load). Radiation boxes may then differ from the flow boxes, which needs a copy between the two layouts for `UII`, `QR`, `KAPPA` and `TMP`.
  2. An optional ordered sweep, in global box order only.
- Cost basis (my estimate, not a measurement): a 16^3 box at K=1 has a 46-plane dependency chain per angle. An ordered sweep over a 4^3 arrangement of such boxes has about 190. The lagged sweep therefore keeps all boxes busy at once; the ordered sweep serialises about four times the chain.

---

## Part 2. Radiation state at regrid

Legend. **Rebuild**: recompute from the new mesh; nothing to carry. **Transfer**: interpolate or restrict from old to new data. **Carry**: keep as is (not mesh-shaped). **Fill**: new cells need an initial value by rule 1.4.

| State | Defined at | What it is | At regrid |
|---|---|---|---|
| `UII` | `mesh.f90:51`, alloc `init.f90:680`, init `:850`, restart `dump.f90:3920, 4109` | Total intensity, cell array `(0:IBP1,0:JBP1,0:KBP1)`, read by the gas and particle source terms | **Transfer** (conservative or linear prolong to new fine boxes; restrict to coarse). It is the lagged quantity of 1.3, so it must be a good guess. Ghost layers are refreshed by the exchange. |
| `UIID` | `mesh.f90:338`, alloc `init.f90:854-855`, restart `dump.f90:3921, 4110` | Per-angle-subset partial sums, 4D `(…,UIIDIM)`, the accumulator for `UII` | **Transfer** with `UII` (it is restarted, so it is state, not scratch), or **Fill** uniformly from `UII/UIIDIM` and take one full-angle update (1.4). The second is the simple one. |
| `QR`, `QR_W`, `RADIATION_EMISSION` | `mesh.f90:47-49`, alloc `init.f90:657, 678-679` | Radiation source terms (gas, particle/droplet, emission) | **Rebuild** each radiation update (`radi.f90` ~4949-4983 overwrite them). Between updates `QR` feeds the energy equation, so a regridded level needs a **Transfer** of the last value (conservative restriction/prolong of the divergence, as for any source) so the step after regrid is not missing it. |
| `IL_S`, `IL_R`, `IL_R_OLD` (`OMESH`) | `type.f90:1035`; alloc `main.f90:2392-2400` (`INITIALIZE_RADIATION_EXCHANGE`, `:2384-2404`) | Per-neighbour face intensities `(NIC,NRA,NBANDS)` | **Rebuild** the layout (they are shaped by the neighbour list). **Fill** values by 1.4. These are the double-buffered face buffers of 1.1. Ghost fill from `IL_R`: `radi.f90:4354-4359`; `IL_S` fill: `:4824-4831`. |
| `OMESH` index lists (`IIO_R`, `JJO_R`, `KKO_R`, `IOR_R`, `IIO_S`, ..., `NIC_R`, `NIC_S`, `NIC_MIN`, `NIC_MAX`) | `type.f90:1037-1039`; built in `INITIALIZE_MESH_EXCHANGE_1` (`main.f90:2132-`, definitions `:2163-2172`, `NIC_R` count `:2174-2186`) | Which cells abut which mesh | **Rebuild** from the new BoxArray and level connectivity. In the AMR code this is the FillBoundary/ParallelCopy metadata; there are no `OMESH` pointers to keep. |
| `ALLOCATE_RADIATION_RECV_PKG`, `ALLOCATE_RADIATION_SEND_PKG` | `main.f90:4053-4092`, `:4097-4132`; called `:370-371`, `:561-562`, `:1102-1103` | MPI message buffers sized by `NIC` and the current angle subset | **Rebuild** after regrid (and, as in FDS, at each new angle cycle: the subset size can change). In the AMR code, the communication metadata is cached per BoxArray and invalidated by regrid. |
| `BOUNDARY_RADIA` records (`BR%BAND(n)%ILW(angle)`, `BR%IL(n)`) | type `type.f90:381-384` (`BAND_TYPE` `:183-185`); alloc `func.f90:4343-4382` (zero `ILW`, `IL=SIGMA*TMPA4/PI`); pack/unpack `func.f90:5128-5153` | Per wall, cface or particle: wall-face intensity per angle and band | **Transfer** with the wall (the walls move with their mesh cells). A wall on a new fine cell: **Fill** from the parent coarse wall, not zero (1.4). Slot bookkeeping: `BOUNDARY_RADIA_OCCUPANCY`, `NEXT_AVAILABLE_BOUNDARY_RADIA_SLOT`, `N_BOUNDARY_RADIA_DIM` (`func.f90:4343-4365`) are **Rebuild**. |
| `BR_INDEX` | wall `type.f90:395`, particle `:441`, cface `:1367` | Slot index of the above | **Rebuild**; the Phase 3 wall table must carry it (`02-…` and `00-r1-signoff.md`). Every wall has a nonzero slot (`INCLUDE_BOUNDARY_RADIA_TYPE` default true, `type.f90:1001`). |
| `CELL_ILW` | `mesh.f90:227`; alloc `func.f90:3879, 3903-3904`; use `radi.f90` ~4323, 4365, 4457-4459, 4482-4483, 4545 | Per solid cell, its wall intensities (thin-obstruction wall cells) | **Rebuild** from the new OBST/cell data, then **Fill** like wall `ILW`. |
| Wall `Q_RAD_IN` | `type.f90:312`; init `init.f90:1806`; reset `radi.f90:3887-3893`; accumulate `:4902-4920`, `:4965-4974` | Incoming flux, rebuilt from `ILW` at each update | **Rebuild** (zeroed and re-summed each update). The previous value is read by the solid phase between updates, so **Carry** it onto the new wall records until the next update. |
| `RAD_CALL_COUNTER`, `ANGLE_INC_COUNTER` | `mesh.f90:339`; init `init.f90:851-852`; incremented `radi.f90:3866` and ~4291, 4295, 4302, 4309, 4785, 4811, 4820; start-of-cycle test `main.f90:557` and `:1097`; restart `dump.f90:3949, 4138` | Where we are in the angle cycle | **Carry** unchanged, one global value (all meshes agree). Regrid must not reset it. The angle subset and the random rotation follow from it. |
| Direction coefficients (`RSA`, `DLX`, `DLY`, `DLZ`, `DLN`, `DLB`) | alloc `radi.f90:2827-2847` (`INIT_RADIATION`, `:2788-3296`); recompute `CALCULATE_DIRECTION_COEFFICIENTS` `:3431-3596` | Angle weights and direction cosines | **Carry**. Mesh-independent. Rotation is drawn on rank 0 and broadcast (`radi.f90:3443-3447`), only when `ALLOW_RANDOM_RADIATION_ROTATION` (`radi.f90:2862`: not cylindrical, not 2-D). With random rotation, `INTERPOLATE_IL` (`radi.f90:3651-3709`) re-bases `ILW` and `IL_R` at each new cycle (`main.f90:560, :1101`). It must run on every level and every box, and it needs the pre-regrid or post-regrid state to be consistent (see open item below). |
| `RTE_SOURCE` | `radi.f90` ~4220 alloc, ~4942 dealloc | Per-band work array | **Rebuild** each call; nothing to transfer. |
| `RTE_SOURCE_CORRECTION_FACTOR`, `RAD_Q_SUM`, `KFST4_SUM` | `cons.f90:513-518`; `main.f90:1823-1848`; `dump.f90:3950` | Global energy-balance correction | **Carry** (global scalars). Recomputed each update. |
| `UPDATE_ALL_ANGLES` | `main.f90:1121, 1868` | Flag for a full-angle update | **Set** by regrid (1.4). |
| Particle `ILW` (`BR` of a particle) | `func.f90:4465, 4528, 4641` (pack with particle) | Per-particle intensities | **Carry** with the particle (particles move between levels with the particle container). Not part of any mesh field. |

### Order of operations at a regrid

1. Tag and regrid flow fields (other roles).
2. Rebuild `OMESH`-equivalent connectivity and message metadata (the `_PKG` routines).
3. Transfer `UII`, `UIID`, `QR`, `QR_W`, `RADIATION_EMISSION` as cell data.
4. Transfer or rebuild wall records and `CELL_ILW` (after the wall table is rebuilt), and fill new walls.
5. Fill new box and level faces by 1.4 (upwind interior, then coarse), in a fixed order independent of rank.
6. Set `UPDATE_ALL_ANGLES`. Do not touch `RAD_CALL_COUNTER`.
7. The next radiation update does K passes with an exchange after each.

### Open items for the radiation ADR

- **Wall history across levels (Q3, Q4).** Which level owns a wall's incoming flux, how a thin OBST on a box or coarse/fine face behaves, and what happens to wall `ILW` when its cell is covered by a finer level, are still open. The table above gives the default (transfer with the wall; fill new walls from the parent).
- **Random rotation vs regrid.** `INTERPOLATE_IL` assumes the old angle set and the new one, then rewrites `ILW` and `IL_R`. A regrid in the same step as a new cycle must run the interpolation on the old data before the transfer, or fill after it. Choose one and test it.
- **Periodic self-coupling (Q7).** Whether `IL_S` self-coupling across a periodic boundary works through the generic exchange.
- **Output (Q5).** `RADF` and Smokeview/VTK on refined levels are outputs of `QR` and `UII`; no regrid state of their own.
