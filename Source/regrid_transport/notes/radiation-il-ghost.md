# L1251: radiation intensity ghost table (IlGhost.H)

The AMR replacement of FDS's neighbour-mesh intensities `IL_S`/`IL_R` (docs radiation 06, section 17.2). Entry point: `fdsrt::build_il_ghost_table`
(`IlGhost.H`/`IlGhost.cpp`, host C++, run once per regrid). It builds index tables only. No intensity is interpolated or averaged here.

## Contract with the sweep (accepted by the Radiation Lead)
- The layout is one flat array `ILR(NICR,NRA,NSB)` with the position index fastest. `ILS` has the same layout and the same positions.
- **Pack.** During the sweep of angle `N` and band `IBND`, `s5rad_ils` does `ILS(LL,N,IBND) = IL(IIO_S(LL),JJO_S(LL),KKO_S(LL),ILS_BOX(LL))`, as at radi.f90:4831. The tables give `ILS_BOX`, `IIO_S`, `JJO_S`, `KKO_S`, with `NICS = NICR`.
- **Lag.** After each radiation update, the caller copies all of `ILS` into `ILR`. This is the same-process `IL_R = IL_S` of `MESH_EXCHANGE`. A pointer swap is not equivalent: with `ANGLE_INCREMENT > 1`, the angles not swept in an update keep their old values.
- **Read.** At an `INTERPOLATED_BOUNDARY` wall, `s5rad_wall1` reads the mean of `ILR(W_NIC_MIN:W_NIC_MAX,N,IBND)`, verbatim radi.f90:4357-4361. The tables give `W_NIC_MIN` and `W_NIC_MAX`. That mean is the only sum; the fill writes raw values.
- **Angle rotation.** `ILR` is re-interpolated by the existing `s5rad_ilw_old_copy` / `s5rad_ilw_interp` kernels with an identity gather (`NBR = NICR`). Nothing is added here.

## Table
- **Input.**
  - The radiation boxes: level, and the index box in that level's space.
  - The per-level refinement ratios. Per direction; use 1 in a single-cell direction.
  - The walls: receiver box, FDS `IOR`, and ghost cell, in the caller's wall order.
- **Sender cells for a wall with ghost cell g.**
  - Equal level: g itself.
  - Coarser sender: the cell that contains g. That is one entry, so `W_NIC_MIN = W_NIC_MAX`.
  - Finer sender: g refined and cut to the layer next to the face. The entry count is the product of the tangential ratios, and the sweep's mean is the area mean of `ASSIGN_GHOST_VALUE`.
  - These are the cells the FDS probe search finds (init.f90:3139-3176).
- **Level search order.** The receiver's own level first, then the coarser levels (nearest first), then the finer levels.
  - A disjoint composite tiling (FDS-like meshes) has exactly one covering box.
  - With full levels (coarse under fine), a box takes its ghosts from its own level, then from the coarser level.
  - A face that spans two boxes aborts, the analogue of FDS ERROR(431).
  - If no box covers g, the wall gets `W_NIC_MAX = W_NIC_MIN - 1` and is not an interpolated wall.
- **Positions.**
  - Grouped by (receiver box, sender box) in ascending order, so each neighbour pair is one contiguous range (`pair_lo`).
  - Inside a pair, the walls keep the caller's order, and the cells of one wall are consecutive: K, then J, then I fastest. This is the FDS `IIO_R` order (main.f90:2191-2254).

## Tests (`test_il_ghost`, ctest `regrid_transport_il_ghost`)
- Every box-face ghost cell is a wall. The table is compared against a coordinate-based reimplementation of the FDS `IIO_R` search: two probe points at ±0.475 of the receiver cell, at the mesh separation distance outside the face.
- Checks per wall:
  - the count, the sender cells and their order;
  - the `s5rad_ils` pack, then the copy `ILS → ILR`, then the `s5rad_wall1` mean, equal to the FDS mean of `IL_R` bitwise;
  - equal level: the raw neighbour value;
  - the value is non-negative;
  - pair ranges are contiguous and ascending, and the walls of a pair are consecutive.
- Cases:
  - equal level, three 3D boxes, with one face shared by two senders;
  - coarse | fine | coarse at ratio 2 (3D);
  - 2D at ratio (4,1,4);
  - three levels at ratios 2 and 3;
  - the full-level search order.
- Result: 7370 checks, 0 failed. Negative control: cutting the wrong fine layer gives 660 failures.

## Not done here (interface)
- **Multi-rank exchange.** Each `pair_lo` range, over all (`N`, `IBND`), is the one-way message from the sender box's rank to the receiver box's rank. It is packed as `(range, NRA, NSB)`, contiguous.
  - The MPI send and receive are not written. They need the sweep host's box-to-rank layout of `ILS`/`ILR` (all boxes on every rank, or rank-local rows), which is the Radiation Lead's choice.
  - On one rank, the whole-array copy above is the exchange.
- **Wall list.** The full wall list and its `W_BTYPE = INTERPOLATED_BOUNDARY` marking stay with the sweep host's wall table. A wall is interpolated exactly when its range is non-empty.
- **Which boxes are swept** (a composite tiling or full levels) is the Radiation Lead's call. The table supports both.
