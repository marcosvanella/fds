# 04 · Solid Phase sign-off on blocked-loop families (SP1-SP4)

Owner: AMR Solid Phase Lead · Reviews `amrex/blocked-loop-families.md`, Solid Phase group · Source line numbers refer to FireX `36975d765f` (the working tree has since moved and its wall.f90 numbers are 9 to 12 lines higher).

| Family | Verdict |
|---|---|
| SP1 `CELL_COUNTER` running average (velo.f90:306-351) | **Accept** with two conditions |
| SP2 ghost mirror (velo.f90:355-363) | **Change**: the uniqueness assertion is wrong |
| SP3 `DP` accumulate (divg.f90:532-554) | **Accept** with two conditions |
| SP4 `WALL_BC` loops and wall heat transfer (wall.f90:109-194, 888-947, 3635-3717) | **Change**: the "not races" statement is wrong |

## SP1 · accept
The gather over each gas cell's solid walls in ascending wall index reproduces the serial recurrence exactly. The running average `MU=(1-WGT)*MU+WGT*x` with `WGT=1/n` depends on order, and bit-for-bit equality requires the same recurrence in the same order, not a closed-form mean. Conditions:
1. The list holds only `SOLID_BOUNDARY` walls (velo.f90:310), and in the DNS branch (velo.f90:347) the store `MU=MU_DNS` is applied to each such cell. The kernel must also not touch cells with no solid wall: they keep their earlier `MU`.
2. `CELL_COUNTER` (`IWORK1`) is only scratch here. With a gather it can be a register. Check that no later routine reads `IWORK1` expecting these counts.
No physics change; the order is the one the serial code uses.

## SP2 · change
The proposed assertion that the `(II,JJ,KK)` targets are distinct will fail. `II,JJ,KK` is the solid cell on the far side of the face, and a solid cell has one wall for every exposed face. Corner cells of any obstruction, and every cell of a one-cell-thick obstruction (gas on both sides), are the target of two or more walls with different source gas cells. The serial loop is last-writer-wins, so the highest wall index sets `MU` and `KRES` in that cell.
Rewrite: gather per target cell and take the source of the highest wall index that passes the `SOLID .OR. EXTERIOR` test (velo.f90:358). Do not assert uniqueness. This keeps the result bitwise. Written in any other order, solid-cell `MU` and `KRES` change, and so do the wall-adjacent stress terms that read them. The proposed check "no target equals a gas read cell" is correct and cheap; keep it. Thin walls, where `II,JJ,KK` is a gas cell, are skipped by the test and should stay skipped.

## SP3 · accept
Per-cell gather in ascending wall index, starting from the `DP` already in the cell, reproduces the serial subtraction order. Conditions:
1. The list excludes `NULL_BOUNDARY`, `INTERPOLATED_BOUNDARY` and `OPEN_BOUNDARY` walls (the open-boundary wall only stores `K_G` and cycles, divg.f90:536-540).
2. The thin-wall `CYCLE` on `IOR<0` guards only the `KDTD*` zero stores, not the `DP` term (divg.f90:544-545). Both sides of a thin wall still subtract from their own gas cells. Keep that distinction in the kernel. The `KDTD*` stores write a constant, so a duplicate store is harmless, but keep the `CYCLE` so results do not depend on thread timing.

## SP4 · change
The family text says the walls are independent and nothing races. That holds for the thermal solve, but the wall loops also add into shared state, and upstream protects each with `!$OMP CRITICAL`, so the order is already nondeterministic under OpenMP and serial order is the reference:
- `OBSTRUCTION(OBST_INDEX)%MASS` is reduced by every consuming wall of an obstruction (wall.f90:1390-1391, corrector only), and for external walls it is copied from or posted to the neighbour mesh (wall.f90:1379, `REAL_SEND_PKG8`); a third writer sets it to -1 at burn-through (wall.f90:2695). Several walls share one obstruction, so this is a many-to-one sum, and the result decides burn-away.
- `D_SOURCE` and `M_DOT_PPP` in the adjacent gas cell receive the layer-removed and thin-obstruction gas production (wall.f90:1412-1415 and 1420-1422) and deposition (wall.f90:1730-1733). A cell with several walls is a many-to-one sum, and it feeds the divergence directly.
- The same pattern appears for solid particles (wall.f90:1519-1556: `D_SOURCE`, `M_DOT_PPP`, `Q_DOT(3,4,8)`, `M_DOT`), and `HT3D_TEMPERATURE_EXCHANGE` stores into `M%TMP` (wall.f90:3676).
Required rewrite: the kernel writes each wall's term to a per-wall scratch slot (mass loss per wall, gas production per species per wall). A second pass sums them in ascending wall index per obstruction (a CSR list wall to obstruction) and per gas cell (the CSR list the document already proposes for SP1 and SP3). Do not use device atomics. Atomics would add in an arbitrary order and break FR-005 (i), and `OBSTRUCTION%MASS` drives burn-away, so one different last digit can change the step on which a cell is removed. The external-wall copy from the neighbour mesh is a data dependency, not a sum, so keep it as an explicit exchange step after the sums. Budget terms (`Q_DOT`, `M_DOT`) are global sums and follow the existing rule for the divergence sums: terms on the device, adds in order.
Further notes:
- Solid-phase random draw: `MASS_FLUX_VAR` (wall.f90:1324-1337) uses the box-layout-dependent random stream. On the device it needs the counter-based generator keyed by face key and step (SP-R9). It is statistically equivalent to the CPU draw, not bitwise, so cases that use it cannot take the bitwise test.
- Pointer aliasing: `B1_BACK`, `BC_BACK`, `ONE_D2` and `OMESH` pointers (wall.f90:1880-1884 and 1913-1927; thin-wall pass near 3670) read other walls' records. These are read-only and fine, but the back-wall record must be read before it is written in the same step. The two-sided thin-wall solve (wall.f90:3635-3717) writes `ONE_D%TMP` of the wall and reads the other side's `DELTA_TMP`, so keep it as a separate pass after the 1-D solve, as the CPU code does.
- The staged port order in the document (leaf callees, ragged per-wall tables, then the solve) is fine. Stage 3 (the solve) needs A1-A6 of FR-041b records only once AMR is on; on a static mesh the port is one record per wall.
