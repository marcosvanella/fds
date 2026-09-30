# Patch 0004: `wall.f90`, `ASSIGN_GHOST_VALUE` skips box-boundary cells when they are filled externally (guarded, `WITH_AMREX`)

Apply after 0001 and 0002 (same module `FDS_AMREX_HOOKS` as 0003). Touches only `Source/wall.f90` (+6 lines, all inside `#ifdef WITH_AMREX`).

## What it does
`USE FDS_AMREX_HOOKS, ONLY: EXTERNAL_GHOSTS_FILLED` in `WALL_ROUTINES`, and in `ASSIGN_GHOST_VALUE`, right after `IF (EWC%NOM==0) RETURN`:
`IF (EXTERNAL_GHOSTS_FILLED) RETURN`. With the flag set the `RHOP`, `ZZP`, `RSUM`, `TMP` ghost cells of a cell across a box boundary (first and second layer)
are the ones the AMReX level fill wrote from the neighbour box, not the OMESH average. `WALL_BC` then still does everything for `NOM==0` walls (solid,
open, non-periodic domain sides): only the `NOM>0` walls are skipped, the routine is otherwise unchanged. The `INTERPOLATED_BC` species-flux branch
(`wall.f90` ~821-1005, `EWC%NIC>1` only) is not changed: it needs different-level neighbours, which M2a does not have.

Value note (S4 finding, README "Plain ghost modes"): for a same-level neighbour of equal cell size `ARO=1` and the OMESH average is the plain neighbour value, so the
AMReX fill equals it bitwise for RHO, RHOS and ZZ; TMP, RSUM are derived (`PBAR_P/(RSUM*RHOP)`), so the driver has to recompute them on the AMReX ghost
layer when the flag is used (S5), or keep the ghost fill of TMP as is.

## Evidence
Same as 0003 (one run for 0003+0004+0005): `git apply --check` passes in order; `USE_AMREX=OFF` bitwise check `shunn3_32` 1 rank PASS (16 files),
`shunn3_4mesh_32` 4 ranks PASS (47 files); `USE_AMREX=ON` build with the flag `.FALSE.`: kernel checks and driver tests pass as in 0003. gfortran 14.2 only; oneAPI not checked.

## Not done here
The flag TRUE path is not exercised (S5). `CC_IBM` is not touched.
