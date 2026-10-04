# FDS pressure-area loops on the CPU (Role 2)

Files: `FdsPressureLoops.H/.cpp` (namespace `pressure_backend::fdsloops`), `harness/loops_modes.cpp` (executable `pb_fds_loops`, no AMReX
dependency), `tests/fds_loops_test.py`, `tests/m5_fds_loops.cmake`, `tests/fds_loops/` (reference generator). Source of the loops: FireX
36975d7 `Source/pres.f90`.

## What is there

| Loop | Upstream lines | Function | Notes |
|---|---|---|---|
| L1211 | pres.f90:250-260 | `pres_compute_rhs_div` | PRHS for IPS 1/4/7, non-cylindrical. IPS 2/3/5/6 (L1212-L1214, transposed store) and the cylindrical branch (L1210) are not covered. |
| L1207 | pres.f90:758-764 | `pres_p_from_h` | `P = RHOP*(HP-KRES)` over the full box 0:IBP1 in every direction, ghost cells included. |
| L1220/21/22 | pres.f90:450-462, 466-477, 481-492 | `pres_h_bc_x/y/z` | H boundary fill from BXS..BZF with the LBC/MBC/NBC codes. Ordered IF sequence; codes 5 and 6 exist only in x. |
| L1209 | pres.f90:65-228 | `pres_poisson_boundary_arrays` | BXS..BZF from wall data. The caller supplies neighbour widths, wind arrays and an `evaluate_ramp` callback. |
| L1215 | pres.f90:394-400 | none (retired) | `pressure_backend` writes phi (H) directly (FR-037); the driver copies phi to H on valid cells. |

Each function is tagged `HAND-WRITTEN <loop> pres.f90:a-b (FireX 36975d7)`. The arrays are Fortran-ordered views (`F1`, `F2`, `F3`) with
explicit lower bounds and leading dimensions, so the same code serves FDS-shaped arrays (PRHS(ITRN,JTRN,KTRN)) and ghosted AMReX arrays.
A view that does not cover the loop's index range throws `std::invalid_argument`.

The floating-point order is the Fortran order, and the library file is compiled with `-ffp-contract=off`. This is why the comparison is
bitwise and not a tolerance.

## Tests

`ctest -R pb_fdsloops` (21 tests):
- `pb_fdsloops_ref`: builds the verbatim-Fortran reference. `tests/fds_loops/fds_loops_gen.py` inserts the upstream loop text, pinned by
  SHA-256 of each line range, and refuses to run if `Source/pres.f90` has changed.
- `pb_fdsloops_bitwise`: 93 case files, every output element compared bitwise (max abs 0, max rel 0). L1211: 10 (five sizes including
  JBAR=1 and 1x1x1, PRHS leading dimension IBAR+1 and IBAR); L1207: 5; L1220-22: 70 (every x code 0-6 against every y/z code 0-4);
  L1209: 8 (all six IORs, Neumann/Dirichlet, solid/open/interpolated/null walls, wind on/off, T_IGN window, exact-zero velocities).
- `pb_fdsloops_mutant_1..18`: each deliberately breaks one detail in the C++ and must fail the bitwise test (expect-fail tests).
- `pb_fdsloops_drift`: fails if the upstream line ranges no longer match the pinned hashes.

## Traps found

1. PRHS leading dimensions are ITRN/JTRN/KTRN, not IBAR+1: ITRN is IBAR when x is periodic, JTRN is 1 when JBAR=1.
2. BXS/BXF have J extent JDIM=1 when JBAR=1 (2-D).
3. L1207 runs over the full box including ghost cells; RHOP is allocated -1:IBP1+1.
4. The H boundary fill must run before P is formed from H.
5. In L1209 the VEL_EDDY choice is made by the vent's IOR, not the wall's; V_WIND(K) is indexed by K in all three directions; the TSI
   selection depends on T_IGN within 20 epsilon and on a ramp index of at least 1.
6. NIC>1 branches of L1209 are not reached in the AMR route (D-072).

## Check against real FDS output

`frozen/fds_loops_dump_hook.py` edits a SCRATCH copy of `pres.f90` (never the reference tree) so that a normal FDS run writes, when the
environment variable `FDS_LDUMP_STEPS` is set (for example `1,2,3`, single OpenMP thread), the inputs and outputs of the loops as stream
records: `L1211_fds_n<ICYC>_<P|C>_m<NM>.bin`, `L1207_...`, `HBC_...` (the H boundary fill) and `L1209_...`. With the variable unset every
inserted call returns at once; the hook only reads FDS arrays. Inputs: `frozen/fds_cases/ldump_closed.fds` (one mesh, solid walls with an
inflow velocity ramp, open top, fire) and `ldump_2mesh.fds` (two meshes, interpolated faces, open vents with wind and a pressure ramp). The
scratch FDS was the reference-tree object files with only `pres.f90` recompiled (gfortran, -O1) and relinked, FireX 36975d7.

`pb_fds_loops dir=<directory>` runs the C++ functions on these files and compares every output element bitwise (`verbose=2` also prints the
first differing values). Results over steps 1 to 3, predictor and corrector, both meshes:

| Loop | Real-FDS files | Result |
|---|---|---|
| L1211 | 12 + 4 | bitwise equal, max abs 0, max rel 0 |
| L1207 | 12 + 4 | bitwise equal |
| L1220-L1222 | 12 + 4 | bitwise equal |
| L1209 | 12 + 4 | Neumann, solid-Dirichlet and interpolated walls: no difference. OPEN-boundary walls: **not equal**, see below |

The committed archive `frozen/fds_loops_cases/fds_loops_cases.tar` (18 files, ctest `pb_fdsloops_fds`) holds L1211, L1207 and the H
boundary fill. L1209 files are not in it.

### Open issue: L1209 on OPEN walls against a real FDS dump

For walls with `BOUNDARY_TYPE==OPEN_BOUNDARY` the C++ function and the verbatim Fortran loop (synthetic cases, 8 bitwise-equal files, all six
IORs, wind, ramp) agree, but the arrays dumped from a real FDS run differ. Closed run: 43 of 99 BZF elements, exactly those with
`WW(I,J,KBAR)>0` (the C++ takes the `KRES` branch as the code says; FDS's value is the `H0` branch value). Two-mesh run (wind, pressure ramp): all
BXS values on the wind-inflow face and most BZF values differ by 1 to 2 percent (C++ gives `HP(1,J,K)`, because the dumped `UU(0,J,K)` equals
`U_WIND(K)`; FDS's value is larger). Both are what the loop would produce if `UU` and `WW` at the open-boundary faces had other values while the
loop ran than in the dump taken at the end of the routine. The same dumps agree bitwise for all other wall kinds, and the ramp argument (TSI)
agrees bitwise with FDS's own. A print inserted in the scratch copy at the top of the wall loop (and at the end of the open branch) never fired for
the wall in question, so the cause is not yet understood. Treat the open-wall branch of `pres_poisson_boundary_arrays` as verified against the Fortran
text only, not against FDS, until this is resolved.
