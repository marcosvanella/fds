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
7. FDS velocity arrays start at -1 in their own direction (`U(-1:IBP1,..)`, `V(..,-1:JBP1,..)`, `W(..,..,-1:KBP1)`); a reader of the shape-only dump records must use that lower bound (this was the cause of the former open-wall finding).

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
| L1209 | 12 + 4 | bitwise equal for every wall kind, OPEN-boundary walls included (after the reader fix, see "L1209 on OPEN walls") |

The committed archive `frozen/fds_loops_cases/fds_loops_cases.tar` (23 files, ctest `pb_fdsloops_fds`) holds L1211, L1207 and the H
boundary fill, and five L1209 files: the closed run (step 1, predictor and corrector) and the two-mesh run (step 3, mesh 1 predictor and
corrector, mesh 2 predictor), chosen because they exercise the OPEN-wall branches (inflow ramp, wind, pressure ramp). The negative control
`pb_fdsloops_fds_uu_lb0` reads the same files with the former wrong lower bound and must fail (314 elements differ).

### L1209 on OPEN walls: cause found and closed

The former finding "OPEN-wall BXS..BZF are not reproducible from the dumped inputs" was a defect of the comparison, not of `pres_poisson_boundary_arrays`
and not of the dump point.

1. **Cause.** FDS allocates the velocity arrays with one extra cell in their own direction: `U/US(-1:IBP1,0:JBP1,0:KBP1)`, `V/VS(0:IBP1,-1:JBP1,0:KBP1)`,
   `W/WS(0:IBP1,0:JBP1,-1:KBP1)` (init.f90). The dump hook writes shape-only records (assumed-shape dummy arguments lose the bounds), and the
   reader (`harness/loops_modes.cpp`) used lower bound 0 for UU, VV and WW. Every velocity was therefore read one cell too low in its own direction:
   `WW(I,J,KBAR)` was really `WW(I,J,KBAR-1)`, and `UU(0,J,K)` was `UU(-1,J,K)`. The open-wall branch compares exactly those normal velocities with zero,
   so it picked the wrong branch wherever the normal velocity changes sign or is zero next to a non-zero one (closed run, predictor at the start-up
   pass: 43 of 99 BZF values; two-mesh run: 80 of 99 BXS and up to 69 of 143 BZF values). The synthetic cases of `ref_loops.f90.in` use 0:IBP1 for all
   three arrays, which is why the verbatim-Fortran test never saw it. The wall-based arrays and the other loops (HP, KRES, FV*, BXS..) are 0-based in FDS
   and were read correctly, which is why only OPEN walls differed.
2. **How it was found.** Scratch FDS (-O0, gfortran, only `pres.f90` recompiled with the dump hook) with prints at the end of the open branch
   (`ICYC`, `PREDICTOR`, I, J, WW(I,J,KBAR), KRES(I,J,KBAR), H0, BZF) and one at each call of `PRESSURE_SOLVER_COMPUTE_RHS`. Findings: (a) FDS
   redirects unit 6 into `<CHID>.out`, so the prints written earlier were not missing, they were in the .out file (the "print never fired" observation was
   wrong); (b) the routine is called twice at ICYC 1 for both predictor and corrector (the start-up pass at T = 0 and the first step), the second call overwrites the
   first dump file, which is expected; (c) at the dumped call FDS's own branch choice and its BZF values are consistent with `WW(I,J,KBAR)` as printed
   inside the loop: there is no overwrite between the loop and the dump, no ramp or `H0` override, and the dump point is right.
3. **Fix.** The reader uses lower bound -1 for UU/VV/WW in files written by the hook (files that carry `WV_RAMP`); `uu_lb0=1` restores the old
   reading as a negative control. The function itself is unchanged.
4. **Result.** L1209 bitwise equal to FDS output for all walls, closed run steps 1 and 2 (4 files, 2552 elements) and two-mesh run steps 1 to 3 (12 files,
   8616 elements), predictor and corrector, 0 differing elements. The C++ function may replace the FDS loop on open boundaries; the NIC>1 and synthetic-eddy
   branches remain untested against FDS (the dump hook stops on synthetic eddies).
