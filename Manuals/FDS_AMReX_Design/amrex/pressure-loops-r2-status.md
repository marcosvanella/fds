# Pressure-area loop translation: status (pressure backend)

Scope: the five loops of the Legacy Mapper's ranked list in the pressure area, FireX 36975d7 `pres.f90`. Code and tests are in
`Source/pressure_backend` (`FdsPressureLoops.H/.cpp`, `harness/loops_modes.cpp`, `tests/fds_loops_test.py`, `tests/m5_fds_loops.cmake`,
`tests/fds_loops/`, `frozen/fds-loops-notes.md`). Claims are in `loop_claims.csv`. The generator worktree was not edited.

## Per loop

| Loop | pres.f90 | Result | Tests |
|---|---|---|---|
| L1211 | 250-260 | CPU reference `pres_compute_rhs_div` (IPS 1/4/7, non-cylindrical), same operation order as the Fortran. Generator entry pending (generator owners). | Bitwise vs the verbatim loop: 10 cases (five sizes incl. JBAR=1 and 1x1x1, PRHS leading dimension IBAR+1 and IBAR). Bitwise vs arrays from a real FDS run: 16 files (max abs 0, max rel 0). |
| L1207 | 758-764 | CPU reference `pres_p_from_h`, full box 0:IBP1 in every direction. Generator probe accepts the entry as written. | 5 cases bitwise; 16 real-FDS files bitwise. |
| L1215 | 394-400 | Retired: the backend writes H directly (FR-037); the driver's valid-cell phi to H copy is a plain MultiFab copy. L1216-L1218 would be retired by the same reasoning (not claimed). | none needed |
| L1220-L1222 | 450-492 | Host code `pres_h_bc_x/y/z`, ordered IF sequence, codes 5/6 in x only. | 70 cases (every x code 0-6 against every y/z code 0-4) bitwise; 16 real-FDS files bitwise. |
| L1209 | 65-228 | Host code `pres_poisson_boundary_arrays`. The caller supplies neighbour widths, wind arrays and the ramp callback. NIC>1 not reached (D-072). | 8 cases bitwise (all six IORs, Neumann/Dirichlet, solid/open/interpolated/null, wind on/off, T_IGN window, exact-zero velocities). Real FDS: bitwise for Neumann, solid-Dirichlet and interpolated walls; **not** for OPEN-boundary walls (open issue). |

Every function is tagged `HAND-WRITTEN <loop> pres.f90:a-b (FireX 36975d7)`. The shared reference file is compiled with
`-ffp-contract=off`.

Tests: `ctest -R pb_fdsloops`: the reference fixture, the bitwise test (93 case files, max abs 0, max rel 0), 18 mutants that must fail
(each changes one operand, sign or index; two survived at first because the stand-in ramp for index 0 was constant and no velocity was exactly
0; both fixed), a drift check (the verbatim ranges are pinned by SHA-256 and the generator step refuses if `pres.f90` moves), and the real-FDS
archive test `pb_fdsloops_fds`.

## Real FDS dumps

`frozen/fds_loops_dump_hook.py` patches a scratch copy of `pres.f90` to write the loop inputs and outputs (write-only; off unless
`FDS_LDUMP_STEPS` is set). Two small inputs (one mesh with an inflow ramp and an open top; two meshes with wind and a pressure ramp), steps 1 to
3, predictor and corrector. The patch is not an upstream patch (it is a test aid and touches only a scratch copy), so it has no number in
`upstream-patches`.

## Traps found

1. PRHS leading dimensions are ITRN/JTRN/KTRN, not IBAR+1: ITRN is IBAR when x is periodic, JTRN is 1 when JBAR=1.
2. BXS/BXF have J extent 1 when JBAR=1.
3. L1207 runs over the full box including ghost cells; RHOP is allocated -1:IBP1+1.
4. The H boundary fill must run before P is formed from H.
5. L1209: VEL_EDDY is selected by the vent's IOR, not the wall's; V_WIND(K) is indexed by K in all three directions; the TSI selection depends on T_IGN within 20 epsilon and on a ramp index of at least 1.
6. Role 1's driver already calls the Fortran `fds_p_h_ghost` for physical-face H ghosts; the C++ version serves the fine levels.

## Generator owners

See `pressure-loops-generator-proposal.md`: the L1207 entry is accepted as written; L1211 needs the extents of an `exact` allocation (ITRN,
JTRN, KTRN) to become integer arguments for cell-loop kernels (today only the rank-2 wall-kernel path does that). Same "layout contract" as
L1212-L1214. Neither entry was added to a sidecar by this work.

## Left

- L1209 on OPEN walls against real FDS (above). The verbatim-Fortran test passes, so either the dump (arrays taken at the end of the routine) is not what the loop saw, or the loop text and the build differ in a way not yet found.
- L1210 (cylindrical PRHS), L1212-L1214 (transposed PRHS), L1216-L1218 (transposed copies; retire like L1215).
- Register: the table in `pressure-velocity-h-loops.md` and the work list need regenerating by the Legacy Mapper (`tools/inventory/loop_work_list.py`); the statuses above were written directly in `loop_claims.csv`.
