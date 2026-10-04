# Sign-off note: patch UP-0011 (opt-in write-only dump of the pressure right-hand side and solved H)

Patch file: `UP-0011-pressure-rhs-dump.patch` (Source/pres.f90: new routine `PRESSURE_SOLVER_DUMP` and helper `POISSON_FACE_TYPE`, 1 line in the PUBLIC list;
Source/main.f90: two 6-line call blocks). Target: FireX (tested on 36975d765f); applies with `patch -p1` to master ce1f659cd4 (offset -131 lines in main.f90 and -4 lines for the second `Source/pres.f90` hunk; the first pres.f90 hunk applies without offset)
and to the merged tip bee11f0329 (offset +73 lines in main.f90, none in pres.f90). Author: Role 2 (pressure backend). State: proposed; the patch is only meant to produce check data and may
stay local to the AMR branch if the project owner does not want it upstream.

## What changes (exact hunks)
1. `Source/main.f90`, in `PRESSURE_ITERATION_SCHEME` (routine start main.f90:1674 at the FDS-AMReX tree), inside the `PRESSURE_ITERATION_LOOP` construct (a loop label, `DO` at main.f90:1696; not a routine), directly after the `PRESSURE_SOLVER_COMPUTE_RHS` loop and before `SELECT CASE(PRES_FLAG)`:
   `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX ; CALL PRESSURE_SOLVER_DUMP(NM,T,DT,.FALSE.) ; ENDDO`
2. `Source/main.f90`, directly after the `END SELECT` of the solver choice and before the residual check:
   the same loop with `.TRUE.` (writes H or HS).
3. `Source/pres.f90`: `PRESSURE_SOLVER_DUMP` added to the PUBLIC list; the routine and `POISSON_FACE_TYPE` are added at the end of module PRES (about 190 lines).
Both calls are outside OpenMP regions, once per mesh of the rank.

## Behaviour with the switch off (default)
`FDS_PDUMP_STEPS` is read once (`GET_ENVIRONMENT_VARIABLE`, first call), stored in SAVEd variables, and every later call returns at its first executable
statement. No FDS variable is read for the decision or modified, no MPI call, no file. With the switch on the routine only reads FDS arrays (`PRHS`,
`H`/`HS`, `RHO`/`RHOS`, `KRES`, `PRESSURE_ZONE`, `BXS..BZF`, mesh extents) and writes files in the run directory; `POINT_TO_MESH(NM)` is the same pointer
association the neighbouring routines do. Nothing is written to unit numbers other than ones taken from `GET_FILE_NUMBER`.

## Behaviour-unchanged check
Builds (scratch copies of the reference tarball `firex-36975d7-src.tar.gz`, same CMake flags, `-DUSE_HYPRE=OFF -DUSE_SUNDIALS=OFF`, GNU 14.2 Release, OpenMP):
control (unpatched) and patched. Comparison with `vv-runs/tools/cmp_runs.py` (bitwise except the documented timing strip), `&DUMP SIG_FIGS=17`.

| Case | Ranks | Compared | Result |
|---|---|---|---|
| `csmag_32` with `FISHPAK_BC=0,0,0` (periodic, FFT), control vs patched, switch off | 1 | 40 output files | BITWISE (31 identical, 9 identical after the timing strip) |
| same, switch on (`FDS_PDUMP_STEPS=1,3,10`) | 1 | 40 output files + 42 dump files (EXTRA) | the 40 identical, only the dump files are extra |
| unpatched control vs the V&V reference binary `refbin/gnu_ompi_firex-36975d7/fds` | 1 | 40 | identical except the two banner lines "Hypre library" (the control build has no HYPRE) |
| `shunn3_4mesh_32` (4 meshes, MPI), control vs patched, switch on (steps 1,2,3) | 4 | 90 output files + 168 dump files (EXTRA) | the 90 identical |
| `ns2d_16` (Dirichlet/Neumann faces, y one cell), control vs patched, switch on (steps 1,2,3) | 1 | 22 output files + 42 dump files (EXTRA) | 21 identical; `ns2d_16_1.restart` differs in 12 of 1569876 bytes (values near 1e-310). The same file differs by 9 to 12 bytes between two runs of the unpatched control, and between the control and the patched binary with the switch off (restart is run-to-run unstable in this case), so it is not an effect of the patch |
Not run: the full Tier 1 set, the `-fcheck=all` build (the routine uses only array sections of the mesh's own bounds; a debug run of one case is still advisable before upstreaming), GLMAT/ULMAT cases, cylindrical and tunnel cases (the routine writes a PRHS that is not the plain FFT layout for those; documented, not exercised).

## Enable and format
Set `FDS_PDUMP_STEPS` to the step counts (ICYC) wanted, e.g. `FDS_PDUMP_STEPS=3` or `1,10`, or `ALL`, before starting FDS. Files are written in the run directory, per mesh `m<NM>`, stage `P`
(predictor, H) and `C` (corrector, HS): `<CHID>_pdump_n<ICYC>_<P|C>_m<NM>_{meta.txt,rhs.bin,phi.bin,rho.bin,kres.bin,zone.bin,bc.bin}`. The binary files are float64 (zone as float64),
valid cells, I fastest, the same layout as Role 1's `FDSTL_PDUMP` files (`pdump_<icyc>_<P|C>_{rhs,phi}.bin`, written by the AMR driver, which carries no metadata, no rho, no KRES and no
boundary data); the text file `meta.txt` (format `fds-pdump-1`) carries mesh extents, IBAR/JBAR/KBAR, the cell faces, LBC/MBC/NBC and the decoded face types, `PRES_FLAG`, `IPS`, `FISHPAK_BC`,
step, stage, T, DT, the largest Poisson boundary datum and the number of pressure zones. With the pressure iteration repeated in one step (`ITERATE_PRESSURE`) the last iteration's files remain.
Reader: `Source/pressure_backend/harness/m4_modes.cpp` (`pb_fds_frozen`); committed case and results: `Source/pressure_backend/frozen/fds_csmag32_periodic/README.md`.
