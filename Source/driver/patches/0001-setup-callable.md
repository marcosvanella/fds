# Patch 0001: FDS set-up callable from the C++ main (Source/main.f90)

Apply first. Touches only `Source/main.f90` (+52 lines, no deleted line); every change is inside `#ifdef WITH_AMREX`.

## What it does (WITH_AMREX only)
- `PROGRAM FDS` becomes `SUBROUTINE FDS_SETUP(MODE,FNAME,DT_OUT) BIND(C,NAME='fds_setup')`. `SAVE` keeps its locals alive like a main
  program, because the internal procedures (`MPI_INITIALIZATION_CHORES`, `STOP_CHECK`, `END_FDS`, ...) use them.
  The internal procedures are not moved or edited, so `MPI_INITIALIZATION_CHORES(1..6)` run exactly as before.
- `MPI_INIT_THREAD` is replaced by `MPI_QUERY_THREAD`: the C++ main has already called `MPI_Init_thread(FUNNELED)`.
- The input file name comes from `FNAME` (a C++ main has no Fortran command line for `GET_COMMAND_ARGUMENT`).
- Set-up runs up to the line before `MAIN_LOOP`, then `DT_OUT = DT` and `RETURN`.
- `MODE=1` (finish normally) and `MODE=2` (finish, marked set-up only) call `END_FDS` (output close, `MPI_FINALIZE`, `STOP`). The driver
  calls it after `amrex::Finalize()`, because AMReX must not outlive `MPI_Finalize`. `MODE=2` is used by the S1 skeleton, which has no time loop.
- Without `WITH_AMREX` the preprocessor removes every added line and the file is the original `PROGRAM FDS`.

## Evidence for USE_AMREX=OFF (IR-006)
Toolchain: `tests/env.sh` (GNU 14.2.0, Open MPI 5.0.7, threads 1). Build of the patched tree outside the reference tree
with `USE_AMREX` unset, same options as the reference build except `USE_SYSTEM_HYPRE/SUNDIALS=ON` (the offline copies of the same
HYPRE commit and SUNDIALS 7.5.0) and fixed date/version strings.
1. Control: the unpatched tree, same options, run on `shunn3_32` gives output bitwise identical to the baseline
   (16/16 files: `_hrr.csv`, `_mass.csv`, `_mms.csv`, `_1.restart`, 12 `.sf`). This shows that the scratch build reproduces the baseline.
2. Patched tree, `USE_AMREX` OFF, `tests/check_off_bitwise.sh`:
   - `shunn3_32`, 1 rank: PASS, 16 output files bitwise identical to the baseline.
   - `shunn3_4mesh_32`, 4 ranks: PASS, 47 output files bitwise identical to the baseline.
3. Machine code: `objdump -d` of the control and the patched executables (same options, same source path) has the same number of
   instructions; 67 of 1.86M lines differ, all `mov $imm,mem` whose immediate is a source line number (`+48`/`+41`, the size of the patch
   above the call site); no opcode, register or address differs. The preprocessed `main.f90` (`gfortran -E -P`) is identical to the original
   apart from blank lines.

## Not done here
The patch does not run the FDS time loop in the WITH_AMREX build; the driver replaces it (S5).
