# Patch 0006: `main.f90`, end-of-step FDS outputs for the driver (`FDS_SETUP(MODE=3)`)

Apply after 0001 and 0002 (0003 to 0005 are independent). Touches only `Source/main.f90` (+25 lines, all inside `#ifdef WITH_AMREX`;
the preprocessor output with the macro undefined is unchanged). `git apply --check` passes on the tree at HEAD and in order after 0001 to 0005.

## What it does
- `FDS_SETUP(MODE,...)` gets a fourth mode. `MODE=3` runs the end-of-step output sequence of `MAIN_LOOP` (main.f90, "Apply velocity boundary conditions,
  and update values of HRR, DEVC, etc." to "Dump out diagnostics") and returns: `UPDATE_GLOBAL_OUTPUTS(T,DT,NM)` for the meshes of the process,
  `EXCHANGE_GLOBAL_OUTPUTS`, `UPDATE_CONTROLS`, `DUMP_GLOBAL_OUTPUTS`, `WRITE_STRINGS`, and `WRITE_DIAGNOSTICS` with the same `DIAGNOSTICS` rule as FDS
  (`MOD(ICYC,10**LO10)==0 .OR. MOD(ICYC,DIAGNOSTICS_INTERVAL)==0 .OR. T>=T_END`, `EXCHANGE_DIAGNOSTICS` first when more than one process).
  The slice, boundary-file, particle and restart dumps are not called (outside the M2a scope).
- `T`, `DT`, `ICYC` come from the module `FDS_AMREX_HOOKS` (`OUT_T`, `OUT_DT`, `OUT_ICYC`, set by `fds_hook_set_step`, `Source/driver/fds_amrex_hooks.f90`).
- Set-up sets `STEP_OUTPUTS_PATCHED = .TRUE.` before it returns. The driver calls `MODE=3` only when `fds_hook_step_outputs()` returns 1, so a tree without
  this patch still builds and runs (no FDS-format writers then).
- With the patch the driver writes through FDS's own routines: `CHID_devc.csv`, `CHID_hrr.csv`, `CHID_mass.csv`, `CHID_steps.csv`, `CHID_cpu.csv`, and the
  `.out` step blocks (plus the set-up and `END_FDS` content FDS writes anyway). The driver zeroes `Q_DOT`/`M_DOT` after `T = T + DT` (`fds_p_zero_dot`), as `MAIN_LOOP` does.

## Evidence
- `USE_AMREX=OFF`, patch 0006 applied on top of 0001 to 0005 (scratch tree `role1-s4-work/s4tree`, `tests/check_off_bitwise.sh`):
  `shunn3_32` 1 rank PASS (16 files bitwise identical to the baseline), `shunn3_4mesh_32` 4 ranks PASS (47 files).
- `USE_AMREX=ON`, gfortran 14.2: builds; `tests/run_outputs_check.sh` PASS (shunn3_32 `_mass.csv` bitwise to the baseline, `_hrr.csv` within 8e-10 of the column scale,
  `.out` pressure-iteration lines equal to the baseline); driver tests, decomposition check, kernel checks and the csmag_32 FISHPAK_BC comparison pass as before.

## Not done here
No oneAPI validation. Slice/boundary/particle files and restart files are not written by the driver.
