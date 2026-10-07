# Two-level plotfile ParaView smoke (Integration)

Show-and-tell helper for R=2 cases `ns2d_16_l0` and `blob_16_l0`. Role 1 owns `WriteMultiLevelPlotfile` in the driver; Integration does not edit `src/`. Prefer AMReX plotfile (CLI `--plotfile`, default off, RHO first) over HDF5; FDS VTKHDF deferred.

Scripts and scout table: `prototypes/s4_cuda_mass/stage1/plotfile_smoke/` (`README.md`, `check_plotfile.sh`, `paraview_smoke.py`). Write hook for Role 1: `Source/driver/TwoLevelRun.cpp` after `loop.advance()` or end-of-run — not inside `TimeLoop::advance`. Field: `reg.fields(lev)["RHO"]`.
