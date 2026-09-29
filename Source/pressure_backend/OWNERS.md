# Source/pressure_backend: owner Role 2 (pressure backend interface implementer)

Solver-agnostic pressure interface and its backends (ADR-002 v1.1): MLMG and assembled-matrix HYPRE PCG+BoomerAMG,
the eps_H agreement check, masks, pins and the mean-removed right-hand side. Only Role 2 edits this directory.
Role 2 never edits `Source/driver/`. The interface header that the driver includes is `pressure_backend/PressureIface.H`
(Role 2 owns it; Role 1 reviews signature changes). Edits to existing FDS files go through the Chief Architect as a
patch list, behind `#ifdef WITH_AMREX`.
