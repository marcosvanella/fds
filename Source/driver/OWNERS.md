# Source/driver: owner Role 1 (data layout implementer)

C++ AmrCore driver and the Fortran `bind(C)` shim (ADR-001 Option A/C, NFR-049: C++ only for driver and glue).
Only Role 1 edits this directory. Role 2 (pressure backend) never edits it; the driver calls the backend through
`PressureIface.H`, whose signature is agreed between Role 1 and Role 2 and recorded in `README.md` here.
Edits to existing FDS files (`CMakeLists.txt`, `Source/main.f90`, `mesh.f90`, `velo.f90`, `wall.f90`, `pres.f90`,
`init.f90`, `read.f90`) are proposed as a patch list and applied by the Chief Architect, all behind
`#ifdef WITH_AMREX` so that `USE_AMREX=OFF` stays identical (IR-006).
