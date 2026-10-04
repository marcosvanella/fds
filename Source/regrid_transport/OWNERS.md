# Source/regrid_transport: owner Role 3 (regrid and multi-level transport implementer)

`&AMR` input and its checks, the multi-level hierarchy (`AmrCore` subclass), tagging, regrid with conservative
transfer, coarse-to-fine ghost fill, average-down and the interface flux overwrite (D-050) for mass and species (Phase 3). Only Role 3 edits
this directory. Role 3 never edits `Source/driver/` or `Source/pressure_backend/`; changes needed there are
requests to their owners (Role 1, Role 2). C++ glue only (D-043, NFR-049); the only physics kernels allowed here
are tagging kernels, and a K1 tagging kernel records its reason. Edits to existing FDS files (`read.f90` for
`&AMR`, top-level CMake) go through the Chief Architect as a patch list, behind `#ifdef WITH_AMREX`.
Plan: `notes/plan.md`.
