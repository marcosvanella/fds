# regrid_transport (Role 3): regrid and multi-level transport, Phase 3

Stub (first commit by the Chief Architect, 2026-10-02). Milestones R0 to R6 are in `notes/plan.md`.

| path | content |
|---|---|
| `OWNERS.md` | ownership and edit rules |
| `notes/` | plan and design notes |
| `tests/` | test scripts, run at 1 and 4 ranks, one thread |
| `CMakeLists.txt` | library `fds_regrid_transport_core` (parser) behind the interface target `fds_regrid_transport`; CTest registration; not yet added to the FDS build |
| `AmrInput.H/.cpp` | `&AMR` and `&AMR_REGION` parser and range checks (R0), no AMReX needed |
| `tests/test_amr_input.cpp` | parser unit tests (`ctest`) |
| `notes/readf90-cmake-patch-list.md` | patch list for `read.f90` and top-level CMake (Architect) |
