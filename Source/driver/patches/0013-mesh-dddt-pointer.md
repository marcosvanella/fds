# Patch 0013 (DRAFT): `mesh.f90`, rebind DDDT in POINT_TO_BOX

**Status: DRAFT for the Architect. Applied in the Role 1 prove build; not committed (mesh.f90 is outside Source/driver).**

One line, inside the existing `#ifdef WITH_AMREX` body of `POINT_TO_BOX`, immediately after the `DS` rebind. Without it a fine box keeps the private `M%DDDT` from `POINT_TO_MESH_OBJECT` (allocated and left 0 in `fds_fine_level`), while `composite_fine_pressure_rhs` reads the AMReX `DDDT` MultiFab. `DIVERGENCE_PART_2` then writes the private array and the fine pressure right-hand side misses `(2*DP-DIV)/DT`.

The driver side (view code 23, `BOX_VIEW_TYPE%DDDT`, `FDS_HOOK_BIND_VIEW` CASE(23)) is in `Source/driver` and does not need this patch. The MultiFab is already allocated (`Fields` default list, `setVal(0)`); the Fortran kernel fills it once this pointer is live.

```
IF (ASSOCIATED(BV%DDDT)) DDDT => BV%DDDT
```

Apply on top of 0007. `USE_AMREX=OFF` does not preprocess this line.
