# Patch 0014 (DRAFT): `init.f90`, zero wall-normal velocity on INERT/MIRROR in `UVW_INIT`

**Status: DRAFT for the Architect. Applied in the Role 1 prove build; not committed (`init.f90` is outside Source/driver).**

`UVW_INIT` (around the `W(:,:,0) = W(:,:,KBAR)` copy) loads the UVW file through `K=KBAR` and, for `PERIODIC_TEST==2`, copies that top face onto `W(:,:,0)`. Neither INERT nor MIRROR is set to 0, so the mesh has a nonzero net volume flux. The singular Neumann projection removes the mean (`removed_rel ~ 1`) and the continuity update then grows mass at a constant rate (0.1059 kg/s on the 48³ vortex ring; the uniform div floor 0.0909729 is the mean top-face W).

The added loop, under `#ifdef WITH_AMREX`, runs after that copy and before `US=U` / `U_NORMAL`, and sets the normal velocity to 0 on external `SOLID_BOUNDARY` and `MIRROR_BOUNDARY` faces. Divergence, wall `U_NORMAL` and the initial pressure are then built from a closed box. Open and periodic faces are unchanged.

The driver also zeros those faces again before the initial projection (`TimeLoop::zero_closed_wall_normal`) and shifts `D`/`DS` by the divergence change, so a face that is still nonzero after init does not reintroduce the mean.
