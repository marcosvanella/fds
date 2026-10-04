# Fine-mesh VELOCITY_BC overwrote the level-0 periodic U/W ghosts: cause, reproducer, fix (S14.4)

Hypothesis (confirmed by the reproducer below): the FDS routine was not misbehaving on the fine box; it was not run on the fine box. `TimeLoop::bind_level` builds `BcStep(Lv, F)` from the registry `Level` of the
fine level. That `Level` has `fds_mesh_offset = 0` (the number is only known after `fds_fine_level_create` returns NM0), and `BcStep` calls `fds_p_save_uvw`, `fds_g_fill_om`, `fds_g_velocity_bc`,
`fds_g_viscosity_bc` ... with `box index + 1 + m_l0.fds_mesh_offset`. Fine box i therefore reached FDS mesh number i + 1, which is a LEVEL-0 mesh object (`MESHES(i+1)`), not `FINE_LEVEL(1)%BOX(i+1)`
(mesh number NM0 + i + 1). With one fine box and one level-0 mesh, `VELOCITY_BC` of "mesh 1" ran on the level-0 mesh after the fine level's `fill_omesh` had written fine-level data into the exchange buffers: it
rewrote the periodic ghost layers of level 0 U (z rows) and W (x columns), not V (one cell in y, no periodic neighbour in y). With several ranks the same slip reaches a mesh object that is not local to the
rank (`MESHES(3)` on rank 1): `fds_p_save_uvw` segfaulted (4 ranks, 16 fine boxes). With 4 fine boxes on 4 ranks the number happened to be local, so it only corrupted data.

Reproducer (small, no Fortran change): `tests/run_stage_boundary_fine_check.sh <build>`: `fds_amr ns2d_16_l0.fds --two-level-run --stage-boundary-test` calls `stage_boundary(1,3)` and `(1,6)` on the bound
fine level and compares level 0 U/V/W with ghost layers before and after (`FDSTL_SKIPAFT=0`). With `FDSTL_LEGACY_BCSTEP_OFFSET=1` (offset 0, the old behaviour) level 0 U and W change (6.7e-16 after the fine `fill_omesh` was dropped, see below; the original corruption, with the fine level's `fill_omesh` also run on
level-0 objects, was 1.57 in U and W, V untouched); with the fix all three are unchanged bit for bit. The run in the same script with `legacy` restores the pre-S14.1 guard of `fds_p_save_uvw` (level-0 numbers only) and aborts on the fine box with the guard message:
that is the case that triggers the `stage_boundary(1,3|6)` abort the S14.1 guard removed.

Fix: `BcStep::set_mesh_offset(nm0)` called by `bind_level`. `FDSTL_SKIPAFT` now defaults to 0 (nothing skipped on a fine level); 40-step results with it are identical to the skip (mass change 1.3e-14, max |div u - D| 1.8e-12).

Second finding (S14.4): `fill_omesh` must not run on a fine level at all. It copies the FAB of every box into `MESHES(NM)%OMESH(NOM)`, and a fine box has no OMESH (no neighbour mesh objects: its same-level ghosts
come from the AMReX fill, the coarse-fine ghosts from the composite fill). With the correct fine mesh numbers it indexed `MESHES` beyond its size (4 ranks, 16 fine boxes: segfault in `fds_g_fill_om`).
`BcStep::after_exchange` runs it on level 0 only; `fds_g_fill_om` returns 1 for a mesh number above NMESHES. `VELOCITY_BC` (patch 0008: `POINT_TO_MESH` is `POINT_TO_BOX` in velo.f90) runs on the fine box
through `FINE_LEVEL(L)%BOX`; a fine box has no wall cells, so its wall loops are empty.
