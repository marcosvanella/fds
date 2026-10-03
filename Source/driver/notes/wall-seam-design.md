# Wall-state seam (ADR-001 rulings W1 and W2): host-side part (S10.3)

Basis: ADR-001 "Wall tables and wall-state ownership" (v0.7), `drafts/gpu-staging-plan.md` stage 2, `drafts/gpu-generator-design.md` (wall lists).
There is no separate WP4 design note in the tree; this note is the Role 1 design of the host side and the list of what needs the GPU laptop.

## What is built (host only, `fds_wall_seam.f90`, `TimeLoop.cpp`)
**W1, wall lists.** One table set per mesh object in upstream order with the global `IW`; per box `WLIST_EXT = 1..N_EXTERNAL_WALL_CELLS` and
`WLIST_INT = N_EXTERNAL_WALL_CELLS+1..N_WALL_CELLS` (FDS orders the external walls first, `init.f90`). The driver has one FDS mesh object per box, so the lists are the identity ranges
and `NWL_EXT`, `NWL_INT` are the counts a generated wall kernel gets (`DO IWI=1,NWL; IW=WLIST(IWI)`). Fine-level boxes (patch 0007, `BOX_OBJ`) have zero walls: empty lists.
- C entries: `fds_wseam_counts(nm, &next, &nint)`, `fds_wseam_refresh(nm)` (rebuilds when the wall counts changed, returns 1 when it did: the signal for a device shim to re-gather its tables after
  `REASSIGN_WALL_CELLS` or a reallocation of `WALL`; the driver has no obstruction events, so today it rebuilds only on the first call).

**W2, `UVW_SAVE`, `U_GHOST`, `V_GHOST`, `W_GHOST`.** Producers (`MATCH_VELOCITY`, `ccib.f90`) stay on the host. `fds_wseam_upload(nm)` gathers the four arrays over `WLIST_EXT` into a staging array
`STAGE(4,NWL_EXT)` (the stand-in of the host-to-device copy) and asserts that the checksum of the staged copy equals the checksum of the host arrays (order-sensitive bitwise XOR/rotate over the
IEEE bit patterns). `fds_wseam_check(nm)` asserts later that the host arrays still match the staged copy, i.e. nothing wrote them behind the upload.

**Wiring** (`FDSTL_WSEAM=1`, off by default, no effect on results): staging + assertion at the start of the DENSITY stage and of WALL_BC of every level-0 stage, and the check at the end of both.
`FDSTL_WSEAM=2` is a negative control (a host value is perturbed after the upload; the run must stop with "checksum assertion of ADR-001 W2").
Finding: the DENSITY and WALL_BC kernels do not write the four arrays (the check passes in every stage of the test cases), which is what ruling W2 assumes.

## Tests
`tests/run_wall_seam_check.sh <build> [steps] [work]` on `dec1` (1 rank, 2176 external wall cells) and `dec4_np4` (4 ranks): `FDSTL_WSEAM=1` bitwise equal to the plain run (final fields and step log),
`FDSTL_WSEAM=2` stops with the message. The staging happens once per stage here; ADR-001 asks for once per step before the wall kernels. The per-stage form is the safe superset (MATCH_VELOCITY runs in
both stages); the shim may upload less often as long as `fds_wseam_check` stays true.

## Not covered here (needs the GPU laptop or the Intel/NVIDIA toolchain)
1. The real host-to-device transfer of the four arrays (`omp target update to(...)`, `has_device_addr` / `is_device_ptr` forms of ADR-001 K2 rule 7) and its timing; today the "device" is host memory.
2. The device-side wall kernels that take `WLIST`/`NWL` (generator work, owner: kernel generator) and the measured cost of the extra index gather (`[VERIFY]` on nvfortran in ADR-001).
3. The scatter/gather of host-written tables (`B1_RHO_D_DZDN_F` after a host `WALL_BC`, ADR-001 "Host-written tables") and the table refresh after `REASSIGN_WALL_CELLS` on a device copy.
4. The checksum assertion on the device side (the same function over the device copy); the host function here is the reference value.
5. Multi-box meshes per level and rank (one mesh object with several boxes, per-box `WLIST_*` that are not the identity): the driver has one mesh object per box, so the non-identity lists
   are written by the same code path but untested.

## Open question for the Architect
ADR-001 says "once per step before the wall kernels". The array contents change between the predictor and the corrector stage (MATCH_VELOCITY runs in both), so one upload per step would be
wrong unless the shim re-uploads after the second MATCH_VELOCITY; the host code uploads once per stage. Confirm "once per stage" as the ruling wording.
