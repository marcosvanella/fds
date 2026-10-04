# Which routine rewrites FVX/FVY/FVZ between VELOCITY_FLUX and VELOCITY_CORRECTOR (corrector stage)

Dump brackets (`FDSTL_STAGEG=<cycle>`, tags in README S13.1), case `dec4_np4` (periodic, 4 boxes, no obstructions), cycle 2, corrector:

| step | FVX changed | FVY changed | FVZ changed |
|---|---|---|---|
| after VELOCITY_FLUX + WALL_BC -> after DIVERGENCE_PART_1 -> after DIVERGENCE_PART_2 -> start of the pressure scheme | 0 | 0 | 0 |
| pressure scheme: `MATCH_VELOCITY_FLUX` (after `MESH_EXCHANGE(5)`, first iteration) | 288 values, max change 1.32 | 0 | 288 values, max change 1.32 |
| pressure scheme: `NO_FLUX` (+ `PRESSURE_SOLVER_COMPUTE_RHS`) | 0 | 0 | 0 |
| pressure solve, `VELOCITY_ERROR` etc. -> `c_prevcorr` | 0 | 0 | 0 |

So in this case FVX and FVZ are rewritten only by `MATCH_VELOCITY_FLUX` (velo.f90, driver call `fds_g_match_flux`) at the faces on the box interfaces (each such face gets the mean of its own value and the area-weighted FVX of the neighbour mesh's faces; with one mesh the routine returns at once, NMESHES==1). The other writers of FVX/FVY/FVZ inside the pressure scheme exist in the FDS source but are inert here: `NO_FLUX` (velo.f90, prescribes FVX/FVY/FVZ on faces
of obstructions, at solid boundaries toward the wall normal velocity, and on interface faces, every pressure iteration), `BAROCLINIC_CORRECTION` (baroclinic torque, only with BAROCLINIC) and `CC_NO_FLUX` (cut cells). A case with
obstructions or prescribed normal velocities changes FVX/FVY/FVZ in `NO_FLUX` as well; the brackets `c_fv_matched` -> `c_fv_noflux` show that.
