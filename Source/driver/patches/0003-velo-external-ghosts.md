# Patch 0003: `velo.f90`, box-boundary ghost cells "filled externally" (guarded, `WITH_AMREX`)

Apply after 0001 and 0002 (independent of 0004 and 0005 except that the module `FDS_AMREX_HOOKS`, `Source/driver/fds_amrex_hooks.f90`, is part of the
`USE_AMREX=ON` build). Touches only `Source/velo.f90` (+18 lines, every added line inside `#ifdef WITH_AMREX`, so with the macro undefined the
preprocessor output is byte-identical to the old source).

## What it does
`USE FDS_AMREX_HOOKS, ONLY: EXTERNAL_GHOSTS_FILLED` (a module logical, default `.FALSE.`, set by `fds_hook_set_flag`). When it is `.TRUE.` the
`EXTERNAL_WALL(IW)%NOM>0` branches that write ghost cells of a box boundary from `OMESH(NOM)` are skipped, because the AMReX level ghost fill has written those
cells (ADR-001 Option C pitfall table):

| Routine (line in the unpatched file) | Skipped when the flag is set | Ghost values that are then taken from AMReX |
|---|---|---|
| `VISCOSITY_BC` (~517, `WALL_LOOP`) | the whole `NOM>0` body | MU, KRES, D (DS for the estimated stage) |
| `NO_FLUX` (~1379) | the `HP` fill from `OMESH(NOM)%H/HS` | H/HS ghost (pressure step, Role 2) |
| `VELOCITY_BC` (~1861, first `WALL_LOOP`) | the normal UU/VV/WW ghost face from `OMESH` | U/V/W, US/VS/WS normal ghost face |
| `MATCH_VELOCITY` (~2620, after the CC_IBM branch) | return | shared-face match: done by the driver on AMReX data |
| `MATCH_VELOCITY_FLUX` (~2865, after the CC_IBM branch) | return | same for the flux (multi-box pressure iteration only) |

Not changed: the `EDGE_LOOP` interpolation branch of `VELOCITY_BC` (`INTERPOLATION_IF ... ELSE`, ~2404-2514). Its `OM%U/V/W` reads at edges are interface
edges of the periodic / same-level case, which the driver fills the same way as the faces (S5: edge strip of the level fill, `NOM(ICD)` from the same table);
it is left unguarded on purpose until the S5 step-loop test shows the edge data of the AMReX fill is enough (the BC-chain check of S4 reproduces
FDS's ghost values bitwise with OMESH filled, so the branch is correct today with the flag `.FALSE.`).

## How S4 uses it
S4 keeps the flag `.FALSE.`: the driver (`GhostExchange.cpp`, `fds_ghost_bc.f90`) fills `OMESH(NOM)` of each local box with the block copy that
`MESH_EXCHANGE` does and runs the UNMODIFIED `VISCOSITY_BC`, `VELOCITY_BC`, `MATCH_VELOCITY`; the values are FDS's own by construction (`BCCHAIN_P/C` are
bitwise). The flag TRUE path (no OMESH fill, no O(domain) broadcast) is exercised from S5; until then this patch only adds the hook and is inert.
Cost of the S4 route: one broadcast of each box's FAB to every rank that owns a neighbour, and the OMESH arrays kept allocated by FDS set-up.

## Evidence
- `git apply --check` passes on the tree at HEAD and with 0001, 0002, 0004, 0005 applied in order (scratch tree
  `role1-s4-work/s4tree`, built out of tree).
- `tests/check_off_bitwise.sh` with 0003+0004+0005 applied, `USE_AMREX=OFF`: `shunn3_32` 1 rank PASS (16 files), `shunn3_4mesh_32` 4 ranks PASS (47 files),
  all bitwise identical to the baseline (same run for the three patches; the OFF build compiles the unchanged text).
- `USE_AMREX=ON` build of the patched scratch tree with the flag left `.FALSE.`: builds; `tests/run_kernelcheck.sh` (dump, full+bc, face+bc bitwise; plain
  full and face: BCCHAIN bitwise, other tags within the gate) and `tests/run_driver_tests.sh` ALL DRIVER TESTS PASS (threads=1, 1 and 4 ranks).
- gfortran 14.2 only. `ifx`/`ifort` are not available on this box, so the oneAPI build of the patched text is NOT checked (the added Fortran is plain
  `USE ... ONLY`, `IF (...) CYCLE/RETURN`, nothing compiler specific).

## Not done here
- The flag TRUE path is not exercised (S5). The `EDGE_LOOP` branch and the `CC_IBM` routines are untouched. `MATCH_VELOCITY` with the flag set assumes the driver
  does the shared-face match; the driver does not yet (S5).

Note (S6b): "shared-face match: done by the driver on AMReX data" now covers the periodic domain faces as well (`BcStep::match_periodic_faces`, `0.5*(a + ((b*dA1)*dA2)/(dA1*dA2))`); before S6b only box interfaces were handled, which left the two copies of a periodic flow face unaveraged in a fully periodic 3-D case (csmag_32). No change to the patch itself.
