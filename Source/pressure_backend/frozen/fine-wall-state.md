# Fine burner face: parent wall state (fds_fine_walls_state)

Accessor in `Source/driver/fds_fine_level.f90`, for the fine-level burner (MASS_FLUX / HRRPUA OBST face) work: D terms (Role 2) and rho*Z / enthalpy boundary fluxes (Species).

```c
int fds_fine_walls_state(int lev, int ib, int pred, int nz, int nmax, int* n, int* ijkd /*[nmax][4]*/, double* s /*[nmax][5+2*nz]*/);
```

- `lev` fine level (>= 1), `ib` = `mfi.index()+1` of the level's box (same numbering as `fds_fine_walls_build`).
- `pred` 1 predictor (U_NORMAL_S), 0 corrector (U_NORMAL). `nz` = N_TRACKED_SPECIES (else return 1). `nmax` capacity (return 1 if exceeded).
- Returns 0 on success; `n` = faces listed.

Per face `w` (Fortran column order, C row `w`):

| field | meaning | units |
|---|---|---|
| `ijkd[4w+0..2]` | fine gas cell I,J,K, box-local, 1-based (global = box smallEnd + I-1) | - |
| `ijkd[4w+3]` | D = +-1..3: the face of that gas cell that is the wall (IOR = -D). Face index along \|D\|: i for D<0, i+1 for D>0 (AMReX nodal) | - |
| `s[0]` | face normal velocity in +axis direction (-SIGN(IOR)*U_NORMAL(_S)); D<0 burner inflow => >0 | m/s |
| `s[1]` | fine face area | m2 |
| `s[2]` | fine area / parent wall AREA (area weight) | - |
| `s[3]` | RHO_F | kg/m3 |
| `s[4]` | TMP_F (wall T; enthalpy: caller evaluates h_s(TMP_F, ZZ_F), as FDS) | K |
| `s[5 .. 4+nz]` | ZZ_F(1:nz) | - |
| `s[5+nz .. 4+2nz]` | M_DOT_G_PP_ADJUST(1:nz), > 0 into the gas | kg/m2/s |

Selection: hidden fine walls whose level-0 parent is an internal (OBST) SOLID wall with nonzero U_NORMAL(_S) or M_DOT_G_PP_ADJUST.
Domain-boundary walls and inert OBST faces are not listed.

Timing: nothing is stored; the values are the parent L0 wall's at call time. Call per stage after the level-0 wall BCs of that stage
(`s_wall_bc(pred)` in TimeLoop), i.e. before DIVERGENCE_PART_1 of the predictor (pred=1) and of the corrector (pred=0). Both stages.

Coverage: fine faces only. L0 wall faces under fine boxes are not listed; L0 keeps FDS's own wall loops (level 0 is unchanged), and the
covered L0 cells are overwritten by average-down.

No WALL_TYPE / ONE_D / BOUNDARY_* changes: read-only on the parent's BOUNDARY_PROP1.
