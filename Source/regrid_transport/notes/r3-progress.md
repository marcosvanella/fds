# Role 3 progress notes (one line per step)

- Driver rebuilt from the branch (gfortran): E1-E5 on 1 and 4 ranks pass; E3b ON 6e-16 (1 rank) / 1e-16 (4 ranks), E1 RHO and ZZ 3e-16; DriverModes injects the parent stage arrays on a new fine level (inject_parent_stage_arrays) and checks D and DS in the FR-016(a) ghost gate.
- Tag-buffer / proper-nesting assertions on 64^3 in 8 boxes (3 levels ratio 2, 2 levels ratio 4): ctest regrid_transport_tag_buffer_large (1 and 4 ranks) pass; buffer complete inside every box and along axes, 1993/28285152 and 279/1332000 cube edge/corner cells across level-l box boundaries missing (one box: 0); nesting and grid sizes clean.
