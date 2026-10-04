# Role 3 progress notes (one line per step)

- Driver rebuilt from the branch (gfortran): E1-E5 on 1 and 4 ranks pass; E3b ON 6e-16 (1 rank) / 1e-16 (4 ranks), E1 RHO and ZZ 3e-16; DriverModes injects the parent stage arrays on a new fine level (inject_parent_stage_arrays) and checks D and DS in the FR-016(a) ghost gate.
- Tag-buffer / proper-nesting assertions on 64^3 in 8 boxes (3 levels ratio 2, 2 levels ratio 4): ctest regrid_transport_tag_buffer_large (1 and 4 ranks) pass; buffer complete inside every box and along axes, 1993/28285152 and 279/1332000 cube edge/corner cells across level-l box boundaries missing (one box: 0); nesting and grid sizes clean.
- run_e2e_driver.sh: MASS_TOL and E1_ZZ_TOL 1e-13 (observed 6e-16 / 3e-16), E1 negative control (another blob smoothing fails the same comparison), E3b OFF control 4.4e-5, E4 at steps 2 and 3 incl. D and DS (D-059 hole-corner cells listed: KRES 4, D 2, DS 2 per step), E6 pending patch 0010 guard.
- Converter refuses &TRNX/&TRNY/&TRNZ and TRN?_ID in AMR mode (an &AMR line present); N-5b input and the stand-alone groups refused, 69 checks pass; N-5b was accepted before.
