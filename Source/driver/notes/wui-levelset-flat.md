# Flat WUI level-set cases: stock FDS vs fds_amr (unrefined and refined)

Tree: FDS-AMReX 8c16d56d50 (level set on L0 through 4909f15828). Role 1's 87d8513887 (PHI on every level) and 05a1624507 (--tag-phi) are NOT in this
tree (not on the box), so refined runs were not possible here; the table says what this tree does. Stock FDS: work-species/build_fds/fds.
Script: `tests/run_wui_levelset_compare.sh STOCK_FDS FDS_AMR OUTDIR [CASES]` (REFINED=1 PATCH="i0 i1 j0 j1" adds the refined variant, which needs a build with --tag-phi);
`tests/wui_levelset_compare.py` prints final HRR and the max relative difference of every device. The LS* cases have 50x50x5 or 20^3 cells, so the script adds
`&AMR BLOCKING_FACTOR=1` to the fds_amr copy (default 8 does not divide the domain).

| case (Verification/WUI) | mode | unrefined fds_amr --run vs stock | refined |
|---|---|---|---|
| level_set_fuel_model_1 (40x40x16, 36000 s, 3906 steps) | 4 | runs; device `mass` final 132599.6 vs 132603.6 (3e-5); front: mass>10% at 7997 vs 7993 s, HRR>1 kW at 40 s both; HRR mean 55045.7 vs 55056.4 (0.02%), max 140250 vs 138892 (+1%); final sample 159.7 vs 165.6 (tail, noise level) | static patch 14..26 (no --tag-phi): ERROR STOP in POINT_TO_MESH (mesh 2 is a fine box; level-set code on a fine level, not in this tree) |
| LS_ellipse_{0,5}ms_{00,30}deg, LS_wind_ramp_lin/quad (mode 1, no REAC) | 1 | refused: "case outside the S5 scope: SOLID_PHASE_ONLY" (rc 6) | not run |
| LS4_ember_yield (mode 4 + particles) | 4 | refused: "particle drag Lagrangian particles" (rc 6) | not run |
| LS4_ember_ignition | 4 | rc 1 (not investigated, particle case) | not run |

So the only flat, driver-supported level-set case is level_set_fuel_model_1; the others need the SOLID_PHASE_ONLY / Lagrangian-particle scope (S6).
Next step once Role 1's commits are on the box: `REFINED=1 PATCH="14 26 14 26" tests/run_wui_levelset_compare.sh ... level_set_fuel_model_1` (expect the known mass-gate failure until fine-level coupling).
