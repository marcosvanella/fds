# Upstream FireX 1fb5c300 merge: impact on the loop work list and the generator markers

Scope: the merge commit `8ad685f435` ("Merge FireX (1fb5c300fa) into FDS-AMReX") of the FDS-AMReX tree. Its first parent `9dc218d55b` is the pre-merge FDS-AMReX commit and its second parent is FireX `1fb5c300fa`. Every range below was produced from `git show` / `git diff -U0` of `9dc218d55b` and `8ad685f435` (Source/ only). The upstream merge report for 1fb5c300fa in `inventory/upstream_merge_reports/` ran its trial merge on `e8a9182835`, which is not an ancestor of `8ad685f435`; the real pre-merge commit is `9dc218d55b`.

Files changed by the merge (`git diff --stat 9dc218d55b 8ad685f435 -- Source/`): `dump.f90`, `func.f90`, `init.f90`, `main.f90`, `type.f90`, `wall.f90` (197 insertions, 41 deletions). `pres.f90`, `radi.f90`, `divg.f90`, `mass.f90`, `velo.f90`, `turb.f90`, `soot.f90`, `fire.f90`, `part.f90`, `vege.f90` are unchanged by the merge.

## Line coordinates in use

| Item | Line numbers refer to | Where stated |
|---|---|---|
| `amrex/loop_work_list.csv` `file_lines` | survey source FireX `36975d7` | `tools/inventory/loop_work_list.py:5`, `amrex/loop-work-list.md` section 1 |
| `s5_markers.toml` `lines`, `anchor_line` (s5-gen) | `36975d7` | `[policy] upstream_commit = "36975d7"` (`s5_markers.toml:12` on s5-gen) |
| `markers/rad_kernels.toml` (s5-gen) | FDS-AMReX `bee11f0329` | `markers/rad_kernels.toml:1` on s5-gen |
| `markers/golden_signatures.json` (s5-gen) | no line numbers: `source` is `file:ROUTINE`, plus `body_sha` | |

So "old range" below is the range at the pre-merge commit `9dc218d55b`, "new range" the range at `8ad685f435`, and the survey range (`36975d7`) is given as well because the work list and the markers still use it.

## 1. Affected work-list loops

A work-list loop counts as affected when its range moves between `9dc218d55b` and `8ad685f435`, when its text changes (whitespace-insensitive compare of the whole range), or when a hunk falls inside it. Five of the 231 loops qualify; all five are in `func.f90` or `init.f90`. No loop of `wall.f90` moves: the only `wall.f90` hunk is the one-line change at `wall.f90:2181` in `SOLID_HEAT_TRANSFER` (`B1%Q_CONDENSE` becomes `B1_BACK%Q_CONDENSE`), which no surveyed loop spans.

| Key | Routine | Survey range (36975d7) | Old range (9dc218d55b) | New range (8ad685f435) | Body changed | Translated | Owner (loop_claims.csv) | Action |
|---|---|---|---|---|---|---|---|---|
| L0630 | BLOCK_CELL | func.f90:5509-5524 | func.f90:5509-5524 | func.f90:5581-5627 | y | n | none | hold |
| L0823 | REASSIGN_WALL_CELLS | init.f90:4940-5054 | init.f90:4940-5054 | init.f90:4951-5065 | n | n | AMR Solid Phase Lead | re-pin only |
| L0822 | REASSIGN_WALL_CELLS | init.f90:4917-4936 | init.f90:4917-4936 | init.f90:4928-4947 | n | n | none | re-pin only |
| L0817 | CREATE_OR_REMOVE_OBST | init.f90:4878-4889 | init.f90:4878-4889 | init.f90:4889-4900 | n | n | none | re-pin only |
| L0677 | PACK_CELL | func.f90:5183-5201 | func.f90:5183-5201 | func.f90:5248-5266 | n | n | none | re-pin only |

Notes per row: L0630: body changed by the new single-cell branch (`func.f90:5589-5619`). L0823, L0822, L0817: +11 lines from the `init.f90` hunks above them (`init.f90:1625`, `1817-1818`, `3711-3724`, `4550`, `4698-4712` at `8ad685f435`). L0677: +65 lines from `func.f90:4036-4044`, `4523-4525` and the new subroutine `PACK_WALL_LITE` (`func.f90:4549-4595`).

Counts: **re-pin only 4** (L0677, L0817, L0822, L0823), **re-translate 0** (no translated loop or kernel is in a changed range), **hold 1** (L0630).

Claimed loops whose owner should be told: **L0823 (AMR Solid Phase Lead)**: its lines move by +11 (`init.f90:4951-5065` at `8ad685f435`); the text is unchanged, so this is a re-pin only. No claimed loop has a changed body, so no claim was touched.

## 2. BLOCK_CELL (L0630): what changed, and why it is on hold

`func.f90:5568-5629` at `8ad685f435` (`func.f90:5503-5526` at `9dc218d55b`). The loop L0630 is `func.f90:5581-5627` (was `5509-5524`).
- New declarations: `USE GLOBAL_CONSTANTS,ONLY: RSUM0`, `INTENT(IN)` on all arguments, locals `COUNT`, `TMP_AVG`, `SINGLE_CELL` (`func.f90:5570-5574`).
- `SINGLE_CELL` is set before the nest when the development machine is one cell and `MESHES(NM)%OBSTRUCTION(OBST_INDEX)%CONSUMABLE` (`func.f90:5577-5578`).
- Inside `CASE(0)` a new branch (`func.f90:5589-5619`) reads `CELL(IC)%WALL_INDEX(-3..3)` and the six `TMP` neighbours, then writes `TMP`, `RHO` (from `PBAR_S(K,PRESSURE_ZONE(I,J,K))` and `RSUM0`) and `RHOS`.
- Points for whoever ports it later:
  - `OBSTRUCTION(0)%CONSUMABLE` is read for the single-cell calls with `OBST_INDEX=0` (`read.f90:15089`, `15101`, `15113`). That is legal because `OBSTRUCTION` is allocated from 0 (`read.f90:10860`).
  - `TMP_AVG/REAL(COUNT,EB)` (`func.f90:5616`) divides by zero when all six `WALL_INDEX` entries are non-zero.
  - In `CREATE_OR_REMOVE_OBST`, `BLOCK_CELL` is called at `init.f90:4885`, and loop L0817 (`init.f90:4889-4900`) then sets `RHOS` and `RHO` of the same cells to `RHO_0(K)` (`init.f90:4892-4893`). On that path only the new `TMP` value survives, and with `SOLID_PHASE_ONLY` not even that one (`init.f90:4894`). L0630 and L0817 should be translated together, in this order.
- The upstream merge report classes two of the three BLOCK_CELL hunks as surface-changing. The kernel argument list grows by `WALL_INDEX`, `TMP`, `RHO`, `RHOS`, `PBAR_S`, `PRESSURE_ZONE`, `RSUM0` and the per-obstruction `CONSUMABLE` flag.

## 3. Surveyed loops outside the work list that the merge touches

| Survey key | Routine | Survey class | Old range (9dc218d55b) | New range (8ad685f435) | Body changed |
|---|---|---|---|---|---|
| L0821 | OPEN_AND_CLOSE (`VENT_LOOP`) | host I/O (not on the work list) | init.f90:4731-4835 | init.f90:4742-4846 | n |
| L0848 | MESH_EXCHANGE | exchange, retired (not on the work list) | main.f90:3217-3687 | main.f90:3230-3713 | y |
| L0849 | MESH_EXCHANGE | exchange, retired (not on the work list) | main.f90:3767-4043 | main.f90:3793-4077 | y |

The `OPEN_AND_CLOSE` hunks are in `OBST_LOOP` (`init.f90:4550`, `4698-4712`), not in the surveyed `VENT_LOOP`. `GAS_PHASE_OUTPUT` (`dump.f90:8604-10073`) and the `WALL_TYPE` counters (`type.f90:454-456`) have no surveyed loop. `MESH_EXCHANGE` (`main.f90:3203-4082`) is retired from the work list (exchange class), so its new codes 66/67 do not reach any key.

## 4. New loops added upstream in the changed routines

None. No `DO` statement is added or removed anywhere in the merge diff (`git diff -U0 9dc218d55b 8ad685f435 -- Source/`, added and removed lines searched for `DO`). The new subroutine `PACK_WALL_LITE` (`func.f90:4549-4595`) has only array-section `EQUATE` calls and no loop. In any case `loop_work_list.py` assigns no keys: keys come from the survey CSV `inventory/gpu_generator_loop_classes.csv` (`tools/inventory/loop_work_list.py:275-277`). A new loop would be listed as "new, unkeyed"; there are none.

## 5. Generator markers and golden: re-pin list

**Caused by this merge: nothing.** None of the 113 marker entries that carry `lines` moves between `9dc218d55b` and `8ad685f435`. The `func.f90` kernels (`gsfv_*`, `func.f90:1348-1453`) lie above the first `func.f90` hunk (`func.f90:4036` at `9dc218d55b`). `radi.f90` is identical at `bee11f0329` and `8ad685f435` (no hunks). The hand-written sweep `s5_rad_sweep/s5rad_sweep.toml` (pinned at `bee11f0329`) is therefore unaffected too. The golden (`markers/golden_signatures.json`) has no line numbers, and every marked range has the same text (whitespace-insensitive) at its pin and at `8ad685f435`, so no golden entry changes.

**For the planned re-pin of the markers to `8ad685f435`** (done on s5-gen after the merge lands), the exact ranges are listed below. 83 of the 113 entries move, all by drift from earlier FDS-AMReX commits (branch patches and the `bee11f0329` merge), none by this merge. The other 30 keep their numbers. The body text is identical in all 113. The 14 wall-list variants (`wlist_of`) carry no lines and need no change. The entries copied in `markers/gsfv_field_kernels.toml`, `nx_kernels.toml`, `pvf_kernels.toml`, `soot_kernels.toml` and `species_kernels.toml` have the same lines as their `s5_markers.toml` copies and have to be changed with them.

Caution for the work list: `loop_work_list.py` matches kernel `lines` against the survey ranges at `36975d7` (`tools/inventory/loop_work_list.py:313-329`), with one hand shift for `rad_wall_qin_zero` (`:304`). Re-pinning the markers alone would detach kernels from their loops. Either the survey ranges are re-pinned in the same step, or the script gets a pin-to-survey line map (section 6 of `loop-handout-round-next.md` gives the work-list ranges at `8ad685f435`).

| Kernel | Sidecar | File | Pinned lines (pin) | Anchor line | Lines at 9dc218d55b | **Lines at 8ad685f435** | **Anchor at 8ad685f435** | Moves |
|---|---|---|---|---|---|---|---|---|
| rho_d_interp | s5_markers.toml | divg.f90 | 133-140 (36975d7) | 133 | 142-149 | **142-149** | **142** | +9 |
| d_z_max | s5_markers.toml | divg.f90 | 159-165 (36975d7) | 159 | 168-174 | **168-174** | **168** | +9 |
| rho_d_dzd | s5_markers.toml | divg.f90 | 171-182 (36975d7) | 174 | 180-191 | **180-191** | **183** | +9 |
| wall_rho_d_dzdn | s5_markers.toml | divg.f90 | 192-232 (36975d7) | 192 | 201-241 | **201-241** | **201** | +9 |
| rho_d_maxloc_fix | s5_markers.toml | divg.f90 | 245-258 (36975d7) | 245 | 254-267 | **254-267** | **254** | +9 |
| h_rho_d_dzd | s5_markers.toml | divg.f90 | 298-312 (36975d7) | 301 | 313-327 | **313-327** | **316** | +15 |
| dp_div_heat | s5_markers.toml | divg.f90 | 390-402 (36975d7) | 390 | 405-417 | **405-417** | **405** | +15 |
| del_rho_d_del_z | s5_markers.toml | divg.f90 | 408-416 (36975d7) | 408 | 423-431 | **423-431** | **423** | +15 |
| cp_rhg | s5_markers.toml | divg.f90 | 435-443 (36975d7) | 438 | 450-458 | **450-458** | **453** | +15 |
| conductivity | s5_markers.toml | divg.f90 | 462-471 (36975d7) | 462 | 477-486 | **477-486** | **477** | +15 |
| wall_kp_ghost | s5_markers.toml | divg.f90 | 482-486 (36975d7) | 482 | 497-501 | **497-501** | **497** | +15 |
| kdtd | s5_markers.toml | divg.f90 | 512-523 (36975d7) | 515 | 527-538 | **527-538** | **530** | +15 |
| wall_corr_kdtd | s5_markers.toml | divg.f90 | 532-554 (36975d7) | 532 | 547-569 | **547-569** | **547** | +15 |
| dp_kdtd | s5_markers.toml | divg.f90 | 561-570 (36975d7) | 564 | 576-585 | **576-585** | **579** | +15 |
| dp_species | s5_markers.toml | divg.f90 | 646-655 (36975d7) | 649 | 661-670 | **661-670** | **664** | +15 |
| rho_h_s_p | s5_markers.toml | divg.f90 | 814-822 (36975d7) | 814 | 829-837 | **829-837** | **829** | +15 |
| cell_adv_hs | s5_markers.toml | divg.f90 | 946-966 (36975d7) | 946 | 961-981 | **961-981** | **961** | +15 |
| rho_z_p_divg | s5_markers.toml | divg.f90 | 998-1004 (36975d7) | 998 | 1013-1019 | **1013-1019** | **1013** | +15 |
| spec_wall_zz | s5_markers.toml | divg.f90 | 1021-1084 (36975d7) | 1021 | 1036-1099 | **1036-1099** | **1036** | +15 |
| rho_rmw_divg | s5_markers.toml | divg.f90 | 1097-1105 (36975d7) | 1097 | 1112-1120 | **1112-1120** | **1112** | +15 |
| spec_wall_rmw | s5_markers.toml | divg.f90 | 1120-1183 (36975d7) | 1120 | 1135-1198 | **1135-1198** | **1135** | +15 |
| flux_mw_fix_zz | s5_markers.toml | divg.f90 | 1191-1213 (36975d7) | 1191 | 1206-1228 | **1206-1228** | **1206** | +15 |
| wall_spec_adv2 | s5_markers.toml | divg.f90 | 1237-1266 (36975d7) | 1237 | 1252-1281 | **1252-1281** | **1252** | +15 |
| cell_adv_zz | s5_markers.toml | divg.f90 | 1269-1289 (36975d7) | 1269 | 1284-1304 | **1284-1304** | **1284** | +15 |
| wall_connect_zones | s5_markers.toml | divg.f90 | 1301-1319 (36975d7) | 1301 | 1316-1334 | **1316-1334** | **1316** | +15 |
| wall_bc_dp | s5_markers.toml | divg.f90 | 1574-1604 (36975d7) | 1574 | 1589-1619 | **1589-1619** | **1589** | +15 |
| div2_pred | s5_markers.toml | divg.f90 | 1611-1618 (36975d7) | 1611 | 1626-1633 | **1626-1633** | **1626** | +15 |
| div2_corr | s5_markers.toml | divg.f90 | 1621-1630 (36975d7) | 1621 | 1636-1645 | **1636-1645** | **1636** | +15 |
| div_extrema | s5_markers.toml | divg.f90 | 1675-1711 (36975d7) | 1675 | 1690-1726 | **1690-1726** | **1690** | +15 |
| wall_q_condense_zero | s5_markers.toml | fire.f90 | 1957-1961 (36975d7) | 1957 | 1957-1961 | **1957-1961** | **1957** | no |
| gsfv_central | s5_markers.toml | func.f90 | 1348-1354 (36975d7) | 1348 | 1348-1354 | **1348-1354** | **1348** | no |
| gsfv_godunov | s5_markers.toml | func.f90 | 1358-1368 (36975d7) | 1358 | 1358-1368 | **1358-1368** | **1358** | no |
| gsfv_superbee | s5_markers.toml | func.f90 | 1372-1389 (36975d7) | 1372 | 1372-1389 | **1372-1389** | **1372** | no |
| gsfv_minmod | s5_markers.toml | func.f90 | 1393-1410 (36975d7) | 1393 | 1393-1410 | **1393-1410** | **1393** | no |
| gsfv_charm | s5_markers.toml | func.f90 | 1414-1435 (36975d7) | 1414 | 1414-1435 | **1414-1435** | **1414** | no |
| gsfv_mp5 | s5_markers.toml | func.f90 | 1439-1453 (36975d7) | 1439 | 1439-1453 | **1439-1453** | **1439** | no |
| rho_z_p_mass | s5_markers.toml | mass.f90 | 70-76 (36975d7) | 70 | 76-82 | **76-82** | **76** | +6 |
| mass_wall_zz | s5_markers.toml | mass.f90 | 93-188 (36975d7) | 93 | 99-194 | **99-194** | **99** | +6 |
| rho_rmw_mass | s5_markers.toml | mass.f90 | 201-209 (36975d7) | 201 | 207-215 | **207-215** | **207** | +6 |
| mass_wall_rmw | s5_markers.toml | mass.f90 | 224-320 (36975d7) | 224 | 230-326 | **230-326** | **230** | +6 |
| flux_mw_fix | s5_markers.toml | mass.f90 | 328-350 (36975d7) | 328 | 334-356 | **334-356** | **334** | +6 |
| wall_uvw_interp | s5_markers.toml | mass.f90 | 424-436 (36975d7) | 424 | 430-442 | **430-442** | **430** | +6 |
| zzs_pred | s5_markers.toml | mass.f90 | 442-455 (36975d7) | 442 | 448-461 | **448-461** | **448** | +6 |
| adv_flux_store | s5_markers.toml | mass.f90 | 462-472 (36975d7) | 462 | 468-478 | **468-478** | **468** | +6 |
| rhos_sum | s5_markers.toml | mass.f90 | 505-512 (36975d7) | 505 | 511-518 | **511-518** | **511** | +6 |
| rsum_pred | s5_markers.toml | mass.f90 | 527-534 (36975d7) | 527 | 533-540 | **533-540** | **533** | +6 |
| rsum_ideal_pred | s5_markers.toml | mass.f90 | 556-564 (36975d7) | 556 | 562-570 | **562-570** | **562** | +6 |
| tmp_pbar_pred | s5_markers.toml | mass.f90 | 570-577 (36975d7) | 570 | 576-583 | **576-583** | **576** | +6 |
| wall_uvw_interp_corr | s5_markers.toml | mass.f90 | 595-608 (36975d7) | 595 | 601-614 | **601-614** | **601** | +6 |
| zz_corr | s5_markers.toml | mass.f90 | 619-632 (36975d7) | 619 | 625-638 | **625-638** | **625** | +6 |
| adv_flux_avg | s5_markers.toml | mass.f90 | 650-660 (36975d7) | 650 | 656-666 | **656-666** | **656** | +6 |
| rho_sum | s5_markers.toml | mass.f90 | 687-694 (36975d7) | 687 | 693-700 | **693-700** | **693** | +6 |
| zz_div_rho | s5_markers.toml | mass.f90 | 709-716 (36975d7) | 709 | 715-722 | **715-722** | **715** | +6 |
| rsum_ideal_corr | s5_markers.toml | mass.f90 | 738-746 (36975d7) | 738 | 744-752 | **744-752** | **744** | +6 |
| tmp_pbar_corr | s5_markers.toml | mass.f90 | 752-759 (36975d7) | 752 | 758-765 | **758-765** | **758** | +6 |
| delta_rho_zz_zero | s5_markers.toml | mass.f90 | 871-871 (36975d7) | 871 | 877-877 | **877-877** | **877** | +6 |
| rho_zz_clip_assign | s5_markers.toml | mass.f90 | 931-937 (36975d7) | 931 | 937-943 | **937-943** | **937** | +6 |
| rho_zz_final_check | s5_markers.toml | mass.f90 | 947-961 (36975d7) | 947 | 953-967 | **953-967** | **953** | +6 |
| zeta_clip | s5_markers.toml | mass.f90 | 980-987 (36975d7) | 980 | 986-993 | **986-993** | **986** | +6 |
| wall_b2_work1 | s5_markers.toml | part.f90 | 3654-3659 (36975d7) | 3654 | 3654-3659 | **3654-3659** | **3654** | no |
| wall_hs_bt | s5_markers.toml | pres.f90 | 4351-4356 (36975d7) | 4351 | 4351-4356 | **4351-4356** | **4351** | no |
| rad_wall_qin_zero | rad_kernels.toml | radi.f90 | 3887-3893 (bee11f0329) | 3887 | 3887-3893 | **3887-3893** | **3887** | no |
| rad_gray_kfst4 | rad_kernels.toml | radi.f90 | 4099-4117 (bee11f0329) | 4099 | 4099-4117 | **4099-4117** | **4099** | no |
| rad_thin_kfst4 | rad_kernels.toml | radi.f90 | 4131-4141 (bee11f0329) | 4131 | 4131-4141 | **4131-4141** | **4131** | no |
| rad_corr_kfst4 | rad_kernels.toml | radi.f90 | 4151-4161 (bee11f0329) | 4151 | 4151-4161 | **4151-4161** | **4151** | no |
| rad_extcoe | rad_kernels.toml | radi.f90 | 4206-4206 (bee11f0329) | 4206 | 4206-4206 | **4206-4206** | **4206** | no |
| rad_uiiold_wb | rad_kernels.toml | radi.f90 | 4213-4213 (bee11f0329) | 4213 | 4213-4213 | **4213-4213** | **4213** | no |
| rad_uiiold_gray | rad_kernels.toml | radi.f90 | 4215-4215 (bee11f0329) | 4215 | 4215-4215 | **4215-4215** | **4215** | no |
| rad_rte_source | rad_kernels.toml | radi.f90 | 4221-4222 (bee11f0329) | 4221 | 4221-4222 | **4221-4222** | **4221** | no |
| rad_uiid_wb | rad_kernels.toml | radi.f90 | 4818-4818 (bee11f0329) | 4818 | 4818-4818 | **4818-4818** | **4818** | no |
| rad_uiid_gray | rad_kernels.toml | radi.f90 | 4820-4820 (bee11f0329) | 4820 | 4820-4820 | **4820-4820** | **4820** | no |
| rad_qr_wb | rad_kernels.toml | radi.f90 | 4949-4949 (bee11f0329) | 4949 | 4949-4949 | **4949-4949** | **4949** | no |
| rad_emis_wb | rad_kernels.toml | radi.f90 | 4951-4951 (bee11f0329) | 4951 | 4951-4951 | **4951-4951** | **4951** | no |
| rad_abs_wb | rad_kernels.toml | radi.f90 | 4952-4952 (bee11f0329) | 4952 | 4952-4952 | **4952-4952** | **4952** | no |
| rad_qrw_wb | rad_kernels.toml | radi.f90 | 4954-4954 (bee11f0329) | 4954 | 4954-4954 | **4954-4954** | **4954** | no |
| wall_open_qin | s5_markers.toml | radi.f90 | 4961-4970 (36975d7) | 4961 | 4965-4974 | **4965-4974** | **4965** | +4 |
| rad_uii_sum | rad_kernels.toml | radi.f90 | 4963-4963 (bee11f0329) | 4963 | 4963-4963 | **4963-4963** | **4963** | no |
| rad_qr_gray | rad_kernels.toml | radi.f90 | 4982-4982 (bee11f0329) | 4982 | 4982-4982 | **4982-4982** | **4982** | no |
| rad_qrw_gray | rad_kernels.toml | radi.f90 | 4983-4983 (bee11f0329) | 4983 | 4983-4983 | **4983-4983** | **4983** | no |
| settle_zero | s5_markers.toml | soot.f90 | 56-58 (36975d7) | 56 | 56-58 | **56-58** | **56** | no |
| settle_vel | s5_markers.toml | soot.f90 | 60-132 (36975d7) | 60 | 60-132 | **60-132** | **60** | no |
| settle_wall_zero | s5_markers.toml | soot.f90 | 136-160 (36975d7) | 136 | 136-160 | **136-160** | **136** | no |
| settle_div | s5_markers.toml | soot.f90 | 162-174 (36975d7) | 162 | 162-174 | **162-174** | **162** | no |
| mu_dns | s5_markers.toml | velo.f90 | 80-88 (36975d7) | 80 | 89-97 | **89-97** | **89** | +9 |
| up_deardorff | s5_markers.toml | velo.f90 | 137-145 (36975d7) | 137 | 146-154 | **146-154** | **146** | +9 |
| wall_up_ghost | s5_markers.toml | velo.f90 | 151-164 (36975d7) | 151 | 160-173 | **160-173** | **160** | +9 |
| kres | s5_markers.toml | velo.f90 | 285-294 (36975d7) | 285 | 294-303 | **294-303** | **294** | +9 |
| wall_visc_mirror | s5_markers.toml | velo.f90 | 355-363 (36975d7) | 355 | 364-372 | **364-372** | **364** | +9 |
| wall_strain_rate | s5_markers.toml | velo.f90 | 440-483 (36975d7) | 440 | 449-492 | **449-492** | **449** | +9 |
| vflux_vort_tau | s5_markers.toml | velo.f90 | 619-639 (36975d7) | 619 | 631-651 | **631-651** | **631** | +12 |
| vflux_fvx | s5_markers.toml | velo.f90 | 662-714 (36975d7) | 662 | 674-726 | **674-726** | **674** | +12 |
| vflux_fvy | s5_markers.toml | velo.f90 | 720-772 (36975d7) | 720 | 732-784 | **732-784** | **732** | +12 |
| vflux_fvz | s5_markers.toml | velo.f90 | 778-830 (36975d7) | 778 | 790-842 | **790-842** | **790** | +12 |
| wall_coriolis_ghost | s5_markers.toml | velo.f90 | 981-987 (36975d7) | 981 | 993-999 | **993-999** | **993** | +12 |
| pvf_fvx | s5_markers.toml | velo.f90 | 1105-1127 (36975d7) | 1110 | 1117-1139 | **1117-1139** | **1122** | +12 |
| pvf_fvy | s5_markers.toml | velo.f90 | 1138-1161 (36975d7) | 1143 | 1150-1173 | **1150-1173** | **1155** | +12 |
| pvf_fvz | s5_markers.toml | velo.f90 | 1172-1194 (36975d7) | 1177 | 1184-1206 | **1184-1206** | **1189** | +12 |
| vpred_us | s5_markers.toml | velo.f90 | 1603-1609 (36975d7) | 1603 | 1622-1628 | **1622-1628** | **1622** | +19 |
| vpred_vs | s5_markers.toml | velo.f90 | 1613-1619 (36975d7) | 1613 | 1632-1638 | **1632-1638** | **1632** | +19 |
| vpred_ws | s5_markers.toml | velo.f90 | 1623-1629 (36975d7) | 1623 | 1642-1648 | **1642-1648** | **1642** | +19 |
| vcorr_u | s5_markers.toml | velo.f90 | 1724-1730 (36975d7) | 1724 | 1743-1749 | **1743-1749** | **1743** | +19 |
| vcorr_v | s5_markers.toml | velo.f90 | 1734-1740 (36975d7) | 1734 | 1753-1759 | **1753-1759** | **1753** | +19 |
| vcorr_w | s5_markers.toml | velo.f90 | 1744-1750 (36975d7) | 1744 | 1763-1769 | **1763-1769** | **1763** | +19 |
| cfl_max | s5_markers.toml | velo.f90 | 3059-3080 (36975d7) | 3059 | 3093-3114 | **3093-3114** | **3093** | +34 |
| cfl_wall_max | s5_markers.toml | velo.f90 | 3093-3108 (36975d7) | 3093 | 3127-3142 | **3127-3142** | **3127** | +34 |
| vn_max | s5_markers.toml | velo.f90 | 3121-3134 (36975d7) | 3121 | 3155-3168 | **3155-3168** | **3155** | +34 |
| baro_p_rrho | s5_markers.toml | velo.f90 | 3254-3261 (36975d7) | 3254 | 3288-3295 | **3288-3295** | **3288** | +34 |
| baro_fvx | s5_markers.toml | velo.f90 | 3267-3275 (36975d7) | 3267 | 3301-3309 | **3301-3309** | **3301** | +34 |
| baro_fvy | s5_markers.toml | velo.f90 | 3282-3290 (36975d7) | 3282 | 3316-3324 | **3316-3324** | **3316** | +34 |
| baro_fvz | s5_markers.toml | velo.f90 | 3297-3305 (36975d7) | 3297 | 3331-3339 | **3331-3339** | **3331** | +34 |
| wall_un_store | s5_markers.toml | velo.f90 | 3347-3368 (36975d7) | 3347 | 3381-3402 | **3381-3402** | **3381** | +34 |
| wall_us_pred | s5_markers.toml | velo.f90 | 3378-3408 (36975d7) | 3378 | 3412-3442 | **3412-3442** | **3412** | +34 |
| wall_uvw_nograd | s5_markers.toml | velo.f90 | 3414-3450 (36975d7) | 3414 | 3448-3484 | **3448-3484** | **3448** | +34 |

Wall-list variants (no `lines`, follow their base kernel): `wall_uvw_interp_wl` (wall_uvw_interp), `wall_uvw_interp_corr_wl` (wall_uvw_interp_corr), `wall_up_ghost_wl` (wall_up_ghost), `wall_kp_ghost_wl` (wall_kp_ghost), `wall_rho_d_dzdn_wl` (wall_rho_d_dzdn), `wall_un_store_wl` (wall_un_store), `wall_coriolis_ghost_wl` (wall_coriolis_ghost), `wall_hs_bt_wl` (wall_hs_bt), `wall_us_pred_wl` (wall_us_pred), `wall_strain_rate_wl` (wall_strain_rate), `wall_b2_work1_wl` (wall_b2_work1), `wall_uvw_nograd_wl` (wall_uvw_nograd), `wall_visc_mirror_wl` (wall_visc_mirror), `wall_open_qin_wl` (wall_open_qin).

## 6. Pre-existing drift of the work list (not caused by this merge)

169 of the 231 work-list loops already have different line numbers in the FDS-AMReX tree than at `36975d7`; the shifts come from the branch patches and the `bee11f0329` merge. Three loops also differ in text between `36975d7` and `9dc218d55b`: L1364, L1367 and L1399 (`velo.f90`), each because of a branch-patch guard `IF (EXTERNAL_GHOSTS_FILLED) CYCLE` under `WITH_AMREX` (for example `velo.f90:1392-1394` in L1364 at `8ad685f435`). All three are claimed by the GPU Mesh Data Loops Engineer. The ranges at `8ad685f435` for every loop handed out next round are in `loop-handout-round-next.md`.

## 7. Claims register changes made with this report (`amrex/loop_claims.csv`)

- L0630: new row, no owner, status `hold-upstream`, note "FireX 1fb5c300 changed body; see upstream-1fb5c300-loop-impact.md".
- MP5 (owner ruling: MP5 is not implemented in the generator):
  - L0401 and L0398 are claimed by the Legacy Mapper. The note "MP5 not implemented in the generator (owner ruling)" was appended; status unchanged.
  - L0657 (`GET_SCALAR_FACE_VALUE` MP5 nest, `func.f90:1439-1453`) is unclaimed, but the committed whole-field kernel `gsfv_mp5` on s5-gen already translates it. A row with the ruling note was added and the status left as `translated`, not `closed-owner-ruling`, until the owner decides whether `gsfv_mp5` is withdrawn.
- No body-changed loop is claimed, so no active claim was touched by the hold.

## Mapper tip re-check (2026-10-06)

Re-checked against FDS-AMReX tip `315ced11efc5` (merge `8ad685f435` still the FireX land). Survey pins and the §1 tip ranges for L0630/L0677/L0817/L0822/L0823 still hold. Additional tip notes (no inventory rewrite):

- **L1356** WALE nest is `velo.f90:265-288` at tip/8ad685 (survey `256-279`); work-list survey pin unchanged.
- Species/radiation register loops (L0876/L0877/L1272/L0866/L0875/L0397/L1242 and R1 walls) match the tip DO starts already cited in claims (`mass` +6/+6, `divg` +15, `radi` +1/+4).
- UP-0011 Intel Rel+Dbg PASS; GNU Debug CLOSED; tip `main.f90` has further MAIN_LOOP drift (RHS call ~1740) and already carries the 0011 MODE==3 dump loop. Details: `upstream-patches/line-check-firex-1fb5c300.md` § Mapper tip re-check.
