# WP9b part B: host pieces and the host/device boundary of the divergence chain

Scope: `DIVERGENCE_PART_1` in the FDS-AMReX tree (`Source/divg.f90`, hook call `FDS_HOOK_DIF_FLUX` at divg.f90:279), split at the DIF hook into part A (before) and part B (after). This note fixes, per statement group, what runs as a device kernel and what stays on the host, which arrays cross, who owns each host piece, and how the split is checked. It extends section 3 of `wp3-generator-coverage.md` (the first draft of the split, item 9b-1) and WP9 / WP9b of `stage1-gpu-spike-plan.md`.

Labels: **read** = from source text or documents; **run** = measured in the spike (cc 8.9 test-machine GPU, host); **ESTIMATE** = arithmetic or judgement, not measured. No GPU was run for this note. Nothing in the driver tree, the generator tree or the AMReX clone was written.

Sources read: `Source/divg.f90` (lines 28-800), `Source/driver/fds_flux_hooks.f90`, `Source/driver/TimeLoop.cpp` (`flux_readout_dif`, `flux_apply_dif`, `FluxStages`), the generator's `s5_markers.toml`, the plan sections 4, 7.8, 7.11, 7.12, 7.13a, 10 and 11, `flux-hooks-kernel-review.md`, `wp3-generator-coverage.md`, the Legacy Mapper's table of call order (`work-review-results.md` item 3; that item audits the writers of `FVX/FVY/FVZ`, it holds no statement table of `DIVERGENCE_PART_1`, so only its call-order table was used) and its NIC analysis (item 2).

## 0. Findings that change or sharpen earlier statements

1. **The whole species block is skipped for one species.** Lines 124-438 (`SPECIES_GT_1_IF: IF (N_TOTAL_SCALARS>1)`) contain part A's flux loop, the species-sum fix, the hook call and the species loop of part B. `tg128` and `csmag_32` have one species, so on them the hook is never reached and no flux array is filled. Only `shunn3_32` and `dec2_obst` (two species) execute the split region. Any two-level test of WP9/WP9b needs at least two species. (read)
2. **Line numbers.** The generator's `s5_markers.toml` quotes the unpatched `divg.f90`. The FDS-AMReX file differs from it by a 6-line header (patch 0008) plus 3 `USE` lines (patch 0009) at the top, so every marker line is **+9** in this file up to the hook, and **+15** after the hook call (the call block adds 6 more). Example: marker `rho_d_maxloc_fix` 245-258 is 254-267 here; `dp_div_heat` 390-402 is 405-417. The statement in the review and in item 9b-10 that patch 0009 "shifts lines after 267 by 9" is therefore not what the files show (checked by a diff of the two `divg.f90`; the shift is 9 everywhere before the hook and 15 after it). Line numbers below are the FDS-AMReX file unless marked "marker". (read)
3. **`WALL_LOOP_2` (333-401) sits in part B, after the hook, inside the species loop**, and overwrites `RHO_D_DZDX/Y/Z` and `H_RHO_D_DZDX/Y/Z` at the faces of solid-type walls before the divergence loops use them. It has no kernel and contains a second `MAXLOC` (line 341, over `B1%ZZ_F`). The earlier note listed it only as "host". It is the largest open piece of part B. (read)
4. **The DIF hook can only touch faces that part B does not overwrite**, provided the override list holds interface faces only: `WALL_LOOP_2` skips `NULL`, `OPEN` and `INTERPOLATED` walls (lines 335-337) and the hook writes only listed interior-range faces. A list face that coincides with a solid-wall face would be overwritten silently; no check for this exists (open question 4).
5. **`PREDICT_NORMAL_VELOCITY` runs at the start of the chain** (lines 101-105, all walls) and writes `B1%U_NORMAL_S` / `B1%U_NORMAL`, which `WALL_LOOP_2` (373-377) and the zone-volume flux (780-781) read. A wall-table dump or upload taken before line 105 holds stale normal velocities. This is the same kind of input-timing trap as the corrector-time rule of plan 7.9c. (read)
6. **Real FDS truth for the part A outputs already exists in the driver**: the DIF read-out (hook mode 1) copies exactly the species-sum-fixed `RHO_D_DZDX/Y/Z` (fds_flux_hooks.f90:136-143). Using it in the scratch driver gives truth for `rho_d_dzd` + `wall_rho_d_dzdn` + `rho_d_maxloc_fix` together without a new dump patch. What is still not dumped is `RHO_D` and the `D_Z` table (inputs). (read)

## 1. (a) Part A and part B per statement group

Part A = everything before the hook call; hook = lines 276-280; part B = everything after. "Marker" = line in the generator's file. Kernel names are those of `s5_markers.toml`. "K2" = generated Fortran OpenMP-offload kernel (plan 7.12: blocking launch). "K1" = C++ `ParallelFor` allowed for data-movement only (`adr001-data-movement-kernel-entry.md`).

### Part A (before the hook)

| # | Lines | Statement group | Kernel today | Side in the split | Notes |
|---|---|---|---|---|---|
| A0 | 61-97 | `SOLID_PHASE_ONLY` return, `POINT_TO_MESH`, stage pointers (`DP=>DS` or `D`, `ZZP`, `RHOP`, `UU`...), `R_PBAR = 1/PBAR_P`, `DP = 0`, `MERGE_PRESSURE_ZONES` | none | host control; `DP = 0` becomes a device fill (first device statement of the chain); `R_PBAR` is a 1-D table (uploaded per stage or per step); `MERGE_PRESSURE_ZONES` host, geometry event only (kernel `wall_connect_zones` exists for its wall loop) | the predictor/corrector pointer choice is host logic that selects which device array is passed |
| A1 | 101-105 | `WALL_LOOP3`: `PREDICT_NORMAL_VELOCITY(IW,T,DT)` over external + internal walls | none | **host** (wall routine with ramps); writes `B1%U_NORMAL_S/U_NORMAL` | must finish before B's `WALL_LOOP_2` and the zone-flux loop; wall-state upload (WP10) fires after it (finding 5) |
| A2 | 107-110 | `CC_IBM` branch | none | out of scope (cut cells off) | assert `.NOT.CC_IBM` on the device chain |
| A3 | 112-135 | `D_Z_MAX = 0` (if `CHECK_VN`), `DEL_RHO_D_DEL_Z = 0`, `RHO_D = MAX(0,MU)*RSC_T` or `RHO_D_TURB` (whole-array forms) | none (whole-array forms are the DP1 generator blocker) | device fills; the `MU`-based forms need a kernel or K1 fill; only needed for LES/VLES (`SIM_MODE`) | `csmag_32`-type LES cases reach lines 131 and 133 |
| A4 | 137-150 | per species N: `RHO_D = 0`, `D_Z_N = D_Z(:,N)`, `RHO_D = RHOP*interp(D_Z_N, TMP)` | `rho_d_interp` (marker 133-140) | A, device | the column alias is handled by the generator (`policy.column_alias`); table `D_Z(0:I_MAX_TEMP, NS)`, `I_MAX_TEMP = 5000` (cons.f90:481), uploaded once per run |
| A5 | 152-159 | LES turbulent add `RHO_D = RHO_D + RHO_D_TURB*...` | none (DP1) | A, device (needs the whole-array rewrite) | |
| A6 | 163 | `PERIODIC_TEST==7` manufactured `RHO_D` | none | host-only test path; assert off on the device | |
| A7 | 167-175 | `D_Z_MAX = MAX(D_Z_MAX, RHO_D/(RHOP+eps))` (if `CHECK_VN`) | `d_z_max` (marker 159-165) | A, device | `D_Z_MAX` also gets a second update at 512-522 (kernel not found in the marker list) and updates in `turb.f90` (lines 1348-1437) |
| A8 | 179-192 | `RHO_D_DZDX/Y/Z(I,J,K,N) = .5*(RHO_D(+1)+RHO_D)*dZ/dx` over faces 0..IBAR, 0..JBAR, 0..KBAR | `rho_d_dzd` (marker 171-182) | A, device | |
| A9 | 196 | `TENSOR_DIFFUSIVITY_MODEL(NM,N)` | none | host-only model; assert off on the device chain | |
| A10 | 200-242 | `WALL_LOOP`: zero (default walls), or store/overwrite the flux at `OPEN` / `INTERPOLATED` walls (`NIC>1` branch reads `B1%RHO_D_DZDN_F`, else stores it) | `wall_rho_d_dzdn` (marker 192-232; `wall_rho_d_dzdn_wl` is the list variant for W1) | A, device with the wall-list indirection (W1/`WLIST`) | the `NIC>1` branch (line 219) is not exercised by any case in the harness (plan 7.10a); the thin-wall `CYCLE` (line 205) must be kept |
| A11 | 248-270 | species-sum fix: per face `N = MAXLOC(ZZP(..)+ZZP(neighbour))`, `R(N) = -(SUM(R) - R(N))` (DNS, LES or tensor diffusivity only) | `rho_d_maxloc_fix` (marker 245-258) | A, device | D-051 sign-off pending (open question 1); NaN limit (9b-3); sign-of-zero rule applies (plan 7.13a) |
| A12 | 274 | `SET_EXIMDIFFLX_3D` (`CC_IBM`) | none | out of scope | |

### The hook (lines 276-280, inside `WITH_AMREX`)

`CALL FDS_HOOK_DIF_FLUX(NM, LBOUND(RHO_D_DZDX,4), RHO_D_DZDX, RHO_D_DZDY, RHO_D_DZDZ)`. Mode 0: return. Mode 1: copy the faces `I=0..IBAR, J=1..JBAR, K=1..KBAR` (X; the analogous ranges for Y and Z), all scalars, into the registered arrays. Mode 2: write the listed `(dir, i, j, k)` entries (`VAL(1:NS)`) into the flux arrays. In the split the hook becomes: **device gather of the listed (or interface) faces, copy to the host, [host work of the regrid/transport role], copy of the list values to the device, device scatter**. Mode and list state stay on the host (review finding 1).

### Part B (after the hook)

| # | Lines | Statement group | Kernel today | Side in the split | Notes |
|---|---|---|---|---|---|
| B1 | 282-298 | `STORE_SPECIES_FLUX` output copies `DIF_FX = 0.5*(DIF_FXS - RHO_D_DZDX)` (predictor), `DIF_FXS = -RHO_D_DZDX` (corrector) | none | **host or unsupported on the first device chain**; only active for the `DIFFUSIVE/TOTAL MASS FLUX` output quantities (read.f90:17010-17014), default off | arithmetic, so a K2 kernel rather than a K1 data-movement kernel if ever moved; `DIF_F*(0:IBP1,..,1:NS)` carry state across stages |
| B2 | 302-327 | per species N: `H_RHO_D_DZDX/Y/Z = h_s(N,T_face)*RHO_D_DZD*` | `h_rho_d_dzd` (marker 298-312) | B, device | `WORK5-7` are reused for each N |
| B3 | 333-401 | `WALL_LOOP_2`: wall `RHO_D_DZDN`, its species-sum fix (`MAXLOC(B1%ZZ_F)`), `B1%RHO_D_DZDN_F(N)`, `DIF_F*` at the wall, wall-face overwrite of `RHO_D_DZD*` and `H_RHO_D_DZD*` | none (DP1: array sections, wall species tables, `MAXLOC`, `STORE_SPECIES_FLUX`, enthalpy call) | **host path** with a wall-list gather/scatter (section 3, P4), or a new kernel (generator, P5) | writes `B1%RHO_D_DZDN_F`, which the next `WALL_BC` reads (wall.f90:1193) |
| B4 | 404-418 | `DP += (R(I)*H_RHO_D_DZDX(I) - R(I-1)*H_RHO_D_DZDX(I-1))*RDX*RRN + ...` | `dp_div_heat` (marker 390-402) | B, device | bitwise tests: NO-TEST on the host tier (ported review) |
| B5 | 422-432 | `DEL_RHO_D_DEL_Z(I,J,K,N) = (R(I)*RHO_D_DZDX(I) - R(I-1)*RHO_D_DZDX(I-1))*RDX*RRN + ...` | `del_rho_d_del_z` (marker 408-416) | B, device | |
| B6 | 442-463 | specific heat `CP`, `R_H_G` | `cp_rhg` (marker 435-443) | B, device | only if `.NOT.CONSTANT_SPECIFIC_HEAT_RATIO` |
| B7 | 474-508 | thermal conductivity `KP`, LES add, ghost copy at external walls | `conductivity` (marker 462-471), `wall_kp_ghost` (marker 482-486) | B, device (wall list) | |
| B8 | 512-522 | `D_Z_MAX = MAX(.., KP/(CP*RHOP))` | not found in the marker list | B, device, kernel missing | |
| B9 | 527-543 | `KDTDX/Y/Z` | `kdtd` (marker 512-523); `TENSOR_DIFFUSIVITY_MODEL(NM)` host | B, device | |
| B10 | 547-569 | `CORRECTION_LOOP`: `K_G`, `DP` at walls (`Q_CON_F`, `Q_LEAK`), zero of `KDTD*` | `wall_corr_kdtd` (marker 532-554, gather `DP`) | B, device (wall list) | writes `B1%K_G` (host reader) |
| B11 | 573-598 | `DP += div(k grad T) + Q + QR` | `dp_kdtd` (marker 561-570) | B, device | |
| B12 | 602-616, 651-677 | `ENTHALPY_ADVECTION_NEW`, `SPECIES_ADVECTION_PART_1_NEW/2`, `DP` updates | `cell_adv_hs`, `gsfv_*`, `rho_z_p_divg`, `cell_adv_zz`, `wall_spec_adv2`, `dp_species` (marker 646-655) and others | B, device where kernels exist | not examined statement by statement here; outside the hook split (`wp3-generator-coverage.md` section 2) |
| B13 | 620-647, 681-704 | `RTRM`, `DP = RTRM*DP`; `D_SOURCE`; stratification | none found in the marker list | B, device, kernels missing | |
| B14 | 708-732 | `PERIODIC_TEST==7` manufactured source | none | host test path; assert off | |
| B15 | 740-796 | zone sums `DSUM`, `PSUM`, `USUM` (only if `N_ZONE>0`) | none (D-053: serial add in FDS order) | per-cell terms on the device, serial addition in FDS order, wall-flux sum on the host or ordered | `zone_save/zone_restore` of the re-run design (TimeLoop.cpp:577-581) disappear with the split, unless `N_ZONE>0` forces a snapshot |

Order constraint inside part B, per species N: `h_rho_d_dzd(N)` -> wall fix(N) (B3) -> `dp_div_heat(N)` -> `del_rho_d_del_z(N)`. The wall fix must follow the cell loop B2 and precede B4 and B5 (it overwrites the arrays they read). Per species the part B sequence is therefore three kernels plus the wall fix; part A is four kernels per species plus one species-sum kernel (ESTIMATE of the launch count at two species: A 9, B 6 plus 2 for the wall fix, about 15-17 launches; at the measured 7.7 us per tiny blocking K2 launch (run, plan 7.12) this is about 0.12-0.13 ms of launch floor per stage; ESTIMATE, not measured on this chain).

## 2. (b) Arrays crossing the host/device boundary

Sizes: double precision. `ZZ/ZZS/TMP/RHOS` have two ghost layers `(-1:N+2)`; `MU`, `D`, `DS`, `SWORK*`, `WORK*`, `DEL_RHO_D_DEL_Z` have one `(0:N+1)` (init.f90:524-528, 560-562, 590, 694-703). One component: 128^3 box with two ghost layers 18.4 MB, with one ghost layer 17.6 MB; 16^3 box 64.0 kB / 46.7 kB. Hook face range per direction and scalar (what the read-out copies): 128^3 129x128x128 = 16.9 MB; 16^3 17x16x16 = 34.8 kB. Walls of a box with all six faces external: 6 N^2 = 98,304 at 128^3 and 1,536 at 16^3 (ESTIMATE of the count; real counts depend on the case). **NS = number of scalars (`N_TOTAL_SCALARS`). `tg128` has NS = 1 and does not execute the split region (finding 1): its actual DIF crossing is zero; the "if it ran" columns use the same formulas for NS scalars so that they can be scaled to a 2-species 128^3 case.** All copy times are ESTIMATES from the measured rates of plan 7.8 (12 GB/s large copies, 5 GB/s at 64 KB, 10-25 us latency per copy).

### Device-resident (no crossing in steady state)

| Array | Shape (FDS name) | Direction | When | 128^3 (tg128) | 16^3 box |
|---|---|---|---|---|---|
| `ZZ`, `ZZS` | `(-1:IBP1+1,.., NS)` | produced by density kernels | per stage | 18.4 MB x NS, no copy | 64 kB x NS, no copy |
| `TMP`, `RHOS`/`RHO`, `MU`, `D`/`DS` | cell arrays | produced earlier in the stage | per stage | 18.4 or 17.6 MB each, no copy | 64 or 46.7 kB each, no copy |
| `RHO_D`, `WORK4-7`, `WORK9`, `KP`, `KDTD*`, `CP`, `R_H_G`, `RTRM` | `(0:IBP1,..)` scratch | scratch, written and read on the device | per stage | no copy | no copy |
| `RHO_D_DZDX/Y/Z` (`SWORK1-3`) | `(0:IBP1,..,N_LOWER:NS)`, `N_LOWER` is 0 or 1 (init.f90:606-608) | device, crossing only through the hook below | per stage | 17.6 MB x components x 3 | 46.7 kB x components x 3 |
| `DEL_RHO_D_DEL_Z` | `(0:IBP1,..,NS)` | device; also read by the next predictor `DENSITY` (`DEL_RHO_D_DEL_Z__0`, mass.f90:418) | per stage | 17.6 MB x NS | 46.7 kB x NS |
| `D_Z_MAX` | `(0:IBP1,..)` | device; consumed by `vn_max`; `turb.f90` also updates it (host side unknown, question 6) | per stage | 17.6 MB | 46.7 kB |
| `DP` | `(0:IBP1,..)` | device, output of the chain (input of `DIVERGENCE_PART_2` and the pressure RHS) | per stage | 17.6 MB | 46.7 kB |
| `DIF_F*`, `DIF_F*S` | `(0:IBP1,..,NS)` | device only if B1 is supported; else absent | per stage | 17.6 MB x NS x 3 x 2 | 46.7 kB x NS x 3 x 2 |

### Uploaded tables (host to device)

| Array | Shape | When | Size |
|---|---|---|---|
| `D_Z`, `H_SENS_Z`, `CP_Z`, `K_RSQMW_Z`, `RSQ_MW_Z`, `MW` (species and property tables) | `(0:5000, NS)` or `(NS)` | once per run (WP8 proposal: table ownership by the driver shim) | 40 kB x NS per `(0:5000)` table, the same at both box sizes |
| metrics `RDX, RDY, RDZ, RDXN, RDYN, RDZN, R, RRN` | 1-D | per regrid | under 10 kB at 128^3, under 1 kB at 16^3 |
| `R_PBAR` | `(KBP1, N_ZONE)` | per stage (changes with `PBAR`) | KBP1 x zones x 8 B |
| `SOLID` mask, `PRESSURE_ZONE` | cell integer arrays | at geometry events, per regrid | about 8.8 MB / 23 kB each as 32-bit integers (ESTIMATE of the type) |

### Crossings of part A, the hook and part B

| Array | Shape | Direction | When | Size at 128^3 (tg128, if it ran with NS scalars) | Size at a 16^3 box | Note |
|---|---|---|---|---|---|---|
| hook read-out, **whole-box form** (first build) | 3 directions x hook range x NS | device to host | per stage and level, only with a registered box | 50.7 MB x NS, 4.2 ms at 12 GB/s (the review quotes 52 MB: that is the full `(0:IBP1)^3` array, 52.7 MB) | 104 kB x NS, about 9 us at the large-copy rate, 25-40 us with latency (ESTIMATE) | works with no new kernel; the cost of the first version |
| hook read-out, **gather form** | listed interface faces x NS | device to host | per stage and level | faces x 8 x NS: 6 sides of 128^2 faces = 98,304 faces, 0.79 MB x NS; the 64^3-fine-patch example of the review (6 x 32^2 = 6,144 faces) 49 kB x NS | one side 256 faces = 2.0 kB x NS; all six sides 1,536 faces = 12.3 kB x NS | gather is a K1 `ParallelFor` (data movement, ADR-001 entry) on the development machine stream, so a stream sync is needed before the copy (plan 7.12) |
| face list (dir, i, j, k) | 16 B per face | host to device | when the list changes (regrid), not per stage | 1.6 MB for 98,304 faces | 25 kB for 1,536 faces | uploaded once per regrid |
| override values | NS doubles per face | host to device | per stage and level (values change every stage) | faces x 8 x NS (same as the gather) | same as the gather | scatter is K1 data movement; the device scatter replaces `RX(I,J,K,1:NS)` |
| `B1%U_NORMAL_S`, `B1%U_NORMAL` (written by A1) | per wall | host to device | per stage, after A1 | 8 B x walls: 0.79 MB | 12 kB | part of the wall-state refresh (WP10: 0.32-0.39 ms at 98,304 walls (run)); B3 and B15 read them |
| wall species inputs `B1%ZZ_F`, `B1%RHO_D_F`, `B1%RDN`, `B1%TMP_F` | (2 NS + 2) doubles per wall | host to device if B3 is a kernel; unused by the host path | per stage (set by `WALL_BC`) | about 8 B x (2 NS + 2) x walls | same formula | |
| `B1%RHO_D_DZDN_F(1:NS)` (written by A10 and B3, `NIC>1` read by A10) | NS per wall | device to host if A10 and B3 are kernels; host if the host path is used for B3 | per stage | 8 B x NS x walls = 0.79 MB x NS | 12 kB x NS | read by the next `WALL_BC` (wall.f90:1193, host), by output (dump.f90:10562-10631) and by the MPI pack (func.f90:5045); never written for NS = 1 (stays 0, func.f90:4946) |
| B3 host path, gather | `ZZP(IIG,JJG,KKG,1:NS)` and `TMP(IIG,JJG,KKG)` per wall | device to host | once per stage at the start of part B (inputs are not modified inside the chain) | 8 B x (NS+1) x walls: 2.4 MB at NS = 2 | 37 kB at NS = 2 | ESTIMATE: 0.2 ms at 12 GB/s plus the host loop |
| B3 host path, scatter | `RHO_D_DZD*` and `H_RHO_D_DZD*` value per wall and species + face index | host to device | once per stage, applied between B2 and B4 of each species | 8 B x 2 x NS x walls + 4-12 B index: 3.1 MB at NS = 2 | 49 kB at NS = 2 | K1 scatter (pure overwrite, bitwise); a face is written by exactly one wall except at thin obstructions, where FDS skips one side (line 349) |
| `DIF_F*` / `DIF_F*S` | `(0:IBP1,..,NS)` x 3 | device to host if `STORE_SPECIES_FLUX` and the host outputs read them | per output interval | 50.7 MB x NS if copied as whole arrays | 104 kB x NS | default off (question 3) |
| `DEL_RHO_D_DEL_Z` | cell, NS | device to host | restart write only (dump.f90:3929) and `soot.f90:171` if deposition | 17.6 MB x NS | 46.7 kB x NS | rare |
| zone sums `DSUM`, `PSUM`, `USUM` | `N_ZONE` doubles | device to host | per stage, only if `N_ZONE>0` | 8 B x zones | same | the serial addition order is the constraint (D-053), not the size |
| `D_Z_MAX` | cell | device to host for the output quantity (dump.f90:8978) and for the stability scalar (a reduction, so 8 B) | per output interval / per stage | 8 B (reduction) | 8 B | |

## 3. (c) Part B host pieces, owners and dependencies

Roles only. "Integration Lead" = the kernel-call and measurement owner of this spike; "Generator Engineer" = the GPU Generator Engineer; "Data Layout Implementer" = the owner of the driver tree and the hook module.

| Piece | What | Owner | Depends on | Acceptance |
|---|---|---|---|---|
| P1 Launch sequence A - hook - B | per stage and level: part A kernels, gather, [host], scatter, part B kernels; `readout_dif(level)` becomes "part A + gather", `run_divergence_part1(level)` becomes "scatter + part B"; the double run (`flux_apply_dif`) and `zone_save/zone_restore` go away | Integration Lead (kernel calls), Data Layout Implementer (the `FluxStages` phase wiring and the device hook module of WP9) | WP9 device lists and gather; WP10 wall-state upload before B3; WP8 table upload | one `DIVERGENCE_PART_1` per stage per level; tests D-T1 and D-T2 |
| P2 Hook control state | mode, list length, range check stay host; list upload and value scatter as small device operations | Data Layout Implementer | WP9 (5-6 work-days ESTIMATE, already planned) | D-T1: empty list equals plain run bitwise |
| P3 `PREDICT_NORMAL_VELOCITY` wall loop (A1) | stays host; its output goes into the wall-state refresh | Data Layout Implementer (wall refresh), Integration Lead (order check) | WP10 | `U_NORMAL(_S)` on the device equal the host values (checksum, plan 7.10b) |
| P4 `WALL_LOOP_2` host path (B3) | gather `ZZP`/`TMP` at wall-adjacent cells, run the FDS loop body on the host over the wall tables for all species, scatter the face values; keep the thin-wall guard | Integration Lead (decision and hand-off measurement), Data Layout Implementer (host code in the derived copy) | wall tables on the host (today); a K1 scatter | wall-face values bit-equal to the host chain; hand-off bytes per stage measured against the chain time (item 9b-7) |
| P5 `WALL_LOOP_2` as a kernel (alternative to P4) | needs array sections with `lo:hi` in the emitter, wall species tables, `MAXLOC(B1%ZZ_F)` (second site for the D-051 note), `STORE_SPECIES_FLUX`, the enthalpy call | Generator Engineer | the W1/`WLIST` indirection; wall species tables on the device (WP8); D-051 sign-off for the second `MAXLOC` | bitwise on `dec2_obst`; mutants for the tie rule |
| P6 Output copies `STORE_SPECIES_FLUX` (B1) | recommended: not supported in the first device chain (assert), enable later as a K2 elementwise kernel or a host copy | Integration Lead (decision), Generator Engineer (kernel if wanted) | question 3 | output quantity test, if enabled |
| P7 `MAXLOC` kernels sign-off | `rho_d_maxloc_fix` (A11); plus B3's `MAXLOC` if P5 | Species and Combustion Lead (sign-off), Generator Engineer (NaN limit, tie test) | D-051 note; tie test with 3 or more species | mutants `maxloc_x_last_wins`, `maxloc_z_last_wins` fail the test |
| P8 Real inputs for the held-out kernels | dump `RHO_D` and the `D_Z` table in the scratch driver; rerun `WP3_DIFF` on the test machine's GPU; add the hook read-out as FDS truth for the part A outputs | Integration Lead | scratch driver patch (not applied to the driver tree) | `rho_d_interp`, `rho_d_dzd`, `wall_rho_d_dzdn`, `rho_d_maxloc_fix` outputs equal the read-out (±0 rule) |
| P9 Generator gaps in the chain | whole-array `RHO_D` forms (A3, A5), the second `D_Z_MAX` update (B8), `RTRM` scaling, `D_SOURCE`, stratification (B13) | Generator Engineer | generator classifier (DP1) | per-kernel bitwise tests |
| P10 Zone sums in the chain (B15) | per-cell terms on the device, serial addition in FDS order | Generator Engineer (terms), Integration Lead (order and the wall-flux sum) | D-053; which cases have `N_ZONE>0` (question 5) | per-cell terms bitwise; no atomics |
| P11 Timing and controls | D-T3 (scaled override differs), D-T4 (gather equals the same faces of the whole-box read-out), D-T5 (read-out, upload, scatter, syncs at 64^3 and 128^3) on a two-level case with 1 and 4 boxes | V&V Lead (controls), Integration Lead (timing) | a two-level case with at least 2 species (finding 1) | as in the flux-hooks review |
| P12 Wording | D-061 amended: the override point is inside part 1 of the chain | Chief Architect | | decision log entry |
| P13 Line references | the review and 9b-10 quote a shift of 9 after line 267; the files show 9 before the hook and 15 after (finding 2) | Legacy Mapper | | `port_merge_check` passes on both trees |

Dependencies in one line: P1 needs WP9, WP10 and (for P5) W1/`WLIST`; P4 needs nothing the generator owes; P7 gates counting `rho_d_maxloc_fix` as ported for the device chain.

## 4. Launch sequence per stage (with plan 7.12 in mind)

K2 launches block and do not wait for AMReX streams, so: (1) part A kernels (K2) complete before the next statement; (2) the gather is a K1 kernel on the development machine stream, so **one stream sync before it reads data a K2 kernel wrote is not needed (K2 is complete), and one stream sync after it is needed before the host copy** (a `cudaMemcpy` after a stream kernel is stale without a sync); (3) the list values copy to the device and the K1 scatter run on the stream; **a stream sync before the first part B K2 launch** (a K2 kernel does not wait for the stream); (4) part B K2 kernels. That is two syncs per stage and level, as in the review (finding 4); with several levels all part A runs first, then the fine-to-coarse override on the host, then all part B (D-061 phase order, `TimeLoop.H` comment of `FluxStages`). Overlap of boxes needs one host thread per box (plan 7.12, "streams per box"); not measured with the K1 gather active at the same time.

## 5. (d) Checks

Rules carried over: corrector-time inputs for all dump-based comparisons (plan 7.9c); equality is bit-identical after mapping -0 to +0 and any NaN fails (plan 7.10c, `zero_equal.h`, `zero_equal.py`); tolerance class per kernel from `kernel_registry.txt` (D-070: default bitwise; the registry lists only `cfl_wall_max` in the `libm` class, so every kernel of this chain is bitwise); device = host serial = host 4 threads, each device job once; no-FMA flags (K2 `nofma`, `--fmad=false`, host `-ffp-contract=off`) and, for the Intel host reference, `-O2 -fp-model=precise` (D-079).

| Piece | Check | Existing data that covers it (read from the plan and the harness) | Not covered |
|---|---|---|---|
| A4 `rho_d_interp`, A7 `d_z_max`, B4 `dp_div_heat`, B5 `del_rho_d_del_z`, A8 `rho_d_dzd`, B2 `h_rho_d_dzd`, A11 `rho_d_maxloc_fix` | device = host serial = host 4 threads, bitwise (±0 rule) | stage `WP3_DIFF` and `WP3_STAB` on `csmag_32` and four `dec2_obst` boxes, 155 of 155 outputs (plan 7.13a); inputs real except `RHO_D` (2e-5 x `RHO`) and the `D_Z` table, which are **synthetic**; no FDS truth for these outputs | FDS truth; `csmag_32` has uniform `ZZ` (all gradients zero), so only `dec2_obst` is a real test; corrector-pass values (the dumps are the predictor pass, step 2) |
| End of chain: `DEL_RHO_D_DEL_Z`, `DP` | compare with FDS | `WP1_DIV1` dump record (`p1_div1`): `DEL_RHO_D_DEL_Z` and `DS` after the divergence part, on `shunn3_32` and `dec2_obst` (plan 7.9, `extract_real.py`) | corrector pass; the contribution of B3 and the missing B8/B13 kernels is inside this number, so a difference cannot be attributed to one piece |
| A8/A10/A11 outputs at the hook | compare `RHO_D_DZDX/Y/Z` with FDS | **the DIF read-out of the existing hook** (mode 1) in a scratch driver run on `shunn3_32` / `dec2_obst`; not run yet (P8) | `RHO_D` and the `D_Z` table are not dumped, so the kernel chain cannot be started from real inputs yet |
| A10 `wall_rho_d_dzdn` | `B1%RHO_D_DZDN_F` for `NIC<=1` interpolated walls | wall-table dump `walltab` carries `RHO_D_DZDN_F(1:NS)` (WP1b draft patch); `dec2_obst` 576 external + 8 internal walls per box (plan 7.10a) | `NIC>1` branch: no harness case (the Legacy Mapper found only `race_test_1` and `race_test_4` and derived 2:1 inputs with 2 or more species); the stored walls' `ZZ_F`, `RHO_D_F`, `RDN`, `U_NORMAL` are not in the draft dump |
| B3 host path (P4) | wall-face values and `B1%RHO_D_DZDN_F` bit-equal to the host chain | none yet | all; needs the wall species inputs above |
| Whole split (P1) | D-T1 empty override = plain device run = host; D-T2 no-op override (list = read-out values) through the split chain; D-T3 scaled override differs; D-T4 gather = same faces of the whole-box read-out | the host hook tests `tests/run_flux_hook_check.sh` (empty, no-op and scaled lists on the host path) | any two-level device run; none exists |
| Sign-off P7 | tie test with 3 or more species and exact ties; mutants `maxloc_x_last_wins`, `maxloc_z_last_wins` | generator test `test/make_div1_tests.py`, `test/mutation_div1.sh`; I have not read whether the inputs contain exact ties | the real-data tie frequency |

Dump-time notes for any new dump of this chain: (1) wall tables must be taken after line 105 (finding 5); (2) `ZZ`, `TMP`, `RHO`, `MU` are not modified inside the chain, so an entry-time dump equals what the hook sees; (3) `DEL_RHO_D_DEL_Z` at entry is the previous stage's value and is zeroed at line 126 for NS>1, so the entry value is not an input of the chain.

## 6. (e) Open questions and effort

### Open questions (not determinable from the sources)

1. Is the D-051 sign-off for the species-sum kernel on file? I found no note naming `rho_d_maxloc_fix` (the later documents still list it as open). If B3 becomes a kernel, a second `MAXLOC` site (line 341, over a wall species table) needs the same note. Owner: Species and Combustion Lead.
2. Is the NaN limit of `rho_d_maxloc_fix` to be closed or documented (9b-3)? Owner: Generator Engineer.
3. Is `STORE_SPECIES_FLUX` required on the device chain? If not, an assert is enough and B1 drops out (recommendation). If yes, who accepts K2 for the copies? Owner: Chief Architect.
4. Can a hook override face coincide with a solid-wall face of `WALL_LOOP_2`? The transport role's lists are interface faces, which suggests no, but no source states it. A one-line assertion at list upload would settle it. Owner: transport role.
5. Which test cases have `N_ZONE>0` (zone-sum block 740-796)? Not read; the harness cases carry none of it in the listed data.
6. Where does `turb.f90` (lines 1348-1437) update `D_Z_MAX` relative to `DIVERGENCE_PART_1` and `CHECK_STABILITY`, and is that a host or a device loop in the AMR driver? Not determined.
7. Is the `NIC>1` branch (divg.f90:219, wall.f90:885) reachable in the AMR driver, or does the hook list replace it for every interface? If it is reachable it needs a case; none exists in the harness.
8. Does an existing two-level driver case (`TwoLevelRun`) have 2 or more species? Not read; needed for P11.
9. The second `D_Z_MAX` update (B8), `RTRM` scaling, `D_SOURCE` and stratification (B13): the marker list has no kernels for them by my search; the loop tracker may list them under other loop ids. Not checked.
10. Predictor and corrector: the dumps are the predictor pass of step 2. Corrector-pass inputs for this chain (`ZZ`, `RHO`, `DS`/`D`, wall `U_NORMAL`) are not dumped (WP1b list).
11. Does patch 0008's `POINT_TO_MESH -> POINT_TO_BOX` redefinition in `divg.f90` matter for the device chain's pointer setup (A0)? The device chain uses its own box pointers; whether `POINT_TO_BOX` is what the derived copy should call is a Data Layout Implementer decision.

### Effort (work-days, ESTIMATE; no measurement behind these figures)

| Item | Owner | Work-days |
|---|---|---|
| P1 launch sequence on the 2-species single-box cases (kernels already exist for A4, A7, A8, A10, A11, B2, B4, B5) | Integration Lead | 1.5 to 2 |
| P8 real `RHO_D`/`D_Z` dump, `WP3_DIFF` rerun with the read-out truth | Integration Lead | 1 to 1.5 |
| P4 host `WALL_LOOP_2` path and the hand-off measurement (9b-7) | Integration Lead with Data Layout Implementer | 1.5 to 2.5 |
| P11 controls and timing on a two-level case (shared with WP9) | Integration Lead, V&V Lead | 1 to 1.5 |
| `FluxStages` phase wiring for A-gather / scatter-B | Data Layout Implementer | 1 to 1.5 (device hook module itself is WP9) |
| P7 sign-off and NaN limit, host tests for the three NO-TEST kernels | Species and Combustion Lead 0.5; Generator Engineer 1 to 1.5 | 1.5 to 2 |
| Alternative to P4: P5 `WALL_LOOP_2` kernel | Generator Engineer | 3 to 4, plus the W1 dependency |
| P9, P10 generator gaps and zone sums | Generator Engineer | not estimated here (WP3 and the zone-sum package) |

Integration Lead share with the host path: about **5 to 7.5 work-days** (ESTIMATE), against the plan's 3 to 4 work-days for the whole of WP9b. The plan figure holds only if P8 is booked with the held-out kernel reviews, P11 with WP9, and P4 stays a host path of one to two work-days. Counting everything above in WP9b raises the device-AMR hooks total (10 to 12 work-days in the plan) by about 2 to 3.5 work-days (ESTIMATE; P11 is shared with WP9 and is not double counted in that range).
