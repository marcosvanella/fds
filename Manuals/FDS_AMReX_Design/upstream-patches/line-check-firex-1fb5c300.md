# Line-number re-anchor after FireX 1fb5c300 (src 8ad685f435)

Scope: dry-run apply of driver patch 0010, UP-0011, and the drafts of driver 0011/0012 against FDS-AMReX `Source/` at merge commit `8ad685f435` (FireX `1fb5c300fa` merged; `Source/` blobs identical at tip `ee38f45b6c`). Method as in `line-check-0010-0011.md` and `line-check-output-plan.md`: `git show 8ad685f435:Source/<file>` into a scratch tree, `patch -p1 --dry-run` (and one real apply in scratch for the patched call sites), then remove the scratch. Pre-merge first parent `9dc218d55b` is the "old" column. Nothing was built. `pres.f90` and the 0012 sites in `dump.f90` are unchanged by the merge.

Patches: `docs/upstream-patches/0010-main-amr-level0-nic-guard.patch`, `docs/upstream-patches/UP-0011-pressure-rhs-dump.patch`, `Source/driver/patches/0011-main-mesh-dumps.patch`, `Source/driver/patches/0012-dump-output-mesh-tables.patch` (drafts from src `9442836334`).

## Apply status

| Patch | File | Old line (patch hunk / 9dc218d55b) | New line (8ad685f435, unpatched) | Apply status | Notes |
|---|---|---|---|---|---|
| 0010 hunk 1 | main.f90 | 112 (`N_WRITTEN` decl) | 112 | clean | `AMR_NIC_*` decl; before the FireX shift |
| 0010 hunk 2 | main.f90 | 327 (after `INITIALIZE_MESH_EXCHANGE_1` loop) | 327 | clean | `STOP_CHECK(1)` + `NIC_CHECK` block; patched insert at 330-343 |
| 0010 hunk 3 | main.f90 | 2132 (`SUBROUTINE INITIALIZE_MESH_EXCHANGE_1`) | 2135 | offset +3 | `USE COMP_FUNCTIONS` / `AMR_N_BAD` locals; patched routine starts at 2151 |
| 0010 hunk 4 | main.f90 | 2371 (before `END SUBROUTINE`) | 2374 | offset +3 | NIC>1 guard; patched `END SUBROUTINE` at 2429 |
| UP-0011 main hunk 1 | main.f90 | 1650 (`PRESSURE_SOLVER_COMPUTE_RHS` loop) | 1726 | offset +76 | was +73 on pre-merge FDS-AMReX; FireX adds +3 in `MAIN_LOOP`; patched call `.FALSE.` at 1732 |
| UP-0011 main hunk 2 | main.f90 | 1669 (after solve `END SELECT`) | 1745 context / residual comment 1748 | offset +76 (+3 vs pre-merge) | patched call `.TRUE.` at 1757; residual comment then 1760 |
| UP-0011 pres hunk 1 | pres.f90 | 9 (`PUBLIC` list) | 9 | clean | merge does not touch `pres.f90` |
| UP-0011 pres hunk 2 | pres.f90 | 1076 (before `END MODULE`) | 1076 | clean | routine body unchanged by merge |
| 0011 hunk 1 | main.f90 | 12 (`USE FDS_AMREX_HOOKS`) | 12 | clean | add `OUT_MESH_DUMPS` |
| 0011 hunk 2 | main.f90 | 120 (`MODE==3` comment) | 120 | clean | comment only |
| 0011 hunk 3 | main.f90 | 129 (after `UPDATE_CONTROLS`) | 129 | clean | insert between 131 and 132 (see Role 1) |
| 0012 hunk 1 | dump.f90 | 2 (file-scope) | 2 | clean | `NMESHES_OUT` macro |
| 0012 hunk 2 | dump.f90 | 30 (`USE COMP_FUNCTIONS`) | 30 | clean | `USE FDS_AMREX_HOOKS` under `WITH_AMREX` |
| 0012 hunk 3 | dump.f90 | 402 (`FN_XYZ` tables) | 405 | clean | FireX does not move these lines; header uses post-hunk-1/2 line |
| 0012 hunk 4 | dump.f90 | 2053 (`XLEVEL` block) | 2053 | clean | `XLEVEL/YLEVEL/ZLEVEL` still at 2056-2060 |

UP-0011 on FireX tip `1fb5c300fa` alone: main.f90 hunks at offset +3 (1653, 1678); `pres.f90` clean. So the patch still applies to current FireX; the large +76 is the accumulated FDS-AMReX drift plus this merge's +3.

## Role 1: 0011 / 0012 insertion points at 8ad685f435

| What | File:line | Content |
|---|---|---|
| 0011 insert | main.f90:131-132 | between `CALL UPDATE_CONTROLS(T,DT,CTRL_STOP_STATUS,.FALSE.)` (131) and `CALL DUMP_GLOBAL_OUTPUTS` (132) inside `MODE==3` |
| 0011 USE | main.f90:15 | `USE FDS_AMREX_HOOKS, ONLY: ...` (add `OUT_MESH_DUMPS`) |
| MAIN_LOOP reference dump loop | main.f90:1219-1228 | `N_WRITTEN=0` (1219), five `WROTE_*` (1220-1224), `DO NM=...` / `DUMP_MESH_OUTPUTS` (1225-1228); was 1216-1225 before the merge |
| 0012 table allocs | dump.f90:405-451 | 35 `ALLOCATE` lines; unchanged by FireX |
| 0012 smokeview levels | dump.f90:2056-2060 | `XLEVEL/YLEVEL/ZLEVEL(0:2*NMESHES)`; unchanged by FireX |

## Build Chiefs: 0010 / UP-0011 at 8ad685f435

| Patch | Key anchors (unpatched 8ad685f435) | Patched call / insert (scratch apply) |
|---|---|---|
| 0010 | exchange call loop 326-328; `INITIALIZE_MESH_EXCHANGE_1` 2135-2377; `EWC%NIC` 2187; `STOP_CHECK` 1972-1996 | `STOP_CHECK(1)` + `NIC_CHECK` after loop (330-343); guard before `END SUBROUTINE` (2397-2427) |
| UP-0011 | `PRESSURE_ITERATION_SCHEME` 1677; `PRESSURE_ITERATION_LOOP` 1699; `PRESSURE_SOLVER_COMPUTE_RHS` 1726; solve `SELECT CASE` 1731; residual comment 1748; `pres.f90` PUBLIC 11-12, `END MODULE` region 1076 | `PRESSURE_SOLVER_DUMP(...,.FALSE.)` 1732; `PRESSURE_SOLVER_DUMP(...,.TRUE.)` 1757 |

## Loop ranges (confirm vs next-work r6 / impact note)

Checked at `8ad685f435`; agree with `docs/next-work.md` (r6) and `docs/amrex/upstream-1fb5c300-loop-impact.md`:

| Key | Range at 8ad685f435 |
|---|---|
| L0630 `BLOCK_CELL` nest | func.f90:5581-5627 |
| L0677 `PACK_CELL` | func.f90:5248-5266 |
| L0817 `CREATE_OR_REMOVE_OBST` | init.f90:4889-4900 |
| L0822 `REASSIGN_WALL_CELLS` | init.f90:4928-4947 |
| L0823 `REASSIGN_WALL_CELLS` | init.f90:4951-5065 |

## Cross-refs

- Prior checks: `line-check-0010-0011.md` (pre-FireX), `line-check-output-plan.md` (FireX delta section added there).
- README target lines for 0010 and UP-0011 updated for this merge.
