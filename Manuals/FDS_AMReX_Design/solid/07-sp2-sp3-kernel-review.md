# 07 · Review of the SP2 and SP3 kernels (wall engineer's branch `s5-wall`)

Owner: AMR Solid Phase Lead · Status: **v0.1, review of committed code** · Reference tree for upstream line numbers: FireX `36975d765f` (the sidecar and the tests cite this revision; the working tree is 9 to 12 lines higher).
Reviewed: branch `s5-wall`, tip `cf9f1b6399`, read only (worktree `src-s5wall`, nothing run, nothing edited). Files: `amrex/s4_mass/s5_gen/s5_markers.toml` (entries from line 331), `generated/s5gen_k2.F90`, `wall_checks.py`, and the tests `test_wall_viscmirror.py`, `test_wall_nograd.py`, `test_wall_gather.py`, `test_wall_bcdp.py`, `test_wall_specadv.py`.
Rules reviewed against: SP2 to SP4 of `04-blocked-loop-signoff.md` and the family text of `docs/amrex/blocked-loop-families.md`.

**What was not done.** The shared generator lock is busy, so no test was run. Statements below about test content come from reading the test source; statements about pass counts are not made. The sign-off is therefore "accept on reading", and each kernel needs its own full run (six flag sets, both callee switches, 4 and 8 threads, mutants) from the Wall engineer before the loop counts as translated (claim protocol, `loop-work-list.md` section 4).

## 1. Verdicts

| Kernel | Loop | Rule | Verdict |
|---|---|---|---|
| `wall_visc_mirror` (+ `_wl`) | L1359, velo.f90:355-363 | SP2 | **Accept**, with three follow-ups (F1 to F3) |
| `wall_uvw_nograd` (+ `_wl`) | L1402, velo.f90:3414-3450 | SP2 pattern (O3) | **Accept** |
| `wall_corr_kdtd` | L0375, divg.f90:532-554 | SP3 | **Accept**, one test gap (F4) |
| `wall_bc_dp` | L0394, divg.f90:1574-1604 | P3 sign-off (pressure lead), SP3 analogue | **Accept**, one note (F5) |
| `wall_spec_adv2` | L0405, divg.f90:1237-1266 | S3 sign-off (species lead), SP3 analogue | **Accept**, one caller condition (F6) |

Nothing blocks. Note on scope: the SP3 rule of `04` (list excludes NULL, INTERPOLATED and OPEN walls; thin-wall `CYCLE` guards only the `KDTD*` stores) is the rule for **L0375**, which is `wall_corr_kdtd`. L0394 and L0405 are different loops (not SP3); they are checked in sections 4 and 5 against their own sign-offs and against the same ordering principle.

## 2. SP2: `wall_visc_mirror` (L1359) and `wall_uvw_nograd` (L1402)

### 2.1 Does the kernel implement per-target highest-wall-index gather?

Yes, in the "flag" form, not a literal gather. Evidence:

- Sidecar entry `last_writer = "W_LAST_VISC"` (`s5_markers.toml:344-355`, comment lines 341-343). The generated kernel starts every wall with `IF (W_LAST_VISC(IW) == 0) CYCLE` (`generated/s5gen_k2.F90:1024`; twin at 2550, after the `WLIST_EXT/INT` indirection). A wall keeps the flag 1 only if no later wall stores to the same target, so for each target exactly the highest wall index stores. This is the same result as "gather per target, take the highest wall index" and is race-free without atomics: the flagged walls store to distinct elements.
- The flag builder is `wall_checks.last_writer_flags` (`wall_checks.py:236-246`): `top[k] = iw` over ascending `iw`, so the last assignment wins, and flag 0 for every wall that stores and is not the top. Target of L1359 is `visc_mirror_target` (`wall_checks.py:257-261`): the wall cell `(II,JJ,KK)` if it is in the live SOLID or EXTERIOR mask, else none. This is the guard of velo.f90:358. No distinct-target assertion is made: `store_conflicts` is informational only (`wall_checks.py:248-255`). That matches the SP2 rule "do not assert uniqueness".
- Both stores of the loop (`KRES` and `MU`) sit under the same flag and the same guard (kernel lines 1025-1028), so they cannot take different winners.
- The generator derives `unique` from `last_writer` (`s5gen.py:1854`); the generated code has no atomic or critical construct (checked by text search of `s5gen_k2.F90`).
- A residual hazard the flags do not cover is a cell that one wall stores and another wall reads as its gas cell (the serial loop would feed the stored value on). `wall_checks.visc_mirror_hazards` finds it (`wall_checks.py:263-273`). In a valid box mesh it cannot occur (a gas cell is never a SOLID/EXTERIOR wall cell). The check exists and is tested, but see F1.

L1402 follows the same pattern: `last_writer = "W_LAST_UVW"` (`s5_markers.toml:331-343`), kernel line `IF (W_LAST_UVW(IW) == 0) CYCLE` first (generated lines 969 on), builder target `uvw_nograd_target` (`wall_checks.py:221-234`) that applies the loop's own guards (SOLID, NULL, MIRROR; NULL only if the gas or wall cell is solid; IOR in plus or minus 1 to 3), and the loop range `NWE + N_INTERNAL_WALL_CELLS_AUX` is passed as `n_exec` so walls beyond the executed range neither store nor overwrite. Correct: the value stored depends on the wall (`VEL_N = UN_WALLS(IW)`), so the idempotent form would be wrong and the last-writer form is needed.

### 2.2 Does the test detect a wrong winner?

For L1359 (`test_wall_viscmirror.py`):

- Shared targets occur in the inputs: the Python meshes are required to contain solid cells stored by several walls in at least a third of the meshes (line 208) and walls that fail the guard in at least a quarter (210). The Fortran meshes use the same generator (`make_meshes`, seed 29, 56 meshes; 16 in quick mode).
- `KRES` and `MU` are filled with independent random values per cell (`mesh_text`, line 48), so two walls with different gas cells store different values to a shared target; a wrong winner changes the bits. The driver counts a mesh as vacuous if the reference loop changed nothing (lines 121-122) and the test requires zero vacuous meshes (line 283).
- Reference is the verbatim upstream text lines 355-363 compiled into the harness on derived types (lines 27-29, 99).
- Order test, not vacuous: the loop run backwards with the flags equals the serial loop (line 289); run backwards with all flags 1 it must differ (line 293).
- Mutants of the winner: builder mutants "lowest wall index wins", "all flags 1", "every wall counts (thin walls too)", "only SOLID counts", "only EXTERIOR counts", "target is the gas cell" (lines 188-195, run at 304-308); kernel mutants include flag test inverted or dropped, `.AND.` for `.OR.`, EXTERIOR or SOLID test dropped, wrong source array, swapped gas indices (163-173). All are expected to be caught; whether they are caught in a run is for the Wall engineer's run to show.
- Python-level model of the hazard: injected hazards are found and refused, and an injected hazard really makes the two forms differ, so the check is not vacuous (lines 241-251).

L1402 has the same structure in `test_wall_nograd.py` (builder mutants lines 214-220 include "lowest wall index wins", NULL guard ignored, wrong face, MIRROR not counted, loop range ignored; line 196-197 `flags_lowest`).

### 2.3 Claim-protocol coverage (by reading)

- Flag sets: the bitwise loop in `test_wall_viscmirror.py` runs `WL.FLAGS`, which is the six sets `O0`, `O2`, `O0omp`, `O2omp`, `O2omp_off`, `O2omp_dpd` (`test_wall_wlist.py:22-23`), in the full run (line 277); the quick run uses `O0omp` only.
- Callee switches: `dpd` and `bind` in the full run (line 278). The loop has no callee, so the switch only exercises the kernel entry form.
- Threads: 4 and 8 for the `omp` sets, 1 for the others (line 279), as `generator-howto.md` section 4 states.
- Flags: `-ffp-contract=off` is in `BASE` (`test_wall_wlist.py:24`).
- Mutants: the full run only (line 295). Negative checks for the `last_writer` attribute are in `test_generator_checks.py` and the gather/wlist negatives; not re-read here.
- Gap, the wall-list twin: the twin `wall_visc_mirror_wl` is only checked **structurally** (build, manifest, header text, lines 264-266). No bitwise run feeds it a list. See F2.

### 2.4 Follow-ups for SP2

- **F1 (wire the hazard check).** `assert_visc_mirror_ok` is called only in the test; `wall_checks.check_all` (lines 210-216) runs ranges, `WALL_INDEX` and the face-write check but not the viscosity hazard. The driver's table check should call it after every build and every refresh. Same for `gather_foreign_hazards` (L0394), which `check_all` also omits. Cheap; no kernel change.
- **F2 (list subset).** The flags are built over all walls `1..NWE+NWI`. A wall-list twin that runs over a subset (owner faces only per ruling W1, or a level's walls) is correct only if the highest wall of each target is in the list, or if the flags are built over the list. State the contract: `last_writer_flags` must be built over exactly the walls the kernel visits. Add a bitwise test of the `_wl` twin with a non-identity list (a list that omits the top wall of a target must give the serial result of the *listed* walls). Until then the twin is "structure only".
- **F3 (flags depend on the live masks).** The target test uses the live SOLID and EXTERIOR masks, so the flags change whenever a wall cell changes type (obstruction creation or removal, `REASSIGN_WALL_CELLS`, open/close). The flag tables belong to the table-refresh trigger of `wall-bc-translation-plan.md` section 11.3. Make it explicit in the driver notes, and add a refresh test: after a mask change, flags rebuilt from scratch equal flags from the old build only where nothing changed.

Also for the driver: L1359 reads the `MU` that L1358 (SP1) writes in the gas cell, so the order SP1 then SP2 must be kept. Nothing in the kernel enforces this; it belongs in the pass list.

## 3. L0375: `wall_corr_kdtd` against SP3

SP3 conditions from `04`: (1) the list excludes NULL, INTERPOLATED and OPEN walls; (2) the thin-wall `CYCLE` on `IOR<0` guards only the `KDTD*` stores, not the `DP` term; the gather starts from the `DP` already in the cell.

Evidence, `s5_markers.toml:358-371` and generated lines 1033-1092:

- Two regions, as designed. Region 1 (lines 1046-1069) is the per-wall part: NULL and INTERPOLATED `CYCLE`; OPEN stores `K_G` (average of the two cells) and cycles; others store `K_G = KP(gas cell)`; then the thin guard `IF ((W_THIN(IW)/=0) .AND. BC_IOR(IW) < 0) CYCLE` and the six `KDTD*` zero stores. The `DP` accumulate is not in this region. Same order of statements as divg.f90:534-552.
- Region 2 (lines 1070-1092): per gas cell `ICELL = 1, IBAR*JBAR*KBAR`, walls `W_CSR_GAS_PTR(ICELL)..PTR(ICELL+1)-1` taken from `W_CSR_GAS_LIST` in ascending wall index (builder `csr_cell_walls`, `wall_checks.py`, key "gas"). In the loop: NULL and INTERPOLATED `CYCLE` (condition 1), OPEN `CYCLE` (condition 1: the open wall stored only `K_G` upstream), then `DP(...) = DP(...) - (AREA_ADJUST*Q_CON_F*RDN - Q_LEAK)`, which starts from the `DP` already in the cell (condition "gather starts from the DP in the cell"; mutant "starts from zero" exists, `test_wall_gather.py:179`).
- Condition 2 holds by statement order: the thin guard in region 2 comes **after** the `DP` statement (line 1080 against 1079). Both sides of a thin wall therefore subtract from their own gas cells.
- Equality with the serial loop: the `DP` chain of a cell receives its terms from the walls of that cell in ascending `IW` in both forms; chains of different cells are independent; the `KDTD` zero stores write a constant, so the `idempotent` marker is right.
- The list is not pre-filtered; the exclusions are runtime `CYCLE`s inside region 2. Equivalent to "the list excludes" and simpler to build. Accepted.

Mutants (`test_wall_gather.py:174-190`) cover: OPEN accumulates, NULL/INTERPOLATED test dropped, INTERPOLATED not skipped, sign, start-from-zero, `Q_LEAK` sign, walls visited backwards, first/last cell skipped, wrong list entry, OPEN factor of `K_G`, `K_G` from the wall cell, thin walls processed on both sides (region 1), wrong `KDTD` faces, NULL stores in region 1. Builder mutants: descending order, order by IOR, external walls only, shifted pointers.

- **F4 (test gap).** No mutant moves the thin guard of region 2 **above** the `DP` statement (the exact SP3 condition 2). Today a kernel with that mutation would still pass if no mesh has a thin wall with `IOR<0` and a non-trivial `DP`; the random meshes do include thin walls (`test_wall_gather.py` header), so it would probably be caught, but the mutant should exist. Also, region 2 carries the dead statements `IF (thin) CYCLE` and an empty `SELECT CASE` (lines 1080-1088): harmless, but they show that the generator keeps statements whose stores moved to region 1. A cleanup in the generator is optional.

## 4. L0394: `wall_bc_dp`

Sign-off text (P3, Pressure Solver Lead): cell-to-wall list as a CSR gather in ascending `IW`; the subtract pass before the copy pass; `BC_INDEX` against `B1_INDEX`. My SP3 rule does not apply to this loop (its walls are keyed by the wall cell, and the subtract list is "SOLID_BOUNDARY with SOLID wall cell", not the SP3 exclusion list). Review against the P3 text and the same ordering principle:

- Key is the wall cell: `W_CSR_WC_PTR/LIST` (`wall_checks.csr_cell_walls(key="wc")`), loop `ICELL = 1,(IBAR+2)*(JBAR+2)*(KBAR+2)` (generated lines 1107-1109). Each cell visits its walls in ascending `IW`.
- NULL walls skipped first (line 1110), SOLID case: `IF (.NOT. SOLID(II,JJ,KK)) CYCLE` (1113), `UN_P` from `U_NORMAL_S` when `PREDICTOR`, else `U_NORMAL` (1114-1118), six `IOR` cases with the same factors as divg.f90:1589-1599 (`R(II-1)` for IOR -1, `RDY` for plus or minus 2, `RDZ` for plus or minus 3), copy case for OPEN, MIRROR and INTERPOLATED (1133-1134).
- Single pass instead of "subtract pass, then copy pass". This is equivalent to the serial loop, and in one respect closer to it than two passes: a wall cell that is reached by both a subtracting wall and a copying wall keeps the ascending-wall order of the serial loop. What must hold is that a copy reads a gas cell that no wall stores (otherwise the copy could read a partly updated `DP`); `gather_foreign_hazards` / `bc_dp_stores` check this, and the test builds a mesh with a thin MIRROR wall where the check fires and the kernel really differs from the serial loop (`test_wall_bcdp.py:210`, 227-228, 262). Not vacuous.
- `B1 => BOUNDARY_PROP1(WC%BC_INDEX)` is read through the flat `B1_*` tables; `check_b1_bc_index` refuses `BC_INDEX /= B1_INDEX` (test line 225-226). Matches the family text (patch 0004 note).
- Mutants (`test_wall_bcdp.py:170-181`): backwards, planes skipped, SOLID test dropped, PREDICTOR branches swapped, sign, `R(II)` for IOR -1, `RDX` for IOR 2, `RDY` for IOR -3, copy from the wall cell, MIRROR/INTERPOLATED not copied; plus builder mutants and negative checks (line 284).
- Coverage by reading: 12 flag/callee combinations, 1/4/8 threads (header lines 11-12).
- **F5.** Same as F1: `gather_foreign_hazards` is not part of `check_all`; the driver must run it after each table build. Also state in the kernel notes that the single-pass form replaces the "two passes" wording of the family text, so the family text is updated by its owner.

Verdict: accept.

## 5. L0405: `wall_spec_adv2`

S3 sign-off condition (Species & Combustion Lead): the per-cell wall order is a fixed key (ascending wall index) independent of the box decomposition; the result array is zeroed first.

- Key is the gas cell: `W_CSR_GAS_PTR/LIST`, loop `ICELL = 1, IBAR*JBAR*KBAR` (generated lines 1155-1157), walls ascending (builder `csr_cell_walls`, key "gas"; test builder mutants "descending", "by IOR", "internal before external", "pointers shifted": `test_wall_specadv.py:173`).
- Only NULL walls are skipped (line 1158), as upstream (divg.f90:1243). `UN` per boundary type as upstream (1246-1262): default case from `UU`, `VV`, `WW`, SOLID from `U_NORMAL_S` or `U_NORMAL` by `PREDICTOR`/`CORRECTOR`, INTERPOLATED from `UVW_SAVE(IW)`. The accumulate is `U_DOT_DEL_RHO_Z(gas) - SIGN(1,IOR)*DU*RDN` in the upstream expression order. `UVW_SAVE` has extent `NWE`: INTERPOLATED walls are external only, so the read is in range.
- `UN` is a private scalar. Upstream leaves `UN` stale if neither `PREDICTOR` nor `CORRECTOR` is true at a SOLID wall; main.f90:406-407, 748-749, 928-929 show exactly one of them is always true at this call, so no difference.
- Mutants (`test_wall_specadv.py:141-157`): order backwards, first and last cell skipped, NULL not skipped, sign, `SIGN` argument, `RDN` dropped, species index of `ZZ_F` and `ZZP`, `RHO_F` replaced by `RHOP`, predictor/corrector branches, INTERPOLATED source, wrong `IOR` reads.
- **F6 (caller condition).** The zero fill `U_DOT_DEL_RHO_Z = 0` (divg.f90:1235, `U_DOT_DEL_RHO_Z=>WORK7`) is outside the kernel range 1237-1266. The kernel is bitwise equal to the loop only if the caller zeroes the array first. This is the S3 sign-off condition; put it in the pass list and in the kernel report text, and make the integration test start from a non-zero array to show that the zeroing is the caller's job.

Verdict: accept.

## 6. Requests, in order of weight

1. Full runs of the five test files under the lock (six flag sets, both callee switches, 4 and 8 threads, mutants), reported with counts; I will check the pass counts when they arrive.
2. F2: bitwise test of the `_wl` twins with a non-identity list, and the contract that flags are built over the visited walls.
3. F1/F5: put `assert_visc_mirror_ok` and `gather_foreign_hazards` into the driver's table checks.
4. F3/F6: write the mask-change refresh and the zero-fill precondition into the driver notes.
5. F4: one extra mutant for the thin guard of region 2.

None of these needs a decision; all are for the Wall engineer. The family text for P3 (two passes) should be updated by the Legacy Mapper once F5 is agreed.
