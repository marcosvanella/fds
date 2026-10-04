# 08 · Review of the newly landed wall kernels (four mass-flux wall nests, open-boundary `Q_RAD_IN`, `Q_RAD_IN` zero, `BR_ILW` table)

Owner: AMR Solid Phase Lead · Status: **v0.1, review of committed code, read only** · Upstream line numbers refer to the pinned FireX `36975d765f`; the radiation text calls the same lines 4965-4974 (L1243) after its own merge.
Not repeated here (unchanged since `07`): the verdicts for L1359, L1402, L0375, L0394, L0405. Section 6 lists only what changed in their tests since `07`.

**What was read, what was run.**
- Branches, read only through `git show` (no worktree edited): `s5-four` tip `c6cfffe762` (holds the four nests; 90 kernels), `s5-wall` tip `4fe0ed1cbb` (holds `wall_open_qin`, the `BR_ILW` table checks, changed wall tests; 88 kernels), `s5-gen` tip `bc36fdfc8d` (holds `rad_wall_qin_zero`; 86 kernels). None of the three is an ancestor of another for these commits, so the four nests and `wall_open_qin` are **not yet on one branch** (92 kernels once merged).
- Nothing was run: the shared lock is busy, the disk is tight, and the instructions were to prefer reading. One live log of another run was observed (`test_legacy_walltab`, first part, 82 PASS and 0 FAIL when read, mutant part not reached). Statements about pass counts are therefore not made; "covered" below means "present in the test source".
- Rules judged against: `04` (SP1 to SP4), the claim protocol in `amrex/loop-work-list.md` section 4 (items 6 and 7), `amrex/massflux-wall-table-spec.md`, `radiation/05-br-ilw-table-spec.md`.

## 1. Verdicts

| Kernel | Loop | Branch | Rule | Verdict |
|---|---|---|---|---|
| `mass_wall_zz` | L0880 wall part, mass.f90:93-188 | `s5-four` | table-checked face scatter (not SP2 to SP4) | **Accept** with G2, W1, W2, W5 |
| `mass_wall_rmw` | L0882 wall part, mass.f90:224-320 | `s5-four` | same | **Accept** with the same follow-ups |
| `spec_wall_zz` | L0401, divg.f90:1021-1084 | `s5-four` | same | **Accept** with G2, W1, W2, W5; same limiter caveat (M1 below) |
| `spec_wall_rmw` | L0403, divg.f90:1120-1183 | `s5-four` | same | **Accept** with G2, W1, W2, W5; one open limiter caveat (M1 below) |
| `wall_open_qin` (+ `_wl`) | L1243, radi.f90:4961-4970 | `s5-wall` | SP4 ordered sum, per-band temporary | **Accept** with R1, R2 |
| `BR_ILW` slot table and `check_br_index` | L1243, L1248 support | `s5-wall` | table contract | **Accept** with W3 |
| `rad_wall_qin_zero` | L1239, radi.f90:3887-3893 | `s5-gen` | each wall writes its own row | **Accept on the numerics; one defect in the device list (G1)** |
| `wall_connect_zones` | L0400, divg.f90:1301-1319 | `s5-wall` and `s5-gen` | idempotent constant store | **Accept** (not solid phase; added because it landed unreviewed) |
| `cell_adv_hs`, `cell_adv_zz` | L0399, L0406 cell loops | `s5-wall` and `s5-gen` | cell loops, no wall store | Not reviewed here (no wall store; Legacy Mapper's) |

Nothing blocks. SP1 (`CELL_COUNTER`), SP2 (ghost mirror) and SP3 (`DP` accumulate) have no new kernel in this batch. SP4 applies to `wall_open_qin` only (an ordered sum into one per-wall slot); none of the new kernels touches `OBSTRUCTION%MASS`, `D_SOURCE`, `M_DOT_PPP` or `MASS_FLUX_VAR`, so the SP4 reference package `wall_sums.py` remains the only holder of that rule.

## 2. The four mass-flux wall nests

### 2.1 What they store, and why no SP2 flag is needed

The nests write species or density face fluxes (`FX/FY/FZ`, `FX_ZZ/FY_ZZ/FZ_ZZ`), not solid-cell state. Each wall makes up to two stores:
- the wall-face store (`mass_wall_*` only): `F(IIG-1,..)` for `IOR>0`, `F(IIG,..)` for `IOR<0`, value `RHO_F*ZZ_F(N)` (or `RHO_F/MW_F` for the `rmw` nest), or `0` on the thin-obstruction branch (generated lines 3267-3299 of `s5gen_k2.F90` on `s5-four`);
- the off-wall overwrite: `F(II+1,..)` for `IOR=+1`, `F(II-2,..)` for `IOR=-1` (same on y and z), under the velocity sign and the `WALL_INDEX` guard (lines 3300-3363).

Upstream runs the loop under `!$OMP DO`. The generator accepts the scatter only on that authority ("wall scatter into FX accepted: upstream !$OMP DO above the loop (author's contract)", `generated/s5gen_report.md`, entries of the four nests). So the kernels carry **no flag, no gather and no atomic**; the safety of the launch rests on the table check `wall_checks.face_write_check` (`s5-wall` `wall_checks.py:166-200`, spec section 4 of `massflux-wall-table-spec.md`). Its rules, which I checked against the kernel text:
- (a) at most one non-thin wall-face store per face; several thin-path stores are allowed because every one writes the constant `0`; a thin and a non-thin store on one face is refused (order-dependent);
- (b) at most two off-wall stores per face, and two only with opposite `IOR` on the axis. This is sound: wall A (`IOR=+1`, cell `II_A`) reads `UU(II_A+1)` and wall B (`IOR=-1`, cell `II_B`) reads `UU(II_B-2)`; both name the same face only when they read the same element, and the tests `>0` and `<0` exclude each other (generated lines 3303 and 3313);
- (c) a face with both a wall-face store and an off-wall store is refused (the `WALL_INDEX` guard should have removed the off-wall store).

This is the "uniqueness proof at table-build time" form of claim-protocol item 7; it differs from SP2, where uniqueness cannot be assumed and the highest wall index must win. Here an input that breaks the rule is refused instead of resolved, which matches upstream (an order-dependent result under `OMP DO` is a race there too). The consequence is that the refusal must be a hard stop in the driver (W1).

### 2.2 Line-by-line against upstream

Compared the generated text with mass.f90:93-188, 224-320 and divg.f90:1021-1084, 1120-1183 of the pin:
- Skip of `NULL_BOUNDARY` first (3259, 3385). Thin test `SOLID_BOUNDARY .AND. .NOT.SOLID(II,JJ,KK) .AND. .NOT.EXTERIOR(II,JJ,KK)` with the live masks (3267, 3393). `INTERPOLATED` skipped in the wall-face store (3283). `INTERPOLATED` and `OPEN` skipped in the overwrite (3300).
- The six `IOR` cases use the right far cell and guard slot: `IOR=+1` reads `UU(II+1)` and `WALL_INDEX(II+1,..,+1)`; `IOR=-1` reads `UU(II-2)` and `WALL_INDEX(II-1,..,-1)`; same pattern on y and z (3303-3361).
- `ZS(0:3)` is `(/RHO_Z_P(II+1), RHO_Z_P(II+1:II+2), 0/)` for positive `IOR` and `(/0, RHO_Z_P(II-2:II-1), RHO_Z_P(II-1)/)` for negative `IOR`, as upstream with `DUMMY=0`.
- The `rmw` nest builds `ZZ_GET(M)=B1_ZZ_F(IW,M)` over `1..NS` where upstream copies `1:N_TRACKED_SPECIES`; equal today because `NS = N_TRACKED_SPECIES` (`N_PASSIVE_SCALARS` is never set). The driver rule "abort if `NS > N_TRACKED_SPECIES`" is in the spec and tested in `test_wall_tables.py:70`.
- `spec_wall_zz` (L0401) and `spec_wall_rmw` (L0403): upstream fills only three values of `Z_TEMP` in the divergence-side nests: `Z_TEMP(0:2,..)` for positive `IOR` (divg.f90:1037 and 1136, with the y and z twins) and `Z_TEMP(1:3,..)` for negative `IOR` (divg.f90:1047 and 1146, with twins), so the far element (`Z_TEMP(3)`, respectively `Z_TEMP(0)`) is stale. The mass nests do not have this: they fill four values with `DUMMY` (mass.f90:30, 141, 151). The kernels pad the far element with `0`. That is invisible to every limiter that does not read the far element; the MP5 limiter does. The Legacy Mapper already records that the MP5 bit tests wait for upstream patches UP-0001 and UP-0002 (`loop-work-list.md`, S2 line). **M1** below asks for that status to be stated in the sidecar comment of both entries.

### 2.3 Solid-phase interaction

- The nests read `B1_RHO_F` and `B1_ZZ_F` of the wall row. Under refinement the row is the owner face's row; both values come from the owner remainder of `WALL_BC` (`05`, section 3, step 6, and section 4: area-mean `TMP_F`, aggregated mass flux). The nests therefore need a place in the pass list **after** the owner remainder and **before** the gas advection that reads `FX/FY/FZ`. Nothing in the kernel enforces it (W5).
- The thin-obstruction branch decides on the live `SOLID` and `EXTERIOR` masks (3267). The masks change with burn-away, obstruction creation and `REASSIGN_WALL_CELLS` (L0823). The face-write check and `WALL_INDEX` must be rebuilt together with the masks, in the same refresh trigger as the `last_writer` flags of `07` F3. `check_all` does call `face_write_check` (`wall_checks.py:210-219`); the open point is only that the driver calls `check_all` on every refresh (W1).
- Under refinement a coarse thin obstruction needs the matching fine faces to carry the zero store too. That needs the fine-level wall tables of D-064 (4); until then the nests are for the single-level and per-level-box case, and `face_write_check` has only been exercised on box meshes (W2).

### 2.4 Test coverage against the claim protocol (by reading)

`test_legacy_walltab.py` (`s5-four`): registration checks (pinned arguments equal derived, golden present, no pointer scratch left, six one-face callee calls, wall tables flat, `WALL_INDEX` is the integer cell table, 90 kernels); bitwise against the verbatim upstream loop on derived types with real tables (limiters 0..5, species 1..3, six meshes of 6x5x4 with random boundary types, solid cells and zero-thickness faces, every mesh accepted by `check_all`); mutants `guard_ge`, `wall_index_slot`, `write_dropped`, `limiter_fixed`, and for the mass nests `thin_exterior_dropped` and `interpolated_test` (lines 132-138).
- Flag sets, threads: all six sets, 4 and 8 threads for the `omp` sets (line 108). Full run only (the quick run uses two sets).
- **Callee switch:** `bind` runs with the `O2omp` set only (line 102). The synthetic-table test `test_legacy_scratch.py` runs the full set-by-switch matrix (line 592), so the two tests together cover it. This is a documented split, not a gap, but the protocol wording is "all six flag sets, both callee switches"; accepted with the request to quote both result lines (W2).
- Counts printed in the live log: guard true and false per `IOR`, 648 cases, thin-obstruction walls 13 against 896 other walls, 757 wall-face value writes. The `IOR=0` column is empty as expected. Not printed: how many meshes contain the two situations that make the check non-trivial, a **two-cell gas gap** (rule (b) pair) and a **thin face written from both sides** (rule (a) shared thin). With six small meshes these may be absent. (W2)
- Vacuity of the order claim: the order tests live in `test_wall_tables.py` (`s5-wall`): the upstream loop run forward and reversed for every table the check accepts must agree, and for illegal tables of outcomes (a) and (c) must differ (counters `illegal-not-order-dependent` and `vacuous` are printed, lines 357-393); the legal suite contains the two-cell gap (`gap2`, expecting one opposite pair), the one-cell gap whose off-wall stores the guard removes (`gap1`), and a thin face written from both sides. That chain is sound. What is missing is the same forward-versus-reversed run **on the real-table meshes of `test_legacy_walltab`** for the generated kernel (W2), and a statement of how many of those meshes have a shared face.
- Negative checks: refusals of tables are in `test_wall_tables.py` (outcomes a, b, c, mixed thin and non-thin, ranges, inconsistent `WALL_INDEX`) and checker mutants. No generator-side negative is needed (no new attribute).
- Branch status: the four nests exist on `s5-four` only; `loop_claims.csv` still shows them "in progress" for the Legacy Mapper. They count as translated when the entries are on the generator branch head and the full run is reported.

## 3. L1243: `wall_open_qin` and the `BR_ILW` table

Kernel (`s5-wall` generated lines 1186-1206), against radi.f90:4961-4970 and spec section 7:
- Loop range `IW = 1, NWE` (external walls only), `CYCLE` unless `OPEN_BOUNDARY`; `B1_Q_RAD_IN(IW) = 0`; per band a temporary `DOT1 = 0` summed over the angles in ascending order and added to `Q`. This is the per-band grouping `Q + (sum over angles)`; a flat chain differs for more than one band (51 mismatches quoted in the table spec). The generated order equals the gfortran `SUM(...)` order: start at `+0`, left to right.
- Gather `BR_ILW(W_BR_INDEX(IW), M, IBND)`, table read only, `IBND` private (it is read by a later loop upstream, report line 886).
- SP4 reading: one writer per wall, one slot per wall, no shared target, no atomic: race-free by construction. `unique = true` in the sidecar (`s5_markers.toml:365-374`) is right for the store `B1_Q_RAD_IN(IW)`.
- Precondition `check_br_index` (`wall_checks.py:390-400`): every wall has `0 < BR_INDEX <= NBR` and slots are unique across all walls including `NULL` walls. Needed because the kernel launches over walls and two threads must never share a record.

Tests (`test_wall_brilw.py`, by reading): structure, builder checks (zero, negative, beyond dimension, duplicate, unoccupied slot, growth keeps old slots), bitwise against real `BOUNDARY_RADIA_TYPE` records with slots a random permutation inside `NBR >= NWE+NWI`, values of every class (wide magnitudes, `±0`, ties, cancellation, denormals), `NSB` 1/2/5/6, `NRA` 1/3/17/100/104, six flag sets, both callee switches, 1/4/8 threads, table never written, wall list in random order, mutants (flat chain, sum start, gather by `IW`, slot and angle swapped, wrong band, dropped condition or range, dropped private) and builder mutants. Not vacuous: the flat-chain mutant must fail for `NSB>1`, and the permuted slots make "gather by `IW`" fail.

- **R1 (compiler scope).** The bitwise claim is against the gfortran `SUM` order. A different compiler may vectorise a `SUM` reduction and change the grouping. Say in the kernel notes that the claim is per compiler, and let the V&V Lead say whether the oneAPI reference build of the open-boundary cases needs a tolerance class instead.
- **R2 (row ownership under refinement).** `OPEN` walls live on the domain boundary only. A fine level that touches the domain boundary has its own `OPEN` walls with their own slots (`radiation/01` L-2: each level's walls get their own `Q_RAD_IN`); the coarse walls it covers are not owner rows and must not be in the list that `wall_open_qin_wl` receives. The twin visits `WLIST_EXT` only (`s5-wall` generated line 2735), which is what the owner list needs. State this in the table spec; it needs no code change.
- Not checked by the kernel test: the scatter of `B1_Q_RAD_IN(IW)` back to `BOUNDARY_PROP1(B1_INDEX)`. It relies on unique `B1_INDEX` (func.f90:4226-4246, as the Radiation Lead's sign-off states), but no builder or check for it exists in `wall_checks.py` (W3).

## 4. L1239: `rad_wall_qin_zero` and the flag `B1_PRESENT`

Kernel (`s5-gen` `test/rad_kernels.golden:271-283`), against radi.f90:3887-3893:
```
DO IW = 1, NWI + NWE
   IF (W_B1_PRESENT(IW) == 0 .OR. W_BOUNDARY_TYPE(IW) == NULL_BOUNDARY) CYCLE
   IF (SF_TMP_GAS_FRONT(W_SURF_INDEX(IW)) <= 0._EB) B1_Q_RAD_IN(IW) = 0._EB
```
- **Flag definition.** The kernel keeps the `NULL_BOUNDARY` test as a separate term, so `W_B1_PRESENT` must mean only `B1_INDEX /= 0` (the table spec section 8 says exactly that, and the test builds it as `merge(1,0,B1_INDEX /= 0)`). `B1_INDEX` is never negative, so `/= 0` and `> 0` are the same. If the builder instead folded the `NULL` test into the flag (`B1_INDEX > 0 and not NULL`), the kernel result would not change, but the flag would then stop meaning "this wall has a `BOUNDARY_PROP1` row" for any later consumer. Keep the pure definition. Under refinement the flag stays a property of an owner row; records have no `B1_INDEX`.
- **Open boundary.** The loop zeroes every non-null wall with a row, so it also zeroes `OPEN` walls when `TMP_GAS_FRONT <= 0`; L1243 later overwrites those (`Q_RAD_IN = 0`, then the band sums), so the order L1239 before the sweeps before L1243 must be kept; both are inside `IF (UPDATE_INTENSITY)` upstream.
- **Surfaces with `TMP_GAS_FRONT > 0`.** They are neither zeroed here nor accumulated at radi.f90:4905 (`IF (SF%TMP_GAS_FRONT>0._EB) CYCLE`), and keep the value `EMISSIVITY*SIGMA*TMP_G**4` set in `NEAR_SURFACE_GAS_VARIABLES` (wall.f90:440, 451). The per-surface table gather by `W_SURF_INDEX` is correct for them because `SURF_INDEX` is a property of the face and all records of a face share it.
- **Test (`make_rad_tests.py`, `case_rad_wall_qin_zero`, by reading):** five wall counts from (0 internal, 1 external) to (37 internal, 300 external), 40 repetitions, 20 percent of the walls without a row, boundary type random 0..4, `B1_INDEX` a random permutation (each row used once, which is the real situation), `TMP_GAS_FRONT` drawn from `-1, -0, 0, 1e-300, 5, -1e-300`, coverage floors (rows changed, no-row walls, `NULL` with row, zeroed walls, each at least 100; the wall-row check that a wall without a row is not written), 1/4/8 threads, plus the mutant `.AND.` for `.OR.` in `test_mesh_rad.py:162`. Not vacuous.
- **G1: defect, device list.** The generated directive lists `S4_DEV((B1_Q_RAD_IN,SF_TMP_GAS_FRONT,W_B1_PRESENT,W_BOUNDARY_TYPE))` (`rad_kernels.golden:279`) but the kernel also reads `W_SURF_INDEX` (golden lines 273, 276, 282). On the device that array would not be mapped. The host tests cannot see this. Cause, from `s5gen.py:2363`: the device list is the set of "used" arrays, and a table used **only as the subscript of another table** is not counted. `wall_open_qin` does not have the problem (its list contains `W_BR_INDEX`). I checked the 17 radiation goldens and the 90 + 88 generated kernels for an array argument missing from the first `S4_DEV` list: only this kernel shows it (the few other hits are kernels with two directives, which a single-list check cannot judge). The Generator Engineer should fix the list rule and add a structural check "every array argument appears in the device list of its directive(s)" for all kernels.

## 5. Interaction with the solid-phase reads of `Q_RAD_IN` (owner-face rule)

Where the solid phase reads `Q_RAD_IN` (pin): `SOLID_HEAT_TRANSFER` front boundary (wall.f90:2280, 2311, 2877), `TMP_F` iteration (771, 800), the mass-flux block of `CALCULATE_ZZ_F` (1194, 1202: it divides by `B1%EMISSIVITY` to recover the incoming flux), the back side from the snapshot (2167), the thin-wall average of the two sides (487). Fine-level plan `05`: records take the owner's lagged `Q_RAD_IN` (OQ-S6, Radiation Lead to rule). Checks against the new kernels:
1. **Owner rows only.** L1239 and L1243 write the owner rows; fan-out copies to records after the radiation update. The record-level `CALCULATE_ZZ_F` block (1194/1202) uses the copy, so it sees the same value the owner has.
2. **The copy is exact only while the emissivity is the same for the owner and its records.** `Q_RAD_IN` is the **absorbed** flux: radi.f90:4917 adds `B1%EMISSIVITY*(INRAD_W+BBFA*EFLUX)`. In the pinned source `B1%EMISSIVITY` is set once per wall from the surface (`func.f90:4920`) and for particles (`part.f90:1590`); no wall routine assigns it later. All records of a face share the surface, so the copy is exact. If a later feature makes the emissivity vary per record (layer mixture, burn state), the owner would have to hold the **incoming** flux and the record would absorb with its own emissivity. I propose OQ-S6 = owner, with that condition written in; it needs the Radiation Lead's agreement.
3. **Pass order.** The thin-wall average at wall.f90:487 reads `B1M` and `B1P` `Q_RAD_IN` of the two sides; both owner rows must carry the post-radiation value before that pass (the thin-wall branch of `NEAR_SURFACE_GAS_VARIABLES` stays a separate pass after the wall pass, as already told to the Wall engineer).
4. **Refresh.** `W_B1_PRESENT`, `W_BR_INDEX` and the slot occupancy change when obstructions are created or removed; they belong to the refresh trigger of `wall-bc-translation-plan.md` section 11.3, with `check_br_index` re-run.
5. **Open boundaries** carry no solid phase (no records), so nothing in the solid phase reads the L1243 value; the refinement rule is the row-ownership note R2.

## 6. Changes to the earlier reviewed tests since `07` (branch `s5-wall`, `cf9f1b6399` to `4fe0ed1cbb`)

Closed: F1 (`check_all` now runs `assert_visc_mirror_ok` and the `bc_dp` copy-hazard check, `wall_checks.py:210-219`), F2 (wall-list twin of `wall_visc_mirror` run bitwise on random subsets with flags built over the visited walls, `test_wall_viscmirror.py`: new block "wall-list twin" and the mutant "flags built over all walls"), F4 (mutant "thin-wall guard moved above the DP statement" added, `test_wall_gather.py:184`), and mutants now also run with the loop reversed, so a dropped flag is observable (an earlier version could pass a flag mutant at one thread).
Open: F3 is stated only in the `check_all` docstring (flags rebuilt with the checks), with no refresh test; F6 (zero fill of `U_DOT_DEL_RHO_Z` is the caller's job) I did not find written.

Removed or weakened without a stated replacement (W4):
- `test_wall_nograd.py`: the four generator refusals for `last_writer` (not on a cell kernel, `true` instead of a name, lower-case name, name that collides with another wall table). The messages still exist at `s5gen.py:1188`, `:1379` and the collision check, and I found no remaining test that triggers them (`git grep` over `s5-wall`).
- `test_wall_nograd.py`: builder mutant "loop range `NWE+AUX` ignored". The loop-range contract (walls beyond `NWE + N_INTERNAL_WALL_CELLS_AUX` neither store nor overwrite) is what this mutant guarded.
- `test_wall_bcdp.py`: kernel mutant "NULL walls not skipped" and builder mutant "pointers shifted by one cell"; first and last cell mutants replaced by first and last plane mutants.
- `test_wall_viscmirror.py`: builder mutants "only SOLID wall cells count" and "every wall counts".
The test headers say these were "equivalent on box meshes". That may be true (an exterior ghost cell is reached by one wall only; edge ghost cells have no wall), but it is a property of box meshes. Under refinement fine-level wall tables may share targets differently. Ask for a one-line proof per removed mutant, or restore it.

## 7. Follow-ups per owner

**Wall Loops Engineer**
- **W1.** Make the refusal of `face_write_check` a hard stop in the driver table build and in every refresh (after obstruction create or remove, burn-away, `REASSIGN_WALL_CELLS`); define in the driver notes what happens then (abort, or a serial host loop; I recommend abort with the report, because the serial order of an invalid table is also undefined upstream).
- **W2.** `test_legacy_walltab`: print the number of meshes with a two-cell gap pair (rule b), with a thin face written from both sides, and with a face where the guard removed an off-wall store; require each at least 3 (raise the mesh count above six if needed). Run the generated kernels once with the wall loop reversed on those real-table meshes (equal to the forward run), and report both the `walltab` and the `scratch` result lines so the "both callee switches" requirement is shown by the pair. The quick `check_all` run on a fine-level table (two-level box mesh with a coarse thin face) is the first test for section 2.3, last item.
- **W3.** Add the missing builders and checks to `wall_checks.py`: `W_B1_PRESENT` (`B1_INDEX /= 0`, pure definition), a range check for `W_SURF_INDEX` against `0..N_SURF+N_SURF_RESERVED`, a uniqueness check for `B1_INDEX` over the walls that are written through `B1_Q_RAD_IN`, and the gather and scatter helper for `B1_Q_RAD_IN` (gather by `B1_INDEX`, scatter back by `B1_INDEX`).
- **W4.** Answer section 6: equivalent (one line) or restore, for each removed negative and mutant.
- **W5.** Pass-list notes: the four nests run after the owner remainder of `WALL_BC` and before the gas advection; `W_B1_PRESENT`, `W_BR_INDEX`, `WALL_INDEX` and the face-write check are rebuilt in the same refresh as the masks; `check_br_index` re-run on growth of `NBR`.

**Generator Engineer**
- **G1.** Device list rule: count a table used only as an index of another table (`W_SURF_INDEX` in `rad_wall_qin_zero`); add the structural test "every array argument is in the device list" for all kernels, including the radiation goldens.
- **G2.** For scatters accepted on "upstream `!$OMP DO` (author's contract)", let the sidecar entry name the justification that makes the acceptance checkable (for example a comment or attribute citing `face_write_check`), so claim-protocol item 7 can be audited from the sidecar. Today the four nest entries (`s5_markers.toml`, from "Wall nests with scratch pointers") do not mention the table check.
- **G3.** Merge `s5-four` and `s5-wall` onto one branch so the 92 kernels are counted together, then regenerate once (the golden and generated files differ between the branches).

**Legacy Mapper**
- **M1.** State in the sidecar comments of `spec_wall_zz` and `spec_wall_rmw` (and in the S2 note) that the divergence-side nests leave one far `Z_TEMP` element unset upstream, the kernels pad `0`, and MP5 equality depends on patches UP-0001 and UP-0002. Merge note: `loop_claims.csv` should show the four nests as on `s5-four` until the merge (G3).

**Radiation Lead**
- **R1.** Kernel note: the L1243 bitwise claim is against gfortran `SUM` order; ask the V&V Lead whether the oneAPI reference needs a tolerance class for the open-boundary cases.
- **R2.** Add the row-ownership sentence to `05-br-ilw-table-spec.md`: covered coarse `OPEN` walls are not in the owner list; `W_BR_INDEX` slots are unique over owner rows.
- **R3.** Rule OQ-S6: `Q_RAD_IN` belongs to the owner row, copied to the records after each radiation update; valid while `B1%EMISSIVITY` is per surface (`func.f90:4920`). If the radiation design ever varies emissivity per record, hold the incoming flux in the owner row instead.
- **R4.** Keep `W_B1_PRESENT` as `B1_INDEX /= 0` only (section 4); the kernel owns the `NULL` test.

## 8. Decisions needed

None blocks this week. Two are worth recording:
1. What the driver does when `face_write_check` refuses a table (W1). Proposed: abort with the report. Owner: Chief Architect with the Wall Loops Engineer.
2. OQ-S6, owner versus record for `Q_RAD_IN` (R3). Proposed: owner, condition on emissivity as in section 5. Owner: Radiation Lead.
