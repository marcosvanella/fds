# O2 edge and exchange loops: what is translatable without the neighbour-mesh branches

Scope: the loops L1363, L1364, L1366, L1367, L1399 of the family O2 (`blocked-loop-families.md`). Ruling D-055 applies: no GPU kernels for the neighbour-mesh branches (`NOM>0`) of
`VISCOSITY_BC` (L1399), `VELOCITY_BC` (L1367, normal ghost face) and `NO_FLUX` (L1364, fill of `H`). This note reads each loop in the pinned upstream sources (`36975d7`) and states which part is
still translatable, with effort and risk. It is analysis only; the only code behind it is the L1366 work listed in section 3.

## 1. Result in one table

| loop | upstream lines | modelled time (%) | iterations that do work on a single mesh | translatable part without `NOM>0` | verdict |
|---|---|---|---|---|---|
| L1363 `MATCH_VELOCITY_FLUX` | velo.f90:2891-3015 | 0.994 | none: the loop skips every wall that is not `INTERPOLATED_BOUNDARY`, and an interpolated wall is a wall shared with another mesh | nothing | host; no kernel (D-055) |
| L1364 `NO_FLUX`, fill of `HP` | velo.f90:1376-1400 | 0.426 | none: `IF (NOM==0) CYCLE` | nothing | host; no kernel (D-055) |
| L1367 `VELOCITY_BC`, normal ghost face | velo.f90:1858-1897 | 0.795 | none: `IF (EWC%NOM==0) CYCLE` | nothing | host; no kernel (D-055) |
| L1399 `VISCOSITY_BC` | velo.f90:514-550 | 0.710 | none: `IF (EWC%NOM==0) CYCLE` | nothing | host; no kernel (D-055) |
| L1366 `NO_FLUX`, face values at walls | velo.f90:1463-1559 | 0.056 | all non-null walls | all of it; `NOM` is only a flag here | translatable; kernel built, waits for the merge and for acceptance |

The four loops L1363, L1364, L1367, L1399 consist entirely of the neighbour-mesh branch: on a single mesh they do no arithmetic, and their modelled time (2.925 % together) exists
only in multi-mesh cases. The default of the driver is the AMReX ghost fill (`EXTERNAL_GHOSTS_FILLED`, upstream patches 0003 and 0004), which also skips the `NOM>0` branches of
VISCOSITY_BC, NO_FLUX (the `H` fill) and VELOCITY_BC. So no part of these four needs a kernel while that default holds; what stays on the host is the loop test itself, which costs
one pass over the external walls with a table lookup. Nothing needs to be translated to keep results bitwise, and no non-`NOM>0` remainder exists to carve out.

### 1a. Status re-check against the current work list (docs commit 4f28cbd) and the pinned sources

The claim of section 1 was checked again. The five loops were re-read in the pinned upstream file (`36975d7`, `Source/velo.f90`) at the cited lines, and the work list and the claim register were compared with the note. Nothing changed:

- L1399 (514-550), L1364 (1376-1400) and L1367 (1858-1897) each begin the wall iteration with `IF (EWC%NOM==0) CYCLE` (L1364: `IF (NOM==0) CYCLE`); everything after it reads `OMESH(NOM)`. On a single mesh no iteration does any work.
- L1363 (2891-3015) skips every wall that is not `INTERPOLATED_BOUNDARY` and then uses `OMESH(NOM)` and `MESHES(NOM)` with `NOM = EWC%NOM`; an interpolated wall always has a neighbour mesh, so the loop is the neighbour-mesh branch as a whole. In the AMReX build the routine also returns early under `EXTERNAL_GHOSTS_FILLED` (the shared-face flux match is done by the driver), and for one mesh (`NMESHES==1`).
- L1366 (1463-1559) still reads no neighbour-mesh array: `NOM` only selects the branch and skips internal null walls.
- The work list (`loop-work-list.md`, generated tables) and the claim register (`loop_claims.csv`) still show L1363, L1364, L1367, L1399 as `claimed` (note "O2 neighbour-mesh loop") and L1366 as `in progress` (the `EXTERNAL_WALL(IW)%NOM` designator is on the s5-wall branch, the loop is not translated). In the generator worktree no committed marker file holds a `wall_no_flux` entry yet; the work is being done there by another worker. This note does not change the claim register.

Status rows (to copy into the work list when its owner next refreshes the hand-written part; the generated tables are not edited by hand):

| loop | routine | translatable? | why | status now |
|---|---|---|---|---|
| L1363 | `MATCH_VELOCITY_FLUX` | no | only `INTERPOLATED_BOUNDARY` walls, all neighbour-mesh (`NOM>0`); no kernels for `NOM>0` under D-055 | host; stays claimed, no kernel |
| L1364 | `NO_FLUX`, fill of `HP` | no | `IF (NOM==0) CYCLE`: only the `NOM>0` branch remains (D-055) | host; stays claimed, no kernel |
| L1367 | `VELOCITY_BC`, normal ghost face | no | `IF (EWC%NOM==0) CYCLE`: only the `NOM>0` branch remains (D-055) | host; stays claimed, no kernel |
| L1399 | `VISCOSITY_BC` | no | `IF (EWC%NOM==0) CYCLE`: only the `NOM>0` branch remains (D-055) | host; stays claimed, no kernel |
| L1366 | `NO_FLUX`, face values at walls | yes | no neighbour-mesh array is read; `NOM` is a flag (table column `EW_NOM`); needs the thin-wall uniqueness check | in progress in the generator worktree; not yet in a committed marker file |

## 2. The four neighbour-mesh loops: what each would be, if the ruling is ever lifted

Common shape of all four: one iteration per external wall with `NOM>0`; an overlap range `IIO_MIN..IIO_MAX, JJO_MIN..JJO_MAX, KKO_MIN..KKO_MAX` of the neighbour mesh; an inner K,J,I sum over
that range with no intrinsic and no reordering; a division by the cell count (L1364, L1367, L1399) or a weighted sum divided by the face area inside each term (L1363); a store into the
own mesh at a place fixed by the wall's cell and `IOR`. They would be one-thread-per-wall kernels with the inner triple loop sequential in the thread, which keeps the sum order
and is bitwise equal to the source. They need the neighbour data as a flat receive buffer in the upstream K,J,I order (exchange-buffer layout ADR-001) and three per-wall columns (`EW_NOM`,
overlap count, offset into the buffer).

| loop | reads from the neighbour | writes (own mesh) | notes |
|---|---|---|---|
| L1363 | `FVX/FVY/FVZ` and the cell sizes `DX, DY, DZ` of the neighbour; `DA_OTHER` is itself a sum over the overlap | `FVX(0 or IBAR,JJ,KK)` etc. with `0.5*(own + other)` | each term is `F*DY*DZ/DA_OTHER`, evaluated left to right inside the sum; the two cases per axis differ in the index shift (`IIO-1`) |
| L1364 | `H` or `HS` (by `PREDICTOR`) | `HP(II,JJ,KK)` (ghost cell of the wall) | `HP` is an alias of `H` or `HS` |
| L1367 | `U/V/W` or `US/VS/WS` (by `APPLY_TO_ESTIMATED_VARIABLES`), one component by `IOR` | `UU/VV/WW` at the face next to the ghost cell | followed in the source by `CC_RESCALE_OMESH_STAGGERED_TO_CUTFACE` for `CC_IBM`, outside this loop |
| L1399 | `MU`, `KRES`, and `D` or `DS` | `MU`, `KRES`, `D` or `DS` at the ghost cell | three sums in one pass |

Writes are one cell (or face) per wall; the targets of different external walls are different cells, so the kernels would be race-free. The table builder should assert it. If
and when a kernel is wanted: about 1 day each with the table, 3 to 4 days for the first (it sets the exchange pattern), and a test with overlap counts above 1 (coarse-fine) for the sums. The
risk is only the buffer layout and its refresh each step; the arithmetic is simple.

## 3. L1366: the part that is translatable, and the state of the work

L1366 loops over all external and internal wall cells, skips `INTERPOLATED_BOUNDARY` and `OPEN_BOUNDARY`, and stores one face value of `FVX/FVY/FVZ` per wall with the specified-normal-velocity rule
(`-RDXN*(HP(II+1)-HP(II))*DHFCT - DUUDT`), or zero for a mirror wall. `NOM` (from `EXTERNAL_WALL(IW)%NOM`) is used only to select the branch (`NOM/=0 .OR. SOLID .OR. NULL`) and to skip internal null walls; no
neighbour-mesh array is read. So the neighbour-mesh data do not enter, and the loop is fully translatable once the designator is a table column.

What is built (in the working tree of the wall-loop branch and in `(local workspace)` patches, not in the shared files of s5-gen):

* kernel `wall_no_flux`: the designator becomes a one-column table `EW_NOM`; `PRES_FLAG` is an argument; the branch on the six values of `IOR` and `PREDICTOR` is kept; bitwise tests
  against the verbatim loop (random meshes with all boundary types, `NOM`, `PRES_FLAG` and `PREDICTOR`);
* the race: a thin wall (zero-thickness obstruction) is seen from both sides, and both walls store into the same face, the later wall index winning. A threaded launch of the kernel
  is then a race. The sidecar says `unique = true` and the generator does not verify it. The answer is a host check of the wall tables (`s5_noflux_unique_check`) and a launcher that runs the
  kernel only when the check finds no face written twice, and otherwise leaves the walls to the host in upstream order; one small generator change records the requirement in the report, the
  header, the manifest and the golden signature (a sidecar field `unique_check`);
* effort remaining after the merge: the shared changes below, the regeneration, and the acceptance by the owner of the generator functions touched. About half a day.

Risks: the thin-wall race (handled by the check, not by the kernel); `DHFCT` depends on `PRES_FLAG` and on whether the wall is external, a table or a scalar argument must carry both; stale
`UN` cannot occur since each branch assigns it. The mirror branch zero store and the specified-velocity store can name the same face from different walls; the check counts any second store.

## 4. How the work lands on s5-gen after the merge

1. The wall-branch merge lands on s5-gen (done by the Legacy Mapper). It is not touched from here. Until then the patches stay in `(local workspace)`, unapplied.
2. After the merge: regenerate under the lock, re-run the fast checks, commit the patch for the PVF kernels (the three `pvf_*` entries and one hunk in the field-test generator) as one commit
   whose subject starts with `shared:`.
3. The L1366 changes to `s5gen.py` (two changes: the rule for the `EW_NOM` designator and the gather change) and the `unique_check` change are committed only after the owner of those
   functions (Legacy Mapper) accepts them, as separate `shared:` commits, each with the regenerated generated files and the golden signature.
4. My own files for L1366 (`s5_noflux_check.F90`, its two test files, the README) are new files and are committed after the merge by path; they need nothing from the owner.
5. Coordination: the Wall Loops Engineer owns the wall tables (`EW_NOM` is one more column of the same builder, the cell-to-wall list is theirs) and the face-write check of the wall-table
   spec (section 4 there) which `s5_noflux_unique_check` repeats for this loop; the Legacy Mapper owns the review of the `s5gen.py` hunks. Requests go through the coordinator; no file of the
   other two owners is edited from here.

## 5. Order of work

1. L1366 after the merge and the acceptance (about half a day).
2. Nothing for L1363, L1364, L1367, L1399 under D-055; revisit only if the ghost fill is dropped as default (then 4 to 6 days in total, buffer layout first).
