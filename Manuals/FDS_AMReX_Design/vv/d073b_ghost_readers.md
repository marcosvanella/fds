# D-073(b) / D-077(6): ghost readers of UII / UIID / QR

Owner: AMR V&V Lead (co-sign for Radiation Lead).
Pin: FireX `bee11f0329` (full `bee11f032907e53585dbb33ee76c30b0d2667a66`).
Script: `amrex/s4_mass/s5_gen/s5_rad_sweep/test/sweep_ghost_reads.py` (tree copies under `src-s5*`).
Docs cited: `docs/radiation/06-fr062-sweep-kernel-design.md` §6 (commit `dc44e30`); `vv-runs/gpu_gate/README.md` § "Rule D-073 (2)" (D-077(6) interior-only restart).

## Classification re-run

```
python3 .../sweep_ghost_reads.py --repo (repo root) --commit bee11f0329 --list
```

Result (saved `scratch/d073b/sweep_ghost_reads_bee11f0329.txt`):

- **77** Source lines mention `UII` / `UIID` / `QR`
- Classes: `N=32`, `W=14`, `I=11`, `G=1`, `V=9`, comments `-=10`
- **0 UNCLASSIFIED**; exit `GHOSTREAD OK`
- Nine `V` lines (all `dump.f90`): restart WRITE `3918/3920/3921`, READ `4107/4109/4110`; `CASE(18)` `8767`; `CASE(19)` `8769`; thermocouple `CASE(110)` `9184`

Matches the claim in radiation/06 §6.

## Verdicts

| # | Claim | Verdict |
|---|---|---|
| 1 | Physics: no reader of ghost `UII`/`UIID`/`QR` in the solve; only `IL` ghosts must stay bit-reproduced for radiation | **CONFIRM** |
| 2 | Restart `dump.f90` WRITE ~3918/3920/3921 and READ ~4107/4109/4110 I/O full arrays; compare stays interior-only; no other restart consumer needs ghost bit-identity for these three | **CONFIRM** |
| 3 | Output `CASE(18)`/`CASE(19)` do not clamp unlike `CASE(16)`; device/slice never address ghost **or** name clamp/refuse/keep-ghost-bits | **CONFIRM WITH CONDITION** |
| 4 | Thermocouple `dump.f90` ~9184: same as (3) | **CONFIRM** |

### (1) Physics — CONFIRM

From the classified reader list at the pin:

- Interior loops (`I`): `divg.f90` `QR(I,J,K)` in `1..*BAR`; `radi.f90` `RAD_Q_SUM` / `KFST4_GAS` interior; `dump.f90` `Q_DOT(3)` HRR sum interior.
- Wall gas cell (`G`): `radi.f90:5035` `UII(BC%IIG,BC%JJG,BC%KKG)` only.
- Whole-array (`W`): `UIIOLD`/`UIID`/`QR`/`UII=SUM(UIID)` stay inside `RADIATION_FVM`; downstream physics reads of those results are the interior/`G` sites above.
- No unclassified reader.

Radiation ghost fill that must stay bit-reproduced remains **`IL` only** (`WALL_LOOP1` / face exchange), as in radiation/06.

### (2) Restart — CONFIRM

Spot-check on pin (`git show bee11f0329:Source/dump.f90`):

- WRITE: `QR`, `UII`, conditional `UIID` at 3918/3920/3921 (complete allocatable arrays, ghost layer included; alloc is `0:IBP1` in `init.f90`).
- READ: same three at 4107/4109/4110 into the same full arrays.
- Gate rule D-077(6) / README "Rule D-073 (2)": restart-file compare of these three is **interior-only**; ghost-only difference → NOTE, not FAIL. Mask implementation waits on a restart compare path (none yet).
- In-tree restart consumer after READ is the FDS solve itself. Physics readers of `QR`/`UII`/`UIID` are interior or wall-gas (`I`/`G`); no other process in Source requires ghost bit-identity of these three. External post-processors of `.restart` are out of AMR V&V scope for this clause.

### (3) Output INTEGRATED INTENSITY / RADIATION LOSS — CONFIRM WITH CONDITION

Spot-check on pin:

- `CASE(16)` ABSORPTION COEFFICIENT clamps: `III=MAX(1,MIN(II,IBAR))` (and J/K).
- `CASE(18)` `8767`: `UII(II,JJ,KK)` — **no clamp**.
- `CASE(19)` `8769`: `QR(II,JJ,KK)` — **no clamp**.
- Both quantities default `CELL_POSITION=CELL_CENTER`, `SLCF_APPROPRIATE` and `ISOF_APPROPRIATE` true (`data.f90`).

**Device path — interior:**

- Point DEVC (`SPATIAL_STATISTIC='null'`): `read.f90` after `GET_IJK` clamps cell-center indices to `1..IBAR/JBAR/KBAR` (15376–15379). Runtime then does `I=MIN(IBP1,MAX(0,DV%I(1)))` (`dump.f90:7963`), which cannot pull a ghost once init clamped.
- Volume DEVC: `SDV%I1 = NINT(...)+1` with XB clipped to mesh; cell-center path does not allow zero as a cell index (`read.f90` ~14508).

**Slice / isosurface / smoke3d / plot3d node fill — can address ghosts:**

- Default (node-centered) SLCF: `SL%I1 = NINT(GINV(...)*RDXI)` may be `0`; dump expands `II1=I1`, `II2=I2+1` so `GAS_PHASE_OUTPUT` is called at `0` and `IBP1` when averaging to corners (`dump.f90` ~6920–6944).
- CELL_CENTERED SLCF spanning from the first interior cell can set `SL%I1 = I-1` → `0` (`read.f90` ~16242).
- ISOF fill loops `I,J,K = 0..IBP1` calling `GAS_PHASE_OUTPUT` (`dump.f90` ~4943–4944).
- SMOKE3D / related full-array fills likewise loop `0..IBP1` (`dump.f90` ~5078).

So the bare claim "device/slice never address the ghost layer" is **false for SLCF/ISOF (and any full-mesh node fill) of qty 18/19**. Devices alone would be CONFIRM.

**Named rule (R-D073b-out)** — pick one; V&V co-sign uses (A) as the default closing rule unless Radiation Lead elects (B) or (C):

| Option | Rule |
|---|---|
| **(A) keep-ghost-bits (default)** | If any active SLCF, ISOF, PLOT3D, or user SMOKE3D requests `INTEGRATED INTENSITY` or `RADIATION LOSS`, the box-split path must keep ghost `UII`/`QR` (and thus `UIID` feeding them) bit-reproduced for that run. Otherwise ghosts of `UII`/`UIID`/`QR` need not be bit-reproduced. |
| **(B) clamp** | Change `CASE(18)`/`CASE(19)` to clamp like `CASE(16)`. Upstream behaviour change at domain-boundary corners; needs FireX patch agreement. |
| **(C) refuse** | AMR mode refuses SLCF/ISOF/PLOT3D/SM3D of qty 18/19 until (A) or (B). |

Single-box kernel path already reproduces ghosts (sweep `GHOST` line); this rule binds the **box-split** path only when those outputs are present.

### (4) Thermocouple — CONFIRM

- `CASE(110)` `9184`: `UII(II,JJ,KK)` via `DV` (device state); same cell-center DEVC init clamp as (3).
- Not a slice/iso consumer of ghost `UII` in normal use (device indices stay in `1..*BAR`).
- Therefore thermocouple does **not** require ghost bit-identity of `UII`.

## Co-sign statement (D-073(b))

V&V Lead accepts radiation/06 §6 reader list at pin `bee11f0329` with the amendment that rows 7–8 are not unconditionally interior: DEVC is interior; SLCF/ISOF/full-mesh node fills can read ghosts. Closing clause for box-split:

1. Physics ghosts of `UII`/`UIID`/`QR`: **not required** bit-reproduced (`IL` only).
2. Restart ghosts: **interior-only compare** (D-077(6)); no other restart consumer needs ghost bit-identity.
3. Output: **R-D073b-out (A)** unless Radiation Lead chooses (B) or (C).
4. Thermocouple: **no ghost requirement**.

With (A)–(C) named, D-073(b) / D-077(6) / A-50 reader confirmation is ready for Radiation Lead co-sign.

## Evidence index

| Item | Location |
|---|---|
| Classification stdout | `scratch/d073b/sweep_ghost_reads_bee11f0329.txt` |
| Script | `src-s5wall-dev/.../s5_rad_sweep/test/sweep_ghost_reads.py` (and siblings) |
| Pin tree | this repository @ `bee11f0329` |
| Restart / output / TC lines | `Source/dump.f90` as cited |
| DEVC clamp | `Source/read.f90` 15374–15379 |
| SLCF index | `Source/read.f90` 16229–16253; dump expand ~6920 |
| ISOF fill | `Source/dump.f90` ~4943 |
| Gate restart rule | `vv-runs/gpu_gate/README.md` Rule D-073 (2) |
| Design draft | `docs/radiation/06` §6 (`dc44e30`) |
