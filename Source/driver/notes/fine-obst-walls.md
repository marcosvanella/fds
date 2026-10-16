# Fine-grid-conforming OBST wall cells on fine AMR boxes

Design note. Co-authored: Role 1 (lead: fine wall build, fluid coupling, composite sums) and Solid Phase Lead (wall state, upstream fit, review).

## Role 1 sections

Base: HEAD with patch 0026 (`20f8a91b62`), hidden fine walls (`27ba12cdbe`, `d5575b2c70`), and Role 2's burner work (`ffdcfa8203`,
`4e30615ef7`, `c27dcd7ed9`; `ce1eb49619` on `species-stage2`). Wall state follows S1; the upstream fit follows S2. This part covers the wall
inventory (build), the faces that have no parent, fluid-side coupling, composite sums, tests and effort.

### R1. What is already fine, what is not

- **Fine solid cells are already fine-snapped:** `SideData` comp 0 comes from `pb::obst_list()` (ObstMask.H, FDS `NINT` at the level's dx on
  the input XB exported by `fds_obst_xb` = `OB%X1..Z2`), and `FDS_FINE_SET_SOLID` writes it into the fine `CELL_INDEX`.
- **The wall inventory is still parent-derived:** `FDS_FINE_WALLS_BUILD` puts a wall only on a fine face that lies on a level-0 SOLID/OPEN wall
  face and takes the parent's SURF/VENT. A fine OBST face that does not coincide with a coarse wall face (sub-coarse OBST, a fine step, a VENT
  edge that snaps differently) gets no wall (S3 findings 5b, 6). Role 2's `fds_fine_walls_state`/`_un` select by the same parent rule.

### R2. Critique of the proposed direction

1. *Run FDS's OBST/VENT setup on the fine box.* Goal yes, means no. FDS's setup spans `READ_OBST`/`PROC_OBST`/`READ_VENT` (per-mesh snapping
   into `MESHES(NM)%OBSTRUCTION/VENTS`), `INITIALIZE_MESH_VARIABLES_1/2` and `INIT_WALL_CELL` (`M=>MESHES(NM)`, `M%VENTS`, `M%OBSTRUCTION`,
   `EXTERNAL_WALL`, `SEARCH_OTHER_MESHES`). Re-running it needs `MESH_OBJ` in a dozen read/init routines plus per-box OBST/VENT tables, against
   D-084 and S2 ("do not add more hunks"). Instead: one driver builder doing the same arithmetic on the global input lists (below).
2. *Coarse keeps a snapped OBST only under covered cells.* Rejected. Level 0 stays FDS's own (L0 bit-identity; `aeb2ebf6de` showed that masking
   L0 loses the MASS_FLUX burner velocity). Covered L0 cells are overwritten by average-down. The real risk is an interface that crosses an
   OBST whose coarse and fine snaps differ; handle it with tagging (3).
3. *Tagging keeps sub-coarse OBSTs fine.* Yes, with the rule "tag every OBST/VENT whose coarse snap differs from its fine snap, grown by NBUF
   coarse cells", `--tag-obst-snap NBUF` (owner 2026-10-15: default **on** for two or more levels; single-level baselines stay off).
4. *Fine wall state at regrid by injection from neighbour/parent.* Superseded by S1: nothing is carried on the fine level; the few history
   quantities are read from a level-0 source wall each step. R3 defines that source for faces without a parent. No injection path exists
   that a thick wall could fall into (thick/ONE_D stays an error under a moving box).

### R3. Fine OBST faces with no parent face

Build rule (replaces the parent-face test): a fine gas cell gets a wall on face D if the neighbour across D is a **fine solid cell** or the
**domain boundary**. Coarse-fine faces never get one. Each wall gets a classification:

| case | example | PIW | SURF / VENT | history source (S1 rows T_IGN, M_DOT adjust, Q_RAD_IN) |
|---|---|---|---|---|
| P: on a parent wall face | coarse-aligned OBST, domain wall | L0 wall | parent's | the parent (today) |
| N: no parent, OBST face | sub-coarse OBST, fine step, OBST face inside a coarse gas cell | 0 | `OB%SURF_INDEX(IOR)` of the owning OBST (input list); VENT override if a fine-snapped VENT plane covers it (INIT_WALL_CELL VENT_SEARCH rule: plane, ALLOW_VENT, RADIUS) | **source wall** = the L0 wall of the same OBST and IOR nearest to the face centre (ties: lowest L0 wall index), if one exists; else none |
| D: no parent, domain face | VENT on the domain boundary that snaps differently | 0 | DEFAULT or the fine-snapped VENT | nearest L0 wall on the same domain face with the same SURF |

Rules for case N/D without any source wall (an OBST that vanished on L0, e.g. thinner than a coarse cell and not THIN):
- inert / fixed-T: no history needed (S1 column 1 is a pure function of SURF, T and gas). Supported.
- HRRPUA/MLRPUA burner: T_IGN = `SF%T_IGN` (input), M_DOT_G_PP from SURF ramp, `M_DOT_G_PP_ADJUST` = 1 (no parent total to conserve). The
  L0 mass source is then missing for that burner (L0 has no face): the composite mass gate must count the fine burner inflow; the coarse
  level does not see it except through average-down. Flag it in the log; allowed only with `--tag-obst-snap` so the burner is always fine.
- Q_RAD_IN: the source wall's, else `EMISSIVITY*SIGMA*TMPA4` (the INITIALIZE_BOUNDARY_PROP1 value). Exact once fine radiation exists.
- NET_FLUX / ADIABATIC with no source: error stop (needs Q_RAD_IN that we cannot supply).

Composite sums (S3 finding 6): for case P the parent share stays `FR*parent`. A covered parent face that is partly behind fine solid cells
gets `FR` from the fine gas faces only, so the solid-covered part is no longer counted; the N faces add their own fine terms with no parent
share. That closes findings 5b and 6.

Fluid side: Role 2's accessors keep their C signatures (`pressure_backend/frozen/fine-wall-state.md`). The selection becomes "fine wall whose
SURF has a mass flux" (lists the same faces for coarse-aligned burners, plus case N/D). `s[2]` (area weight to the parent) is 0 for PIW=0;
consumers must not divide by it. The values switch from parent B1 to the fine wall's B1 (built per S1) behind one switch in
`fds_fine_level.f90`, after `ce1eb49619` is merged.

### R4. Files and estimated lines

| file | change | lines |
|---|---|---|
| `driver/fds_mesh_query.f90` | export OBST `SURF_INDEX(-3:3)`, ALLOW_VENT, THIN; VENT list (XB, IOR, SURF, RADIUS, X0..) | +60 |
| `driver/fds_fine_level.f90` | builder v2 (solid-neighbour scan, fine NINT, VENT match, P/N/D, source wall), refuse rules | +150 / -40 |
| `driver/fds_fine_level.f90` | state for PIW=0 (source wall or SURF input), Role 2 accessor source switch | +60 / -20 |
| `driver/fds_fine_level.f90`, `TimeLoop.cpp` | sums: FR from fine gas faces, N faces without parent share | +10 |
| `regrid_transport/RegridAmrCore.*`, `pressure_backend/ObstMask.H` | `--tag-obst-snap`, `snap_differs()` | +50 |
| tests | `run_fine_obst_walls.sh`, two inputs | +80 |
| upstream | none (S2: no new hunks) | 0 |

### R5. L0 bit-identity, D-084, rebase

- Level 0 untouched; `MESH_OBJ` null outside the fine routines; exports read-only. Proof as for 0026: single-level `--run` byte compare, and
  OMP=1 two-level TWO-LEVEL lines/devc/steps against the previous HEAD (OMP=2 runs are not reproducible run to run beyond ~step 90).
- D-084: one face rule, one snap rule (shared with the masks), no copy of READ_OBST/INIT_WALL_CELL, no new upstream code.
- Rebase: all wall access stays in three routines (build, state, sums) via `B1_INDEX/B2_INDEX` and named B1 fields; S2's fallback (allocate
  BOUNDARY_* directly if ALLOCATE_STORAGE changes) fits the builder unchanged.

### R6. Test plan

0. **Owner's primary acceptance test (run first).**
   - **Setup:** an isothermal case (no fire, no heat sources, `TMPA` everywhere) with one or more closed compartments, separated by thin OBST
     walls thinner than a coarse cell (e.g. `XB=1.80,1.83,...` at coarse dx 0.1 m). Each compartment is its own `&ZONE`.
   - **Hierarchy:** driven through refinement and derefinement/re-snap with forced tags that move every few steps (`--regrid-every 2`), so
     each wall is alternately fine (solid cells or a fine thin plane) and coarse (a thin plane). Variant: the wall is shifted so that its coarse
     plane moves at re-snap. Run with 2 and 3 levels.
   - **Pass:**
     - T stays isothermal: uniform and unchanged to round-off on every level after every regrid.
     - Z unchanged to round-off.
     - Background pressure PBAR of every zone unchanged to round-off.
     - Zones stay separate: zone IDs match the L0 IDs cell by cell, and no flow or mass crosses a wall.
   - **Diagnostic (not a pass criterion):** the total and per-zone mass change at each regrid (`REGRID ZONE-MASS`, R10 item 4), expected to be
     the re-snapped volume times rho.
   - L0-only run bit-identical.
1. **simple_test burner fine** (`--tag-obst 1 --nofine -1 --patch 0 7 0 7`): coarse-aligned burner, so every wall is case P; inventory = L0
   inventory refined (4 fine faces per coarse face, same SURF); Role 2's `run_fine_burner_check.sh` numbers unchanged; L0 bit-identical.
2. **simple_test, burner shifted by half a fine cell**: case N faces exist; burner mass flux on fine faces = HRRPUA*A/dH within 1e-12 relative;
   composite mass gate PASS with the fine burner counted; `--tag-obst-snap` keeps the interface off the OBST.
3. **pyramid with fine-conforming steps** (`pyramid_48_4m.fds` variant, steps of half a coarse cell): per-IOR fine wall counts equal a
   single-mesh FDS run at fine dx over the refined region (exact); Q_CON_F sum within resolution of FDS fine; closed-domain mass PASS.
4. Regrid (`--regrid-every 5`): inventory rebuilt, nothing carried, L0 unchanged.

### R7. Effort

| step | content | estimate |
|---|---|---|
| A | exports, builder v2 (P/N/D), snap tagging, tests 1 and 3 (counts) | 1.5 days |
| B | PIW=0 state, accessor switch, composite FR fix, tests 2 and 4 | 1.5 days |
| review | Solid Lead questions, Role 2 merge | 0.5 day |
| **total** |  | **~3.5 days** |
| later | walls live in the fine kernels (WALL_BC, velocity BC, radiation) | 2 weeks or more, separate note |

### R9. Owner rule: every OBST on every level (thin OBSTs)

Rule (owner): an OBST resolved on a fine level but thinner than a cell of a coarser level must not vanish there; on that level it becomes a
zero-thickness (thin) obstruction, as FDS does. Level 0 stays FDS's own.

**What FDS already does on L0 (read.f90, READ_OBST):**
- Snap: `OB%I1 = NINT(GINV(XB1-XS)*RDXI)`, `OB%I2 = NINT(GINV(XB2-XS)*RDXI)` (read.f90:11157-11162), same for J, K.
- Without THICKEN, an OBST thinner than `0.25*DX` that still snapped to two different faces is collapsed onto the nearer face, `I1=I2`
  (read.f90:11175-11181, and :11192, :11209 for J, K).
- Collapsed in **one** direction: kept as a zero-thickness OBST; `OB%THIN = .TRUE.` when the input length is below `0.5*DX`
  (read.f90:11315-11317). Its two faces become wall cells on both gas sides of that face plane.
- Collapsed in **two or three** directions (a rod or a point at that resolution): thrown out (read.f90:11220-11225, "Throw out obstructions
  that are too small"). Such an OBST does vanish on L0. That conflicts with the owner rule, and we cannot change it without touching L0
  (owner decision D1 below).
- An OBST thinner than a cell that straddles a face can snap to a full cell (`I1/=I2`, thickness `>=0.25*DX`): FDS thickens it. Same rule at
  every level.
- `THICKEN=.TRUE.` (read.f90:11167-11172) makes a collapsed OBST one cell thick instead.

So on L0 FDS keeps thin OBSTs as zero-thickness walls in the one-direction case, and drops them only in the rod/point case.

**What the driver does today on levels >= 1 (a gap):** `pb::obst_cells` (ObstMask.H) returns an empty box when any direction snaps to zero
thickness, so a thin OBST marks no cell on a fine or intermediate level. It is also absent from the fine face mask and from the hidden walls.
A sheet OBST that is thin on L1 is lost on L1 today.

**Proposed per-level rule (each level uses FDS's rule at its own dx):** for level l, snap the input XB with FDS NINT at dx_l, apply the
`0.25*dx` collapse rule:
- 0 collapsed directions: solid cells (today).
- 1 collapsed direction: a **thin OBST**. No solid cells. In `SideData`, both faces at that plane get the wall code on the two gas cells
  (comps 1..6), so the pressure/transport masks block the face. The wall builder creates two walls (IOR +-d), one per gas side, with the
  OBST's `SURF_INDEX(+-d)`, exactly as FDS does for THIN.
- 2-3 collapsed directions: dropped, as on L0 (D1).
- THICKEN: one cell thick, as FDS.

**3-level runs:** an OBST can be solid on L2, thin on L1, and thin or gone on L0. Each level carries its own representation. The covered L1
cells under L2 are overwritten by average-down, so the L1 thin wall matters only where L1 is uncovered, and at the L1-L2 interface.
- **Interface rule:** the interface must not cut an OBST whose representation differs between the two levels (solid vs thin, or different
  snapped extent). Otherwise the two levels disagree on which faces are blocked.
- **Tagging:** any OBST whose snap result (class or extent) differs between level l and level l+1 is a `--tag-obst-snap` candidate at level l,
  grown by NBUF cells of level l. That includes every OBST that is thin on a coarse level and solid on a finer one, and every OBST dropped on
  L0 but present on a finer level.
- **Composite sums:** a thin L0 OBST under fine cells is a parent wall pair; its fine solid faces are case P only where the planes coincide,
  otherwise case N (R3).

**Code impact** (adds to R4): `ObstMask.H` per-level class (solid / thin / dropped) and face-mask marking for thin planes, about 40 lines.
The builder's thin-plane walls are about 30 lines. `snap_differs()` compares class and extent, about 10 lines. Effort: +0.5 day (total about
4 days).

**Test:** a sheet OBST of 1/4 coarse cell in simple_test (e.g. a baffle `XB=1.80,1.83,...`), three levels. On every level the face count of
the thin plane equals a single-mesh FDS run at that level's dx (thin walls both sides). L0 is bit-identical. Mass PASS.

### R10. Pressure zones and volume

Owner rules: (a) closed compartments (pressure zones) stay separate on every level, and a thin or zero-thickness coarse OBST still seals
them. (b) Re-snapping an OBST to a different coarse grid plane changes the room volume. Whatever we do must not change T or Z.

**How FDS assigns and seals zones (L0, unchanged):**
- `INITIALIZE_PRESSURE_ZONES` (main.f90:2713) calls `ASSIGN_PRESSURE_ZONE` for each `&ZONE` XYZ (main.f90:2742) and again for the mesh
  interfaces (main.f90:2787, 2988, 3011).
- `ASSIGN_PRESSURE_ZONE` (func.f90:5661-5832) is a flood fill over gas cells:
  - It does not march into non-removable solid cells (func.f90:5721-5722).
  - It **does not cross a thin OBST face**: for each direction it checks every OBST with `I1==I2` (or J, K) on that face and cycles
    (func.f90:5763-5778). A zero-thickness OBST therefore seals two zones on L0, provided it survived the snap (R9: it collapses in one
    direction only).
  - Overlap of two ZONEs is a setup error (main.f90:2743-2745).
- So on L0, FDS's own rounding already gives sealed zones with thin walls. The volume of each zone is the L0 gas-cell volume. Re-snapping a
  wall to a different coarse plane changes it by up to one coarse cell layer of wall area.

**Today in the driver (a gap):** `FDS_FINE_LEVEL_CREATE` refuses a rank whose level-0 meshes hold more than one pressure zone (IERR 3,
fds_fine_level.f90) and gives every fine box the single zone. Multi-zone cases do not run with fine levels yet.

**Proposed rules:**
1. **Zone IDs come from L0 and are inherited, never flood-filled again on a fine level.** A fine gas cell takes the zone of the covering L0
   gas cell. Where the fine snap differs (a fine gas cell under an L0 solid cell, or on the other side of an L0 thin plane), it takes the
   zone of the L0 gas cell it is connected to through fine gas faces that do not cross a fine solid or thin face. This is a flood fill
   restricted to the fine box, seeded from L0, run at build/regrid. A fine cell that connects to two different L0 zones is a build error:
   the fine walls would open two compartments to each other.
   - Same IDs on every level and at every regrid. No new zone numbers. The thin-OBST seal applies on every level through R9's per-level
     thin faces (blocked face mask + two walls).
   - Intermediate levels (3-level runs) inherit from their parent level the same way.
2. **Zone volume is defined from the finest covering level.** `V_z` = the sum of uncovered gas-cell volumes over all levels, each cell
   counted at its finest covering level (composite volume, the same mask as the composite energy sums in TimeLoop `energy_sums`). This is
   used for the zone's `D_PBAR_DT` (the composite pressure-zone integrals). On L0-only runs it is FDS's volume (bit-identical).
3. **Coarse covered cells are bookkeeping only:** filled by average-down of the fine conserved fields (rho, rho*Z, rho*h_s). An L0 cell that
   is solid on L0 but gas on the fine level gets the fine average and is not part of L0's FDS update. An L0 gas cell that is partly solid on
   the fine level gets the gas-volume-weighted fine average. Uncovered L0 stays FDS's own.
4. **Derefinement / OBST re-snap (owner correction, supersedes conservation):** preserve the **intensive** state: Z (species mass
   fractions), T, background pressure PBAR (per zone) and the divergence D. Do **not** enforce zone mass or energy conservation here.
   - A coarse cell that becomes gas (it was solid on the coarse level, gas on the fine level) takes T and Z as the gas-volume-weighted average of
     the fine gas cells it covered. rho then follows from the equation of state at the zone's unchanged PBAR: `rho = PBAR/(R(Z) T)`. D is the
     covered fine average (average-down, as for any covered cell).
   - Fine gas that falls inside a coarse solid cell simply disappears with its volume. Neighbouring cells are not changed.
   - Refinement (the mirror case): fine cells under a coarse solid that become gas get T and Z from the adjacent same-zone gas cells, and rho
     from the EOS at PBAR.
   - Result: no jump in T, Z, PBAR or D. The total and per-zone mass (and energy) change by the volume difference times rho. That change is
     accepted and **reported as a diagnostic**: a `REGRID ZONE-MASS` log line per regrid and zone with dm and dm/m. It is excluded from the
     mass gate the same way as the boundary inflow (the gate subtracts it).
   - The gate checks that max|dT|, max|dZ| and max|dD| over unchanged cells are exactly 0 across a regrid, and that PBAR is unchanged.
5. Thin coarse OBSTs never open: a coarse thin plane under a fine box is also a fine thin plane or fine solid cells at the same position
   (or the interface is kept off it by `--tag-obst-snap`, R9). So a zone can never leak through a level change.

**Code impact** (adds to R4): zone inheritance flood fill in the builder (~60 lines, fds_fine_level.f90), lifting IERR 3 (~10), composite
zone volume and zone integrals in the pressure-zone path (~40, TimeLoop.cpp/CompositeOps), regrid fill for status-changing cells at fixed T/Z/PBAR (~60,
regrid_transport CellTransfer), per-zone gate lines (~30). About +1.5 days (total about 5.5 days).

**Test:** a closed compartment whose dividing wall is thinner than a coarse cell. In simple_test geometry, a `&ZONE` on each side of a wall
`XB=1.80,1.83,0,2.4,0,2.4`, with no leak and a burner on one side, 2 and 3 levels, `--regrid-every 5`.
- L0: the wall is a thin OBST (read.f90:11315), and FDS gives two zones (func.f90:5763-5778).
- Each level: zone IDs match L0 cell by cell (comparing the covering L0 cell).
- Across every regrid: T, Z and D are unchanged in cells whose status did not change, and PBAR is unchanged. The per-zone mass change is
  reported (diagnostic) and equals the re-snapped volume times rho.
- The zone without the burner has `D_PBAR_DT = 0` to round-off (no leak through the wall).
- Variant: shift the wall by half a coarse cell, so the coarse plane changes. Same checks, and the zone volume reported from the finest level
  matches FDS fine within one fine cell layer.
- L0-only run bit-identical.

### R8. Questions for the Solid Phase Lead and the owner

1. [SP] Case N source wall "same OBST and IOR, nearest L0 wall": acceptable, or wait for the FR-046 face key?
2. [SP] Case N burner with no L0 face: `M_DOT_G_PP_ADJUST = 1` and input T_IGN. Acceptable, given the L0 mass source is then absent?
3. [SP] D-010 keeps stateful VENTs on L0 (no-refine buffer) while Role 2 runs the burner on L1. Confirm this is within the static-hierarchy
   clarification.
4. Owner: level 0 never re-snapped (R2 point 2); `--tag-obst-snap` as the sanctioned interface rule.
5. Owner: live fine walls are a separate decision.

---

## S1. Wall-state initialization at build/regrid (Solid Phase Lead)

Rule: a fine wall owns only what is a pure function of (surface, geometry, current gas). Everything with time history stays on the level-0 parent wall and is read, never transferred.

| Quantity | Inert / fixed-T (SPECIFIED_TEMPERATURE) | NET_FLUX / ADIABATIC | Thermally thick, burner with ONE_D |
|---|---|---|---|
| TMP_F | recomputed each call (wall.f90:737-747: TMP_0 or TMP_FRONT ramp; TMP_G if U_NORMAL>0) | recomputed from the face balance (wall.f90:781-827, uses Q_RAD_IN); parent value is only the iteration start | stays on L0 (D-010) |
| TMP_B, ONE_D (TMP, layer thickness, mass, char) | none | none | L0 only; fine box must not carry the face (FR-041, D-010) |
| T_IGN, ramp clocks (TIME_TEMP, TIME_HEAT) | copy from parent (needed by ramps, wall.f90:739,764) | copy | L0 |
| HRRPUA / MLRPUA burner VENT (no ONE_D) | n/a | n/a | M_DOT_G_PP_* is a function of T, T_IGN, ramp and area: recompute from parent T_IGN; adjust factor (M_DOT_G_PP_ADJUST) copied from parent so the parent face total is conserved |
| U_NORMAL, ZZ_F, RHO_F | from the parent (burner/leak) or the fine balance | same | L0 |
| Q_RAD_IN | parent value (no fine radiation) | parent value | L0 |
| HTC, Q_CON_F, U_TANG, TMP_G, RHO_G, ZZ_G | recomputed from fine gas | recomputed | L0 |
| B2: U_TAU, Y_PLUS | must be recomputed at the fine DN, not copied (see R3) | same | L0 |

Carried across regrid: nothing on the fine level. The fine wall set is rebuilt at every bind/regrid from the parent wall index plus the surface; the state it needs (T_IGN, M_DOT adjust, Q_RAD_IN) is read from the parent each step. This keeps D-010 trivially satisfied for thin walls and needs no prolongation/restriction of wall data. Thermally thick or ONE_D-stateful faces under a moving fine box are an error (current behaviour); for a static hierarchy D-010 allows them to stay coarse-owned.

Simplicity: one copy of B1 from parent per wall per step (already done), no wall arrays saved across regrids, no per-wall flags beyond the parent index (PNM, PIW).

## S2. Fit with the upcoming upstream WALL/ONE_D/BOUNDARY layout change

No wall refactor on our side. Patch 0026 touch points and their risk:

1. mesh.f90 `MESH_OBJ` pointer: low risk (new line after MESHES declaration).
2. func.f90 HEAT_TRANSFER_COEFFICIENT `M =>` (func.f90:3156): medium. Upstream may pass B1/B2 directly instead of an index; then the hunk disappears (better for us).
3. func.f90 ALLOCATE_STORAGE `M =>` (func.f90:3968): high. This is where WALL/BOUNDARY_PROP1/ONE_D storage is laid out; any upstream split of ONE_D will rewrite it.
4. func.f90 INITIALIZE_BOUNDARY_PROP1 `B1 =>` (func.f90:4981): high, same reason.
5. wall.f90 `PUBLIC SURFACE_HEAT_TRANSFER`: low.

Keep the surface minimal: the three hunks are one identical selector line each; do not add more. If upstream changes ALLOCATE_STORAGE's signature, prefer to drop hunks 3-4 and allocate the hidden walls' BOUNDARY_COORD/PROP1/PROP2 directly (they need no ONE_D for the supported surfaces), which removes our dependency on the storage layout entirely. Do not depend on ONE_D field names in driver code.

## S3. Review of Role 1 milestones 1-3

Commits: 20f8a91b62 (patch 0026), 27ba12cdbe (hidden fine walls, fds_fine_level.f90), d5575b2c70 (composite Q_CONV/Q_COND/Q_DIFF, TimeLoop.cpp).

**Verdict: Accept-with-notes.**

Findings:

1. Q_CONV/Q_COND/Q_DIFF definitions and signs match UPDATE_HRR (dump.f90:11738-11753): `-U_N rho_f h_s A`, `-Q_CON_F A`, `-2 sum h_s,a rho D_f (Z_g-Z_f) RDN A`, AREA_F = B1%AREA. U_N selection per IOR matches dump.f90:11724-11735. OK.
2. HTC distance: with MESH_OBJ associated and WALL_INDEX_IN=IW, HEAT_TRANSFER_COEFFICIENT reads the hidden wall's B1, and `DN = 1/B1%RDN` (func.f90:3185) uses the fine RDN set in FDS_FINE_WALLS_SUMS. Default model then includes the resolved limit `2 CONV_LENGTH/DN` (func.f90:3253,3264) at the fine spacing: correct. Correction note: the review question "is fine DELTA used" — yes for DEFAULT/natural/forced.
3. LOGLAW (and other B2-based models, func.f90:3266-3271): B2 is copied from the parent (U_TAU, Y_PLUS at the coarse first-cell distance). With LOGLAW_HTC_MODEL on a fine box this gives the coarse y+. Either recompute U_TAU/Y_PLUS from fine U_TANG and DN (one WALL_MODEL call) or error out for LOGLAW like VEL_T. Recommend the error (fewer lines) until needed.
4. Inert walls: FDS INERT is SPECIFIED_TEMPERATURE with TMP_FRONT = ambient; SURFACE_HEAT_TRANSFER overwrites TMP_F (wall.f90:737-747), so copying the parent TMP_F is harmless and consistent. Adiabatic / NET_FLUX walls recompute TMP_F from the fine balance with the parent Q_RAD_IN (wall.f90:781-827): consistent given no fine radiation; the copied TMP_F is only the iteration start. The doc string "fine walls copy parent TMP_F" should say "as initial value; recomputed for all supported thermal BCs" (except NO_CONVECTION, set to TMP_G, and OPEN, set by inflow/outflow).
5. Missing fine radiation does not explain the Q_COND gap directly: for SPECIFIED_TEMPERATURE walls Q_CON_F = HTC*(TMP_G-TMP_F) has no radiation term (wall.f90:760). Radiation enters only indirectly via the fine gas TMP_G (fine cells exchange no radiative source) and via TMP_F of net-flux walls. More likely causes of -7.99 vs -11.69 kW: (a) only the covered fraction of the wall area is fine (2-12% share), so the composite can at most move that share toward fine FDS; (b) fine-conforming OBST faces with no parent face get no wall (face test, FDS_FINE_WALLS_BUILD) and are missing entirely; (c) 5 s window transients. Suggested check: report the fine-covered fraction of the total wall area and the Q_COND per unit area on fine vs FDS fine on the same faces.
6. Area consistency: composite = L0 - sum FR*parent + fine, FR = A_fine/A_parent. No double count where the parent face is covered by fine gas faces. Where a covered parent face is partly faced by fine SOLID cells (fine-conforming OBST), the parent share for that part stays counted although that area is now inside a fine obstacle, and the fine obstacle's own faces are absent (finding 5b). Net error is of order the non-conforming area; acceptable for milestone 3, to be closed by the Role 1 fine-OBST-face section.
7. Identity check (Q(9)) runs a second SURFACE_HEAT_TRANSFER per wall every step, while it is printed every 50 steps. Simplicity: gate it on the print cadence (pass a flag) to halve the wall cost.
8. Identity check uses BC0 for IOR/indices but HTC looks up `M%BOUNDARY_COORD(WC%BC_INDEX)` of the hidden wall (MESH_OBJ still associated): same IOR, so it passes, but it is not a pure parent evaluation. Fine as a path check; say so in the comment.
9. K_G and RHO_D_F for OPEN and Q_DIFF are parent values (coarse film state). Small error; acceptable, note it.
10. Parent mesh search in FDS_FINE_WALLS_BUILD is O(cells x meshes); acceptable at build/regrid only.
11. Patch 0026 is L0 bit-identical when MESH_OBJ is unassociated; the risk is upstream drift (S2), not correctness.

## Owner decisions

- **D1** OBSTs that collapse in two or three directions on L0 (rods/points) are dropped by FDS (read.f90:11220-11225). Keep L0 as FDS's own
  (they stay absent on L0 and are always tagged to the finest level), or treat this as an accepted exception to "every OBST on every level"?
- **D2** Level 0 is never re-snapped or masked by the driver (R2 point 2).
- **D3** Levels >= 1 follow FDS's snap rule at their own dx, including thin (zero-thickness) OBSTs as blocked faces with two walls (R9).
- **D4** `--tag-obst-snap NBUF`: the coarse-fine interface must not cut an OBST whose representation differs between the two levels. Opt-in
  (default off) or default on?
- **D5** Faces with no parent face (R3 case N/D): history source = the nearest L0 wall of the same OBST and IOR, else SURF input. Needs the
  Solid Lead's sign-off (R8 items 1 and 2).
- **D6** Live fine walls (in the fine kernels) are a separate decision, not part of this note.
- **D7** Pressure-zone IDs come from L0 and are inherited by every level (seeded flood fill inside fine boxes, error if a fine cell joins two
  L0 zones), and are never renumbered at regrid (R10).
- **D8** Zone volume and zone integrals are composite: each gas cell is counted at its finest covering level. Coarse covered cells are
  average-down bookkeeping only.
- **D9** At derefinement or OBST re-snap, preserve the intensive state (Z, T, PBAR, D) and do not enforce zone mass/energy
  conservation. The mass change from the volume difference is accepted and reported as a diagnostic per regrid and zone, and excluded
  from the mass gate (owner correction).

---

# N-level extension (Role 1 design, 2026-10-15)

Status: M0/M1 DONE (2026-10-16, see R20); the text below is the design, corrected where the code showed otherwise. Owner lifted the hold on
fine OBST walls and asked for N levels, not a level-1 patch. The level>=2 segfault in `fds_fine_walls_build` was diagnosed: see R20.1. Section numbers continue R1-R10 above; this part is R11-R19.
Solid Phase side (R16): written from S1/S2 and this note only. The Solid Lead has not been consulted; items marked [SP] need sign-off.

## R11. Where the code is level-1 only today (the inventory of special cases)

Read from `fds_fine_level.f90` and `TimeLoop.cpp` at HEAD:

1. `FDS_FINE_WALLS_BUILD(LEV,IB)` took every fine wall from a level-0 parent wall (found by scanning `MESHES`). **Correction:** an earlier version of
   this item blamed "no L0 wall for a face interior to a refined L1 region" for the level-2 crash. That was wrong, see R20.1: the real cause is
   that the parent L0 mesh belongs to another MPI rank and its `WALL` array is not stored there. Replaced by the per-level builder (R20.2).
2. `FDS_FINE_WALLS_SUMS(LEV,...)` takes `FR = A_fine/A_parent` against an L0 parent wall `PNM/PIW`.
3. `TimeLoop::composite_energy` calls `fds_fine_walls_sums(1, ...)` only when `lc.size()==2`; with 3+ levels it prints a warning and the
   three columns stay level-0 sums (not leaf).
4. `FINE_DEVC`, `LEAF_WALL_COUNT/LEAF_WALL/LEAF_WALL_COVERED` hooks are bound from level 1 builds; they iterate `FWALL(:)` (already
   level-indexed: `FWALL(LEV)%B(IB)`, good) but the covered test is "under a level-1 box".
5. `fds_fine_walls_state/_un` (burner faces, Role 2 accessors) take `lev` but select by the L0 parent rule; the callers pass `l` already.
6. The OBST snap (`ObstMask.H`, R9) and `SideData` are already per-level; the thin-plane marking (R9) is still to be implemented.
7. Create/destroy ordering: `FDS_FINE_LEVEL_CREATE` only lets the top level be (re)made (`ANY(N_BOXES(LEV+1:)/=0)` -> IERR 1), so a level's
   walls can be rebuilt only when all finer levels are gone. Compatible with a level-looped build, but see R14 (regrid order).
8. IERR 3: multi-zone L0 rank refused (R10).

None of these is design-essential; they all come from "parent = L0". The fix is one idea: **the parent of a level-l wall is whatever
level l-1 holds at that face, and the source wall is found by walking down the chain.**

## R12. One level-looped wall path (D-N1)

**Rule.** For every level l >= 1 and every box, the builder runs the same code. There is no `if (lev==1)`. Level 0 keeps FDS's own walls.

Inputs per level l (all already per-level or exportable): the box's fine `CELL_INDEX` solid flags (from `SideData`, snapped by FDS NINT at dx_l,
R1), the per-level thin-plane list (R9), the global OBST/VENT lists (R4 exports), the domain boundary.

Per fine gas cell and direction D (the R3 rule, unchanged): a wall if the neighbour is a level-l solid cell, a thin plane on that face, or the
domain boundary; never on a coarse-fine face. Classification P/N/D as R3, with "parent" redefined:

- **Parent wall** of a level-l face = the level-(l-1) wall on the covering face, if any. It is found by the *geometric covering chain*: face
  centre -> covering cell on level l-1 (box containing it, searched in `FINE_LEVEL(l-1)%BOX`, or in `MESHES` when l-1 = 0) -> the wall of the
  same IOR on that cell. The parent level-(l-1) wall may itself be hidden (l-1 >= 1) or real (l-1 = 0). Storage: `PNM` = mesh index
  (global numbering already used by `NM0+IB`), `PIW` = wall index in that mesh, `PLEV` = parent level (new, 1 integer).
- **Source wall** (history, R3 case N/D): walk up the chain `l-1, l-2, ... 0` until a wall of the same OBST/IOR (or same domain face and SURF)
  exists; take the nearest one at the *lowest* level that has it (history lives only on L0 by S1, so the walk ends at L0; a chain without an
  L0 wall is "no source" as in R3). `SRCNM/SRCIW` is stored once at build, so the per-step read is O(1) and level independent.
- **Burner VENT stateful on L0** (D-010): the history is read from the L0 source wall at any level. No wall state is stored at any level
  above 0 (S1 stands unchanged for N levels).

**Search cost.** Replace the O(cells x meshes) scan by a per-level box lookup (AMReX BoxArray on level l-1 via a small C++ callback, or a
Fortran loop over `FINE_LEVEL(l-1)%BOX`, then `MESHES` for l-1=0). Build/regrid only. Fortran loop first (about 15 lines).

**Where the code goes.** Same three routines (build, state, sums); `FWALL(LEV)%B(IB)` already generalizes; the level-1 specifics are removed:
the `MESHES` scan becomes `FWALL_PARENT(LEV,...)`, a helper that takes level, point, IOR.

## R13. Leaf-only sums for N levels (D-N2)

Current: composite = L0 - sum FR*parent + fine, FR = A_fine/A_parent, two levels. Generalize with a **leaf** definition:

- A wall is *leaf* if no finer-level fluid face covers any part of it. L0 real walls: leaf when not covered by an L1 box face.
  Hidden level-l walls: leaf when not covered by a level-(l+1) box (their covered fraction goes to the child).
- For each level l >= 1 hidden wall w with parent p (level l-1): the contribution
  `S(w) = sum_over_leaf_faces(w) value*A` and the parent carries **`-FR*p`** where FR = A_covered_by_children / A_p is evaluated per
  level pair. Chain invariant: the terms telescope. Total = sum over leaf walls only. Implementation: a single loop over levels from the finest
  to level 1; for each level add `+ sum leaf-part(w)` and subtract the covered share of the parent. Because each level's covered fraction
  is taken against its own parent, a 3-level pyramid gives L0 - FR1*p0 + (L1 - FR2*p1) + L2 = leaf set, with no double count.
- N faces (no parent): add their value, no parent share (R3). Thin plane walls: both IORs, as for any wall.
- Covered fractions use fine **gas** faces only (R3 closes findings 5b/6), now from level l+1 box faces under the level-l wall.
- Quantities: Q_CONV, Q_COND, Q_DIFF (UPDATE_HRR signs, S3 finding 1). The same leaf loop feeds: **wall devices** (`FINE_DEVC` iterates
  all `FWALL(LEV)` levels; covered L(l) walls are excluded by `LEAF_WALL_COVERED` on any finer box), **MLR** (burner mass loss from leaf walls;
  parent M_DOT_G_PP_ADJUST conserved per R3/S1), **gate inflow** (the composite mass gate counts leaf-wall mass flux; covered coarse faces
  contribute none, case N burners count at their level) and **wall-statistics** hooks (`LEAF_WALL_COUNT/LEAF_WALL/LEAF_WALL_COVERED` take a
  level argument and apply the "no finer box covers it" test).
- `TimeLoop::composite_energy`: remove the `lc.size()==2` branch and the three-level warning; call `fds_fine_walls_sums(l, ...)` for
  l = 1..finest and add the vector; `Q(12)` diagnostic remains L0-based (covered fraction per level printed as `nwall_l`).
- Reduction: one `ReduceRealSum` over the summed vector, unchanged.

## R14. Regrid, bind order, derefinement (D-N3)

1. **Order.** Walls are a function of (solid cells at each level, parent chain). Rebuild a level's walls whenever its box set or its parent's
   box set changes. Processing order: finest destroyed first (already enforced by IERR 1), rebuilt coarse to fine. A 1-level change at l
   forces re-bind of l+1.. because their parent walls (PNM/PIW) point into level l's wall arrays. Do not carry parent indices across a regrid.
2. **Derefinement of level l+1 with OBST re-snap.** The level-l walls of the former parent become leaf; their covered fraction goes to 0.
   Nothing is carried. Fluid state follows R10 item 4 (intensive state preserved, mass changes by volume, reported).
3. **Intermediate-level thin planes.** A thin plane on L1 under L2 is bookkeeping; the L2 plane/solid is the active one. When L2 is dropped,
   the L1 plane must exist (R9 rule: every OBST on every level), so no gap appears. Test below.
4. **Moving boxes** still forbid thick/ONE_D faces (S1). For static hierarchies they stay on L0 (D-010).

## R15. Zones, volume, level-0 identity (D-N4)

- Zone IDs: R10 rule 1 applies per level through the chain: level l inherits from level l-1 (not directly from L0), seeded flood fill inside
  the box, error on a fine cell connecting two zones. IERR 3 is replaced by this.
- Volume: composite, each gas cell at its finest covering level (R10 rule 2). Generalizes by "finest level covering this cell" without
  change; use the same covering mask as the leaf mask in R13, so zone volume and wall sums share one definition.
- Level 0 bit-identity: unchanged (R5). L0 code paths are not touched; `MESH_OBJ` null outside fine routines. Single-level run byte compare
  and OMP=1 two-level comparison against the previous HEAD still hold for the lev==1 cases; the loop must reproduce today's L1 inventory
  exactly (test M1).

## R16. Solid Phase: ONE_D/solid state and history for new fine faces (D-N5) [SP]

Based on S1/S2 and this note; the Solid Lead has not confirmed. Proposed, all level independent:

- **Nothing is stored above L0.** TMP_F, Q_CON_F, HTC, U_TANG, TMP_G, RHO_G, ZZ_G, B2 are recomputed from fine gas each step (S1 column 1/2).
  History (T_IGN, ramp clocks, M_DOT_G_PP_ADJUST, Q_RAD_IN) is read from the **L0 source wall** found through the chain (R12). For a new fine
  face (no source): `SF%T_IGN`, adjust 1, `EMISSIVITY*SIGMA*TMPA^4` (R3 rules, unchanged).
- **ONE_D / thermally thick on fine faces:** unsupported at any level above 0 (error), as S1. A thick OBST face refined to level l stays
  coarse-owned: its L0 wall provides the thermal BC; the fine face reads TMP_F from the L0 source. This is the "static hierarchy" reading of D-010;
  moving boxes over thick walls are an error. If the Solid Lead later wants fine ONE_D, the history is the problem (new fine faces have no
  profile): the only well-defined init is "copy and interpolate the L0 profile", which is a separate decision (not proposed here).
- **Upstream layout fit (S2):** the new work adds fields only to the driver-side `FWALL` records (PLEV, SRCNM, SRCIW). No new upstream hunks.
  The existing three hunks stay one line each; if upstream splits ONE_D and ALLOCATE_STORAGE changes, apply S2's fallback (allocate BOUNDARY_*
  directly). The chain walk reads only `B1_INDEX` and named B1 fields of the L0 wall (S2 rule).
- **Questions for the Solid Lead** (add to R8): (6) is "read L0 source wall at any level" acceptable for 3+ levels? (7) For thin OBST on a coarse
  level that has no burner state, is `SURF` input enough? (8) Confirm that no ONE_D on levels >= 1 is acceptable for the owner's sealed-compartment
  tests (these are inert/isothermal, so the answer should be yes).
- Coordination: the Solid Lead may not respond; the design above is complete without new demands on solid code. If no answer arrives, the
  default is as above (nothing stored above L0, thick = error).

## R17. Sizing, files, milestones, tests

**Estimate (builds on R7, whose A+B steps are about 5.5 days with R9/R10):** the N-level generalization itself adds the parent-chain helper (about
60 lines Fortran), the loop in `composite_energy` and leaf hooks (about 40 lines), the zone chain (already in R10), and per-level wall counts.
**About +2 days**, so about **7.5 days total** (R7 5.5 + 2), and 1 day of that is test infrastructure. Review and Solid Lead round trip +0.5 day,
done in parallel. Fits on the upstream layout change because only the three driver routines read wall storage.

| Milestone | Content | Test (pass) | Days |
|---|---|---|---|
| M0 | wait for the crash fix; add the chain helper without behavior change (lev=1 uses it) | L1 inventory identical to HEAD (wall list compare); L0 bit-identical; `--levels 2` Q lines unchanged | 0.5 |
| M1 | level-looped builder + exports (R4), no thin yet, P/N/D at any level | **3-level pyramid** (`pyramid_48_4m.fds` variant, half-cell steps at each level): per-IOR wall counts per level equal a single-mesh FDS at that dx; no segfault at L2; closed-domain mass PASS | 1.5 |
| M2 | leaf sums for N levels (R13) incl. devices, MLR, gate | 3-level static box: composite Q_CONV/COND/DIFF identity: sum over leaf = L0 - covered + fine (exact to round-off); area check: sum of leaf wall areas = same case with levels collapsed; remove 3-level warning | 1.5 |
| M3 | thin planes per level (R9) + zones through the chain (R10) | thin baffle `XB=1.80,1.83,...`: face count per level = FDS at that dx; zone IDs match covering parent cell by cell; zone without burner D_PBAR_DT = 0 | 1.5 |
| M4 | regrid: derefinement and OBST re-snap (R14, R10 item 4) | **Owner primary test (R6.0)** with 3 levels, `--regrid-every 2` moving tags: T, Z, PBAR unchanged to round-off every regrid, zones separate, REGRID ZONE-MASS diagnostic equals volume x rho; shifted-wall variant | 1.5 |
| M5 | burner / N faces (R3) at levels >=2, `--tag-obst-snap` | burner shifted half fine cell at L2: mass flux HRRPUA*A/dH within 1e-12; gate PASS | 1 |
| review | Solid Lead and Role 2 | questions R16 | 0.5 |

Always-on gates per milestone: L0-only run byte identical; OMP=1 two-level lines/devc against previous HEAD; unit tests (`driver_unit_tests`).

## R18. Risks

1. **Parent-chain across regrid** (high): stale PNM/PIW if a level l-1 is rebuilt without rebuilding level l. Mitigation: rebuild l..finest together;
   assert parent mesh/wall index in range at build; never cache across regrid.
2. **Leaf double count at 3 levels** (high): wall covered partly by L2 but wall of L1 itself covered partly L1->L0. Mitigation: telescoping per level
   pair (R13), tested by the exact identity in M2 against a collapsed hierarchy.
3. **Crash fix collision** (closed): the crash cause is the other-rank parent mesh (R20.1), not a missing L0 parent; the per-level builder does not read
   parent walls to decide which walls exist, so it cannot recur.
4. **Thin planes on intermediate levels** (medium): L1 thin but L2 solid, interface cuts. Mitigation: `--tag-obst-snap` default on for
   hierarchies with two or more levels (see D4; owner 2026-10-15).
5. **Upstream WALL/ONE_D change** (medium): see S2; hidden fine-wall hunks are 3 lines; fallback documented.
6. **Zone inheritance error cases** (medium): fine cell joining two zones, or rank holding multiple zones (IERR 3 lifted) with MPI: zones across
   ranks need the same IDs; start single-rank, MPI after.
7. **Rod/point OBSTs dropped on L0** (low/known): D1.
8. **Derefinement physics**: intensive state kept, mass changes (D9); repeated regrids drift total mass by design; diagnostics only.
9. **Performance**: per-step source-wall reads O(1); builder O(boxes) per level; negligible.
10. **No ONE_D above L0**: limits thermally thick fine walls to the L0 owner; an owner who wants fine thick walls needs a new design.

## R19. Status and recommendation for each open decision

| D | Status | Recommendation |
|---|---|---|
| D1 rod/point OBST dropped on L0 | Owner decided (2026-10-15) | Accept exception: FDS drops them on L0; always tag them to the finest level they exist on and log. Cannot be fixed without touching L0. |
| D2 L0 never re-snapped/masked | Owner decided (implied by L0 bit-identity; confirmed with 2026-10-15 wall decisions) | Yes: L0 stays FDS's own. |
| D3 levels >=1 follow FDS snap, thin as blocked faces with two walls | Agreed in note, not yet implemented | Adopt; N-level: each level uses its own dx (R9). |
| D4 `--tag-obst-snap` | Owner decided (2026-10-15) | Default **on** for two or more levels (mismatch can appear with a single refinement jump); off for single-level. |
| D5 no-parent faces history source = nearest L0 wall of same OBST/IOR, else SURF input | Needs Solid Lead sign-off (R8 1-2, R16) | Adopt with the chain walk; fall back to SURF input; no response = proceed as designed. |
| D6 live fine walls in fine kernels | Out of scope, separate | Keep hidden walls; design leaves a clean upgrade path (leaf walls are already the right set). |
| D7 zone IDs from L0, inherited per level, never renumbered | Agreed in note | Adopt; implement via level-(l-1) chain (R15), lift IERR 3. |
| D8 composite zone volume from finest level | Agreed | Adopt; same mask as the R13 leaf definition. |
| D9 derefinement keeps T/Z/PBAR/D, mass change diagnostic only | Owner decided | Adopt as written; gate test M4. |

N-level follow-ons (owner 2026-10-15 and prior rules):
- **D-N1** one level-looped builder with chain parents (R12): **adopt** (uniform wall build on every level; no level-1-only path).
- **D-N2** leaf sums by telescoping per level pair (R13): **adopt** (same leaf-integral rule as owner 2026-10-14).
- **D-N3** rebuild l..finest together at regrid, nothing cached (R14): **adopt as written**.
- **D-N4** zone inheritance via level l-1 (R15): **adopt** (zones stay separate on every level).
- **D-N5** nothing stored above L0, thick/ONE_D unsupported above L0 (R16, [SP]): **owner decided** (2026-10-15); rods/points on L0 only.

Still open: **D5** Solid Lead sign-off on no-parent history source (proceed as designed if no response).

## R20. M0/M1 as built (Role 1, 2026-10-16)

### R20.1 Real root cause of the level-2 crash

Reproduced on the box with the HEAD tree minus the rank guard of `ca61ac9f01` (`IF (PROCESS(NM0)/=MY_RANK) CYCLE`): `simple_test` split into two L0
meshes, `mpirun -np 2 ... --two-level-run --levels 3` (static patches, no regrid needed): rank 1 segfaults inside `fds_fine_walls_build` (SIGSEGV,
rc 11). A single-mesh 3-level run does not crash (it ran to the end at HEAD). Cause: the builder looked up the fine cell's L0 parent mesh by
coordinates (`MESHES(L)%XS..XF`, known on every rank) and then read `M0%WALL(...)` / `M0%CELL(...)%WALL_INDEX`. On a rank that does not own `M0`
the `WALL` array is not allocated (only the mesh extents and `CELL` exist). Fine boxes follow the AMReX distribution, not the L0 owner, so any fine
box (level 1 too, with several L0 meshes and ranks) can sit over a mesh of another rank. It looked "level-2 only" because `--levels 3` is the first
configuration in which the fine distribution stops following the L0 owners. It has nothing to do with a missing L0 wall for interior faces.

### R20.2 Per-level builder (no level-1 special case)

`FDS_FINE_WALLS_BUILD(LEV,IB)` is the same code for every level >= 1 and does not copy the parent's wall inventory:

1. `FDS_FINE_SET_LISTS` (called once from `load_obst_list`, `ObstList.H`) stores the global OBST list (XB, `SURF_INDEX(-3:3)`, from `fds_obst_xb/fds_obst_surf`)
   and the global VENT list (XB, IOR, SURF_INDEX, BOUNDARY_TYPE, `fds_vent_list`, includes the MIRROR/OPEN mesh-boundary vents), gathered over ranks.
2. `SNAP_SETUP(box)` snaps both lists at the box's own dx with FDS rounding `NINT((X-X0)/DX)` (same as `pb::obst_cells`), and sets the box offset `G0` and the
   domain size `NG` at that level.
3. For every gas cell and face D, `FACE_SURFACE`: a wall exists if the cell across is an OBST cell at this level or lies outside the domain; box-box and
   coarse-fine faces get none. Surface = last covering VENT (plane and extent snapped at this dx, IOR = -D), else the OBST face surface (last OBST containing
   the cell across), else `SURF_DEFAULT`. Boundary types other than SOLID/OPEN are skipped as before.
4. Thick (ONE_D) walls, `VEL_T`, `TMP_GAS_FRONT` and the log-law HTC model above L0 stop with a clear error (owner decision D-N5).
5. Parent link: `PNM/PIW` is the L0 wall with the same surface and boundary type whose face contains the fine face, only if that L0 mesh is on this rank; otherwise 0.
   It is a state source only (TMP_F, radiation, burner state), never an inventory. Consumers that need it skip `PNM==0` (D5 history source still open).
6. Rebuild at regrid is unchanged: `bind_level` rebuilds the level and everything above it; nothing is cached across regrids.

`--tag-obst-snap` (default ON for two or more levels, `--tag-obst-snap 0` off): on each level below the top, every OBST whose snap at level l+1 is not the
refinement of its snap at level l (including OBSTs thin at level l) is force-tagged (rounded outward, grown by 1 cell), so a coarse-fine interface does not cut it.

Files (M1): `fds_fine_level.f90` (builder, +165/-33 lines, consumers guard `PNM==0`), `fds_mesh_query.f90` (+38, two exports), `ObstList.H` (+30),
`RegridAmrCore.H/.cpp` (+18), `TwoLevelRun.H/.cpp`, `main.cpp` (option). No change to `init.f90`, `WALL_TYPE`, `ONE_D`; patch 0026 unchanged. Env `FDS_AMR_WALL_DIAG=1`
prints one `WALLDIAG level L box B g0 ... nwall N perD(...) withparent K` line per box build.

### R20.3 Evidence (box, small runs)

- L0 `--run 20`: all 15 `run_final_*.bin` plus mass/steps csv byte-identical to the HEAD binary.
- Two-level `simple_test` 100 steps: hrr and devc csv byte-identical to HEAD, 2304 fine walls in both.
- 3-level single mesh and 3-level two meshes with 2 ranks: run without the crash. Per-box wall counts by direction equal an independent Python snap of the same OBST list at
  that dx for 100/100 boxes (levels 1 and 2), also with an extra OBST (`XB=1.0,1.43,...`) that snaps differently on each level (91/91, `--tag-obst-snap 0` and 1).

### R20.4 Owner decisions applied (status of D1-D9)

D1 rods/points dropped on L0: accepted. D2 L0 never re-snapped: yes. D3: implemented for solid OBSTs at every level; **thin (zero-thickness) OBSTs are not yet built** (M3).
D4 default on for 2+ levels: implemented. D5 open (Solid Lead). D6 hidden walls only. D7/D8 zones and composite volume: not yet. D9 intensive state kept at derefinement: not yet (M4).
Thick and ONE_D walls above L0: unsupported, refused with an error. Leaf-only sums over the finest level covering each location: M2.
