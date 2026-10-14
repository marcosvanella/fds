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
   coarse cells", opt-in `--tag-obst-snap NBUF` (default off, so existing baselines hold).
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
