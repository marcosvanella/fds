# Fine-grid-conforming OBST wall cells on fine AMR boxes

Design note. Co-authored: Role 1 (lead: fine wall build, fluid coupling, composite sums) and Solid Phase Lead (wall state, upstream fit, review).

## Role 1 sections (placeholder)

_To be written by Role 1: wall-cell build on fine boxes (FDS_FINE_WALLS_BUILD), fine-conforming OBST faces that have no parent face, fluid-side BCs, composite HRR terms, test plan._

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
