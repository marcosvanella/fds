// GhostExchange.H: ghost fill at the MESH_EXCHANGE positions of MAIN_LOOP and the boundary-condition step for the periodic and
// same-level box-to-box case (M2a, S4).
//
// Kernel-facing rules (M2a):
//  (a) Passive scalars: N_TOTAL_SCALARS beyond the tracked species are handled by Fields.cpp; ZZ/ZZS carry ncomp = N_TOTAL_SCALARS and every
//      exchange here moves all components of a registered field, so a passive scalar needs no special case.
//  (b) Only uniform Cartesian metrics are used; R(I)/RRN(I) are dropped, CYLINDRICAL and TRN* meshes are rejected at level-0 assembly (IR-002).
//
// Exchange codes of MAIN_LOOP (main.f90 MESH_EXCHANGE, verified at S4; README "Ghost fill at the MESH_EXCHANGE positions"):
//   1  after the predictor density update     RHOS, ZZS (+ MU, KRES, D, Q of the sender)
//   4  after the corrector density update     RHO,  ZZ  (+ MU, KRES, DS, Q)
//   3  after the predictor pressure solve     US, VS, WS, HS
//   6  after the corrector                    U,  V,  W,  H
//   5  inside the pressure iteration          FVX, FVY, FVZ, H (predictor) or HS (corrector)
// The other codes (0, 2, 7-11, 14-20) are radiation, obstruction, HVAC, particle and level-set exchanges, outside M2a.
//
// Boundary-condition step (bc_step): FDS gets the values across a mesh boundary from OMESH(NOM) (init: main.f90 INITIALIZE_MESH_EXCHANGE_1) and the
// unmodified boundary routines VISCOSITY_BC and VELOCITY_BC write them into the ghost cells (edge and corner cells, boundary-face strips, MU/KRES
// clamps). The driver fills OMESH(NOM) of every local box from the box data (block copy of MESH_EXCHANGE) and calls those routines: the ghost values are FDS's own by
// construction. Data of boxes on other ranks is broadcast (M2a: O(domain) traffic; the guarded patches 0003/0004 are the hooks to replace this by the
// AMReX ghost data without OMESH storage).
#ifndef FDS_GHOSTEXCHANGE_H_
#define FDS_GHOSTEXCHANGE_H_

#include <functional>
#include <string>
#include <vector>

#include "Fields.H"

namespace fdsamr {

// Fields moved by MESH_EXCHANGE code `code` (1, 3, 4, 5 or 6), as registered names. Empty for any other code. `predictor` selects H or HS for code 5.
std::vector<std::string> exchange_fields(int code, bool predictor = true);

// FillBoundary of exchange_fields(code): faces, edges, corners, periodic images. Physical-boundary ghosts are not touched.
// H/HS are exchanged like FDS does (their ghosts hold the pressure solver's boundary values, which belong to the pressure backend: S6).
void ghost_exchange(Fields& F, int code, bool predictor = true);

// Boundary-condition step that FDS runs right after exchange `code` (1, 3, 4, 6; other codes: nothing). Fills OMESH(NOM) of every local box from the
// box data (ghost fill of `code` must have been done), then calls the unmodified routine of MAIN_LOOP:
//   code 3 (predictor end):  MATCH_VELOCITY, VELOCITY_BC(estimated)     code 6 (corrector end): MATCH_VELOCITY, VELOCITY_BC(final)
//   code 1 (predictor):      VISCOSITY_BC(final velocities)             code 4 (corrector):     VISCOSITY_BC(estimated)
// followed, for codes 1 and 4, by WALL_BC (main.f90 886 / 1057) at the caller's convenience through wall_bc(). t, dt: time and step of the stage.
class BcStep {
public:
    BcStep(const Level0& l0, Fields& F);
    // Optional hook called with true/false around the velocity boundary routines: the driver uses it to turn box-interface walls into no-wall faces while
    // VELOCITY_BC runs (its edge loop would otherwise use the OMESH interpolation branch instead of the regular edge formulas of the single mesh).
    bool ext_ghost = false;   // EXTERNAL_GHOSTS_FILLED mode (patches 0003/0004): MATCH_VELOCITY is skipped by FDS, the driver saves the face velocities instead
    std::function<void(bool)> iface_hook;
    void after_exchange(int code, double t, double dt);
    void wall_bc(int predictor, double t, double dt);   // WALL_BC of the stage (all local boxes)
    // Kernel-check replay: the box data is ALREADY matched (a frozen dump state), so MATCH_VELOCITY is not run again; the ghost-producing routines only:
    // OMESH fill, VISCOSITY_BC (both stages), VELOCITY_BC (both stages), WALL_BC (predictor, then corrector).
    void replay(double t, double dt, bool predictor);          // + WALL_BC of the record's stage
    void replay_velocity(double t, double dt, bool predictor);  // without WALL_BC
    // OMESH(NOM) of every local box from the current box data (the block copy of MESH_EXCHANGE, all registered fields). The time loop calls it after each
    // ghost_exchange(5): NO_FLUX, MATCH_VELOCITY_FLUX and COMPUTE_VELOCITY_ERROR read H/HS and FVX.. of the neighbour from OMESH.
    void exchange_om() { fill_omesh(); }

private:
    const Level0& m_l0;
    Fields& m_F;
    bool local(int i) const;
    void fill_omesh();
};

}  // namespace fdsamr

#endif
