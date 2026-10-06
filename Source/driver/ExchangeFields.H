// ExchangeFields.H: which registered fields each FDS MESH_EXCHANGE code moves. Header only and pure (std only: no AMReX, no FDS), so a test can use it without linking the driver
// (Role 3's stand-alone ghost test). Moved out of GhostExchange.cpp in S10; behaviour unchanged. Kernel-facing rules (M2a): (a) passive scalars are extra components of ZZ/ZZS
// (Fields.cpp), all components move; (b) only uniform Cartesian metrics are used.
#ifndef FDS_EXCHANGEFIELDS_H_
#define FDS_EXCHANGEFIELDS_H_

#include <string>
#include <vector>

namespace fdsamr {

// Fields moved by MESH_EXCHANGE code `code` (1, 3, 4, 5 or 6), as registered names. Empty for any other code (including the back-wall codes 66 and 67).
// `predictor` selects H or HS for code 5. Mapping table (FireX merge 1fb5c300 / FDS-AMReX tip after 8ad685f435):
//   1  RHOS, ZZS, MU, KRES, D, TMP, RSUM     4  RHO, ZZ, MU, KRES, DS, TMP, RSUM
//   3  US, VS, WS, HS                        6  U, V, W, H   (velocity/pressure ghosts only; no back-wall data)
//   5  FVX, FVY, FVZ, H|HS
//   66  full PACK_WALL / PACK_THIN_WALL (once at set-up, INITIALIZE_BACK_WALL_EXCHANGE) — FDS host path only
//   67  PACK_WALL_LITE (every step after MESH_EXCHANGE(6), and HT3D) — FDS host path only
// Codes 66 and 67 must never be routed through the driver ghost fill (ghost_exchange / BcStep::exchange / after_exchange): they return empty here and
// after_exchange returns immediately for any code other than 1, 3, 4, 6. See notes/exchange-codes.md.
inline std::vector<std::string> exchange_fields(int code, bool predictor = true)
{
    switch (code) {
    // TMP and RSUM: the ghost layer of a cell across a box interface. FDS gets them from WALL_BC/ASSIGN_GHOST_VALUE (OMESH average, patch 0004: skipped when the ghosts
    // are filled externally); the driver fills them from the neighbouring box like the other scalars and keeps the interface walls out of WALL_BC (TimeLoop::iface).
    case 1: return {"RHOS", "ZZS", "MU", "KRES", "D", "TMP", "RSUM"};
    case 4: return {"RHO", "ZZ", "MU", "KRES", "DS", "TMP", "RSUM"};
    case 3: return {"US", "VS", "WS", "HS"};
    case 6: return {"U", "V", "W", "H"};   // after FireX: back-wall left code 6; only U, V, W, H ghosts
    case 5: return {"FVX", "FVY", "FVZ", predictor ? "H" : "HS"};
    // 66, 67 and every other code: empty (host path; never through the driver fill)
    default: return {};
    }
}

}  // namespace fdsamr

#endif
