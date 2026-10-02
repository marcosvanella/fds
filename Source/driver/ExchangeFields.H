// ExchangeFields.H: which registered fields each FDS MESH_EXCHANGE code moves. Header only and pure (std only: no AMReX, no FDS), so a test can use it without linking the driver
// (Role 3's stand-alone ghost test). Moved out of GhostExchange.cpp in S10; behaviour unchanged. Kernel-facing rules (M2a): (a) passive scalars are extra components of ZZ/ZZS
// (Fields.cpp), all components move; (b) only uniform Cartesian metrics are used.
#ifndef FDS_EXCHANGEFIELDS_H_
#define FDS_EXCHANGEFIELDS_H_

#include <string>
#include <vector>

namespace fdsamr {

// Fields moved by MESH_EXCHANGE code `code` (1, 3, 4, 5 or 6), as registered names. Empty for any other code. `predictor` selects H or HS for code 5.
inline std::vector<std::string> exchange_fields(int code, bool predictor = true)
{
    switch (code) {
    // TMP and RSUM: the ghost layer of a cell across a box interface. FDS gets them from WALL_BC/ASSIGN_GHOST_VALUE (OMESH average, patch 0004: skipped when the ghosts
    // are filled externally); the driver fills them from the neighbouring box like the other scalars and keeps the interface walls out of WALL_BC (TimeLoop::iface).
    case 1: return {"RHOS", "ZZS", "MU", "KRES", "D", "TMP", "RSUM"};
    case 4: return {"RHO", "ZZ", "MU", "KRES", "DS", "TMP", "RSUM"};
    case 3: return {"US", "VS", "WS", "HS"};
    case 6: return {"U", "V", "W", "H"};
    case 5: return {"FVX", "FVY", "FVZ", predictor ? "H" : "HS"};
    default: return {};
    }
}

}  // namespace fdsamr

#endif
