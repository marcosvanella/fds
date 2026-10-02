// PressureBcMap.H: FDS pressure boundary codes of one direction -> pb::BC of the two faces, as the driver passes them to pb::solve_pressure (S9, thin-direction check).
// Header only, pure (no FDS, no AMReX objects), so tests/test_units.cpp can test it. Kernel-facing rules (M2a): (a) passive scalars are not involved; (b) only uniform Cartesian metrics.
//
// FDS codes LBC/MBC/NBC (init.f90 ~2494, pres.f90 452-487): 0 periodic, 1 Dirichlet-Dirichlet, 2 Dirichlet(low)-Neumann(high), 3 Neumann-Neumann, 4 Neumann(low)-Dirichlet(high); 5, 6 (cylindrical axis) are outside the scope.
// amrex::FFT::Poisson ignores a direction of ONE cell (the cell contributes no coupling term): a Dirichlet face in it would be dropped, while the full discrete operator
// (FDS H3CZSS, or any other backend) has an extra diagonal term. FDS itself solves a one-cell y direction as 2-D (JBAR == 1 sets TWO_D, read.f90 703; H2CZSS has no y), so there the
// ignored direction IS the FDS operator and the faces are mapped to Neumann (even, no coupling, correct for every backend). A one-cell direction x or z is solved by FDS with its Dirichlet
// term; no backend here reproduces that through the FFT string, so a Dirichlet face there aborts with a message (a Neumann or periodic one-cell direction is exact).
#ifndef FDS_PRESSUREBCMAP_H_
#define FDS_PRESSUREBCMAP_H_

#include <string>

#include "PressureIface.H"

namespace fdsamr {

struct DirBc {
    pb::BC lo = pb::BC::Neumann, hi = pb::BC::Neumann;
    std::string error;   // non-empty: the case must be refused with this message
    std::string note;    // non-empty: the mapping changed the FDS code (documented, not an error)
};

inline char bc_char(pb::BC b) { return b == pb::BC::Periodic ? 'P' : b == pb::BC::Dirichlet ? 'D' : 'N'; }
inline std::string bc_string(const DirBc& d) { return std::string() + bc_char(d.lo) + bc_char(d.hi); }

// dir 0..2, code_lo/code_hi: the FDS code at the low/high domain face (as reduced over the boxes), ncell: cells of the level in this direction, level_periodic: the level is periodic here,
// fds_ignores: FDS does not solve this direction (TWO_D, y, one cell).
inline DirBc map_pressure_bc_direction(int dir, int code_lo, int code_hi, int ncell, bool level_periodic, bool fds_ignores)
{
    static const char* const dn[3] = {"x", "y", "z"};
    DirBc r;
    const std::string d = dn[dir];
    if (code_lo < 0 || code_hi < 0) { r.error = "no box found at a domain face (pressure boundary map)"; return r; }
    if (level_periodic) { r.lo = r.hi = pb::BC::Periodic; return r; }   // level-wide periodic whatever the per-mesh code (several meshes: the mesh-to-mesh interpolated boundaries the level solve replaces)
    if (code_lo == 0 || code_hi == 0) {
        if (ncell != 1) { r.error = "FDS pressure code 0 (periodic) on a non-periodic domain direction " + d + " with more than one cell (inconsistent, S6)"; return r; }
        r.lo = r.hi = pb::BC::Neumann;   // one cell: the Neumann and periodic operators coincide
        return r;
    }
    if (code_lo > 4 || code_hi > 4) { r.error = "pressure boundary code 5/6 (cylindrical axis) is outside the S5 scope"; return r; }
    const bool dlo = (code_lo == 1 || code_lo == 2), dhi = (code_hi == 1 || code_hi == 4);
    if (ncell == 1) {
        if (fds_ignores) {   // y of a TWO_D case: FDS has no y operator; Neumann adds nothing in any backend
            r.lo = r.hi = pb::BC::Neumann;
            if (dlo || dhi) r.note = "direction " + d + " has one cell and FDS does not solve it (TWO_D): Dirichlet face(s) mapped to Neumann (no coupling)";
            return r;
        }
        if (dlo || dhi) {
            r.error = "Dirichlet (open) pressure boundary on the one-cell direction " + d + ": amrex::FFT::Poisson ignores a one-cell direction but the FDS operator has the Dirichlet term there; not supported (use at least 2 cells or a closed/periodic direction)";
            return r;
        }
    }
    r.lo = dlo ? pb::BC::Dirichlet : pb::BC::Neumann;
    r.hi = dhi ? pb::BC::Dirichlet : pb::BC::Neumann;
    return r;
}

// The direction FDS does not solve (TWO_D y) does not enter the operator, so its face types are free. When every face of the two solved directions is Dirichlet (a fully open 2-D box) it takes
// Dirichlet as well, so that pb::select_backend sees six uniform open faces instead of refusing a "mixed open/closed" set that is not mixed in the operator; otherwise it stays Neumann.
// Only a Neumann ignored direction is changed (a periodic one is the level's property). Returns true when it changed.
inline bool ignored_direction_follows_open(DirBc& ign, const DirBc& a, const DirBc& b)
{
    if (ign.lo != pb::BC::Neumann || ign.hi != pb::BC::Neumann) return false;
    const bool open = a.lo == pb::BC::Dirichlet && a.hi == pb::BC::Dirichlet && b.lo == pb::BC::Dirichlet && b.hi == pb::BC::Dirichlet;
    if (!open) return false;
    ign.lo = ign.hi = pb::BC::Dirichlet;
    ign.note = "the ignored (one-cell, TWO_D) direction takes Dirichlet because both solved directions are open: the backend selector then sees six open faces (no effect on the operator: a one-cell direction has no coupling term)";
    return true;
}

}  // namespace fdsamr

#endif
