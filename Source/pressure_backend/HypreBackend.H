// HypreBackend.H: assembled-matrix HYPRE backend (ADR-002 backend (H)): an IJ matrix of the same discrete operator as the
// MLMG backend, solved with PCG (or GMRES / BiCGSTAB) preconditioned by BoomerAMG. Internal header (not part of PressureIface.H).
//
// Operator. Unknowns are the uncovered cells of all levels; the equations are the composite equations of the MLMG backend
// (frozen/composite-notes.md), built here cell by cell from the layout alone:
//  - same level, both cells unknowns: (phi_n - phi_c)/h^2 per face (7-point, second order);
//  - physical faces: Neumann ghost = phi_c, homogeneous Dirichlet ghost = -phi_c, periodic wraps; a direction of one cell has
//    no term (D-057: the y term of a one-cell y is dropped whatever its BC);
//  - fine cell next to a coarse/fine (C/F) face (refinement ratio r in the face-normal direction): MLMG's ghost value is the
//    linear (maxorder 2) extrapolation through the boundary datum b located at the centre of the coarse cell, r/2 fine cells outside
//    the face (MLMGBndry bloc), and the first interior fine cell: ghost = a*b + (1-a)*phi_c with a = 2/(r+1)  (a = 2/3 for r = 2,
//    2/5 for r = 4). b is interpolated from the coarser level with AMReX's order-3 sliding-parabola stencil
//    (interpbndrydata_*_o3: centre + first and second differences along the tangential neighbours that are in the domain and not
//    covered by the fine level, + the mixed term), so every C/F ghost is a fixed combination of coarse unknowns;
//  - coarse cell next to a covered cell: the coarse face flux is replaced by the average of the fine face fluxes
//    a*(phi_f - b_f)/h_f (fine-flux consistency, MLMG reflux) over the r^(active tangential directions) fine faces.
// The matrix is the volume-scaled operator A = -V*L (positive definite, symmetric on one level, FDS ULMAT/UGLMAT scaling);
// the right-hand side is F = -V*b. Covered cells keep a decoupled diagonal row (rhs 0). For a singular problem (no Dirichlet face
// in a direction that carries an operator term) one unknown (the lowest-index uncovered cell of the coarsest level that has one)
// is pinned: its row is a (diagonal-scaled) identity with rhs 0 and its column is dropped from all other rows, as FDS ULMAT
// does. The mean removal before and the gauge after the solve are the common layer's, not done here.
// Host buffers only (HYPRE_MEMORY_HOST); the set-up (matrix, BoomerAMG hierarchy) is kept for repeated solves with the same layout.
#ifndef PB_HYPRE_BACKEND_H_
#define PB_HYPRE_BACKEND_H_

#include "PressureIface.H"

#include <memory>
#include <string>
#include <vector>

namespace pb {

struct HypreLayoutLevel {
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::Geometry geom;
    amrex::IntVect ratio{1};                      // to the next coarser level (ignored on level 0)
};

class HypreSystem {
public:
    // `ebc` are the effective boundary types (effective_bc()). pin_singular = false builds the bare operator (tests).
    // Collective. ok() is false with message() set when the layout cannot be built (e.g. more rows than HYPRE's integer range).
    HypreSystem (std::vector<HypreLayoutLevel> const& levels, std::array<BC,6> const& ebc, HypreOptions const& opt,
                 bool singular, bool pin_singular = true);
    ~HypreSystem ();
    HypreSystem (HypreSystem const&) = delete;
    HypreSystem& operator= (HypreSystem const&) = delete;

    bool ok () const;
    std::string const& message () const;
    bool matches (std::vector<HypreLayoutLevel> const& levels, std::array<BC,6> const& ebc, HypreOptions const& opt, bool singular) const;

    // Solve L phi = b for the unit-operator right-hand sides b (valid cells; covered cells ignored). phi: valid cells, initial guess
    // if use_guess, solution on return (covered cells untouched; the pinned unknown is 0 before the common layer's gauge).
    BackendStatus solve (std::vector<amrex::MultiFab*> const& phi, std::vector<amrex::MultiFab const*> const& b,
                         double tol_rel, int max_iter, bool use_guess);
    // y = A x on valid cells (volume-scaled operator A = -V*L; rows of covered cells: the diagonal scale times x). Test hook.
    void apply (std::vector<amrex::MultiFab const*> const& x, std::vector<amrex::MultiFab*> const& y) const;

    long long rows () const;
    long long nonzeros () const;
    double setup_seconds () const;
    amrex::IntVect pin_cell () const;             // level, not used for pins on level 0: (i,j,k) and pin_level()
    int pin_level () const;
    std::string method () const;

    struct Impl;
private:
    std::unique_ptr<Impl> m_impl;
};

} // namespace pb
#endif
