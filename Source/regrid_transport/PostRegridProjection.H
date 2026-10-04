// PostRegridProjection.H: the post-regrid composite projection (ruling D-063), behind an abstract composite Poisson solver. AMReX only.
//
// After a regrid that creates new fine faces next to retained ones, the face velocity of the new hierarchy does not satisfy div u = D at the seams between
// retained and prolonged faces (D-060 keeps the divergence of every PARENT cell, not of the fine cell next to a retained fine face). D-063: the driver then solves, over the
// whole composite hierarchy,
//        Lap phi = div u - D          (uncovered cells of all levels; covered cells carry no equation)
// subtracts grad phi from the face velocities, and lets the fine face values overwrite the coarse faces under them (the interface flux overwrite of ADR-002 v1.2).
// Acceptance: max |div u - D| over all uncovered cells is at most the pressure tolerance. The band re-prolongation is NOT used.
//
// project_after_regrid does the bookkeeping (average-down of faces, divergence, right-hand side, update, check) and calls a CompositePoissonSolver for the linear solve.
// The solver is the plug-in point of Role 2's composite MLMG: its entry points are PressureWorkspace::rebuild (hierarchy changed, keeps no pointers to old
// MultiFabs) and a PressureLevel vector. Those names are not in Source/pressure_backend yet (only the single-level PressureIface exists); ProjectionLevel and
// CompositePoissonSolver::rebuild/solve are the corresponding names on this side and will be renamed or wrapped when the headers land. Tests use the mock in
// tests/MlmgCompositeSolver.H (amrex::MLMG, periodic).
//
// Sign and discretisation contract of a solver: solve() returns, on every valid face of every level, the gradient of ITS OWN solution, g = (phi_hi - phi_lo) / dx in the
// normal direction, evaluated with the solver's own coarse-fine treatment (a fine face on a coarse-fine interface uses the solver's interpolated ghost value), so that
// the cell divergence of g equals Lap phi of its operator on every uncovered cell, where the coarse faces under a fine face are the area mean of the fine faces. That is
// exactly what makes div (u - g) = D hold to the solver tolerance after the fine faces have overwritten the coarse ones.
#ifndef FDSRT_POSTREGRIDPROJECTION_H_
#define FDSRT_POSTREGRIDPROJECTION_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace fdsrt {

// One level of the composite hierarchy (non-owning view). The face MultiFabs are nodal in their direction, one component, any ghost width (valid faces only are used).
struct ProjectionLevel {
    amrex::Geometry geom;
    amrex::IntVect ref_ratio = amrex::IntVect(1);          // ratio of this level to the next COARSER one (unused on level 0)
    std::array<amrex::MultiFab*, 3> vel{nullptr, nullptr, nullptr};   // face velocity, modified by the projection
    const amrex::MultiFab* D = nullptr;                    // target divergence of the cells (one component); nullptr = 0
    const amrex::iMultiFab* covered = nullptr;             // nonzero = covered by the next finer level; nullptr = nothing covered (finest level)
};

struct SolveReport {
    bool ok = false;
    int iterations = 0;
    double residual = -1.0;          // final residual norm in the solver's own measure (-1: none)
    std::string message;
};

class CompositePoissonSolver {
public:
    virtual ~CompositePoissonSolver() = default;
    virtual const char* name() const = 0;
    // The hierarchy changed: drop everything that refers to the old grids and set up for `levels` (collective). Keeps no pointers to the old MultiFabs.
    virtual void rebuild(const std::vector<ProjectionLevel>& levels) = 0;
    // Solve Lap phi = rhs[l] on the uncovered cells (rhs[l] is zero on covered cells; one component, cell-centred on the box arrays of the levels given to rebuild) and fill
    // grad[l][d] (nodal in d, one component, valid faces) as described in the contract above. Periodic or Neumann domain faces: the composite mean of rhs is zero.
    virtual SolveReport solve(const std::vector<const amrex::MultiFab*>& rhs, std::vector<std::array<amrex::MultiFab*, 3>>& grad, double tol_rel, double tol_abs, int max_iter) = 0;
};

struct ProjectionOptions {
    bool enabled = true;          // false: only measure (the "projection off" negative control, and the Phase 3 default of ruling (b) 5 until the real solver lands)
    double tol_rel = 1.0e-12;     // relative tolerance given to the solver (residual against the right-hand side)
    double tol_abs = 0.0;         // absolute residual tolerance given to the solver (0: none); set it below accept_abs so that a right-hand side of round-off size converges
    int max_iter = 200;
    double accept_abs = 0.0;      // acceptance bound on max|div u - D| after the projection: accept_abs + accept_rel * (max|div u - D| before)
    double accept_rel = 1.0e-8;   // (the pressure tolerance of the run enters here)
};

struct ProjectionReport {
    bool ran = false;                       // the solve was done (false: disabled, or max|div u - D| was already within accept_abs and the solve was skipped)
    bool solved = false;                    // the solver reported success
    std::string solver, message;
    double div_before = 0.0, div_after = 0.0;               // max |div u - D| over the uncovered cells of all levels (reduced over ranks)
    std::vector<double> div_before_level, div_after_level;  // by level
    std::vector<double> max_change_level;                   // max |grad phi| over the faces of the level (the velocity change), by level
    double rhs_sum_rel = 0.0;               // |sum over uncovered cells of volume * rhs| / (sum of volume * |rhs|): zero for a solvable periodic problem
    int iterations = 0;
    bool accepted = false;                  // div_after <= accept_abs + accept_rel * div_before  (and the solve happened)
};

// Composite check alone (after the fine faces have overwritten the coarse ones this does not modify the velocity): max |div u - D| over the uncovered cells, by level and overall.
// Collective. The caller must have done the face average-down (project_after_regrid does it).
double composite_divergence_error(const std::vector<ProjectionLevel>& levels, std::vector<double>* by_level = nullptr);

// Fine faces overwrite the coarse faces under them, finest level first (area mean). Collective.
void average_down_velocity(const std::vector<ProjectionLevel>& levels);

// The post-regrid projection. `solver` may be null only with opts.enabled == false. Collective. The velocity is modified only if enabled and the solve succeeded.
ProjectionReport project_after_regrid(std::vector<ProjectionLevel>& levels, CompositePoissonSolver* solver, const ProjectionOptions& opts);

}  // namespace fdsrt

#endif
