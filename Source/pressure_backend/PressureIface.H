// PressureIface.H: solver-agnostic pressure (Poisson for H) interface, single-level subset (milestone M1).
//
// The driver supplies the problem; a backend (FFT::Poisson or single-level MLMG) only performs the linear
// solve. Everything else (component labelling, pins, mean removal, gauge, true-residual check, removed-mean
// diagnostic) is the common layer and wraps every backend from outside. The interface, not the backend,
// defines the discrete operator: the 7-point cell-centred Laplacian (see CommonLayer.H, apply_operator).
// Backend fluxes are never used.
//
// Not built (returned as Status::NotBuilt, never a silent fallback): masked cells, composite (more than one
// level), covered cells, and an unmasked single level with mixed open/closed domain faces.
#ifndef PB_PRESSURE_IFACE_H_
#define PB_PRESSURE_IFACE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace pb {

enum class BC { Neumann, Periodic, Dirichlet };   // Dirichlet: homogeneous in M1 (values folded into the RHS)
enum class BackendKind { Auto, FFT, MLMG };
enum class Status { Ok, NotBuilt, InvalidInput, NotConverged };

// Domain face index: side*3 + dir (side 0 = low, 1 = high; dir 0..2 = x,y,z).
inline int face_index (int dir, int side) { return side*3 + dir; }

struct PressureProblem {
    int nlevels = 1;                              // > 1 is a composite request: NotBuilt
    amrex::BoxArray ba;                           // level 0
    amrex::DistributionMapping dm;
    amrex::Geometry geom;                         // dx from geom; periodicity must match `bc`
    std::array<BC,6> bc = {BC::Neumann, BC::Neumann, BC::Neumann, BC::Neumann, BC::Neumann, BC::Neumann};

    // Coefficients: M1 supports the unit operator only (nullptr). Non-null is NotBuilt.
    amrex::MultiFab const* cell_coef_a = nullptr;
    std::array<amrex::MultiFab const*,3> face_coef_b = {nullptr, nullptr, nullptr};

    // Masks (same BoxArray/DistributionMapping as `ba`). cell_class: nonzero = special (unknown, OBST, gap,
    // known-value or pinned) cell. uncovered: zero = covered by a finer level. Either one non-trivial: NotBuilt.
    amrex::iMultiFab const* cell_class = nullptr;
    amrex::iMultiFab const* uncovered = nullptr;
    // Driver-supplied component ids: M1 labels components itself; non-null is NotBuilt.
    amrex::iMultiFab const* component_id = nullptr;

    amrex::MultiFab const* rhs = nullptr;         // one component, valid cells only are read
    amrex::MultiFab* phi = nullptr;               // out: H on valid cells (gauge-fixed); in: initial guess if used
};

struct PressureOptions {
    BackendKind backend = BackendKind::Auto;      // Auto = selector rule; FFT/MLMG = request, still rule-checked
    double tol_rel = 1.0e-12;                     // backend relative tolerance
    int max_iter = 200;                           // backend iteration cap
    bool use_initial_guess = false;               // false: start from zero
    bool remove_mean = true;                      // test hook: false disables mean removal (singular solves then warn)
    bool check_residual = true;                   // true composite residual b - L phi of the mean-removed problem
    double residual_tol = 1.0e-12;                // warn above this relative 2-norm
    double removed_mean_warn = 1.0e-10;           // warn when removed mean (rel. to rms b) exceeds this
    int verbose = 1;                              // 0 silent, 1 prints WARNING and error lines
};

struct ComponentInfo {
    int id = 0;
    bool singular = false;
    amrex::Long ncells = 0;
    amrex::IntVect pin{0};                        // lowest global index cell; recorded, not applied by MLMG/FFT
    double removed_mean = 0.0;                    // mean actually subtracted from the RHS (0 when at round-off)
    double removed_rel = 0.0;                     // |measured mean| relative to rms of the RHS (component); drives the diagnostic
    double gauge_shift = 0.0;                     // mean removed from H after the solve
};

struct BackendStatus {
    bool converged = false;
    int iterations = 0;                           // 0 for direct backends
    double own_residual = -1.0;                   // backend-reported relative residual (-1: none)
};

struct PressureResult {
    Status status = Status::Ok;
    std::string message;                          // reason for any non-Ok status
    std::string backend;                          // "FFT", "MLMG" (empty if not built)
    BackendStatus backend_status;
    std::vector<ComponentInfo> components;
    bool residual_checked = false;
    double residual_rel2 = 0.0;                   // ||b - L H||_2 / ||b||_2 (mean-removed b)
    double residual_relmax = 0.0;                 // max|b - L H| / max|b|
    bool residual_ok = true;
    std::vector<std::string> warnings;
};

// Selector: FFT only on a single level with no masked cells and a uniformly closed (Neumann/periodic) or
// uniformly open (Dirichlet) set of domain faces; everything else is NotBuilt. Collective.
struct Selection { bool ok = false; BackendKind kind = BackendKind::Auto; std::string message; };
Selection select_backend (PressureProblem const& p, BackendKind requested);

// Full solve: selector, components and pins, mean removal, backend, gauge, checks. Collective. On any
// status other than Ok/NotConverged the phi field is not modified.
PressureResult solve_pressure (PressureProblem const& p, PressureOptions const& o);

const char* to_string (Status s);
const char* to_string (BC b);

} // namespace pb
#endif
