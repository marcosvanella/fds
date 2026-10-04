// PressureIface.H: solver-agnostic pressure (Poisson for H) interface, single-level subset (milestone M1).
//
// The driver supplies the problem; a backend (FFT::Poisson or single-level MLMG) only performs the linear
// solve. Everything else (component labelling, pins, mean removal, gauge, true-residual check, removed-mean
// diagnostic) is the common layer and wraps every backend from outside. The interface, not the backend,
// defines the discrete operator: the 7-point cell-centred Laplacian (see CommonLayer.H, apply_operator).
// Backend fluxes are never used.
//
// Two APIs share one entry point:
//  - single level (milestone M1): the fields ba/dm/geom/rhs/phi of PressureProblem, FFT or MLMG backend;
//  - composite (milestone M2, first version): `PressureProblem::levels` non-empty. Then the list is authoritative
//    and the single-level ba/dm/geom/rhs/phi/nlevels are ignored. Only the MLMG backend does composite
//    (amrex::MLMG over MLPoisson on the whole hierarchy). See frozen/composite-notes.md.
//
// Not built (returned as Status::NotBuilt, never a silent fallback): masked cells (obstructed/solid, known-value
// or pinned cells), covered cells on the single-level API, non-uniform cell widths, non-Cartesian
// (cylindrical, RZ) geometry, variable coefficients, mixed open/closed domain faces (single level and
// composite), anisotropic refinement ratios (except ratio 1 in a one-cell direction), ratios other than 2 and 4,
// composite on the FFT backend, composite gauge weight/offset, and a multi-level request through the
// single-level fields (`nlevels` > 1 with an empty `levels`).
#ifndef PB_PRESSURE_IFACE_H_
#define PB_PRESSURE_IFACE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace pb {

enum class BC { Neumann, Periodic, Dirichlet };   // Dirichlet: homogeneous in M1 (values folded into the RHS)
enum class BackendKind { Auto, FFT, MLMG };
enum class Status { Ok, NotBuilt, InvalidInput, NotConverged };

// Domain face index: side*3 + dir (side 0 = low, 1 = high; dir 0..2 = x,y,z).
inline int face_index (int dir, int side) { return side*3 + dir; }

// One level of a composite request. All fields are mandatory. One component, valid cells only.
// Levels are ordered coarse to fine. Level 0 covers its domain exactly (as the single-level API requires).
// Levels l > 0 are a union of boxes on the refined domain of level l-1 (domain(l) = refine(domain(l-1), ref_ratio)),
// aligned to the ratio (boxes coarsenable by ref_ratio), properly nested in level l-1 with a buffer of one
// coarse cell wherever the coarsened level is not at a non-periodic domain face. No mask is passed: the cells of
// level l covered by level l+1 are those under the BoxArray of level l+1 coarsened by levels[l+1].ref_ratio.
struct PressureLevel {
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::Geometry geom;                         // this level's Cartesian geometry; periodicity must match `bc`
    amrex::IntVect ref_ratio{1};                  // ratio to the next COARSER level; ignored on level 0
    amrex::MultiFab const* rhs = nullptr;         // in: right-hand side, valid cells (covered cells are ignored)
    amrex::MultiFab* phi = nullptr;               // out: valid cells; in: initial guess if use_initial_guess
};

class PressureWorkspace;                          // opaque, below

struct PressureProblem {
    int nlevels = 1;                              // single-level API: > 1 is a composite request through these fields: NotBuilt
    // Composite hierarchy. Non-empty = authoritative (the single-level ba, dm, geom, rhs, phi and nlevels are
    // ignored); the flags below (bc, cylindrical, cell_width, coefficients, masks, gauge fields) are still checked.
    std::vector<PressureLevel> levels;
    amrex::BoxArray ba;                           // level 0
    amrex::DistributionMapping dm;
    amrex::Geometry geom;                         // dx from geom; periodicity must match `bc`. Cartesian only.
    // Optional cell widths per direction (one entry per cell of the domain). Empty = uniform, taken from geom.
    // Widths that are not all equal are NotBuilt: FFT::Poisson and MLPoisson assume uniform spacing. Widths that
    // are all equal must agree with geom (else InvalidInput).
    std::array<std::vector<amrex::Real>,3> cell_width;
    // FDS &MESH CYLINDRICAL=T (x is the radius; matrix rows and RHS are scaled by the radius factor RC). The 3-D
    // AMReX geometry cannot express it, so the driver sets this flag; true is NotBuilt.
    bool cylindrical = false;
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

    // Optional gauge inputs for singular components (FDS SYMM_INDEFINITE gauge, pres.f90 ~1830-1878): after the
    // solve, phi -= sum(w*(phi - g)) / sum(w) with w = cell volume * gauge_weight and g = gauge_offset. FDS passes
    // gauge_weight = RHOP (RHO predictor, RHOS corrector) and gauge_offset = KRES. Null = weight 1 / offset 0,
    // which is the plain volume-weighted mean removal. Same BoxArray/DistributionMapping as `ba`, one component.
    amrex::MultiFab const* gauge_weight = nullptr;
    amrex::MultiFab const* gauge_offset = nullptr;

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
    int verbose = 1;                              // 0 silent, 1 prints WARNING and error lines; >= 3 (composite) also the MLMG iteration history
};

struct ComponentInfo {
    int id = 0;
    bool singular = false;
    amrex::Long ncells = 0;
    amrex::IntVect pin{0};                        // lowest global index cell; recorded, not applied by MLMG/FFT
    double removed_mean = 0.0;                    // mean actually subtracted from the RHS (0 when at round-off)
    double removed_rel = 0.0;                     // |measured mean| relative to rms of the RHS (component); drives the diagnostic
    double gauge_shift = 0.0;                     // weighted mean (of H - gauge_offset) removed from H after the solve
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
    double residual_rel2 = 0.0;                   // ||b - L H||_2 / ||b||_2 (mean-removed b); composite: volume-weighted over uncovered cells
    double residual_relmax = 0.0;                 // max|b - L H| / max|b|; composite: over uncovered cells
    bool residual_ok = true;
    std::vector<std::string> warnings;
    bool workspace_rebuilt = false;               // composite: the workspace passed in did not match and was rebuilt
    int nlevels = 1;                              // composite: number of levels solved
    amrex::Long ncells_uncovered = 0;             // composite: number of uncovered cells over all levels
};

// Selector: FFT only on a single level with no masked cells and a uniformly closed (Neumann/periodic) or
// uniformly open (Dirichlet) set of domain faces; everything else is NotBuilt. Collective.
struct Selection { bool ok = false; BackendKind kind = BackendKind::Auto; std::string message; };
Selection select_backend (PressureProblem const& p, BackendKind requested);

// Hierarchy-dependent solver setup for the composite path (D-058): the covered-cell masks, the cell volumes and the
// MLMG hierarchy (MLPoisson with its multigrid levels and bottom solver). Build it once per hierarchy and pass it to
// solve_pressure() for every solve until the hierarchy (BoxArrays, DistributionMappings, geometries, ratios or BC)
// changes; then call rebuild(). It stores copies of the layout only: no pointer or reference into any rhs/phi
// MultiFab of the problem survives a call, so the old MultiFabs may be destroyed at once after a regrid.
// If a solve is given a workspace that does not match the problem it is rebuilt automatically
// (PressureResult::workspace_rebuilt = true), so a stale workspace is never used.
class PressureWorkspace {
public:
    PressureWorkspace ();
    ~PressureWorkspace ();
    PressureWorkspace (PressureWorkspace&&) noexcept;
    PressureWorkspace& operator= (PressureWorkspace&&) noexcept;
    PressureWorkspace (PressureWorkspace const&) = delete;
    PressureWorkspace& operator= (PressureWorkspace const&) = delete;

    // Discard the old setup and build it for `p.levels`. Collective. Returns Ok, or NotBuilt/InvalidInput with
    // `message` (the same checks as the solve). The workspace is empty after a non-Ok return.
    Status rebuild (PressureProblem const& p, std::string* message = nullptr);
    bool built () const;
    int num_levels () const;
    // True when the workspace was built for the layout of `p` (ba, dm, geom, ratios, bc).
    bool matches (PressureProblem const& p) const;

    struct Impl;
    Impl* impl () { return m_impl.get(); }
private:
    std::unique_ptr<Impl> m_impl;
};

// Full solve: selector, components and pins, mean removal, backend, gauge, checks. Collective. On any
// status other than Ok/NotConverged the phi field is not modified. `ws` is used only by the composite path
// (null = build the hierarchy setup and discard it after the call).
PressureResult solve_pressure (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws = nullptr);

// Face-centred dphi/dn of the composite solution, one component, valid faces, for every level, in the sign
// convention grad = (phi(+) - phi(-)) / dx in the direction of increasing index. grad_out[l][d] must be nodal in d
// on level l's BoxArray (amrex::convert(ba, IntVect::TheDimensionVector(d))), one component, no ghost cells.
// Uses p.levels[l].phi (valid cells; the hierarchy is as in the solve). Domain faces: Neumann gives 0, Dirichlet
// (homogeneous) 2*phi_in/dx outward-signed as in the operator, periodic wraps. Same-level faces and the C/F
// interpolation of the fine ghost cells are those of the solve (AMReX MLMG, maxorder 2). Afterwards every coarse
// face that is covered by fine faces (the C/F faces and the coarse faces under the fine level) is overwritten by
// amrex::average_down_faces of the fine gradients (D-032 (3)). For a one-cell direction the gradient in that
// direction is 0. Collective. Returns the same status values as the solve (grad_out untouched when not Ok).
PressureResult face_gradient_composite (PressureProblem const& p,
                                        std::vector<std::array<amrex::MultiFab*,3>>& grad_out);

const char* to_string (Status s);
const char* to_string (BC b);

} // namespace pb
#endif
