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
//    and the single-level ba/dm/geom/rhs/phi/nlevels are ignored. The MLMG backend (amrex::MLMG over MLPoisson on the
//    whole hierarchy) and the HYPRE backend do composite. See frozen/composite-notes.md.
//
// Backends: FFT (single uniform level), MLMG (default for composite; single level on request), HYPRE (assembled matrix, PCG or
// GMRES/BiCGSTAB with BoomerAMG, ADR-002 backend (H); on request only through `PressureOptions::backend =
// BackendKind::HYPRE`, never selected by Auto). HYPRE assembles the same discrete operator as MLMG (including the
// coarse/fine interface treatment) into an IJ matrix, see HypreBackend.H and frozen/hypre-notes.md.
//
// Not built (returned as Status::NotBuilt, never a silent fallback): masked cells on any path but the single-level MLMG (the masked
// branch, see PressureProblem::cell_class and frozen/masked-notes.md), covered cells on the single-level API, non-uniform cell widths, non-Cartesian
// (cylindrical, RZ) geometry, variable coefficients, a Dirichlet face in a one-cell x or z direction, anisotropic refinement ratios (except ratio 1 in a one-cell direction), ratios other than 2 and 4,
// composite on the FFT backend, and a multi-level request through the single-level fields (`nlevels` > 1 with an
// empty `levels`).
//
// Mean removal and gauge (D-067). Singular components (all faces Neumann or periodic) are made solvable by the common
// layer before the backend runs, and the solution constant is fixed after it, for every backend (FFT, MLMG, composite).
//  - Mean removal, default `MeanKind::Volume`: the composite volume-weighted mean of the right-hand side over the
//    uncovered cells of the component is subtracted (exact sum, one fixed-point scale, bitwise independent of
//    decomposition). `MeanKind::ScaledArithmetic` is the FDS parity switch: the arithmetic mean of the volume-scaled
//    right-hand side F = v*b is removed from F (b_k -= mean(F)/v_k), which is what FDS ULMAT/UGLMAT do. On uniform
//    cells both are the same operation up to rounding; they differ on a hierarchy (cell volumes differ between levels)
//    and on stretched cells, and only when the right-hand side is not already compatible.
//  - Gauge, always applied to a singular component: sum(rho*V*(KRES - H)) = 0 over the uncovered cells of the
//    component, with rho = gauge_weight (default 1) and KRES = gauge_offset (default 0; then the plain volume-weighted
//    mean of H is zero). A component that is not singular (some Dirichlet face) has a unique solution; nothing is
//    shifted. Exact sums, so the constant does not depend on the decomposition up to the round-off of H itself.
//  - Components: without masks the whole domain (all levels) is one component. With a masked single-level problem the components are
//    the flood-filled gas regions (gas-gas faces); mean removal and gauge are per component over its gas cells.
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
enum class BackendKind { Auto, FFT, MLMG, HYPRE };   // Auto never selects HYPRE
// Message of Status::NotBuilt when BackendKind::HYPRE is requested in a build without PB_WITH_HYPRE (README.md, "Build options").
inline constexpr const char* kHypreNotBuilt = "built without HYPRE";
enum class Status { Ok, NotBuilt, InvalidInput, NotConverged };

// PressureProblem::cell_class codes (frozen/masked-notes.md).
enum CellClass : int { CellGas = 0, CellSolid = 1, CellKnown = 2 };

// Which zero mode of the right-hand side is removed from a singular component (see the header comment and
// CommonLayer.H). Default Volume (D-067); ScaledArithmetic is the FDS parity switch.
//   Volume:           b -= sum(v*b)/sum(v)        (constant shift of b; the scaled F = v*b changes by c*v)
//   ScaledArithmetic: b_k -= mean(F)/v_k          (F = v*b, arithmetic mean over the uncovered cells; FDS ULMAT/UGLMAT)
enum class MeanKind { Volume, ScaledArithmetic };

// Domain face index: side*3 + dir (side 0 = low, 1 = high; dir 0..2 = x,y,z).
inline int face_index (int dir, int side) { return side*3 + dir; }

// Effective boundary types (D-057). A one-cell y direction is the FDS TWO_D case: the y term of the operator is dropped
// (FFT::Poisson sets its y factor to 0 for a length-1 direction whatever the BC, and a Neumann or periodic MLPoisson
// direction of one cell contributes 0 too). A Dirichlet face of a one-cell y therefore does not make the problem
// non-singular and does not change the operator: it is treated as Neumann by every backend, by the reference
// operator and by the singular-component test. Any other BC is unchanged. `dom0` is the (coarsest) domain.
// A Dirichlet face of a one-cell x or z direction is NOT equivalent (the term -2 phi/dx^2 stays) and is NotBuilt.
inline std::array<BC,6> effective_bc (std::array<BC,6> bc, amrex::Box const& dom0)
{
    if (dom0.length(1) == 1) {
        for (int side = 0; side < 2; ++side) {
            if (bc[face_index(1,side)] == BC::Dirichlet) { bc[face_index(1,side)] = BC::Neumann; }
        }
    }
    return bc;
}

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
    // Optional gauge inputs of this level (D-067): rho (FDS RHOP: RHO in the predictor, RHOS in the corrector) and KRES.
    // Same BoxArray/DistributionMapping as `ba`, one component, valid cells; only the uncovered cells are read.
    // Either both given on every level, or one of them on every level, or none (null = rho 1 / KRES 0). The gauge
    // constant is sum(V_l*rho*(phi - KRES)) / sum(V_l*rho) over the uncovered cells of all levels (one exact sum).
    amrex::MultiFab const* gauge_weight = nullptr;
    amrex::MultiFab const* gauge_offset = nullptr;
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

    // Masks (same BoxArray/DistributionMapping as `ba`; frozen/masked-notes.md). cell_class codes (CellClass): 0 Gas (unknown), 1 Solid
    // (OBST or gap cell: excluded, no flux through any of its faces), 2 Known (known-value cell behind a face-Dirichlet condition: not an
    // unknown, holds the face value `known_value`; the gas neighbour gets the flux 2 (g - phi)/dx^2). Null or all zero: unmasked problem
    // (unchanged code). Masked single-level problems are built for MLMG only (the selector picks MLMG for Auto; FFT, HYPRE and the
    // composite paths: NotBuilt). Any other code: InvalidInput. uncovered: zero = covered by a finer level, NotBuilt.
    amrex::iMultiFab const* cell_class = nullptr;
    amrex::iMultiFab const* uncovered = nullptr;
    // Known values g of the Known cells (same BoxArray/DistributionMapping as `ba`, one component; read on Known cells only; null = 0).
    // On return phi holds g on Known cells and 0 on Solid cells.
    amrex::MultiFab const* known_value = nullptr;
    // Driver-supplied component ids: the common layer labels components itself (flood fill over gas-gas faces); non-null is NotBuilt.
    amrex::iMultiFab const* component_id = nullptr;

    // Mean removal of the right-hand side of singular components (D-067): Volume (default) or the FDS parity switch
    // ScaledArithmetic. See the header comment at the top of this file.
    MeanKind mean_kind = MeanKind::Volume;

    // Optional gauge inputs for singular components (FDS SYMM_INDEFINITE gauge, pres.f90 ~1830-1878): after the
    // solve, phi -= sum(w*(phi - g)) / sum(w) with w = cell volume * gauge_weight and g = gauge_offset. FDS passes
    // gauge_weight = RHOP (RHO predictor, RHOS corrector) and gauge_offset = KRES. Null = weight 1 / offset 0,
    // which is the plain volume-weighted mean removal. Same BoxArray/DistributionMapping as `ba`, one component.
    // Single-level API only; for a composite hierarchy use PressureLevel::gauge_weight/gauge_offset (these two
    // fields set together with `levels` are NotBuilt, to avoid an ambiguous layout).
    amrex::MultiFab const* gauge_weight = nullptr;
    amrex::MultiFab const* gauge_offset = nullptr;

    amrex::MultiFab const* rhs = nullptr;         // one component, valid cells only are read
    amrex::MultiFab* phi = nullptr;               // out: H on valid cells (gauge-fixed); in: initial guess if used
};

// FR-039 trigger points for the expensive checks. The caller (driver) says what kind of solve this is by OR-ing flags into
// PressureOptions::trigger; a PressureWorkspace adds FirstSolve (its first solve) and FirstAfterRegrid (first solve after
// rebuild(), explicit or automatic) by itself. If any flag of the solve is in PressureOptions::full_checks_on the solve runs
// the FULL checks, otherwise only the cheap path.
//   Full checks: structural validation of the layout, component bookkeeping (pin), compatibility diagnostic (rms of the
//                right-hand side, removed-mean warning), true residual b - L H (if check_residual).
//   Cheap path:  selector, layout consistency of rhs/phi, mean removal, backend, gauge, convergence status. The solution is
//                bitwise the same as on the full path; removed_mean and gauge_shift are still reported.
enum SolveTrigger : unsigned { SolveRoutine = 0u, SolveFirst = 1u, SolveFirstAfterRegrid = 2u, SolveDebug = 4u };

// HYPRE backend settings (ADR-002 backend (H)); defaults follow FDS (pres.f90 3020-3035, imkl.f90 237-263) except where noted.
// The tolerance and iteration cap are PressureOptions::tol_rel and max_iter.
//  krylov:         PCG (symmetric positive definite matrices: every single-level problem) or GMRES / BiCGSTAB for the composite
//                  matrix, which is not symmetric (the C/F interpolation is one-sided); Auto = PCG on one level, GMRES on a
//                  hierarchy.
//  coarsen_type:   BoomerAMG coarsening (FDS 8 = PMIS).     relax_type: smoother (FDS 18 = L1-Jacobi, symmetric).
//  num_sweeps:     sweeps per level.                        strong_threshold: BoomerAMG strength threshold, < 0 = HYPRE default.
//  interp_type:    interpolation, < 0 = HYPRE default.      agg_levels: aggressive coarsening levels.
// Tolerance: the Krylov solver is run to 0.1*tol_rel (its recurrence residual is slightly optimistic); `converged` is judged against
// tol_rel. Singular problems are pinned by an identity row (FDS ULMAT), which concentrates the residual of all other rows in
// the pin row (the residual of a conservative operator with a mean-free right-hand side sums to zero), so the raw relative
// residual of a large singular problem (about 1e6 unknowns) exceeds 1e-12 although the solution is accurate to 1e-13. The check
// therefore warns on the residual with the pin cell excluded, against the size-scaled limit max(residual_tol, residual_floor)
// (PressureResult::residual_check / residual_limit, frozen/hypre-notes.md, section "Residual check"); the full residual and the pin row are reported. Typical counts at tol_rel 1e-12: 20-27 PCG
// iterations single level, 18-28 GMRES iterations composite (MLMG: 9-12); a tuned set (coarsen 10, relax 6, strong 0.5, agg 1)
// needs 19-22 and is about 2x faster.
// Round-off factor c of the residual_floor that bounds the residual limit from below (FR-031, A-67): about the number of floating-point operations that enter one row of
// b - L H (7-point stencil with the boundary ghost, 9 flops), times the unit round-off 2^-53.
inline constexpr double kResidualRoundoff = 10.0;

enum class HypreKrylov { Auto, PCG, GMRES, BiCGSTAB };
struct HypreOptions {
    HypreKrylov krylov = HypreKrylov::Auto;
    int coarsen_type = 8;
    int relax_type = 18;
    int num_sweeps = 1;
    double strong_threshold = -1.0;
    int interp_type = -1;
    int agg_levels = 0;
    int gmres_kdim = 30;
    bool operator== (HypreOptions const& o) const
    {
        return krylov == o.krylov && coarsen_type == o.coarsen_type && relax_type == o.relax_type && num_sweeps == o.num_sweeps &&
               strong_threshold == o.strong_threshold && interp_type == o.interp_type && agg_levels == o.agg_levels && gmres_kdim == o.gmres_kdim;
    }
};

struct PressureOptions {
    BackendKind backend = BackendKind::Auto;      // Auto = selector rule; FFT/MLMG/HYPRE = request, still rule-checked
    HypreOptions hypre;                           // used when backend == HYPRE
    double tol_rel = 1.0e-12;                     // backend relative tolerance
    int max_iter = 200;                           // backend iteration cap
    bool use_initial_guess = false;               // false: start from zero
    bool remove_mean = true;                      // test hook: false disables mean removal (singular solves then warn)
    bool check_residual = true;                   // true composite residual b - L phi of the mean-removed problem
    double residual_tol = 1.0e-12;                // base of the residual limit: warn when the checked residual (pin cell excluded if a pin was applied) exceeds max(residual_tol, residual_floor)
    double removed_mean_warn = 1.0e-10;           // warn when removed mean (rel. to rms b) exceeds this
    // Trigger points (FR-039). Default: every solve is a Debug solve, i.e. full checks always (the behaviour before this option).
    // A driver sets SolveRoutine for ordinary solves (FDS solves up to 10 times per step) and passes SolveFirst /
    // SolveFirstAfterRegrid / SolveDebug when it wants the full checks.
    unsigned trigger = SolveDebug;
    unsigned full_checks_on = SolveFirst | SolveFirstAfterRegrid | SolveDebug;
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
    bool plan_reused = false;                     // FFT/HYPRE: the set-up cached in the PressureWorkspace was used (not built for this solve)
    // HYPRE only (0 otherwise): matrix size, set-up time (assembly + BoomerAMG setup; 0 when the cached set-up was reused),
    // solver time of this solve, Krylov method actually used ("PCG", "GMRES", "BiCGSTAB").
    long long hypre_rows = 0, hypre_nnz = 0;
    double hypre_setup_seconds = 0.0, hypre_solve_seconds = 0.0;
    std::string hypre_method;
    // The solve pinned one unknown by an identity row (HYPRE on a singular problem): level and cell. Used by the common layer's
    // residual check (the pin row is not an independent equation, see PressureResult::residual_check).
    bool pin_applied = false;
    int pin_level = 0;
    amrex::IntVect pin_cell{0};
};

struct PressureResult {
    Status status = Status::Ok;
    std::string message;                          // reason for any non-Ok status
    std::string backend;                          // "FFT", "MLMG", "HYPRE" (empty if not built)
    BackendStatus backend_status;
    std::vector<ComponentInfo> components;
    bool residual_checked = false;
    double residual_rel2 = 0.0;                   // ||b - L H||_2 / ||b||_2 (mean-removed b); composite: volume-weighted over uncovered cells
    double residual_relmax = 0.0;                 // max|b - L H| / max|b|; composite: over uncovered cells
    bool residual_ok = true;                      // residual_check <= residual_limit
    // Pin-aware residual check (definition: README "Residual check", frozen/hypre-notes.md). r = b - L H on the uncovered cells; the
    // norms below are scaled exactly like residual_rel2 / residual_relmax (weighted 2-norm / ||b||_2, max / max|b|; weights are the cell
    // volumes on a hierarchy). residual_rel2 and residual_relmax stay the FULL residual (pin row included) and are always reported.
    //   residual_rel2_nopin    ||r with the pin cell zeroed||_2 / ||b||_2       (= residual_rel2 when no pin was applied)
    //   residual_relmax_nopin  max over the non-pin cells of |r| / max|b|        (= residual_relmax when no pin was applied)
    //   residual_pin_rel       |r_pin| (weighted like the norm) / ||b||_2: the residual of the identity pin row, reported separately
    //   residual_pin_abs       |r_pin| unweighted, absolute
    //   residual_check         the value warned on: residual_rel2_nopin for a singular component, residual_rel2 otherwise
    //   residual_limit         max(residual_tol, residual_floor): the value residual_check is compared with (FR-031, A-67; equals residual_tol
    //                          while the round-off floor is below it, i.e. on small problems)
    // Extra (residual_floor enters the limit; the others are informational):
    //   residual_rel2_mr       ||r - mean(r)||_2 / ||b||_2, exact weighted mean of r (equals the raw value for singular components
    //                          with a mean-free b: the mean of r is round-off; equals residual_rel2 when not singular)
    //   residual_floor         10 * 2^-53 * ||A|| ||H||_2 / ||b||_2 with ||A|| = 4 sum_d 1/dx_d^2 (finest level, directions of more than
    //                          one cell): the round-off size of a floating-point residual; set for singular and non-singular components
    //                          alike and used as the lower bound of residual_limit
    //   residual_backward      ||r_nopin||_2 / (||b||_2 + ||A|| ||H||_2): normwise backward error
    double residual_rel2_nopin = 0.0;
    double residual_relmax_nopin = 0.0;
    double residual_pin_rel = 0.0;
    double residual_pin_abs = 0.0;
    double residual_check = 0.0;
    double residual_limit = 0.0;
    double residual_rel2_mr = 0.0;
    double residual_floor = 0.0;
    double residual_backward = 0.0;
    std::vector<std::string> warnings;
    unsigned triggers = 0;                        // effective trigger flags of this solve (options OR workspace-detected)
    bool full_checks = true;                      // the full checks ran (false: cheap path; residual_checked is then false, removed_rel is 0)
    bool workspace_rebuilt = false;               // composite: the workspace passed in did not match and was rebuilt
    int nlevels = 1;                              // composite: number of levels solved
    amrex::Long ncells_uncovered = 0;             // composite: number of uncovered cells over all levels
};

// Selector: FFT only on a single level with no masked cells and a uniformly closed (Neumann/periodic) or
// uniformly open (Dirichlet) set of domain faces; everything else is NotBuilt. Collective.
struct Selection { bool ok = false; BackendKind kind = BackendKind::Auto; std::string message; };
Selection select_backend (PressureProblem const& p, BackendKind requested);

// Per-layout solver set-up that outlives one solve (D-058). Two independent parts, both keyed by the layout and discarded by
// rebuild():
//  - composite path: the covered-cell masks, the cell volumes and the MLMG hierarchy (MLPoisson with its multigrid levels and
//    bottom solver), for a problem with `levels`;
//  - single-level FFT path: the FFT::Poisson plan and its spectral work arrays, keyed by BoxArray, DistributionMapping,
//    Geometry and the six boundary types of the single-level fields (FDS solves up to ten times per step on one layout).
// Build it once per layout and pass it to solve_pressure() for every solve until the layout (BoxArrays, DistributionMappings,
// geometries, ratios or BC) changes; then call rebuild(). A workspace stores copies of the layout only: no pointer or
// reference into any rhs/phi MultiFab of the problem survives a call, so the old MultiFabs may be destroyed at once after a
// regrid. A solve given a workspace that does not match the problem rebuilds the stale part automatically
// (PressureResult::workspace_rebuilt = true), so a stale workspace is never used. Results are bitwise the same with or
// without a workspace. The workspace also counts solves and pending regrids, which supply the FR-039 trigger flags
// SolveFirst and SolveFirstAfterRegrid (see PressureOptions::trigger).
class PressureWorkspace {
public:
    PressureWorkspace ();
    ~PressureWorkspace ();
    PressureWorkspace (PressureWorkspace&&) noexcept;
    PressureWorkspace& operator= (PressureWorkspace&&) noexcept;
    PressureWorkspace (PressureWorkspace const&) = delete;
    PressureWorkspace& operator= (PressureWorkspace const&) = delete;

    // Discard the old setup (composite and FFT) and build it for `p`: `p.levels` non-empty = composite hierarchy, otherwise
    // the single-level fields (the FFT plan is built now if the selector picks FFT; MLMG single level has nothing to cache).
    // Collective. Returns Ok, or NotBuilt/InvalidInput with `message` (the same checks as the solve). The workspace is empty
    // after a non-Ok return. Marks a regrid as pending (FirstAfterRegrid for the next solve).
    Status rebuild (PressureProblem const& p, std::string* message = nullptr);
    bool built () const;                          // composite set-up or FFT plan present
    int num_levels () const;                      // composite levels, 0 when only an FFT plan (or nothing) is built
    // True when the workspace was built for the layout of `p` (composite: ba, dm, geom, ratios, bc; single level: the FFT
    // plan key).
    bool matches (PressureProblem const& p) const;

    // FR-039 bookkeeping. Counts completed solves that used this workspace; regrid_pending() is true after rebuild() until the
    // next solve completes; auto_triggers() = the flags this workspace adds to the next solve.
    unsigned long num_solves () const { return m_nsolves; }
    bool regrid_pending () const { return m_regrid; }
    unsigned auto_triggers () const { return (m_nsolves == 0 ? unsigned(SolveFirst) : 0u) | (m_regrid ? unsigned(SolveFirstAfterRegrid) : 0u); }
    // FFT plan counters (single-level path): plans constructed, solves that reused the cached plan.
    int fft_plan_builds () const { return m_fft_builds; }
    int fft_plan_reuses () const { return m_fft_reuses; }
    bool fft_plan_built () const;
    bool hypre_built () const;                    // a HYPRE set-up (matrix + BoomerAMG) is cached for the single-level path

    // Internal (used by solve_pressure).
    void note_solve_done () { ++m_nsolves; m_regrid = false; }
    void note_fft (bool reused) { if (reused) { ++m_fft_reuses; } else { ++m_fft_builds; } }
    void mark_regrid () { m_regrid = true; }
    void discard ();                              // drop all set-up (counters are kept)
    struct Impl;
    Impl* impl () { return m_impl.get(); }
    Impl& ensure_impl ();
private:
    std::unique_ptr<Impl> m_impl;
    unsigned long m_nsolves = 0;
    bool m_regrid = false;
    int m_fft_builds = 0, m_fft_reuses = 0;
};

// Full solve: selector, components and pins, mean removal, backend, gauge, checks. Collective. On any
// status other than Ok/NotConverged the phi field is not modified. `ws` holds the composite hierarchy set-up and the FFT plan
// (null = build the set-up for this call and discard it after the call). Which checks run: see SolveTrigger.
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

// ---- Inhomogeneous boundary data: fold into the right-hand side (M2 item d) ------------------------------------------------
// The solve is for the HOMOGENEOUS boundary problem (zero-gradient / zero-value / periodic). Non-zero boundary data enter by moving
// the known ghost-cell terms to the right-hand side, once, before solve_pressure.
//
// Convention (FDS pres.f90: header at lines 57-58, boundary application at 452-454, ULMAT F_H terms 1620-1679):
//   Neumann face, data g = dH/dx_d, the derivative along the INCREASING coordinate at both the low and the high face
//     (so g is the outward derivative at a high face and minus the outward derivative at a low face);
//   Dirichlet face, data H_b = the value of H at the wall (not at the ghost cell).
// Ghost values used by the operator L (CommonLayer.H) with data, h = cell width in direction d, phi_1 / phi_n the first / last
// interior cell of the line:
//   low Neumann    phi_g = phi_1 - h*g       high Neumann  phi_g = phi_n + h*g       Dirichlet (either side)  phi_g = 2*H_b - phi_1
// Moving the data terms to the right-hand side (b is the right-hand side of L_hom phi = b, L_hom the homogeneous operator):
//   low Neumann:   b += g / h          high Neumann:  b -= g / h          Dirichlet (either side):  b -= 2*H_b / h^2
// applied to the boundary-adjacent layer of cells of the face only (corner and edge cells collect one term per face). No other factor
// (no cell volume): b is the unit-operator right-hand side of this interface (volume-scaled F = V*b is formed by the caller if needed).
// Periodic faces take no data. A one-cell y direction drops its term (D-057), so data on its faces are ignored.
// A Dirichlet face in a one-cell x or z direction is NotBuilt, as in the solve.
// Data per face: `value[f]` (cell-centred MultiFab, one component, any BoxArray / DistributionMapping; it is read in the boundary-adjacent
// layer of cells, cell indices as the domain, cells of the layer that it does not cover are read as 0) or, if null, the uniform `constant[f]`.
struct BoundaryData {
    std::array<amrex::MultiFab const*,6> value{};     // index face_index(dir, side)
    std::array<amrex::Real,6> constant{};
};
// Adds the fold terms to `rhs` (BoxArray / DistributionMapping of p.ba / p.dm, valid cells). Collective. Single-level problems only.
// InvalidInput: data on a periodic face, rhs layout different from the problem, bad data MultiFab. rhs untouched on error.
Status fold_boundary_data (PressureProblem const& p, BoundaryData const& bd, amrex::MultiFab& rhs, std::string* message = nullptr);

const char* to_string (Status s);
const char* to_string (BC b);

} // namespace pb
#endif
