// CommonLayer.H: parts applied to every backend, from outside (mean removal, components and pins, gauge,
// reference operator and true-residual check).
//
// Discrete operator (defined here, not by a backend): for each valid cell
//   (L phi) = sum_d ( phi(+e_d) + phi(-e_d) - 2 phi ) / dx_d^2 ,
// with the out-of-domain neighbour of a non-periodic face taken as phi (Neumann, zero gradient) or -phi
// (homogeneous Dirichlet: ghost 2*phi_b - phi_1 with phi_b = 0).
#ifndef PB_COMMON_LAYER_H_
#define PB_COMMON_LAYER_H_

#include "PressureIface.H"

#include <cstdint>
#include <memory>

namespace pb {

struct ComponentMap {
    amrex::iMultiFab label;                       // component id per valid cell (>= 0)
    std::vector<ComponentInfo> comps;
    int rounds = 0;                               // masked flood fill: ghost-exchange rounds (incl. the final no-change round); 0 = not a masked fill
    double seconds = 0.0;                         // masked flood fill: wall time of the labelling (max over ranks)
};

// ---- Masked branch (frozen/masked-notes.md) ---------------------------------------------------------------------------------
// True when p.cell_class is set and has a non-zero cell (collective). An all-zero or null cell_class is the unmasked problem.
bool is_masked (PressureProblem const& p);
// cell_class / known_value with one ghost layer: periodic images are filled; cells outside a non-periodic domain face are
// CellSolid (class) and 0 (value). Collective. known_ghost is all zero when p.known_value is null.
amrex::iMultiFab class_ghost (PressureProblem const& p);
amrex::MultiFab known_ghost (PressureProblem const& p);

// Content hash of p.cell_class (valid cells; 0 for a null cell_class). Collective. Built from per-box FNV-1a hashes of (box, data) that are
// summed over ranks in exact integer arithmetic, so it does not depend on which rank owns a box. Key of the workspace's label cache.
std::uint64_t class_content_hash (PressureProblem const& p);

// M1: the whole box is one component (no masks); singular unless some face is Dirichlet. Records the pin
// (lowest global index cell, x fastest) per component (skipped with record_pin = false: cheap path). With a masked problem
// (is_masked) the components are the flood-filled gas regions (gas-gas faces, periodic wrap included): label = component id on gas
// cells and -1 on Solid and Known cells; ids are ordered by the lowest global index of the component, which is its pin cell; a
// component is singular unless a gas cell of it touches a domain Dirichlet face or a Known cell. The pin is always recorded.
ComponentMap label_components (PressureProblem const& p, bool record_pin = true);
// Composite masked components (D4).
struct CompositePin { int level = 0; amrex::IntVect cell{0}; };
struct CompositeComponentMap {
    std::vector<std::unique_ptr<amrex::iMultiFab>> label;
    std::vector<ComponentInfo> comps;
    std::vector<CompositePin> pins;   // one per singular component, ComponentInfo::pin mirrored on comps
    int rounds = 0;
    double seconds = 0.0;
};
bool composite_is_masked (PressureProblem const& p);
// W_unc[l] = uncovered mask (1 = uncovered). Collective. Non-empty return = InvalidInput message (e.g. Known on C/F).
std::string label_components_composite_masked (PressureProblem const& p,
    std::vector<std::unique_ptr<amrex::iMultiFab>> const& W_unc,
    std::vector<amrex::Geometry> const& geom,
    std::vector<amrex::IntVect> const& ratio,
    CompositeComponentMap& cm);

// D2: zone id (PressureProblem::component_id) of every component; a component whose gas cells carry two zone ids is an error (returned
// message; empty = ok). Sets ComponentInfo::zone. Collective.
std::string assign_zones (amrex::iMultiFab const& zone, ComponentMap& cm);

// Which zero mode of the right-hand side is removed when cell volumes vary: `MeanKind` (PressureIface.H). The solved
// system is the unit operator L (cell-centred 7-point, uniform cells), or, for a future non-uniform mesh, the
// volume-scaled finite-volume system K x = F with F = v*b (rows scaled by the cell volume v, entries area/distance,
// as FDS ULMAT/UGLMAT build it, pres.f90 ULMAT_SOLVE_ZONE and GET_FH_FROM_PRHS_AND_BCS). Its compatibility condition
// is the plain sum of F being zero.
//   Volume (default, D-067):  b -= sum(v*b)/sum(v)   (constant shift of b; F changes by c*v)
//   ScaledArithmetic (FDS):   b_k -= mean(F)/v_k      (arithmetic mean of F removed from F, as FDS does)
// Both make sum(F) zero. They give the same system when v is uniform and differ by a multiple of
// (v_k - mean v) otherwise; ScaledArithmetic is the FDS-parity switch (see frozen/mean-removal-vs-fds.md).

// Mean removal per singular component, over uncovered cells, with exact sums. Idempotent: a mean at or below
// 2^-52 of the largest term is left alone, so a second application is a bitwise no-op. Fills removed_mean and
// removed_rel of the components (for ScaledArithmetic they refer to F = v*b). Bitwise independent of
// decomposition.
//  - cell_volume == nullptr: uniform cells of volume `cell_volume_uniform`; both kinds are the same operation up to
//    rounding (volume-weighted mean of b, i.e. the arithmetic mean of v*b divided by v).
//  - cell_volume != nullptr: per-cell volumes (same BoxArray/DistributionMapping), `kind` selects the zero mode.
//  - measure = false skips the rms of the right-hand side: removed_rel is then 0 (cheap path); the removed constant, and so the
//    result, is bitwise the same.
void remove_mean (amrex::MultiFab& rhs, ComponentMap& cm, amrex::iMultiFab const* uncovered,
                  amrex::Real cell_volume_uniform, amrex::MultiFab const* cell_volume = nullptr,
                  MeanKind kind = MeanKind::Volume, bool measure = true);

// Gauge after the solve, per singular component: phi -= sum(W*(phi - g)) / sum(W) over uncovered cells, with
// W = v*rho. v is the cell volume (per-cell field or the uniform value), rho is `gauge_weight` (null = 1) and
// g is `gauge_offset` (null = 0). With both null this is the exact volume-weighted mean removed from phi.
// FDS: rho = RHOP (RHO in the predictor, RHOS in the corrector), g = KRES, H = phi (pres.f90 ~1830-1878).
// Components that are not singular are not shifted.
void apply_gauge (amrex::MultiFab& phi, ComponentMap& cm, amrex::iMultiFab const* uncovered,
                  amrex::Real cell_volume_uniform, amrex::MultiFab const* cell_volume = nullptr,
                  amrex::MultiFab const* gauge_weight = nullptr, amrex::MultiFab const* gauge_offset = nullptr);

// E2 (docs/pressure/08 section 4.2): out := phi minus the volume-weighted mean (exact sum, uncovered cells) of each SINGULAR component;
// non-singular components and cells with label < 0 are copied unchanged. Used on the solver's own solution (before the D-067 gauge) so the
// residual verdict and the round-off floor are those of the system actually solved and cannot depend on the gauge constant.
// `out` is (re)defined on phi's layout with one ghost layer. No side effect on cm.
void mean_removed_copy (amrex::MultiFab const& phi, ComponentMap const& cm, amrex::iMultiFab const* uncovered, amrex::Real cell_volume_uniform,
                        amrex::MultiFab& out);

// out = L phi on valid cells (reference implementation). phi needs >= 1 ghost layer; its ghost cells are
// overwritten. out has the same BoxArray and DistributionMapping.
void apply_operator (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab& out);

// Relative 2-norm and max-norm of b - L phi.
struct ResidualNorms { double rel2 = 0.0; double relmax = 0.0; };
ResidualNorms true_residual (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab const& rhs);

// Bound of the discrete operator norm on one level: 4 sum_d 1/dx_d^2 over directions with more than one cell (D-084: shared by
// residual_sums and the composite E3 floor). dx is CellSizeArray(); dom is the level Domain().
inline double operator_anorm (amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx, amrex::Box const& dom)
{
    double a = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (dom.length(d) > 1) { a += 4.0 / (double(dx[d]) * double(dx[d])); }
    }
    return a;
}

// ---- Scaled, pin-aware residual check (see PressureResult::residual_check and frozen/hypre-notes.md) -----------------------
// The sums behind it, all with the same cell weights w (1 on a single level, the cell volume on a hierarchy) over the uncovered cells:
//   r2 = sum w r^2, r2_nopin = same with every pin cell's r set to 0, pin2 = sum_pins w_pin r_pin^2, sumr = sum w r, sumw = sum w,
//   b2 = sum w b^2, phi2 = sum w H^2, anorm = bound of the operator norm (operator_anorm; single level: that grid; composite:
//   min(2*anorm_eff, A_finest) with anorm_eff level-wise A_l weighted like phi2, E3 / D-087), rmax = max |r| over all cells,
//   rmax_nopin = same without the pin cells, bmax = max |b|, pin_abs = max |r_pin| over pins.
struct ResidualSums {
    double r2 = 0, r2_nopin = 0, pin2 = 0, sumr = 0, sumw = 0, b2 = 0, phi2 = 0, anorm = 0, rmax = 0, rmax_nopin = 0, bmax = 0, pin_abs = 0;
    bool singular = false;
    bool pin_excluded = false;
    const char* what = "";                        // wording in the warning ("", "composite ")
};
// Fills the residual_* fields of R and appends the warning if the check fails. residual_check = residual_rel2_nopin for a singular
// component (the full residual when no pin was applied), residual_rel2 otherwise; it is compared with
// residual_limit = max(residual_tol, residual_floor), residual_floor = kResidualRoundoff 2^-53 anorm sqrt(phi2)/sqrt(b2) for every
// component, singular or not (FR-031, A-67). residual_rel2_mr / residual_backward are informational.
void evaluate_residual (PressureResult& R, PressureOptions const& o, ResidualSums const& S);

// Single-level: the sums for phi (gauge-fixed H) and the mean-removed right-hand side rhs; `pins` are the pinned cells of the
// solve (empty = none). Collective. phi needs one ghost layer (overwritten).
ResidualSums residual_sums (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab const& rhs, bool singular,
                            std::vector<amrex::IntVect> const& pins = {});

} // namespace pb
#endif
