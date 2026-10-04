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

namespace pb {

struct ComponentMap {
    amrex::iMultiFab label;                       // component id per valid cell (>= 0)
    std::vector<ComponentInfo> comps;
};

// M1: the whole box is one component (no masks); singular unless some face is Dirichlet. Records the pin
// (lowest global index cell, x fastest) per component (skipped with record_pin = false: cheap path). Flood fill over beta != 0 faces is M2.
ComponentMap label_components (PressureProblem const& p, bool record_pin = true);

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

// out = L phi on valid cells (reference implementation). phi needs >= 1 ghost layer; its ghost cells are
// overwritten. out has the same BoxArray and DistributionMapping.
void apply_operator (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab& out);

// Relative 2-norm and max-norm of b - L phi.
struct ResidualNorms { double rel2 = 0.0; double relmax = 0.0; };
ResidualNorms true_residual (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab const& rhs);

// ---- Scaled, pin-aware residual check (see PressureResult::residual_check and frozen/hypre-notes.md) -----------------------
// The sums behind it, all with the same cell weights w (1 on a single level, the cell volume on a hierarchy) over the uncovered cells:
//   r2 = sum w r^2, r2_nopin = same with the pin cell's r set to 0, pin2 = w_pin r_pin^2, sumr = sum w r, sumw = sum w,
//   b2 = sum w b^2, phi2 = sum w H^2, anorm = bound of the operator norm (4 sum_d 1/dx_d^2 over directions with more than one cell,
//   finest level), rmax = max |r| over all cells, rmax_nopin = same without the pin cell, bmax = max |b|, pin_abs = |r_pin|.
struct ResidualSums {
    double r2 = 0, r2_nopin = 0, pin2 = 0, sumr = 0, sumw = 0, b2 = 0, phi2 = 0, anorm = 0, rmax = 0, rmax_nopin = 0, bmax = 0, pin_abs = 0;
    bool singular = false;
    bool pin_excluded = false;
    const char* what = "";                        // wording in the warning ("", "composite ")
};
// Fills the residual_* fields of R and appends the warning if the check fails. Non-singular: residual_check = residual_rel2 against
// residual_tol, exactly as before. Singular: residual_check = residual_rel2_nopin (pin cell excluded when a pin was applied) against
// the same residual_tol; residual_floor / residual_rel2_mr / residual_backward are informational.
void evaluate_residual (PressureResult& R, PressureOptions const& o, ResidualSums const& S);

// Single-level: the sums for phi (gauge-fixed H) and the mean-removed right-hand side rhs; `pin` is the pinned cell of the solve
// or nullptr. Collective. phi needs one ghost layer (overwritten).
ResidualSums residual_sums (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab const& rhs, bool singular,
                            amrex::IntVect const* pin);

} // namespace pb
#endif
