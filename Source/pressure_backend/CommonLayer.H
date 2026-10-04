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
// (lowest global index cell, x fastest) per component. Flood fill over beta != 0 faces is M2.
ComponentMap label_components (PressureProblem const& p);

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
void remove_mean (amrex::MultiFab& rhs, ComponentMap& cm, amrex::iMultiFab const* uncovered,
                  amrex::Real cell_volume_uniform, amrex::MultiFab const* cell_volume = nullptr,
                  MeanKind kind = MeanKind::Volume);

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

} // namespace pb
#endif
