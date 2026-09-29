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

// Volume-weighted (cell volume) exact-sum mean removal per singular component, over uncovered cells.
// Idempotent: a mean at or below 2^-52 max|rhs| is left alone, so a second application is a bitwise no-op.
// Fills removed_mean and removed_rel of the components. Bitwise independent of decomposition.
void remove_mean (amrex::MultiFab& rhs, ComponentMap& cm, amrex::iMultiFab const* uncovered,
                  amrex::Real cell_volume);

// Gauge after the solve: subtract the exact volume-weighted mean per singular component.
void apply_gauge (amrex::MultiFab& phi, ComponentMap& cm, amrex::iMultiFab const* uncovered,
                  amrex::Real cell_volume);

// out = L phi on valid cells (reference implementation). phi needs >= 1 ghost layer; its ghost cells are
// overwritten. out has the same BoxArray and DistributionMapping.
void apply_operator (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab& out);

// Relative 2-norm and max-norm of b - L phi.
struct ResidualNorms { double rel2 = 0.0; double relmax = 0.0; };
ResidualNorms true_residual (PressureProblem const& p, amrex::MultiFab& phi, amrex::MultiFab const& rhs);

} // namespace pb
#endif
