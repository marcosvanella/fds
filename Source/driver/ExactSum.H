// ExactSum.H: decomposition-independent (exact) fixed-point sums for the driver (D-028, FR-005 (v)). M2a, S5.
//
// Kernel-facing rules (M2a): (a) passive scalars are not involved (scalar sums are taken one component at a time); (b) only uniform Cartesian metrics are used (the
// cell volume terms are the products DX*RC*DY*DZ that the callers build, RC=1).
//
// Method. Every term is converted to an integer multiple of 2^-s, with ONE scale s per group that is derived from the global maximum |term| of the group
// (max-reduction: exact, order independent), so that every integer is below 2^62. The integers are added in 128-bit arithmetic (exact, associative), the rank
// partial sums are combined limb-wise with an integer sum reduction, and the result is converted back once (one rounding). The result therefore depends
// on the multiset of terms only: not on the order, the box split, the tiling, or the rank count. The conversion of a term rounds it to 2^-(62-e) with
// e = exponent of the group maximum: terms smaller than 2^-62 of the maximum are kept to that absolute resolution (documented, deterministic).
//
// This file is the driver's own implementation; it does not include Source/pressure_backend/ExactSum.H (Role 2's private header).
#ifndef FDS_EXACTSUM_H_
#define FDS_EXACTSUM_H_

#include <AMReX_MultiFab.H>

#include <vector>

namespace fdsamr {

// Sums of `ngroups` groups of terms. Local terms: group[i] (0-based) and val[i]. Returns the global exact sums (the same on every rank).
std::vector<double> exact_group_sums(int ngroups, const std::vector<int>& group, const std::vector<double>& val);

// Global exact sum of comp `comp` of the valid cells of `mf`, each cell weighted by `weight` (constant cell volume).
double exact_sum(const amrex::MultiFab& mf, int comp, double weight);

// Global exact sum of (a(i,j,k,ca) * b(i,j,k,cb)) * weight over valid cells (for example mass = RHO * Z * cell volume).
double exact_sum_product(const amrex::MultiFab& a, int ca, const amrex::MultiFab& b, int cb, double weight);

}  // namespace fdsamr

#endif
