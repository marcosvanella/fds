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
#include <AMReX_iMultiFab.H>

#include <vector>

namespace fdsamr {

// Sums of `ngroups` groups of terms. Local terms: group[i] (0-based) and val[i]. Returns the global exact sums (the same on every rank).
std::vector<double> exact_group_sums(int ngroups, const std::vector<int>& group, const std::vector<double>& val);

// Global exact sum of comp `comp` of the valid cells of `mf`, each cell weighted by `weight` (constant cell volume).
double exact_sum(const amrex::MultiFab& mf, int comp, double weight);

// Global exact sum of (a(i,j,k,ca) * b(i,j,k,cb)) * weight over valid cells (for example mass = RHO * Z * cell volume).
double exact_sum_product(const amrex::MultiFab& a, int ca, const amrex::MultiFab& b, int cb, double weight);

// Uncovered-cell variants (S9, D-028 under AMR, Role 3 plan section 4 item 2): `covered` is the fine-mask of the level (same BoxArray/DistributionMapping as the
// data, any ghost width; nonzero = the valid cell is covered by a finer level, see LevelRegistry::covered_mask). Covered cells contribute nothing, so a
// level-by-level sum of the uncovered cells counts every point of the domain once. nullptr = nothing covered: bitwise equal to exact_sum / exact_sum_product.
// The result still depends on the multiset of terms only (not on box split or rank count), the multiset being the uncovered cells.
double exact_sum_uncovered(const amrex::MultiFab& mf, int comp, double weight, const amrex::iMultiFab* covered);
double exact_sum_product_uncovered(const amrex::MultiFab& a, int ca, const amrex::MultiFab& b, int cb, double weight, const amrex::iMultiFab* covered);

// Several levels in one fixed-point sum (one scale for the whole hierarchy, so the total is independent of how the cells are split between levels and ranks):
// level l carries its own cell weight (volume of a level-l cell). Vectors have one entry per level; covered[l] may be nullptr (finest level).
double exact_sum_hierarchy(const std::vector<const amrex::MultiFab*>& mf, int comp, const std::vector<double>& weight, const std::vector<const amrex::iMultiFab*>& covered);

}  // namespace fdsamr

#endif
