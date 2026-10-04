// ExactSum.H: private, mask-aware, volume-weighted exact fixed-point sum (D-028 rule applied to the
// pressure path). Every term w*x is rounded to a 128-bit integer on a common power-of-two grid; integer
// addition is associative, so the result is bitwise independent of box layout, rank count, thread count
// and summation order. The grid is set by the global maximum |w*x| (an order-independent reduction), so
// each term is kept to 2^-88 of that maximum; up to 2^38 selected cells cannot overflow the accumulator.
//
// Weights: each term is `weight * w(i,j,k) * x(i,j,k)`, where `w` is the optional per-cell weight field `wfield`
// (same BoxArray and DistributionMapping as `mf`, first component; nullptr = 1). Per-cell volumes of a
// non-uniform mesh, or rho*volume for the gauge, are passed this way. Without `wfield` the term is
// `weight * x`, exactly as before.
//
// Cell selection: `uncovered` (nonzero = counted; nullptr = all cells) and `label` (component id >= 0;
// negative = excluded; nullptr = component 0). Both must share the BoxArray and DistributionMapping of `mf`.
#ifndef PB_EXACT_SUM_H_
#define PB_EXACT_SUM_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelContext.H>
#include <cmath>
#include <cstdint>
#include <vector>

namespace pb {

struct ExactSumResult {
    std::vector<double> sum;        // sum over selected cells of weight*x, per component
    std::vector<amrex::Long> count; // number of selected cells, per component
    bool finite = true;             // false if any selected term was NaN or Inf (sum is then NaN)
};

namespace detail {
using i128 = __int128;
constexpr int kFracBits = 88;

inline i128 to_i128 (double v) { return static_cast<i128>(v); }   // v is integer valued and |v| <= 2^88
inline double from_i128 (i128 a, int e)
{
    const auto lo = static_cast<std::uint64_t>(static_cast<unsigned __int128>(a));
    const auto hi = static_cast<std::int64_t>(a >> 64);
    long double v = static_cast<long double>(hi) * 18446744073709551616.0L + static_cast<long double>(lo);
    return static_cast<double>(std::ldexp(v, e - kFracBits));
}
}

inline ExactSumResult exact_sum (amrex::MultiFab const& mf, int comp, double weight, int ncomp,
                                 amrex::iMultiFab const* uncovered = nullptr,
                                 amrex::iMultiFab const* label = nullptr,
                                 amrex::MultiFab const* wfield = nullptr)
{
    static_assert(sizeof(amrex::Real) == 8, "ExactSum assumes double precision");
    using namespace amrex;
    auto selected = [&] (Array4<int const> const& u, Array4<int const> const& l, int i, int j, int k, int& c) {
        if (u.contains(i,j,k) && u(i,j,k) == 0) { return false; }
        c = l.contains(i,j,k) ? l(i,j,k) : 0;
        return c >= 0 && c < ncomp;
    };
    // Pass 1: global max |w*x| and finiteness (both order independent).
    double red[2] = {0.0, 0.0};
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.const_array(mfi, comp);
        Array4<int const> u, l;
        Array4<Real const> wf;
        if (uncovered) { u = uncovered->const_array(mfi); }
        if (label) { l = label->const_array(mfi); }
        if (wfield) { wf = wfield->const_array(mfi); }
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            int c;
            if (!selected(u, l, i, j, k, c)) { return; }
            const double t = wfield ? weight * wf(i,j,k) * a(i,j,k) : weight * a(i,j,k);
            if (!std::isfinite(t)) { red[1] = 1.0; } else { red[0] = std::max(red[0], std::abs(t)); }
        });
    }
    ParallelDescriptor::ReduceRealMax(red, 2);
    ExactSumResult R;
    R.sum.assign(ncomp, 0.0);
    R.count.assign(ncomp, 0);
    R.finite = (red[1] == 0.0);
    int e = 0;
    if (red[0] > 0.0) { std::frexp(red[0], &e); }   // |t| < 2^e
    std::vector<detail::i128> acc(ncomp, 0);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.const_array(mfi, comp);
        Array4<int const> u, l;
        Array4<Real const> wf;
        if (uncovered) { u = uncovered->const_array(mfi); }
        if (label) { l = label->const_array(mfi); }
        if (wfield) { wf = wfield->const_array(mfi); }
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            int c;
            if (!selected(u, l, i, j, k, c)) { return; }
            R.count[c] += 1;
            const double t = wfield ? weight * wf(i,j,k) * a(i,j,k) : weight * a(i,j,k);
            if (std::isfinite(t)) {
                acc[c] += detail::to_i128(std::nearbyint(std::ldexp(t, detail::kFracBits - e)));
            }
        });
    }
    // Integer reduction across ranks (order independent).
    std::vector<std::int64_t> mine(3*ncomp);
    for (int c = 0; c < ncomp; ++c) {
        mine[3*c]   = static_cast<std::int64_t>(static_cast<std::uint64_t>(static_cast<unsigned __int128>(acc[c])));
        mine[3*c+1] = static_cast<std::int64_t>(acc[c] >> 64);
        mine[3*c+2] = R.count[c];
    }
    const int np = ParallelContext::NProcsSub();
    std::vector<std::int64_t> all(mine.size() * np);
#ifdef AMREX_USE_MPI
    MPI_Allgather(mine.data(), static_cast<int>(mine.size()*sizeof(std::int64_t)), MPI_BYTE,
                  all.data(), static_cast<int>(mine.size()*sizeof(std::int64_t)), MPI_BYTE,
                  ParallelContext::CommunicatorSub());
#else
    all = mine;
#endif
    for (int c = 0; c < ncomp; ++c) {
        detail::i128 s = 0; Long n = 0;
        for (int r = 0; r < np; ++r) {
            const std::int64_t* p = &all[static_cast<std::size_t>(r)*mine.size() + 3*c];
            s += (static_cast<detail::i128>(p[1]) << 64) + static_cast<detail::i128>(static_cast<std::uint64_t>(p[0]));
            n += p[2];
        }
        R.count[c] = n;
        R.sum[c] = R.finite ? detail::from_i128(s, e) : std::nan("");
    }
    return R;
}

} // namespace pb
#endif
