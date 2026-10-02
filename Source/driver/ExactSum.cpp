// ExactSum.cpp: see ExactSum.H. Kernel-facing rules (M2a): (a) passive scalars are summed one component at a time by the callers; (b) only uniform Cartesian
// metrics are used.
#include "ExactSum.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace fdsamr {

namespace {
using i128 = __int128;

// limbs of a signed 128-bit value: v = l0 + l1 2^32 + l2 2^64 + l3 2^96, l0..l2 in [0,2^32), l3 signed
inline void split(i128 v, long l[4])
{
    l[0] = static_cast<long>(static_cast<std::uint64_t>(v) & 0xffffffffULL);
    l[1] = static_cast<long>(static_cast<std::uint64_t>(v >> 32) & 0xffffffffULL);
    l[2] = static_cast<long>(static_cast<std::uint64_t>(v >> 64) & 0xffffffffULL);
    l[3] = static_cast<long>(v >> 96);
}
inline i128 join(const long l[4])
{
    return static_cast<i128>(l[0]) + (static_cast<i128>(l[1]) << 32) + (static_cast<i128>(l[2]) << 64) + (static_cast<i128>(l[3]) << 96);
}
}  // namespace

std::vector<double> exact_group_sums(int ngroups, const std::vector<int>& group, const std::vector<double>& val)
{
    std::vector<double> mx(ngroups, 0.0);
    for (std::size_t i = 0; i < val.size(); ++i) mx[group[i]] = std::max(mx[group[i]], std::abs(val[i]));
    if (ngroups > 0) amrex::ParallelAllReduce::Max(mx.data(), ngroups, amrex::ParallelContext::CommunicatorSub());
    std::vector<int> sc(ngroups, 0);
    for (int g = 0; g < ngroups; ++g) {
        int e = 0;
        if (mx[g] > 0.0) std::frexp(mx[g], &e);   // mx = m 2^e, m in [0.5,1)
        sc[g] = 62 - e;                            // |x| 2^sc < 2^62
    }
    std::vector<i128> acc(ngroups, 0);
    for (std::size_t i = 0; i < val.size(); ++i) {
        const int g = group[i];
        acc[g] += static_cast<i128>(std::llround(std::ldexp(val[i], sc[g])));
    }
    std::vector<long> limbs(4 * static_cast<std::size_t>(ngroups));
    for (int g = 0; g < ngroups; ++g) split(acc[g], &limbs[4 * static_cast<std::size_t>(g)]);
    if (ngroups > 0) amrex::ParallelAllReduce::Sum(limbs.data(), 4 * ngroups, amrex::ParallelContext::CommunicatorSub());
    std::vector<double> out(ngroups, 0.0);
    for (int g = 0; g < ngroups; ++g) {
        const i128 v = join(&limbs[4 * static_cast<std::size_t>(g)]);
        out[g] = std::ldexp(static_cast<double>(v), -sc[g]);
    }
    return out;
}

double exact_sum(const amrex::MultiFab& mf, int comp, double weight)
{
    std::vector<int> g;
    std::vector<double> v;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.const_array(mfi);
        const amrex::Box b = mfi.validbox();
        for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k)
            for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j)
                for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) { g.push_back(0); v.push_back(a(i, j, k, comp) * weight); }
    }
    return exact_group_sums(1, g, v)[0];
}

double exact_sum_product(const amrex::MultiFab& ma, int ca, const amrex::MultiFab& mb, int cb, double weight)
{
    std::vector<int> g;
    std::vector<double> v;
    for (amrex::MFIter mfi(ma); mfi.isValid(); ++mfi) {
        const auto a = ma.const_array(mfi);
        const auto b2 = mb.const_array(mfi);
        const amrex::Box b = mfi.validbox();
        for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k)
            for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j)
                for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) { g.push_back(0); v.push_back(a(i, j, k, ca) * b2(i, j, k, cb) * weight); }
    }
    return exact_group_sums(1, g, v)[0];
}

namespace {
inline bool is_covered(const amrex::iMultiFab* cov, amrex::MFIter& mfi, int i, int j, int k)
{
    return cov != nullptr && (*cov)[mfi](amrex::IntVect(i, j, k)) != 0;
}
}  // namespace

double exact_sum_uncovered(const amrex::MultiFab& mf, int comp, double weight, const amrex::iMultiFab* covered)
{
    std::vector<int> g;
    std::vector<double> v;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const auto a = mf.const_array(mfi);
        const amrex::Box b = mfi.validbox();
        for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k)
            for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j)
                for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i)
                    if (!is_covered(covered, mfi, i, j, k)) { g.push_back(0); v.push_back(a(i, j, k, comp) * weight); }
    }
    return exact_group_sums(1, g, v)[0];
}

double exact_sum_product_uncovered(const amrex::MultiFab& ma, int ca, const amrex::MultiFab& mb, int cb, double weight, const amrex::iMultiFab* covered)
{
    std::vector<int> g;
    std::vector<double> v;
    for (amrex::MFIter mfi(ma); mfi.isValid(); ++mfi) {
        const auto a = ma.const_array(mfi);
        const auto b2 = mb.const_array(mfi);
        const amrex::Box b = mfi.validbox();
        for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k)
            for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j)
                for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i)
                    if (!is_covered(covered, mfi, i, j, k)) { g.push_back(0); v.push_back(a(i, j, k, ca) * b2(i, j, k, cb) * weight); }
    }
    return exact_group_sums(1, g, v)[0];
}

double exact_sum_hierarchy(const std::vector<const amrex::MultiFab*>& mf, int comp, const std::vector<double>& weight, const std::vector<const amrex::iMultiFab*>& covered)
{
    AMREX_ALWAYS_ASSERT(mf.size() == weight.size() && mf.size() == covered.size());
    std::vector<int> g;
    std::vector<double> v;
    for (std::size_t l = 0; l < mf.size(); ++l)
        for (amrex::MFIter mfi(*mf[l]); mfi.isValid(); ++mfi) {
            const auto a = mf[l]->const_array(mfi);
            const amrex::Box b = mfi.validbox();
            for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k)
                for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j)
                    for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i)
                        if (!is_covered(covered[l], mfi, i, j, k)) { g.push_back(0); v.push_back(a(i, j, k, comp) * weight[l]); }
        }
    return exact_group_sums(1, g, v)[0];
}

}  // namespace fdsamr
