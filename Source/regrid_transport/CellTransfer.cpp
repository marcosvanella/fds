// CellTransfer.cpp: see CellTransfer.H.
#include "CellTransfer.H"

#include <AMReX_BoxArray.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>

namespace fdsrt {

int clip_children(double* ch, int n, double parent, double floor)
{
    bool below = false;
    for (int m = 0; m < n; ++m) below = below || ch[m] < floor;
    if (!below) return 0;
    if (parent < floor) {
        for (int m = 0; m < n; ++m) ch[m] = parent;
        return 2;
    }
    double s = 0.0;
    for (int m = 0; m < n; ++m) {
        ch[m] = std::max(ch[m], floor) - floor;
        s += ch[m];
    }
    const double target = n * (parent - floor);
    const double f = s > 0.0 ? target / s : 0.0;
    for (int m = 0; m < n; ++m) ch[m] = floor + ch[m] * f;
    if (s <= 0.0) for (int m = 0; m < n; ++m) ch[m] = parent;   // all children were clipped to the floor: parent == floor (target 0), keep it
    return 1;
}

namespace {

// monotonised-central limited slope (change over one coarse cell) from the two one-sided differences
inline double mc_slope(double dm, double dp)
{
    if (dm * dp <= 0.0) return 0.0;
    const double a = std::abs(dm), b = std::abs(dp), c = 0.5 * std::abs(dm + dp);
    const double s = std::min(c, std::min(2.0 * a, 2.0 * b));
    return dm > 0.0 ? s : -s;
}

}  // namespace

ProlongStats prolong_conserved(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& crse_geom, const amrex::IntVect& ratio,
                               const ProlongOpts& opt, const amrex::MultiFab* old_fine)
{
    using namespace amrex;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(crse.nGrow() >= 1, "prolong_conserved: the coarse data need one valid ghost layer");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine.nComp() == crse.nComp(), "prolong_conserved: component counts differ");
    const int nc = fine.nComp();
    BoxArray cba = fine.boxArray();
    for (int i = 0; i < static_cast<int>(cba.size()); ++i) {
        const Box b = cba[i];
        for (int d = 0; d < 3; ++d)
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(b.smallEnd(d) % ratio[d] == 0 && (b.bigEnd(d) + 1) % ratio[d] == 0, "prolong_conserved: fine box is not a union of whole coarse cells");
    }
    cba.coarsen(ratio);
    MultiFab tmp(cba, fine.DistributionMap(), nc, 1);
    tmp.setVal(0.0);
    tmp.ParallelCopy(crse, 0, 0, nc, 1, 1, crse_geom.periodicity());   // valid and ghost cells of the coarse data, periodic images included

    int rd[3] = {ratio[0], ratio[1], ratio[2]};
    bool act[3] = {rd[0] > 1, rd[1] > 1, rd[2] > 1};
    double xmax[3], xoff[3][8];
    for (int d = 0; d < 3; ++d) {
        xmax[d] = 0.5 - 0.5 / rd[d];
        for (int m = 0; m < rd[d]; ++m) xoff[d][m] = (m + 0.5) / rd[d] - 0.5;
    }
    const int nchild = rd[0] * rd[1] * rd[2];
    ProlongStats st;
    std::vector<double> ch(nchild);
    for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
        const Box cb = tmp.box(mfi.index());   // valid coarse box = coarsened fine box
        auto q = tmp.const_array(mfi);
        auto f = fine.array(mfi);
        LoopOnCpu(cb, [&](int i, int j, int k) {
            ++st.parents;
            bool any_limited = false;
            for (int n = 0; n < nc; ++n) {
                const double qc = q(i, j, k, n);
                double s[3] = {0, 0, 0}, qmin = qc, qmax = qc;
                for (int d = 0; d < 3; ++d) {
                    if (!act[d]) continue;
                    const int e[3] = {d == 0, d == 1, d == 2};
                    const double qm = q(i - e[0], j - e[1], k - e[2], n), qp = q(i + e[0], j + e[1], k + e[2], n);
                    s[d] = mc_slope(qc - qm, qp - qc);
                    if (s[d] == 0.0 && (qc - qm) * (qp - qc) <= 0.0) any_limited = true;
                }
                // extremum bounds over the 3x3x3 neighbourhood (active directions only)
                for (int kk = -(act[2] ? 1 : 0); kk <= (act[2] ? 1 : 0); ++kk)
                    for (int jj = -(act[1] ? 1 : 0); jj <= (act[1] ? 1 : 0); ++jj)
                        for (int ii = -(act[0] ? 1 : 0); ii <= (act[0] ? 1 : 0); ++ii) {
                            const double v = q(i + ii, j + jj, k + kk, n);
                            qmin = std::min(qmin, v);
                            qmax = std::max(qmax, v);
                        }
                double lo = qmin;
                if (opt.use_floor) lo = std::max(lo, opt.floor);
                double dev = 0.0;
                for (int d = 0; d < 3; ++d) dev += std::abs(s[d]) * xmax[d];
                double alpha = 1.0;
                if (dev > 0.0) {
                    const double room = std::min(qc - lo, qmax - qc);
                    if (room < dev * (1.0 + 1e-12)) {
                        alpha = std::min(1.0, std::max(0.0, room) / dev) * (1.0 - 1e-12);   // tiny safety factor: the bound must hold after floating-point rounding
                        any_limited = true;
                    }
                }
                int m = 0;
                for (int kc = 0; kc < rd[2]; ++kc)
                    for (int jc = 0; jc < rd[1]; ++jc)
                        for (int ic = 0; ic < rd[0]; ++ic)
                            ch[m++] = qc + alpha * (s[0] * xoff[0][ic] + s[1] * xoff[1][jc] + s[2] * xoff[2][kc]);
                if (opt.use_floor) {
                    const int r = clip_children(ch.data(), nchild, qc, opt.floor);
                    if (r) { ++st.clips; if (r == 2) ++st.unfixable; }
                }
                m = 0;
                for (int kc = 0; kc < rd[2]; ++kc)
                    for (int jc = 0; jc < rd[1]; ++jc)
                        for (int ic = 0; ic < rd[0]; ++ic)
                            f(i * rd[0] + ic, j * rd[1] + jc, k * rd[2] + kc, n) = ch[m++];
            }
            if (any_limited) ++st.limited;
        });
    }
    if (old_fine) fine.ParallelCopy(*old_fine, 0, 0, nc, 0, 0);   // bitwise copy of the cells that already had fine data
    long v[4] = {st.parents, st.limited, st.clips, st.unfixable};
    ParallelAllReduce::Sum(v, 4, ParallelContext::CommunicatorSub());
    st.parents = v[0]; st.limited = v[1]; st.clips = v[2]; st.unfixable = v[3];
    return st;
}

}  // namespace fdsrt
