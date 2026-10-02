// test_celltransfer.cpp (R4 step 3): conservative limited prolongation of rho*Z_n (positivity-preserving, parent-sum clip) and the regrid budget.
// Usage: test_celltransfer (run on 1 and 4 ranks). Fine cell volume is 1/(r_x r_y r_z), coarse cell volume 1.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

#include "CellTransfer.H"
#include "LevelOps.H"
#include "TagOps.H"
#include "check.H"

using namespace fdsrt;
using Fn = std::function<double(int, int, int, int)>;   // value of component n at coarse cell (i,j,k)

namespace {

constexpr int NC = 3;   // tracked species

double hash01(int i, int j, int k, int c)
{
    unsigned long long h = 1469598103934665603ULL;
    for (long long v : {static_cast<long long>(i), static_cast<long long>(j), static_cast<long long>(k), static_cast<long long>(c)}) { h ^= static_cast<unsigned long long>(v + 7919); h *= 1099511628211ULL; h ^= h >> 31; }
    return static_cast<double>(h % 1000003ULL) / 1000003.0;
}
double smooth(int i, int j, int k, int n) { return 0.35 + 0.2 * std::sin(0.45 * i + 0.3 * n) * std::cos(0.35 * j) + 0.1 * std::sin(0.25 * k + n); }
double rough(int i, int j, int k, int n) { double v = 0.02 + hash01(i, j, k, n); if (i >= 7 && i <= 9) v *= 40.0; return v; }   // positive, steep, with 40x fronts

std::string sci(double v) { char b[32]; std::snprintf(b, sizeof b, "%.2e", v); return b; }

struct Setup {
    amrex::Geometry cg, fg;
    amrex::BoxArray cba;
    amrex::DistributionMapping cdm;
    amrex::IntVect ratio;
};

Setup make_setup(const amrex::IntVect& ratio, bool hidden)
{
    Setup s;
    s.ratio = ratio;
    const amrex::IntVect n = hidden ? amrex::IntVect(16, 1, 16) : amrex::IntVect(16, 16, 16);
    amrex::Box dom(amrex::IntVect(0), n - 1);
    amrex::RealBox rb({0., 0., 0.}, {1., 1., 1.});
    s.cg = amrex::Geometry(dom, rb, 0, {1, 1, 1});
    s.fg = amrex::refine(s.cg, ratio);
    s.cba = amrex::BoxArray(dom);
    s.cba.maxSize(hidden ? amrex::IntVect(8, 1, 8) : amrex::IntVect(8));
    s.cdm = amrex::DistributionMapping(s.cba);
    return s;
}

amrex::MultiFab make_coarse(const Setup& s, const Fn& f)
{
    amrex::MultiFab mf(s.cba, s.cdm, NC, 1);
    mf.setVal(0.0);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf.array(mfi);
        for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, n) = f(i, j, k, n); });
    }
    mf.FillBoundary(s.cg.periodicity());
    return mf;
}

amrex::BoxArray fine_ba_from_coarse_box(const amrex::Box& cb, const amrex::IntVect& ratio, int maxsize)
{
    amrex::BoxArray ba(amrex::refine(cb, ratio));
    ba.maxSize(maxsize);
    return ba;
}

// composite integral of every component: uncovered coarse cells (volume 1) plus fine cells (volume 1/r^3)
std::array<double, NC> composite(const Setup& s, const amrex::MultiFab& C, const amrex::MultiFab* F)
{
    std::array<long double, NC> sum{};
    amrex::iMultiFab cov = make_covered_mask(s.cba, s.cdm, F ? F->boxArray() : amrex::BoxArray(), s.ratio);
    for (amrex::MFIter mfi(C); mfi.isValid(); ++mfi) {
        auto a = C.const_array(mfi); auto m = cov.const_array(mfi);
        for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (!m(i, j, k)) sum[n] += a(i, j, k, n); });
    }
    if (F) {
        const long double vf = 1.0L / (s.ratio[0] * s.ratio[1] * s.ratio[2]);
        for (amrex::MFIter mfi(*F); mfi.isValid(); ++mfi) {
            auto a = F->const_array(mfi);
            for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { sum[n] += vf * a(i, j, k, n); });
        }
    }
    std::array<double, NC> out{};
    for (int n = 0; n < NC; ++n) { double v = static_cast<double>(sum[n]); amrex::ParallelAllReduce::Sum(v, amrex::ParallelContext::CommunicatorSub()); out[n] = v; }
    return out;
}

double min_value(const amrex::MultiFab& mf) { return mf.min(0, 0, false) < mf.min(1, 0, false) ? (mf.min(0, 0, false) < mf.min(2, 0, false) ? mf.min(0, 0, false) : mf.min(2, 0, false)) : (mf.min(1, 0, false) < mf.min(2, 0, false) ? mf.min(1, 0, false) : mf.min(2, 0, false)); }

// Y_n = rho*Z_n / rho with rho = sum_n: checks the range [0,1] and the sum; returns the number of violations
long check_mass_fractions(const amrex::MultiFab& F)
{
    long bad = 0;
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            double rho = 0; for (int n = 0; n < NC; ++n) rho += a(i, j, k, n);
            double ys = 0; bool ok = rho > 0;
            for (int n = 0; n < NC && ok; ++n) { const double y = a(i, j, k, n) / rho; ys += y; if (y < 0.0 || y > 1.0) ok = false; }
            if (!ok || std::abs(ys - 1.0) > 1e-14) ++bad;
        });
    }
    amrex::ParallelAllReduce::Sum(bad, amrex::ParallelContext::CommunicatorSub());
    return bad;
}

// Regrid scenario: old fine patch -> new fine patch (partial overlap, one part dropped, one part new). Returns max relative composite change.
double regrid_budget(const Setup& s, const Fn& f, const std::string& tag, bool perturb, ProlongStats* stats_out = nullptr, long* overlap_bad = nullptr, long* overlap_n = nullptr)
{
    const amrex::IntVect r = s.ratio;
    const bool hid = s.cg.Domain().length(1) == 1;
    const amrex::Box old_cb(amrex::IntVect(3, 0, 3), hid ? amrex::IntVect(8, 0, 8) : amrex::IntVect(8, 8, 8));
    const amrex::Box new_cb(amrex::IntVect(6, 0, 6), hid ? amrex::IntVect(13, 0, 13) : amrex::IntVect(13, 7, 13));
    amrex::MultiFab C = make_coarse(s, f);
    amrex::BoxArray oba = fine_ba_from_coarse_box(old_cb, r, 8);
    amrex::MultiFab Fold(oba, amrex::DistributionMapping(oba), NC, 0);
    ProlongOpts opt;
    fdsrt::prolong_conserved(Fold, C, s.cg, r, opt);
    for (amrex::MFIter mfi(Fold); mfi.isValid(); ++mfi) {   // emulate evolved fine data: multiplicative noise, positive
        auto a = Fold.array(mfi);
        for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, n) *= 0.8 + 0.4 * hash01(i, j, k, 50 + n); });
    }
    fdsrt::average_down_cells(Fold, C, s.fg, s.cg, r, 0, NC);   // covered coarse cells now hold the fine average
    C.FillBoundary(s.cg.periodicity());
    const auto before = composite(s, C, &Fold);

    amrex::BoxArray nba = fine_ba_from_coarse_box(new_cb, r, 8);
    amrex::MultiFab Fnew(nba, amrex::DistributionMapping(nba), NC, 0);
    const ProlongStats st = fdsrt::prolong_conserved(Fnew, C, s.cg, r, opt, &Fold);
    if (stats_out) *stats_out = st;
    if (overlap_bad) {   // cells present in both: bitwise equal to the old data
        amrex::MultiFab X(nba, Fnew.DistributionMap(), NC, 0);
        X.setVal(-1.0);
        X.ParallelCopy(Fold, 0, 0, NC, 0, 0);
        long bad = 0, n = 0;
        for (amrex::MFIter mfi(X); mfi.isValid(); ++mfi) {
            auto x = X.const_array(mfi); auto y = Fnew.const_array(mfi);
            for (int c = 0; c < NC; ++c) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (x(i, j, k, c) != -1.0) { ++n; if (x(i, j, k, c) != y(i, j, k, c)) ++bad; } });
        }
        amrex::ParallelAllReduce::Sum(bad, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Sum(n, amrex::ParallelContext::CommunicatorSub());
        *overlap_bad = bad; *overlap_n = n;
    }
    if (perturb) {   // negative control: one fine cell (inside the new patch) is changed by 1e-3 relative
        const amrex::Box b0 = nba[0]; const amrex::IntVect pc = b0.smallEnd() + amrex::IntVect(1, 0, 1);
        for (amrex::MFIter mfi(Fnew); mfi.isValid(); ++mfi)
            if (mfi.validbox().contains(pc)) Fnew[mfi].plus<amrex::RunOn::Host>(1e-3 * Fnew[mfi](pc, 0), amrex::Box(pc, pc), 0, 1);
    }
    fdsrt::average_down_cells(Fnew, C, s.fg, s.cg, r, 0, NC);
    const auto after = composite(s, C, &Fnew);
    double worst = 0;
    for (int n = 0; n < NC; ++n) worst = std::max(worst, std::abs(after[n] - before[n]) / std::abs(before[n]));
    const long yviol = perturb ? 0 : check_mass_fractions(Fnew);   // collective: outside the IOProcessor branch
    const double minv = perturb ? 0.0 : min_value(Fnew);
    if (!perturb && amrex::ParallelDescriptor::IOProcessor())
        std::printf("  %-44s relative composite change %.2e, parents %ld, limited %ld, clips %ld, Y violations %ld, min value %.3g\n", tag.c_str(), worst, st.parents, st.limited, st.clips, yviol, minv);
    return worst;
}

}  // namespace

int main(int, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    {
        // (1) regrid budget, ratios 2 and 4, smooth and rough positive data, and the hidden-direction (2-D) case
        struct Case { std::string name; amrex::IntVect ratio; bool hidden; Fn f; };
        const std::vector<Case> cases = {{"ratio 2 smooth", amrex::IntVect(2), false, smooth}, {"ratio 2 rough", amrex::IntVect(2), false, rough},
                                         {"ratio 4 smooth", amrex::IntVect(4), false, smooth}, {"ratio 4 rough", amrex::IntVect(4), false, rough},
                                         {"ratio (2,1,2) 2-D smooth", amrex::IntVect(2, 1, 2), true, smooth}, {"ratio (4,1,4) 2-D rough", amrex::IntVect(4, 1, 4), true, rough}};
        for (const Case& c : cases) {
            Setup s = make_setup(c.ratio, c.hidden);
            ProlongStats st; long ob = -1, on = 0;
            const double w = regrid_budget(s, c.f, c.name, false, &st, &ob, &on);
            CHECK_MSG(w <= 1e-12, c.name + ": composite mass and species change per regrid <= 1e-12 relative, got " + std::to_string(w));
            CHECK_MSG(st.clips == 0, c.name + ": zero clips on positive data, got " + std::to_string(st.clips));
            CHECK_MSG(ob == 0 && on > 0, c.name + ": overlap with the old fine data is copied bitwise, cells " + std::to_string(on) + " bad " + std::to_string(ob));
            // negative control: the same scenario with one fine cell perturbed by 1e-6 must violate the budget
            const double wp = regrid_budget(s, c.f, c.name, true);
            CHECK_MSG(wp > 1e-11, c.name + ": negative control (one fine cell changed by 1e-3) is detected, relative change " + sci(wp));
        }
        // (2) the operator on its own: per-parent conservation, bounds, mass fractions, DMP; unlimited central slopes would go negative on the rough data
        for (int rr : {2, 4}) {
            Setup s = make_setup(amrex::IntVect(rr), false);
            amrex::MultiFab C = make_coarse(s, rough);
            amrex::BoxArray fba = fine_ba_from_coarse_box(amrex::Box(amrex::IntVect(4, 4, 4), amrex::IntVect(11, 11, 11)), s.ratio, 16);
            amrex::MultiFab F(fba, amrex::DistributionMapping(fba), NC, 0);
            ProlongOpts opt;
            const ProlongStats st = fdsrt::prolong_conserved(F, C, s.cg, s.ratio, opt);
            long bad_sum = 0, bad_dmp = 0, bad_unlim = 0, n_dmp = 0;
            for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
                auto a = F.const_array(mfi);
                amrex::LoopOnCpu(amrex::coarsen(mfi.validbox(), s.ratio), [&](int i, int j, int k) {
                    for (int n = 0; n < NC; ++n) {
                        double sum = 0, lo = 1e300, hi = -1e300;
                        for (int di = -1; di <= 1; ++di) for (int dj = -1; dj <= 1; ++dj) for (int dk = -1; dk <= 1; ++dk) { const double v = rough(i + di, j + dj, k + dk, n); lo = std::min(lo, v); hi = std::max(hi, v); }
                        for (int kc = 0; kc < rr; ++kc) for (int jc = 0; jc < rr; ++jc) for (int ic = 0; ic < rr; ++ic) {
                            const double v = a(i * rr + ic, j * rr + jc, k * rr + kc, n);
                            sum += v; ++n_dmp; if (v < lo * (1 - 1e-13) || v > hi * (1 + 1e-13)) ++bad_dmp;
                        }
                        if (std::abs(sum - rr * rr * rr * rough(i, j, k, n)) > 1e-13 * rr * rr * rr * rough(i, j, k, n)) ++bad_sum;
                        // what unlimited central slopes would give at the corner child: negative values appear on rough data
                        double dev = 0; for (int d = 0; d < 3; ++d) { const int e[3] = {d == 0, d == 1, d == 2}; dev += std::abs(0.5 * (rough(i + e[0], j + e[1], k + e[2], n) - rough(i - e[0], j - e[1], k - e[2], n))) * (0.5 - 0.5 / rr); }
                        if (rough(i, j, k, n) - dev < 0) ++bad_unlim;
                    }
                });
            }
            long v3[4] = {bad_sum, bad_dmp, bad_unlim, n_dmp}; amrex::ParallelAllReduce::Sum(v3, 4, amrex::ParallelContext::CommunicatorSub());
            const std::string t = "operator ratio " + std::to_string(rr) + " rough data";
            CHECK_MSG(v3[0] == 0, t + ": children average to the parent in every component, bad parents " + std::to_string(v3[0]));
            CHECK_MSG(v3[1] == 0 && v3[3] > 0, t + ": children stay inside the min/max of the 3x3x3 neighbourhood (no new extremes), bad " + std::to_string(v3[1]));
            CHECK_MSG(v3[2] > 0, t + ": negative control: unlimited central slopes would give negative values at " + std::to_string(v3[2]) + " parents");
            CHECK_MSG(st.clips == 0 && check_mass_fractions(F) == 0, t + ": zero clips, mass fractions in [0,1] and summing to 1");
            CHECK_MSG(st.limited > 0, t + ": the limiter was active");
        }
        // (3) exact for a linear field, bitwise for a constant, independent of the box split
        {
            Setup s = make_setup(amrex::IntVect(4, 2, 2), false);
            const Fn lin = [](int i, int j, int k, int n) { return 5.0 + 0.3 * i + 0.2 * j + 0.1 * k + n; };
            amrex::MultiFab C = make_coarse(s, lin);
            amrex::BoxArray fba = fine_ba_from_coarse_box(amrex::Box(amrex::IntVect(5, 5, 5), amrex::IntVect(10, 10, 10)), s.ratio, 16);
            amrex::MultiFab F(fba, amrex::DistributionMapping(fba), NC, 0);
            ProlongOpts opt;
            fdsrt::prolong_conserved(F, C, s.cg, s.ratio, opt);
            double err = 0;
            for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
                auto a = F.const_array(mfi);
                for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    const double x = (i + 0.5) / s.ratio[0] - 0.5, y = (j + 0.5) / s.ratio[1] - 0.5, z = (k + 0.5) / s.ratio[2] - 0.5;   // child centre in coarse-cell units
                    err = std::max(err, std::abs(a(i, j, k, n) - lin(0, 0, 0, n) - 0.3 * x - 0.2 * y - 0.1 * z));
                });
            }
            amrex::ParallelAllReduce::Max(err, amrex::ParallelContext::CommunicatorSub());
            CHECK_MSG(err < 1e-13, "linear field is reproduced at the child centres, max error " + std::to_string(err));

            amrex::MultiFab Cc = make_coarse(s, [](int, int, int, int n) { return 0.37 + n; });
            fdsrt::prolong_conserved(F, Cc, s.cg, s.ratio, opt);
            double dev = 0; for (int n = 0; n < NC; ++n) { amrex::MultiFab t(F.boxArray(), F.DistributionMap(), 1, 0); amrex::MultiFab::Copy(t, F, n, 0, 1, 0); t.plus(-(0.37 + n), 0, 1); dev = std::max(dev, t.norminf(0)); }
            CHECK_MSG(dev == 0.0, "constant field stays bitwise constant, max deviation " + std::to_string(dev));

            amrex::MultiFab C2 = make_coarse(s, rough);
            amrex::BoxArray fa = fine_ba_from_coarse_box(amrex::Box(amrex::IntVect(4, 4, 4), amrex::IntVect(11, 11, 11)), s.ratio, 32);
            amrex::BoxArray fb = fine_ba_from_coarse_box(amrex::Box(amrex::IntVect(4, 4, 4), amrex::IntVect(11, 11, 11)), s.ratio, 8);
            amrex::MultiFab Fa(fa, amrex::DistributionMapping(fa), NC, 0), Fb(fb, amrex::DistributionMapping(fb), NC, 0), Fc(fa, Fa.DistributionMap(), NC, 0);
            fdsrt::prolong_conserved(Fa, C2, s.cg, s.ratio, opt);
            fdsrt::prolong_conserved(Fb, C2, s.cg, s.ratio, opt);
            Fc.ParallelCopy(Fb, 0, 0, NC, 0, 0);
            amrex::MultiFab::Subtract(Fc, Fa, 0, 0, NC, 0);
            CHECK_MSG(Fc.norminf(0) == 0.0 && Fc.norminf(1) == 0.0 && Fc.norminf(2) == 0.0, "result does not depend on how the fine level is cut into boxes (bitwise)");
        }
        // (4) clip: children below the floor are raised, the others rescaled, the parent sum is kept, every clip is counted; negative parent is counted as unfixable
        {
            double ch[4] = {-0.1, 0.5, 0.7, 0.9};   // sum 2.0 = 4 * 0.5
            const int r = clip_children(ch, 4, 0.5, 0.0);
            double sum = 0, mn = 1e9; for (double v : ch) { sum += v; mn = std::min(mn, v); }
            CHECK_MSG(r == 1 && std::abs(sum - 2.0) < 1e-15 && mn >= 0.0, "clip_children: returns 1, children sum kept, none below the floor");
            double ok[4] = {0.1, 0.2, 0.3, 0.4};
            CHECK_MSG(clip_children(ok, 4, 0.25, 0.0) == 0 && ok[0] == 0.1 && ok[3] == 0.4, "clip_children: untouched when no child is below the floor");
            double neg[4] = {-0.1, -0.1, -0.1, -0.1};
            CHECK_MSG(clip_children(neg, 4, -0.1, 0.0) == 2 && neg[2] == -0.1, "clip_children: negative parent is reported (2) and kept piecewise constant");
            Setup s = make_setup(amrex::IntVect(2), false);
            amrex::MultiFab C = make_coarse(s, [](int i, int j, int k, int n) { return (i == 7 && j == 7 && k == 7 && n == 1) ? -0.2 : 0.5; });
            amrex::BoxArray fba = fine_ba_from_coarse_box(amrex::Box(amrex::IntVect(4, 4, 4), amrex::IntVect(11, 11, 11)), s.ratio, 16);
            amrex::MultiFab F(fba, amrex::DistributionMapping(fba), NC, 0);
            ProlongOpts opt;
            const ProlongStats st = fdsrt::prolong_conserved(F, C, s.cg, s.ratio, opt);
            CHECK_MSG(st.clips >= 1 && st.unfixable == 1, "negative control: a negative parent is counted as a clip (" + std::to_string(st.clips) + ") and as unfixable (" + std::to_string(st.unfixable) + ")");
            // and with the floor switched off nothing is counted
            opt.use_floor = false;
            const ProlongStats st2 = fdsrt::prolong_conserved(F, C, s.cg, s.ratio, opt);
            CHECK_MSG(st2.clips == 0, "use_floor = false: no clip counted");
        }
    }
    const long nfail = fdstest::report("test_celltransfer");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
