// test_tagging.cpp (R4 step 2): tagging kernels against an independent cell-by-cell evaluation, FR-011 style.
// The reference is a plain function of the cell index (periodic wrap), no arrays; the kernels run on a MultiFab split into boxes over the ranks,
// so box edges and rank boundaries are exercised through the ghost exchange. Usage: test_tagging (run on 1 and 4 ranks).
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_TagBox.H>

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

#include "TagOps.H"
#include "check.H"

using namespace fdsrt;

namespace {

constexpr int N = 24;   // periodic cube; boxes of 8^3

int wrap(int i) { return ((i % N) + N) % N; }

// smooth part plus a sharp front (a jump of 40 K) so that every criterion has both tagged and untagged cells
double temp(int i, int j, int k)
{
    i = wrap(i); j = wrap(j); k = wrap(k);
    double t = 300.0 + 2.0 * std::sin(0.3 * i + 0.2 * j) + 1.5 * std::cos(0.25 * k);
    if (i >= 10 && i < 14 && j >= 6 && k < 18) t += 40.0;
    if (i >= 3 && i < 5 && j >= 20 && k >= 4 && k < 6) t += 25.0;
    if (i >= 8 && i < 10 && j >= 8 && j < 12 && k < 6) t += 10.0;      // 10 K steps: tagged only with the relaxed (kept) threshold of the hysteresis test
    if (i >= 18 && i < 20 && j >= 8 && j < 12 && k >= 18) t += 10.0;   // same step outside the covered region
    return t;
}
double rho(int i, int j, int k) { return 1.2 * 300.0 / temp(i, j, k); }
double zz(int i, int j, int k) { i = wrap(i); j = wrap(j); k = wrap(k); return 0.1 + 0.05 * std::sin(0.5 * i) * std::cos(0.4 * j) + (k >= 12 && k < 14 ? 0.3 : 0.0); }   // rho*Z stored
double yval(int i, int j, int k) { return zz(i, j, k) / rho(i, j, k); }

amrex::MultiFab make_field(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::Geometry& geom, const std::function<double(int, int, int)>& f, int ng)
{
    amrex::MultiFab mf(ba, dm, 1, ng);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = f(i, j, k); });
    }
    mf.FillBoundary(geom.periodicity());   // ghost cells come from the neighbour boxes / other ranks / periodic wrap
    return mf;
}

// number of cells where the tag differs from ref(i,j,k) (summed over ranks)
long mismatches(const amrex::TagBoxArray& tags, const std::function<bool(int, int, int)>& ref)
{
    long bad = 0;
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
        auto a = tags.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if ((a(i, j, k) != 0) != ref(i, j, k)) ++bad; });
    }
    amrex::ParallelAllReduce::Sum(bad, amrex::ParallelContext::CommunicatorSub());
    return bad;
}

double maxdiff(const std::function<double(int, int, int)>& v, int i, int j, int k, bool rel, const bool dirs[3])
{
    double best = 0.0;
    const int off[3][2][3] = {{{-1, 0, 0}, {1, 0, 0}}, {{0, -1, 0}, {0, 1, 0}}, {{0, 0, -1}, {0, 0, 1}}};
    for (int d = 0; d < 3; ++d) {
        if (!dirs[d]) continue;
        for (int s = 0; s < 2; ++s) {
            const double a = v(i, j, k), b = v(i + off[d][s][0], j + off[d][s][1], k + off[d][s][2]);
            double df = std::abs(a - b);
            if (rel) { const double m = std::max(std::abs(a), std::abs(b)); df = m > 0 ? df / m : 0.0; }
            best = std::max(best, df);
        }
    }
    return best;
}

}  // namespace

int main(int argc, char** argv)
{
    int one = 1;
    (void)argc;
    amrex::Initialize(one, argv);
    {
        amrex::Box dom(amrex::IntVect(0), amrex::IntVect(N - 1));
        amrex::RealBox rb({0., 0., 0.}, {1., 1., 1.});
        amrex::Geometry geom(dom, rb, 0, {1, 1, 1});
        amrex::BoxArray ba(dom);
        ba.maxSize(8);
        amrex::DistributionMapping dm(ba);
        const bool all[3] = {true, true, true};

        auto T = make_field(ba, dm, geom, temp, 1);
        auto R = make_field(ba, dm, geom, rho, 1);
        auto Zr = make_field(ba, dm, geom, zz, 1);
        auto ymf = make_field(ba, dm, geom, yval, 1);

        auto fresh = [&]() { return amrex::TagBoxArray(ba, dm, 0); };
        auto temp_fn = std::function<double(int, int, int)>(temp);

        {   // (1) temperature rise above ambient
            auto tags = fresh();
            TagCriterion c; c.mode = TagMode::Above; c.base = 300.0; c.thr = 20.0;
            tag_cells(tags, T, nullptr, 0, nullptr, c);
            const long bad = mismatches(tags, [&](int i, int j, int k) { return temp(i, j, k) - 300.0 > 20.0; });
            const long n = count_tags(tags);
            CHECK_MSG(bad == 0, "temperature above ambient: tags equal the independent evaluation, mismatches " + std::to_string(bad));
            CHECK_MSG(n > 0 && n < N * N * N, "temperature above ambient: some but not all cells tagged, " + std::to_string(n));
            // negative control: a reference with a slightly different threshold must disagree (the comparison can fail)
            const long ctl = mismatches(tags, [&](int i, int j, int k) { return temp(i, j, k) - 300.0 > 27.0; });
            CHECK_MSG(ctl > 0, "negative control: a different threshold gives mismatches, " + std::to_string(ctl));
        }
        {   // (2) temperature difference, absolute and relative; field is the same MultiFab
            for (int rel = 0; rel < 2; ++rel) {
                auto tags = fresh();
                TagCriterion c; c.mode = TagMode::Diff; c.thr = rel ? 0.05 : 15.0; c.relative = rel != 0;
                tag_cells(tags, T, nullptr, 0, nullptr, c);
                const long bad = mismatches(tags, [&](int i, int j, int k) { return maxdiff(temp_fn, i, j, k, rel != 0, all) > c.thr; });
                CHECK_MSG(bad == 0, std::string("temperature ") + (rel ? "relative" : "absolute") + " undivided difference: mismatches " + std::to_string(bad));
                CHECK_MSG(count_tags(tags) > 0, "difference criterion tags something");
            }
        }
        {   // (3) density (relative difference)
            auto tags = fresh();
            TagCriterion c; c.mode = TagMode::Diff; c.thr = 0.08; c.relative = true;
            tag_cells(tags, R, nullptr, 0, nullptr, c);
            const long bad = mismatches(tags, [&](int i, int j, int k) { return maxdiff(rho, i, j, k, true, all) > 0.08; });
            CHECK_MSG(bad == 0, "density relative difference: mismatches " + std::to_string(bad));
            CHECK_MSG(count_tags(tags) > 0, "density criterion tags something");
        }
        {   // (4) species mass fraction Y = (rho*Z)/rho from two fields
            auto tags = fresh();
            TagCriterion c; c.mode = TagMode::Diff; c.thr = 0.15;
            tag_cells(tags, Zr, &R, 0, nullptr, c);
            const long bad = mismatches(tags, [&](int i, int j, int k) { return maxdiff(yval, i, j, k, false, all) > 0.15; });
            CHECK_MSG(bad == 0, "species mass fraction difference (rho*Z / rho): mismatches " + std::to_string(bad));
            CHECK_MSG(count_tags(tags) > 0, "species criterion tags something");
        }
        {   // (5) heat release rate per unit volume threshold (a cell field above a level; ambient 0)
            auto q = make_field(ba, dm, geom, [](int i, int j, int k) { return (wrap(i) < 6 && wrap(j) > 12 && wrap(k) < 9) ? 150.0 + wrap(i) : 3.0; }, 0);
            auto tags = fresh();
            TagCriterion c; c.mode = TagMode::Above; c.base = 0.0; c.thr = 100.0;
            tag_cells(tags, q, nullptr, 0, nullptr, c);
            const long bad = mismatches(tags, [&](int i, int j, int k) { return wrap(i) < 6 && wrap(j) > 12 && wrap(k) < 9; });
            CHECK_MSG(bad == 0, "HRRPUV threshold: mismatches " + std::to_string(bad));
        }
        {   // (6) hysteresis: threshold 15 K difference, keep factor 0.5 inside the region covered by the finer level
            const double thr = 15.0, keep = 0.5;
            auto in_cov = [](int i, int j, int k) { return i >= 6 && i < 18 && j >= 4 && j < 20 && k >= 0 && k < 12; };
            amrex::BoxArray fine_ba(amrex::Box(amrex::IntVect(12, 8, 0), amrex::IntVect(35, 39, 23)));   // refined (ratio 2) footprint = coarse cells 6..17, 4..19, 0..11
            amrex::iMultiFab cov = make_covered_mask(ba, dm, fine_ba, amrex::IntVect(2));
            long cn = 0;
            for (amrex::MFIter mfi(cov); mfi.isValid(); ++mfi) { auto a = cov.const_array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (a(i, j, k)) ++cn; }); }
            amrex::ParallelAllReduce::Sum(cn, amrex::ParallelContext::CommunicatorSub());
            CHECK_MSG(cn == 12 * 16 * 12, "covered mask has 12x16x12 cells, got " + std::to_string(cn));
            TagCriterion c; c.mode = TagMode::Diff; c.thr = thr; c.keepfac = keep;
            auto tags = fresh();
            tag_cells(tags, T, nullptr, 0, &cov, c);
            const long bad = mismatches(tags, [&](int i, int j, int k) { return maxdiff(temp_fn, i, j, k, false, all) > (in_cov(i, j, k) ? thr * keep : thr); });
            CHECK_MSG(bad == 0, "hysteresis: covered cells use thr*keep, others thr, mismatches " + std::to_string(bad));
            // negative control: without the mask the result differs from the hysteresis reference
            auto tags2 = fresh();
            tag_cells(tags2, T, nullptr, 0, nullptr, c);
            const long ctl = mismatches(tags2, [&](int i, int j, int k) { return maxdiff(temp_fn, i, j, k, false, all) > (in_cov(i, j, k) ? thr * keep : thr); });
            CHECK_MSG(ctl > 0, "negative control: no mask gives mismatches against the hysteresis reference, " + std::to_string(ctl));
            CHECK_MSG(count_tags(tags) > count_tags(tags2), "hysteresis keeps more tags than the plain threshold");
        }
        {   // (7) single-cell direction: dirs[1] = 0, the y neighbours are not read (poisoned ghost values in y must not matter)
            auto q = make_field(ba, dm, geom, temp, 1);
            for (amrex::MFIter mfi(q); mfi.isValid(); ++mfi) {
                auto a = q.array(mfi);
                const amrex::Box vb = mfi.validbox();
                const amrex::Box gb = mfi.growntilebox();
                amrex::LoopOnCpu(gb, [&](int i, int j, int k) { if (j < vb.smallEnd(1) || j > vb.bigEnd(1)) a(i, j, k) = 1.0e30; });
            }
            auto tags = fresh();
            TagCriterion c; c.mode = TagMode::Diff; c.thr = 15.0; c.dirs = {1, 0, 1};
            // the y ghosts of interior boxes are poisoned too, so compare against a reference that reads neither y neighbours nor poisoned cells:
            // cells whose x and z neighbours are valid or good x/z ghosts. Poison sits only in y ghost layers, so x/z neighbours of valid cells are untouched.
            tag_cells(tags, q, nullptr, 0, nullptr, c);
            const bool xz[3] = {true, false, true};
            const long bad = mismatches(tags, [&](int i, int j, int k) { return maxdiff(temp_fn, i, j, k, false, xz) > 15.0; });
            CHECK_MSG(bad == 0, "single-cell direction: y neighbours are not read, mismatches " + std::to_string(bad));
        }
        {   // (8) force boxes, clip to the refinable region, count
            auto tags = fresh();
            amrex::Box ub1(amrex::IntVect(2, 2, 2), amrex::IntVect(9, 5, 3));          // crosses the box boundary at i = 8
            amrex::Box ub2(amrex::IntVect(16, 16, 16), amrex::IntVect(17, 17, 17));
            tag_boxes(tags, {ub1, ub2});
            CHECK_MSG(count_tags(tags) == ub1.numPts() + ub2.numPts(), "forced boxes: " + std::to_string(count_tags(tags)) + " cells tagged");
            const long bad = mismatches(tags, [&](int i, int j, int k) { return ub1.contains(amrex::IntVect(i, j, k)) || ub2.contains(amrex::IntVect(i, j, k)); });
            CHECK_MSG(bad == 0, "forced boxes: exactly the cells of the boxes");
            amrex::Box allow(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 7, 7));      // keeps part of ub1 only
            const long removed = clip_tags_to_boxes(tags, {allow});
            amrex::Box kept = ub1 & allow;
            CHECK_MSG(count_tags(tags) == kept.numPts(), "clip: remaining tags = tags inside the allowed boxes, " + std::to_string(count_tags(tags)));
            CHECK_MSG(removed == ub1.numPts() - kept.numPts() + ub2.numPts(), "clip: removed count " + std::to_string(removed));
            const long bad2 = mismatches(tags, [&](int i, int j, int k) { return kept.contains(amrex::IntVect(i, j, k)); });
            CHECK_MSG(bad2 == 0, "clip: remaining tags are exactly ub1 inside the allowed box");
            // several allowed boxes, one empty intersection
            auto tags3 = fresh();
            tag_boxes(tags3, {ub1});
            const long removed3 = clip_tags_to_boxes(tags3, {amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(4, 7, 7)), amrex::Box(amrex::IntVect(6, 0, 0), amrex::IntVect(9, 3, 7))});
            amrex::Box k1(amrex::IntVect(2, 2, 2), amrex::IntVect(4, 5, 3)), k2(amrex::IntVect(6, 2, 2), amrex::IntVect(9, 3, 3));
            CHECK_MSG(count_tags(tags3) == k1.numPts() + k2.numPts() && removed3 == ub1.numPts() - k1.numPts() - k2.numPts(), "clip to two boxes: counts");
        }
    }
    const long nfail = fdstest::report("test_tagging");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
