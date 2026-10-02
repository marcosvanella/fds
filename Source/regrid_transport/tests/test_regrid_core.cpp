// test_regrid_core.cpp (R4 step 4): dynamic hierarchy through RegridAmrCore with a store that owns rho*Z_n per level.
//  A  t = 0 hierarchy from tags (initial conditions evaluated directly on each new level), then regrids of a moving blob (prescribed translation,
//     data rebuilt at every cycle as the composite state of the exact blob): conservation across each regrid, the feature stays on the finest level, tags
//     outside the refinable region are discarded and counted, the tag buffer does not reach outside the region, hysteresis.
//  B  the same run with a different distribution of boxes over ranks gives the same grids and bitwise the same data.
//  C  negative controls: no tag buffer loses the feature; skipping the average-down before the regrid breaks the budget.
// The line "HASH hier=... data=..." is printed on rank 0; tests/run_regrid_rank_check.sh runs the test at 1 and 4 ranks and compares the lines.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <memory>
#include <string>

#include "CellTransfer.H"
#include "DriverAdapter.H"
#include "LevelOps.H"
#include "RegridAmrCore.H"
#include "TagOps.H"
#include "check.H"

using namespace fdsrt;

namespace {

constexpr int NC = 3;

// the blob moves with the prescribed velocity: centre x = 0.30 + V*t (V = 0.75 coarse cells = 1.5 level-1 cells per unit time, below the 2-cell tag buffer)
constexpr double V = 0.75 / 32.0;
double blob(double x, double y, double z, double t)
{
    const double xc = 0.30 + V * t, r2 = (x - xc) * (x - xc) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
    return 0.1 + 0.8 * std::exp(-r2 / (2.0 * 0.06 * 0.06));
}
double ic(int comp, double x, double y, double z, double t)
{
    const double b = blob(x, y, z, t);
    if (comp == 0) return b;
    if (comp == 1) return 0.3 * (1.0 + 0.3 * std::sin(6.0 * x + 2.0 * y));
    return 0.2 + 0.1 * b;
}

class Store : public LevelListener, public LevelDataTransfer {
public:
    explicit Store(RegridAmrCore& c) : core(c) {}
    RegridAmrCore& core;
    std::vector<std::unique_ptr<amrex::MultiFab>> q, old;
    double t = 0.0;
    int n_make = 0, n_remake = 0, n_clear = 0, n_begin = 0, n_end = 0, n_done = 0;
    long clips = 0;

    void ensure(int l) { if (static_cast<int>(q.size()) <= l) { q.resize(l + 1); old.resize(l + 1); } }
    void fill_ic(int l)
    {
        const amrex::Geometry& g = core.Geom(l);
        const double dx = g.CellSize(0), dy = g.CellSize(1), dz = g.CellSize(2);
        for (amrex::MFIter mfi(*q[l]); mfi.isValid(); ++mfi) {
            auto a = q[l]->array(mfi);
            for (int n = 0; n < NC; ++n)
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, n) = ic(n, (i + 0.5) * dx, (j + 0.5) * dy, (k + 0.5) * dz, t); });
        }
        q[l]->FillBoundary(g.periodicity());
    }
    void make_level(const LevelLayout& l) override
    {
        ensure(l.level); ++n_make;
        q[l.level] = std::make_unique<amrex::MultiFab>(l.ba, l.dm, NC, 1);
        q[l.level]->setVal(0.0);
        if (l.level == 0) fill_ic(0);
    }
    void remake_level(const LevelLayout& l) override
    {
        ensure(l.level); ++n_remake;
        old[l.level] = std::move(q[l.level]);
        q[l.level] = std::make_unique<amrex::MultiFab>(l.ba, l.dm, NC, 1);
        q[l.level]->setVal(0.0);
    }
    void clear_level(int level) override { ++n_clear; q[level].reset(); old[level].reset(); }
    void begin_regrid() override { ++n_begin; }
    void end_regrid() override { ++n_end; }

    void fill_initial_level(const LevelLayout& l) override { fill_ic(l.level); }
    void prolong(const LevelLayout& l, const amrex::MultiFab* old_fine)
    {
        q[l.level - 1]->FillBoundary(core.Geom(l.level - 1).periodicity());
        const ProlongStats st = prolong_conserved(*q[l.level], *q[l.level - 1], core.Geom(l.level - 1), core.refRatio(l.level - 1), ProlongOpts(), old_fine);
        clips += st.clips;
        q[l.level]->FillBoundary(core.Geom(l.level).periodicity());
    }
    void fill_new_level(const LevelLayout& l) override { prolong(l, nullptr); }
    void fill_remade_level(const LevelLayout& l) override { prolong(l, old[l.level].get()); old[l.level].reset(); }
    void hierarchy_done(bool) override { ++n_done; average_down_all(); }
    void average_down_all()
    {
        for (int l = core.finestLevel(); l >= 1; --l) {
            average_down_cells(*q[l], *q[l - 1], core.Geom(l), core.Geom(l - 1), core.refRatio(l - 1), 0, NC);
            q[l - 1]->FillBoundary(core.Geom(l - 1).periodicity());
        }
        for (int l = 0; l <= core.finestLevel(); ++l) q[l]->FillBoundary(core.Geom(l).periodicity());
    }
    // the exact composite state of the blob at time t on every level, consistent between levels
    void set_state(double tt)
    {
        t = tt;
        for (int l = 0; l <= core.finestLevel(); ++l) fill_ic(l);
        average_down_all();
    }
    std::array<double, NC> composite() const
    {
        std::array<long double, NC> s{};
        for (int l = 0; l <= core.finestLevel(); ++l) {
            amrex::iMultiFab cov = core.covered_mask(l);
            const long double vol = std::pow(0.125L, l);
            for (amrex::MFIter mfi(*q[l]); mfi.isValid(); ++mfi) {
                auto a = q[l]->const_array(mfi); auto m = cov.const_array(mfi);
                for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (!m(i, j, k)) s[n] += vol * a(i, j, k, n); });
            }
        }
        std::array<double, NC> out{};
        for (int n = 0; n < NC; ++n) { double v = static_cast<double>(s[n]); amrex::ParallelAllReduce::Sum(v, amrex::ParallelContext::CommunicatorSub()); out[n] = v; }
        return out;
    }
};

struct Run {
    double worst_budget = 0.0;       // largest relative composite change over the regrids
    long feature_uncovered = 0;      // feature cells at the next regrid time that are not on the finest level, summed over cycles
    int finest = 0;
    long clips = 0;
    uint64_t hier = 14695981039346656037ULL, data = 0;
    RegridStats stats;
    int n_changed = 0;
    long outside_region_cells = 0;   // cells of levels >= 1 outside the refinable region (must be 0)
    int n_make = 0, n_remake = 0, n_clear = 0, n_begin = 0, n_end = 0;
};

uint64_t mix(uint64_t h, uint64_t v) { h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); return h; }

Run run(int nbuf, bool rotated_dm, bool skip_avg, bool hysteresis)
{
    Run r;
    Report rep;
    std::string text = "&AMR MAX_LEVEL=2, REF_RATIO=2, BLOCKING_FACTOR=4, MAX_GRID_SIZE=16, N_ERROR_BUF=" + std::to_string(nbuf) + ", N_PROPER=1 /\n"
                       "&AMR_REGION XB=0.125,0.875,0.125,0.875,0.125,0.875 /\n";
    AmrParams p = parse_amr_params(text, rep);
    MeshInput m;
    m.ijk[0] = m.ijk[1] = m.ijk[2] = 32;
    m.xb[0] = m.xb[2] = m.xb[4] = 0.0; m.xb[1] = m.xb[3] = m.xb[5] = 1.0;
    Hierarchy h;
    const bool ok = build_hierarchy_from_meshes({m}, p, {true, true, true}, h, rep);
    CHECK_MSG(ok && rep.ok(), "hierarchy for the dynamic run");
    if (!ok) return r;
    RegridAmrCore core(h, p);
    Store st(core);
    core.set_data_transfer(&st);
    RegridAmrCore::DmFn dm_fn;
    amrex::DistributionMapping dm0;
    if (rotated_dm) {   // a different box-to-rank assignment on every level (the grids must not depend on it)
        dm_fn = [](int, const amrex::BoxArray& ba) {
            amrex::Vector<int> pm(ba.size());
            for (int i = 0; i < static_cast<int>(pm.size()); ++i) pm[i] = (3 * i + 1) % amrex::ParallelDescriptor::NProcs();
            return amrex::DistributionMapping(pm);
        };
    }
    core.set_tag_function([&](int lev, amrex::TagBoxArray& tags, amrex::Real) {
        TagCriterion c; c.mode = TagMode::Above; c.base = 0.1; c.thr = 0.05; c.keepfac = hysteresis ? 0.8 : 1.0;
        amrex::iMultiFab cov = core.covered_mask(lev);
        tag_cells(tags, *st.q[lev], nullptr, 0, hysteresis ? &cov : nullptr, c);
    });
    st.t = 0.0;
    core.init_from_tags(st, 0.0, true, dm_fn);
    r.finest = core.finestLevel();
    CHECK_MSG(core.finestLevel() == 2, "initial hierarchy reaches MAX_LEVEL = 2 (finest level " + std::to_string(core.finestLevel()) + ")");

    const double dt = 1.0;
    for (int cyc = 0; cyc < 6; ++cyc) {
        const double t = cyc * dt;
        st.set_state(t);
        if (skip_avg) {   // negative control: stale coarse values under the fine patches (as if the coarse level had not seen the fine data)
            for (int l = 0; l < core.finestLevel(); ++l) {
                amrex::iMultiFab cov = core.covered_mask(l);
                for (amrex::MFIter mfi(*st.q[l]); mfi.isValid(); ++mfi) {
                    auto qa = st.q[l]->array(mfi); auto m = cov.const_array(mfi);
                    for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (m(i, j, k)) qa(i, j, k, n) = 0.123; });
                }
            }
        }
        const auto before = st.composite();
        core.regrid_dynamic(t);
        const auto after = st.composite();
        for (int n = 0; n < NC; ++n) r.worst_budget = std::max(r.worst_budget, std::abs(after[n] - before[n]) / std::abs(before[n]));
        // the feature at the next regrid time must lie on the finest level of this hierarchy
        const int lf = core.finestLevel();
        const amrex::Geometry& g = core.Geom(2);
        long unc = 0;
        for (int k = 0; k < g.Domain().length(2); ++k) for (int j = 0; j < g.Domain().length(1); ++j) {
            if (std::abs((j + 0.5) * g.CellSize(1) - 0.5) > 0.2 || std::abs((k + 0.5) * g.CellSize(2) - 0.5) > 0.2) continue;
            for (int i = 0; i < g.Domain().length(0); ++i)
                if (blob((i + 0.5) * g.CellSize(0), (j + 0.5) * g.CellSize(1), (k + 0.5) * g.CellSize(2), t + dt) > 0.15 && !(lf == 2 && core.boxArray(2).contains(amrex::IntVect(i, j, k)))) ++unc;
        }
        r.feature_uncovered += unc;
        for (int l = 0; l <= lf; ++l) for (int b = 0; b < static_cast<int>(core.boxArray(l).size()); ++b) {
            const amrex::Box bx = core.boxArray(l)[b];
            r.hier = mix(r.hier, static_cast<uint64_t>(l * 1000003 + bx.smallEnd(0) * 7919 + bx.smallEnd(1) * 104729 + bx.smallEnd(2) * 1299709 + bx.bigEnd(0) * 15485863 + bx.bigEnd(1) * 32452843 + bx.bigEnd(2) * 49979687));
        }
        // no cell of the finest possible level outside the refinable region (region boxes of the hierarchy, refined); a level below MAX_LEVEL may exceed it
        // by the proper-nesting cover of the next finer level (N_PROPER cells, snapped to the blocking factor of the level), which is allowed here
        for (int l = 1; l <= lf; ++l) {
            amrex::BoxArray allowed;
            amrex::BoxList bl;
            for (const RegionBox& rb : h.levels[l - 1].taggable) {
                amrex::Box rbx = amrex::refine(to_amrex(rb.box), core.refRatio(l - 1));
                if (l < h.max_level) rbx.grow(h.n_proper + h.levels[l].blocking_factor[0]);
                bl.push_back(rbx);
            }
            allowed = amrex::BoxArray(bl);
            long inside = 0;   // the single region box is disjoint from itself, so the intersections do not double count
            for (int b = 0; b < static_cast<int>(core.boxArray(l).size()); ++b)
                for (const auto& is : allowed.intersections(core.boxArray(l)[b])) inside += is.second.numPts();
            r.outside_region_cells += core.boxArray(l).numPts() - inside;
        }
    }
    uint64_t d = 0;   // order-independent checksum of all valid data
    for (int l = 0; l <= core.finestLevel(); ++l)
        for (amrex::MFIter mfi(*st.q[l]); mfi.isValid(); ++mfi) {
            auto a = st.q[l]->const_array(mfi);
            for (int n = 0; n < NC; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                double v = a(i, j, k, n); uint64_t bits; std::memcpy(&bits, &v, 8);
                d += mix(mix(mix(mix(static_cast<uint64_t>(l), static_cast<uint64_t>(i + 1000)), static_cast<uint64_t>(j + 1000)), static_cast<uint64_t>(k + 1000) * 31 + n), bits);
            });
        }
    unsigned long long dd = d;
    amrex::ParallelAllReduce::Sum(dd, amrex::ParallelContext::CommunicatorSub());
    r.data = dd;
    r.clips = st.clips;
    r.stats = core.stats();
    r.n_changed = core.stats().n_changed;
    r.n_make = st.n_make; r.n_remake = st.n_remake; r.n_clear = st.n_clear; r.n_begin = st.n_begin; r.n_end = st.n_end;
    CHECK_MSG(st.n_begin == st.n_end && st.n_begin == 6, "begin_regrid / end_regrid bracket every regrid (6)");
    return r;
}

}  // namespace

int main(int, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    {
        const Run a = run(2, false, false, true);
        if (amrex::ParallelDescriptor::IOProcessor())
            std::printf("  A moving blob: worst composite change per regrid %.2e, feature cells off the finest level %ld, clips %ld, regrids that changed grids %d of 6\n"
                        "    tags set level0/level1: %ld/%ld, discarded outside region %ld/%ld, buffer blocks removed %ld/%ld, make/remake/clear %d/%d/%d\n",
                        a.worst_budget, a.feature_uncovered, a.clips, a.n_changed, a.stats.tags_set[0], a.stats.tags_set[1], a.stats.discarded_outside_region[0],
                        a.stats.discarded_outside_region[1], a.stats.buffer_blocks_discarded[0], a.stats.buffer_blocks_discarded[1], a.n_make, a.n_remake, a.n_clear);
        CHECK_MSG(a.worst_budget <= 1e-12, "A: composite mass and species change per regrid <= 1e-12 relative, got " + std::to_string(a.worst_budget));
        CHECK_MSG(a.clips == 0, "A: zero clips on positive data");
        CHECK_MSG(a.feature_uncovered == 0, "A: with a 2-cell tag buffer every feature cell of the next regrid time is on the finest level, missing " + std::to_string(a.feature_uncovered));
        CHECK_MSG(a.n_changed >= 3, "A: the hierarchy follows the blob (grids changed in " + std::to_string(a.n_changed) + " of 6 regrids)");
        CHECK_MSG(a.n_remake > 0 && a.n_clear >= 0 && a.n_make >= 3, "A: levels were made and remade");
        CHECK_MSG(a.outside_region_cells == 0, "A: no cell of the finest level outside the refinable region, none of a coarser refined level beyond the nesting cover");

        const Run b = run(2, true, false, true);
        CHECK_MSG(a.hier == b.hier && a.data == b.data, "B: a different box-to-rank assignment gives the same grids and bitwise the same data");
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("HASH hier=%016llx data=%016llx\n", static_cast<unsigned long long>(a.hier), static_cast<unsigned long long>(a.data));

        const Run c1 = run(0, false, false, true);
        CHECK_MSG(c1.feature_uncovered > 0, "C: negative control, no tag buffer: feature cells leave the finest level (" + std::to_string(c1.feature_uncovered) + ")");
        const Run c2 = run(2, false, true, true);
        CHECK_MSG(c2.worst_budget > 1e-6, "C: negative control, stale coarse data under the fine patches breaks the budget (change " + std::to_string(c2.worst_budget) + ")");
        const Run c3 = run(2, false, false, false);
        CHECK_MSG(c3.hier != a.hier || c3.stats.tags_set[1] != a.stats.tags_set[1], "C: hysteresis changes the tags (TAG_KEEP is active)");
        CHECK_MSG(a.stats.tags_set[1] >= c3.stats.tags_set[1], "C: with hysteresis at least as many cells stay tagged as without");
    }
    const long nfail = fdstest::report("test_regrid_core");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
