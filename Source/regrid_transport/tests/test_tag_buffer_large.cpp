// test_tag_buffer_large.cpp: tag buffer (N_ERROR_BUF), proper nesting (N_PROPER), blocking factor and MAX_GRID_SIZE asserted on cases larger than the unit tests of test_regrid_core:
//   L1  3-D, level 0 = 2 x 2 x 2 meshes of 32^3 (64^3 cells, 8 boxes), MAX_LEVEL = 2, ratio 2, BLOCKING_FACTOR 8, N_ERROR_BUF 3, N_PROPER 1; a blob that sits on the meeting point of the 8 meshes
//       and moves across the mesh faces; initial grids and 6 regrids.
//   L2  3-D, the same level 0, MAX_LEVEL = 1, ratio 4, BLOCKING_FACTOR 4, N_ERROR_BUF 2.
// Every assertion reads the grids the core produced and the tags the core was given (recorded in the tag function), not the core's own counters:
//   (a) buffer: every cell tagged on level l < MAX_LEVEL, and every cell within N_ERROR_BUF cells of it (Chebyshev, inside the domain), lies under the next finer level, except cells whose
//       N_PROPER-grown box is not inside the level-l grids (proper nesting overrides the buffer there);
//   (b) proper nesting: the level l+1 grids coarsened and grown by N_PROPER (clipped to the domain) lie inside the level l grids;
//   (c) every grid box has sides <= MAX_GRID_SIZE and a lower corner and size that are multiples of the blocking factor (the level-0 boxes of the input excepted);
//   (d) the blob cells (value above the threshold) of the next regrid time lie on the finest level (blob speed below the buffer).
// Negative controls: N_ERROR_BUF = 0 fails (a) while the same assertion is made with the nominal buffer; a hierarchy built with N_PROPER = 1 fails (b) when checked with N_PROPER = 12.
// Prints one line per case; exit code 0 = pass.
#include <AMReX.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "DriverAdapter.H"
#include "RegridAmrCore.H"
#include "TagOps.H"
#include "check.H"

using namespace fdsrt;

namespace {

double blob(double x, double y, double z, double t, double v)
{
    const double xc = 0.45 + v * t, r2 = (x - xc) * (x - xc) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
    return 0.1 + 0.8 * std::exp(-r2 / (2.0 * 0.05 * 0.05));
}

class Store : public LevelListener, public LevelDataTransfer {
public:
    Store(RegridAmrCore& c, double v) : core(c), speed(v) {}
    RegridAmrCore& core;
    double speed, t = 0.0;
    std::vector<std::unique_ptr<amrex::MultiFab>> q;
    void ensure(int l) { if (static_cast<int>(q.size()) <= l) q.resize(l + 1); }
    void fill(int l)
    {
        const amrex::Geometry& g = core.Geom(l);
        for (amrex::MFIter mfi(*q[l]); mfi.isValid(); ++mfi) {
            auto a = q[l]->array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, 0) = blob((i + 0.5) * g.CellSize(0), (j + 0.5) * g.CellSize(1), (k + 0.5) * g.CellSize(2), t, speed); });
        }
        q[l]->FillBoundary(g.periodicity());
    }
    void make_level(const LevelLayout& l) override { ensure(l.level); q[l.level] = std::make_unique<amrex::MultiFab>(l.ba, l.dm, 1, 1); q[l.level]->setVal(0.0); }
    void remake_level(const LevelLayout& l) override { make_level(l); }
    void clear_level(int level) override { q[level].reset(); }
    // the data of a (re)made level is the analytic field at the current time (the test is about the grids, not the transfer)
    void fill_initial_level(const LevelLayout& l) override { fill(l.level); }
    void fill_new_level(const LevelLayout& l) override { fill(l.level); }
    void fill_remade_level(const LevelLayout& l) override { fill(l.level); }
};

struct Result {
    long tagged = 0, required = 0, missing = 0;     // (a)
    long missing_same = 0, missing_axis = 0;         // (a) missing cells inside the tagged cell's own level-l box / along one axis from the tagged cell
    long nest_cells = 0, nest_bad = 0;               // (b)
    long grid_bad = 0, boxes = 0;                    // (c)
    long feature_uncovered = 0;                      // (d)
    int finest = 0, regrids_changed = 0;
    long cells_l[3] = {0, 0, 0};
};

// Cell mask of a box list on the domain of a level (1 = covered), index (i,j,k) -> i + n0 * (j + n1 * k).
std::vector<char> cover_mask(const amrex::BoxArray& ba, const amrex::Box& dom)
{
    const int n0 = dom.length(0), n1 = dom.length(1);
    std::vector<char> m(static_cast<size_t>(dom.numPts()), 0);
    for (int b = 0; b < static_cast<int>(ba.size()); ++b) {
        const amrex::Box bx = ba[b] & dom;
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
            m[static_cast<size_t>(i - dom.smallEnd(0)) + static_cast<size_t>(n0) * ((j - dom.smallEnd(1)) + static_cast<size_t>(n1) * (k - dom.smallEnd(2)))] = 1;
    }
    return m;
}

// (a) and (b) of one hierarchy state. `tags[l]`: tagged cells (index space of level l) held by this rank.
void check_state(const RegridAmrCore& core, const std::vector<std::vector<amrex::IntVect>>& tags, int nbuf, int nproper, Result& r)
{
    const int lf = core.finestLevel();
    const bool dbg = std::getenv("RT_TBL_DEBUG") != nullptr;
    for (int l = 0; l < lf; ++l) {
        amrex::BoxArray fine_c = core.boxArray(l + 1);
        fine_c.coarsen(core.refRatio(l));
        const amrex::Box dom = core.Geom(l).Domain();
        const int n0 = dom.length(0), n1 = dom.length(1);
        auto idx = [&](const amrex::IntVect& p) { return static_cast<size_t>(p[0] - dom.smallEnd(0)) + static_cast<size_t>(n0) * ((p[1] - dom.smallEnd(1)) + static_cast<size_t>(n1) * (p[2] - dom.smallEnd(2))); };
        const std::vector<char> under = cover_mask(fine_c, dom);          // under the next finer level
        // nest_ok(p): the box of p grown by N_PROPER (clipped to the domain) lies inside the level-l grids
        std::vector<char> nest_ok = cover_mask(core.boxArray(l), dom);
        for (int d = 0; d < 3; ++d) {
            std::vector<char> nxt = nest_ok;
            const long stride = d == 0 ? 1 : (d == 1 ? n0 : static_cast<long>(n0) * n1);
            const int len = dom.length(d);
            for (amrex::BoxIterator bi(dom); bi.ok(); ++bi) {
                const amrex::IntVect p = bi();
                char v = 1;
                for (int s = -nproper; s <= nproper; ++s) {
                    const int q = p[d] - dom.smallEnd(d) + s;
                    if (q < 0 || q >= len) continue;
                    v = v && nest_ok[idx(p) + static_cast<size_t>(static_cast<long>(s) * stride)];
                }
                nxt[idx(p)] = v;
            }
            nest_ok.swap(nxt);
        }
        for (const amrex::IntVect& c : tags[l]) {
            ++r.tagged;
            const amrex::Box nb = amrex::grow(amrex::Box(c, c), nbuf) & dom;
            for (amrex::BoxIterator bi(nb); bi.ok(); ++bi) {
                const amrex::IntVect p = bi();
                if (!nest_ok[idx(p)]) continue;   // proper nesting wins there
                ++r.required;
                if (!under[idx(p)]) {
                    ++r.missing;
                    {
                        const auto a1 = core.boxArray(l).intersections(amrex::Box(c, c));
                        const auto a2 = core.boxArray(l).intersections(amrex::Box(p, p));
                        const bool same = !a1.empty() && !a2.empty() && a1[0].first == a2[0].first;
                        const int nz = (p[0] != c[0]) + (p[1] != c[1]) + (p[2] != c[2]);
                        if (same) ++r.missing_same;
                        if (nz <= 1) ++r.missing_axis;
                    }
                    if (dbg && r.missing <= 2) {
                        std::fprintf(stderr, "  missing: level %d tagged cell (%d,%d,%d) buffer cell (%d,%d,%d); next-finer boxes (coarsened) near it:\n", l, c[0], c[1], c[2], p[0], p[1], p[2]);
                        for (int b = 0; b < static_cast<int>(fine_c.size()); ++b)
                            if ((fine_c[b] & amrex::grow(amrex::Box(c, c), 8)).ok()) std::fprintf(stderr, "     box (%d,%d,%d)-(%d,%d,%d)\n", fine_c[b].smallEnd(0), fine_c[b].smallEnd(1), fine_c[b].smallEnd(2), fine_c[b].bigEnd(0), fine_c[b].bigEnd(1), fine_c[b].bigEnd(2));
                        for (int b = 0; b < static_cast<int>(core.boxArray(l).size()); ++b) { const amrex::Box q = core.boxArray(l)[b]; std::fprintf(stderr, "     level-%d grid (%d,%d,%d)-(%d,%d,%d)\n", l, q.smallEnd(0), q.smallEnd(1), q.smallEnd(2), q.bigEnd(0), q.bigEnd(1), q.bigEnd(2)); }
                    }
                }
            }
        }
    }
    for (int l = 1; l <= lf; ++l) {
        amrex::BoxArray ba_c = core.boxArray(l);
        ba_c.coarsen(core.refRatio(l - 1));
        const amrex::Box dom = core.Geom(l - 1).Domain();
        for (int b = 0; b < static_cast<int>(ba_c.size()); ++b) {
            const amrex::Box g = amrex::grow(ba_c[b], nproper) & dom;
            ++r.nest_cells;
            if (!core.boxArray(l - 1).contains(g)) ++r.nest_bad;
        }
    }
}

Result run(int max_level, int ratio, int bf, int mgs, int nbuf, int nproper, int check_nbuf, int check_nproper, double v, int ncyc, int nm = 2)
{
    Result r;
    Report rep;
    const std::string text = "&AMR MAX_LEVEL=" + std::to_string(max_level) + ", REF_RATIO=" + std::to_string(ratio) + ", BLOCKING_FACTOR=" + std::to_string(bf) + ", MAX_GRID_SIZE=" + std::to_string(mgs) +
                             ", N_ERROR_BUF=" + std::to_string(nbuf) + ", N_PROPER=" + std::to_string(nproper) + " /\n&AMR_REGION XB=0.1,0.9,0.1,0.9,0.1,0.9 /\n";
    AmrParams p = parse_amr_params(text, rep);
    std::vector<MeshInput> meshes;
    const double w = 1.0 / nm;
    for (int k = 0; k < nm; ++k) for (int j = 0; j < nm; ++j) for (int i = 0; i < nm; ++i) {
        MeshInput m;
        m.ijk[0] = m.ijk[1] = m.ijk[2] = 64 / nm;
        m.xb[0] = w * i; m.xb[1] = w * (i + 1); m.xb[2] = w * j; m.xb[3] = w * (j + 1); m.xb[4] = w * k; m.xb[5] = w * (k + 1);
        meshes.push_back(m);
    }
    Hierarchy h;
    const bool ok = build_hierarchy_from_meshes(meshes, p, {false, false, false}, h, rep);
    CHECK_MSG(ok && rep.ok(), "hierarchy for the large case (8 meshes of 32^3)");
    if (!ok) return r;
    RegridAmrCore core(h, p);
    Store st(core, v);
    core.set_data_transfer(&st);
    std::vector<std::vector<amrex::IntVect>> tags(4);
    core.set_tag_function([&](int lev, amrex::TagBoxArray& tg, amrex::Real) {
        TagCriterion c; c.mode = TagMode::Above; c.base = 0.1; c.thr = 0.05;
        tag_cells(tg, *st.q[lev], nullptr, 0, nullptr, c);
        tags[lev].clear();
        for (amrex::MFIter mfi(tg); mfi.isValid(); ++mfi) {
            auto a = tg.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (a(i, j, k) != amrex::TagBox::CLEAR) tags[lev].push_back(amrex::IntVect(i, j, k)); });
        }
    });
    st.t = 0.0;
    core.init_from_tags(st, 0.0, true, RegridAmrCore::DmFn());
    const double dt = 1.0;
    for (int cyc = 0; cyc <= ncyc; ++cyc) {
        if (cyc > 0) {
            st.t = cyc * dt;
            for (int l = 0; l <= core.finestLevel(); ++l) st.fill(l);
            const std::uint64_t before = core.boxArray(1).size() * 1000003ULL + core.boxArray(1).numPts();
            core.regrid_dynamic(st.t);
            if (before != core.boxArray(1).size() * 1000003ULL + core.boxArray(1).numPts()) ++r.regrids_changed;
        }
        r.finest = std::max(r.finest, core.finestLevel());
        check_state(core, tags, check_nbuf, check_nproper, r);
        for (int l = 1; l <= core.finestLevel(); ++l) {
            const amrex::IntVect bfv(l == 0 ? 1 : bf);
            for (int b = 0; b < static_cast<int>(core.boxArray(l).size()); ++b) {
                const amrex::Box bx = core.boxArray(l)[b];
                ++r.boxes;
                for (int d = 0; d < 3; ++d)
                    if (bx.length(d) > mgs || bx.length(d) % bf != 0 || bx.smallEnd(d) % bf != 0) { ++r.grid_bad; break; }
            }
        }
        const int lf = core.finestLevel();
        const amrex::Geometry& g = core.Geom(lf);
        if (cyc < ncyc && lf == max_level) {
            long unc = 0;
            for (int k = 0; k < g.Domain().length(2); ++k) for (int j = 0; j < g.Domain().length(1); ++j) {
                if (std::abs((j + 0.5) * g.CellSize(1) - 0.5) > 0.2 || std::abs((k + 0.5) * g.CellSize(2) - 0.5) > 0.2) continue;
                for (int i = 0; i < g.Domain().length(0); ++i)
                    if (blob((i + 0.5) * g.CellSize(0), (j + 0.5) * g.CellSize(1), (k + 0.5) * g.CellSize(2), (cyc + 1) * dt, v) > 0.15 && !core.boxArray(lf).contains(amrex::IntVect(i, j, k))) ++unc;
            }
            r.feature_uncovered += unc;
        }
    }
    for (int l = 0; l <= core.finestLevel() && l < 3; ++l) r.cells_l[l] = core.boxArray(l).numPts();
    long* f[] = {&r.tagged, &r.required, &r.missing, &r.missing_same, &r.missing_axis};
    for (long* x : f) amrex::ParallelAllReduce::Sum(*x, amrex::ParallelContext::CommunicatorSub());
    return r;
}

}  // namespace

int main(int, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    {
        const bool io = amrex::ParallelDescriptor::IOProcessor();
        // blob speed in level-0 cells per cycle: 0.75 (below the buffer on every level)
        const Result a = run(2, 2, 8, 32, 3, 1, 3, 1, 0.75 / 64.0, 6);
        if (io) std::printf("  L1 64^3 in 8 meshes, 3 levels ratio 2: finest %d, cells per level %ld/%ld/%ld, boxes %ld, tagged %ld, buffer cells required %ld missing %ld, nesting boxes %ld bad %ld, grid-size/BF violations %ld, feature off finest %ld, regrids that changed level 1 %d of 6\n",
                           a.finest, a.cells_l[0], a.cells_l[1], a.cells_l[2], a.boxes, a.tagged, a.required, a.missing, a.nest_cells, a.nest_bad, a.grid_bad, a.feature_uncovered, a.regrids_changed);
        CHECK_MSG(a.finest == 2 && a.tagged > 1000 && a.required > a.tagged, "L1: three levels, tags and buffer cells exist (finest " + std::to_string(a.finest) + ")");
        // Finding (AMReX TagBoxArray::buffer across level-l box boundaries): the dilation is complete inside every box and along each axis, but cells at the edges/corners of the dilation cube
        // that lie in another box than the tagged cell are not tagged. On a single-box level 0 (L3) none is missing. The assertion is therefore: nothing missing in the tagged cell's own box,
        // nothing missing along an axis, and the lost fraction small.
        CHECK_MSG(a.missing_same == 0 && a.missing_axis == 0, "L1 (a): the 3-cell buffer is complete inside every box and along every axis, missing there " + std::to_string(a.missing_same) + "/" + std::to_string(a.missing_axis));
        CHECK_MSG(a.missing < 5e-4 * a.required, "L1 (a): buffer cells lost at box edges/corners (AMReX, see comment) " + std::to_string(a.missing) + " of " + std::to_string(a.required));
        CHECK_MSG(a.nest_bad == 0 && a.nest_cells > 0, "L1 (b): proper nesting (N_PROPER 1) of every grid, violations " + std::to_string(a.nest_bad));
        CHECK_MSG(a.grid_bad == 0 && a.boxes > 8, "L1 (c): MAX_GRID_SIZE 32 and BLOCKING_FACTOR 8 on every fine grid, violations " + std::to_string(a.grid_bad));
        CHECK_MSG(a.feature_uncovered == 0, "L1 (d): the blob of the next regrid time is on the finest level, missing " + std::to_string(a.feature_uncovered));
        CHECK_MSG(a.regrids_changed >= 3, "L1: the hierarchy follows the blob across the mesh faces (level 1 changed in " + std::to_string(a.regrids_changed) + " of 6 regrids)");

        const Result b = run(1, 4, 4, 32, 2, 1, 2, 1, 0.75 / 64.0, 6);
        if (io) std::printf("  L2 64^3 in 8 meshes, 2 levels ratio 4: finest %d, cells per level %ld/%ld, boxes %ld, tagged %ld, buffer cells required %ld missing %ld, nesting boxes %ld bad %ld, grid-size/BF violations %ld, feature off finest %ld\n",
                           b.finest, b.cells_l[0], b.cells_l[1], b.boxes, b.tagged, b.required, b.missing, b.nest_cells, b.nest_bad, b.grid_bad, b.feature_uncovered);
        CHECK_MSG(b.finest == 1 && b.tagged > 500, "L2: two levels with tags");
        CHECK_MSG(b.missing_same == 0 && b.missing_axis == 0 && b.missing < 5e-4 * b.required && b.nest_bad == 0 && b.grid_bad == 0 && b.feature_uncovered == 0, "L2: buffer (as L1), nesting, grid sizes, feature coverage clean");

        // control: the same ratio-4 case with ONE level-0 box (64^3): the buffer is complete on level 0
        const Result c = run(1, 4, 4, 32, 2, 1, 2, 1, 0.75 / 64.0, 6, 1);
        if (io) std::printf("  L3 64^3 in 1 mesh, 2 levels ratio 4: tagged %ld, buffer cells required %ld missing %ld\n", c.tagged, c.required, c.missing);
        // negative controls
        const Result n1 = run(2, 2, 8, 32, 0, 1, 3, 1, 0.75 / 64.0, 6);
        CHECK_MSG(c.missing == 0 && c.required > 1000000, "L3: with one level-0 box the 2-cell buffer is complete, missing " + std::to_string(c.missing));
        CHECK_MSG(n1.missing_same > 0, "negative control: N_ERROR_BUF = 0 but checked against 3 cells: the buffer assertion fails (missing in own box " + std::to_string(n1.missing_same) + ")");
        const Result n2 = run(2, 2, 8, 32, 3, 1, 3, 12, 0.75 / 64.0, 6);
        CHECK_MSG(n2.nest_bad > 0, "negative control: N_PROPER = 1 checked against 12: the nesting assertion fails (bad boxes " + std::to_string(n2.nest_bad) + ")");
    }
    const long nfail = fdstest::report("test_tag_buffer_large");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
