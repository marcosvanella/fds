// test_amrcore.cpp: R2a checks with AMReX (no FDS Fortran, no GPU): the AmrCore subclass builds the static two-level hierarchy of an
// ns2d_16_int_1to2_refinement-style input and notifies the driver's LevelRegistry; average-down conserves the integral to round-off; a uniform state
// stays uniform across the coarse-fine face after the ghost fill; ghost layers 1 and 2 hold the coarse value (also with ratio 4, periodic wrap and a
// neighbouring fine box). Usage: test_amrcore <cases dir>. Run on 1 rank (and 4 ranks by hand).
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <cstdio>
#include <string>

#include "DriverAdapter.H"
#include "LevelOps.H"
#include "RegridAmrCore.H"
#include "check.H"
#include "mesh_text.H"

using namespace fdsrt;

// Test-only mirror of fdsamr::exchange_fields (Source/driver/GhostExchange.cpp, which also holds the FDS-linked BcStep and cannot be linked without the Fortran
// objects). If the driver list changes, this copy must follow; the sentinel checks below do not depend on the exact list for H/HS and RHO/ZZ.
namespace fdsamr {
std::vector<std::string> exchange_fields(int code, bool predictor)
{
    switch (code) {
    case 1: return {"RHOS", "ZZS", "MU", "KRES", "D", "TMP", "RSUM"};
    case 4: return {"RHO", "ZZ", "MU", "KRES", "DS", "TMP", "RSUM"};
    case 3: return {"US", "VS", "WS", "HS"};
    case 6: return {"U", "V", "W", "H"};
    case 5: return {"FVX", "FVY", "FVZ", predictor ? "H" : "HS"};
    default: return {};
    }
}
}  // namespace fdsamr

namespace {

double fval(int i, int j, int k, int n) { return 1.0 + 0.001 * i + 0.0173 * j + 0.31 * k + 7.0 * n + std::sin(0.37 * i + 0.11 * k); }
double gval(int i, int j, int k, int n) { return 100.0 + 0.5 * i - 0.25 * j + 0.125 * k + 3.0 * n; }   // coarse values: distinct per cell

void fill_all(amrex::MultiFab& mf, double (*f)(int, int, int, int), bool with_ghost)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf.array(mfi);
        const amrex::Box b = with_ghost ? mfi.growntilebox() : mfi.validbox();
        for (int n = 0; n < mf.nComp(); ++n) amrex::LoopOnCpu(b, [&](int i, int j, int k) { a(i, j, k, n) = f(i, j, k, n); });
    }
}

// Direct test of the operators (no registry): coarse domain 16^3 periodic, fine boxes with given ratio.
void test_ops(int rr, bool two_boxes)
{
    const std::string tag = "ops ratio " + std::to_string(rr) + (two_boxes ? " two boxes" : " one box");
    amrex::Box cdom(amrex::IntVect(0), amrex::IntVect(15));
    amrex::RealBox rb({0, 0, 0}, {1, 1, 1});
    amrex::Array<int, 3> per{1, 1, 1};
    amrex::Geometry cg(cdom, rb, 0, per);
    amrex::IntVect ratio(rr);
    amrex::Geometry fg = amrex::refine(cg, ratio);
    amrex::BoxArray cba(cdom);
    cba.maxSize(8);
    amrex::DistributionMapping cdm(cba);
    // fine region = coarse cells (0..7, 4..11, 4..11): touches the low x edge, so the x ghost layers wrap through the periodic boundary
    amrex::BoxList bl;
    amrex::Box fregion = amrex::refine(amrex::Box(amrex::IntVect(0, 4, 4), amrex::IntVect(7, 11, 11)), ratio);
    if (two_boxes) {
        amrex::Box a = fregion, b = fregion;
        const int mid = fregion.smallEnd(0) + fregion.length(0) / 2 - 1;
        a.setBig(0, mid);
        b.setSmall(0, mid + 1);
        bl.push_back(a);
        bl.push_back(b);
    } else {
        bl.push_back(fregion);
    }
    amrex::BoxArray fba(bl);
    amrex::DistributionMapping fdm(fba);
    const int ncomp = 2, ng = 3;
    amrex::MultiFab C(cba, cdm, ncomp, 0), F(fba, fdm, ncomp, ng);
    fill_all(C, gval, false);
    F.setVal(-777.0);
    fill_all(F, fval, false);
    F.FillBoundary(fg.periodicity());   // same-level fill first (as the driver does)
    amrex::MultiFab F_before(fba, fdm, ncomp, ng);
    amrex::MultiFab::Copy(F_before, F, 0, 0, ncomp, ng);
    const long nf = fill_cf_ghosts_pc(F, C, fg, cg, ratio, 2, 0, ncomp);
    long bad_c = 0, bad_same = 0, bad_beyond = 0, n2 = 0;
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        auto b0 = F_before.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (vb.contains(iv)) return;
            const int layer = std::max({vb.smallEnd(0) - i, i - vb.bigEnd(0), vb.smallEnd(1) - j, j - vb.bigEnd(1), vb.smallEnd(2) - k, k - vb.bigEnd(2)});
            amrex::IntVect w = iv;   // wrapped position
            for (int d = 0; d < 3; ++d) { if (w[d] < 0) w[d] += fg.Domain().length(d); if (w[d] >= fg.Domain().length(d)) w[d] -= fg.Domain().length(d); }
            const bool in_fine = fba.contains(w);
            for (int n = 0; n < ncomp; ++n) {
                if (in_fine) { if (a(i, j, k, n) != b0(i, j, k, n)) ++bad_same; }
                else if (layer <= 2) {
                    const amrex::IntVect c = amrex::coarsen(iv, ratio);
                    amrex::IntVect cw = c;
                    for (int d = 0; d < 3; ++d) { if (cw[d] < 0) cw[d] += 16; if (cw[d] >= 16) cw[d] -= 16; }
                    if (a(i, j, k, n) != gval(cw[0], cw[1], cw[2], n)) ++bad_c;
                    if (n == 0) ++n2;
                } else if (a(i, j, k, n) != b0(i, j, k, n)) ++bad_beyond;
            }
        });
    }
    // layer 2 equals layer 1: the cell one step further out in a face-normal direction sits in the same coarse cell
    long bad_l12 = 0;
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        const amrex::IntVect probes[3] = {{vb.smallEnd(0), vb.smallEnd(1) + 1, vb.smallEnd(2) + 1}, {vb.smallEnd(0) + 1, vb.smallEnd(1), vb.smallEnd(2) + 1},
                                          {vb.smallEnd(0) + 1, vb.smallEnd(1) + 1, vb.smallEnd(2)}};
        for (int d = 0; d < 3; ++d) {
            amrex::IntVect p1 = probes[d], p2 = probes[d];
            p1[d] = vb.smallEnd(d) - 1; p2[d] = vb.smallEnd(d) - 2;
            amrex::IntVect w = p1; for (int e = 0; e < 3; ++e) { if (w[e] < 0) w[e] += fg.Domain().length(e); }
            if (fba.contains(w)) continue;   // a neighbouring fine box holds it
            for (int n = 0; n < ncomp; ++n) if (a(p1[0], p1[1], p1[2], n) != a(p2[0], p2[1], p2[2], n)) ++bad_l12;
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n2);   // ranks without a fine box have none
    CHECK_MSG(nf >= 0, tag);
    CHECK_MSG(n2 > 0, tag + ": some ghost cells were filled");
    CHECK_MSG(bad_c == 0, tag + ": ghost cells (layers 1,2) hold the value of the coarse cell that contains them, bad=" + std::to_string(bad_c));
    CHECK_MSG(bad_same == 0, tag + ": cells held by a neighbouring fine box are not overwritten");
    CHECK_MSG(bad_beyond == 0, tag + ": layer 3 is untouched");
    CHECK_MSG(bad_l12 == 0, tag + ": layer 2 equals layer 1");

    // average-down: covered coarse cells = mean of the fine cells, the rest untouched, integral conserved
    amrex::MultiFab C2(cba, cdm, ncomp, 0);
    fill_all(C2, gval, false);
    average_down_cells(F, C2, fg, cg, ratio, 0, ncomp);
    double sum_all = 0, sum_unc = 0, sum_fine = 0, maxdev = 0;
    long untouched_bad = 0, covered_n = 0;
    amrex::BoxArray covered = fba;
    covered.coarsen(ratio);
    for (amrex::MFIter mfi(C2); mfi.isValid(); ++mfi) {
        auto a = C2.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const bool cov = covered.contains(amrex::IntVect(i, j, k));
            for (int n = 0; n < ncomp; ++n) {
                sum_all += a(i, j, k, n);
                if (!cov) { sum_unc += a(i, j, k, n); if (a(i, j, k, n) != gval(i, j, k, n)) ++untouched_bad; }
                else {
                    double m = 0;
                    for (int c = 0; c < rr * rr * rr; ++c) m += fval(i * rr + c % rr, j * rr + (c / rr) % rr, k * rr + c / (rr * rr), n);
                    m /= rr * rr * rr;
                    maxdev = std::max(maxdev, std::abs(a(i, j, k, n) - m) / std::abs(m));
                }
            }
            if (cov) ++covered_n;
        });
    }
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < ncomp; ++n) sum_fine += a(i, j, k, n); });
    }
    amrex::ParallelDescriptor::ReduceRealSum(sum_all);
    amrex::ParallelDescriptor::ReduceRealSum(sum_unc);
    amrex::ParallelDescriptor::ReduceRealSum(sum_fine);
    amrex::ParallelDescriptor::ReduceRealMax(maxdev);
    amrex::ParallelDescriptor::ReduceLongSum(untouched_bad);
    amrex::ParallelDescriptor::ReduceLongSum(covered_n);
    const double vc = 1.0, vf = 1.0 / (rr * rr * rr);   // cell volumes in units of a coarse cell
    const double total_before = sum_unc * vc + sum_fine * vf, total_after = sum_all * vc;
    CHECK_MSG(covered_n == 8 * 8 * 8, tag + ": 512 covered coarse cells, got " + std::to_string(covered_n));
    CHECK_MSG(untouched_bad == 0, tag + ": uncovered coarse cells unchanged");
    CHECK_MSG(maxdev < 1e-14, tag + ": covered cell = mean of its fine cells, max rel dev " + std::to_string(maxdev));
    CHECK_MSG(std::abs(total_after - total_before) < 1e-13 * std::abs(total_before), tag + ": integral conserved, rel diff " + std::to_string(std::abs(total_after - total_before) / total_before));
}

}  // namespace

int main(int argc, char** argv)
{
    int one = 1;   // our only argument is the cases directory: do not let AMReX read it as an inputs file
    amrex::Initialize(one, argv);
    const std::string dir = argc > 1 ? argv[1] : ".";
    {
        test_ops(2, false);
        test_ops(4, false);
        test_ops(2, true);
        test_ops(4, true);

        // ---- the ns2d_16_int_1to2_refinement-style hierarchy through AmrCore and the driver's LevelRegistry ----
        const std::string text = fdsrt_test::read_file(dir + "/ns2d_16_int_1to2.fds");
        auto meshes = fdsrt_test::meshes_from_text(text);
        Report rep;
        AmrParams p = parse_amr_params(text, rep);
        Hierarchy h;
        const bool ok = build_hierarchy_from_meshes(meshes, p, {true, false, true}, h, rep);
        CHECK(ok && rep.ok());
        if (ok) {
            // level 0 as the driver would hold it: the level-0 boxes in order, FDS-style objects adopted by the registry
            amrex::Geometry g0 = make_level0_geometry(h);
            amrex::BoxList bl0;
            for (const GridBox& g : h.levels[0].grids) bl0.push_back(to_amrex(g.box));
            amrex::BoxArray ba0(bl0);
            amrex::DistributionMapping dm0(ba0);
            fdsamr::DomainInfo dom{};
            dom.periodic[0] = 1; dom.periodic[2] = 1;
            dom.n_tracked = 2; dom.n_total = 2; dom.nranks = amrex::ParallelDescriptor::NProcs();
            fdsamr::Level l0 = fdsamr::make_layout_level(0, g0, ba0, dm0, amrex::IntVect(1), dom);
            fdsamr::Fields F0(l0, 2);
            fdsamr::SideData sd0(l0, fdsamr::layout_cell_walls(l0));
            fdsamr::LevelRegistry reg(dom, 2);
            reg.adopt_level0(l0, F0, sd0);

            RegridAmrCore core(h, p);
            core.init_static(reg, false, RegridAmrCore::DmFn(), &dm0);
            CHECK(core.finestLevel() == 1);
            CHECK(reg.num_levels() == 2 && reg.has_level(1) && reg.n_make == 1);
            CHECK(core.boxArray(0) == ba0 && core.boxArray(1).size() == 1);
            CHECK(core.boxArray(1)[0] == amrex::Box(amrex::IntVect(8, 0, 8), amrex::IntVect(23, 0, 23)));
            CHECK(core.refRatio(0) == amrex::IntVect(2, 1, 2));
            CHECK(reg.level(1).ref_ratio_from_parent == amrex::IntVect(2, 1, 2));
            CHECK(std::abs(core.Geom(1).CellSize(0) - 0.5 * core.Geom(0).CellSize(0)) < 1e-15);
            CHECK(std::abs(core.Geom(1).CellSize(1) - core.Geom(0).CellSize(1)) < 1e-15);   // hidden direction is not refined
            CHECK(reg.fields(1)["RHO"].boxArray() == core.boxArray(1));
            const amrex::iMultiFab* cov = reg.covered_mask(0);
            long ncov = 0;
            if (cov) for (amrex::MFIter mfi(*cov); mfi.isValid(); ++mfi) { auto a = cov->const_array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { ncov += a(i, j, k); }); }
            amrex::ParallelDescriptor::ReduceLongSum(ncov);
            CHECK_MSG(ncov == 64, "covered coarse cells under the fine patch: 8 x 8, got " + std::to_string(ncov));

            const std::vector<std::string> names = {"RHO", "ZZ", "TMP", "MU"};
            // (1) uniform state stays uniform across the interface; fine level: valid cells 2.5, ghost cells start as the sentinel -777 (layer 3, H/HS must keep it)
            for (int l = 0; l < 2; ++l)
                for (const auto& n : reg.fields(l).names()) reg.fields(l)[n].setVal(-777.0);
            for (const auto& n : names) if (reg.fields(0).has(n)) reg.fields(0)[n].setVal(2.5);   // coarse level: everywhere (valid and ghost cells)
            for (const auto& n : names) if (reg.fields(1).has(n)) { amrex::MultiFab& m = reg.fields(1)[n]; for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) m[mfi].setVal<amrex::RunOn::Host>(2.5, mfi.validbox()); }
            fdsamr::CfGhostHook hook = make_cf_ghost_hook(reg);
            // code 4 (RHO, ZZ, ...) and code 3 (US.. and HS, which are not transferred)
            for (const auto& n : names) if (reg.fields(1).has(n)) reg.fields(1).fill_ghosts(n);
            fdsamr::CfGhostRequest rq{1, 4, true, fdsamr::exchange_fields(4, true)};
            hook(rq);
            fdsamr::CfGhostRequest rq3{1, 3, true, fdsamr::exchange_fields(3, true)};
            hook(rq3);
            long nonuni = 0, nlayer12 = 0, beyond_untouched = 1;
            for (const auto& n : {std::string("RHO"), std::string("ZZ")}) {
                const amrex::MultiFab& m = reg.fields(1)[n];
                const amrex::Box fdom = reg.level(1).geom.Domain();
                for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
                    auto a = m.const_array(mfi);
                    const amrex::Box vb = mfi.validbox();
                    amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                        const amrex::IntVect iv(i, j, k);
                        if (vb.contains(iv)) return;
                        const int layer = std::max({vb.smallEnd(0) - i, i - vb.bigEnd(0), vb.smallEnd(1) - j, j - vb.bigEnd(1), vb.smallEnd(2) - k, k - vb.bigEnd(2)});
                        const bool inside = fdom.contains(iv);   // periodic wrap: the patch is away from the x/z domain edges, y is not periodic
                        for (int c = 0; c < m.nComp(); ++c) {
                            if (inside && layer <= 2) { ++nlayer12; if (a(i, j, k, c) != 2.5) ++nonuni; }
                            else if (inside && layer == 3 && a(i, j, k, c) != -777.0) beyond_untouched = 0;
                            else if (!inside && a(i, j, k, c) != -777.0) beyond_untouched = 0;
                        }
                    });
                }
            }
            amrex::ParallelDescriptor::ReduceLongSum(nonuni);
            amrex::ParallelDescriptor::ReduceLongSum(nlayer12);
            CHECK_MSG(nlayer12 > 0 && nonuni == 0, "uniform state stays uniform in ghost layers 1 and 2 of RHO and ZZ across the coarse-fine faces, bad=" + std::to_string(nonuni));
            CHECK_MSG(beyond_untouched == 1, "layer 3 of RHO and cells outside the non-periodic y edge are not written");
            long hs_bad = 0;
            { const amrex::MultiFab& m = reg.fields(1)["HS"]; for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) { auto a = m.const_array(mfi); amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) { if (a(i, j, k) != -777.0) ++hs_bad; }); } }
            CHECK_MSG(hs_bad == 0, "H/HS ghost cells are not transferred by this hook");

            // (2) average-down through the registry: distinct fine values, coarse covered cells become the 2x1x2 mean
            for (int l = 0; l < 2; ++l) fill_all(reg.fields(l)["RHO"], l == 0 ? gval : fval, false);
            average_down_registry(reg, {"RHO"});
            double maxdev = 0;
            {
                const amrex::MultiFab& c = reg.fields(0)["RHO"];
                for (amrex::MFIter mfi(c); mfi.isValid(); ++mfi) {
                    auto a = c.const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox() & amrex::Box(amrex::IntVect(4, 0, 4), amrex::IntVect(11, 0, 11)), [&](int i, int j, int k) {
                        double m = 0;
                        for (int di = 0; di < 2; ++di) for (int dk = 0; dk < 2; ++dk) m += fval(2 * i + di, j, 2 * k + dk, 0);
                        m /= 4;
                        maxdev = std::max(maxdev, std::abs(a(i, j, k) - m) / std::abs(m));
                    });
                }
            }
            amrex::ParallelDescriptor::ReduceRealMax(maxdev);
            CHECK_MSG(maxdev < 1e-14, "registry average-down: covered cell = mean of its 4 fine cells, max rel dev " + std::to_string(maxdev));
        }
    }
    const long fails = fdstest::report("regrid_transport amrcore (R2a)");
    amrex::Finalize();
    return fails == 0 ? 0 : 1;
}
