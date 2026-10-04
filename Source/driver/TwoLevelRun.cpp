// TwoLevelRun.cpp: see TwoLevelRun.H.
#include "TwoLevelRun.H"

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <cmath>
#include <iomanip>
#include <cstdio>
#include <sstream>

#include "DriverAdapter.H"
#include "Hierarchy.H"
#include "PostRegridProjection.H"
#include "PressureBackendSolver.H"
#include "RegistryTransfer.H"
#include "RegridAmrCore.H"
#include "TimeLoopWiring.H"

extern "C" void fds_p_set_legacy_save_guard(int on);   // fds_step.f90

namespace fdsamr {

namespace {

// Composite sums over the uncovered cells of all levels: volume integral of RHO and of RHO*ZZ_n.
std::vector<double> composite_mass(LevelRegistry& reg, int nl, int ns)
{
    std::vector<double> m(1 + ns, 0.0);
    for (int l = 0; l < nl; ++l) {
        const Level& lv = reg.level(l);
        Fields& F = reg.fields(l);
        const amrex::iMultiFab* cov = reg.covered_mask(l);
        const double vol = lv.dx[0] * lv.dx[1] * lv.dx[2];
        for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
            auto r = F["RHO"].const_array(mfi);
            auto z = F["ZZ"].const_array(mfi);
            amrex::Array4<const int> cv;
            if (cov) cv = cov->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (cov && cv(i, j, k) != 0) return;
                m[0] += r(i, j, k) * vol;
                for (int n = 0; n < ns; ++n) m[1 + n] += r(i, j, k) * z(i, j, k, n) * vol;
            });
        }
    }
    amrex::ParallelAllReduce::Sum(m.data(), static_cast<int>(m.size()), amrex::ParallelContext::CommunicatorSub());
    return m;
}

std::vector<fdsrt::ProjectionLevel> projection_levels(LevelRegistry& reg, int nl)
{
    std::vector<fdsrt::ProjectionLevel> lv;
    for (int l = 0; l < nl; ++l) {
        Fields& F = reg.fields(l);
        fdsrt::ProjectionLevel pl;
        pl.geom = reg.level(l).geom;
        pl.ref_ratio = l > 0 ? reg.level(l).ref_ratio_from_parent : amrex::IntVect(1);
        pl.vel = {&F["U"], &F["V"], &F["W"]};
        pl.D = &F["D"];
        pl.covered = reg.covered_mask(l);
        lv.push_back(pl);
    }
    return lv;
}

}  // namespace

int two_level_run(TimeLoop& loop, const Level0& l0, double dt_setup, const TwoLevelOptions& o)
{
    int fails = 0;
    LevelRegistry& reg = loop.registry();
    const int ns = l0.dom.n_total;
    auto probe = [&](const char* tag) {   // FDSTL_TLDIAG: periodic z ghost of U on level 0 at x face 5 (must equal the opposite valid row)
        if (!std::getenv("FDSTL_TLDIAG")) return;
        const amrex::MultiFab& u = reg.fields(0)["U"];
        for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
            auto a = u.const_array(mfi);
            amrex::Print() << "PROBE " << tag << ": U(5,0,-1) " << a(5, 0, -1) << " U(5,0,15) " << a(5, 0, 15) << " U(5,0,16) " << a(5, 0, 16) << " U(5,0,0) " << a(5, 0, 0) << "\n";
        }
    };
    probe("entry");
    const amrex::Box dom = l0.geom.Domain();

    // ---- the hierarchy: &AMR and &AMR_REGION text -> AmrParams, the level-0 meshes -> Hierarchy -> RegridAmrCore; level 1 is made by the first regrid_dynamic (the tag function tags the patch)
    amrex::IntVect rr(o.ratio, o.ratio, o.ratio);
    for (int d = 0; d < 3; ++d) if (dom.length(d) == 1) rr[d] = 1;
    const int plo[3] = {o.patch[0], dom.smallEnd(1), o.patch[2]}, phi_[3] = {o.patch[1], dom.bigEnd(1), o.patch[3]};
    amrex::Box patch(amrex::IntVect(plo[0], plo[1], plo[2]), amrex::IntVect(phi_[0], phi_[1], phi_[2]));
    std::ostringstream amr;
    amr << "&AMR MAX_LEVEL=1, REF_RATIO=" << o.ratio << ", BLOCKING_FACTOR=" << o.blocking << ", MAX_GRID_SIZE=" << o.maxsize << ", N_PROPER=1, N_ERROR_BUF=1, REGRID_INTERVAL=0, POST_REGRID_PROJECTION='" << o.projection << "' /\n";
    amr << std::setprecision(17) << "&AMR_REGION XB=";
    for (int d = 0; d < 3; ++d) amr << l0.geom.ProbLo(d) + plo[d] * l0.dx[d] << "," << l0.geom.ProbLo(d) + (phi_[d] + 1) * l0.dx[d] << (d < 2 ? "," : "");
    amr << ", LEVEL=1 /\n";
    fdsrt::Report rep;
    fdsrt::AmrParams ap = fdsrt::parse_amr_params(amr.str(), rep);
    std::vector<fdsrt::MeshInput> meshes;
    for (const MeshInfo& m : l0.mesh) {
        fdsrt::MeshInput mi;
        for (int d = 0; d < 3; ++d) mi.ijk[d] = m.ijk[d];
        for (int q = 0; q < 6; ++q) mi.xb[q] = m.xb[q];
        mi.rank = m.rank;
        meshes.push_back(mi);
    }
    std::array<bool, 3> per{l0.dom.periodic[0] != 0, l0.dom.periodic[1] != 0, l0.dom.periodic[2] != 0};
    fdsrt::Hierarchy h;
    if (!rep.ok() || !fdsrt::build_hierarchy_from_meshes(meshes, ap, per, h, rep) || !rep.ok()) {
        for (const auto& e : rep.errors) amrex::Print() << "TWO-LEVEL hierarchy error: " << e << "\n";
        return 1;
    }
    for (const auto& w : rep.warnings) amrex::Print() << "TWO-LEVEL hierarchy warning: " << w << "\n";
    fdsrt::RegridAmrCore core(h, ap);
    RegistryTransfer tr(reg, l0.dom.n_tracked);
    // The stage arrays a regrid does not transfer (D, DS, RSUM, MU, KRES, H, HS) of a new fine level: the parent cell value (piecewise constant). D is the projection target of the first regrid; RSUM, MU,
    // KRES are the equation-of-state and transport terms the first DIVERGENCE_PART_1 reads before the density update (which does not run in the first cycle, ICYC<=1) rebuilds them.
    tr.derive = [&](int level) {
        if (level < 1) return;
        const Level& lf = reg.level(level);
        amrex::BoxArray cba = lf.ba; cba.coarsen(lf.ref_ratio_from_parent);
        const amrex::IntVect ratio = lf.ref_ratio_from_parent;
        for (const char* nm : {"D", "DS", "RSUM", "MU", "KRES", "H", "HS"}) {
            if (!reg.fields(level).has(nm) || !reg.fields(level - 1).has(nm)) continue;
            amrex::MultiFab& Df = reg.fields(level)[nm];
            amrex::MultiFab ct(cba, lf.dm, Df.nComp(), 0);
            ct.setVal(0.0);
            ct.ParallelCopy(reg.fields(level - 1)[nm], 0, 0, Df.nComp(), 0, 0, reg.level(level - 1).geom.periodicity());
            for (amrex::MFIter mfi(Df); mfi.isValid(); ++mfi) {
                auto d = Df.array(mfi); auto c = ct.const_array(mfi); const int nc = Df.nComp();
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < nc; ++n) d(i, j, k, n) = c(amrex::coarsen(i, ratio[0]), amrex::coarsen(j, ratio[1]), amrex::coarsen(k, ratio[2]), n); });
            }
        }
    };
    core.set_data_transfer(&tr);
    core.set_tag_function([&](int lev, amrex::TagBoxArray& tags, amrex::Real) {
        if (lev != 0) return;
        for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
            auto a = tags.array(mfi);
            const amrex::Box b = mfi.validbox() & patch;
            if (b.ok()) amrex::LoopOnCpu(b, [&](int i, int j, int k) { a(i, j, k) = amrex::TagBox::SET; });
        }
    });
    core.init_static(reg, false, fdsrt::RegridAmrCore::DmFn(), &l0.dm);
    probe("after init_static");

    // ---- level 0 alone: one corrector-form pass for the D of the initial state (the target of the projection and the parent value of the new level's D)
    loop.set_state(0.0, dt_setup, 0);
    probe("before first prime");
    loop.prime_levels();
    probe("after first prime");
    fdsrt::PressureBackendSolver solver;
    fdsrt::PostRegridProjectionState pst;
    pst.solver = &solver;
    pst.D = [&](int l) { return &reg.fields(l)["D"]; };
    pst.options.accept_abs = 1.0e-9;
    fdsrt::install_post_regrid_projection(core, reg, pst);   // D-063: the hook of RegridAmrCore::regrid_dynamic
    core.regrid_dynamic(0.0);
    probe("after regrid_dynamic");
    const int finest = core.finestLevel();
    if (finest != 1) { amrex::Print() << "TWO-LEVEL: expected one fine level, got finest level " << finest << "\n"; return 1; }
    amrex::Print() << "TWO-LEVEL regrid_dynamic: changed " << core.last_outcome().changed << ", new cells " << core.last_outcome().new_cells << ", retained " << core.last_outcome().retained_cells
                   << ", post-regrid hook calls " << pst.calls;
    if (pst.calls > 0) amrex::Print() << " (max|div u - D| before " << pst.last.div_before << ", after " << pst.last.div_after << ", accepted " << pst.last.accepted << ")";
    amrex::Print() << "\n";
    for (int l = 1; l <= finest; ++l) loop.bind_level(l);
    fdsrt::CfHookStats cfs;
    fdsrt::install_cf_ghost_hooks(loop, finest, fdsrt::ThermoProvider(), &cfs);
    loop.set_flux_overwrite(o.overwrite);
    probe("after bind/hooks");
    amrex::Print() << "TWO-LEVEL: level 1 = " << reg.level(1).ba.size() << " box(es), " << reg.level(1).ba.numPts() << " cells, ratio (" << rr[0] << "," << rr[1] << "," << rr[2] << "), blocking factor "
                   << o.blocking << ", flux overwrite " << (o.overwrite ? "on" : "off") << ", projection " << o.projection << "\n";
    const int nl = finest + 1;
    if (o.boundary_test != 0) {
        // stage_boundary(1,3) and (1,6): the boundary routines FDS runs after the exchanges of the predictor and the corrector, on the bound fine level (fds_p_save_uvw on a fine box number, and with
        // FDSTL_SKIPAFT=0 also fill_omesh / VELOCITY_BC). Test 2 restores the pre-S14.1 guard (level-0 numbers only): the call must then stop the run. Fine-level calls must leave level 0 alone.
        if (o.boundary_test == 2) fds_p_set_legacy_save_guard(1);
        std::vector<amrex::MultiFab> keep;
        const char* vn[3] = {"U", "V", "W"};
        for (int d = 0; d < 3; ++d) { const amrex::MultiFab& m = reg.fields(0)[vn[d]]; keep.emplace_back(m.boxArray(), m.DistributionMap(), m.nComp(), m.nGrow()); amrex::MultiFab::Copy(keep.back(), m, 0, 0, m.nComp(), m.nGrow()); }
        amrex::Print() << "STAGE-BOUNDARY-FINE: calling stage_boundary(1,3) and stage_boundary(1,6)" << (o.boundary_test == 2 ? " with the pre-S14.1 level-0-only guard" : "") << "\n";
        loop.stage_boundary(1, 3);
        loop.stage_boundary(1, 6);
        double worst = 0.0;
        for (int d = 0; d < 3; ++d) {
            const amrex::MultiFab& m = reg.fields(0)[vn[d]];
            amrex::MultiFab diff(m.boxArray(), m.DistributionMap(), m.nComp(), m.nGrow());
            amrex::MultiFab::Copy(diff, m, 0, 0, m.nComp(), m.nGrow());
            amrex::MultiFab::Subtract(diff, keep[d], 0, 0, m.nComp(), m.nGrow());
            const double w = diff.norminf(0, m.nGrow());
            amrex::Print() << "STAGE-BOUNDARY-FINE: level 0 " << vn[d] << " (ghost layers included) changed by max |diff| " << w << "\n";
            worst = std::max(worst, w);
        }
        amrex::Print() << "STAGE-BOUNDARY-FINE: returned, level 0 U/V/W max change " << worst << (worst == 0.0 ? " (PASS)" : " (FAIL)") << "\n";
        return worst == 0.0 ? 0 : 1;
    }

    // ---- the fine level's own D, KRES, MU, ... by one corrector-form pass, then the post-regrid projection of the two-level velocity with that D (the body of the hook, Role 3 PostRegridProjection)
    loop.set_state(0.0, dt_setup, 0);
    loop.prime_levels();
    probe("after second prime");
    if (std::getenv("FDSTL_TLDIAG")) {
        for (int l = 0; l < nl; ++l)
            for (const char* nm : {"RHO", "TMP", "ZZ", "D", "DS", "DDDT", "KRES", "MU", "U", "V", "W"}) {
                if (!reg.fields(l).has(nm)) continue;
                const amrex::MultiFab& mf = reg.fields(l)[nm];
                amrex::Print() << "TLDIAG prime L" << l << " " << nm << " min " << mf.min(0, 0) << " max " << mf.max(0, 0) << (mf.contains_nan() ? " NAN" : "") << "\n";
            }
    }
    {
        std::vector<fdsrt::ProjectionLevel> lv = projection_levels(reg, nl);
        fdsrt::ProjectionOptions po = pst.options;
        po.enabled = (o.projection != "OFF");
        const fdsrt::ProjectionReport pr = fdsrt::project_after_regrid(lv, &solver, po);
        amrex::Print() << "TWO-LEVEL initial projection (" << o.projection << "): ran " << pr.ran << " solved " << pr.solved << " max|div u - D| before " << pr.div_before << " after " << pr.div_after
                       << " (level 0 " << (pr.div_after_level.empty() ? 0.0 : pr.div_after_level[0]) << ", level 1 " << (pr.div_after_level.size() > 1 ? pr.div_after_level[1] : 0.0) << "), accepted " << pr.accepted
                       << ", rhs mean/|rhs| " << pr.rhs_sum_rel << ", solver " << pr.solver << "\n";
        if (po.enabled && pr.ran && pr.solved) {
            for (int l = 0; l < nl; ++l) {
                Fields& F = reg.fields(l);
                const char* st[3] = {"US", "VS", "WS"}; const char* vn[3] = {"U", "V", "W"};
                for (int d = 0; d < 3; ++d) { amrex::MultiFab::Copy(F[st[d]], F[vn[d]], 0, 0, 1, 0); F[st[d]].FillBoundary(reg.level(l).geom.periodicity()); F[vn[d]].FillBoundary(reg.level(l).geom.periodicity()); }
            }
            loop.prime_levels();   // ghost faces of the projected velocity, D of the projected state
        }
        if (po.enabled && pr.ran && !pr.accepted) { ++fails; amrex::Print() << "TWO-LEVEL: the initial projection was not accepted\n"; }
    }

    // ---- the run
    const std::vector<double> m0 = composite_mass(reg, nl, ns);
    std::string logname = o.outdir + "/" + o.chid + "_two_level.csv";
    std::FILE* lf = amrex::ParallelDescriptor::IOProcessor() ? std::fopen(logname.c_str(), "w") : nullptr;
    if (lf) std::fprintf(lf, "step,t,dt,mass,mass_rel_change,divmax,div_l0,div_l1,removed_rel\n");
    double worst_mass = 0.0, worst_div = 0.0;
    std::vector<double> bylev;
    for (int s = 1; s <= o.steps; ++s) {
        if (!loop.advance()) { amrex::Print() << "TWO-LEVEL: STOP at step " << s << " (instability or non-finite DT)\n"; ++fails; break; }
        const std::vector<double> m = composite_mass(reg, nl, ns);
        const double dm = (m[0] - m0[0]) / m0[0];
        std::vector<fdsrt::ProjectionLevel> lv = projection_levels(reg, nl);
        const double dv = fdsrt::composite_divergence_error(lv, &bylev);
        worst_mass = std::max(worst_mass, std::abs(dm));
        worst_div = std::max(worst_div, dv);
        const double dts = loop.steps().empty() ? 0.0 : loop.steps().back().dt;
        if (lf) std::fprintf(lf, "%d,%.9g,%.9g,%.17g,%.6e,%.6e,%.6e,%.6e,%.3e\n", s, loop.time(), dts, m[0], dm, dv, bylev[0], bylev.size() > 1 ? bylev[1] : 0.0, loop.composite_report().removed_rel);
        if (o.log_every > 0 && (s % o.log_every == 0 || s == o.steps))
            amrex::Print() << "TWO-LEVEL step " << s << " t " << loop.time() << " dt " << dts << " composite mass change " << dm << " max|div u - D| " << dv << " (L0 " << bylev[0] << ", L1 " << (bylev.size() > 1 ? bylev[1] : 0.0) << ")\n";
    }
    if (lf) std::fclose(lf);
    const std::vector<double> m1 = composite_mass(reg, nl, ns);
    const fdsamr::TimeLoop::CompositeReport& cr = loop.composite_report();
    amrex::Print() << "TWO-LEVEL RESULT steps " << loop.steps().size() << " composite mass " << m0[0] << " -> " << m1[0] << " relative change " << (m1[0] - m0[0]) / m0[0] << " (worst over steps " << worst_mass << ")\n";
    for (int n = 0; n < ns; ++n)
        if (std::abs(m0[1 + n]) > 0.0) amrex::Print() << "TWO-LEVEL RESULT composite rho*Z" << n + 1 << " " << m0[1 + n] << " -> " << m1[1 + n] << " relative change " << (m1[1 + n] - m0[1 + n]) / m0[1 + n] << "\n";
    amrex::Print() << "TWO-LEVEL RESULT max|div u - D| over uncovered cells, worst over steps " << worst_div << "\n";
    amrex::Print() << "TWO-LEVEL RESULT composite pressure: " << cr.solves << " solves, backend " << cr.backend << ", removed mean (largest) " << cr.removed_mean << ", relative to rms " << cr.removed_rel
                   << ", fine H ghost cells set " << cr.h_ghost_cells << "\n";
    if (cr.res_checked > 0) amrex::Print() << "TWO-LEVEL RESULT pressure residual check: " << cr.res_checked << " solves checked, largest residual " << cr.res_check_max << ", largest limit " << cr.res_limit_max << " (round-off floor " << cr.res_floor_max << "), " << cr.res_failed << " above the limit\n";
    amrex::Print() << "TWO-LEVEL RESULT cf scalar ghost hook: calls " << cfs.calls << ", fine ghost cells " << cfs.fine_ghost_cells << ", covered cells " << cfs.covered_cells << ", conflicts " << cfs.conflicts << "\n";
    return fails;
}

}  // namespace fdsamr
