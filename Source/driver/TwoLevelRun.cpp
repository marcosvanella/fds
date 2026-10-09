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
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include "TimeLoopWiring.H"

extern "C" void fds_p_set_legacy_save_guard(int on);   // fds_step.f90
extern "C" int fds_find_tracked_species(const char* name);  // fds_mesh_query.f90; 0-based ZZ index or -1

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

namespace {

// Cell-centred plot Multifab: RHO, U, V, W + undivided OMEGA_Y (ω_y = Δ_z U − Δ_x W).
// Owner viz: ω_y on mid-y xz + dark-grey mesh/AMR boxes. VORT (|ω|) stays TAG_VORT-only via ensure().
void write_rho_vel_plotfile(const std::string& path, LevelRegistry& reg, int nl, double t, int step)
{
    amrex::Vector<amrex::MultiFab> mf(nl);
    amrex::Vector<const amrex::MultiFab*> pmf(nl);
    amrex::Vector<amrex::Geometry> geom(nl);
    amrex::Vector<int> level_steps(nl, step);
    amrex::Vector<amrex::IntVect> ref_ratio(std::max(nl - 1, 0));
    for (int l = 0; l < nl; ++l) {
        const Level& Lv = reg.level(l);
        Fields& F = reg.fields(l);
        const amrex::Periodicity per = Lv.geom.periodicity();
        F["U"].FillBoundary(per);
        F["W"].FillBoundary(per);
        mf[l].define(Lv.ba, Lv.dm, 5, 0);
        amrex::MultiFab::Copy(mf[l], F["RHO"], 0, 0, 1, 0);
        amrex::Array<const amrex::MultiFab*, 3> faces{&F["U"], &F["V"], &F["W"]};
        amrex::average_face_to_cellcenter(mf[l], 1, faces);
        const amrex::Box dom = Lv.geom.Domain();
        const int dx_on = (dom.length(0) > 1) ? 1 : 0;
        const int dz_on = (dom.length(2) > 1) ? 1 : 0;
        const amrex::MultiFab& U = F["U"];
        const amrex::MultiFab& W = F["W"];
        for (amrex::MFIter mfi(mf[l]); mfi.isValid(); ++mfi) {
            auto a = mf[l].array(mfi);
            auto u = U.const_array(mfi);
            auto w = W.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                // undivided ω_y (FDS dump.f90 vorticity without /dx)
                a(i, j, k, 4) = dz_on * (u(i, j, k + 1) - u(i, j, k)) - dx_on * (w(i + 1, j, k) - w(i, j, k));
            });
        }
        pmf[l] = &mf[l];
        geom[l] = Lv.geom;
        if (l + 1 < nl) ref_ratio[l] = reg.level(l + 1).ref_ratio_from_parent;
    }
    amrex::Vector<std::string> names{"RHO", "U", "V", "W", "OMEGA_Y"};
    amrex::WriteMultiLevelPlotfile(path, nl, pmf, names, geom, t, level_steps, ref_ratio);
    amrex::Print() << "PLOTFILE wrote " << path << " (RHO,U,V,W,OMEGA_Y cell-centred; OMEGA_Y undivided ω_y, " << nl << " level(s))\n";
}

}  // namespace



namespace {

// Cell-centred undivided |ω| from face U,V,W (FDS dump.f90 vorticity without /dx).
// ωx = (W_{j+1}-W_j) - (V_{k+1}-V_k), cyclic. Thin directions contribute 0.
void fill_undivided_vorticity(LevelRegistry& reg, int nl)
{
    for (int l = 0; l < nl; ++l) {
        Fields& F = reg.fields(l);
        F.ensure("VORT");
        const amrex::Periodicity per = reg.level(l).geom.periodicity();
        F["U"].FillBoundary(per);
        F["V"].FillBoundary(per);
        F["W"].FillBoundary(per);
        const amrex::Box dom = reg.level(l).geom.Domain();
        const int dx_on = (dom.length(0) > 1) ? 1 : 0;
        const int dy_on = (dom.length(1) > 1) ? 1 : 0;
        const int dz_on = (dom.length(2) > 1) ? 1 : 0;
        amrex::MultiFab& Vo = F["VORT"];
        const amrex::MultiFab& U = F["U"];
        const amrex::MultiFab& V = F["V"];
        const amrex::MultiFab& W = F["W"];
        for (amrex::MFIter mfi(Vo); mfi.isValid(); ++mfi) {
            auto vo = Vo.array(mfi);
            auto u = U.const_array(mfi);
            auto v = V.const_array(mfi);
            auto w = W.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                const double ox = dy_on * (w(i, j + 1, k) - w(i, j, k)) - dz_on * (v(i, j, k + 1) - v(i, j, k));
                const double oy = dz_on * (u(i, j, k + 1) - u(i, j, k)) - dx_on * (w(i + 1, j, k) - w(i, j, k));
                const double oz = dx_on * (v(i + 1, j, k) - v(i, j, k)) - dy_on * (u(i, j + 1, k) - u(i, j, k));
                vo(i, j, k) = std::sqrt(ox * ox + oy * oy + oz * oz);
            });
        }
    }
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

    // ---- the hierarchy: &AMR / &AMR_REGION -> AmrParams, level-0 meshes -> Hierarchy -> RegridAmrCore; levels 1..max_level by regrid_dynamic (nested tags)
    if (o.max_level < 1 || o.max_level > 8) { amrex::Print() << "MULTI-LEVEL: max_level must be 1..8, got " << o.max_level << "\n"; return 1; }
    const char* tag = (o.max_level == 1) ? "TWO-LEVEL" : "MULTI-LEVEL";
    amrex::IntVect rr(o.ratio, o.ratio, o.ratio);
    for (int d = 0; d < 3; ++d) if (dom.length(d) == 1) rr[d] = 1;
    const int plo[3] = {o.patch[0], dom.smallEnd(1), o.patch[2]}, phi_[3] = {o.patch[1], dom.bigEnd(1), o.patch[3]};
    amrex::Box patch(amrex::IntVect(plo[0], plo[1], plo[2]), amrex::IntVect(phi_[0], phi_[1], phi_[2]));
    const int p2lo[3] = {o.patch2[0], dom.smallEnd(1), o.patch2[2]}, p2hi[3] = {o.patch2[1], dom.bigEnd(1), o.patch2[3]};
    amrex::Box patch2(amrex::IntVect(p2lo[0], p2lo[1], p2lo[2]), amrex::IntVect(p2hi[0], p2hi[1], p2hi[2]));
    if (o.max_level >= 2 && !patch.contains(patch2)) {
        amrex::Print() << tag << ": patch2 must lie inside patch (level-0 cells)\n";
        return 1;
    }
    std::ostringstream amr;
    amr << "&AMR MAX_LEVEL=" << o.max_level << ", REF_RATIO=" << o.ratio << ", BLOCKING_FACTOR=" << o.blocking << ", MAX_GRID_SIZE=" << o.maxsize
        << ", N_PROPER=1, N_ERROR_BUF=1, REGRID_INTERVAL=" << o.regrid_every << ", POST_REGRID_PROJECTION='" << o.projection << "'";
    // Adaptive TAG_*: vortex-ring TAG_VORT path, else blob TAG_SPEC_*. Static (regrid_every==0) uses nested set_tag_function only.
    const bool use_tag_vort = o.tag_vort_auto || !std::isnan(o.tag_vort);
    double tag_vort_thr = o.tag_vort;
    if (o.regrid_every > 0 && use_tag_vort) {
        fill_undivided_vorticity(reg, 1);   // L0 |ω| before AmrParams (auto thr needs max)
        if (o.tag_vort_auto) {
            // Owner overnight: 0.45*max hugs ring core tighter than 0.35 (less fat L1).
            tag_vort_thr = 0.45 * reg.fields(0)["VORT"].max(0, 0);   // MultiFab::max is global
            amrex::Print() << tag << ": TAG_VORT auto thr = 0.45 * max|VORT|_L0 = " << tag_vort_thr << "\n";
        }
        amr << ", TAG_VORT=" << std::setprecision(17) << tag_vort_thr << ", TAG_KEEP=0.85";
    } else if (o.regrid_every > 0) {
        amr << ", TAG_SPEC_ID='TRACER', TAG_SPEC_DIFF=0.05, TAG_KEEP=0.8";
    }
    amr << " /\n";
    // Vortex adaptive: whole-domain refinable region so the ring can leave the static patch.
    amr << std::setprecision(17) << "&AMR_REGION XB=";
    if (use_tag_vort) {
        for (int d = 0; d < 3; ++d) amr << l0.geom.ProbLo(d) << "," << l0.geom.ProbHi(d) << (d < 2 ? "," : "");
    } else {
        for (int d = 0; d < 3; ++d) amr << l0.geom.ProbLo(d) + plo[d] * l0.dx[d] << "," << l0.geom.ProbLo(d) + (phi_[d] + 1) * l0.dx[d] << (d < 2 ? "," : "");
    }
    amr << ", LEVEL=" << o.max_level << " /\n";
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
        for (const auto& e : rep.errors) amrex::Print() << tag << " hierarchy error: " << e << "\n";
        return 1;
    }
    for (const auto& w : rep.warnings) amrex::Print() << tag << " hierarchy warning: " << w << "\n";
    fdsrt::RegridAmrCore core(h, ap);
    RegistryTransfer tr(reg, l0.dom.n_tracked);
    // Stage arrays a regrid does not transfer (D, DS, RSUM, MU, KRES, H, HS) of a new fine level: parent cell value (piecewise constant), any level >= 1 (D-053).
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
    // D-090: map TAG_SPEC_ID -> ZZ component before set_tag_fields (no hardcode). Abort if the name is unknown.
    int zz_comp_tag = 0;
    if (!ap.tag_spec_id.empty()) {
        const int idx = fds_find_tracked_species(ap.tag_spec_id.c_str());
        if (idx < 0) {
            amrex::Abort("D-090: TAG_SPEC_ID='" + ap.tag_spec_id + "' is not a tracked species or tracer (SPECIES_MIXTURE%ID)");
        }
        zz_comp_tag = idx;
        amrex::Print() << tag << ": TAG_SPEC_ID='" << ap.tag_spec_id << "' -> ZZ component " << zz_comp_tag << "\n";
    }
    // Registry fields for ErrorEst TAG_* (preferred when AmrParams has a criterion). Fallback nested set_tag_function for static.
    core.set_tag_fields([&](int lev) {
        Fields& F = reg.fields(lev);
        const amrex::Periodicity per = reg.level(lev).geom.periodicity();
        fdsrt::TagFields f;
        if (F.has("TMP")) { F["TMP"].FillBoundary(per); f.tmp = &F["TMP"]; }
        if (F.has("RHO")) { F["RHO"].FillBoundary(per); f.rho = &F["RHO"]; }
        if (F.has("ZZ")) {
            F["ZZ"].FillBoundary(per);
            f.zz = &F["ZZ"];
            f.zz_comp = zz_comp_tag;   // D-090: from TAG_SPEC_ID via fds_find_tracked_species
        }
        // D-090: TAG_HRRPUV field — HRRPUV if registered, else FDS Q (W/m3; thr×1000 in apply_tag_criteria)
        if (F.has("HRRPUV")) f.hrrpuv = &F["HRRPUV"];
        else if (F.has("Q")) f.hrrpuv = &F["Q"];
        if (F.has("VORT")) { F["VORT"].FillBoundary(per); f.vort = &F["VORT"]; }
        for (int d = 0; d < 3; ++d) f.dirs[d] = (dom.length(d) == 1) ? 0 : 1;
        return f;
    });
    // Nested tags: level 0 tags `patch` (creates level 1); levels 1..max_level-1 tag the refined image of `patch2` (creates levels 2..max_level).
    core.set_tag_function([&](int lev, amrex::TagBoxArray& tags, amrex::Real) {
        if (lev >= o.max_level) return;
        amrex::IntVect cum(1, 1, 1);
        for (int l = 0; l < lev; ++l) for (int d = 0; d < 3; ++d) cum[d] *= rr[d];
        const amrex::Box src = (lev == 0) ? patch : patch2;
        amrex::Box tagbox(amrex::IntVect(src.smallEnd(0) * cum[0], src.smallEnd(1) * cum[1], src.smallEnd(2) * cum[2]),
                          amrex::IntVect(src.bigEnd(0) * cum[0] + cum[0] - 1, src.bigEnd(1) * cum[1] + cum[1] - 1, src.bigEnd(2) * cum[2] + cum[2] - 1));
        for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
            auto a = tags.array(mfi);
            const amrex::Box b = mfi.validbox() & tagbox;
            if (b.ok()) amrex::LoopOnCpu(b, [&](int i, int j, int k) { a(i, j, k) = amrex::TagBox::SET; });
        }
    });
    core.init_static(reg, false, fdsrt::RegridAmrCore::DmFn(), &l0.dm);
    probe("after init_static");

    // ---- level 0 alone: one corrector-form pass for the D of the initial state (projection target and parent value of new levels' D)
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
    loop.zero_closed_wall_normal();   // IC: Neumann faces carry no net volume flux into the first projection
    // One AmrCore::regrid grows at most one new fine level; repeat until max_level (or stuck).
    for (int pass = 0; pass < o.max_level; ++pass) {
        if (use_tag_vort) fill_undivided_vorticity(reg, reg.num_levels());
        core.regrid_dynamic(0.0);
        if (core.finestLevel() >= o.max_level) break;
        if (pass > 0 && !core.last_outcome().changed) break;
    }
    probe("after regrid_dynamic");
    const int finest = core.finestLevel();
    if (finest != o.max_level) {
        amrex::Print() << tag << ": expected finest level " << o.max_level << ", got " << finest << " (tagging / blocking / region)\n";
        return 1;
    }
    amrex::Print() << tag << " regrid_dynamic: changed " << core.last_outcome().changed << ", new cells " << core.last_outcome().new_cells << ", retained " << core.last_outcome().retained_cells
                   << ", post-regrid hook calls " << pst.calls;
    if (pst.calls > 0) amrex::Print() << " (max|div u - D| before " << pst.last.div_before << ", after " << pst.last.div_after << ", accepted " << pst.last.accepted << ")";
    amrex::Print() << "\n";
    for (int l = 1; l <= finest; ++l) loop.bind_level(l);   // mesh offset + derive of parent stage arrays already done by RegistryTransfer::derive at make_level
    fdsrt::CfHookStats cfs;
    fdsrt::install_cf_ghost_hooks(loop, finest, fdsrt::ThermoProvider(), &cfs);
    loop.set_flux_overwrite(o.overwrite);
    probe("after bind/hooks");
    for (int l = 1; l <= finest; ++l)
        amrex::Print() << tag << ": level " << l << " = " << reg.level(l).ba.size() << " box(es), " << reg.level(l).ba.numPts() << " cells, ratio from parent ("
                       << reg.level(l).ref_ratio_from_parent[0] << "," << reg.level(l).ref_ratio_from_parent[1] << "," << reg.level(l).ref_ratio_from_parent[2] << ")\n";
    amrex::Print() << tag << ": levels 0.." << finest << ", blocking " << o.blocking << ", flux overwrite " << (o.overwrite ? "on" : "off") << ", projection " << o.projection
                   << ", bind_only " << (o.bind_only ? "yes" : "no") << "\n";
    const int nl = finest + 1;
    if (o.bind_only) {
        amrex::Print() << tag << " BIND-ONLY PASS: registry has " << reg.num_levels() << " level(s), TimeLoop bound " << loop.num_levels() << "\n";
        return (reg.num_levels() == nl && loop.num_levels() == nl) ? 0 : 1;
    }
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
        loop.zero_closed_wall_normal();   // fine faces prolonged onto INERT/MIRROR, then the projection
        const fdsrt::ProjectionReport pr = fdsrt::project_after_regrid(lv, &solver, po);
        amrex::Print() << tag << " initial projection (" << o.projection << "): ran " << pr.ran << " solved " << pr.solved << " max|div u - D| before " << pr.div_before << " after " << pr.div_after
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
        if (po.enabled && pr.ran && !pr.accepted) { ++fails; amrex::Print() << tag << ": the initial projection was not accepted\n"; }
    }

    // ---- the run
    const std::vector<double> m0 = composite_mass(reg, nl, ns);
    std::string logname = o.outdir + "/" + o.chid + "_multi_level.csv";
    std::FILE* lf = amrex::ParallelDescriptor::IOProcessor() ? std::fopen(logname.c_str(), "w") : nullptr;
    if (lf) { std::fprintf(lf, "step,t,dt,mass,mass_rel_change,divmax"); for (int l = 0; l < nl; ++l) std::fprintf(lf, ",div_l%d", l); std::fprintf(lf, ",removed_rel\n"); }
    double worst_mass = 0.0, worst_div = 0.0;
    std::vector<double> bylev;
    int nl_now = nl;
    for (int s = 1; s <= o.steps; ++s) {
        if (!loop.advance()) { amrex::Print() << tag << ": STOP at step " << s << " (instability or non-finite DT)\n"; ++fails; break; }
        // Mid-run adaptive regrid (Role 3 call shape): gate on AmrParams.regrid_interval; unbind top-down via finest_bound;
        // regrid_dynamic brackets LevelListener — do NOT call reg.begin_regrid/end_regrid here.
        {
            const int ri = core.amr_params().regrid_interval;
            if (ri > 0 && s % ri == 0) {
                const int fb = loop.finest_bound();
                for (int l = fb; l >= 1; --l) loop.unbind_level(l);
                if (use_tag_vort) fill_undivided_vorticity(reg, reg.num_levels());
                // auto: thr stays 0.45*max|VORT|_L0 of now (AmrParams thr fixed, so rescale VORT by thr/(0.45*max))
                if (o.tag_vort_auto) { const double sc = tag_vort_thr / (0.45 * reg.fields(0)["VORT"].max(0, 0)); for (int l = 0; l < reg.num_levels(); ++l) reg.fields(l)["VORT"].mult(sc); }
                const bool changed = core.regrid_dynamic(loop.time());
                const int finest = core.finestLevel();
                bool rebound = false;
                for (int l = 1; l <= finest; ++l) { loop.bind_level(l); rebound = true; }
                if (finest >= 1) {
                    fdsrt::install_cf_ghost_hooks(loop, finest, fdsrt::ThermoProvider(), &cfs);
                    loop.set_flux_overwrite(o.overwrite);
                }
                if (changed || rebound) loop.prime_levels();
                nl_now = finest + 1;
                amrex::Print() << tag << " mid-regrid after step " << s << " t " << loop.time()
                               << " changed " << changed << " finest " << finest
                               << " boxes L1 " << (finest >= 1 ? reg.level(1).ba.size() : 0) << "\n";
                // Mid-run dump shares the same write helper (Architect: regrid dump path); only without a --plot-every cadence (else it dumped every regrid).
                if (!o.plotfile.empty() && o.plot_every <= 0) {
                    const std::string plt = o.plotfile + "/plt_" + o.chid + "_s" + std::to_string(s);
                    write_rho_vel_plotfile(plt, reg, nl_now, loop.time(), loop.icyc());
                }
            }
        }
        // Owner-show cadence (D-089 pack): --plot-every N with --plotfile; dump when s % N == 0 (1-based step). Final dump unchanged; mid-regrid dump off.
        if (o.plot_every > 0 && !o.plotfile.empty() && (s % o.plot_every) == 0) {
            const std::string plt = o.plotfile + "/plt_" + o.chid + "_s" + std::to_string(s);
            write_rho_vel_plotfile(plt, reg, nl_now, loop.time(), loop.icyc());
        }
        const std::vector<double> m = composite_mass(reg, nl_now, ns);
        const double dm = (m[0] - m0[0]) / m0[0];
        std::vector<fdsrt::ProjectionLevel> lv = projection_levels(reg, nl_now);
        const double dv = fdsrt::composite_divergence_error(lv, &bylev);
        worst_mass = std::max(worst_mass, std::abs(dm));
        worst_div = std::max(worst_div, dv);
        const double dts = loop.steps().empty() ? 0.0 : loop.steps().back().dt;
        // CLIP_MASS / CLIP_CELLS (P3-C11/C12): every dump step when P3-C01 CSV path is active (emit zeros too).
        if (lf) {
            const int ntr = l0.dom.n_tracked;
            for (int l = 0; l < nl_now; ++l)
                for (int n = 1; n <= ntr; ++n) {
                    amrex::Print() << "CLIP_MASS N=" << n << " L=" << l << " mass=" << std::setprecision(17) << loop.clip_mass(n, l) << "\n";
                    amrex::Print() << "CLIP_CELLS N=" << n << " L=" << l << " count=" << loop.clip_cells(n, l) << "\n";
                }
            std::fprintf(lf, "%d,%.9g,%.9g,%.17g,%.6e,%.6e", s, loop.time(), dts, m[0], dm, dv);
            for (int l = 0; l < nl_now; ++l) std::fprintf(lf, ",%.6e", l < static_cast<int>(bylev.size()) ? bylev[l] : 0.0);
            std::fprintf(lf, ",%.3e\n", loop.composite_report().removed_rel);
        }
        if (o.log_every > 0 && (s % o.log_every == 0 || s == o.steps)) {
            amrex::Print() << tag << " step " << s << " t " << loop.time() << " dt " << dts << " composite mass change " << dm << " max|div u - D| " << dv;
            for (int l = 0; l < nl_now; ++l) amrex::Print() << " L" << l << " " << (l < static_cast<int>(bylev.size()) ? bylev[l] : 0.0);
            amrex::Print() << "\n";
        }
    }
    if (lf) std::fclose(lf);
    loop.write_perf();
    if (!o.plotfile.empty() && fails == 0) {
        const std::string plt = o.plotfile + "/plt_" + o.chid + "_final";
        write_rho_vel_plotfile(plt, reg, loop.num_levels(), loop.time(), loop.icyc());
    }
    const std::vector<double> m1 = composite_mass(reg, loop.num_levels(), ns);
    const fdsamr::TimeLoop::CompositeReport& cr = loop.composite_report();
    amrex::Print() << tag << " RESULT steps " << loop.steps().size() << " composite mass " << m0[0] << " -> " << m1[0] << " relative change " << (m1[0] - m0[0]) / m0[0] << " (worst over steps " << worst_mass << ")\n";
    for (int n = 0; n < ns; ++n)
        if (std::abs(m0[1 + n]) > 0.0) amrex::Print() << tag << " RESULT composite rho*Z" << n + 1 << " " << m0[1 + n] << " -> " << m1[1 + n] << " relative change " << (m1[1 + n] - m0[1 + n]) / m0[1 + n] << "\n";
    amrex::Print() << tag << " RESULT max|div u - D| over uncovered cells, worst over steps " << worst_div << "\n";
    amrex::Print() << tag << " RESULT composite pressure: " << cr.solves << " solves, backend " << cr.backend << ", removed mean (largest) " << cr.removed_mean << ", relative to rms " << cr.removed_rel
                   << ", fine H ghost cells set " << cr.h_ghost_cells << "\n";
    if (cr.res_checked > 0) amrex::Print() << tag << " RESULT pressure residual check: " << cr.res_checked << " solves checked, largest residual " << cr.res_check_max << ", largest limit " << cr.res_limit_max << " (round-off floor " << cr.res_floor_max << "), " << cr.res_failed << " above the limit\n";
    amrex::Print() << tag << " RESULT cf scalar ghost hook: calls " << cfs.calls << ", fine ghost cells " << cfs.fine_ghost_cells << ", covered cells " << cfs.covered_cells << ", conflicts " << cfs.conflicts << "\n";
    return fails;
}

}  // namespace fdsamr
