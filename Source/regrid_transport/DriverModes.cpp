// DriverModes.cpp: Role 3 driver-level end-to-end modes (R2b); see DriverModes.H. Built and run with gfortran and with oneAPI (ifx), patches 0005-0009 as validated.
//
// Mode "transport": prescribed constant velocity (no pressure solve, no velocity update), two species of equal molecular weight (the second one is a passive tracer blob), uniform density
// and temperature. One step is the FDS stage sequence of TimeLoop::advance() with the pressure and velocity parts left out, run through the per-level stage entry points, with Role 3's
// FluxStageRunner between the pieces (D-061 order, notes/flux-stage-wiring.md):
//   predictor: state, viscosity (all levels) | ADV read-out (all levels), runner ADV | density (finest first) | exchange 1 + boundary 1 (the coarse-fine ghost hook runs in exchange) |
//              velocity flux, init divergence, WALL_BC | DIF read-out (all levels), runner DIF | DIVERGENCE_PART_1 with the DIF overrides (finest first)
//   corrector: the same with exchange / boundary code 4. US/VS/WS are set equal to U/V/W (the velocity is prescribed).
// Mode "ghost": the FR-016 ghost check through a real stage_exchange (see mode_ghost).
// Output lines start with "RTE2E" and are read by tests/run_e2e_driver.sh.
#include "DriverModes.H"

#include <AMReX.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "DriverAdapter.H"
#include "FluxStageRunner.H"
#include "GhostShare.H"
#include "RegistryTransfer.H"
#include "TimeLoopWiring.H"
#include "tests/dump_reader.H"

namespace fdsrt {
namespace {

using fdsamr::Fields;
using fdsamr::Level;
using fdsamr::LevelRegistry;
using fdsamr::RegistryTransfer;
using fdsamr::TimeLoop;

struct Opt {
    std::string test = "transport";
    std::string mode = "patch";        // patch | full | none
    int p[4] = {4, 11, 4, 11};         // coarse cell range x0 x1 z0 z1 of the level-1 patch
    int py[2] = {-1, -1};              // coarse cell range y0 y1 of the patch (default: the whole y extent; a 3D case with an interior range keeps the fine box off every domain edge)
    int steps = 6;
    bool overwrite = true;
    double u0 = 1.0, w0 = 0.5;
    double dt = 0.0;                   // 0: CFL 0.24 on the finest level
    std::vector<std::array<double, 4>> blobs;   // physical x0 x1 z0 z1 of tracer blobs (aligned with coarse cell faces)
    std::string dump, label = "run", dumps_prefix;
    int maxsize = 16, ratio = 2, step = 3;
    bool no_blob = false;
    double smooth = 0.0;               // > 0: smooth tracer profile 0.5 + smooth*sin(2 pi x) cos(2 pi z) instead of blobs (no clipping, no limiter undershoot)
    int kinds = 3;                     // bit 0: overwrite ADV fluxes, bit 1: overwrite DIF fluxes (diagnostic; default both)
};


Opt parse(int argc, char** argv)
{
    Opt o;
    o.test = argc > 3 ? argv[3] : "list";
    for (int i = 4; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](int n) { if (i + n >= argc) amrex::Abort("rt-e2e: missing value after " + a); };
        if (a == "--patch") { need(4); for (int q = 0; q < 4; ++q) o.p[q] = std::atoi(argv[++i]); o.mode = "patch"; }
        else if (a == "--py") { need(2); o.py[0] = std::atoi(argv[++i]); o.py[1] = std::atoi(argv[++i]); }
        else if (a == "--full") o.mode = "full";
        else if (a == "--none") o.mode = "none";
        else if (a == "--steps") { need(1); o.steps = std::atoi(argv[++i]); }
        else if (a == "--overwrite") { need(1); o.overwrite = std::atoi(argv[++i]) != 0; }
        else if (a == "--velocity") { need(2); o.u0 = std::atof(argv[++i]); o.w0 = std::atof(argv[++i]); }
        else if (a == "--dt") { need(1); o.dt = std::atof(argv[++i]); }
        else if (a == "--blob") { need(4); std::array<double, 4> b; for (int q = 0; q < 4; ++q) b[q] = std::atof(argv[++i]); o.blobs.push_back(b); }
        else if (a == "--no-blob") o.no_blob = true;
        else if (a == "--smooth") { need(1); o.smooth = std::atof(argv[++i]); }
        else if (a == "--kinds") { need(1); o.kinds = std::atoi(argv[++i]); }
        else if (a == "--dump") { need(1); o.dump = argv[++i]; }
        else if (a == "--label") { need(1); o.label = argv[++i]; }
        else if (a == "--maxsize") { need(1); o.maxsize = std::atoi(argv[++i]); }
        else if (a == "--ratio") { need(1); o.ratio = std::atoi(argv[++i]); }
        else if (a == "--dumps") { need(1); o.dumps_prefix = argv[++i]; }
        else if (a == "--step") { need(1); o.step = std::atoi(argv[++i]); }
        else amrex::Abort("rt-e2e: unknown option " + a);
    }
    return o;
}

// Everything one transport run needs.
struct Env {
    TimeLoop& loop;
    const fdsamr::Level0& l0;
    Opt o;
    LevelRegistry& reg;
    int nl = 1, ns = 0;
    std::unique_ptr<RegistryTransfer> tr;
    std::unique_ptr<fdsamr::FluxStages> fs;
    std::unique_ptr<FluxStageRunner> runner;
    CfHookStats cf_stats;
    double max_flux[2] = {0.0, 0.0};   // largest |face flux| seen at the override positions: ADV, DIF (shows that the DIF overrides carry nonzero values)
    double t = 0.0, dt = 0.0;
    bool in_pred = true;
    int icyc = 1;   // FDS skips the density update while ICYC<=1 (DENSITY_PRE_CLIP), so the first step is cycle 2
    Env(TimeLoop& l, const fdsamr::Level0& z, const Opt& op) : loop(l), l0(z), o(op), reg(l.registry()) {}
};

void set_velocity(Env& e, int lev, double u, double w)
{
    Fields& F = e.reg.fields(lev);
    const char* nm[6] = {"U", "V", "W", "US", "VS", "WS"};
    const double val[6] = {u, 0.0, w, u, 0.0, w};
    for (int q = 0; q < 6; ++q) if (F.has(nm[q])) F[nm[q]].setVal(val[q], 0, F[nm[q]].nComp(), F[nm[q]].nGrowVect());
}

bool in_blob(const Opt& o, double x, double z)
{
    for (const auto& b : o.blobs) if (x > b[0] && x < b[1] && z > b[2] && z < b[3]) return true;
    return false;
}

// tracer (species index 1): blobs or a smooth profile on `lev`, background elsewhere; RHO and TMP stay as the FDS set-up made them (equal molecular weights: same RSUM, same density)
void init_profile(Env& e, int lev)
{
    Fields& F = e.reg.fields(lev);
    const Level& L = e.reg.level(lev);
    const double plo[3] = {L.geom.ProbLo(0), L.geom.ProbLo(1), L.geom.ProbLo(2)};
    const double TWO_PI = 6.283185307179586;
    for (const char* zn : {"ZZ", "ZZS"}) {
        for (amrex::MFIter mfi(F[zn]); mfi.isValid(); ++mfi) {
            auto z = F[zn].array(mfi);
            const amrex::Box gb = amrex::grow(mfi.validbox(), 0);
            amrex::LoopOnCpu(gb, [&](int i, int j, int k) {
                const double x = plo[0] + (i + 0.5) * L.dx[0], zc = plo[2] + (k + 0.5) * L.dx[2];
                double t = (!e.o.no_blob && in_blob(e.o, x, zc)) ? 1.0 : 0.0;
                if (e.o.smooth > 0.0 && !e.o.no_blob) t = 0.5 + e.o.smooth * std::sin(TWO_PI * x) * std::cos(TWO_PI * zc);
                for (int n = 0; n < e.ns; ++n) z(i, j, k, n) = (n == 1) ? t : (n == 0 ? 1.0 - t : 0.0);
            });
        }
        F.fill_ghosts(zn);
    }
}

void make_level1(Env& e, bool transfer = true)
{
    const Level& L0 = e.reg.level(0);
    const amrex::Box dom = e.l0.geom.Domain();
    amrex::Box cb = dom;
    if (e.o.mode == "patch") cb = amrex::Box(amrex::IntVect(e.o.p[0], e.o.py[0] >= 0 ? e.o.py[0] : dom.smallEnd(1), e.o.p[2]), amrex::IntVect(e.o.p[1], e.o.py[0] >= 0 ? e.o.py[1] : dom.bigEnd(1), e.o.p[3]));
    const amrex::IntVect rr(e.o.ratio, dom.length(1) == 1 ? 1 : e.o.ratio, e.o.ratio);
    fdsrt::LevelLayout fl;
    fl.level = 1;
    fl.ref_ratio_from_parent = rr;
    fl.geom = amrex::refine(e.l0.geom, rr);
    fl.ba = amrex::BoxArray(amrex::refine(cb, rr));
    fl.ba.maxSize(amrex::IntVect(e.o.maxsize, 100000, e.o.maxsize));
    fl.dm.define(fl.ba);
    (void)L0;
    e.reg.begin_regrid();
    e.reg.make_level(fl);
    if (transfer) e.tr->fill_new_level(fl);
    e.reg.end_regrid();
    if (transfer) e.tr->hierarchy_done(false);
}

// ---- stage pieces over all levels
bool trace() { static const bool t = std::getenv("RTE2E_TRACE") != nullptr; return t; }
#define TR(msg) do { if (trace()) amrex::Print() << "RTE2E-TRACE " << msg << "\n"; } while (0)
template <class F> void up(Env& e, F&& f) { for (int l = 0; l < e.nl; ++l) { TR("  level " << l); f(l); } }
template <class F> void down(Env& e, F&& f) { for (int l = e.nl - 1; l >= 0; --l) { TR("  level " << l); f(l); } }

double flux_max(Env& e, FluxKind k)
{
    double m = 0.0;
    for (int l = 0; l < e.nl; ++l) for (int d = 0; d < 3; d += 2) {
        const double v = static_cast<double>(e.fs->stage_flux(l, k, d).norminf(0, e.ns, amrex::IntVect(0), true));
        TR("flux " << (k == FluxKind::Adv ? "ADV" : "DIF") << " level " << l << " dir " << d << " max " << v << "; U max " << e.reg.fields(l)["U"].norminf(0, 1, amrex::IntVect(0), true) << " US max " << e.reg.fields(l)["US"].norminf(0, 1, amrex::IntVect(0), true) << " FX max " << (e.reg.fields(l).has("FX") ? e.reg.fields(l)["FX"].norminf(0, 1, amrex::IntVect(0), true) : -1.0));
        if (trace()) { const auto& mf = e.fs->stage_flux(l, k, d); double sm = 0, sa = 0; long n = 0; for (amrex::MFIter mi(mf); mi.isValid(); ++mi) { const auto a = mf.const_array(mi); amrex::Box bx = mi.fabbox(); amrex::LoopOnCpu(bx, [&](int i, int j, int kk) { for (int q = 0; q < e.ns; ++q) { sa += std::abs(a(i, j, kk, q)); sm = std::max(sm, a(i, j, kk, q)); } ++n; }); } TR("  scan fabbox: n " << n << " sumabs " << sa << " max " << sm << " ncomp " << mf.nComp() << " ngrow " << mf.nGrow() << " type " << mf.ixType().toIntVect()); }
        m = std::max(m, v);
    }
    amrex::ParallelAllReduce::Max(m, amrex::ParallelContext::CommunicatorSub());
    return m;
}

void set_ov(Env& e, FluxKind k)
{
    TR("set_overrides");
    e.max_flux[k == FluxKind::Adv ? 0 : 1] = std::max(e.max_flux[k == FluxKind::Adv ? 0 : 1], flux_max(e, k));
    const bool want = e.o.overwrite && (e.o.kinds & (k == FluxKind::Adv ? 1 : 2)) && !(std::getenv("RTE2E_SKIP_STAGE") && std::atoi(std::getenv("RTE2E_SKIP_STAGE")) == (e.in_pred ? 1 : 2));
    const bool was = e.runner->overwrite();
    e.runner->set_overwrite(want);
    e.runner->set_overrides(*e.fs, k);
    e.runner->set_overwrite(was);
}

void minmax(Env& e, const char* tag)
{
    if (!trace()) return;
    for (int l = 0; l < e.nl; ++l) {
        Fields& F = e.reg.fields(l);
        for (const char* nm : {"RHO", "RHOS", "TMP", "ZZ", "ZZS", "D", "DS", "MU"}) {
            if (!F.has(nm)) continue;
            const auto& mf = F[nm];
            for (int c = 0; c < mf.nComp(); ++c) amrex::Print() << "RTE2E-TRACE " << tag << " L" << l << " " << nm << "[" << c << "] min " << mf.min(c, 0) << " max " << mf.max(c, 0) << "\n";
        }
    }
}

// One stage (predictor or corrector) of the transport-only step. `prime` = the part that only refreshes the lagged DIVERGENCE_PART_1 state (no density update).
void stage(Env& e, bool pred, bool with_density, bool first_pass = true)
{
    TimeLoop& lp = e.loop;
    e.in_pred = pred;
    lp.stage_state(pred, first_pass);
    if (with_density) {
        TR("stage_viscosity");
        up(e, [&](int l) { lp.stage_viscosity(l, pred); });
    TR("compute_stage_fluxes");
        up(e, [&](int l) { e.fs->compute_stage_fluxes(l, pred); });
    TR("set_overrides");
        set_ov(e, FluxKind::Adv);
    TR("apply_flux_divergence");
        down(e, [&](int l) { e.fs->apply_flux_divergence(l, pred); });
        minmax(e, "after-density");
        TR("stage_exchange");
        up(e, [&](int l) { lp.stage_exchange(l, pred ? 1 : 4, pred); });
        TR("stage_boundary");
        up(e, [&](int l) { lp.stage_boundary(l, pred ? 1 : 4); });
    } else {
        TR("stage_viscosity");
        up(e, [&](int l) { lp.stage_viscosity(l, pred); });
        TR("stage_exchange");
        up(e, [&](int l) { lp.stage_exchange(l, pred ? 1 : 4, pred); });
        TR("stage_boundary");
        up(e, [&](int l) { lp.stage_boundary(l, pred ? 1 : 4); });
    }
    TR("stage_velocity_flux");
    up(e, [&](int l) { lp.stage_velocity_flux(l, pred); });
    TR("init_div");
    lp.stage_init_divergence();
    TR("stage_wall_bc");
    up(e, [&](int l) { lp.stage_wall_bc(l, pred); });
    TR("readout_dif");
    up(e, [&](int l) { e.fs->readout_dif(l); });
    TR("set_overrides");
    set_ov(e, FluxKind::Dif);
    TR("run_divergence_part1");
    down(e, [&](int l) { if (l > 0 && std::getenv("RTE2E_NO_REDIV1_FINE")) return; e.fs->run_divergence_part1(l); });
}

// velocity positions of the step (exchange + MATCH_VELOCITY / VELOCITY_BC), as TimeLoop::advance() has them; the velocity itself is prescribed and not updated
void velocity_match(Env& e, int code, bool pred)
{
    TR("velocity_match " << code);
    // Level 0 only: the velocity boundary routines of a fine level (MATCH_VELOCITY bookkeeping, fds_p_save_uvw) are not fine-ready in patches 0007-0009 (found by this test: the call
    // aborts with "mesh number 2 is not a level-0 FDS mesh"). The fine velocity is constant everywhere including its ghost faces here, so nothing is lost for this prescribed-velocity test.
    e.loop.stage_exchange(0, code, pred);
    e.loop.stage_boundary(0, code);
}

// A new level (and level 0 at the start) has no DEL_RHO_D_DEL_Z / D / KRES / MU of a previous stage: one corrector-form pass without the density update (notes/level-binding.md).
void prime(Env& e)
{
    e.loop.set_state(e.t, e.dt, 0);
    velocity_match(e, 6, false);
    stage(e, false, false, false);
    minmax(e, "after-prime");
}

// ---- diagnostics
struct Totals { std::vector<double> mass; double rho = 0.0, zmin = 1e300, zmax = -1e300, drho = 0.0, dz = 0.0, rmin = 1e300, rmax = -1e300, zsum = 0.0; };

Totals totals(Env& e, const std::vector<double>& rho_ref = {}, const char* rn = "RHO", const char* zn = "ZZ")
{
    Totals T;
    T.mass.assign(e.ns, 0.0);
    double rho_tot = 0.0, zmin = 1e300, zmax = -1e300, drho = 0.0, dz = 0.0, rmin = 1e300, rmax = -1e300, zsum = 0.0;
    for (int l = 0; l < e.nl; ++l) {
        const Level& L = e.reg.level(l);
        Fields& F = e.reg.fields(l);
        const amrex::iMultiFab* cov = e.reg.covered_mask(l);
        const double vol = L.dx[0] * L.dx[1] * L.dx[2];
        for (amrex::MFIter mfi(F[rn]); mfi.isValid(); ++mfi) {
            auto r = F[rn].const_array(mfi); auto z = F[zn].const_array(mfi);
            amrex::Array4<const int> cv; if (cov) cv = cov->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (cov && cv(i, j, k) != 0) return;
                rho_tot += r(i, j, k) * vol;
                rmin = std::min(rmin, r(i, j, k)); rmax = std::max(rmax, r(i, j, k));
                { double sz = 0.0; for (int n = 0; n < e.ns; ++n) sz += z(i, j, k, n); zsum = std::max(zsum, std::abs(sz - 1.0)); }
                for (int n = 0; n < e.ns; ++n) T.mass[n] += r(i, j, k) * z(i, j, k, n) * vol;
                zmin = std::min(zmin, z(i, j, k, 1)); zmax = std::max(zmax, z(i, j, k, 1));
                if (!rho_ref.empty()) { drho = std::max(drho, std::abs(r(i, j, k) - rho_ref[0])); dz = std::max(dz, std::abs(z(i, j, k, 1) - rho_ref[1])); }
            });
        }
    }
    const auto cm = amrex::ParallelContext::CommunicatorSub();
    amrex::ParallelAllReduce::Sum(rho_tot, cm);
    amrex::ParallelAllReduce::Sum(T.mass.data(), e.ns, cm);
    amrex::ParallelAllReduce::Min(zmin, cm); amrex::ParallelAllReduce::Max(zmax, cm); amrex::ParallelAllReduce::Min(rmin, cm); amrex::ParallelAllReduce::Max(rmax, cm); amrex::ParallelAllReduce::Max(zsum, cm);
    amrex::ParallelAllReduce::Max(drho, cm); amrex::ParallelAllReduce::Max(dz, cm);
    T.rho = rho_tot; T.rmin = rmin; T.rmax = rmax; T.zsum = zsum; T.zmin = zmin; T.zmax = zmax; T.drho = drho; T.dz = dz;
    return T;
}

void step(Env& e)
{
    ++e.icyc;
    e.loop.set_state(e.t, e.dt, e.icyc);
    stage(e, true, true);
    minmax(e, "after-pred");
    if (trace()) { const Totals ts = totals(e, {}, "RHOS", "ZZS"); amrex::Print() << "RTE2E-TRACE predictor composite rho*Z (star) " << ts.mass[0] << " " << ts.mass[1] << "\n"; }
    velocity_match(e, 3, true);
    e.t += e.dt;
    e.loop.set_state(e.t, e.dt, e.icyc);
    stage(e, false, true);
    velocity_match(e, 6, false);
    minmax(e, "after-corr");
}

// Valid cells of level `lev` (RHO and ZZ) to <prefix>.L<lev>.r<rank>.bin: int32 nboxes, per box 6 int32 (lo, hi) then RHO, ZZ_1..ZZ_ns as doubles, i fastest.
void dump_level(Env& e, int lev, const std::string& prefix)
{
    Fields& F = e.reg.fields(lev);
    const std::string fn = prefix + ".L" + std::to_string(lev) + ".r" + std::to_string(amrex::ParallelDescriptor::MyProc()) + ".bin";
    std::FILE* f = std::fopen(fn.c_str(), "wb");
    if (!f) amrex::Abort("cannot write " + fn);
    int nb = 0;
    for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) ++nb;
    std::fwrite(&nb, 4, 1, f);
    for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        int h[6] = {b.smallEnd(0), b.smallEnd(1), b.smallEnd(2), b.bigEnd(0), b.bigEnd(1), b.bigEnd(2)};
        std::fwrite(h, 4, 6, f);
        auto r = F["RHO"].const_array(mfi); auto z = F["ZZ"].const_array(mfi);
        std::vector<double> buf;
        for (int c = -1; c < e.ns; ++c)
            for (int k = h[2]; k <= h[5]; ++k) for (int j = h[1]; j <= h[4]; ++j) for (int i = h[0]; i <= h[3]; ++i) buf.push_back(c < 0 ? r(i, j, k) : z(i, j, k, c));
        std::fwrite(buf.data(), 8, buf.size(), f);
    }
    std::fclose(f);
}

int mode_transport(TimeLoop& loop, const fdsamr::Level0& l0, const Opt& o)
{
    Env e(loop, l0, o);
    int fails = 0;
    auto check = [&](bool ok, const std::string& what) {
        amrex::Print() << "RTE2E " << o.label << (ok ? " CHECK-OK " : " CHECK-FAIL ") << what << "\n";
        if (!ok) ++fails;
    };
    Fields& F0 = e.reg.fields(0);
    e.ns = F0["ZZ"].nComp();
    if (e.ns < 2) amrex::Abort("rt-e2e transport: the case needs two species (ZZ has " + std::to_string(e.ns) + " component)");
    e.tr.reset(new RegistryTransfer(e.reg, l0.dom.n_tracked));
    set_velocity(e, 0, o.u0, o.w0);
    init_profile(e, 0);
    e.nl = 1;
    if (o.mode != "none") { make_level1(e); e.nl = 2; set_velocity(e, 1, o.u0, o.w0); if (o.smooth > 0.0) { init_profile(e, 1); e.tr->hierarchy_done(false); } }
    for (int l = 0; l < e.nl; ++l) { e.reg.fields(l)["ZZ"].FillBoundary(e.reg.level(l).geom.periodicity()); }
    if (e.nl > 1) {
        loop.bind_level(1);
        install_cf_ghost_hooks(loop, e.nl - 1, ThermoProvider(), &e.cf_stats);
    }
    e.fs.reset(new fdsamr::FluxStages(loop));
    {
        std::vector<StageLevel> sl;
        for (int l = 0; l < e.nl; ++l) {
            const Level& L = e.reg.level(l);
            StageLevel s; s.geom = L.geom; s.ba = L.ba; s.dm = L.dm; s.ratio_from_parent = L.ref_ratio_from_parent;
            sl.push_back(s);
        }
        e.runner.reset(new FluxStageRunner(sl, o.overwrite));
    }
    const Level& Lt = e.reg.level(e.nl - 1);
    const double umax = std::abs(o.u0) / Lt.dx[0] + std::abs(o.w0) / Lt.dx[2];
    e.dt = o.dt > 0.0 ? o.dt : 0.24 / umax;
    e.t = 0.0;
    amrex::Print() << "RTE2E " << o.label << " levels = " << e.nl << ", level-1 boxes = " << (e.nl > 1 ? static_cast<int>(e.reg.level(1).ba.size()) : 0)
                   << ", dt = " << e.dt << ", steps = " << o.steps << ", overwrite = " << (o.overwrite ? "on" : "off") << ", velocity = (" << o.u0 << ", " << o.w0 << ")\n";
    if (e.nl > 1) {
        const GhostShareCount gc = count_shared_ghost_cells(e.reg.level(1).ba, e.reg.level(1).ref_ratio_from_parent, e.reg.level(0).geom);
        amrex::Print() << "RTE2E " << o.label << " D-059 shared ghost cells: " << gc.to_string() << "\n";
    }
    prime(e);
    if (std::getenv("RTE2E_PRIME2")) prime(e);
    const Totals t0 = totals(e);
    const double r0 = t0.rho / (e.reg.level(0).geom.ProbLength(0) * e.reg.level(0).geom.ProbLength(1) * e.reg.level(0).geom.ProbLength(2));
    std::vector<double> ref{0.0, 0.0};
    {   // reference values for the uniform check: the state at the start, first cell of the first level-0 box
        double rr = 0.0, zr = 0.0;
        for (amrex::MFIter mfi(F0["RHO"]); mfi.isValid(); ++mfi) { rr = F0["RHO"].const_array(mfi)(mfi.validbox().smallEnd(0), mfi.validbox().smallEnd(1), mfi.validbox().smallEnd(2)); zr = F0["ZZ"].const_array(mfi)(mfi.validbox().smallEnd(0), mfi.validbox().smallEnd(1), mfi.validbox().smallEnd(2), 1); break; }
        amrex::ParallelAllReduce::Max(rr, amrex::ParallelContext::CommunicatorSub());
        ref = {rr, zr};
    }
    (void)r0;
    for (int s = 0; s < o.steps; ++s) step(e);
    const Totals t1 = totals(e, ref);
    for (int n = 0; n < e.ns; ++n) {
        const double d = (t1.mass[n] - t0.mass[n]) / std::max(1e-300, std::abs(t0.mass[n]));
        amrex::Print() << "RTE2E " << o.label << " composite rho*Z" << n + 1 << " start " << t0.mass[n] << " end " << t1.mass[n] << " relative change " << d << "\n";
    }
    const double drho = (t1.rho - t0.rho) / t0.rho;
    amrex::Print() << "RTE2E " << o.label << " composite mass start " << t0.rho << " end " << t1.rho << " relative change " << drho << "\n";
    double worst = std::abs(drho);
    for (int n = 0; n < e.ns; ++n) worst = std::max(worst, std::abs((t1.mass[n] - t0.mass[n]) / std::max(1e-300, std::abs(t0.mass[n]))));
    amrex::Print() << "RTE2E " << o.label << " WORST-MASS-DRIFT " << worst << "\n";
    amrex::Print() << "RTE2E " << o.label << " RHO range [" << t1.rmin << ", " << t1.rmax << "], max |sum Z - 1| " << t1.zsum << "\n";
    amrex::Print() << "RTE2E " << o.label << " tracer range [" << t1.zmin << ", " << t1.zmax << "], deviation from the start values: max |rho - rho0| " << t1.drho << ", max |Z_tracer - Z0| " << t1.dz << "\n";
    amrex::Print() << "RTE2E " << o.label << " largest face flux: ADV " << e.max_flux[0] << ", DIF " << e.max_flux[1] << " (kg/m2/s of rho*Z, all components)\n";
    if (e.nl > 1) {
        amrex::Print() << "RTE2E " << o.label << " runner: ADV entries " << e.runner->stats().adv_entries << ", DIF entries " << e.runner->stats().dif_entries << ", set_overrides calls " << e.runner->stats().calls
                       << "; cf hook: calls " << e.cf_stats.calls << ", covered coarse cells written " << e.cf_stats.covered_cells << ", fine ghost cells written " << e.cf_stats.fine_ghost_cells
                       << ", conflicts " << e.cf_stats.conflicts << "\n";
    }
    if (o.no_blob) { check(t1.dz < 1.0e-12 && t1.drho < 1.0e-12, "uniform state stays uniform (max |dRHO| " + std::to_string(t1.drho) + ", max |dZ| " + std::to_string(t1.dz) + ")"); }
    if (!o.dump.empty()) {
        if (e.nl > 1) e.tr->hierarchy_done(false);
        for (int l = 0; l < e.nl; ++l) dump_level(e, l, o.dump);
    }
    if (e.nl > 1) loop.unbind_level(1);
    return fails;
}

// ---- E4: FR-016 ghost cells through a real stage. The reference dumps of the instrumented, unmodified FDS (ns2d_16_int_1to2_refinement: 12 coarse meshes of 4x1x4 cells + the 16x1x16 fine
// mesh 13 over coarse cells 4..11) are loaded into the registry fields of level 0 (this case's single 16x1x16 box; the 8x8 cells under the fine mesh are poison) and level 1 (the fine box),
// level 1 is bound, and TimeLoop::stage_exchange (same-level fill + Role 3's coarse-fine ghost hook, exactly what a step calls) is run. The ghost layers 1 and 2 that FDS itself wrote
// (coarse meshes inside the hole, fine mesh outside its box) are then compared with the registry values. Cells that two coarse meshes share (hole corners) are listed separately (D-059).
int mode_ghost(TimeLoop& loop, const fdsamr::Level0& l0, const Opt& o0)
{
    Opt o = o0;
    o.mode = "patch"; o.p[0] = 4; o.p[1] = 11; o.p[2] = 4; o.p[3] = 11; o.maxsize = 100000;
    Env e(loop, l0, o);
    int rc = 0;
    if (o.dumps_prefix.empty()) { amrex::Print() << "RTE2E ghost SKIP: no --dumps <prefix>\n"; return 0; }
    if (amrex::ParallelDescriptor::NProcs() != 1 || l0.ba.size() != 1 || l0.geom.Domain().length(0) != 16 || l0.geom.Domain().length(2) != 16) {
        amrex::Print() << "RTE2E ghost: needs 1 rank and the single 16x1x16 mesh of ns2d_16_l0 (got " << l0.ba.size() << " boxes, " << amrex::ParallelDescriptor::NProcs() << " ranks)\n";
        return 1;
    }
    e.ns = e.reg.fields(0)["ZZ"].nComp();
    e.tr.reset(new RegistryTransfer(e.reg, l0.dom.n_tracked));
    make_level1(e, false);
    e.nl = 2;
    loop.bind_level(1);
    const amrex::Box hole(amrex::IntVect(4, 0, 4), amrex::IntVect(11, 0, 11));
    const amrex::Box fbox = amrex::refine(hole, amrex::IntVect(2, 1, 2));
    const amrex::IntVect ratio(2, 1, 2);
    // mesh m (1..12): block (bi, bk) of 4x4 coarse cells, I fastest, the middle 2x2 blocks skipped
    std::vector<amrex::Box> mbox(14);
    {
        int m = 0;
        for (int bk = 0; bk < 4; ++bk) for (int bi = 0; bi < 4; ++bi) {
            if ((bi == 1 || bi == 2) && (bk == 1 || bk == 2)) continue;
            mbox[++m] = amrex::Box(amrex::IntVect(4 * bi, 0, 4 * bk), amrex::IntVect(4 * bi + 3, 0, 4 * bk + 3));
        }
        mbox[13] = fbox;
    }
    std::map<int, std::vector<fdsrt_test::DumpRec>> dumps;
    for (int m = 1; m <= 13; ++m) dumps[m] = fdsrt_test::read_dump(m == 1 ? o.dumps_prefix : o.dumps_prefix + "." + std::to_string(m));
    struct Case { const char* label; const char* rec; bool pred; std::vector<std::string> rho_zz; std::vector<std::string> extra; };
    const std::vector<Case> cases = {
        {"corrector RHO/ZZ/TMP/RSUM", "DENS_P", false, {"RHO", "ZZ"}, {"TMP", "RSUM"}},
        {"predictor RHOS/ZZS", "DENS_C", true, {"RHOS", "ZZS"}, {}},
        {"viscosity MU/KRES", "VISC_P", false, {}, {"MU", "KRES"}},
    };
    const std::vector<std::string> all_cell = {"RHO", "RHOS", "ZZ", "ZZS", "TMP", "RSUM", "MU", "KRES", "D", "DS"};
    for (const Case& cs : cases) {
        std::map<int, const fdsrt_test::DumpRec*> rec;
        for (int m = 1; m <= 13; ++m)
            for (const auto& r : dumps[m]) if (r.name == cs.rec && r.icyc == o.step) { rec[m] = &r; break; }
        if (rec.size() != 13) { amrex::Print() << "RTE2E ghost SKIP " << cs.label << ": record " << cs.rec << " step " << o.step << " not on all 13 meshes\n"; continue; }
        const bool have_rz = !cs.rho_zz.empty();
        Fields* Fl[2] = {&e.reg.fields(0), &e.reg.fields(1)};
        for (int l = 0; l < 2; ++l) for (const auto& n : all_cell) if (Fl[l]->has(n)) (*Fl[l])[n].setVal(1.0, 0, (*Fl[l])[n].nComp(), (*Fl[l])[n].nGrowVect());
        std::vector<std::string> loaded = cs.extra;
        if (have_rz) { loaded.push_back(cs.rho_zz[0]); loaded.push_back(cs.rho_zz[1]); }
        for (int l = 0; l < 2; ++l) for (const auto& n : loaded) (*Fl[l])[n].setVal(-1.0e30, 0, (*Fl[l])[n].nComp(), (*Fl[l])[n].nGrowVect());
        auto load = [&](const std::string& reg_name, const char* dump_name, int m) {
            const int lev = (m == 13) ? 1 : 0;
            auto it = rec[m]->bef.find(dump_name);
            if (it == rec[m]->bef.end()) return;
            const auto& arr = it->second;
            const amrex::Box vb = mbox[m];
            amrex::MultiFab& mf = (*Fl[lev])[reg_name];
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                auto a = mf.array(mfi);
                amrex::LoopOnCpu(vb, [&](int i, int j, int k) { a(i, j, k, 0) = arr.at(i - vb.smallEnd(0) + 1, j - vb.smallEnd(1) + 1, k - vb.smallEnd(2) + 1); });
            }
        };
        for (int m = 1; m <= 13; ++m) {
            if (have_rz) { load(cs.rho_zz[0], cs.rho_zz[0].c_str(), m); load(cs.rho_zz[1], cs.rho_zz[1].c_str(), m); }
            for (const auto& n : cs.extra) load(n, n.c_str(), m);
        }
        const double rsum0 = rec[13]->bef.count("RSUM") ? rec[13]->bef.at("RSUM").at(1, 1, 1) : 0.0;
        const double pbar0 = rec[13]->bef.count("PBAR") ? rec[13]->bef.at("PBAR").at(1, 1, 1) : 0.0;
        ThermoProvider th;
        th.n_tracked = 1; th.rsum = [rsum0](const double*) { return rsum0; }; th.pbar = [pbar0](int, const amrex::IntVect&) { return pbar0; };
        const bool use_th = have_rz && rsum0 > 0 && pbar0 > 0 && !cs.pred;
        CfHookStats st;
        install_cf_ghost_hooks(loop, 1, use_th ? th : ThermoProvider(), &st);
        const int code = cs.pred ? 1 : 4;
        loop.set_state(0.0, 1.0e-3, o.step);
        for (int l = 0; l < 2; ++l) loop.stage_exchange(l, code, cs.pred);
        amrex::Print() << "RTE2E ghost " << cs.label << " (record " << cs.rec << ", step " << o.step << ", stage_exchange code " << code << "): covered cells written " << st.covered_cells << ", fine ghost cells written "
                       << st.fine_ghost_cells << ", conflict cells " << st.conflicts << (use_th ? "" : " (no EOS: TMP/RSUM not rebuilt)") << "\n";
        struct Stat { long n = 0, bad = 0; double maxrel = 0; };
        std::map<std::string, Stat> stat;
        auto in_corner = [&](const amrex::IntVect& c) {
            int nr = 0;
            for (int d = 0; d < 3; ++d) { if (hole.length(d) <= 1) continue; if (c[d] - hole.smallEnd(d) <= 1 || hole.bigEnd(d) - c[d] <= 1) ++nr; }
            return nr >= 2;
        };
        for (int m = 1; m <= 13; ++m) {
            const int lev = (m == 13) ? 1 : 0;
            const amrex::Box vb = mbox[m];
            const auto& r = *rec[m];
            const amrex::Box shell = amrex::grow(vb, amrex::IntVect(2, 0, 2));
            amrex::LoopOnCpu(shell, [&](int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                if (vb.contains(iv)) return;
                int nout = 0, layer = 0;
                for (int d = 0; d < 3; ++d) { const int dd = std::max({vb.smallEnd(d) - iv[d], iv[d] - vb.bigEnd(d), 0}); if (dd > 0) { ++nout; layer = dd; } }
                if (nout != 1) return;
                if (lev == 0 && !hole.contains(iv)) return;
                if (lev == 1 && !e.reg.level(1).geom.Domain().contains(iv)) return;
                const bool cz = (lev == 0) && in_corner(iv);
                const std::string side = lev ? "fine" : "coarse";
                auto one = [&](const std::string& rname, const char* dname, int comp, int maxlayer) {
                    const amrex::MultiFab& mf = (*Fl[lev])[rname];
                    if (layer > maxlayer || layer > mf.nGrow() || !r.bef.count(dname)) return;
                    const auto& arr = r.bef.at(dname);
                    const int fi = i - vb.smallEnd(0) + 1, fj = j - vb.smallEnd(1) + 1, fk = k - vb.smallEnd(2) + 1;
                    if (!arr.has(fi, fj, fk)) return;
                    double ours = std::nan("");
                    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) if (mf.fabbox(mfi.index()).contains(iv)) ours = mf.const_array(mfi)(i, j, k, comp);
                    const double fds = arr.at(fi, fj, fk);
                    const double rel = std::abs(ours - fds) / std::max(1e-300, std::abs(fds));
                    Stat& s = stat[side + " " + dname + " layer " + std::to_string(layer) + (cz ? " [hole edge/corner]" : "")];
                    ++s.n; s.maxrel = std::max(s.maxrel, rel);
                    if (!(rel <= 1e-12)) ++s.bad;
                };
                if (have_rz) { one(cs.rho_zz[0], cs.rho_zz[0].c_str(), 0, 2); one(cs.rho_zz[1], cs.rho_zz[1].c_str(), 0, 2); }
                for (const auto& n : cs.extra) one(n, n.c_str(), 0, (n == "TMP") ? 2 : 1);
            });
        }
        for (auto& kv : stat) {
            amrex::Print() << "RTE2E ghost   " << kv.first << ": cells " << kv.second.n << ", mismatches " << kv.second.bad << ", max rel diff " << kv.second.maxrel << "\n";
            if (kv.first.find("[hole edge/corner]") == std::string::npos && kv.second.bad > 0) ++rc;
        }
    }
    loop.unbind_level(1);
    return rc;
}

}  // namespace

int driver_mode(int argc, char** argv, TimeLoop& loop, const fdsamr::Level0& l0, double /*dt_setup*/)
{
    const Opt o = parse(argc, argv);
    if (o.test == "list") { amrex::Print() << "RTE2E modes: transport ghost\n"; return 0; }
    int fails = 0;
    if (o.test == "transport") fails = mode_transport(loop, l0, o);
    else if (o.test == "ghost") fails = mode_ghost(loop, l0, o);
    else amrex::Abort("rt-e2e: unknown test " + o.test);
    amrex::Print() << "RTE2E " << o.label << (fails == 0 ? " PASS" : " FAIL") << "\n";
    return fails;
}

}  // namespace fdsrt
