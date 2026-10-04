// DriverAdapter.cpp: see DriverAdapter.H.
#include "DriverAdapter.H"

#include <algorithm>

#include "LevelOps.H"

namespace fdsrt {

bool is_transfer_scalar(const fdsamr::Fields& F, const std::string& name)
{
    if (!F.has(name)) return false;
    if (name == "H" || name == "HS") return false;
    return F.spec(name).stag == fdsamr::Stag::Cell;
}

// Restriction of the named fields of level l to level l-1: a density with its mass fractions (RHO+ZZ, RHOS+ZZS) is restricted by average_down_species (rho and rho*Z
// averaged, Z rebuilt); a mass-fraction field named without its density gets the density of the same stage averaged as well. Other fields: plain volume average.
static void average_down_fields(fdsamr::LevelRegistry& reg, int l, const std::vector<std::string>& names)
{
    const fdsamr::Level& lf = reg.level(l);
    const fdsamr::Level& lc = reg.level(l - 1);
    fdsamr::Fields& Ff = reg.fields(l);
    fdsamr::Fields& Fc = reg.fields(l - 1);
    auto has = [&](const std::string& n) { return std::find(names.begin(), names.end(), n) != names.end(); };
    const amrex::iMultiFab* cov = reg.covered_mask(l - 1);
    for (const char* pr : {"RHO", "RHOS"}) {
        const std::string rho = pr, zz = (rho == "RHO") ? "ZZ" : "ZZS";
        const bool want_zz = has(zz) && Ff.has(zz) && Fc.has(zz) && Ff.has(rho) && Fc.has(rho);
        if (want_zz) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cov != nullptr, "average_down_fields: no covered mask");
            average_down_species(Ff[rho], Ff[zz], Fc[rho], Fc[zz], *cov, lf.geom, lc.geom, lf.ref_ratio_from_parent);
        } else if (has(rho) && Ff.has(rho) && Fc.has(rho)) {
            average_down_cells(Ff[rho], Fc[rho], lf.geom, lc.geom, lf.ref_ratio_from_parent, 0, 1);
        }
    }
    for (const std::string& n : names) {
        if (n == "RHO" || n == "RHOS" || n == "ZZ" || n == "ZZS") continue;
        if (!Ff.has(n) || !Fc.has(n)) continue;
        average_down_cells(Ff[n], Fc[n], lf.geom, lc.geom, lf.ref_ratio_from_parent, 0, Fc[n].nComp());
    }
}

fdsamr::CfGhostHook make_cf_ghost_hook(fdsamr::LevelRegistry& reg, const ThermoProvider& th, CfHookStats* stats)
{
    return [&reg, th, stats](const fdsamr::CfGhostRequest& rq) {
        const int l = rq.level;
        if (l < 1 || !reg.has_level(l) || !reg.has_level(l - 1)) return;
        fdsamr::Fields& Ff = reg.fields(l);
        fdsamr::Fields& Fc = reg.fields(l - 1);
        const fdsamr::Level& lf = reg.level(l);
        const fdsamr::Level& lc = reg.level(l - 1);
        auto in_req = [&](const char* n) { return std::find(rq.fields.begin(), rq.fields.end(), n) != rq.fields.end(); };
        auto pick = [&](fdsamr::Fields& F, const char* n) -> amrex::MultiFab* { return (in_req(n) && is_transfer_scalar(F, n)) ? &F[n] : nullptr; };
        auto stage = [&](fdsamr::Fields& F, std::vector<std::string>& names) {
            ScalarStage s;
            const bool pred = in_req("RHOS");
            s.rho = pick(F, pred ? "RHOS" : "RHO");
            s.zz = pick(F, pred ? "ZZS" : "ZZ");
            s.tmp = pick(F, "TMP");
            s.rsum = pick(F, "RSUM");
            if (s.rho) names.push_back(pred ? "RHOS" : "RHO");
            if (s.zz) names.push_back(pred ? "ZZS" : "ZZ");
            if (s.tmp) names.push_back("TMP");
            if (s.rsum) names.push_back("RSUM");
            for (const char* n : {"MU", "KRES", "D", "DS"})
                if (auto* m = pick(F, n)) { s.mean.push_back(m); names.push_back(n); }
            return s;
        };
        std::vector<std::string> nf, nc;
        ScalarStage fs = stage(Ff, nf), cs = stage(Fc, nc);
        if (!fs.rho || !cs.rho || nf != nc) return;
        average_down_fields(reg, l, nf);
        Thermo tc, tf;
        if (th.rsum && th.pbar) {
            tc.n_tracked = tf.n_tracked = th.n_tracked;
            tc.rsum = tf.rsum = th.rsum;
            tc.pbar = [&th, l](const amrex::IntVect& c) { return th.pbar(l - 1, c); };
            tf.pbar = [&th, l](const amrex::IntVect& c) { return th.pbar(l, c); };
        }
        const bool have_th = th.rsum && th.pbar;
        const CfStats sc = fill_covered_ghosts_fds(cs, fs, lf.geom, lc.geom, lf.ref_ratio_from_parent, have_th ? &tc : nullptr);
        for (const std::string& n : nc) Fc.fill_ghosts(n);   // the coarse boxes next to the covered region read these cells as their ghost cells
        const CfStats sf = fill_fine_ghosts_fds(fs, cs, lf.geom, lc.geom, lf.ref_ratio_from_parent, have_th ? &tf : nullptr);
        if (stats) { stats->fine_ghost_cells += sf.ghost_cells; stats->covered_cells += sc.ghost_cells; stats->conflicts += sc.conflicts; ++stats->calls; }
    };
}

void average_down_registry(fdsamr::LevelRegistry& reg, const std::vector<std::string>& names)
{
    for (int l = reg.num_levels() - 2; l >= 0; --l) {
        const fdsamr::Fields& Fc = reg.fields(l);
        const fdsamr::Fields& Ff = reg.fields(l + 1);
        const std::vector<std::string> use = names.empty() ? Fc.names() : names;
        std::vector<std::string> sel;
        for (const std::string& n : use)
            if (is_transfer_scalar(Fc, n) && Ff.has(n)) sel.push_back(n);
        average_down_fields(reg, l + 1, sel);
    }
}

std::vector<ProjectionLevel> registry_projection_levels(fdsamr::LevelRegistry& reg, const RegridAmrCore& core, const std::function<const amrex::MultiFab*(int level)>& D)
{
    std::vector<ProjectionLevel> lv;
    for (int l = 0; l <= core.finestLevel(); ++l) {
        fdsamr::Fields& F = reg.fields(l);
        ProjectionLevel pl;
        pl.geom = core.Geom(l);
        pl.ref_ratio = l > 0 ? core.refRatio(l - 1) : amrex::IntVect(1);
        pl.vel = {&F["U"], &F["V"], &F["W"]};
        pl.D = D ? D(l) : nullptr;
        pl.covered = reg.covered_mask(l);
        lv.push_back(pl);
    }
    return lv;
}

void install_post_regrid_projection(RegridAmrCore& core, fdsamr::LevelRegistry& reg, PostRegridProjectionState& state)
{
    RegridAmrCore* cp = &core;
    fdsamr::LevelRegistry* rp = &reg;
    PostRegridProjectionState* sp = &state;
    core.set_post_regrid_hook([cp, rp, sp](const RegridOutcome&) {
        std::vector<ProjectionLevel> lv = registry_projection_levels(*rp, *cp, sp->D);
        if (sp->before) sp->before(lv);
        ProjectionOptions o = sp->options;
        o.enabled = true;
        sp->last = project_after_regrid(lv, sp->solver, o);
        ++sp->calls;
        if (sp->last.ran && sp->last.solved) {
            const char* st[3] = {"US", "VS", "WS"};
            const char* vn[3] = {"U", "V", "W"};
            for (int l = 0; l < static_cast<int>(lv.size()); ++l) {
                fdsamr::Fields& F = rp->fields(l);
                for (int d = 0; d < 3; ++d)
                    if (F.has(st[d])) amrex::MultiFab::Copy(F[st[d]], F[vn[d]], 0, 0, 1, 0);
            }
        }
        if (sp->after) sp->after(sp->last, lv);
    });
}

}  // namespace fdsrt
