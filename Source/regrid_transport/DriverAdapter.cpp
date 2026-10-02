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
        for (const std::string& n : nf) average_down_cells(Ff[n], Fc[n], lf.geom, lc.geom, lf.ref_ratio_from_parent, 0, Fc[n].nComp());
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
        const fdsamr::Level& lc = reg.level(l);
        const fdsamr::Level& lf = reg.level(l + 1);
        fdsamr::Fields& Fc = reg.fields(l);
        const fdsamr::Fields& Ff = reg.fields(l + 1);
        const std::vector<std::string> use = names.empty() ? Fc.names() : names;
        for (const std::string& n : use) {
            if (!is_transfer_scalar(Fc, n) || !Ff.has(n)) continue;
            average_down_cells(Ff[n], Fc[n], lf.geom, lc.geom, lf.ref_ratio_from_parent, 0, Fc[n].nComp());
        }
    }
}

}  // namespace fdsrt
