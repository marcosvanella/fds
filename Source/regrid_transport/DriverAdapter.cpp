// DriverAdapter.cpp: see DriverAdapter.H.
#include "DriverAdapter.H"

#include "LevelOps.H"

namespace fdsrt {

bool is_transfer_scalar(const fdsamr::Fields& F, const std::string& name)
{
    if (!F.has(name)) return false;
    if (name == "H" || name == "HS") return false;
    return F.spec(name).stag == fdsamr::Stag::Cell;
}

fdsamr::CfGhostHook make_cf_ghost_hook(fdsamr::LevelRegistry& reg, int nlayers)
{
    return [&reg, nlayers](const fdsamr::CfGhostRequest& rq) {
        const int l = rq.level;
        if (l < 1 || !reg.has_level(l) || !reg.has_level(l - 1)) return;
        fdsamr::Fields& Ff = reg.fields(l);
        const fdsamr::Fields& Fc = reg.fields(l - 1);
        const fdsamr::Level& lf = reg.level(l);
        const fdsamr::Level& lc = reg.level(l - 1);
        for (const std::string& n : rq.fields) {
            if (!is_transfer_scalar(Ff, n) || !Fc.has(n)) continue;
            amrex::MultiFab& mf = Ff[n];
            fill_cf_ghosts_pc(mf, Fc[n], lf.geom, lc.geom, lf.ref_ratio_from_parent, std::min(nlayers, mf.nGrow()), 0, mf.nComp());
        }
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
