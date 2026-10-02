// RegridAmrCore.cpp: see RegridAmrCore.H.
#include "RegridAmrCore.H"

#include <AMReX_RealBox.H>

namespace fdsrt {

amrex::Box to_amrex(const IBox& b)
{
    return amrex::Box(amrex::IntVect(b.lo[0], b.lo[1], b.lo[2]), amrex::IntVect(b.hi[0], b.hi[1], b.hi[2]));
}

amrex::Geometry make_level0_geometry(const Hierarchy& h)
{
    amrex::RealBox rb({h.dom_lo[0], h.dom_lo[1], h.dom_lo[2]}, {h.dom_hi[0], h.dom_hi[1], h.dom_hi[2]});
    amrex::Array<int, 3> per{h.periodic[0] ? 1 : 0, h.periodic[1] ? 1 : 0, h.periodic[2] ? 1 : 0};
    return amrex::Geometry(to_amrex(h.levels[0].domain), rb, 0, per);
}

amrex::AmrInfo make_amr_info(const Hierarchy& h, const AmrParams& p)
{
    amrex::AmrInfo info;
    info.max_level = h.max_level;
    info.ref_ratio.clear();
    info.blocking_factor.clear();
    info.max_grid_size.clear();
    info.n_error_buf.clear();
    for (int l = 0; l <= h.max_level; ++l) {
        const Level& L = h.levels[l];
        info.blocking_factor.push_back(amrex::IntVect(L.blocking_factor[0], L.blocking_factor[1], L.blocking_factor[2]));
        info.max_grid_size.push_back(amrex::IntVect(L.max_grid_size[0], L.max_grid_size[1], L.max_grid_size[2]));
        info.n_error_buf.push_back(amrex::IntVect(p.n_error_buf));
        if (l < h.max_level) {
            const IVec& r = h.levels[l + 1].ref_from_parent;
            info.ref_ratio.push_back(amrex::IntVect(r[0], r[1], r[2]));
        }
    }
    if (info.ref_ratio.empty()) info.ref_ratio.push_back(amrex::IntVect(2));  // AmrInfo needs at least one entry
    info.grid_eff = p.grid_eff;
    info.n_proper = p.n_proper;
    info.check_input = false;  // level-0 grids are the meshes as given (not blocking-factor multiples); R1 already ran the checks
    return info;
}

RegridAmrCore::RegridAmrCore(const Hierarchy& h, const AmrParams& p)
    : amrex::AmrCore(make_level0_geometry(h), make_amr_info(h, p)), h_(h)
{
}

LevelLayout RegridAmrCore::layout(int level) const
{
    LevelLayout l;
    l.level = level;
    l.geom = Geom(level);
    l.ba = boxArray(level);
    l.dm = DistributionMap(level);
    if (level > 0) {
        const IVec& r = h_.levels[level].ref_from_parent;
        l.ref_ratio_from_parent = amrex::IntVect(r[0], r[1], r[2]);
    }
    return l;
}

void RegridAmrCore::init_static(LevelListener& listener, bool notify_level0, const DmFn& dm_fn, const amrex::DistributionMapping* level0_dm)
{
    listener_ = &listener;
    const int top = h_.top_level();
    for (int l = 0; l <= top; ++l) {
        amrex::BoxList bl;
        for (const GridBox& g : h_.levels[l].grids) bl.push_back(to_amrex(g.box));
        amrex::BoxArray ba(bl);
        amrex::DistributionMapping dm;
        if (l == 0 && level0_dm) dm = *level0_dm;
        else if (dm_fn) dm = dm_fn(l, ba);
        else dm = amrex::DistributionMapping(ba);
        SetBoxArray(l, ba);
        SetDistributionMap(l, dm);
    }
    SetFinestLevel(top);
    for (int l = notify_level0 ? 0 : 1; l <= top; ++l) listener_->make_level(layout(l));
}

void RegridAmrCore::MakeNewLevelFromScratch(int lev, amrex::Real, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    if (listener_) listener_->make_level(layout(lev));
}
void RegridAmrCore::MakeNewLevelFromCoarse(int lev, amrex::Real t, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    MakeNewLevelFromScratch(lev, t, ba, dm);  // R4 adds the interpolation of the data
}
void RegridAmrCore::RemakeLevel(int lev, amrex::Real, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    if (listener_) listener_->remake_level(layout(lev));
}
void RegridAmrCore::ClearLevel(int lev)
{
    if (listener_) listener_->clear_level(lev);
}
void RegridAmrCore::ErrorEst(int, amrex::TagBoxArray&, amrex::Real, int) {}

}  // namespace fdsrt
