// RegridAmrCore.cpp: see RegridAmrCore.H.
#include "RegridAmrCore.H"

#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include "TagOps.H"

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
    if (transfer_ && lev > 0) transfer_->fill_initial_level(layout(lev));
}
void RegridAmrCore::MakeNewLevelFromCoarse(int lev, amrex::Real, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    if (listener_) listener_->make_level(layout(lev));
    if (transfer_) transfer_->fill_new_level(layout(lev));
}
void RegridAmrCore::RemakeLevel(int lev, amrex::Real, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    if (listener_) listener_->remake_level(layout(lev));
    if (transfer_) transfer_->fill_remade_level(layout(lev));
}
void RegridAmrCore::ClearLevel(int lev)
{
    if (listener_) listener_->clear_level(lev);
    if (transfer_) transfer_->level_cleared(lev);
}

amrex::DistributionMapping RegridAmrCore::MakeDistributionMap(int lev, amrex::BoxArray const& ba)
{
    if (lev == 0 && have_level0_dm_) return level0_dm_;
    if (dm_fn_) return dm_fn_(lev, ba);
    return amrex::AmrCore::MakeDistributionMap(lev, ba);
}

void RegridAmrCore::ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int)
{
    if (stats_.tags_set.size() < h_.levels.size()) {
        stats_.tags_set.assign(h_.levels.size(), 0);
        stats_.discarded_outside_region.assign(h_.levels.size(), 0);
        stats_.buffer_blocks_discarded.assign(h_.levels.size(), 0);
        reported_discard_.assign(h_.levels.size(), false);
    }
    if (tag_fn_) tag_fn_(lev, tags, time);
    // static finer &MESH footprints are always refined (IR-002): force-tag their coarsened boxes on this level
    if (lev + 1 <= h_.top) {
        const amrex::IntVect r = amrex::AmrCore::refRatio(lev);
        std::vector<amrex::Box> forced;
        for (const GridBox& g : h_.levels[lev + 1].grids)
            if (g.mesh >= 0) forced.push_back(amrex::coarsen(to_amrex(g.box), r));
        if (!forced.empty()) tag_boxes(tags, forced);
    }
    // IR-008: tags outside the refinable region are discarded and counted; one message per level
    std::vector<amrex::Box> allowed;
    for (const RegionBox& rb : h_.levels[lev].taggable) allowed.push_back(to_amrex(rb.box));
    const long removed = clip_tags_to_boxes(tags, allowed);
    stats_.discarded_outside_region[lev] += removed;
    if (removed > 0 && !reported_discard_[lev]) {
        reported_discard_[lev] = true;
        amrex::Print() << "AMR: " << removed << " tagged cell(s) on level " << lev << " lie outside the refinable region and were discarded (reported once per level)\n";
    }
    stats_.tags_set[lev] += count_tags(tags);
}

void RegridAmrCore::ManualTagsPlacement(int lev, amrex::TagBoxArray& tags, const amrex::Vector<amrex::IntVect>& bf_lev)
{
    // tags are in cells of level lev coarsened by bf_lev[lev]; the region boxes are aligned with that (snapped to the blocking factor of level+1)
    std::vector<amrex::Box> allowed;
    for (const RegionBox& rb : h_.levels[lev].taggable) allowed.push_back(amrex::coarsen(to_amrex(rb.box), bf_lev[lev]));
    const long removed = clip_tags_to_boxes(tags, allowed);
    if (stats_.buffer_blocks_discarded.size() > static_cast<size_t>(lev)) stats_.buffer_blocks_discarded[lev] += removed;
}

amrex::iMultiFab RegridAmrCore::covered_mask(int lev) const
{
    if (lev >= finestLevel()) {
        amrex::iMultiFab m(boxArray(lev), DistributionMap(lev), 1, 0);
        m.setVal(0);
        return m;
    }
    return make_covered_mask(boxArray(lev), DistributionMap(lev), boxArray(lev + 1), refRatio(lev));
}

void RegridAmrCore::init_from_tags(LevelListener& listener, amrex::Real time, bool notify_level0, const DmFn& dm_fn, const amrex::DistributionMapping* level0_dm)
{
    listener_ = &listener;
    dm_fn_ = dm_fn;
    // level 0 = the meshes as given
    amrex::BoxList bl;
    for (const GridBox& g : h_.levels[0].grids) bl.push_back(to_amrex(g.box));
    amrex::BoxArray ba0(bl);
    if (level0_dm) { level0_dm_ = *level0_dm; have_level0_dm_ = true; }
    amrex::DistributionMapping dm0 = MakeDistributionMap(0, ba0);
    SetBoxArray(0, ba0);
    SetDistributionMap(0, dm0);
    SetFinestLevel(0);
    if (notify_level0) listener_->make_level(layout(0));
    if (max_level > 0) {
        amrex::Vector<amrex::BoxArray> new_grids(max_level + 1);
        new_grids[0] = ba0;
        do {
            int new_finest;
            MakeNewGrids(finest_level, time, new_finest, new_grids);   // calls ErrorEst on the current finest level
            if (new_finest <= finest_level) break;
            finest_level = new_finest;
            amrex::DistributionMapping dm = MakeDistributionMap(new_finest, new_grids[new_finest]);
            MakeNewLevelFromScratch(new_finest, time, new_grids[new_finest], dm);   // make_level + fill_initial_level (new_finest is one more than before)
        } while (finest_level < max_level);
    }
    if (transfer_) transfer_->hierarchy_done(true);
}

bool RegridAmrCore::regrid_dynamic(amrex::Real time)
{
    if (max_level == 0) return false;
    ++stats_.n_regrids;
    const int old_finest = finestLevel();
    std::vector<amrex::BoxArray> old_ba;
    for (int l = 0; l <= old_finest; ++l) old_ba.push_back(boxArray(l));
    if (listener_) listener_->begin_regrid();
    regrid(0, time);
    if (listener_) listener_->end_regrid();
    bool changed = finestLevel() != old_finest;
    for (int l = 1; !changed && l <= old_finest; ++l) changed = !(boxArray(l) == old_ba[l]);
    if (changed) ++stats_.n_changed;
    if (transfer_) transfer_->hierarchy_done(false);
    return changed;
}

}  // namespace fdsrt
