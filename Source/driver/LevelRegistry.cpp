// LevelRegistry.cpp: see LevelRegistry.H. Kernel-facing rules (M2a): (a) passive scalars are the Fields' business; (b) only uniform Cartesian metrics are used.
// No Fortran symbol is referenced, so the unit tests link this file without the FDS objects.
#include "LevelRegistry.H"

#include <AMReX_ParallelDescriptor.H>

namespace fdsamr {

CellWallProvider layout_cell_walls(const Level& lev)
{
    const amrex::Box dom = lev.geom.Domain();
    return [dom](int /*nm*/, const amrex::Box& vb, amrex::Array4<int> const& a) {
        amrex::LoopOnCpu(vb, [&](int i, int j, int k) {
            a(i, j, k, 0) = 0;
            for (int f = 0; f < 6; ++f) {
                amrex::IntVect nb(i, j, k);
                nb[f / 2] += (f % 2 == 0) ? -1 : 1;
                // inside the box: no wall; outside the domain: wall (periodic edges too); else the neighbour is another box of this level or the coarser level: open
                a(i, j, k, 1 + f) = vb.contains(nb) ? kNoWall : !dom.contains(nb) ? kWall : kInterface;
            }
        });
    };
}

amrex::iMultiFab build_covered_mask(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::BoxArray& fine_ba, const amrex::IntVect& ratio)
{
    amrex::iMultiFab m(ba, dm, 1, 0);
    m.setVal(0);
    const amrex::BoxArray cba = amrex::coarsen(fine_ba, ratio);
    for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
        auto a = m.array(mfi);
        const amrex::Box vb = mfi.validbox();
        for (const auto& is : cba.intersections(vb)) {
            amrex::LoopOnCpu(is.second, [&](int i, int j, int k) { a(i, j, k) = 1; });
        }
    }
    return m;
}

LevelRegistry::LevelRegistry(const DomainInfo& dom, int nscalars, const std::vector<std::string>& names) : m_dom(dom), m_ns(nscalars), m_names(names) {}

void LevelRegistry::adopt_level0(const Level& l0, Fields& F, SideData& sd)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_have0, "LevelRegistry: level 0 adopted twice");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l0.level == 0, "LevelRegistry: adopt_level0 needs a level-0 Level");
    ensure(0);
    m_lv[0].lev = &l0; m_lv[0].F = &F; m_lv[0].sd = &sd;
    m_have0 = true;
    refresh_covered();
}

const Level& LevelRegistry::level(int l) const { AMREX_ALWAYS_ASSERT_WITH_MESSAGE(has_level(l), "LevelRegistry: no such level"); return *m_lv[l].lev; }
Fields& LevelRegistry::fields(int l) { AMREX_ALWAYS_ASSERT_WITH_MESSAGE(has_level(l), "LevelRegistry: no such level"); return *m_lv[l].F; }
const Fields& LevelRegistry::fields(int l) const { AMREX_ALWAYS_ASSERT_WITH_MESSAGE(has_level(l), "LevelRegistry: no such level"); return *m_lv[l].F; }
SideData& LevelRegistry::side_data(int l) { AMREX_ALWAYS_ASSERT_WITH_MESSAGE(has_level(l), "LevelRegistry: no such level"); return *m_lv[l].sd; }

void LevelRegistry::build_slot(const fdsrt::LevelLayout& l, Slot& s)
{
    // destroy in dependency order (SideData and Fields hold a reference to the Level)
    s.sd_own.reset(); s.F_own.reset(); s.lev_own.reset();
    s.lev_own.reset(new Level(make_layout_level(l.level, l.geom, l.ba, l.dm, l.ref_ratio_from_parent, m_dom)));
    s.lev = s.lev_own.get();
    s.F_own.reset(new Fields(*s.lev, m_ns, m_names));
    s.F = s.F_own.get();
    s.sd_own.reset();
    s.sd = nullptr;
    // the SideData (collective: FillBoundary) is built by rebuild_side_data
}

void LevelRegistry::make_level(const fdsrt::LevelLayout& l)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_have0, "LevelRegistry: make_level before adopt_level0");
    ++n_make;
    if (l.level == 0) {
        const Level& l0 = *m_lv[0].lev;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l.ba == l0.ba && l.dm.ProcessorMap() == l0.dm.ProcessorMap(),
                                         "LevelRegistry: level 0 is the FDS &MESH set; make_level(0) with a different layout is refused");
        return;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l.level <= num_levels(), "LevelRegistry: make_level skips a level (levels are made coarse to fine)");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!has_level(l.level), "LevelRegistry: make_level on an existing level (use remake_level)");
    ensure(l.level);
    build_slot(l, m_lv[l.level]);
    rebuild_side_data(l.level);
    refresh_covered();
}

void LevelRegistry::remake_level(const fdsrt::LevelLayout& l)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_have0, "LevelRegistry: remake_level before adopt_level0");
    if (l.level == 0) { make_level(l); return; }   // accepted only when unchanged (checked there)
    ++n_remake;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(has_level(l.level), "LevelRegistry: remake_level on a level that does not exist");
    Slot& s = m_lv[l.level];
    // keep the old objects readable for the transfer
    s.old_sd = std::move(s.sd_own); s.old_F = std::move(s.F_own); s.old_lev = std::move(s.lev_own);
    s.lev = nullptr; s.F = nullptr; s.sd = nullptr;
    // build_slot resets (now empty) own pointers and creates the new objects
    build_slot(l, s);
    rebuild_side_data(l.level);
    refresh_covered();
}

void LevelRegistry::clear_level(int level)
{
    if (level == 0) amrex::Abort("LevelRegistry: level 0 cannot be cleared (FDS &MESH set)");
    if (level >= num_levels()) return;
    ++n_clear;
    // a cleared level takes all finer levels with it (AmrCore clears from the finest down; be defensive)
    for (int l = num_levels() - 1; l >= level; --l) {
        Slot& s = m_lv[l];
        s.old_sd.reset(); s.old_F.reset(); s.old_lev.reset(); s.covered.reset();
        s.sd_own.reset(); s.F_own.reset(); s.lev_own.reset();
        s.lev = nullptr; s.F = nullptr; s.sd = nullptr;
    }
    m_lv.resize(level);
    refresh_covered();
}

void LevelRegistry::rebuild_side_data(int l)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(has_level(l), "LevelRegistry: rebuild_side_data on a level that does not exist");
    if (l == 0) amrex::Abort("LevelRegistry: the level-0 SideData is built from the FDS set-up by TimeLoop and is not rebuilt here");
    Slot& s = m_lv[l];
    s.sd_own.reset();
    s.sd_own.reset(new SideData(*s.lev, s.provider ? s.provider : layout_cell_walls(*s.lev)));
    s.sd = s.sd_own.get();
    ++n_side_rebuilds;
}

void LevelRegistry::refresh_covered()
{
    for (int l = 0; l < num_levels(); ++l) {
        m_lv[l].covered.reset();
        if (!has_level(l) || !has_level(l + 1)) continue;
        const Level& c = *m_lv[l].lev; const Level& f = *m_lv[l + 1].lev;
        m_lv[l].covered.reset(new amrex::iMultiFab(build_covered_mask(c.ba, c.dm, f.ba, f.ref_ratio_from_parent)));
    }
}

const amrex::iMultiFab* LevelRegistry::covered_mask(int l) const { return (l >= 0 && l < num_levels()) ? m_lv[l].covered.get() : nullptr; }

const Level* LevelRegistry::retired_level(int l) const { return (l >= 0 && l < num_levels()) ? m_lv[l].old_lev.get() : nullptr; }
Fields* LevelRegistry::retired_fields(int l) { return (l >= 0 && l < num_levels()) ? m_lv[l].old_F.get() : nullptr; }
void LevelRegistry::release_retired(int l)
{
    if (l < 0 || l >= num_levels()) return;
    m_lv[l].old_sd.reset(); m_lv[l].old_F.reset(); m_lv[l].old_lev.reset();
}

}  // namespace fdsamr
