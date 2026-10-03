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
        if (m_in_regrid && s.F_own) {   // D-058: keep the removed level readable until end_regrid
            Cleared& cl = m_cleared[l];
            cl.sd = std::move(s.sd_own); cl.F = std::move(s.F_own); cl.lev = std::move(s.lev_own);
        }
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

const Level* LevelRegistry::retired_level(int l) const
{
    if (l >= 0 && l < num_levels() && m_lv[l].old_lev) return m_lv[l].old_lev.get();
    auto it = m_cleared.find(l);
    return it != m_cleared.end() ? it->second.lev.get() : nullptr;
}
Fields* LevelRegistry::retired_fields(int l)
{
    if (l >= 0 && l < num_levels() && m_lv[l].old_F) return m_lv[l].old_F.get();
    auto it = m_cleared.find(l);
    return it != m_cleared.end() ? it->second.F.get() : nullptr;
}
void LevelRegistry::release_retired(int l)
{
    m_cleared.erase(l);
    if (l < 0 || l >= num_levels()) return;
    m_lv[l].old_sd.reset(); m_lv[l].old_F.reset(); m_lv[l].old_lev.reset();
}

void LevelRegistry::begin_regrid()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_in_regrid, "LevelRegistry: begin_regrid inside a regrid (the brackets do not nest)");
    m_in_regrid = true;
    ++n_begin_regrid;
}

void LevelRegistry::end_regrid()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_in_regrid, "LevelRegistry: end_regrid without begin_regrid");
    for (int l = 0; l < num_levels(); ++l) { m_lv[l].old_sd.reset(); m_lv[l].old_F.reset(); m_lv[l].old_lev.reset(); }
    m_cleared.clear();
    m_in_regrid = false;
    ++n_end_regrid;
}

void LevelRegistry::fill_initial_level(int l, const InitialFill& f)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l >= 1 && has_level(l), "LevelRegistry::fill_initial_level: level must exist and be > 0 (level 0 is the FDS initialisation)");
    Fields& F = fields(l);
    const Level& lev = *m_lv[l].lev;
    const int ns = m_ns;
    const double* plo = lev.geom.ProbLo();
    const double dx[3] = {lev.dx[0], lev.dx[1], lev.dx[2]};
    auto has = [&](const char* n) { return F.has(n); };
    if (f.cell) {
        for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            std::vector<double> zz(ns);
            for (int k = vb.smallEnd(2); k <= vb.bigEnd(2); ++k)
                for (int j = vb.smallEnd(1); j <= vb.bigEnd(1); ++j)
                    for (int i = vb.smallEnd(0); i <= vb.bigEnd(0); ++i) {
                        double rho = 0.0, tmp = 0.0;
                        f.cell(plo[0] + (i + 0.5) * dx[0], plo[1] + (j + 0.5) * dx[1], plo[2] + (k + 0.5) * dx[2], &rho, &tmp, zz.data());
                        auto put = [&](const char* n, double v) { if (has(n)) F[n].array(mfi)(i, j, k, 0) = v; };
                        put("RHO", rho); put("RHOS", rho); put("TMP", tmp);
                        for (int n = 0; n < ns; ++n) {
                            if (has("ZZ"))  F["ZZ"].array(mfi)(i, j, k, n) = zz[n];
                            if (has("ZZS")) F["ZZS"].array(mfi)(i, j, k, n) = zz[n];
                        }
                    }
        }
    }
    if (f.velocity) {
        const char* nm[3][2] = {{"U", "US"}, {"V", "VS"}, {"W", "WS"}};
        for (int d = 0; d < 3; ++d) {
            if (!has(nm[d][0])) continue;
            for (amrex::MFIter mfi(F[nm[d][0]]); mfi.isValid(); ++mfi) {
                amrex::Box fb = mfi.validbox();   // nodal in d for a face field
                auto a = F[nm[d][0]].array(mfi);
                for (int k = fb.smallEnd(2); k <= fb.bigEnd(2); ++k)
                    for (int j = fb.smallEnd(1); j <= fb.bigEnd(1); ++j)
                        for (int i = fb.smallEnd(0); i <= fb.bigEnd(0); ++i) {
                            const int idx[3] = {i, j, k};
                            double x[3];
                            for (int e = 0; e < 3; ++e) x[e] = plo[e] + (idx[e] + (e == d ? 0.0 : 0.5)) * dx[e];
                            a(i, j, k, 0) = f.velocity(d, x[0], x[1], x[2]);
                        }
                if (has(nm[d][1])) {
                    auto as = F[nm[d][1]].array(mfi);
                    for (int k = fb.smallEnd(2); k <= fb.bigEnd(2); ++k)
                        for (int j = fb.smallEnd(1); j <= fb.bigEnd(1); ++j)
                            for (int i = fb.smallEnd(0); i <= fb.bigEnd(0); ++i) as(i, j, k, 0) = a(i, j, k, 0);
                }
            }
        }
    }
}

}  // namespace fdsamr
