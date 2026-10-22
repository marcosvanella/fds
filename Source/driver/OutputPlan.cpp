// OutputPlan.cpp: see OutputPlan.H. Keep the fill path short: COPY from the matching compute level, average_down one finer when present, OverrideSync on faces.
#include "../pressure_backend/FaceSync.H"
#include "OutputPlan.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>

#include "LevelOps.H"

namespace fdsamr {
namespace {

amrex::BoxArray chop_to_max(const amrex::BoxArray& ba, int maxsize)
{
    amrex::BoxList bl;
    for (int i = 0; i < ba.size(); ++i) {
        amrex::BoxArray tmp(ba[i]);
        tmp.maxSize(maxsize);
        for (int j = 0; j < tmp.size(); ++j) bl.push_back(tmp[j]);
    }
    return amrex::BoxArray(std::move(bl));
}

}  // namespace

amrex::DistributionMapping OutputPlan::writers_by_cells(const amrex::BoxArray& ba)
{
    const int nb = ba.size(), np = amrex::ParallelDescriptor::NProcs();
    amrex::Vector<int> pmap(nb, 0);
    if (nb == 0) return amrex::DistributionMapping(std::move(pmap));
    long total = 0;
    for (int i = 0; i < nb; ++i) total += ba[i].numPts();
    const long goal = std::max(1L, (total + np - 1) / np);
    int r = 0;
    long acc = 0;
    for (int i = 0; i < nb; ++i) {
        if (r < np - 1 && acc >= goal && i > 0) { ++r; acc = 0; }
        pmap[i] = r;
        acc += ba[i].numPts();
    }
    return amrex::DistributionMapping(std::move(pmap));
}

void OutputPlan::make(const LevelRegistry& reg, const OutputPlanOptions& o)
{
    m_o = o;
    m_g.clear();
    const Level& l0 = reg.level(0);
    const amrex::Box domain = l0.geom.Domain();
    amrex::IntVect r1(1);
    for (int k = 0; k < o.L_out; ++k) r1 *= o.ratio;

    if (o.refinable.ok()) {
        const amrex::Box fine = amrex::refine(o.refinable & domain, r1);
        OutputGroup g;
        g.level = o.L_out;
        amrex::RealBox rb = l0.geom.ProbDomain();
        g.geom.define(amrex::refine(domain, r1), rb, 0, l0.geom.isPeriodic());
        g.ba = chop_to_max(amrex::BoxArray(fine), o.maxsize);
        g.dm = writers_by_cells(g.ba);
        m_g.push_back(std::move(g));
    }

    {
        amrex::BoxList out;
        const amrex::Box ref = o.refinable.ok() ? (o.refinable & domain) : amrex::Box();
        for (int i = 0; i < l0.ba.size(); ++i) {
            const amrex::Box b = l0.ba[i];
            if (!ref.ok() || !b.intersects(ref)) { out.push_back(b); continue; }
            amrex::BoxList diff = amrex::boxDiff(b, ref);
            for (const auto& p : diff) if (p.ok()) out.push_back(p);
        }
        if (out.isNotEmpty()) {
            OutputGroup g;
            g.level = 0;
            g.geom = l0.geom;
            g.ba = chop_to_max(amrex::BoxArray(std::move(out)), o.maxsize);
            g.dm = writers_by_cells(g.ba);
            m_g.push_back(std::move(g));
        }
    }

    const int nc = std::max(1, static_cast<int>(o.cell_names.size()));
    for (auto& g : m_g) {
        g.cell = std::make_unique<amrex::MultiFab>(g.ba, g.dm, nc, 0);
        g.cell->setVal(0.0);
        g.amr_level = std::make_unique<amrex::MultiFab>(g.ba, g.dm, 1, 0);
        g.amr_level->setVal(static_cast<amrex::Real>(g.level));
        if (!o.face_names.empty()) {
            for (int d = 0; d < 3; ++d) {
                amrex::BoxArray fba = g.ba;
                fba.surroundingNodes(d);
                g.face[d] = std::make_unique<amrex::MultiFab>(fba, g.dm, 1, 0);
                g.face[d]->setVal(0.0);
            }
        }
    }
}

int OutputPlan::writer_of(int g, int m) const { return m_g.at(g).dm[m]; }

void OutputPlan::fill(const LevelRegistry& reg)
{
    const int nl = reg.num_levels();
    const int nc = static_cast<int>(m_o.cell_names.size());
    for (auto& g : m_g) {
        // Matching compute level: min(group level, finest). Coarser output group (level 0) gets average_down from level 1 when present.
        const int src = std::min(g.level, nl - 1);
        const Level& ls = reg.level(src);
        Fields& Fs = const_cast<Fields&>(reg.fields(src));

        amrex::MultiFab tmp(ls.ba, ls.dm, std::max(nc, 1), 0);
        tmp.setVal(0.0);
        for (int c = 0; c < nc; ++c) {
            if (!Fs.has(m_o.cell_names[c])) continue;
            amrex::MultiFab work(Fs[m_o.cell_names[c]], amrex::make_alias, 0, 1);
            if (g.level == 0 && src == 0 && nl > 1 && reg.fields(1).has(m_o.cell_names[c])) {
                // Covered coarse cells: replace with average_down of the fine field (display fill; integrals use D6 / ExactSum).
                amrex::MultiFab dst(ls.ba, ls.dm, 1, 0);
                dst.ParallelCopy(work, 0, 0, 1, 0, 0);
                amrex::MultiFab fine(const_cast<Fields&>(reg.fields(1))[m_o.cell_names[c]], amrex::make_alias, 0, 1);
                fdsrt::average_down_cells(fine, dst, reg.level(1).geom, ls.geom, reg.level(1).ref_ratio_from_parent, 0, 1);
                tmp.ParallelCopy(dst, 0, c, 1, 0, 0);
            } else {
                tmp.ParallelCopy(work, 0, c, 1, 0, 0);
            }
        }
        // Piecewise-constant injection from level 0 where the fine compute does not cover a fine output cell.
        if (g.level > 0 && src == g.level) {
            g.cell->ParallelCopy(tmp, 0, 0, nc, 0, 0);   // fine compute -> fine output (COPY; uncovered fine output cells stay 0 until injection)
            // Inject level 0 into a coarse-aligned temp, then expand into unset fine cells: for the static case the fine BA covers the refinable region, so this is a no-op when BAs match.
            if (nl > 0 && reg.fields(0).has(m_o.cell_names[0])) {
                // Skeleton: COPY already filled every cell under the fine BA; remaining zeros would need pc injection — left for a later pass when output BA exceeds fine BA.
            }
        } else {
            g.cell->ParallelCopy(tmp, 0, 0, nc, 0, 0);
        }
        g.amr_level->setVal(static_cast<amrex::Real>(src));

        if (m_o.face_names.empty()) continue;
        for (int d = 0; d < 3; ++d) {
            if (!Fs.has(m_o.face_names[d])) continue;
            amrex::MultiFab& mf = Fs[m_o.face_names[d]];
            mf.OverrideSync(ls.geom.periodicity());
            amrex::MultiFab work(mf.boxArray(), mf.DistributionMap(), 1, 0);
            work.ParallelCopy(mf, 0, 0, 1, 0, 0);
            if (g.level == 0 && nl > 1 && reg.fields(1).has(m_o.face_names[d])) {
                amrex::MultiFab& ff = const_cast<Fields&>(reg.fields(1))[m_o.face_names[d]];
                pb::sync_faces_low_owner(ff, reg.level(1).geom);   // low-face owner, as before every average-down
                amrex::average_down_faces(ff, work, reg.level(1).ref_ratio_from_parent, ls.geom);
            }
            g.face[d]->ParallelCopy(work, 0, 0, 1, 0, 0);
        }
    }
}

std::vector<double> composite_mass_exact(const LevelRegistry& reg, int nl, int ns, bool skip_fine_solid)
{
    std::vector<double> m(1 + ns, 0.0);
    std::vector<const amrex::MultiFab*> rhos;
    std::vector<double> vols;
    std::vector<amrex::iMultiFab> skip(nl);   // leaf mask: covered by a finer level, or (skip_fine_solid) a fine OBST cell that holds the burner wall state, not gas
    std::vector<const amrex::iMultiFab*> covs;
    for (int l = 0; l < nl; ++l) {
        rhos.push_back(&reg.fields(l)["RHO"]);
        const double* dx = reg.level(l).dx;
        vols.push_back(dx[0] * dx[1] * dx[2]);
        const amrex::iMultiFab* cv = reg.covered_mask(l);
        if (l > 0 && skip_fine_solid) {
            const amrex::iMultiFab& sd = const_cast<LevelRegistry&>(reg).side_data(l).mask();
            skip[l].define(sd.boxArray(), sd.DistributionMap(), 1, 0);
            for (amrex::MFIter mfi(skip[l]); mfi.isValid(); ++mfi) {
                auto o = skip[l].array(mfi);
                auto s = sd.const_array(mfi);
                amrex::Array4<const int> c;
                if (cv) c = cv->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { o(i, j, k) = (s(i, j, k, 0) != 0 || (cv && c(i, j, k) != 0)) ? 1 : 0; });
            }
            cv = &skip[l];
        }
        covs.push_back(cv);
    }
    m[0] = exact_sum_hierarchy(rhos, 0, vols, covs);
    for (int n = 0; n < ns; ++n) {
        double s = 0.0;
        for (int l = 0; l < nl; ++l)
            s += exact_sum_product_uncovered(reg.fields(l)["RHO"], 0, reg.fields(l)["ZZ"], n, vols[l], covs[l]);
        m[1 + n] = s;
    }
    return m;
}

}  // namespace fdsamr
