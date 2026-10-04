// RegistryTransfer.cpp: see RegistryTransfer.H.
#include "RegistryTransfer.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <memory>
#include <cmath>
#include <string>
#include <vector>

#include "CellTransfer.H"
#include "FaceTransfer.H"
#include "LevelOps.H"

namespace fdsamr {

namespace {
const char* const kCell[] = {"RHO", "RHOS", "TMP", "ZZ", "ZZS"};
const char* const kFace[3][2] = {{"U", "US"}, {"V", "VS"}, {"W", "WS"}};

// rho*Z_n on `ba`/`dm`, ncomp = ns, `ng` ghost layers computed from the ghost layers of RHO and ZZ (ng <= 2).
void make_rz(const Fields& F, int ns, int ng, amrex::MultiFab& rz)
{
    const amrex::MultiFab& rho = F["RHO"];
    const amrex::MultiFab& zz = F["ZZ"];
    for (amrex::MFIter mfi(rz); mfi.isValid(); ++mfi) {
        const amrex::Box b = amrex::grow(mfi.validbox(), ng);
        auto r = rho.const_array(mfi); auto z = zz.const_array(mfi); auto o = rz.array(mfi);
        amrex::LoopOnCpu(b, [&](int i, int j, int k) { for (int n = 0; n < ns; ++n) o(i, j, k, n) = r(i, j, k) * z(i, j, k, n); });
    }
}

// Parent value of a coarse cell field injected into the fine cells (piecewise constant).
void inject(const amrex::MultiFab& crse, amrex::MultiFab& fine, const amrex::IntVect& ratio)
{
    amrex::BoxArray cba = fine.boxArray(); cba.coarsen(ratio);
    amrex::MultiFab tmp(cba, fine.DistributionMap(), fine.nComp(), 0);
    tmp.ParallelCopy(crse, 0, 0, fine.nComp(), 0, 0);
    for (amrex::MFIter mfi(fine); mfi.isValid(); ++mfi) {
        auto c = tmp.const_array(mfi); auto f = fine.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const int ic = amrex::coarsen(i, ratio[0]), jc = amrex::coarsen(j, ratio[1]), kc = amrex::coarsen(k, ratio[2]);
            for (int n = 0; n < fine.nComp(); ++n) f(i, j, k, n) = c(ic, jc, kc, n);
        });
    }
}
}  // namespace

void RegistryTransfer::fill_initial_level(const fdsrt::LevelLayout& l)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(initial.cell || initial.velocity, "RegistryTransfer::fill_initial_level: no initial-condition function set (member `initial`)");
    m_reg.fill_initial_level(l.level, initial);
    ++m_st.n_initial;
    if (derive) derive(l.level);
}

void RegistryTransfer::fill_new_level(const fdsrt::LevelLayout& l) { prolong(l.level, false); ++m_st.n_new; }
void RegistryTransfer::fill_remade_level(const fdsrt::LevelLayout& l) { prolong(l.level, true); ++m_st.n_remade; }

void RegistryTransfer::prolong(int level, bool remade)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(level >= 1 && m_reg.has_level(level) && m_reg.has_level(level - 1), "RegistryTransfer: the level and its parent must exist");
    const Level& lf = m_reg.level(level);
    const Level& lc = m_reg.level(level - 1);
    Fields& Ff = m_reg.fields(level);
    Fields& Fc = m_reg.fields(level - 1);
    const amrex::IntVect ratio = lf.ref_ratio_from_parent;
    const int ns = Ff["ZZ"].nComp();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ntracked >= 1 && m_ntracked <= ns, "RegistryTransfer: n_tracked out of range");

    Fields* Fo = remade ? m_reg.retired_fields(level) : nullptr;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!remade || Fo != nullptr, "RegistryTransfer::fill_remade_level: the previous level objects are gone (call it between begin_regrid and end_regrid)");

    // ---- conserved densities rho*Z_n
    amrex::MultiFab rzc(lc.ba, lc.dm, ns, 1);
    rzc.setVal(0.0);
    make_rz(Fc, ns, 1, rzc);
    rzc.FillBoundary(lc.geom.periodicity());
    amrex::MultiFab rzf(lf.ba, lf.dm, ns, 0);
    std::unique_ptr<amrex::MultiFab> rzo;
    if (Fo) {
        const Level* lo = m_reg.retired_level(level);
        rzo.reset(new amrex::MultiFab(lo->ba, lo->dm, ns, 0));
        make_rz(*Fo, ns, 0, *rzo);
    }
    fdsrt::ProlongOpts po; po.use_floor = true; po.floor = 0.0;
    const fdsrt::ProlongStats ps = fdsrt::prolong_conserved(rzf, rzc, lc.geom, ratio, po, rzo.get());
    m_st.parents += ps.parents; m_st.limited += ps.limited; m_st.clips += ps.clips; m_st.unfixable += ps.unfixable;

    // ---- rho = sum of the tracked rho*Z, ZZ = rho*Z / rho
    const int nt = m_ntracked;
    for (amrex::MFIter mfi(rzf); mfi.isValid(); ++mfi) {
        auto q = rzf.const_array(mfi);
        auto rho = Ff["RHO"].array(mfi); auto rhos = Ff["RHOS"].array(mfi);
        auto zz = Ff["ZZ"].array(mfi); auto zzs = Ff["ZZS"].array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            double s = 0.0;
            for (int n = 0; n < nt; ++n) s += q(i, j, k, n);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(s > 0.0, "RegistryTransfer: non-positive density after the prolongation");
            rho(i, j, k) = s; rhos(i, j, k) = s;
            for (int n = 0; n < ns; ++n) { const double z = q(i, j, k, n) / s; zz(i, j, k, n) = z; zzs(i, j, k, n) = z; }
        });
    }
    // ---- TMP: parent value (fallback; `derive` may recompute it from the equation of state)
    inject(Fc["TMP"], Ff["TMP"], ratio);

    // ---- face velocities
    {
        std::array<std::unique_ptr<amrex::MultiFab>, 3> ct;
        for (int d = 0; d < 3; ++d) {
            amrex::BoxArray cba = lf.ba; cba.coarsen(ratio); cba.surroundingNodes(d);
            ct[d].reset(new amrex::MultiFab(cba, lf.dm, 1, 0));
            ct[d]->setVal(0.0);
            ct[d]->ParallelCopy(Fc[kFace[d][0]], 0, 0, 1, 0, 0);
        }
        for (amrex::MFIter mfi(Ff["RHO"]); mfi.isValid(); ++mfi) {
            std::array<const amrex::FArrayBox*, 3> cf{&(*ct[0])[mfi], &(*ct[1])[mfi], &(*ct[2])[mfi]};
            std::array<amrex::FArrayBox*, 3> ff{&Ff["U"][mfi], &Ff["V"][mfi], &Ff["W"][mfi]};
            fdsrt::prolong_faces_normal_linear(cf, ff, mfi.validbox(), ratio);
        }
    }
    for (int d = 0; d < 3; ++d) amrex::MultiFab::Copy(Ff[kFace[d][1]], Ff[kFace[d][0]], 0, 0, 1, 0);

    // ---- remade level: bitwise copy of the previous fine data where it overlaps
    if (Fo) {
        long n = 0;
        for (const char* nm : kCell) { Ff[nm].ParallelCopy(Fo->operator[](nm), 0, 0, Ff[nm].nComp(), 0, 0); }
        for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) Ff[kFace[d][s]].ParallelCopy(Fo->operator[](kFace[d][s]), 0, 0, 1, 0, 0);
        const amrex::BoxArray& ob = m_reg.retired_level(level)->ba;
        for (amrex::MFIter mfi(Ff["RHO"]); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            for (const amrex::Box& o : ob.boxList()) { amrex::Box x = b & o; if (x.ok()) n += x.numPts(); }
        }
        amrex::ParallelAllReduce::Sum(n, amrex::ParallelContext::CommunicatorSub());
        m_st.copied_old_cells += n;
    }
    if (derive) derive(level);
}

void RegistryTransfer::hierarchy_done(bool initial)
{
    ++m_st.n_done;
    for (int l = m_reg.num_levels() - 1; l >= 1; --l) {
        if (!m_reg.has_level(l)) continue;
        const Level& lf = m_reg.level(l);
        const Level& lc = m_reg.level(l - 1);
        Fields& Ff = m_reg.fields(l);
        Fields& Fc = m_reg.fields(l - 1);
        // Conservative average-down of the species: rho and rho*Z_n are volume-averaged, the mass fractions of the covered coarse cells are rebuilt as (rho*Z_n)/rho. Averaging the
        // mass fractions themselves would not conserve the species mass.
        const int ns = Ff["ZZ"].nComp();
        amrex::MultiFab rzf(lf.ba, lf.dm, ns, 0), rzc(lc.ba, lc.dm, ns, 0);
        make_rz(Ff, ns, 0, rzf);
        make_rz(Fc, ns, 0, rzc);
        fdsrt::average_down_cells(rzf, rzc, lf.geom, lc.geom, lf.ref_ratio_from_parent, 0, ns);
        for (const char* nm : {"RHO", "TMP"}) fdsrt::average_down_cells(Ff[nm], Fc[nm], lf.geom, lc.geom, lf.ref_ratio_from_parent, 0, 1);
        const amrex::iMultiFab* cov = m_reg.covered_mask(l - 1);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cov != nullptr, "RegistryTransfer::hierarchy_done: no covered mask");
        for (amrex::MFIter mfi(rzc); mfi.isValid(); ++mfi) {
            auto q = rzc.const_array(mfi); auto m = cov->const_array(mfi);
            auto rho = Fc["RHO"].array(mfi); auto rhos = Fc["RHOS"].array(mfi); auto zz = Fc["ZZ"].array(mfi); auto zzs = Fc["ZZS"].array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (!m(i, j, k)) return;
                rhos(i, j, k) = rho(i, j, k);
                for (int n = 0; n < ns; ++n) { zz(i, j, k, n) = q(i, j, k, n) / rho(i, j, k); zzs(i, j, k, n) = zz(i, j, k, n); }
            });
        }
        for (int d = 0; d < 3; ++d)
            for (int s = 0; s < 2; ++s)
                amrex::average_down_faces(Ff[kFace[d][s]], Fc[kFace[d][s]], lf.ref_ratio_from_parent, lc.geom);
    }
    for (int l = 0; l < m_reg.num_levels(); ++l) {
        if (!m_reg.has_level(l)) continue;
        Fields& F = m_reg.fields(l);
        const amrex::Periodicity per = m_reg.level(l).geom.periodicity();
        for (const char* nm : kCell) F[nm].FillBoundary(per);
        for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) F[kFace[d][s]].FillBoundary(per);
    }
    if (after_hierarchy) after_hierarchy(initial);
}

}  // namespace fdsamr
