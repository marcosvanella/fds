// two_level_transport.H: prescribed-velocity two-level scalar transport on the LevelRegistry fields (S12; AMReX only, no FDS kernels). PRE-VALIDATION: it tests the wiring of
//   Role 3's FluxOverrideOps / FluxStageRunner (interface flux overwrite, finest level first), RegistryTransfer (prolongation, conservative average-down, regrid) and the registry data
// layout, NOT the FDS operators on level 1. Included by test_units.cpp.
//
// The transported quantities are rho*Z_n (2 species, rho = their sum) with a centred advective face flux (face value x prescribed face velocity) and a two-point diffusive flux, one SSPRK2
// step (predictor, corrector) with the same dt on both levels. The velocity is the exactly divergence-free discrete stream-function field of the registry transfer test and does not change.
// Stage order (RegridInterface.H, flux-stage-wiring.md): compute fluxes on all levels, runner sets the ADV and DIF overrides finest first, cell update finest first, average-down after the step.
// This is a simplified single-phase form of the FDS order (the FDS density / DIVERGENCE_PART_1 split and the lagged DEL_RHO_D_DEL_Z are not modelled).
#ifndef FDS_TEST_TWO_LEVEL_TRANSPORT_H_
#define FDS_TEST_TWO_LEVEL_TRANSPORT_H_

#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <cmath>

#include "FluxStageRunner.H"
#include "RegistryTransfer.H"

namespace fdstest {

class RegTransport : public fdsrt::FluxAccess {
public:
    static constexpr int NC = 2;
    RegTransport(fdsamr::LevelRegistry& reg, fdsamr::RegistryTransfer& tr, double dt, double diff) : m_reg(reg), m_tr(tr), m_dt(dt), m_D(diff) { rebuild(); }

    struct LV {
        amrex::MultiFab q[3];                 // 0 = start of step, 1 = predictor result, 2 = corrector result; rho*Z_n, 1 ghost cell
        amrex::MultiFab adv[3], dif[3];       // stage face fluxes, nodal in d
        std::vector<std::vector<fdsrt::FluxOverride>> ov[2];
    };
    std::vector<LV> lv;

    void rebuild()   // after any change of the hierarchy
    {
        const int nl = m_reg.num_levels();
        lv.clear(); lv.resize(nl);
        std::vector<fdsrt::StageLevel> sl;
        for (int l = 0; l < nl; ++l) {
            const fdsamr::Level& L = m_reg.level(l);
            for (int s = 0; s < 3; ++s) { lv[l].q[s].define(L.ba, L.dm, NC, 1); lv[l].q[s].setVal(0.0); }
            for (int d = 0; d < 3; ++d) {
                const amrex::BoxArray nba = amrex::convert(L.ba, amrex::IntVect::TheDimensionVector(d));
                lv[l].adv[d].define(nba, L.dm, NC, 0); lv[l].dif[d].define(nba, L.dm, NC, 0);
                lv[l].adv[d].setVal(0.0); lv[l].dif[d].setVal(0.0);
            }
            fdsrt::StageLevel s;
            s.geom = L.geom; s.ba = L.ba; s.dm = L.dm; s.ratio_from_parent = L.ref_ratio_from_parent;
            sl.push_back(s);
        }
        m_runner.reset(new fdsrt::FluxStageRunner(sl, m_overwrite));
    }
    void set_overwrite(bool on) { m_overwrite = on; m_runner->set_overwrite(on); }

    // registry (RHO, ZZ) -> q[0] (valid cells), and back (RHO = sum, ZZ = q/RHO, RHOS, ZZS), then the conservative average-down of the registry transfer
    void pull()
    {
        for (int l = 0; l < static_cast<int>(lv.size()); ++l) {
            auto& F = m_reg.fields(l);
            for (amrex::MFIter mfi(lv[l].q[0]); mfi.isValid(); ++mfi) {
                auto q = lv[l].q[0].array(mfi); auto r = F["RHO"].const_array(mfi); auto z = F["ZZ"].const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < NC; ++n) q(i, j, k, n) = r(i, j, k) * z(i, j, k, n); });
            }
        }
    }
    void push()
    {
        for (int l = 0; l < static_cast<int>(lv.size()); ++l) {
            auto& F = m_reg.fields(l);
            for (amrex::MFIter mfi(lv[l].q[2]); mfi.isValid(); ++mfi) {
                auto q = lv[l].q[2].const_array(mfi); auto r = F["RHO"].array(mfi); auto rs = F["RHOS"].array(mfi); auto z = F["ZZ"].array(mfi); auto zs = F["ZZS"].array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    double s = 0.0;
                    for (int n = 0; n < NC; ++n) s += q(i, j, k, n);
                    r(i, j, k) = s; rs(i, j, k) = s;
                    for (int n = 0; n < NC; ++n) { z(i, j, k, n) = q(i, j, k, n) / s; zs(i, j, k, n) = z(i, j, k, n); }
                });
            }
        }
        m_tr.hierarchy_done(false);
    }

    void step()
    {
        for (int stage = 0; stage < 2; ++stage) {
            const bool pred = (stage == 0);
            const int nl = static_cast<int>(lv.size());
            for (int l = 0; l < nl; ++l) compute_stage_fluxes(l, pred);
            m_runner->set_overrides(*this, fdsrt::FluxKind::Adv);
            m_runner->set_overrides(*this, fdsrt::FluxKind::Dif);
            for (int l = nl - 1; l >= 0; --l) apply_flux_divergence(l, pred);
        }
        for (auto& L : lv) amrex::MultiFab::Copy(L.q[0], L.q[2], 0, 0, NC, 0);   // for push(): q[2] is the result
        push();
        pull();
    }

    // composite totals over uncovered cells: total[n] = sum rho*Z_n * cell volume
    std::array<double, NC> totals() const
    {
        std::array<double, NC> t{};
        for (int l = 0; l < static_cast<int>(lv.size()); ++l) {
            const fdsamr::Level& L = m_reg.level(l);
            const amrex::iMultiFab* cov = m_reg.covered_mask(l);
            const double vol = L.dx[0] * L.dx[1] * L.dx[2];
            for (amrex::MFIter mfi(lv[l].q[0]); mfi.isValid(); ++mfi) {
                auto q = lv[l].q[0].const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    if (cov && cov->const_array(mfi)(i, j, k) != 0) return;
                    for (int n = 0; n < NC; ++n) t[n] += q(i, j, k, n) * vol;
                });
            }
        }
        for (int n = 0; n < NC; ++n) amrex::ParallelAllReduce::Sum(t[n], amrex::ParallelContext::CommunicatorSub());
        return t;
    }

    // ---- fdsrt::FluxAccess
    void compute_stage_fluxes(int l, bool predictor) override
    {
        amrex::MultiFab& qin = lv[l].q[predictor ? 0 : 1];
        fill_ghost(l, qin);
        const fdsamr::Level& L = m_reg.level(l);
        auto& F = m_reg.fields(l);
        const char* vn[3] = {"U", "V", "W"};
        for (int d = 0; d < 3; ++d) {
            const double h = L.dx[d];
            for (amrex::MFIter mfi(lv[l].adv[d]); mfi.isValid(); ++mfi) {
                auto a = lv[l].adv[d].array(mfi); auto df = lv[l].dif[d].array(mfi); auto q = qin.const_array(mfi); auto u = F[vn[d]].const_array(mfi);
                const amrex::Box fb = mfi.validbox();
                const int e[3] = {d == 0, d == 1, d == 2};
                amrex::LoopOnCpu(fb, [&](int i, int j, int k) {
                    for (int n = 0; n < NC; ++n) {
                        if (d == 1) { a(i, j, k, n) = 0.0; df(i, j, k, n) = 0.0; continue; }   // hidden direction: one cell, no flux
                        const double ql = q(i - e[0], j - e[1], k - e[2], n), qr = q(i, j, k, n);
                        a(i, j, k, n) = 0.5 * (ql + qr) * u(i, j, k);
                        df(i, j, k, n) = -m_D * (qr - ql) / h;
                    }
                });
            }
        }
    }
    void apply_flux_divergence(int l, bool predictor) override
    {
        const fdsamr::Level& L = m_reg.level(l);
        int mi = 0;
        for (amrex::MFIter mfi(lv[l].q[0]); mfi.isValid(); ++mfi, ++mi) {
            for (int kind = 0; kind < 2; ++kind) {
                if (lv[l].ov[kind].size() <= static_cast<std::size_t>(mi)) continue;
                for (const fdsrt::FluxOverride& o : lv[l].ov[kind][mi]) {
                    auto f = (kind == 0 ? lv[l].adv[o.dir] : lv[l].dif[o.dir]).array(mfi);
                    for (std::size_t q = 0; q < o.face.size(); ++q)
                        for (int n = 0; n < o.nscal && n < NC; ++n) f(o.face[q][0], o.face[q][1], o.face[q][2], n) = o.value[q * o.nscal + n];
                }
            }
        }
        amrex::MultiFab& q0 = lv[l].q[0]; amrex::MultiFab& qs = lv[l].q[1]; amrex::MultiFab& qo = lv[l].q[predictor ? 1 : 2];
        for (amrex::MFIter mfi(q0); mfi.isValid(); ++mfi) {
            auto a0 = q0.const_array(mfi); auto as = qs.const_array(mfi); auto ao = qo.array(mfi);
            auto ax = lv[l].adv[0].const_array(mfi); auto az = lv[l].adv[2].const_array(mfi); auto dx_ = lv[l].dif[0].const_array(mfi); auto dz_ = lv[l].dif[2].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                for (int n = 0; n < NC; ++n) {
                    const double div = ((ax(i + 1, j, k, n) + dx_(i + 1, j, k, n)) - (ax(i, j, k, n) + dx_(i, j, k, n))) / L.dx[0] +
                                       ((az(i, j, k + 1, n) + dz_(i, j, k + 1, n)) - (az(i, j, k, n) + dz_(i, j, k, n))) / L.dx[2];
                    ao(i, j, k, n) = predictor ? a0(i, j, k, n) - m_dt * div : 0.5 * (a0(i, j, k, n) + as(i, j, k, n) - m_dt * div);
                }
            });
        }
    }
    const amrex::MultiFab& stage_flux(int l, fdsrt::FluxKind kind, int d) const override { return kind == fdsrt::FluxKind::Adv ? lv[l].adv[d] : lv[l].dif[d]; }
    void set_flux_override(int l, fdsrt::FluxKind kind, const std::vector<std::vector<fdsrt::FluxOverride>>& lists) override { lv[l].ov[kind == fdsrt::FluxKind::Adv ? 0 : 1] = lists; }
    double local_dt(int) const override { return m_dt; }
    void run_divergence_part1(int) override {}
    double max_divergence_error(int l) const override
    {
        const fdsamr::Level& L = m_reg.level(l);
        const auto& F = m_reg.fields(l);
        double m = 0.0;
        for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
            auto u = F["U"].const_array(mfi); auto w = F["W"].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { m = std::max(m, std::abs((u(i + 1, j, k) - u(i, j, k)) / L.dx[0] + (w(i, j, k + 1) - w(i, j, k)) / L.dx[2])); });
        }
        amrex::ParallelAllReduce::Max(m, amrex::ParallelContext::CommunicatorSub());
        return m;
    }
    const fdsrt::FluxStageRunner& runner() const { return *m_runner; }

private:
    // ghost layer of q in x and z: coarse value injected (piecewise constant) on a fine level, then the same-level / periodic copy over it
    void fill_ghost(int l, amrex::MultiFab& q)
    {
        const fdsamr::Level& L = m_reg.level(l);
        if (l > 0) {
            const fdsamr::Level& C = m_reg.level(l - 1);
            const amrex::IntVect r = L.ref_ratio_from_parent;
            amrex::BoxList bl;
            for (int b = 0; b < static_cast<int>(L.ba.size()); ++b) { amrex::Box g = amrex::grow(L.ba[b], amrex::IntVect(1, 0, 1)); bl.push_back(amrex::coarsen(g, r)); }
            amrex::BoxArray cba(bl);
            amrex::MultiFab ct(cba, L.dm, NC, 0);
            ct.setVal(0.0);
            ct.ParallelCopy(lv[l - 1].q[&q == &lv[l].q[0] ? 0 : (&q == &lv[l].q[1] ? 1 : 2)], 0, 0, NC, 0, 0, C.geom.periodicity());
            for (amrex::MFIter mfi(q); mfi.isValid(); ++mfi) {
                auto a = q.array(mfi); auto c = ct.const_array(mfi);
                const amrex::Box vb = mfi.validbox();
                const amrex::Box gb = amrex::grow(vb, amrex::IntVect(1, 0, 1));
                amrex::LoopOnCpu(gb, [&](int i, int j, int k) {
                    if (vb.contains(i, j, k)) return;
                    for (int n = 0; n < NC; ++n) a(i, j, k, n) = c(amrex::coarsen(i, r[0]), amrex::coarsen(j, r[1]), amrex::coarsen(k, r[2]), n);
                });
            }
        }
        q.FillBoundary(L.geom.periodicity());
    }
    fdsamr::LevelRegistry& m_reg;
    fdsamr::RegistryTransfer& m_tr;
    double m_dt, m_D;
    bool m_overwrite = true;
    std::unique_ptr<fdsrt::FluxStageRunner> m_runner;
};

}  // namespace fdstest

#endif
