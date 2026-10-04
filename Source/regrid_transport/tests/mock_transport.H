// mock_transport.H: a small conservative scalar-transport code that speaks the FluxAccess protocol (RegridInterface.H) over a multi-level hierarchy. It is NOT FDS: it
// has the same stage structure as one SSPRK2 step of FDS (predictor, corrector; advective face flux = face value x face velocity; diffusive face flux; DEL_RHO_D_DEL_Z
// lagged by one stage; same time step on all levels), but the kernels are a two-point centred scheme. Its purpose is to test, without the driver's fine-level binding,
// everything Role 3 owns in the interface flux overwrite: the override lists (FluxOverrideOps), the finest-first ordering and clearing (FluxStageRunner), the shared faces
// between boxes, periodic wrap, hidden directions, three levels, and composite conservation at round-off with the overwrite on and its failure with it off.
// The cell update and the divergence of the diffusive fluxes use the override lists exactly as the driver's hooks do (replace the listed face value for all components).
#ifndef FDSRT_TEST_MOCK_TRANSPORT_H_
#define FDSRT_TEST_MOCK_TRANSPORT_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <cmath>
#include <functional>
#include <vector>

#include "FluxStageRunner.H"
#include "RegridInterface.H"

namespace mock {

using fdsrt::FluxAccess;
using fdsrt::FluxKind;
using fdsrt::FluxOverride;

struct Level {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::IntVect ratio{1, 1, 1};
    amrex::MultiFab S[3];                 // 0 = q (start of step), 1 = qs (predictor result), 2 = qn (corrector result); 1 ghost cell
    amrex::MultiFab vel[3];               // face velocity, nodal in its own direction
    amrex::MultiFab adv[3], dif[3];       // stage face fluxes (kg/(m2 s) of the transported quantity)
    amrex::MultiFab ddz;                  // DEL_RHO_D_DEL_Z, lagged
    std::vector<std::vector<FluxOverride>> ov[2];   // [0] ADV, [1] DIF; one vector per local box
    amrex::iMultiFab covered;             // 1 where level+1 covers the cell
    double dx[3] = {1, 1, 1};
};

class Transport : public FluxAccess {
public:
    // levels: geometry/ba/ratio filled by the caller; init() allocates and fills.
    std::vector<Level> lv;
    int nc = 1;
    double D = 0.0;

    void init(int ncomp, double diff,
              const std::function<double(double, double, double, int)>& q0,
              const std::function<std::array<double, 3>(int level, int dir, double x, double y, double z, const double* dx)>& vel_face)
    {
        nc = ncomp;
        D = diff;
        for (int l = 0; l < static_cast<int>(lv.size()); ++l) {
            Level& L = lv[l];
            if (L.dm.empty()) L.dm = amrex::DistributionMapping(L.ba);
            for (int d = 0; d < 3; ++d) L.dx[d] = L.geom.CellSize(d);
            for (int s = 0; s < 3; ++s) L.S[s].define(L.ba, L.dm, nc, 1);
            for (int d = 0; d < 3; ++d) {
                amrex::IntVect nodal(0); nodal[d] = 1;
                amrex::BoxArray nba = amrex::convert(L.ba, nodal);
                L.vel[d].define(nba, L.dm, 1, 0);
                L.adv[d].define(nba, L.dm, nc, 0);
                L.dif[d].define(nba, L.dm, nc, 0);
            }
            L.ddz.define(L.ba, L.dm, nc, 0);
            for (int s = 0; s < 3; ++s) L.S[s].setVal(0.0);
            L.ddz.setVal(0.0);
            L.ov[0].clear(); L.ov[1].clear();
            for (amrex::MFIter mfi(L.S[0]); mfi.isValid(); ++mfi) { L.ov[0].emplace_back(); L.ov[1].emplace_back(); }
            for (amrex::MFIter mfi(L.S[0]); mfi.isValid(); ++mfi) {
                auto a = L.S[0].array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    const double x = L.geom.ProbLo(0) + (i + 0.5) * L.dx[0], y = L.geom.ProbLo(1) + (j + 0.5) * L.dx[1], z = L.geom.ProbLo(2) + (k + 0.5) * L.dx[2];
                    for (int n = 0; n < nc; ++n) a(i, j, k, n) = q0(x, y, z, n);
                });
            }
            for (int d = 0; d < 3; ++d)
                for (amrex::MFIter mfi(L.vel[d]); mfi.isValid(); ++mfi) {
                    auto v = L.vel[d].array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                        // lower corner of the face: the face index is the low face of cell (i,j,k); the position of the face centre is cell centre minus half a cell in d
                        double p[3] = {L.geom.ProbLo(0) + (i + 0.5) * L.dx[0], L.geom.ProbLo(1) + (j + 0.5) * L.dx[1], L.geom.ProbLo(2) + (k + 0.5) * L.dx[2]};
                        p[d] -= 0.5 * L.dx[d];
                        v(i, j, k) = vel_face(l, d, p[0], p[1], p[2], L.dx)[d];
                    });
                }
        }
        update_masks();
        // ghosts and the lagged diffusion term of the start state
        average_down(0);
        fill_ghosts(0);
        for (size_t l = 0; l < lv.size(); ++l) readout_dif(static_cast<int>(l), 0);
        runner->set_overrides(*this, FluxKind::Dif);
        for (size_t l = 0; l < lv.size(); ++l) compute_ddz(static_cast<int>(l));
    }

    fdsrt::FluxStageRunner* runner = nullptr;

    void update_masks()
    {
        for (size_t l = 0; l < lv.size(); ++l) {
            Level& L = lv[l];
            L.covered.define(L.ba, L.dm, 1, 0);
            L.covered.setVal(0);
            if (l + 1 < lv.size()) {
                amrex::BoxArray cf = lv[l + 1].ba;
                cf.coarsen(lv[l + 1].ratio);
                for (amrex::MFIter mfi(L.covered); mfi.isValid(); ++mfi) {
                    auto m = L.covered.array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (cf.contains(amrex::IntVect(i, j, k))) m(i, j, k) = 1; });
                }
            }
        }
    }

    void average_down(int s)
    {
        for (int l = static_cast<int>(lv.size()) - 2; l >= 0; --l)
            amrex::average_down(lv[l + 1].S[s], lv[l].S[s], 0, nc, lv[l + 1].ratio);
    }

    // ghost cells (1 layer): piecewise-constant from the coarser level where no same-level cell exists, same-level copy and periodic images otherwise
    void fill_ghosts(int s)
    {
        for (size_t l = 0; l < lv.size(); ++l) {
            Level& L = lv[l];
            if (l > 0) {
                const Level& C = lv[l - 1];
                amrex::MultiFab cg(C.ba, C.dm, nc, 1);   // coarse level with its own ghosts so that a periodic parent can be read
                amrex::MultiFab::Copy(cg, C.S[s], 0, 0, nc, 0);
                cg.FillBoundary(C.geom.periodicity());
                // value of the parent cell for every ghost cell of every fine box (parent array gathered on the fine boxes' coarsened, grown box)
                amrex::BoxArray cba = L.ba;
                cba.coarsen(L.ratio);
                amrex::MultiFab par(cba, L.dm, nc, 1);
                par.ParallelCopy(cg, 0, 0, nc, 1, 1, C.geom.periodicity());
                for (amrex::MFIter mfi(L.S[s]); mfi.isValid(); ++mfi) {
                    auto f = L.S[s].array(mfi);
                    auto p = par.const_array(mfi);
                    const amrex::Box vb = mfi.validbox();
                    amrex::LoopOnCpu(amrex::grow(vb, 1), [&](int i, int j, int k) {
                        if (vb.contains(i, j, k)) return;
                        const int ci = amrex::coarsen(i, L.ratio[0]), cj = amrex::coarsen(j, L.ratio[1]), ck = amrex::coarsen(k, L.ratio[2]);
                        for (int n = 0; n < nc; ++n) f(i, j, k, n) = p(ci, cj, ck, n);
                    });
                }
            }
            L.S[s].FillBoundary(L.geom.periodicity());   // same-level neighbours and periodic images overwrite the injected values
        }
    }

    // ---- stage pieces ------------------------------------------------------------------------------------------------------------------------------
    void readout_adv(int l, int s)
    {
        Level& L = lv[l];
        for (int d = 0; d < 3; ++d)
            for (amrex::MFIter mfi(L.adv[d]); mfi.isValid(); ++mfi) {
                auto F = L.adv[d].array(mfi);
                auto q = L.S[s].const_array(mfi);
                auto v = L.vel[d].const_array(mfi);
                amrex::IntVect e(0); e[d] = 1;
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    for (int n = 0; n < nc; ++n) F(i, j, k, n) = (0.5 * (q(i - e[0], j - e[1], k - e[2], n) + q(i, j, k, n))) * v(i, j, k);
                });
            }
    }
    void readout_dif(int l, int s)
    {
        Level& L = lv[l];
        for (int d = 0; d < 3; ++d)
            for (amrex::MFIter mfi(L.dif[d]); mfi.isValid(); ++mfi) {
                auto F = L.dif[d].array(mfi);
                auto q = L.S[s].const_array(mfi);
                amrex::IntVect e(0); e[d] = 1;
                const double r = D / L.dx[d];
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    for (int n = 0; n < nc; ++n) F(i, j, k, n) = -r * (q(i, j, k, n) - q(i - e[0], j - e[1], k - e[2], n));
                });
            }
    }
    // face array with the override list of `kind` applied (copy; the read-out stays untouched)
    amrex::MultiFab with_override(const amrex::MultiFab& F, int d, int kind, int l) const
    {
        amrex::MultiFab T(F.boxArray(), F.DistributionMap(), nc, 0);
        amrex::MultiFab::Copy(T, F, 0, 0, nc, 0);
        int ib = 0;
        for (amrex::MFIter mfi(T); mfi.isValid(); ++mfi, ++ib) {
            auto a = T.array(mfi);
            for (const FluxOverride& o : lv[l].ov[kind][ib]) {
                if (o.dir != d) continue;
                for (size_t f = 0; f < o.face.size(); ++f)
                    for (int n = 0; n < nc; ++n) a(o.face[f][0], o.face[f][1], o.face[f][2], n) = o.value[f * nc + n];
            }
        }
        return T;
    }
    void compute_ddz(int l)
    {
        Level& L = lv[l];
        L.ddz.setVal(0.0);
        for (int d = 0; d < 3; ++d) {
            amrex::MultiFab T = with_override(L.dif[d], d, 1, l);
            amrex::IntVect e(0); e[d] = 1;
            for (amrex::MFIter mfi(L.ddz); mfi.isValid(); ++mfi) {
                auto z = L.ddz.array(mfi);
                auto F = T.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    for (int n = 0; n < nc; ++n) z(i, j, k, n) += -(F(i + e[0], j + e[1], k + e[2], n) - F(i, j, k, n)) / L.dx[d];
                });
            }
        }
    }
    // out = a*base + b*(in + dt*(-div ADV + ddz)), cell update of one stage on level l
    void update_cells(int l, int in, int base, int out, double a, double b, double dt)
    {
        Level& L = lv[l];
        amrex::MultiFab R(L.ba, L.dm, nc, 0);
        R.setVal(0.0);
        for (int d = 0; d < 3; ++d) {
            amrex::MultiFab T = with_override(L.adv[d], d, 0, l);
            amrex::IntVect e(0); e[d] = 1;
            for (amrex::MFIter mfi(R); mfi.isValid(); ++mfi) {
                auto r = R.array(mfi);
                auto F = T.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    for (int n = 0; n < nc; ++n) r(i, j, k, n) += -(F(i + e[0], j + e[1], k + e[2], n) - F(i, j, k, n)) / L.dx[d];
                });
            }
        }
        for (amrex::MFIter mfi(R); mfi.isValid(); ++mfi) {
            auto r = R.const_array(mfi);
            auto qi = L.S[in].const_array(mfi);
            auto qb = L.S[base].const_array(mfi);
            auto qo = L.S[out].array(mfi);
            auto z = L.ddz.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                for (int n = 0; n < nc; ++n) qo(i, j, k, n) = a * qb(i, j, k, n) + b * (qi(i, j, k, n) + dt * (r(i, j, k, n) + z(i, j, k, n)));
            });
        }
    }

    // one SSPRK2 step with the stage structure of FDS: after each stage, average-down and ghost fill, then the diffusive read-out, its overrides and DEL_RHO_D_DEL_Z
    void step(double dt)
    {
        const int top = static_cast<int>(lv.size());
        average_down(0);
        fill_ghosts(0);
        for (int l = 0; l < top; ++l) readout_adv(l, 0);
        runner->set_overrides(*this, FluxKind::Adv);
        for (int l = top - 1; l >= 0; --l) update_cells(l, 0, 0, 1, 0.0, 1.0, dt);
        average_down(1);
        fill_ghosts(1);
        for (int l = 0; l < top; ++l) readout_dif(l, 1);
        runner->set_overrides(*this, FluxKind::Dif);
        for (int l = top - 1; l >= 0; --l) compute_ddz(l);
        // corrector
        for (int l = 0; l < top; ++l) readout_adv(l, 1);
        runner->set_overrides(*this, FluxKind::Adv);
        for (int l = top - 1; l >= 0; --l) update_cells(l, 1, 0, 2, 0.5, 0.5, dt);   // qn = 0.5 q + 0.5 (qs + dt L(qs))
        average_down(2);
        fill_ghosts(2);
        for (int l = 0; l < top; ++l) readout_dif(l, 2);
        runner->set_overrides(*this, FluxKind::Dif);
        for (int l = top - 1; l >= 0; --l) compute_ddz(l);
        for (int l = 0; l < top; ++l) amrex::MultiFab::Copy(lv[l].S[0], lv[l].S[2], 0, 0, nc, 1);
    }

    // composite integral of component n: uncovered cells of every level times their volume (long double, summed in a fixed order per rank, then over ranks)
    long double composite(int n, int s = 0) const
    {
        long double sum = 0;
        for (size_t l = 0; l < lv.size(); ++l) {
            const Level& L = lv[l];
            const double vol = L.dx[0] * L.dx[1] * L.dx[2];
            for (amrex::MFIter mfi(L.S[s]); mfi.isValid(); ++mfi) {
                auto q = L.S[s].const_array(mfi);
                auto m = L.covered.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (!m(i, j, k)) sum += static_cast<long double>(q(i, j, k, n)) * vol; });
            }
        }
        // exact enough reduction over ranks: every rank's long double sum is split into hi + lo doubles, gathered, and added in long double in rank order
        const double hi = static_cast<double>(sum), lo = static_cast<double>(sum - static_cast<long double>(hi));
        std::vector<double> ah(amrex::ParallelDescriptor::NProcs()), al(amrex::ParallelDescriptor::NProcs());
        amrex::ParallelAllGather::AllGather(hi, ah.data(), amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllGather::AllGather(lo, al.data(), amrex::ParallelContext::CommunicatorSub());
        long double tot = 0;
        for (size_t i = 0; i < ah.size(); ++i) tot += static_cast<long double>(ah[i]) + static_cast<long double>(al[i]);
        return tot;
    }
    double max_abs_diff_from(double v, int n) const   // max |q - v| over the uncovered cells (all levels)
    {
        double mx = 0;
        for (size_t l = 0; l < lv.size(); ++l) {
            const Level& L = lv[l];
            for (amrex::MFIter mfi(L.S[0]); mfi.isValid(); ++mfi) {
                auto q = L.S[0].const_array(mfi);
                auto m = L.covered.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (!m(i, j, k)) mx = std::max(mx, std::abs(q(i, j, k, n) - v)); });
            }
        }
        amrex::ParallelAllReduce::Max(mx, amrex::ParallelContext::CommunicatorSub());
        return mx;
    }

    // ---- FluxAccess ---------------------------------------------------------------------------------------------------------------------------------
    void compute_stage_fluxes(int, bool) override { amrex::Abort("mock: the stage pieces are called directly (see step())"); }
    void apply_flux_divergence(int, bool) override { amrex::Abort("mock: the stage pieces are called directly (see step())"); }
    const amrex::MultiFab& stage_flux(int l, FluxKind kind, int d) const override { return kind == FluxKind::Adv ? lv[l].adv[d] : lv[l].dif[d]; }
    void set_flux_override(int l, FluxKind kind, const std::vector<std::vector<FluxOverride>>& per_box) override
    {
        auto& dst = lv[l].ov[kind == FluxKind::Adv ? 0 : 1];
        const size_t nb = dst.size();
        dst.assign(nb, {});
        for (size_t i = 0; i < per_box.size() && i < nb; ++i) dst[i] = per_box[i];
    }
    double local_dt(int) const override { return 1.0e300; }
    void run_divergence_part1(int) override { amrex::Abort("mock: not used"); }
    double max_divergence_error(int) const override { return 0.0; }
};

}  // namespace mock

#endif
