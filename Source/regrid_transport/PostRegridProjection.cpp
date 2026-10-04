// PostRegridProjection.cpp: see PostRegridProjection.H.
#include "PostRegridProjection.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>

#include "FaceTransfer.H"

namespace fdsrt {

namespace {

amrex::MultiFab make_cell_mf(const ProjectionLevel& L)
{
    const amrex::BoxArray cba = amrex::convert(L.vel[0]->boxArray(), amrex::IntVect(0));
    return amrex::MultiFab(cba, L.vel[0]->DistributionMap(), 1, 0);
}

// out = div u - D on uncovered cells, 0 on covered ones; returns the local max of |out|; and the local volume-weighted sums of the right-hand side
double divergence_residual(const ProjectionLevel& L, amrex::MultiFab& out, double* sum, double* abssum)
{
    const double dx = L.geom.CellSize(0), dy = L.geom.CellSize(1), dz = L.geom.CellSize(2);
    const double vol = dx * dy * dz;
    double mx = 0.0, s = 0.0, as = 0.0;
    for (amrex::MFIter mfi(out); mfi.isValid(); ++mfi) {
        auto o = out.array(mfi);
        auto u = L.vel[0]->const_array(mfi);
        auto v = L.vel[1]->const_array(mfi);
        auto w = L.vel[2]->const_array(mfi);
        amrex::Array4<const double> d;
        if (L.D) d = L.D->const_array(mfi);
        amrex::Array4<const int> cv;
        if (L.covered) cv = L.covered->const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            if (L.covered && cv(i, j, k) != 0) { o(i, j, k) = 0.0; return; }
            double r = (u(i + 1, j, k) - u(i, j, k)) / dx + (v(i, j + 1, k) - v(i, j, k)) / dy + (w(i, j, k + 1) - w(i, j, k)) / dz;
            if (L.D) r -= d(i, j, k);
            o(i, j, k) = r;
            mx = std::max(mx, std::abs(r));
            s += vol * r;
            as += vol * std::abs(r);
        });
    }
    if (sum) *sum += s;
    if (abssum) *abssum += as;
    return mx;
}

}  // namespace

void average_down_velocity(const std::vector<ProjectionLevel>& levels)
{
    for (int l = static_cast<int>(levels.size()) - 1; l >= 1; --l)
        average_down_faces({levels[l].vel[0], levels[l].vel[1], levels[l].vel[2]}, {levels[l - 1].vel[0], levels[l - 1].vel[1], levels[l - 1].vel[2]}, levels[l].ref_ratio);
}

double composite_divergence_error(const std::vector<ProjectionLevel>& levels, std::vector<double>* by_level)
{
    double all = 0.0;
    if (by_level) by_level->assign(levels.size(), 0.0);
    for (size_t l = 0; l < levels.size(); ++l) {
        amrex::MultiFab r = make_cell_mf(levels[l]);
        double m = divergence_residual(levels[l], r, nullptr, nullptr);
        amrex::ParallelDescriptor::ReduceRealMax(m);
        if (by_level) (*by_level)[l] = m;
        all = std::max(all, m);
    }
    return all;
}

ProjectionReport project_after_regrid(std::vector<ProjectionLevel>& levels, CompositePoissonSolver* solver, const ProjectionOptions& opts)
{
    ProjectionReport rep;
    const int nl = static_cast<int>(levels.size());
    rep.solver = solver ? solver->name() : "none";
    rep.max_change_level.assign(nl, 0.0);
    average_down_velocity(levels);
    std::vector<amrex::MultiFab> rhs;
    rhs.reserve(nl);
    double sum = 0.0, abssum = 0.0;
    for (int l = 0; l < nl; ++l) {
        rhs.push_back(make_cell_mf(levels[l]));
        divergence_residual(levels[l], rhs[l], &sum, &abssum);
    }
    rep.div_before = composite_divergence_error(levels, &rep.div_before_level);
    rep.div_after = rep.div_before;
    rep.div_after_level = rep.div_before_level;
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(abssum);
    rep.rhs_sum_rel = abssum > 0.0 ? std::abs(sum) / abssum : 0.0;
    if (!opts.enabled) {
        rep.accepted = rep.div_after <= opts.accept_abs + opts.accept_rel * rep.div_before;
        return rep;
    }
    if (rep.div_before <= opts.accept_abs) {   // nothing to correct: a solve against a right-hand side of round-off size would only chase round-off
        rep.accepted = true;
        rep.message = "already within the bound, solve skipped";
        return rep;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solver != nullptr, "project_after_regrid: enabled without a solver");
    solver->rebuild(levels);
    std::vector<std::array<amrex::MultiFab, 3>> grad(nl);
    std::vector<std::array<amrex::MultiFab*, 3>> gp(nl);
    std::vector<const amrex::MultiFab*> rp(nl);
    for (int l = 0; l < nl; ++l) {
        rp[l] = &rhs[l];
        for (int d = 0; d < 3; ++d) {
            grad[l][d].define(levels[l].vel[d]->boxArray(), levels[l].vel[d]->DistributionMap(), 1, 0);
            grad[l][d].setVal(0.0);
            gp[l][d] = &grad[l][d];
        }
    }
    const SolveReport sr = solver->solve(rp, gp, opts.tol_rel, opts.tol_abs, opts.max_iter);
    rep.ran = true;
    rep.solved = sr.ok;
    rep.iterations = sr.iterations;
    rep.message = sr.message;
    if (!sr.ok) return rep;   // velocity untouched
    for (int l = 0; l < nl; ++l) {
        for (int d = 0; d < 3; ++d) {
            rep.max_change_level[l] = std::max(rep.max_change_level[l], grad[l][d].norminf(0));
            amrex::MultiFab::Subtract(*levels[l].vel[d], grad[l][d], 0, 0, 1, 0);
        }
    }
    average_down_velocity(levels);
    rep.div_after = composite_divergence_error(levels, &rep.div_after_level);
    rep.accepted = rep.div_after <= opts.accept_abs + opts.accept_rel * rep.div_before;
    return rep;
}

}  // namespace fdsrt
