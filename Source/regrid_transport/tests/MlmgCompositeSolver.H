// MlmgCompositeSolver.H: MOCK composite Poisson solver for the post-regrid projection tests (amrex::MLPoisson + MLMG, periodic domain or homogeneous Neumann).
// This stands in for Role 2's composite MLMG (PressureWorkspace::rebuild + PressureLevel vector) until it lands; it is NOT the production pressure path.
#ifndef FDSRT_TESTS_MLMGCOMPOSITESOLVER_H_
#define FDSRT_TESTS_MLMGCOMPOSITESOLVER_H_

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#include <cstdlib>
#include <exception>
#include <memory>

#include "PostRegridProjection.H"

namespace fdsrt {

class MlmgCompositeSolver : public CompositePoissonSolver {
public:
    const char* name() const override { return "mock amrex::MLMG"; }

    void rebuild(const std::vector<ProjectionLevel>& levels) override
    {
        mlmg_.reset();
        op_.reset();
        phi_.clear();
        const int nl = static_cast<int>(levels.size());
        amrex::Vector<amrex::Geometry> geoms;
        amrex::Vector<amrex::BoxArray> grids;
        amrex::Vector<amrex::DistributionMapping> dms;
        int hidden = -1;
        for (int l = 0; l < nl; ++l) {
            geoms.push_back(levels[l].geom);
            grids.push_back(amrex::convert(levels[l].vel[0]->boxArray(), amrex::IntVect(0)));
            dms.push_back(levels[l].vel[0]->DistributionMap());
        }
        for (int d = 0; d < 3; ++d)
            if (levels[0].geom.Domain().length(d) == 1) hidden = d;
        amrex::LPInfo info;
        info.setMaxCoarseningLevel(30);
        if (hidden >= 0) info.setHiddenDirection(hidden);
        hidden_ = hidden;
        op_ = std::make_unique<amrex::MLPoisson>(geoms, grids, dms, info);
        amrex::Array<amrex::LinOpBCType, 3> lo, hi;
        for (int d = 0; d < 3; ++d) {
            const bool per = levels[0].geom.isPeriodic(d);
            lo[d] = hi[d] = per ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
        }
        op_->setDomainBC(lo, hi);
        phi_.resize(nl);
        for (int l = 0; l < nl; ++l) {
            phi_[l].define(grids[l], dms[l], 1, 1);
            phi_[l].setVal(0.0);
            op_->setLevelBC(l, &phi_[l]);
        }
        mlmg_ = std::make_unique<amrex::MLMG>(*op_);
        mlmg_->setVerbose(std::getenv("RT_MLMG_VERBOSE") ? 2 : 0);
        mlmg_->setThrowException(true);   // a solver that does not converge is reported, not aborted
        mlmg_->setBottomVerbose(std::getenv("RT_MLMG_VERBOSE") ? 1 : 0);
    }

    SolveReport solve(const std::vector<const amrex::MultiFab*>& rhs, std::vector<std::array<amrex::MultiFab*, 3>>& grad, double tol_rel, double tol_abs, int max_iter) override
    {
        SolveReport r;
        amrex::Vector<amrex::MultiFab*> sol;
        amrex::Vector<const amrex::MultiFab*> rh;
        for (size_t l = 0; l < phi_.size(); ++l) { phi_[l].setVal(0.0); sol.push_back(&phi_[l]); rh.push_back(rhs[l]); }
        mlmg_->setMaxIter(max_iter);
        double res = -1.0;
        try {
            res = mlmg_->solve(sol, rh, tol_rel, tol_abs);
        } catch (std::exception const& e) {
            r.ok = false;
            r.iterations = mlmg_->getNumIters();
            r.message = e.what();
            return r;
        }
        amrex::Vector<amrex::Array<amrex::MultiFab*, 3>> g;
        for (auto& a : grad) g.push_back(a);
        mlmg_->getGradSolution(g);
        // the gradient along a hidden (single-cell) direction is zero by definition; MLMG's own value there is not (it reads its ghost cells)
        if (hidden_ >= 0)
            for (auto& a : grad) a[hidden_]->setVal(0.0);
        r.ok = true;
        r.residual = res;
        r.iterations = mlmg_->getNumIters();
        return r;
    }

private:
    std::unique_ptr<amrex::MLPoisson> op_;
    std::unique_ptr<amrex::MLMG> mlmg_;
    std::vector<amrex::MultiFab> phi_;
    int hidden_ = -1;
};

}  // namespace fdsrt

#endif
