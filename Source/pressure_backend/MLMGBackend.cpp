// MLMG backend: single level MLPoisson at setMaxOrder(2), default bottom solver. Fluxes are never requested.
// The common-layer pin is ignored: singular problems rely on mean removal (before) and gauge (after), as P2.
#include "PressureBackend.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

namespace pb {

using namespace amrex;

namespace {
LinOpBCType lin_bc (BC b)
{
    switch (b) {
    case BC::Neumann: return LinOpBCType::Neumann;
    case BC::Periodic: return LinOpBCType::Periodic;
    case BC::Dirichlet: return LinOpBCType::Dirichlet;
    }
    return LinOpBCType::Neumann;
}

class MLMGBackend final : public PressureBackend {
public:
    const char* name () const override { return "MLMG"; }
    BackendStatus solve (PressureProblem const& p, PressureOptions const& o,
                         MultiFab& phi, MultiFab const& rhs) override
    {
        LPInfo info;
        MLPoisson mlp({p.geom}, {p.ba}, {p.dm}, info);
        mlp.setMaxOrder(2);
        Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(p.bc[face_index(0,0)]), lin_bc(p.bc[face_index(1,0)]), lin_bc(p.bc[face_index(2,0)])};
        Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(p.bc[face_index(0,1)]), lin_bc(p.bc[face_index(1,1)]), lin_bc(p.bc[face_index(2,1)])};
        mlp.setDomainBC(lo, hi);
        mlp.setLevelBC(0, &phi);      // homogeneous Dirichlet values; other BCs ignore the ghost data
        MLMG mlmg(mlp);
        mlmg.setVerbose(0);
        mlmg.setMaxIter(o.max_iter);
        mlmg.setConvergenceNormType(MLMGNormType::bnorm);
        mlmg.setThrowException(true);
        BackendStatus s;
        try {
            mlmg.solve({&phi}, {&rhs}, Real(o.tol_rel), Real(0.0));
            s.converged = true;
        } catch (std::exception const&) {
            s.converged = false;
        }
        s.iterations = mlmg.getNumIters();
        const Real b0 = mlmg.getInitRHS();
        s.own_residual = (b0 > 0) ? mlmg.getFinalResidual() / b0 : mlmg.getFinalResidual();
        return s;
    }
};
}

std::unique_ptr<PressureBackend> make_mlmg_backend () { return std::make_unique<MLMGBackend>(); }

} // namespace pb
