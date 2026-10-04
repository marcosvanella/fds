// MLMG backend: single level MLPoisson at setMaxOrder(2), default bottom solver (MLABecLaplacian with face coefficients for a masked
// problem, see solve_masked). Fluxes are never requested.
// The common-layer pin is ignored: singular problems rely on mean removal (before) and gauge (after), as P2.
#include "PressureBackend.H"
#include "CommonLayer.H"

#include <AMReX_MLABecLaplacian.H>
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

// Masked single-level solve (frozen/masked-notes.md): A phi = s with A = alpha - div(beta grad), beta = 1 on gas-gas faces and on the
// domain faces of gas cells, 0 on every face that touches a Solid or Known cell; alpha = sum of 2/dx^2 over the Known faces of a gas
// cell, sum_d 1/dx_d^2 on excluded cells (identity row, decoupled); s = -rhs + sum_known 2 g/dx^2 on gas cells, 0 elsewhere. This is
// -L phi = -rhs of the masked operator of CommonLayer.cpp (apply_operator) on the gas cells. rhs is zero on excluded cells.
BackendStatus solve_masked (PressureProblem const& p, PressureOptions const& o, MultiFab& phi, MultiFab const& rhs)
{
    const Box dom = p.geom.Domain();
    const auto dx = p.geom.CellSizeArray();
    const Real idx2[3] = {Real(1)/(dx[0]*dx[0]), Real(1)/(dx[1]*dx[1]), Real(1)/(dx[2]*dx[2])};
    std::array<BC,6> const eb = effective_bc(p.bc, dom);
    iMultiFab cls = class_ghost(p);
    MultiFab kv = known_ghost(p);

    MultiFab alpha(p.ba, p.dm, 1, 0), s(p.ba, p.dm, 1, 0);
    for (MFIter mfi(alpha); mfi.isValid(); ++mfi) {
        auto const& al = alpha.array(mfi);
        auto const& sr = s.array(mfi);
        auto const& r = rhs.const_array(mfi);
        auto const& cl = cls.const_array(mfi);
        auto const& g = kv.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            if (cl(i,j,k) != CellGas) { al(i,j,k) = idx2[0] + idx2[1] + idx2[2]; sr(i,j,k) = Real(0); return; }
            Real a = 0, rs = -r(i,j,k);
            bool coupled = false;     // any gas neighbour or Dirichlet domain face
            for (int d = 0; d < 3; ++d) {
                for (int side = 0; side < 2; ++side) {
                    int nb[3] = {i, j, k};
                    nb[d] += (side == 0) ? -1 : 1;
                    const int iv[3] = {i, j, k};
                    const bool edge = (side == 0) ? (iv[d] == dom.smallEnd(d)) : (iv[d] == dom.bigEnd(d));
                    if (edge && eb[face_index(d,side)] != BC::Periodic) { if (eb[face_index(d,side)] == BC::Dirichlet) { coupled = true; } continue; }
                    if (cl(nb[0],nb[1],nb[2]) == CellGas) { coupled = true; }
                    if (cl(nb[0],nb[1],nb[2]) == CellKnown) { a += Real(2)*idx2[d]; rs += Real(2)*idx2[d]*g(nb[0],nb[1],nb[2]); }
                }
            }
            if (!coupled && a == Real(0)) { a = idx2[0] + idx2[1] + idx2[2]; rs = Real(0); }   // isolated singular cell: identity row, phi = 0 (rhs is zero after mean removal)
            al(i,j,k) = a; sr(i,j,k) = rs;
        });
    }
    std::array<MultiFab,3> beta;
    for (int d = 0; d < 3; ++d) {
        BoxArray nba = amrex::convert(p.ba, IntVect::TheDimensionVector(d));
        beta[d].define(nba, p.dm, 1, 0);
        for (MFIter mfi(beta[d]); mfi.isValid(); ++mfi) {
            auto const& b = beta[d].array(mfi);
            auto const& cl = cls.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                int hi[3] = {i, j, k}, lo[3] = {i, j, k};
                lo[d] -= 1;                                   // face (i,j,k) separates cell lo (index-1) and cell hi (index)
                const int idx = hi[d];
                const bool lo_edge = (idx == dom.smallEnd(d)), hi_edge = (idx == dom.bigEnd(d) + 1);
                bool gas;
                if ((lo_edge || hi_edge) && eb[face_index(d, lo_edge ? 0 : 1)] != BC::Periodic) {
                    gas = lo_edge ? (cl(hi[0],hi[1],hi[2]) == CellGas) : (cl(lo[0],lo[1],lo[2]) == CellGas);
                } else {
                    gas = (cl(hi[0],hi[1],hi[2]) == CellGas) && (cl(lo[0],lo[1],lo[2]) == CellGas);
                }
                b(i,j,k) = gas ? Real(1) : Real(0);
            });
        }
    }
    LPInfo info;
    MLABecLaplacian mlab({p.geom}, {p.ba}, {p.dm}, info);
    mlab.setMaxOrder(2);
    Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(eb[face_index(0,0)]), lin_bc(eb[face_index(1,0)]), lin_bc(eb[face_index(2,0)])};
    Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(eb[face_index(0,1)]), lin_bc(eb[face_index(1,1)]), lin_bc(eb[face_index(2,1)])};
    mlab.setDomainBC(lo, hi);
    mlab.setScalars(Real(1), Real(1));
    mlab.setACoeffs(0, alpha);
    Array<MultiFab const*,AMREX_SPACEDIM> bptr{&beta[0], &beta[1], &beta[2]};
    mlab.setBCoeffs(0, bptr);
    // Excluded cells solve to zero (identity rows); keep the caller's initial guess on gas cells only.
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const& f = phi.array(mfi);
        auto const& cl = cls.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (cl(i,j,k) != CellGas) { f(i,j,k) = Real(0); } });
    }
    mlab.setLevelBC(0, &phi);
    MLMG mlmg(mlab);
    // Averaged coefficients make the coarse levels of a masked problem poor approximations: smooth harder and use BiCGStab+CG as the
    // bottom solver (plain BiCGStab stalls on sealed multi-component systems whose coarse operator couples the components).
    mlmg.setPreSmooth(8);
    mlmg.setPostSmooth(8);
    mlmg.setBottomSolver(MLMG::BottomSolver::bicgcg);
    mlmg.setVerbose(0);
    mlmg.setMaxIter(o.max_iter);
    mlmg.setConvergenceNormType(MLMGNormType::bnorm);
    mlmg.setThrowException(true);
    BackendStatus st;
    // MLMG measures the residual against |s|_inf, which contains the Known-cell terms 2 g/dx^2; the common layer measures against |rhs|.
    // Scale the relative tolerance so that both mean the same thing.
    Real rel = Real(o.tol_rel);
    {
        const Real ns = s.norm0(), nr = rhs.norm0();
        if (ns > Real(0) && nr > Real(0) && nr < ns) { rel *= nr / ns; }
    }
    try {
        mlmg.solve({&phi}, {&s}, Real(rel), Real(0.0));
        st.converged = true;
    } catch (std::exception const&) {
        st.converged = false;
    }
    st.iterations = mlmg.getNumIters();
    const Real b0 = mlmg.getInitRHS();
    st.own_residual = (b0 > 0) ? mlmg.getFinalResidual() / b0 : mlmg.getFinalResidual();
    return st;
}

class MLMGBackend final : public PressureBackend {
public:
    const char* name () const override { return "MLMG"; }
    BackendStatus solve (PressureProblem const& p, PressureOptions const& o,
                         MultiFab& phi, MultiFab const& rhs) override
    {
        if (is_masked(p)) { return solve_masked(p, o, phi, rhs); }
        LPInfo info;
        MLPoisson mlp({p.geom}, {p.ba}, {p.dm}, info);
        mlp.setMaxOrder(2);
        std::array<BC,6> const eb = effective_bc(p.bc, p.geom.Domain());
        Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(eb[face_index(0,0)]), lin_bc(eb[face_index(1,0)]), lin_bc(eb[face_index(2,0)])};
        Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(eb[face_index(0,1)]), lin_bc(eb[face_index(1,1)]), lin_bc(eb[face_index(2,1)])};
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
