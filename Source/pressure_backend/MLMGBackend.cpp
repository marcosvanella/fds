// MLMG backend: single level MLPoisson at setMaxOrder(2), default bottom solver (MLABecLaplacian with face coefficients for a masked
// problem, see solve_masked). Fluxes are never requested.
// The common-layer pin is ignored: singular problems rely on mean removal (before) and gauge (after), as P2.
#include "PressureBackend.H"
#include "CommonLayer.H"

#include <AMReX_BiCGStab.H>
#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_PCG.H>

namespace pb {

using namespace amrex;

namespace {
// E1 (docs/pressure/08 section 4.2): MLMG stops on a max-norm rule (|r|_inf / |b|_inf) while the common layer checks a 2-norm, and a
// V-cycle count is discrete, so a solve that stops "just in time" can end up to ~1.3 times above tol_rel in the checked norm
// (3 of 24 default solves in docs/pressure/07 section 15.2). Run MLMG to half the requested tolerance, as the HYPRE path runs its
// Krylov solver to 0.1*tol_rel; `converged` is then judged against the stricter value, never the looser one.
constexpr double kMlmgTolFactor = 0.5;
constexpr int kMlmgPolish = 3;          // extra solves of the 2-norm polish (unmasked path)

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
    // D7: averaged coefficients make the coarse levels of a masked problem poor approximations, so a plain V-cycle iteration can need
    // ~200 cycles (Dirichlet face plus a thin slab). GMRES with one V-cycle (8+8 sweeps, smoother bottom: a fixed linear preconditioner)
    // per iteration recovers the convergence. GMRES judges its own 2-norm residual against ||s||_2, so the caller's tol_rel, defined on
    // ||rhs||_2 by the common layer, is retargeted by ||rhs||_2/||s||_2 when the Known terms make s larger.
    mlmg.setPreSmooth(8);
    mlmg.setPostSmooth(8);
    mlmg.setVerbose(0);
    const Real ns = s.norm2(), nr = rhs.norm2();
    const Real rel = (ns > Real(0) && nr > Real(0) && nr < ns) ? Real(o.tol_rel) * nr / ns : Real(o.tol_rel);
    return mlmg_solve(mlmg, {&phi}, {&s}, o.mlmg_krylov == Krylov::Auto ? Krylov::GMRES : o.mlmg_krylov, rel, o.max_iter);
}

class MLMGBackend final : public PressureBackend {
public:
    const char* name () const override { return "MLMG"; }
    // Unmasked: the operator is built once per layout (the workspace keeps this object; P2 "set up once per regrid").
    bool plan_matches (PressureProblem const& p) const override
    {
        return m_mlmg && !is_masked(p) && m_ba == p.ba && m_dm == p.dm && m_bc == p.bc && same_geometry(m_geom, p.geom);
    }
    void prepare (PressureProblem const& p) override
    {
        m_mlmg.reset();
        LPInfo info;
        m_mlp = std::make_unique<MLPoisson>(Vector<Geometry>{p.geom}, Vector<BoxArray>{p.ba}, Vector<DistributionMapping>{p.dm}, info);
        m_mlp->setMaxOrder(2);
        std::array<BC,6> const eb = effective_bc(p.bc, p.geom.Domain());
        Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(eb[face_index(0,0)]), lin_bc(eb[face_index(1,0)]), lin_bc(eb[face_index(2,0)])};
        Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(eb[face_index(0,1)]), lin_bc(eb[face_index(1,1)]), lin_bc(eb[face_index(2,1)])};
        m_mlp->setDomainBC(lo, hi);
        m_mlp->setLevelBC(0, nullptr);   // homogeneous Dirichlet values; other BCs ignore the ghost data
        m_mlmg = std::make_unique<MLMG>(*m_mlp);
        m_mlmg->setVerbose(0);
        m_mlmg->setConvergenceNormType(MLMGNormType::bnorm);
        m_mlmg->setThrowException(true);
        m_ba = p.ba; m_dm = p.dm; m_bc = p.bc; m_geom = p.geom;
    }
    BackendStatus solve (PressureProblem const& p, PressureOptions const& o,
                         MultiFab& phi, MultiFab const& rhs) override
    {
        const double t0 = amrex::second();
        if (is_masked(p)) { return solve_masked(p, o, phi, rhs); }
        const bool reused = plan_matches(p);
        if (!reused) { prepare(p); }
        MLMG& mlmg = *m_mlmg;
        mlmg.setMaxIter(o.max_iter);
        const double t_setup = amrex::second() - t0;
        const Real tgt = Real(kMlmgTolFactor*o.tol_rel);
        if (o.mlmg_krylov != Krylov::Auto && o.mlmg_krylov != Krylov::None) {
            BackendStatus s = mlmg_solve(mlmg, {&phi}, {&rhs}, o.mlmg_krylov, tgt, o.max_iter);
            s.setup_seconds = reused ? 0.0 : t_setup; s.plan_reused = reused;
            return s;
        }
        BackendStatus s;
        s.method = "V-cycle";
        s.setup_seconds = reused ? 0.0 : t_setup; s.plan_reused = reused;
        // Polish toward residual_tol (proxy for residual_limit) when that is looser than half-tol; never chase below 0.5*tol_rel
        // when the caller asked for a loose solve (negative controls with tol_rel >> residual_tol must still warn).
        const Real polish_tgt = std::max(tgt, Real(o.residual_tol > 0.0 ? o.residual_tol : tgt));
        Real tol = tgt;
        int iters = 0;
        bool first_ok = false;
        try {
            mlmg.solve({&phi}, {&rhs}, tol, Real(0.0));
            iters = mlmg.getNumIters();
            first_ok = true;
            // 2-norm polish (E1): max-norm stop vs common-layer 2-norm; hoist residual MultiFab; keep first-solve status if polish fails.
            const Real bn2 = rhs.norm2(0);
            MultiFab res(p.ba, p.dm, 1, 0);
            for (int pass = 0; pass < kMlmgPolish && bn2 > Real(0); ++pass) {
                mlmg.compResidual({&res}, {&phi}, {&rhs});
                const Real r2 = res.norm2(0) / bn2;
                if (r2 <= polish_tgt) { break; }
                const Real rmax = res.norm0(0) / rhs.norm0(0);
                const Real q = (rmax > Real(0)) ? r2 / rmax : Real(1);
                tol = Real(0.8) * polish_tgt / std::max(q, Real(1.0e-3));
                mlmg.solve({&phi}, {&rhs}, tol, Real(0.0));
                iters += mlmg.getNumIters();
            }
            s.converged = true;
        } catch (std::exception const&) {
            s.converged = first_ok;   // polish (or a later pass) failed: keep Ok if the half-tol first solve met its target
            iters += mlmg.getNumIters();
        }
        s.iterations = iters;
        const Real b0 = mlmg.getInitRHS();
        s.own_residual = (b0 > 0) ? mlmg.getFinalResidual() / b0 : mlmg.getFinalResidual();
        s.solve_seconds = amrex::second() - t0 - t_setup;
        return s;
    }
private:
    std::unique_ptr<MLPoisson> m_mlp;
    std::unique_ptr<MLMG> m_mlmg;                 // refers to *m_mlp: declared after it, destroyed first
    BoxArray m_ba;
    DistributionMapping m_dm;
    Geometry m_geom;
    std::array<BC,6> m_bc{};
};
}

namespace {
// -L and -(one V-cycle) for PCG: the MLPoisson operator is negative (semi-)definite, PCG needs a positive one.
struct NegatedMLMG : GMRESMLMG {
    using GMRESMLMG::GMRESMLMG;
    void apply (VEC& lhs, VEC const& rhs) const { GMRESMLMG::apply(lhs, rhs); scale(lhs, RT(-1)); }
    void precond (VEC& lhs, VEC const& rhs) const { GMRESMLMG::precond(lhs, rhs); scale(lhs, RT(-1)); }
};
}

BackendStatus mlmg_solve (MLMG& mlmg, Vector<MultiFab*> const& x, Vector<MultiFab const*> const& b, Krylov k, Real tol, int max_iter, bool keep_bottom)
{
    BackendStatus st;
    const double t0 = amrex::second();
    if (k == Krylov::Auto || k == Krylov::None) {
        st.method = "V-cycle";
        // None (the study's explicit choice) stops at tol/10 like the Krylov solvers: the MLMG max-norm criterion can leave the true relative L2 residual
        // just above tol after a regrid. Auto keeps the default stop of the time loop.
        const Real t = (k == Krylov::None) ? Real(0.1) * tol : tol;
        try { mlmg.solve(x, b, t, Real(0)); st.converged = true; } catch (std::exception const&) {}
        st.iterations = mlmg.getNumIters();
        const Real b0 = mlmg.getInitRHS();
        st.own_residual = (b0 > 0) ? mlmg.getFinalResidual() / b0 : mlmg.getFinalResidual();
        st.solve_seconds = amrex::second() - t0;
        return st;
    }
    // Correction form of GMRESMLMG::solve: L cor = L x - b from cor = 0, then x -= cor. GMRESMLMG supplies the operator, the one-V-cycle
    // preconditioner (smoother bottom solver: a fixed linear map) and the vector algebra that every AMReX Krylov solver needs.
    using V = Vector<MultiFab>;
    GMRESMLMG op(mlmg);
    const BottomSolver bottom = mlmg.getBottomSolver();
    if (!keep_bottom) { mlmg.setBottomSolver(BottomSolver::smoother); }   // keep_bottom: the masked composite keeps its BoomerAMG bottom (pressure/13 ruling 6)
    V res = op.makeVecLHS(), cor = op.makeVecLHS();
    mlmg.apply(GetVecOfPtrs(res), x);
    for (int l = 0; l < int(x.size()); ++l) { MultiFab::Saxpy(res[l], Real(-1), *b[l], 0, 0, 1, 0); }
    const Real r0 = op.norm2(res);
    auto run = [&] (auto& s) {
        s.setMaxIters(max_iter);
        s.solve(cor, res, Real(0.1) * tol, Real(0));   // the recursive Krylov residual drifts from the true one: stop at tol/10, as the HYPRE solvers
        st.converged = (s.getStatus() == 0);
        st.iterations = s.getNumIters();
        st.own_residual = (r0 > 0) ? s.getResidualNorm() / r0 : s.getResidualNorm();
    };
    if (k == Krylov::GMRES) { st.method = "GMRES"; run(op.getGMRES()); }
    else if (k == Krylov::BiCGSTAB) { st.method = "BiCGSTAB"; BiCGStab<V,GMRESMLMG> s; s.define(op); run(s); }
    else { st.method = "PCG"; NegatedMLMG nop(mlmg); GMRESMLMG::scale(res, Real(-1)); PCG<V,NegatedMLMG> s; s.define(nop); run(s); }
    for (int l = 0; l < int(x.size()); ++l) { MultiFab::Saxpy(*x[l], Real(-1), cor[l], 0, 0, 1, 0); }
    mlmg.setBottomSolver(bottom);
    st.solve_seconds = amrex::second() - t0;
    return st;
}

std::unique_ptr<PressureBackend> make_mlmg_backend () { return std::make_unique<MLMGBackend>(); }

} // namespace pb
