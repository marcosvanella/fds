// PressureBackendSolver.cpp: see PressureBackendSolver.H.
#include "PressureBackendSolver.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "PressureIface.H"

namespace fdsrt {

struct PressureBackendSolver::Impl {
    pb::PressureWorkspace ws;
    std::vector<amrex::BoxArray> ba;
    std::vector<amrex::DistributionMapping> dm;
    std::vector<amrex::Geometry> geom;
    std::vector<amrex::IntVect> ratio;
    std::vector<amrex::MultiFab> phi;       // owned solution storage, valid cells + one ghost layer
    std::vector<amrex::MultiFab> rhs0;      // zero right-hand sides: layout carriers for rebuild
    std::array<pb::BC, 6> bc;
    double residual_rel2 = 0.0;
    int rebuilds = 0;

    pb::PressureProblem problem(const std::vector<const amrex::MultiFab*>* rhs)
    {
        pb::PressureProblem p;
        p.bc = bc;
        for (size_t l = 0; l < ba.size(); ++l) {
            pb::PressureLevel L;
            L.ba = ba[l]; L.dm = dm[l]; L.geom = geom[l]; L.ref_ratio = ratio[l];
            L.rhs = rhs ? (*rhs)[l] : &rhs0[l];
            L.phi = &phi[l];
            p.levels.push_back(L);
        }
        return p;
    }
};

PressureBackendSolver::PressureBackendSolver() : impl_(std::make_unique<Impl>()) {}
PressureBackendSolver::~PressureBackendSolver() = default;
double PressureBackendSolver::last_residual_rel2() const { return impl_->residual_rel2; }
int PressureBackendSolver::workspace_rebuilds() const { return impl_->rebuilds; }

void PressureBackendSolver::rebuild(const std::vector<ProjectionLevel>& levels)
{
    Impl& m = *impl_;
    m.ba.clear(); m.dm.clear(); m.geom.clear(); m.ratio.clear(); m.phi.clear(); m.rhs0.clear();
    for (const auto& L : levels) {
        m.ba.push_back(amrex::convert(L.vel[0]->boxArray(), amrex::IntVect(0)));
        m.dm.push_back(L.vel[0]->DistributionMap());
        m.geom.push_back(L.geom);
        m.ratio.push_back(L.ref_ratio);
    }
    for (size_t l = 0; l < m.ba.size(); ++l) {
        m.phi.emplace_back(m.ba[l], m.dm[l], 1, 1);
        m.phi.back().setVal(0.0);
        m.rhs0.emplace_back(m.ba[l], m.dm[l], 1, 0);
        m.rhs0.back().setVal(0.0);
    }
    for (int d = 0; d < 3; ++d) {
        const pb::BC b = levels[0].geom.isPeriodic(d) ? pb::BC::Periodic : pb::BC::Neumann;
        m.bc[pb::face_index(d, 0)] = b;
        m.bc[pb::face_index(d, 1)] = b;
    }
    ++m.rebuilds;
    std::string msg;
    const pb::Status st = m.ws.rebuild(m.problem(nullptr), &msg);
    if (st != pb::Status::Ok) amrex::Print() << "PressureBackendSolver::rebuild: " << pb::to_string(st) << ": " << msg << "\n";   // solve() reports the same reason
}

SolveReport PressureBackendSolver::solve(const std::vector<const amrex::MultiFab*>& rhs, std::vector<std::array<amrex::MultiFab*, 3>>& grad, double tol_rel, double /*tol_abs*/, int max_iter)
{
    Impl& m = *impl_;
    SolveReport r;
    for (auto& f : m.phi) f.setVal(0.0);
    pb::PressureProblem p = m.problem(&rhs);
    pb::PressureOptions o;
    o.backend = pb::BackendKind::MLMG;
    o.tol_rel = tol_rel;
    o.max_iter = max_iter;
    o.use_initial_guess = false;
    o.verbose = 0;
    const pb::PressureResult res = pb::solve_pressure(p, o, &m.ws);
    r.iterations = res.backend_status.iterations;
    m.residual_rel2 = res.residual_rel2;
    r.residual = res.backend_status.own_residual;
    if (res.status != pb::Status::Ok) {
        r.ok = false;
        r.message = std::string(pb::to_string(res.status)) + ": " + res.message;
        return r;
    }
    std::vector<std::array<amrex::MultiFab*, 3>> g = grad;
    const pb::PressureResult gr = pb::face_gradient_composite(p, g);
    if (gr.status != pb::Status::Ok) {
        r.ok = false;
        r.message = std::string("face_gradient_composite ") + pb::to_string(gr.status) + ": " + gr.message;
        return r;
    }
    r.ok = true;
    return r;
}

}  // namespace fdsrt
