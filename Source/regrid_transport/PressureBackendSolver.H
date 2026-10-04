// PressureBackendSolver.H: binds the CompositePoissonSolver of the post-regrid projection (PostRegridProjection.H) to Role 2's composite MLMG
// (Source/pressure_backend: PressureProblem::levels, PressureWorkspace::rebuild, solve_pressure, face_gradient_composite). AMReX only.
//
// rebuild(levels)  -> PressureWorkspace::rebuild for the new hierarchy (layout copies only; the old MultiFabs may be destroyed at once).
// solve(...)       -> solve_pressure(problem, options, &workspace) with the per-level right-hand sides, then face_gradient_composite for the face gradients
//                     (single-cell directions get a zero gradient; the fine fluxes overwrite the coarse C/F faces, D-032 (3)).
// Domain faces: periodic where the level-0 geometry is periodic, Neumann elsewhere (Dirichlet and mixed open/closed faces are not built by the backend).
// Limits of the backend that apply here (frozen/composite-notes.md): refinement ratios 2 and 4 only; no masks or OBST; the gradient on the fine side of
// C/F faces is first order (up to about 12 % at the corners of ratio-4 patches), so gate numbers use ratio 2. A non-Ok status of the backend becomes
// SolveReport::ok == false with the backend's message (the velocity is then not touched).
// tol_abs of the solver interface is not used (the backend takes a relative tolerance only); a right-hand side of round-off size is skipped by the caller.
#ifndef FDSRT_PRESSUREBACKENDSOLVER_H_
#define FDSRT_PRESSUREBACKENDSOLVER_H_

#include <memory>

#include "PostRegridProjection.H"

namespace fdsrt {

class PressureBackendSolver : public CompositePoissonSolver {
public:
    PressureBackendSolver();
    ~PressureBackendSolver() override;
    const char* name() const override { return "pressure_backend composite MLMG"; }
    void rebuild(const std::vector<ProjectionLevel>& levels) override;
    SolveReport solve(const std::vector<const amrex::MultiFab*>& rhs, std::vector<std::array<amrex::MultiFab*, 3>>& grad, double tol_rel, double tol_abs, int max_iter) override;
    // Diagnostics of the last solve (backend's own composite true residual, relative L2).
    double last_residual_rel2() const;
    int workspace_rebuilds() const;      // number of rebuild() calls
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace fdsrt

#endif
