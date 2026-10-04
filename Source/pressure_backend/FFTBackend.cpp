// FFT backend: amrex::FFT::Poisson (single level, box domain, uniform closed or uniform open faces).
#include "PressureBackend.H"

#include <AMReX_FFT_Poisson.H>

namespace pb {

using namespace amrex;

namespace {
FFT::Boundary fft_bc (BC b)
{
    switch (b) {
    case BC::Neumann: return FFT::Boundary::even;
    case BC::Periodic: return FFT::Boundary::periodic;
    case BC::Dirichlet: return FFT::Boundary::odd;
    }
    return FFT::Boundary::even;
}

class FFTBackend final : public PressureBackend {
public:
    const char* name () const override { return "FFT"; }

    // Key of the cached plan: BoxArray, DistributionMapping, Geometry and the six boundary types.
    bool plan_matches (PressureProblem const& p) const override
    {
        return m_plan && m_bc == effective_bc(p.bc, p.geom.Domain()) && m_ba == p.ba && m_dm == p.dm && same_geometry(m_geom, p.geom);
    }
    void prepare (PressureProblem const& p) override
    {
        m_plan.reset();
        m_ba = p.ba; m_dm = p.dm; m_geom = p.geom; m_bc = effective_bc(p.bc, p.geom.Domain());
        std::array<BC,6> const& eb = m_bc;
        Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> fbc{
            std::make_pair(fft_bc(eb[face_index(0,0)]), fft_bc(eb[face_index(0,1)])),
            std::make_pair(fft_bc(eb[face_index(1,0)]), fft_bc(eb[face_index(1,1)])),
            std::make_pair(fft_bc(eb[face_index(2,0)]), fft_bc(eb[face_index(2,1)]))};
        m_plan = std::make_unique<FFT::Poisson<MultiFab>>(p.geom, fbc);
    }
    BackendStatus solve (PressureProblem const& p, PressureOptions const&,
                         MultiFab& phi, MultiFab const& rhs) override
    {
        BackendStatus s;
        s.plan_reused = plan_matches(p);
        if (!s.plan_reused) { prepare(p); }
        m_plan->solve(phi, rhs);
        s.converged = true;     // direct solve
        return s;
    }
private:
    std::unique_ptr<FFT::Poisson<MultiFab>> m_plan;
    BoxArray m_ba;
    DistributionMapping m_dm;
    Geometry m_geom;
    std::array<BC,6> m_bc{};
};
}

std::unique_ptr<PressureBackend> make_fft_backend () { return std::make_unique<FFTBackend>(); }

} // namespace pb
