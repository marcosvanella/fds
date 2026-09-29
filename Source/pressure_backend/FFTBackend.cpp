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
    BackendStatus solve (PressureProblem const& p, PressureOptions const&,
                         MultiFab& phi, MultiFab const& rhs) override
    {
        Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> fbc{
            std::make_pair(fft_bc(p.bc[face_index(0,0)]), fft_bc(p.bc[face_index(0,1)])),
            std::make_pair(fft_bc(p.bc[face_index(1,0)]), fft_bc(p.bc[face_index(1,1)])),
            std::make_pair(fft_bc(p.bc[face_index(2,0)]), fft_bc(p.bc[face_index(2,1)]))};
        FFT::Poisson<MultiFab> fft(p.geom, fbc);
        fft.solve(phi, rhs);
        BackendStatus s;
        s.converged = true;     // direct solve
        return s;
    }
};
}

std::unique_ptr<PressureBackend> make_fft_backend () { return std::make_unique<FFTBackend>(); }

} // namespace pb
