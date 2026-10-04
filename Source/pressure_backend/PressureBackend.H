// PressureBackend.H: contract of a backend. A backend solves L phi = rhs for H on valid cells of a
// single uniform level and reports status. Mean removal, pins, gauge and checks are done by the caller.
#ifndef PB_PRESSURE_BACKEND_H_
#define PB_PRESSURE_BACKEND_H_

#include "PressureIface.H"
#include <cmath>
#include <memory>

namespace pb {

// Same domain, problem extent and periodicity.
inline bool same_geometry (amrex::Geometry const& a, amrex::Geometry const& b)
{
    if (a.Domain() != b.Domain()) { return false; }
    for (int d = 0; d < 3; ++d) {
        if (a.isPeriodic(d) != b.isPeriodic(d)) { return false; }
        if (std::abs(a.ProbLo(d) - b.ProbLo(d)) > 1.0e-12*a.ProbLength(d)) { return false; }
        if (std::abs(a.ProbHi(d) - b.ProbHi(d)) > 1.0e-12*a.ProbLength(d)) { return false; }
    }
    return true;
}

class PressureBackend {
public:
    virtual ~PressureBackend () = default;
    virtual const char* name () const = 0;
    // rhs: valid cells; phi: 1 ghost layer allocated, valid cells hold the initial guess on entry.
    virtual BackendStatus solve (PressureProblem const& p, PressureOptions const& o,
                                 amrex::MultiFab& phi, amrex::MultiFab const& rhs) = 0;
    // Cached set-up (FFT plans). plan_matches: the cached set-up was built for this problem's layout (BoxArray,
    // DistributionMapping, Geometry, boundary types). prepare: build it now (discarding an older one). Backends without a
    // cache keep the defaults.
    virtual bool plan_matches (PressureProblem const&) const { return false; }
    virtual void prepare (PressureProblem const&) {}
};

std::unique_ptr<PressureBackend> make_fft_backend ();
std::unique_ptr<PressureBackend> make_mlmg_backend ();
#ifdef PB_WITH_HYPRE
std::unique_ptr<PressureBackend> make_hypre_backend ();     // HypreBackend.cpp (only with PB_WITH_HYPRE)
#endif

} // namespace pb
#endif
