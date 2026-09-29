// PressureBackend.H: contract of a backend. A backend solves L phi = rhs for H on valid cells of a
// single uniform level and reports status. Mean removal, pins, gauge and checks are done by the caller.
#ifndef PB_PRESSURE_BACKEND_H_
#define PB_PRESSURE_BACKEND_H_

#include "PressureIface.H"
#include <memory>

namespace pb {

class PressureBackend {
public:
    virtual ~PressureBackend () = default;
    virtual const char* name () const = 0;
    // rhs: valid cells; phi: 1 ghost layer allocated, valid cells hold the initial guess on entry.
    virtual BackendStatus solve (PressureProblem const& p, PressureOptions const& o,
                                 amrex::MultiFab& phi, amrex::MultiFab const& rhs) = 0;
};

std::unique_ptr<PressureBackend> make_fft_backend ();
std::unique_ptr<PressureBackend> make_mlmg_backend ();

} // namespace pb
#endif
