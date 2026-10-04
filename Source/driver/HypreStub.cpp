// Link stub for the HYPRE backend (pressure_backend/HypreBackend.cpp) in driver builds whose AMReX has no HYPRE: the pressure backend's
// workspace and composite solve name HypreSystem and make_hypre_backend(), so the symbols must exist; selecting HYPRE then reports
// "not built" instead of failing at link time. Compiled only when HypreBackend.cpp exists and AMReX was built without HYPRE (CMakeLists.txt).
#include "HypreBackend.H"
#include "PressureBackend.H"

namespace pb {

struct HypreSystem::Impl { std::string msg = "the HYPRE backend is not built in this driver (AMReX without HYPRE)"; };

HypreSystem::HypreSystem (std::vector<HypreLayoutLevel> const&, std::array<BC,6> const&, HypreOptions const&, bool, bool) : m_impl(new Impl) {}
HypreSystem::~HypreSystem () = default;
bool HypreSystem::ok () const { return false; }
std::string const& HypreSystem::message () const { return m_impl->msg; }
bool HypreSystem::matches (std::vector<HypreLayoutLevel> const&, std::array<BC,6> const&, HypreOptions const&, bool) const { return false; }
BackendStatus HypreSystem::solve (std::vector<amrex::MultiFab*> const&, std::vector<amrex::MultiFab const*> const&, double, int, bool)
{
    BackendStatus s; s.converged = false; return s;
}
double HypreSystem::setup_seconds () const { return 0.0; }
std::unique_ptr<PressureBackend> make_hypre_backend () { return nullptr; }

} // namespace pb
