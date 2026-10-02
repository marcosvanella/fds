// fds_stubs.cpp: stand-ins for the FDS Fortran entry points that GhostExchange.cpp declares, so the tests can link the real BcStep::exchange (same-level fill and
// the coarse-fine hook) without the FDS objects. None of them may be called by these tests; each aborts with its name if it is.
#include <cstdio>
#include <cstdlib>

#define STUB(name, ...) extern "C" __VA_ARGS__ { std::fprintf(stderr, "fds_stubs: %s called (not available without the FDS objects)\n", #name); std::abort(); }
STUB(fds_g_fill_om, int fds_g_fill_om(int, int, int, const int*, const int*, int, const double*))
STUB(fds_g_phase, void fds_g_phase(int))
STUB(fds_g_match, void fds_g_match(int))
STUB(fds_p_save_uvw, void fds_p_save_uvw(int, int))
STUB(fds_g_wall_bc, void fds_g_wall_bc(double, double, int))
STUB(fds_g_velocity_bc, void fds_g_velocity_bc(double, int, int))
STUB(fds_g_viscosity_bc, void fds_g_viscosity_bc(int, int))
STUB(fds_g_mu_edges, void fds_g_mu_edges(int))
STUB(fds_g_mu_edges_dom, void fds_g_mu_edges_dom(int, int, int))
