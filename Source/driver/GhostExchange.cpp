// GhostExchange.cpp: see GhostExchange.H. Kernel-facing rules (M2a): (a) passive scalars are extra components of ZZ/ZZS (Fields.cpp), all
// components move; (b) only uniform Cartesian metrics are used, R(I)/RRN(I) are dropped.
#include "GhostExchange.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <map>

extern "C" {
int fds_g_fill_om(int nm, int nom, int which, const int* lb, const int* ext, int nc, const double* data);
void fds_g_phase(int pred);
void fds_g_match(int nm);
void fds_p_save_uvw(int nm, int pred);
void fds_g_wall_bc(double t, double dt, int nm);
void fds_g_velocity_bc(double t, int nm, int est);
void fds_g_viscosity_bc(int nm, int est);
void fds_g_mu_edges(int nm);
}

namespace fdsamr {

std::vector<std::string> exchange_fields(int code, bool predictor)
{
    switch (code) {
    // TMP and RSUM: the ghost layer of a cell across a box interface. FDS gets them from WALL_BC/ASSIGN_GHOST_VALUE (OMESH average, patch 0004: skipped when the ghosts
    // are filled externally); the driver fills them from the neighbouring box like the other scalars and keeps the interface walls out of WALL_BC (TimeLoop::iface).
    case 1: return {"RHOS", "ZZS", "MU", "KRES", "D", "TMP", "RSUM"};
    case 4: return {"RHO", "ZZ", "MU", "KRES", "DS", "TMP", "RSUM"};
    case 3: return {"US", "VS", "WS", "HS"};
    case 6: return {"U", "V", "W", "H"};
    case 5: return {"FVX", "FVY", "FVZ", predictor ? "H" : "HS"};
    default: return {};
    }
}

void ghost_exchange(Fields& F, int code, bool predictor)
{
    for (const auto& n : exchange_fields(code, predictor)) {
        if (!F.has(n)) continue;
        F.fill_ghosts(n);
    }
}

BcStep::BcStep(const Level0& l0, Fields& F) : m_l0(l0), m_F(F) {}

bool BcStep::local(int i) const { return m_l0.dm[i] == amrex::ParallelDescriptor::MyProc(); }

namespace {
// OMESH array index of FDS_G_FILL_OM
const std::vector<std::pair<const char*, int>> kOm = {{"MU", 1}, {"RHO", 2}, {"RHOS", 3}, {"U", 4}, {"V", 5}, {"W", 6}, {"US", 7}, {"VS", 8}, {"WS", 9}, {"H", 10}, {"HS", 11},
                                                        {"FVX", 12}, {"FVY", 13}, {"FVZ", 14}, {"D", 15}, {"DS", 16}, {"KRES", 17}, {"ZZ", 19}, {"ZZS", 20}};
}

// OMESH(NOM) of every local box NM from the native FAB of every box NOM (the owner broadcasts its FAB when NOM is on another rank)
void BcStep::fill_omesh()
{
    const int nbox = static_cast<int>(m_l0.ba.size());
    for (const auto& kv : kOm) {
        if (!m_F.has(kv.first)) continue;
        amrex::MultiFab& mf = m_F[kv.first];
        const FieldSpec& sp = m_F.spec(kv.first);
        for (int nom = 0; nom < nbox; ++nom) {
            const FdsBounds nb = fds_bounds(sp, m_l0.ba[nom]);
            const int nc = mf.nComp();
            const std::size_t n = static_cast<std::size_t>(nb.ext[0]) * nb.ext[1] * nb.ext[2] * nc;
            std::vector<double> buf(n);
            const int owner = m_l0.dm[nom];
            if (local(nom)) {
                const amrex::FArrayBox& fab = mf[nom];
                std::copy(fab.dataPtr(), fab.dataPtr() + n, buf.begin());
            }
            if (amrex::ParallelDescriptor::NProcs() > 1) amrex::ParallelDescriptor::Bcast(buf.data(), n, owner);
            for (int nm = 0; nm < nbox; ++nm)
                if (local(nm)) fds_g_fill_om(nm + 1, nom + 1, kv.second, nb.lb, nb.ext, nc, buf.data());
        }
    }
}

void BcStep::after_exchange(int code, double t, double dt)
{
    (void)dt;
    const int nbox = static_cast<int>(m_l0.ba.size());
    if (code != 1 && code != 3 && code != 4 && code != 6) return;
    fill_omesh();
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        if (code == 3 || code == 6) {
            fds_g_phase(code == 3 ? 1 : 0);
            if (ext_ghost) fds_p_save_uvw(nm + 1, code == 3 ? 1 : 0); else fds_g_match(nm + 1);
        }
    }
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        if (code == 3 || code == 6) { fds_g_phase(code == 3 ? 1 : 0); if (iface_hook) iface_hook(true); fds_g_velocity_bc(t, nm + 1, code == 3 ? 1 : 0); if (iface_hook) iface_hook(false); }
        else fds_g_viscosity_bc(nm + 1, code == 4 ? 1 : 0);
    }
}

void BcStep::wall_bc(int predictor, double t, double dt)
{
    fds_g_phase(predictor);
    for (int nm = 0; nm < static_cast<int>(m_l0.ba.size()); ++nm)
        if (local(nm)) fds_g_wall_bc(t, dt, nm + 1);
}

void BcStep::replay_velocity(double t, double dt, bool predictor)
{
    (void)dt;
    // A frozen dump state is a snapshot: the ghost values of the estimated (US,VS,WS) and final (U,V,W) velocities were written at different points of the
    // step, from interiors that changed in between (DENSITY restores the boundary-face values, mass.f90:426/598). So the replay regenerates the ghosts that the
    // kernels of THIS stage read: predictor records U,V,W (VELOCITY_BC after the previous corrector), corrector records US,VS,WS (after the predictor).
    const int nbox = static_cast<int>(m_l0.ba.size());
    fill_omesh();
    fds_g_phase(predictor ? 0 : 1);   // FDS calls VELOCITY_BC(final velocities) with CORRECTOR set (main.f90 1174) and VELOCITY_BC(estimated) with PREDICTOR set (969)
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        fds_g_viscosity_bc(nm + 1, predictor ? 0 : 1);
        fds_g_mu_edges(nm + 1);   // the clamped edge-cell copies of MU, KRES that COMPUTE_VISCOSITY ends with (a frozen state has them from the periodic fill)
        fds_g_velocity_bc(t, nm + 1, predictor ? 0 : 1);
    }
    fds_g_phase(predictor ? 1 : 0);
}

void BcStep::replay(double t, double dt, bool predictor)
{
    replay_velocity(t, dt, predictor);
    // the density/species/temperature ghosts that a kernel of this record reads were set by the WALL_BC of the previous stage: CORRECTOR for the predictor kernels
    // (they read RHO, ZZ, TMP), PREDICTOR for the corrector kernels (RHOS, ZZS)
    wall_bc(predictor ? 0 : 1, t, dt);
    fds_g_phase(predictor ? 1 : 0);
}

}  // namespace fdsamr
