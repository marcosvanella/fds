// GhostExchange.cpp: see GhostExchange.H. Kernel-facing rules (M2a): (a) passive scalars are extra components of ZZ/ZZS (Fields.cpp), all
// components move; (b) only uniform Cartesian metrics are used, R(I)/RRN(I) are dropped.
#include "GhostExchange.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cstdlib>
#include <map>
#include <string>

extern "C" {
int fds_g_fill_om(int nm, int nom, int which, const int* lb, const int* ext, int nc, const double* data);
void fds_g_phase(int pred);
void fds_g_match(int nm);
void fds_p_save_uvw(int nm, int pred);
void fds_g_wall_bc(double t, double dt, int nm);
void fds_g_velocity_bc(double t, int nm, int est);
void fds_g_viscosity_bc(int nm, int est);
void fds_g_mu_edges(int nm);
void fds_g_mu_edges_dom(int nm, int mask, int which);
}

namespace {
// Fault injection for the regression tests (tests/run_periodic_regression.sh): FDSTL_SKIP_FIX=match,mu,kres (any subset, comma separated) switches the named periodic-case fix
// off again. Never set in production runs; the regression test sets each in turn and requires the csmag_32 periodic comparison to FAIL, so a fix cannot silently disappear.
bool skip_fix(const char* name)
{
    const char* e = std::getenv("FDSTL_SKIP_FIX");
    if (!e) return false;
    const std::string s = std::string(",") + e + ",";
    return s.find(std::string(",") + name + ",") != std::string::npos;
}
}  // namespace

namespace fdsamr {

void ghost_exchange(Fields& F, int code, bool predictor)
{
    for (const auto& n : exchange_fields(code, predictor)) {
        if (!F.has(n)) continue;
        F.fill_ghosts(n);
    }
}

BcStep::BcStep(const Level0& l0, Fields& F) : m_l0(l0), m_F(F) {}

bool BcStep::local(int i) const { return m_l0.dm[i] == amrex::ParallelDescriptor::MyProc(); }

extern double g_prof[4];   // Fields.cpp
namespace {
struct ProfScope { int k; double t0; explicit ProfScope(int kk) : k(kk), t0(amrex::second()) {} ~ProfScope() { g_prof[k] += amrex::second() - t0; } };
}

namespace {
// OMESH array index of FDS_G_FILL_OM
const std::vector<std::pair<const char*, int>> kOm = {{"MU", 1}, {"RHO", 2}, {"RHOS", 3}, {"U", 4}, {"V", 5}, {"W", 6}, {"US", 7}, {"VS", 8}, {"WS", 9}, {"H", 10}, {"HS", 11},
                                                        {"FVX", 12}, {"FVY", 13}, {"FVZ", 14}, {"D", 15}, {"DS", 16}, {"KRES", 17}, {"ZZ", 19}, {"ZZS", 20}};
}

// OMESH(NOM) of every local box NM from the native FAB of every box NOM (the owner broadcasts its FAB when NOM is on another rank)
void BcStep::fill_omesh()
{
    const double tp0 = amrex::second();
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
                if (local(nm)) fds_g_fill_om(nm + 1 + m_l0.fds_mesh_offset, nom + 1 + m_l0.fds_mesh_offset, kv.second, nb.lb, nb.ext, nc, buf.data());
        }
    }
    g_prof[1] += amrex::second() - tp0;
}

void BcStep::exchange(int code, bool predictor)
{
    if (m_l0.level > 0 && !ext_ghost) amrex::Abort("BcStep: a level > 0 has coarse-fine faces and needs EXTERNAL_GHOSTS_FILLED (ext_ghost, patches 0003/0004); the OMESH route is same-level only");
    ghost_exchange(m_F, code, predictor);
    if (cf_ghost_hook) {
        CfGhostRequest r{m_l0.level, code, predictor, {}};
        for (const auto& n : exchange_fields(code, predictor)) if (m_F.has(n)) r.fields.push_back(n);
        cf_ghost_hook(r);
    }
}

void BcStep::after_exchange(int code, double t, double dt)
{
    (void)dt;
    const int nbox = static_cast<int>(m_l0.ba.size());
    if (code != 1 && code != 3 && code != 4 && code != 6) return;
    if (ext_ghost && (code == 3 || code == 6))   // UVW_SAVE: the face velocities before the match (DENSITY restores them at the wall faces)
        for (int nm = 0; nm < nbox; ++nm)
            if (local(nm)) { fds_g_phase(code == 3 ? 1 : 0); fds_p_save_uvw(nm + 1 + m_l0.fds_mesh_offset, code == 3 ? 1 : 0); }
    if (ext_ghost && (code == 3 || code == 6)) {
        // the periodic domain faces: FDS's MATCH_VELOCITY averages the two copies of the flow face (patch 0003 skips it, the driver does it on the AMReX data)
        static const char* const pn[2][3] = {{"U", "V", "W"}, {"US", "VS", "WS"}};
        for (int d = 0; d < 3; ++d)
            if (m_l0.dom.periodic[d] && m_F.has(pn[code == 3][d]) && !skip_fix("match")) match_periodic_faces(pn[code == 3][d], d);
        // the ghost layer next to a matched face holds the matched value in FDS (the wall-cell loop of MATCH_VELOCITY runs before VELOCITY_BC reads it): refill
        for (int d = 0; d < 3; ++d)
            if (m_l0.dom.periodic[d] && m_F.has(pn[code == 3][d])) m_F.fill_ghosts(pn[code == 3][d]);
    }
    fill_omesh();   // after the match: FDS's MATCH_VELOCITY also writes the averaged values into OMESH, which VELOCITY_BC reads for the periodic ghosts
    ProfScope prof_bc(3);
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        if (code == 3 || code == 6) {
            fds_g_phase(code == 3 ? 1 : 0);
            if (!ext_ghost) fds_g_match(nm + 1 + m_l0.fds_mesh_offset);
        }
    }
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        if (code == 3 || code == 6) { fds_g_phase(code == 3 ? 1 : 0); if (iface_hook) iface_hook(true); fds_g_velocity_bc(t, nm + 1 + m_l0.fds_mesh_offset, code == 3 ? 1 : 0); if (iface_hook) iface_hook(false); }
        else {
            fds_g_viscosity_bc(nm + 1 + m_l0.fds_mesh_offset, code == 4 ? 1 : 0);
            // FDS ends COMPUTE_VISCOSITY with clamped copies of MU, KRES in the edge cells of the domain; the full ghost fill above replaced them by periodic images
            const amrex::Box& b = m_l0.ba[nm];
            const amrex::Box& d = m_l0.geom.Domain();
            int mask = 0;
            for (int dir = 0; dir < 3; ++dir) {
                if (b.smallEnd(dir) == d.smallEnd(dir)) mask |= 1 << (2 * dir);
                if (b.bigEnd(dir) == d.bigEnd(dir)) mask |= 2 << (2 * dir);
            }
            const int which = (skip_fix("mu") ? 0 : 1) | (skip_fix("kres") ? 0 : 2);   // bit 0: MU, bit 1: KRES
            if (which) fds_g_mu_edges_dom(nm + 1 + m_l0.fds_mesh_offset, mask, which);
        }
    }
}

void BcStep::match_periodic_faces(const std::string& name, int dir)
{
    const double tp0 = amrex::second();
    amrex::MultiFab& mf = m_F[name];
    const amrex::Box dom = m_l0.geom.Domain();
    const int n = dom.length(dir);
    amrex::Box p0 = amrex::surroundingNodes(dom, dir);
    p0.setBig(dir, p0.smallEnd(dir));
    amrex::BoxArray b0(p0);
    b0.maxSize(32);
    amrex::IntVect sh(0); sh[dir] = n;
    amrex::BoxArray bn(b0);
    bn.shift(sh);
    amrex::DistributionMapping dm(b0);
    amrex::MultiFab P0(b0, dm, 1, 0), PN(bn, dm, 1, 0);
    P0.setVal(0.0); PN.setVal(0.0);
    P0.ParallelCopy(mf, 0, 0, 1, 0, 0);
    PN.ParallelCopy(mf, 0, 0, 1, 0, 0);
    // FDS: DA_OTHER = 0 + A1*A2; UU_OTHER = 0 + ((OM*A1)*A2)/DA_OTHER; UU = 0.5*(UU + UU_OTHER). A1, A2: the two other cell sizes (DY,DZ / DX,DZ / DX,DY)
    const double a1 = dir == 0 ? m_l0.dx[1] : m_l0.dx[0], a2 = dir == 2 ? m_l0.dx[1] : m_l0.dx[2];
    const double da = 0.0 + a1 * a2;
    for (amrex::MFIter mfi(P0); mfi.isValid(); ++mfi) {
        auto lo = P0.array(mfi);
        auto hi = PN.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const int ih = i + (dir == 0 ? n : 0), jh = j + (dir == 1 ? n : 0), kh = k + (dir == 2 ? n : 0);
            const double a = lo(i, j, k), b = hi(ih, jh, kh);
            const double other_hi = 0.0 + ((b * a1) * a2) / da;   // what the low face sees of the high face
            const double other_lo = 0.0 + ((a * a1) * a2) / da;   // and the other way round
            lo(i, j, k) = 0.5 * (a + other_hi);
            hi(ih, jh, kh) = 0.5 * (b + other_lo);
        });
    }
    mf.ParallelCopy(P0, 0, 0, 1, 0, 0);
    mf.ParallelCopy(PN, 0, 0, 1, 0, 0);
    g_prof[2] += amrex::second() - tp0;
}

void BcStep::wall_bc(int predictor, double t, double dt)
{
    fds_g_phase(predictor);
    for (int nm = 0; nm < static_cast<int>(m_l0.ba.size()); ++nm)
        if (local(nm)) fds_g_wall_bc(t, dt, nm + 1 + m_l0.fds_mesh_offset);
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
        fds_g_viscosity_bc(nm + 1 + m_l0.fds_mesh_offset, predictor ? 0 : 1);
        fds_g_mu_edges(nm + 1 + m_l0.fds_mesh_offset);   // the clamped edge-cell copies of MU, KRES that COMPUTE_VISCOSITY ends with (a frozen state has them from the periodic fill)
        fds_g_velocity_bc(t, nm + 1 + m_l0.fds_mesh_offset, predictor ? 0 : 1);
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
