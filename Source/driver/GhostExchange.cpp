// GhostExchange.cpp: see GhostExchange.H. Kernel-facing rules (M2a): (a) passive scalars are extra components of ZZ/ZZS (Fields.cpp), all
// components move; (b) only uniform Cartesian metrics are used, R(I)/RRN(I) are dropped.
#include "GhostExchange.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cstdio>
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
void fds_wall_ghost(double u, double rho, double tmp, double muc, double dz, double rough, double* g);
void fds_ambient_state(int k, int nz, double* rho, double* tmp, double* zz);
int fds_n_tracked_species();
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

void mirror_domain_edges(const Level& lev, amrex::MultiFab& mf, int nodal_dir, const FaceBc& bc)
{
    const amrex::Box dom = amrex::convert(lev.geom.Domain(), mf.ixType());   // nodal in nodal_dir: its ends are the wall faces
    const int nc = mf.nComp();
    for (int d = 0; d < 3; ++d) {
        if (lev.geom.isPeriodic(d)) continue;
        const bool normal = d == nodal_dir;
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox(), fb = mf[mfi].box();
            auto a = mf.array(mfi);
            for (int side = 0; side < 2; ++side) {
                const int edge = side == 0 ? dom.smallEnd(d) : dom.bigEnd(d), f = bc[2 * d + side];
                if ((side == 0 ? vb.smallEnd(d) : vb.bigEnd(d)) != edge || (normal && f != kFaceInert && f != kFaceMirror)) continue;
                const double s = normal ? 0.0 : f == kFaceInert ? -1.0 : 1.0;   // normal: 0 on the wall face itself (g = 0 only); tangential: odd (no-slip) or even copy
                const int g1 = normal ? 0 : mf.nGrow(d);
                for (int g = normal ? 0 : 1; g <= g1; ++g) {
                    amrex::Box slab = fb;
                    slab.setSmall(d, side == 0 ? edge - g : edge + g); slab.setBig(d, slab.smallEnd(d));
                    amrex::IntVect sh(0, 0, 0);
                    sh[d] = side == 0 ? g : -g;
                    amrex::LoopOnCpu(slab, nc, [&](int i, int j, int k, int n) { a(i, j, k, n) = s * a(i + sh[0], j + sh[1], k + sh[2], n); });
                }
            }
        }
    }
}

BcStep::BcStep(const Level0& l0, Fields& F) : m_l0(l0), m_mesh_offset(l0.fds_mesh_offset), m_F(F) {}

bool BcStep::local(int i) const { return m_l0.dm[i] == amrex::ParallelDescriptor::MyProc(); }


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
                if (local(nm)) fds_g_fill_om(nm + 1 + m_mesh_offset, nom + 1 + m_mesh_offset, kv.second, nb.lb, nb.ext, nc, buf.data());
        }
    }
    g_prof[1] += amrex::second() - tp0;
}

void BcStep::exchange(int code, bool predictor)
{
    if (m_l0.level > 0 && !ext_ghost) amrex::Abort("BcStep: a level > 0 has coarse-fine faces and needs EXTERNAL_GHOSTS_FILLED (ext_ghost, patches 0003/0004); the OMESH route is same-level only");
    ghost_exchange(m_F, code, predictor);
    if (cf_ghost_hook) {
        ProfScope prof_cf(5);
        CfGhostRequest r{m_l0.level, code, predictor, {}};
        for (const auto& n : exchange_fields(code, predictor)) if (m_F.has(n)) r.fields.push_back(n);
        cf_ghost_hook(r);
    }
    if (m_l0.level > 0) {   // domain-edge faces of a fine box: mirror ghost layers (no FDS wall cell there)
        for (const auto& n : exchange_fields(code, predictor)) {
            if (!m_F.has(n)) continue;
            const FieldSpec& sp = m_F.spec(n);
            int nd = -1;
            for (int d = 0; d < 3; ++d) if (sp.nodal(d)) nd = d;
            const bool vel = n == "U" || n == "V" || n == "W" || n == "US" || n == "VS" || n == "WS";
            mirror_domain_edges(m_l0, m_F[n], nd, vel ? vel_bc : FaceBc{});
            if (!vel) inflow_ambient_ghosts(n, predictor);
            if (vel && (n == "U" || n == "V" || n == "US" || n == "VS")) ground_wall_ghosts(n);
            if (vel && nd >= 0) zero_mixed_wall(m_F[n], nd, false);
        }
    }
}

// An OPEN domain face of a fine box where the gas flows in: the scalar ghost layers hold the ambient state (as the level-0 inflow wall does), not the mirror of the cell inside (fire gas
// there would re-enter through the face and the fuel in it would be counted as inflow).
void BcStep::inflow_ambient_ghosts(const std::string& n, bool predictor)
{
    const bool s = n == "RHOS" || n == "ZZS";
    if (!(s || n == "RHO" || n == "ZZ" || n == "TMP")) return;
    const int nz = fds_n_tracked_species();
    amrex::MultiFab& mf = m_F[n];
    const amrex::Box dom = m_l0.geom.Domain();
    const char* un[3] = {s ? "US" : "U", s ? "VS" : "V", s ? "WS" : "W"};
    (void)predictor;
    for (int d = 0; d < 3; ++d) {
        if (m_l0.geom.isPeriodic(d)) continue;
        for (int side = 0; side < 2; ++side) {
            const int f = vel_bc[2 * d + side];
            if (f == kFaceInert || f == kFaceMirror) continue;
            const int edge = side == 0 ? dom.smallEnd(d) : dom.bigEnd(d) + 1;   // the wall face of the nodal velocity
            const int cell = side == 0 ? dom.smallEnd(d) : dom.bigEnd(d);
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const amrex::Box vb = mfi.validbox();
                if ((side == 0 ? vb.smallEnd(d) : vb.bigEnd(d)) != cell) continue;
                const auto u = m_F[un[d]].const_array(mfi);
                auto a = mf.array(mfi);
                for (int g = 1; g <= mf.nGrow(d); ++g) {
                    amrex::Box slab = mf[mfi].box();
                    slab.setRange(d, side == 0 ? cell - g : cell + g);
                    amrex::LoopOnCpu(slab, [&](int i, int j, int k) {
                        int c[3] = {i, j, k}; c[d] = edge;
                        const double vn = u(c[0], c[1], c[2]);
                        if (!(side == 0 ? vn > 0.0 : vn < 0.0)) return;
                        double rho, tmp, zz[16];
                        fds_ambient_state(amrex::coarsen(k, m_l0.geom.Domain().length(2) / l0_len[2]) + 1, nz, &rho, &tmp, zz);
                        if (n == "RHO" || n == "RHOS") a(i, j, k) = rho;
                        else if (n == "TMP") a(i, j, k) = tmp;
                        else for (int q = 0; q < nz; ++q) a(i, j, k, q) = zz[q];
                    });
                }
            }
        }
    }
}

// A fine box on the ground (INERT floor): the tangential velocity ghost layer below the first cell is the wall-model stress of FDS written as a ghost value (fds_wall_ghost) instead of
// the odd (no-slip) copy, which with the eddy viscosity of the cell is a far larger drag than the log-law stress the level-0 ground has.
void BcStep::ground_wall_ghosts(const std::string& n)
{
    if (vel_bc[4] != kFaceInert || std::getenv("FDSTL_NOWALLMODEL")) return;
    const amrex::Box dom = m_l0.geom.Domain();
    const bool s = n.back() == 'S';
    const double dz = m_l0.geom.CellSize(2);
    const int rx = dom.length(0) / l0_len[0], ry = dom.length(1) / l0_len[1];
    amrex::MultiFab& mf = m_F[n];
    // The face between two fine boxes belongs to both (nodal): both must compute the same ghost for it, so the cell states they read are first made the same in the ghost cells of both
    // (a stale copy gave the two boxes different shear at that face: composite mass off by ~1e-4 once a regrid had made fine-fine box faces next to a fire).
    for (const char* nm : {"MU", "TMP", s ? "RHOS" : "RHO"}) m_F[nm].FillBoundary(m_l0.geom.periodicity());
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        if (mfi.validbox().smallEnd(2) != dom.smallEnd(2)) continue;
        auto a = mf.array(mfi);
        const auto rho = m_F[s ? "RHOS" : "RHO"].const_array(mfi);
        const auto tmp = m_F["TMP"].const_array(mfi);
        const auto mu = m_F["MU"].const_array(mfi);
        amrex::Box slab = mf[mfi].box();
        slab.setRange(2, dom.smallEnd(2) - 1);
        const amrex::Box cells = mfi.validbox();
        amrex::LoopOnCpu(slab, [&](int i, int j, int k) {
            const int fu = n.front() == 'U' ? 1 : 0, fv = n.front() == 'V' ? 1 : 0;   // the face index of the nodal direction: the low face of the first cell is smallEnd - 1
            if (i < cells.smallEnd(0) - fu || i > cells.bigEnd(0) || j < cells.smallEnd(1) - fv || j > cells.bigEnd(1)) return;
            const double u = a(i, j, k + 1);
            if (u == 0.0) return;
            const int ic = std::max(i, cells.smallEnd(0)), jc = std::max(j, cells.smallEnd(1));
            double g;
            fds_wall_ghost(u, rho(ic, jc, k + 1), tmp(ic, jc, k + 1), mu(ic, jc, k + 1), dz, ground_rough[amrex::coarsen(ic, rx) + l0_len[0] * amrex::coarsen(jc, ry)], &g);
            a(i, j, k) = g;
        });
    }
}

void BcStep::zero_wall_flux()
{
    static const char* const fv[3] = {"FVX", "FVY", "FVZ"};
    const amrex::Box& dom = m_l0.geom.Domain();
    for (int d = 0; d < 3; ++d)
        for (amrex::MFIter mfi(m_F[fv[d]]); mfi.isValid(); ++mfi) {
            const amrex::Box fb = mfi.fabbox();
            auto a = m_F[fv[d]].array(mfi);
            for (int side = 0; side < 2; ++side) {
                const int f = vel_bc[2 * d + side];
                if (f != kFaceInert && f != kFaceMirror) continue;
                amrex::Box s = fb;
                s.setRange(d, side == 0 ? dom.smallEnd(d) - 1 : dom.bigEnd(d));
                s &= fb;
                if (s.ok()) amrex::LoopOnCpu(s, [&](int i, int j, int k) { a(i, j, k) = 0.0; });
            }
        }
    for (int d = 0; d < 3; ++d) zero_mixed_wall(m_F[fv[d]], d, true);   // mixed face: FV = 0 over the SOLID cells
}

void BcStep::zero_mixed_wall(amrex::MultiFab& mf, int d, bool cell_shaped)
{
    const amrex::Box dom = m_l0.geom.Domain();
    const int t1 = d == 0 ? 1 : 0, t2 = d == 2 ? 1 : 2;
    const int r1 = dom.length(t1) / l0_len[t1], r2 = dom.length(t2) / l0_len[t2];
    for (int side = 0; side < 2; ++side) {
        const std::vector<int>& w = face_wtype[2 * d + side];
        if (w.empty()) continue;
        const int plane = cell_shaped ? (side == 0 ? dom.smallEnd(d) - 1 : dom.bigEnd(d)) : (side == 0 ? dom.smallEnd(d) : dom.bigEnd(d) + 1);
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            amrex::Box s = mf[mfi].box();
            if (plane < s.smallEnd(d) || plane > s.bigEnd(d)) continue;
            s.setRange(d, plane);
            auto a = mf.array(mfi);
            amrex::LoopOnCpu(s, [&](int i, int j, int k) {
                const int c[3] = {i, j, k};
                const int a1 = amrex::coarsen(c[t1] - dom.smallEnd(t1), r1), a2 = amrex::coarsen(c[t2] - dom.smallEnd(t2), r2);
                if (a1 < 0 || a2 < 0 || a1 >= l0_len[t1] || a2 >= l0_len[t2]) return;
                if (w[a1 + l0_len[t1] * a2] == 1) a(i, j, k) = 0.0;   // SOLID_BOUNDARY
            });
        }
    }
}

void BcStep::after_exchange(int code, double t, double dt)
{
    (void)dt;
    const int nbox = static_cast<int>(m_l0.ba.size());
    if (code != 1 && code != 3 && code != 4 && code != 6) return;
    // A level > 0 has no external wall cell: its domain-edge ghosts are the mirror rule (mirror_domain_edges), its coarse-fine ghosts the composite fill. FDSTL_SKIPAFT (debug) skips steps on a
    // level > 0: bit 1 save_uvw, 4 velocity_bc (bit 2 was fill_omesh: never run on a fine level, a fine box has no OMESH). Default 0 (S14.4). Until S14.3 the default was 6 because VELOCITY_BC of a fine box overwrote the periodic ghosts of LEVEL 0:
    // the cause was the mesh number (this BcStep called the FDS routines with box index + 1 instead of FINE_LEVEL NM0 + box index + 1, so they ran on level-0 mesh objects); see notes/fine-velocity-bc.md.
    static const int dbg_skip = std::getenv("FDSTL_SKIPAFT") ? std::atoi(std::getenv("FDSTL_SKIPAFT")) : 0;
    const bool fine_dbg = m_l0.level > 0;
    if (ext_ghost && (code == 3 || code == 6) && !(fine_dbg && (dbg_skip & 1)))   // UVW_SAVE: the face velocities before the match (DENSITY restores them at the wall faces)
        for (int nm = 0; nm < nbox; ++nm)
            if (local(nm)) { if (std::getenv("FDSTL_SAVEDBG")) std::fprintf(stderr, "SAVEDBG rank %d level %d box %d of %d fds_mesh %d offset %d\n", amrex::ParallelDescriptor::MyProc(), m_l0.level, nm, nbox, nm + 1 + m_mesh_offset, m_mesh_offset); fds_g_phase(code == 3 ? 1 : 0); fds_p_save_uvw(nm + 1 + m_mesh_offset, code == 3 ? 1 : 0); }
    if (ext_ghost && (code == 3 || code == 6)) {
        // the periodic domain faces: FDS's MATCH_VELOCITY averages the two copies of the flow face (patch 0003 skips it, the driver does it on the AMReX data)
        static const char* const pn[2][3] = {{"U", "V", "W"}, {"US", "VS", "WS"}};
        for (int d = 0; d < 3; ++d)
            if (m_l0.dom.periodic[d] && m_F.has(pn[code == 3][d]) && !skip_fix("match")) match_periodic_faces(pn[code == 3][d], d);
        // the ghost layer next to a matched face holds the matched value in FDS (the wall-cell loop of MATCH_VELOCITY runs before VELOCITY_BC reads it): refill
        for (int d = 0; d < 3; ++d)
            if (m_l0.dom.periodic[d] && m_F.has(pn[code == 3][d])) m_F.fill_ghosts(pn[code == 3][d]);
    }
    if (!fine_dbg) fill_omesh();   // a fine box has no OMESH (no neighbour mesh objects): its same-level ghosts come from the AMReX fill, MESHES(NM) must not be indexed with a fine number   // after the match: FDS's MATCH_VELOCITY also writes the averaged values into OMESH, which VELOCITY_BC reads for the periodic ghosts
    ProfScope prof_bc(3);
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        if (code == 3 || code == 6) {
            fds_g_phase(code == 3 ? 1 : 0);
            if (!ext_ghost) fds_g_match(nm + 1 + m_mesh_offset);
        }
    }
    for (int nm = 0; nm < nbox; ++nm) {
        if (!local(nm)) continue;
        if (code == 3 || code == 6) { fds_g_phase(code == 3 ? 1 : 0); if (iface_hook) iface_hook(true); if (!(fine_dbg && (dbg_skip & 4))) { if (std::getenv("FDSTL_SAVEDBG")) std::fprintf(stderr, "SAVEDBG velocity_bc rank %d level %d fds_mesh %d\n", amrex::ParallelDescriptor::MyProc(), m_l0.level, nm + 1 + m_mesh_offset); fds_g_velocity_bc(t, nm + 1 + m_mesh_offset, code == 3 ? 1 : 0); } if (iface_hook) iface_hook(false); }
        else {
            fds_g_viscosity_bc(nm + 1 + m_mesh_offset, code == 4 ? 1 : 0);
            // FDS ends COMPUTE_VISCOSITY with clamped copies of MU, KRES in the edge cells of the domain; the full ghost fill above replaced them by periodic images
            const amrex::Box& b = m_l0.ba[nm];
            const amrex::Box& d = m_l0.geom.Domain();
            int mask = 0;
            for (int dir = 0; dir < 3; ++dir) {
                if (b.smallEnd(dir) == d.smallEnd(dir)) mask |= 1 << (2 * dir);
                if (b.bigEnd(dir) == d.bigEnd(dir)) mask |= 2 << (2 * dir);
            }
            const int which = (skip_fix("mu") ? 0 : 1) | (skip_fix("kres") ? 0 : 2);   // bit 0: MU, bit 1: KRES
            if (which) fds_g_mu_edges_dom(nm + 1 + m_mesh_offset, mask, which);
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
    amrex::MultiFab P0(b0, dm, 1, 0), PN(bn, dm, 1, 0), M0(b0, dm, 1, 0), MN(bn, dm, 1, 0), one(mf.boxArray(), mf.DistributionMap(), 1, 0);
    P0.setVal(0.0); PN.setVal(0.0); M0.setVal(0.0); MN.setVal(0.0); one.setVal(1.0);
    P0.ParallelCopy(mf, 0, 0, 1, 0, 0);
    PN.ParallelCopy(mf, 0, 0, 1, 0, 0);
    M0.ParallelCopy(one, 0, 0, 1, 0, 0);   // 1 where the level has that copy of the face (a fine level may have only one)
    MN.ParallelCopy(one, 0, 0, 1, 0, 0);
    // FDS: DA_OTHER = 0 + A1*A2; UU_OTHER = 0 + ((OM*A1)*A2)/DA_OTHER; UU = 0.5*(UU + UU_OTHER). A1, A2: the two other cell sizes (DY,DZ / DX,DZ / DX,DY)
    const double a1 = dir == 0 ? m_l0.dx[1] : m_l0.dx[0], a2 = dir == 2 ? m_l0.dx[1] : m_l0.dx[2];
    const double da = 0.0 + a1 * a2;
    for (amrex::MFIter mfi(P0); mfi.isValid(); ++mfi) {
        auto lo = P0.array(mfi);
        auto hi = PN.array(mfi);
        auto m0 = M0.const_array(mfi); auto mn = MN.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const int ih = i + (dir == 0 ? n : 0), jh = j + (dir == 1 ? n : 0), kh = k + (dir == 2 ? n : 0);
            if (m0(i, j, k) == 0.0 || mn(ih, jh, kh) == 0.0) return;   // a fine level may hold one copy only: nothing to match
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
        if (local(nm)) fds_g_wall_bc(t, dt, nm + 1 + m_mesh_offset);
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
        fds_g_viscosity_bc(nm + 1 + m_mesh_offset, predictor ? 0 : 1);
        fds_g_mu_edges(nm + 1 + m_mesh_offset);   // the clamped edge-cell copies of MU, KRES that COMPUTE_VISCOSITY ends with (a frozen state has them from the periodic fill)
        fds_g_velocity_bc(t, nm + 1 + m_mesh_offset, predictor ? 0 : 1);
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
