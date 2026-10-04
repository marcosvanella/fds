// test_blob_registry.cpp (R4 part 2): end-to-end dynamic regrid of a moving blob through the REAL pieces of the data path: Role 3's RegridAmrCore (tags -> grids),
// Role 1's LevelRegistry (per-level Fields) and Role 1's RegistryTransfer (the LevelDataTransfer on the registry: conservative prolongation of rho*Z, D-060 face
// prolongation, bitwise copy of old fine data, conservative restriction). The FDS kernels are NOT involved (no fine-level FDS binding yet); the physics is a small
// conservative upwind transport of rho*Z_n (n = 1..3) with a prescribed analytic velocity that is constant in time (uniform translation plus a periodic field), one global
// time step on all levels, and the interface flux overwrite (coarse flux on a covered face = mean of the fine fluxes, no refluxing).
//
// Runs (all periodic, species n=1 is the blob, n=3 is ~1e-9 away from it so that positivity is stressed):
//   3-D  32^3 coarse, 2 levels ratio 2, two level-0 boxes; 2-D 48x1x48 coarse, 3 levels ratio 2; 2-D 32x1x32 coarse, 2 levels ratio 4.
// Checks per run: composite mass and each species change per regrid <= 1e-12 relative, and over the whole run; no negative rho*Z and zero clips of the prolongation;
// every cell above the tag threshold at the next regrid lies on the finest level; max|div u - D| (D = 0 for the prescribed field) reported per level after each regrid;
// D-059 shared ghost cell count printed per regrid; the AMR result against the uniform-fine run (and the uniform-coarse run for reference).
// Corner check: in the first steps a static two-level hierarchy equals the uniform-fine run BITWISE in every fine cell farther than n+1 cells from the patch outline (the
// dependence distance of an n-step upwind scheme), corners of the patch included.
// Negative controls: no flux overwrite, linear Z average-down, no tag buffer.
// Prints HASH (grids + data) for the 1/2/4-rank comparison done by run_regrid_rank_check.sh.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "DriverAdapter.H"
#include "ExactSum.H"
#include "FaceTransfer.H"
#include "GhostShare.H"
#include "Hierarchy.H"
#include "LevelOps.H"
#include "MlmgCompositeSolver.H"
#include "PostRegridProjection.H"
#include "RegistryTransfer.H"
#include "RegridAmrCore.H"
#include "TagOps.H"
#include "check.H"

using namespace fdsrt;

namespace {

constexpr int NS = 3;
const double PI = 3.14159265358979323846;

struct Cfg {
    std::string name;
    int dim = 3;                 // 3 or 2 (hidden direction y)
    int n0 = 32;                 // coarse cells per non-hidden direction
    int max_level = 1;
    int ratio = 2;
    int nbuf = 2;
    int interval = 4;            // steps between regrids; 0 = static hierarchy
    int steps = 40;
    bool overwrite = true;       // interface flux overwrite
    bool linearZ = false;        // negative control: restriction by linear averaging of Z
    double speed = 3.0;
    bool null_solver = false;    // negative control of the acceptance test: a solver that reports success and returns zero gradients
    bool project = false;        // D-063 post-regrid composite projection (mock MLMG solver); false: the divergence is only measured
};

struct Blob { double xc, yc, zc; };

// state at a point: rho, mass fractions (sum 1)
void cell_state(const Cfg& c, double x, double y, double z, double* rho, double* zz)
{
    const double r2 = (x - 0.30) * (x - 0.30) + (c.dim == 3 ? (y - 0.40) * (y - 0.40) : 0.0) + (z - 0.50) * (z - 0.50);
    const double s = std::exp(-r2 / (2.0 * 0.06 * 0.06));
    zz[0] = 0.05 + 0.9 * s;
    const double w = std::exp(-r2 / (2.0 * 0.03 * 0.03));              // ~1e-9 far from the centre
    zz[2] = (1.0 - zz[0]) * 0.9 * w;
    zz[1] = 1.0 - zz[0] - zz[2];
    *rho = 1.2 / (1.0 + 3.0 * s);
}

double vel_comp(const Cfg& c, int dir, double x, double y, double z)
{
    const double U0[3] = {c.speed, c.dim == 3 ? 0.6 * c.speed : 0.0, 0.4 * c.speed};
    const double a = 0.15 * c.speed;
    double u = U0[dir];
    if (c.dim == 3) {   // curl of A = a/(2pi) (sin(2pi y) sin(2pi z), sin(2pi z) sin(2pi x), sin(2pi x) sin(2pi y)), divergence free analytically
        if (dir == 0) u += a * std::sin(2 * PI * x) * (std::cos(2 * PI * y) - std::cos(2 * PI * z));
        if (dir == 1) u += a * std::sin(2 * PI * y) * (std::cos(2 * PI * z) - std::cos(2 * PI * x));
        if (dir == 2) u += a * std::sin(2 * PI * z) * (std::cos(2 * PI * x) - std::cos(2 * PI * y));
    } else {            // stream function psi = a/(2pi) sin(2pi x) sin(2pi z): u_x = -dpsi/dz, u_z = dpsi/dx
        if (dir == 0) u += -a * std::sin(2 * PI * x) * std::cos(2 * PI * z);
        if (dir == 2) u += a * std::cos(2 * PI * x) * std::sin(2 * PI * z);
    }
    return u;
}

// Face value of the prescribed velocity on a grid of cell size dx: the DISCRETE curl of the analytic vector potential sampled on the cell edges (A_y only in 2-D), plus the uniform
// translation. Divergence-free to round-off on every level (D = 0 exactly), so that on a hierarchy any max|div u - D| comes from the transfer between levels, not from truncation.
double vector_potential(const Cfg& c, int comp, double x, double y, double z)
{
    const double a = 0.15 * c.speed / (2 * PI);
    if (c.dim == 3) {
        if (comp == 0) return a * std::sin(2 * PI * y) * std::sin(2 * PI * z);
        if (comp == 1) return a * std::sin(2 * PI * z) * std::sin(2 * PI * x);
        return a * std::sin(2 * PI * x) * std::sin(2 * PI * y);
    }
    return comp == 1 ? a * std::sin(2 * PI * x) * std::sin(2 * PI * z) : 0.0;
}

double vel_face(const Cfg& c, int dir, double x, double y, double z, const double* dx)
{
    const double U0[3] = {c.speed, c.dim == 3 ? 0.6 * c.speed : 0.0, 0.4 * c.speed};
    auto d = [&](int comp, int along) {   // d A_comp / d x_along, centred difference over one cell
        double p[3] = {x, y, z}, m[3] = {x, y, z};
        p[along] += 0.5 * dx[along]; m[along] -= 0.5 * dx[along];
        return (vector_potential(c, comp, p[0], p[1], p[2]) - vector_potential(c, comp, m[0], m[1], m[2])) / dx[along];
    };
    double u = U0[dir];
    if (dir == 0) u += d(2, 1) - d(1, 2);
    if (dir == 1) u += d(0, 2) - d(2, 0);
    if (dir == 2) u += d(1, 0) - d(0, 1);
    return u;
}

// Negative control for the acceptance test: claims success, corrects nothing.
class NullSolver : public fdsrt::CompositePoissonSolver {
public:
    const char* name() const override { return "null solver (control)"; }
    void rebuild(const std::vector<fdsrt::ProjectionLevel>&) override {}
    fdsrt::SolveReport solve(const std::vector<const amrex::MultiFab*>&, std::vector<std::array<amrex::MultiFab*, 3>>& grad, double, double, int) override
    {
        for (auto& a : grad) for (auto* m : a) m->setVal(0.0);
        fdsrt::SolveReport r; r.ok = true; return r;
    }
};

struct Sim {
    Cfg c;
    Hierarchy h;
    AmrParams p;
    fdsamr::Level0 l0;
    std::unique_ptr<fdsamr::Fields> F0;
    std::unique_ptr<fdsamr::SideData> sd0;
    std::unique_ptr<fdsamr::LevelRegistry> reg;
    std::unique_ptr<RegridAmrCore> core;
    std::unique_ptr<fdsamr::RegistryTransfer> rt;
    double t = 0.0, dt = 0.0;
    // diagnostics
    double worst_regrid = 0.0;                  // largest relative composite change over a regrid (rho and the species)
    std::vector<double> mass0;                  // composite at the start (rho, rho*Z_n)
    long feature_late = 0;                      // cells above the threshold on a level below max_level just before a regrid
    double worst_div = 0.0;                     // max |div u| over the levels >= 1 after a regrid (seams to retained fine data included)
    double worst_new_div = 0.0;                 // max |div u| over cells of the levels >= 1 whose six faces were all produced by the prolongation (no retained fine neighbour)
    double worst_parent = 0.0;                  // max |div(child) - div(parent)| over those cells (D-060: zero to round-off)
    long n_new_cells = 0;
    long min_neg = 0;                           // cells with a negative rho*Z_n after any step
    int n_changed = 0;
    uint64_t hier_hash = 1469598103934665603ULL;
    std::vector<std::string> log;
    // post-regrid projection (D-063): statistics over the regrids
    std::unique_ptr<fdsrt::CompositePoissonSolver> psolver;
    int proj_n = 0, proj_notaccepted = 0, proj_iter = 0, proj_failed = 0;
    std::string proj_msg;
    double proj_before = 0, proj_after = 0, proj_rhs_sum = 0;
    std::vector<double> proj_change;            // by level: largest velocity change of the projection
    std::vector<double> prof = std::vector<double>(8, 0.0);   // largest |velocity change| by distance (cells of its level) to the nearest cell with a source: 0,1,..,6, >=7
};

uint64_t mixh(uint64_t h, uint64_t v) { h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); return h; }

double total_cells_dx(const Cfg& c) { return 1.0 / (c.n0 * std::pow(c.ratio, c.max_level)); }

void fill_analytic(Sim& S, int l, bool cells, bool faces)
{
    fdsamr::Fields& F = S.reg->fields(l);
    const amrex::Geometry& g = S.core ? S.core->Geom(l) : S.l0.geom;
    const double* plo = g.ProbLo();
    const double dx[3] = {g.CellSize(0), g.CellSize(1), g.CellSize(2)};
    if (cells)
        for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
            auto rho = F["RHO"].array(mfi), rhos = F["RHOS"].array(mfi), tmp = F["TMP"].array(mfi);
            auto zz = F["ZZ"].array(mfi), zzs = F["ZZS"].array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                double r, z[NS];
                cell_state(S.c, plo[0] + (i + 0.5) * dx[0], plo[1] + (j + 0.5) * dx[1], plo[2] + (k + 0.5) * dx[2], &r, z);
                rho(i, j, k) = r; rhos(i, j, k) = r; tmp(i, j, k) = 300.0;
                for (int n = 0; n < NS; ++n) { zz(i, j, k, n) = z[n]; zzs(i, j, k, n) = z[n]; }
            });
        }
    if (faces) {
        const char* nm[3][2] = {{"U", "US"}, {"V", "VS"}, {"W", "WS"}};
        for (int d = 0; d < 3; ++d)
            for (amrex::MFIter mfi(F[nm[d][0]]); mfi.isValid(); ++mfi) {
                auto a = F[nm[d][0]].array(mfi), as = F[nm[d][1]].array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    const int idx[3] = {i, j, k};
                    double x[3];
                    for (int e = 0; e < 3; ++e) x[e] = plo[e] + (idx[e] + (e == d ? 0.0 : 0.5)) * dx[e];
                    a(i, j, k) = vel_face(S.c, d, x[0], x[1], x[2], dx);
                    as(i, j, k) = a(i, j, k);
                });
            }
    }
}

std::vector<double> composite(Sim& S)   // exact [rho, rho*Z_1..3] over the uncovered cells of all levels
{
    const int L = S.core->finestLevel();
    std::vector<amrex::MultiFab> q(L + 1);
    std::vector<const amrex::MultiFab*> mf;
    std::vector<const amrex::iMultiFab*> cov;
    std::vector<double> w;
    for (int l = 0; l <= L; ++l) {
        fdsamr::Fields& F = S.reg->fields(l);
        q[l].define(F["RHO"].boxArray(), F["RHO"].DistributionMap(), 1 + NS, 0);
        for (amrex::MFIter mfi(q[l]); mfi.isValid(); ++mfi) {
            auto o = q[l].array(mfi); auto r = F["RHO"].const_array(mfi); auto z = F["ZZ"].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { o(i, j, k, 0) = r(i, j, k); for (int n = 0; n < NS; ++n) o(i, j, k, 1 + n) = r(i, j, k) * z(i, j, k, n); });
        }
        mf.push_back(&q[l]);
        cov.push_back(S.reg->covered_mask(l));
        w.push_back(S.core->Geom(l).CellSize(0) * S.core->Geom(l).CellSize(1) * S.core->Geom(l).CellSize(2));
    }
    std::vector<double> out;
    for (int n = 0; n < 1 + NS; ++n) out.push_back(fdsamr::exact_sum_hierarchy(mf, n, w, cov));
    return out;
}

void projection_step(Sim& S);

void build(Sim& S, const Cfg& c)
{
    S.c = c;
    Report rep;
    std::string text = "&AMR MAX_LEVEL=" + std::to_string(c.max_level) + ", REF_RATIO=" + std::to_string(c.ratio) + ", BLOCKING_FACTOR=4, MAX_GRID_SIZE=16, N_ERROR_BUF=" +
                       std::to_string(c.nbuf) + ", N_PROPER=1 /\n" +
                       "&AMR_REGION XB=0.0625,0.9375," + (c.dim == 3 ? "0.0625,0.9375" : "0.0,1.0") + ",0.0625,0.9375 /\n";
    S.p = parse_amr_params(text, rep);
    const int ny = c.dim == 3 ? c.n0 : 1;
    std::vector<MeshInput> meshes;
    amrex::Vector<fdsamr::MeshInfo> mi;
    const int nr = amrex::ParallelDescriptor::NProcs();
    for (int q = 0; q < 2; ++q) {   // two level-0 boxes, split in x
        MeshInput m;
        m.ijk[0] = c.n0 / 2; m.ijk[1] = ny; m.ijk[2] = c.n0;
        m.xb[0] = 0.5 * q; m.xb[1] = 0.5 * (q + 1); m.xb[2] = 0.0; m.xb[3] = ny * (1.0 / c.n0); m.xb[4] = 0.0; m.xb[5] = 1.0;
        m.rank = q % nr;
        meshes.push_back(m);
        fdsamr::MeshInfo f{};
        for (int d = 0; d < 3; ++d) f.ijk[d] = m.ijk[d];
        for (int d = 0; d < 6; ++d) f.xb[d] = m.xb[d];
        f.rank = m.rank;
        mi.push_back(f);
    }
    const bool ok = build_hierarchy_from_meshes(meshes, S.p, {true, true, true}, S.h, rep);
    if (!ok) { for (auto& e : rep.errors) std::fprintf(stderr, "hierarchy error: %s\n", e.c_str()); amrex::Abort("test_blob_registry: hierarchy"); }
    fdsamr::DomainInfo dom{};
    dom.periodic[0] = dom.periodic[1] = dom.periodic[2] = 1;
    dom.n_tracked = NS; dom.n_total = NS; dom.nranks = nr;
    S.l0 = fdsamr::assemble_level0(mi, dom);
    S.F0 = std::make_unique<fdsamr::Fields>(S.l0, NS);
    S.sd0 = std::make_unique<fdsamr::SideData>(S.l0, fdsamr::layout_cell_walls(S.l0));
    S.reg = std::make_unique<fdsamr::LevelRegistry>(S.l0.dom, NS);
    S.reg->adopt_level0(S.l0, *S.F0, *S.sd0);
    S.core = std::make_unique<RegridAmrCore>(S.h, S.p);
    S.rt = std::make_unique<fdsamr::RegistryTransfer>(*S.reg, NS);
    Sim* sp = &S;
    S.rt->initial.cell = [sp](double x, double y, double z, double* rho, double* tmp, double* zz) { cell_state(sp->c, x, y, z, rho, zz); *tmp = 300.0; };
    S.rt->initial.velocity = [sp](int dir, double x, double y, double z) { return vel_comp(sp->c, dir, x, y, z); };
    S.core->set_data_transfer(S.rt.get());
    S.core->set_tag_function([sp](int lev, amrex::TagBoxArray& tags, amrex::Real) {
        fdsamr::Fields& F = sp->reg->fields(lev);

        F["ZZ"].FillBoundary(sp->core->Geom(lev).periodicity());
        TagCriterion cr;
        cr.mode = TagMode::Above; cr.base = 0.05; cr.thr = 0.25; cr.keepfac = 0.8; cr.comp = 0;
        cr.dirs = {1, sp->c.dim == 3 ? 1 : 0, 1};
        amrex::iMultiFab cv = sp->core->covered_mask(lev);
        tag_cells(tags, F["ZZ"], nullptr, 0, &cv, cr);
    });
    // level 0 data, then the hierarchy from the tags at t = 0 (initial conditions evaluated on each new level, one average-down at the end)
    fill_analytic(S, 0, true, true);
    S.core->init_from_tags(*S.reg, 0.0, false, RegridAmrCore::DmFn(), &S.l0.dm);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(S.core->boxArray(0) == S.l0.ba, "level-0 grids of AmrCore and registry differ");
    for (int l = 0; l <= S.core->finestLevel(); ++l) { S.reg->fields(l)["RHO"].FillBoundary(S.core->Geom(l).periodicity()); S.reg->fields(l)["ZZ"].FillBoundary(S.core->Geom(l).periodicity()); }
    // the prescribed velocity of every level is the discrete curl on that level's own grid (divergence free to round-off on each level)
    for (int l = 0; l <= S.core->finestLevel(); ++l) fill_analytic(S, l, false, true);
    if (c.null_solver) S.psolver = std::make_unique<NullSolver>(); else S.psolver = std::make_unique<fdsrt::MlmgCompositeSolver>();
    if (c.project && S.core->finestLevel() > 0) projection_step(S);   // the t = 0 hierarchy: the coarse faces under the fine ones took the fine mean
    S.proj_n = S.proj_notaccepted = S.proj_iter = 0; S.proj_before = S.proj_after = S.proj_rhs_sum = 0; S.proj_change.clear(); std::fill(S.prof.begin(), S.prof.end(), 0.0);
    S.mass0 = composite(S);
}

// one transport step with the global dt: upwind fluxes of q_n = rho*Z_n on every level, interface flux overwrite finest first, update of every cell, rho = sum q_n, Z = q/rho,
// restriction (mass weighted, or linear for the negative control), same-level ghost refill.
void step(Sim& S)
{
    const Cfg& c = S.c;
    const int L = S.core->finestLevel();
    std::vector<amrex::MultiFab> q(L + 1);
    std::vector<std::array<amrex::MultiFab, 3>> flux(L + 1);
    const bool act[3] = {true, c.dim == 3, true};
    for (int l = 0; l <= L; ++l) {
        fdsamr::Fields& F = S.reg->fields(l);
        const amrex::BoxArray& ba = F["RHO"].boxArray();
        q[l].define(ba, F["RHO"].DistributionMap(), NS, 1);
        q[l].setVal(0.0);
        for (amrex::MFIter mfi(q[l]); mfi.isValid(); ++mfi) {
            auto o = q[l].array(mfi); auto r = F["RHO"].const_array(mfi); auto z = F["ZZ"].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < NS; ++n) o(i, j, k, n) = r(i, j, k) * z(i, j, k, n); });
        }
        q[l].FillBoundary(S.core->Geom(l).periodicity());
        if (l >= 1) fill_cf_ghosts_pc(q[l], q[l - 1], S.core->Geom(l), S.core->Geom(l - 1), S.core->refRatio(l - 1), 1, 0, NS);
        const char* vn[3] = {"U", "V", "W"};
        for (int d = 0; d < 3; ++d) {
            amrex::IntVect e(0); e[d] = 1;
            flux[l][d].define(amrex::convert(ba, e), F["RHO"].DistributionMap(), NS, 0);
            flux[l][d].setVal(0.0);
            if (!act[d]) continue;
            for (amrex::MFIter mfi(flux[l][d]); mfi.isValid(); ++mfi) {
                auto fl = flux[l][d].array(mfi); auto u = F[vn[d]].const_array(mfi); auto qa = q[l].const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    const double uf = u(i, j, k);
                    const int im = i - e[0], jm = j - e[1], km = k - e[2];
                    for (int n = 0; n < NS; ++n) fl(i, j, k, n) = uf * (uf > 0.0 ? qa(im, jm, km, n) : qa(i, j, k, n));
                });
            }
        }
    }
    if (c.overwrite)
        for (int l = L; l >= 1; --l)
            average_down_faces({&flux[l][0], &flux[l][1], &flux[l][2]}, {&flux[l - 1][0], &flux[l - 1][1], &flux[l - 1][2]}, S.core->refRatio(l - 1));
    for (int l = 0; l <= L; ++l) {
        fdsamr::Fields& F = S.reg->fields(l);
        const amrex::Geometry& g = S.core->Geom(l);
        for (amrex::MFIter mfi(q[l]); mfi.isValid(); ++mfi) {
            auto qa = q[l].array(mfi);
            auto rho = F["RHO"].array(mfi), rhos = F["RHOS"].array(mfi), zz = F["ZZ"].array(mfi), zzs = F["ZZS"].array(mfi);
            std::array<amrex::Array4<const double>, 3> fl{flux[l][0].const_array(mfi), flux[l][1].const_array(mfi), flux[l][2].const_array(mfi)};
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                double s = 0.0, qn[NS];
                for (int n = 0; n < NS; ++n) {
                    double v = qa(i, j, k, n);
                    for (int d = 0; d < 3; ++d) {
                        if (!act[d]) continue;
                        const int ii = i + (d == 0), jj = j + (d == 1), kk = k + (d == 2);
                        v -= S.dt / g.CellSize(d) * (fl[d](ii, jj, kk, n) - fl[d](i, j, k, n));
                    }
                    qn[n] = v; s += v;
                    if (v < 0.0) ++S.min_neg;
                }
                rho(i, j, k) = s; rhos(i, j, k) = s;
                for (int n = 0; n < NS; ++n) { zz(i, j, k, n) = qn[n] / s; zzs(i, j, k, n) = zz(i, j, k, n); }
            });
        }
    }
    if (c.linearZ) {   // negative control: what the restriction did before the species fix
        for (int l = L; l >= 1; --l) {
            fdsamr::Fields& Ff = S.reg->fields(l); fdsamr::Fields& Fc = S.reg->fields(l - 1);
            for (const char* nm : {"RHO", "ZZ"}) average_down_cells(Ff[nm], Fc[nm], S.core->Geom(l), S.core->Geom(l - 1), S.core->refRatio(l - 1), 0, Fc[nm].nComp());
        }
    } else {
        average_down_registry(*S.reg, {"RHO", "ZZ"});
    }
    for (int l = 0; l <= L; ++l) { S.reg->fields(l).fill_ghosts("RHO"); S.reg->fields(l).fill_ghosts("ZZ"); }
    S.t += S.dt;
}

double max_div(Sim& S, int l)
{
    fdsamr::Fields& F = S.reg->fields(l);
    double d = max_divergence_error_local({&F["U"], &F["V"], &F["W"]}, S.core->Geom(l), nullptr);
    amrex::ParallelDescriptor::ReduceRealMax(d);
    return d;
}

long cells_above_below_max(Sim& S)
{
    long n = 0;
    for (int l = 0; l < S.c.max_level && l <= S.core->finestLevel(); ++l) {
        const amrex::iMultiFab* cv = S.reg->covered_mask(l);
        fdsamr::Fields& F = S.reg->fields(l);
        for (amrex::MFIter mfi(F["ZZ"]); mfi.isValid(); ++mfi) {
            auto z = F["ZZ"].const_array(mfi);
            amrex::Array4<const int> m;
            if (cv) m = cv->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if ((!cv || !m(i, j, k)) && z(i, j, k, 0) - 0.05 > 0.25) ++n; });
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

// divergence of the level-l face velocities, cell by cell (one component)
amrex::MultiFab divergence(Sim& S, int l)
{
    fdsamr::Fields& F = S.reg->fields(l);
    amrex::MultiFab dv(F["RHO"].boxArray(), F["RHO"].DistributionMap(), 1, 0);
    const amrex::Geometry& g = S.core->Geom(l);
    for (amrex::MFIter mfi(dv); mfi.isValid(); ++mfi) {
        auto d = dv.array(mfi); auto u = F["U"].const_array(mfi); auto v = F["V"].const_array(mfi); auto w = F["W"].const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            d(i, j, k) = (u(i + 1, j, k) - u(i, j, k)) / g.CellSize(0) + (v(i, j + 1, k) - v(i, j, k)) / g.CellSize(1) + (w(i, j, k + 1) - w(i, j, k)) / g.CellSize(2);
        });
    }
    return dv;
}

// ---- D-063 post-regrid projection (mock composite MLMG solver) -------------------------------------------------------------------------------------------------
std::vector<fdsrt::ProjectionLevel> projection_levels(Sim& S)
{
    std::vector<fdsrt::ProjectionLevel> lv;
    for (int l = 0; l <= S.core->finestLevel(); ++l) {
        fdsamr::Fields& F = S.reg->fields(l);
        fdsrt::ProjectionLevel pl;
        pl.geom = S.core->Geom(l);
        pl.ref_ratio = l > 0 ? S.core->refRatio(l - 1) : amrex::IntVect(1);
        pl.vel = {&F["U"], &F["V"], &F["W"]};
        pl.covered = S.reg->covered_mask(l);
        lv.push_back(pl);
    }
    return lv;
}

// distance (in cells of the level, along the active directions) to the nearest cell with a source, up to 7; the profile of the velocity change against it
void add_profile(Sim& S, int l, const amrex::MultiFab& r, double rthr, const std::array<amrex::MultiFab, 3>& dv)
{
    const amrex::BoxArray cba = r.boxArray();
    amrex::iMultiFab dist(cba, r.DistributionMap(), 1, 1);
    dist.setVal(99);
    for (amrex::MFIter mfi(dist); mfi.isValid(); ++mfi) {
        auto d = dist.array(mfi); auto rr = r.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (std::abs(rr(i, j, k)) > rthr) d(i, j, k) = 0; });
    }
    const int ny = S.c.dim == 3 ? 1 : 0;
    for (int it = 1; it <= 7; ++it) {
        dist.FillBoundary(S.core->Geom(l).periodicity());
        for (amrex::MFIter mfi(dist); mfi.isValid(); ++mfi) {
            auto d = dist.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (d(i, j, k) <= it) return;
                const bool nb = d(i - 1, j, k) == it - 1 || d(i + 1, j, k) == it - 1 || d(i, j, k - 1) == it - 1 || d(i, j, k + 1) == it - 1 ||
                                (ny && (d(i, j - 1, k) == it - 1 || d(i, j + 1, k) == it - 1));
                if (nb) d(i, j, k) = it;
            });
        }
    }
    dist.FillBoundary(S.core->Geom(l).periodicity());
    for (int dir = 0; dir < 3; ++dir) {
        if (dir == 1 && !ny) continue;
        for (amrex::MFIter mfi(dv[dir]); mfi.isValid(); ++mfi) {
            auto a = dv[dir].const_array(mfi); auto d = dist.const_array(mfi);
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                const int dm = std::min(d(i, j, k), d(i - e[0], j - e[1], k - e[2]));
                const int b = std::min(dm, 7);
                S.prof[b] = std::max(S.prof[b], std::abs(a(i, j, k)));
            });
        }
    }
}

void projection_step(Sim& S)
{
    std::vector<fdsrt::ProjectionLevel> lv = projection_levels(S);
    const int nl = static_cast<int>(lv.size());
    const double scale = S.c.speed / total_cells_dx(S.c);
    fdsrt::ProjectionOptions po;
    po.enabled = S.c.project;
    po.tol_rel = 1.0e-11;
    po.max_iter = 100;
    po.tol_abs = 1.0e-10 * scale;
    po.accept_abs = 1.0e-9 * scale;     // the pressure tolerance of this mock run, against the divergence scale u / dx_fine
    po.accept_rel = 0.0;
    fdsrt::average_down_velocity(lv);
    std::vector<std::array<amrex::MultiFab, 3>> before(nl);
    std::vector<amrex::MultiFab> r(nl);
    std::vector<double> rmax(nl, 0.0);
    for (int l = 0; l < nl; ++l) {
        for (int d = 0; d < 3; ++d) {
            before[l][d].define(lv[l].vel[d]->boxArray(), lv[l].vel[d]->DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(before[l][d], *lv[l].vel[d], 0, 0, 1, 0);
        }
        r[l] = divergence(S, l);
        if (lv[l].covered)
            for (amrex::MFIter mfi(r[l]); mfi.isValid(); ++mfi) {
                auto a = r[l].array(mfi); auto m = lv[l].covered->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (m(i, j, k)) a(i, j, k) = 0.0; });
            }
        rmax[l] = r[l].norminf(0);
    }
    const fdsrt::ProjectionReport rep = fdsrt::project_after_regrid(lv, S.psolver.get(), po);
    ++S.proj_n;
    if (rep.ran && !rep.solved) { ++S.proj_failed; S.proj_msg = rep.message; }
    if (!rep.accepted) ++S.proj_notaccepted;
    S.proj_before = std::max(S.proj_before, rep.div_before);
    S.proj_after = std::max(S.proj_after, rep.div_after);
    S.proj_rhs_sum = std::max(S.proj_rhs_sum, rep.rhs_sum_rel);
    S.proj_iter = std::max(S.proj_iter, rep.iterations);
    if (S.proj_change.size() < static_cast<size_t>(nl)) S.proj_change.resize(nl, 0.0);
    if (rep.ran && rep.solved) {
        const char* st[3] = {"US", "VS", "WS"};
        const char* vn[3] = {"U", "V", "W"};
        for (int l = 0; l < nl; ++l) {
            fdsamr::Fields& F = S.reg->fields(l);
            std::array<amrex::MultiFab, 3> dv;
            for (int d = 0; d < 3; ++d) {
                amrex::MultiFab::Copy(F[st[d]], F[vn[d]], 0, 0, 1, 0);
                dv[d].define(before[l][d].boxArray(), before[l][d].DistributionMap(), 1, 0);
                amrex::MultiFab::Copy(dv[d], *lv[l].vel[d], 0, 0, 1, 0);
                amrex::MultiFab::Subtract(dv[d], before[l][d], 0, 0, 1, 0);
            }
            S.proj_change[l] = std::max(S.proj_change[l], rep.max_change_level[l]);
            if (rmax[l] > 1e-6 * scale) add_profile(S, l, r[l], 1e-3 * rmax[l], dv);
        }
    }
}

// after a regrid: for the cells of level l that have no cell of the OLD level l within one cell (all faces new): divergence against the parent's
void new_cell_divergence(Sim& S, int l, const amrex::BoxArray& old_ba, double* worst_abs, double* worst_parent, long* ncells)
{
    amrex::MultiFab df = divergence(S, l), dc = divergence(S, l - 1);
    const amrex::IntVect R = S.core->refRatio(l - 1);
    amrex::BoxArray cba = S.core->boxArray(l);
    cba.coarsen(R);
    amrex::MultiFab pd(cba, S.core->DistributionMap(l), 1, 0);
    pd.ParallelCopy(dc, 0, 0, 1, 0, 0, S.core->Geom(l - 1).periodicity());
    for (amrex::MFIter mfi(df); mfi.isValid(); ++mfi) {
        auto f = df.const_array(mfi); auto p = pd.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const amrex::Box nb(amrex::IntVect(i - 1, S.c.dim == 3 ? j - 1 : j, k - 1), amrex::IntVect(i + 1, S.c.dim == 3 ? j + 1 : j, k + 1));
            if (!old_ba.empty() && old_ba.intersects(nb)) return;
            *worst_abs = std::max(*worst_abs, std::abs(f(i, j, k)));
            *worst_parent = std::max(*worst_parent, std::abs(f(i, j, k) - p(amrex::coarsen(i, R[0]), amrex::coarsen(j, R[1]), amrex::coarsen(k, R[2]))));
            ++*ncells;
        });
    }
}

void run(Sim& S, bool verbose)
{
    const Cfg& c = S.c;
    // dt from the finest cell size and the largest velocity component sum (CFL-like factor 0.8 on the sum of the per-direction rates: positivity)
    const double dxf = total_cells_dx(c);
    const double umax = 1.25 * c.speed;
    S.dt = 0.8 * dxf / ((c.dim == 3 ? 3.0 : 2.0) * umax);
    int nreg = 0;
    for (int s = 1; s <= c.steps; ++s) {
        step(S);
        if (c.interval > 0 && s % c.interval == 0 && s < c.steps) {
            S.feature_late += cells_above_below_max(S);
            const std::vector<double> before = composite(S);
            std::vector<amrex::BoxArray> old_ba;
            for (int l = 0; l <= S.core->finestLevel(); ++l) old_ba.push_back(S.core->boxArray(l));
            const bool ch = S.core->regrid_dynamic(S.t);
            const std::vector<double> after = composite(S);
            ++nreg;
            if (ch) ++S.n_changed;
            double rel = 0;
            for (int n = 0; n < 1 + NS; ++n) rel = std::max(rel, std::abs(after[n] - before[n]) / std::abs(before[n]));
            S.worst_regrid = std::max(S.worst_regrid, rel);
            std::string dv, gs;
            double nd_abs = 0, nd_par = 0; long nn = 0;
            for (int l = 1; l <= S.core->finestLevel(); ++l) {
                new_cell_divergence(S, l, l < static_cast<int>(old_ba.size()) ? old_ba[l] : amrex::BoxArray(), &nd_abs, &nd_par, &nn);
                const double d = max_div(S, l);
                S.worst_div = std::max(S.worst_div, d);
                dv += " L" + std::to_string(l) + "=" + std::to_string(d);
                const GhostShareCount gc = count_shared_ghost_cells(S.core->boxArray(l), S.core->refRatio(l - 1), S.core->Geom(l - 1));
                gs += " L" + std::to_string(l) + ":" + std::to_string(gc.serving_two_faces) + "/" + std::to_string(gc.serving_several);
                for (int b = 0; b < static_cast<int>(S.core->boxArray(l).size()); ++b) {
                    const amrex::Box bx = S.core->boxArray(l)[b];
                    S.hier_hash = mixh(S.hier_hash, static_cast<uint64_t>(l * 1000003 + bx.smallEnd(0) * 7919 + bx.smallEnd(1) * 104729 + bx.smallEnd(2) * 1299709 + bx.bigEnd(0) * 15485863 + bx.bigEnd(1) * 32452843 + bx.bigEnd(2) * 49979687));
                }
            }
            amrex::ParallelDescriptor::ReduceRealMax(nd_abs); amrex::ParallelDescriptor::ReduceRealMax(nd_par); amrex::ParallelDescriptor::ReduceLongSum(nn);
            S.worst_new_div = std::max(S.worst_new_div, nd_abs); S.worst_parent = std::max(S.worst_parent, nd_par); S.n_new_cells += nn;
            std::string pj;
            if (ch && S.core->finestLevel() > 0) {   // D-063: the projection after a regrid that changed the grids (enabled by the case; otherwise only measured)
                const double b0 = S.proj_before;
                projection_step(S);
                pj = " | projection " + std::string(S.c.project ? "on" : "off (measured)") + ": composite max|div u - D| before " + std::to_string(S.proj_before) + " (run max, was " + std::to_string(b0) + "), after " + std::to_string(S.proj_after) + " (run max), solver iterations " + std::to_string(S.proj_iter);
            }
            if (verbose && amrex::ParallelDescriptor::IOProcessor())
                std::printf("    regrid %2d t=%.4f finest level %d%s: composite change %.1e, max|div u - D| by level (all cells incl. seams to retained fine data):%s, in all-new cells %.1e (child - parent %.1e, %ld cells), D-059 corner/edge cells (two-face/all shared):%s, clips so far %ld%s\n", nreg, S.t, S.core->finestLevel(),
                            ch ? " (grids changed)" : "", rel, dv.c_str(), nd_abs, nd_par, nn, gs.c_str(), S.rt->stats().clips, pj.c_str());
        }
    }
}

uint64_t data_hash(Sim& S)
{
    uint64_t d = 0;
    for (int l = 0; l <= S.core->finestLevel(); ++l)
        for (amrex::MFIter mfi(S.reg->fields(l)["RHO"]); mfi.isValid(); ++mfi) {
            auto r = S.reg->fields(l)["RHO"].const_array(mfi); auto z = S.reg->fields(l)["ZZ"].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                for (int n = -1; n < NS; ++n) {
                    double v = (n < 0) ? r(i, j, k) : z(i, j, k, n); uint64_t bits; std::memcpy(&bits, &v, 8);
                    d += mixh(mixh(mixh(mixh(static_cast<uint64_t>(l), static_cast<uint64_t>(i + 1000)), static_cast<uint64_t>(j + 1000)), static_cast<uint64_t>(k + 1000) * 31 + n + 1), bits);
                }
            });
        }
    unsigned long long dd = d;
    amrex::ParallelAllReduce::Sum(dd, amrex::ParallelContext::CommunicatorSub());
    return dd;
}

// composite rho*Z_1 of the hierarchy projected (piecewise constant) onto the uniform grid of the finest cell size, finer levels overwrite coarser ones
amrex::MultiFab project(Sim& S, const amrex::BoxArray& uba, const amrex::DistributionMapping& udm, int total_ratio)
{
    amrex::MultiFab out(uba, udm, 1, 0);
    out.setVal(0.0);
    const int L = S.core->finestLevel();
    int cum = total_ratio;
    for (int l = 0; l <= L; ++l) {
        // cumulative ratio from level l to the uniform grid
        int rr = total_ratio;
        for (int m = 1; m <= l; ++m) rr /= S.core->refRatio(m - 1)[0];
        (void)cum;
        const amrex::IntVect R(rr, S.c.dim == 3 ? rr : 1, rr);
        fdsamr::Fields& F = S.reg->fields(l);
        const amrex::BoxArray eba = amrex::refine(F["RHO"].boxArray(), R);
        amrex::MultiFab E(eba, F["RHO"].DistributionMap(), 1, 0);
        for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi) {
            auto e = E.array(mfi); auto r = F["RHO"].const_array(mfi); auto z = F["ZZ"].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                const amrex::IntVect c(amrex::coarsen(i, R[0]), amrex::coarsen(j, R[1]), amrex::coarsen(k, R[2]));
                e(i, j, k) = r(c[0], c[1], c[2]) * z(c[0], c[1], c[2], 0);
            });
        }
        out.ParallelCopy(E, 0, 0, 1, 0, 0);
    }
    return out;
}

double l1_diff(const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    amrex::MultiFab bb(ba, dm, 1, 0);
    bb.ParallelCopy(b, 0, 0, 1, 0, 0);
    double d = 0, n = 0;
    for (amrex::MFIter mfi(a); mfi.isValid(); ++mfi) {
        auto x = a.const_array(mfi); auto y = bb.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { d += std::abs(x(i, j, k) - y(i, j, k)); n += std::abs(y(i, j, k)); });
    }
    double v[2] = {d, n}; amrex::ParallelAllReduce::Sum(v, 2, amrex::ParallelContext::CommunicatorSub());
    return v[0] / v[1];
}

struct Outcome { double worst_regrid = 0, total_drift = 0, worst_div = 0, worst_new_div = 0, worst_parent = 0; long n_new = 0; double err_amr = 0, err_coarse = 0; long late = 0, neg = 0, clips = 0; int changed = 0; uint64_t hh = 0, dh = 0; int finest = 0;
    int proj_n = 0, proj_notaccepted = 0, proj_iter = 0, proj_failed = 0; std::string proj_msg; double proj_before = 0, proj_after = 0, proj_rhs_sum = 0; std::vector<double> proj_change, prof; };

Outcome dynamic_run(const Cfg& c, bool verbose, bool compare)
{
    Outcome o;
    Sim S;
    build(S, c);
    if (verbose && amrex::ParallelDescriptor::IOProcessor()) std::printf("  %s: initial hierarchy has %d level(s), boxes per level:", c.name.c_str(), S.core->finestLevel() + 1);
    if (verbose) { for (int l = 0; l <= S.core->finestLevel(); ++l) if (amrex::ParallelDescriptor::IOProcessor()) std::printf(" %zu", (size_t)S.core->boxArray(l).size()); if (amrex::ParallelDescriptor::IOProcessor()) std::printf("\n"); }
    run(S, verbose);
    const std::vector<double> end = composite(S);
    for (int n = 0; n < 1 + NS; ++n) o.total_drift = std::max(o.total_drift, std::abs(end[n] - S.mass0[n]) / std::abs(S.mass0[n]));
    o.worst_regrid = S.worst_regrid; o.worst_div = S.worst_div; o.worst_new_div = S.worst_new_div; o.worst_parent = S.worst_parent; o.n_new = S.n_new_cells; o.late = S.feature_late; o.neg = S.min_neg; o.clips = S.rt->stats().clips; o.changed = S.n_changed;
    o.hh = S.hier_hash; o.dh = data_hash(S); o.finest = S.core->finestLevel();
    o.proj_failed = S.proj_failed; o.proj_msg = S.proj_msg; o.proj_n = S.proj_n; o.proj_notaccepted = S.proj_notaccepted; o.proj_iter = S.proj_iter; o.proj_before = S.proj_before; o.proj_after = S.proj_after; o.proj_rhs_sum = S.proj_rhs_sum; o.proj_change = S.proj_change; o.prof = S.prof; amrex::ParallelDescriptor::ReduceRealMax(o.prof.data(), static_cast<int>(o.prof.size()));
    long neg = o.neg; amrex::ParallelDescriptor::ReduceLongSum(neg); o.neg = neg;
    if (compare) {
        // uniform-fine and uniform-coarse runs with the same dt and number of steps
        int total_ratio = 1; for (int l = 0; l < c.max_level; ++l) total_ratio *= c.ratio;
        Cfg cf = c; cf.project = false; cf.max_level = 0; cf.n0 = c.n0 * total_ratio; cf.interval = 0; cf.name = c.name + " uniform fine";
        Cfg cc = c; cc.project = false; cc.max_level = 0; cc.interval = 0; cc.name = c.name + " uniform coarse";
        Sim F;
        build(F, cf);
        F.dt = S.dt; F.t = 0;
        for (int s = 0; s < c.steps; ++s) step(F);
        Sim C;
        build(C, cc);
        C.dt = S.dt; C.t = 0;
        for (int s = 0; s < c.steps; ++s) step(C);
        const amrex::BoxArray uba = F.core->boxArray(0);
        const amrex::DistributionMapping udm = F.core->DistributionMap(0);
        const amrex::MultiFab ref = project(F, uba, udm, 1);
        const amrex::MultiFab a = project(S, uba, udm, total_ratio);
        const amrex::MultiFab cr = project(C, uba, udm, total_ratio);
        o.err_amr = l1_diff(a, ref, uba, udm);
        o.err_coarse = l1_diff(cr, ref, uba, udm);
    }
    return o;
}

// Corner check: static two-level hierarchy against the uniform-fine run, n steps, bitwise outside the dependence cone of the patch outline.
void corner_check(const Cfg& base, int nsteps)
{
    Cfg c = base; c.interval = 0; c.steps = nsteps; c.max_level = 1;
    Sim S; build(S, c);
    int ratio = c.ratio;
    Cfg cf = c; cf.max_level = 0; cf.n0 = c.n0 * ratio; cf.name = "uniform fine";
    Sim F; build(F, cf);
    const double dxf = total_cells_dx(c);
    S.dt = F.dt = 0.8 * dxf / ((c.dim == 3 ? 3.0 : 2.0) * 1.25 * c.speed);
    for (int s = 0; s < nsteps; ++s) { step(S); step(F); }
    fdsamr::Fields& A = S.reg->fields(1);
    const amrex::BoxArray fba = S.core->boxArray(1);
    amrex::MultiFab ref(fba, S.core->DistributionMap(1), 1 + NS, 0), cur(fba, S.core->DistributionMap(1), 1 + NS, 0), tmp(F.core->boxArray(0), F.core->DistributionMap(0), 1 + NS, 0);
    for (amrex::MFIter mfi(tmp); mfi.isValid(); ++mfi) {
        auto o = tmp.array(mfi); auto r = F.reg->fields(0)["RHO"].const_array(mfi); auto z = F.reg->fields(0)["ZZ"].const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { o(i, j, k, 0) = r(i, j, k); for (int n = 0; n < NS; ++n) o(i, j, k, 1 + n) = z(i, j, k, n); });
    }
    ref.ParallelCopy(tmp, 0, 0, 1 + NS, 0, 0);
    for (amrex::MFIter mfi(cur); mfi.isValid(); ++mfi) {
        auto o = cur.array(mfi); auto r = A["RHO"].const_array(mfi); auto z = A["ZZ"].const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { o(i, j, k, 0) = r(i, j, k); for (int n = 0; n < NS; ++n) o(i, j, k, 1 + n) = z(i, j, k, n); });
    }
    long deep = 0, deep_bad = 0, near = 0, near_diff = 0;
    double near_max = 0;
    const int rad = nsteps + 1;
    for (amrex::MFIter mfi(cur); mfi.isValid(); ++mfi) {
        auto a = cur.const_array(mfi); auto b = ref.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            amrex::Box nb(amrex::IntVect(i - rad, c.dim == 3 ? j - rad : j, k - rad), amrex::IntVect(i + rad, c.dim == 3 ? j + rad : j, k + rad));
            const bool inside = fba.contains(nb);
            double dmax = 0;
            for (int n = 0; n < 1 + NS; ++n) dmax = std::max(dmax, std::abs(a(i, j, k, n) - b(i, j, k, n)));
            if (inside) { ++deep; if (dmax != 0.0) ++deep_bad; }
            else { ++near; if (dmax > 1e-12) ++near_diff; near_max = std::max(near_max, dmax); }
        });
    }
    long v[4] = {deep, deep_bad, near, near_diff}; amrex::ParallelDescriptor::ReduceLongSum(v, 4); amrex::ParallelDescriptor::ReduceRealMax(near_max);
    if (amrex::ParallelDescriptor::IOProcessor())
        std::printf("  corner check (%s, %d steps, static patch): %ld fine cells farther than %d cells from the patch outline (corners included) differ from the uniform-fine run in %ld cells (bitwise test); "
                    "within the cone %ld cells, %ld of them differ by more than 1e-12 (max %.1e, the coarse-ghost effect)\n", base.name.c_str(), nsteps, v[0], rad, v[1], v[2], v[3], near_max);
    CHECK_MSG(v[0] > 0, base.name + ": corner check has cells beyond the dependence cone");
    CHECK_MSG(v[1] == 0, base.name + ": no corner effect: fine cells beyond the dependence cone of the patch outline equal the uniform-fine run bitwise, bad=" + std::to_string(v[1]));
    CHECK_MSG(v[3] > 0, base.name + ": the check can fail: cells inside the cone do feel the coarse ghost values (" + std::to_string(v[3]) + " differ)");
}

Outcome check_run(const Cfg& c, bool compare, bool print_hash, uint64_t* hh, uint64_t* dh)
{
    const Outcome o = dynamic_run(c, true, compare);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::printf("  %s: finest level %d, grids changed in %d regrids, composite change per regrid <= %.1e, whole run %.1e, cells above threshold below max level at a regrid %ld, negative rho*Z %ld, clips %ld, "
                    "max|div u - D| after a regrid: %.2e over all cells (seams to retained fine data), %.2e over %ld all-new cells, child - parent %.1e (u/dx scale %.1f)\n", c.name.c_str(), o.finest, o.changed, o.worst_regrid,
                    o.total_drift, o.late, o.neg, o.clips, o.worst_div, o.worst_new_div, o.n_new, o.worst_parent, c.speed / total_cells_dx(c));
        if (compare) std::printf("    L1 difference of rho*Z_1 to the uniform-fine run: AMR %.3e, uniform coarse %.3e\n", o.err_amr, o.err_coarse);
    }
    CHECK_MSG(o.finest == c.max_level, c.name + ": the hierarchy reaches MAX_LEVEL");
    CHECK_MSG(o.changed >= 3, c.name + ": the hierarchy follows the blob (grids changed in " + std::to_string(o.changed) + " regrids)");
    CHECK_MSG(o.worst_regrid <= 1e-12, c.name + ": composite mass and species change per regrid <= 1e-12 relative, got " + std::to_string(o.worst_regrid));
    CHECK_MSG(o.total_drift <= 1e-12, c.name + ": composite mass and species over the whole run <= 1e-12 relative, got " + std::to_string(o.total_drift));
    CHECK_MSG(o.n_new > 0 && o.worst_parent <= 1e-11 * c.speed / total_cells_dx(c), c.name + ": D-060: every all-new fine cell has the divergence of its parent to round-off, got " + std::to_string(o.worst_parent) + " over " + std::to_string(o.n_new) + " cells");
    CHECK_MSG(o.neg == 0, c.name + ": no negative rho*Z_n");
    CHECK_MSG(o.clips == 0, c.name + ": zero clips of the prolongation on positive data");
    CHECK_MSG(o.late == 0, c.name + ": every cell above the tag threshold lies on the finest level at the next regrid, missing " + std::to_string(o.late));
    if (compare) CHECK_MSG(o.err_amr < o.err_coarse, c.name + ": the AMR result is closer to the uniform-fine run than the uniform-coarse one");
    if (hh) { *hh = o.hh; *dh = o.dh; }
    (void)print_hash;
    return o;
}

// D-063: the same case with the post-regrid composite projection (mock solver). Bound = the pressure tolerance of the mock run, expressed against the divergence scale u/dx_fine.
// g_assert_projection (command line --report-only turns it off) is the single switch between "numbers reported" and "bound asserted"; the composite conservation checks are always asserted.
bool g_assert_projection = true;

void check_projection(const Cfg& c0, const Outcome& off, uint64_t* hier_hash)
{
    Cfg c = c0; c.project = true; c.name = c0.name + " + projection";
    const Outcome o = dynamic_run(c, true, true);
    const double scale = c.speed / total_cells_dx(c);
    const double bound = 1.0e-9 * scale;
    if (o.proj_failed == 0 && amrex::ParallelDescriptor::IOProcessor()) {
        std::printf("  %s [mock solver: NOT a Phase 3 gate until the real composite solver lands]: %d projections, composite max|div u - D| before %.2e, after %.2e (bound %.1e = 1e-9 u/dx), solver iterations <= %d, "
                    "|sum rhs| / sum|rhs| <= %.1e\n", c.name.c_str(), o.proj_n, o.proj_before, o.proj_after, bound, o.proj_iter, o.proj_rhs_sum);
        std::printf("    largest velocity change of the projection by level:");
        for (double v : o.proj_change) std::printf(" %.2e", v);
        std::printf(" (u = %.1f, u/dx_fine = %.1f)\n", c.speed, scale);
        std::printf("    velocity change against the distance to the nearest divergence source (cells of the level; max over regrids and levels):");
        for (int b = 0; b < 8; ++b) std::printf(" d%s%d: %.2e", b == 7 ? ">=" : "=", b, o.prof[b]);
        std::printf("\n");
        std::printf("    conservation: composite change per regrid %.1e (without projection %.1e), whole run %.1e (without %.1e); L1 of rho*Z_1 to the uniform-fine run: AMR+projection %.3e (uniform coarse %.3e)\n",
                    o.worst_regrid, off.worst_regrid, o.total_drift, off.total_drift, o.err_amr, o.err_coarse);
        std::printf("    control, projection off, same case: max|div u - D| before %.2e, accepted in %d of %d regrids\n", off.proj_before, off.proj_n - off.proj_notaccepted, off.proj_n);
    }
    if (o.proj_failed > 0) {   // no projection support for this case in the mock solver: reported, not a gate
        if (amrex::ParallelDescriptor::IOProcessor())
            std::printf("  %s: NOT GATE: the mock composite solver did not converge in %d of %d projections (%s); the case has no projection support yet, only its transport conservation is asserted\n",
                        c.name.c_str(), o.proj_failed, o.proj_n, o.proj_msg.c_str());
        CHECK_MSG(o.worst_regrid <= 1e-12 && o.total_drift <= 1e-12, c.name + ": composite mass and species conserved (projection not applied)");
        return;
    }
    CHECK_MSG(o.proj_n > 0 && off.proj_n > 0, c.name + ": projections were done");
    CHECK_MSG(o.worst_regrid <= 1e-12 && o.total_drift <= 1e-12, c.name + ": the projection does not perturb composite mass and species (per regrid " + std::to_string(o.worst_regrid) + ", run " + std::to_string(o.total_drift) + ")");
    CHECK_MSG(o.neg == 0 && o.clips == 0, c.name + ": no negative rho*Z and no clips with the projection");
    CHECK_MSG(o.err_amr < o.err_coarse, c.name + ": with the projection the AMR result is still closer to the uniform-fine run than the uniform-coarse one");
    CHECK_MSG(off.proj_before > bound, c.name + ": control: without the projection the divergence at the seams exceeds the bound (" + std::to_string(off.proj_before) + " > " + std::to_string(bound) + ")");
    CHECK_MSG(off.proj_notaccepted == off.proj_n, c.name + ": control: the acceptance test fails for every regrid when the projection is off");
    if (g_assert_projection) {
        CHECK_MSG(o.proj_notaccepted == 0, c.name + ": D-063 acceptance: max|div u - D| <= bound on all uncovered cells after every projection (not accepted: " + std::to_string(o.proj_notaccepted) + ")");
        CHECK_MSG(o.proj_after <= bound, c.name + ": D-063 acceptance: max|div u - D| after projection " + std::to_string(o.proj_after) + " <= " + std::to_string(bound));
    }
    if (hier_hash) *hier_hash = o.hh;
}

}  // namespace

int main(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--report-only") g_assert_projection = false;
    int one = 1;
    amrex::Initialize(one, argv);
    uint64_t hh = 0, dh = 0;
    {
        Cfg a; a.name = "3-D 32^3, 2 levels, ratio 2"; a.dim = 3; a.n0 = 32; a.max_level = 1; a.ratio = 2; a.steps = 36;
        const Outcome oa = check_run(a, true, true, &hh, &dh);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("HASH hier=%016llx data=%016llx\n", (unsigned long long)hh, (unsigned long long)dh);
        uint64_t ph = 0;
        check_projection(a, oa, &ph);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("HASH_PROJ hier=%016llx\n", (unsigned long long)ph);
        {
            Cfg nc = a; nc.project = true; nc.null_solver = true; nc.name = "control: projection with a solver that corrects nothing";
            const Outcome on = dynamic_run(nc, false, false);
            if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  control null solver: max|div u - D| after %.2e, accepted in %d of %d regrids\n", on.proj_after, on.proj_n - on.proj_notaccepted, on.proj_n);
            CHECK_MSG(on.proj_n > 0 && on.proj_notaccepted == on.proj_n, "control: the acceptance test fails when the solver corrects nothing");
        }
        Cfg b; b.name = "2-D 48x1x48, 3 levels, ratio 2"; b.dim = 2; b.n0 = 48; b.max_level = 2; b.ratio = 2; b.steps = 36;
        const Outcome ob = check_run(b, true, false, nullptr, nullptr);
        check_projection(b, ob, nullptr);
        Cfg d; d.name = "2-D 32x1x32, 2 levels, ratio 4"; d.dim = 2; d.n0 = 32; d.max_level = 1; d.ratio = 4; d.steps = 36;
        const Outcome od = check_run(d, true, false, nullptr, nullptr);
        check_projection(d, od, nullptr);

        // corner check
        corner_check(a, 5);
        corner_check(d, 5);

        // negative controls (3-D unless noted)
        Cfg n1 = a; n1.name = "control: no flux overwrite"; n1.overwrite = false;
        const Outcome o1 = dynamic_run(n1, false, false);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  control no flux overwrite: composite drift over the run %.2e\n", o1.total_drift);
        CHECK_MSG(o1.total_drift > 1e-9, "control: without the interface flux overwrite the composite mass drifts (" + std::to_string(o1.total_drift) + ")");
        Cfg n2 = a; n2.name = "control: linear Z restriction"; n2.linearZ = true;
        const Outcome o2 = dynamic_run(n2, false, false);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  control linear Z average-down: worst change over a regrid %.2e, over the run %.2e\n", o2.worst_regrid, o2.total_drift);
        CHECK_MSG(std::max(o2.worst_regrid, o2.total_drift) > 1e-9, "control: a linear average of Z breaks the species mass (" + std::to_string(std::max(o2.worst_regrid, o2.total_drift)) + ")");
        Cfg n3 = a; n3.name = "control: no tag buffer"; n3.nbuf = 0;
        const Outcome o3 = dynamic_run(n3, false, false);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  control no tag buffer: cells above threshold off the finest level %ld\n", o3.late);
        CHECK_MSG(o3.late > 0, "control: without a tag buffer the feature leaves the finest level");
    }
    const long nfail = fdstest::report("test_blob_registry");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
