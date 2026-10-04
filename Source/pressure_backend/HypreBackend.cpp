// HYPRE assembled-matrix backend: see HypreBackend.H for the operator and frozen/hypre-notes.md for the derivation, the checks and
// the numbers. The matrix is built cell by cell from the layout (BoxArrays, DistributionMappings, geometries, ratios, boundary types)
// with a formula-based global numbering, so no communication is needed to assemble it.
#include "HypreBackend.H"

// The whole file is compiled only with PB_WITH_HYPRE (README.md, "Build options"); without it the object is empty and no HYPRE symbol is
// referenced anywhere in the pressure backend.
#ifdef PB_WITH_HYPRE
#include "PressureBackend.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_krylov.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>

namespace pb {

using namespace amrex;

struct HypreSystem::Impl {
    std::vector<HypreLayoutLevel> lev;
    std::array<BC,6> ebc;
    HypreOptions opt;
    bool singular = false, pin_enabled = true;
    int nlev = 0;
    bool ok = true;
    std::string message;
    std::array<bool,3> active{{true,true,true}};      // direction carries an operator term (more than one cell)
    std::vector<std::array<double,3>> dx;
    std::vector<double> vol, dscale;                  // cell volume; typical diagonal 2*V*sum_d 1/h_d^2
    std::vector<std::vector<long long>> base;         // first row of each box
    std::vector<BoxArray> cfine;                      // cfine[l]: cells of level l covered by level l+1
    long long ilower = 0, iupper = -1, nglobal = 0;
    bool has_pin = false;
    int pin_lev = 0;
    IntVect pin_iv{0};
    long long pin_id = -1;
    long long nnz_local = 0, nnz_global = 0;
    double setup_s = 0.0;
    int kind = 0;                                     // 0 PCG, 1 GMRES, 2 BiCGSTAB
    int rank = 0;

    HYPRE_IJMatrix A = nullptr;
    HYPRE_ParCSRMatrix parA = nullptr;
    HYPRE_IJVector bvec = nullptr, xvec = nullptr;
    HYPRE_ParVector parb = nullptr, parx = nullptr;
    HYPRE_Solver solver = nullptr, precond = nullptr;

    ~Impl ();
    void build ();
    void assemble ();
    void make_solver ();

    // ---- layout helpers ----------------------------------------------------------------------------------------------
    static Long lin (Box const& b, IntVect const& iv)
    {
        return Long(iv[0]-b.smallEnd(0)) + Long(b.length(0))*(Long(iv[1]-b.smallEnd(1)) + Long(b.length(1))*Long(iv[2]-b.smallEnd(2)));
    }
    bool wrap (int l, IntVect& iv) const
    {
        Box const& dom = lev[l].geom.Domain();
        for (int d = 0; d < 3; ++d) {
            if (iv[d] < dom.smallEnd(d) || iv[d] > dom.bigEnd(d)) {
                if (!lev[l].geom.isPeriodic(d)) { return false; }
                const int len = dom.length(d);
                iv[d] = dom.smallEnd(d) + (((iv[d] - dom.smallEnd(d)) % len) + len) % len;
            }
        }
        return true;
    }
    // Row of the cell (any box on level l), or -1 if the cell is outside the domain (non-periodic) or not on the level.
    long long id_at (int l, IntVect iv) const
    {
        if (!wrap(l, iv)) { return -1; }
        auto isects = lev[l].ba.intersections(Box(iv, iv), true, IntVect(0));
        if (isects.empty()) { return -1; }
        const int b = isects[0].first;
        return base[l][b] + lin(lev[l].ba[b], iv);
    }
    bool covered_wrapped (int l, IntVect const& wiv) const { return l + 1 < nlev && cfine[l].contains(wiv); }

    // Cell of level L-1 usable by the interpolation to level L: on the level, in the domain, not covered by level L.
    bool avail (int L, IntVect c, long long& id) const
    {
        if (!wrap(L-1, c)) { return false; }
        if (cfine[L-1].contains(c)) { return false; }
        id = id_at(L-1, c);
        return id >= 0;
    }

    struct W { long long id; double w; };
    static void addw (std::vector<W>& v, long long id, double w)
    {
        for (auto& e : v) { if (e.id == id) { e.w += w; return; } }
        v.push_back({id, w});
    }
    // Boundary value at the fine ghost cell g (index space of level L, normal direction d) as a combination of level L-1 cells:
    // AMReX interpbndrydata_{x,y,z}_o3. False if the base cell is not available (nesting violated).
    bool bdry_stencil (int L, IntVect const& g, int d, std::vector<W>& out) const
    {
        out.clear();
        IntVect const r = lev[L].ratio;
        IntVect cc;
        for (int k = 0; k < 3; ++k) { cc[k] = amrex::coarsen(g[k], r[k]); }
        long long id0 = -1;
        if (!avail(L, cc, id0)) { return false; }
        addw(out, id0, 1.0);
        int t[2], nt = 0;
        for (int k = 0; k < 3; ++k) { if (k != d) { t[nt++] = k; } }
        int lo[2] = {0,0}, hi[2] = {0,0};
        double pos[2] = {0.0,0.0};
        bool act[2];
        long long idl[2] = {-1,-1}, idh[2] = {-1,-1};
        for (int m = 0; m < 2; ++m) {
            act[m] = active[t[m]];
            if (!act[m]) { continue; }
            IntVect cm = cc, cp = cc; cm[t[m]] -= 1; cp[t[m]] += 1;
            long long a = -1;
            if (avail(L, cm, a)) { lo[m] = -1; idl[m] = a; }
            if (avail(L, cp, a)) { hi[m] = 1; idh[m] = a; }
            pos[m] = -0.5 + (double(g[t[m]] - cc[t[m]]*r[t[m]]) + 0.5)/double(r[t[m]]);
            const double fac = (hi[m] == lo[m] + 1) ? 1.0 : 0.5;
            if (hi[m] != lo[m]) {
                const long long ih = (hi[m] == 1) ? idh[m] : id0, il = (lo[m] == -1) ? idl[m] : id0;
                addw(out, ih, pos[m]*fac);
                addw(out, il, -pos[m]*fac);
            }
            if (hi[m] == lo[m] + 2) {
                const double c2 = pos[m]*pos[m]*0.5;
                addw(out, idh[m], c2); addw(out, id0, -2.0*c2); addw(out, idl[m], c2);
            }
        }
        if (act[0] && act[1]) {
            long long ids[2][2];            // [sign of t0][sign of t1], index 0 = minus, 1 = plus
            bool all = true;
            for (int a = 0; a < 2 && all; ++a) {
                for (int b = 0; b < 2 && all; ++b) {
                    IntVect c = cc; c[t[0]] += (a ? 1 : -1); c[t[1]] += (b ? 1 : -1);
                    long long id = -1;
                    if (!avail(L, c, id)) { all = false; } else { ids[a][b] = id; }
                }
            }
            if (all) {
                const double f = 0.25*pos[0]*pos[1];
                addw(out, ids[1][1], f); addw(out, ids[0][1], -f); addw(out, ids[0][0], f); addw(out, ids[1][0], -f);
            }
        }
        return true;
    }
};

namespace {
void hcheck (int ierr, const char* what, std::string& msg, bool& ok)
{
    if (ierr != 0 && ok) {
        std::ostringstream o; o << "HYPRE call " << what << " returned " << ierr;
        msg = o.str(); ok = false;
        HYPRE_ClearAllErrors();
    }
}
}

HypreSystem::Impl::~Impl ()
{
    if (solver) {
        if (kind == 0) { HYPRE_ParCSRPCGDestroy(solver); }
        else if (kind == 1) { HYPRE_ParCSRGMRESDestroy(solver); }
        else { HYPRE_ParCSRBiCGSTABDestroy(solver); }
    }
    if (precond) { HYPRE_BoomerAMGDestroy(precond); }
    if (A) { HYPRE_IJMatrixDestroy(A); }
    if (bvec) { HYPRE_IJVectorDestroy(bvec); }
    if (xvec) { HYPRE_IJVectorDestroy(xvec); }
}

void HypreSystem::Impl::build ()
{
    const double t0 = amrex::second();
    nlev = static_cast<int>(lev.size());
    rank = ParallelDescriptor::MyProc();
    const int nprocs = ParallelDescriptor::NProcs();
    Box const& dom0 = lev[0].geom.Domain();
    for (int d = 0; d < 3; ++d) { active[d] = (dom0.length(d) > 1); }
    dx.resize(nlev); vol.resize(nlev); dscale.resize(nlev);
    for (int l = 0; l < nlev; ++l) {
        double v = 1.0, s = 0.0;
        for (int d = 0; d < 3; ++d) { dx[l][d] = lev[l].geom.CellSize(d); v *= dx[l][d]; }
        vol[l] = v;
        for (int d = 0; d < 3; ++d) { if (active[d]) { s += 2.0*v/(dx[l][d]*dx[l][d]); } }
        dscale[l] = s > 0.0 ? s : 1.0;
    }
    cfine.resize(nlev);
    for (int l = 0; l + 1 < nlev; ++l) { cfine[l] = lev[l+1].ba; cfine[l].coarsen(lev[l+1].ratio); }
    // Numbering: rank-major, then level, then box index, then x-fastest inside the box. All covered cells keep a row.
    std::vector<long long> count(nprocs, 0);
    for (int l = 0; l < nlev; ++l) {
        for (int b = 0; b < static_cast<int>(lev[l].ba.size()); ++b) { count[lev[l].dm[b]] += lev[l].ba[b].numPts(); }
    }
    std::vector<long long> start(nprocs + 1, 0);
    for (int r = 0; r < nprocs; ++r) { start[r+1] = start[r] + count[r]; }
    nglobal = start[nprocs];
    if (nglobal >= (long long)std::numeric_limits<HYPRE_Int>::max() / 2) {
        ok = false; message = "HYPRE backend: more rows than the HYPRE integer type holds (needs a bigint HYPRE build)"; return;
    }
    ilower = start[rank]; iupper = start[rank+1] - 1;
    std::vector<long long> cur(start.begin(), start.end() - 1);
    base.assign(nlev, {});
    for (int l = 0; l < nlev; ++l) {
        base[l].resize(lev[l].ba.size());
        for (int b = 0; b < static_cast<int>(lev[l].ba.size()); ++b) {
            const int r = lev[l].dm[b];
            base[l][b] = cur[r]; cur[r] += lev[l].ba[b].numPts();
        }
    }
    // Pin (singular problems): lowest index (x fastest) uncovered cell of the coarsest level that has uncovered cells.
    has_pin = singular && pin_enabled;
    if (has_pin) {
        int pl = -1;
        for (int l = 0; l < nlev && pl < 0; ++l) {
            Long ncov = (l + 1 < nlev) ? cfine[l].numPts() : 0;
            if (lev[l].ba.numPts() - ncov > 0) { pl = l; }
        }
        pin_lev = pl;
        Box const& dom = lev[pl].geom.Domain();
        const Long nx = dom.length(0), ny = dom.length(1);
        Long lowest = std::numeric_limits<Long>::max();
        for (int b = 0; b < static_cast<int>(lev[pl].ba.size()); ++b) {
            if (lev[pl].dm[b] != rank) { continue; }
            Box const& bx = lev[pl].ba[b];
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const Long idx = Long(i - dom.smallEnd(0)) + nx*(Long(j - dom.smallEnd(1)) + ny*Long(k - dom.smallEnd(2)));
                if (idx >= lowest) { continue; }
                if (pl + 1 < nlev && cfine[pl].contains(IntVect(i,j,k))) { continue; }
                lowest = idx;
            }
        }
        ParallelDescriptor::ReduceLongMin(lowest);
        pin_iv = IntVect(int(lowest % nx) + dom.smallEnd(0), int((lowest / nx) % ny) + dom.smallEnd(1), int(lowest / (nx*ny)) + dom.smallEnd(2));
        pin_id = id_at(pin_lev, pin_iv);
    }
    assemble();
    if (!ok) { return; }
    make_solver();
    setup_s = amrex::second() - t0;
}

void HypreSystem::Impl::assemble ()
{
    struct Ent { long long col; double val; };
    std::vector<HYPRE_Int> ncols;
    std::vector<HYPRE_BigInt> rowids, cols;
    std::vector<double> vals;
    std::vector<Ent> ent;
    std::vector<W> st;
    const long long nloc = iupper - ilower + 1;
    ncols.reserve(std::max<long long>(nloc, 0)); rowids.reserve(std::max<long long>(nloc, 0));
    cols.reserve(std::max<long long>(nloc, 0)*7); vals.reserve(std::max<long long>(nloc, 0)*7);
    auto bad = [&] (const char* what) { if (ok) { ok = false; message = std::string("HYPRE backend: ") + what; } };

    for (int l = 0; l < nlev; ++l) {
        for (int b = 0; b < static_cast<int>(lev[l].ba.size()); ++b) {
            if (lev[l].dm[b] != rank) { continue; }
            Box const& bx = lev[l].ba[b];
            // uncovered flags of this box
            std::vector<char> unc(bx.numPts(), 1);
            if (l + 1 < nlev) {
                for (auto const& is : cfine[l].intersections(bx)) {
                    Box const& ib = is.second;
                    for (int k = ib.smallEnd(2); k <= ib.bigEnd(2); ++k) for (int j = ib.smallEnd(1); j <= ib.bigEnd(1); ++j) for (int i = ib.smallEnd(0); i <= ib.bigEnd(0); ++i) {
                        unc[lin(bx, IntVect(i,j,k))] = 0;
                    }
                }
            }
            IntVect const rf = (l + 1 < nlev) ? lev[l+1].ratio : IntVect(1);
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const IntVect c(i,j,k);
                const long long myid = base[l][b] + lin(bx, c);
                ent.clear();
                double diag = 0.0;
                if (!unc[lin(bx, c)]) {
                    diag = dscale[l];                                   // covered: decoupled row
                } else if (has_pin && myid == pin_id) {
                    diag = dscale[l];                                   // pin: identity row (scaled by the local diagonal), rhs 0
                } else {
                    for (int d = 0; d < 3; ++d) {
                        if (!active[d]) { continue; }
                        const double h = dx[l][d];
                        const double kk = vol[l]/(h*h);
                        for (int side = 0; side < 2; ++side) {
                            IntVect n = c; n[d] += side ? 1 : -1;
                            // classify the neighbour
                            enum { STD, PHYS, COVERED, CF } cls;
                            long long nid = -1;
                            if (bx.contains(n)) {
                                if (unc[lin(bx, n)]) { cls = STD; nid = base[l][b] + lin(bx, n); } else { cls = COVERED; }
                            } else {
                                IntVect w = n;
                                if (!wrap(l, w)) { cls = PHYS; }
                                else {
                                    nid = id_at(l, w);
                                    if (nid < 0) { cls = CF; }
                                    else if (covered_wrapped(l, w)) { cls = COVERED; }
                                    else { cls = STD; }
                                }
                            }
                            if (cls == STD) {
                                diag += kk; ent.push_back({nid, -kk});
                            } else if (cls == PHYS) {
                                if (ebc[face_index(d, side)] == BC::Dirichlet) { diag += 2.0*kk; }
                            } else if (cls == CF) {
                                if (l == 0 || !bdry_stencil(l, n, d, st)) { bad("coarse/fine interpolation stencil not available (nesting)"); continue; }
                                // MLMG: ghost = a*b + (1-a)*phi_c with a = 2/(r+1): linear extrapolation through (-r/2, b) and (0.5, phi_c)
                                const double ra = 2.0/(double(lev[l].ratio[d]) + 1.0);
                                diag += ra*kk;
                                for (auto const& e : st) { ent.push_back({e.id, -ra*kk*e.w}); }
                            } else {   // COVERED: reflux. Average of the fine face fluxes 2*(phi_f - b_f)/h_f replaces the coarse face flux.
                                IntVect wn = n; wrap(l, wn);
                                double F = 1.0;
                                for (int t = 0; t < 3; ++t) { if (t != d && active[t]) { F *= double(rf[t]); } }
                                const double ra = 2.0/(double(rf[d]) + 1.0);
                                const double coef = ra*vol[l]/(h*F*dx[l+1][d]);
                                // fine index ranges of the face layer
                                const int rd = rf[d];
                                const int f_n = side ? wn[d]*rd : wn[d]*rd + rd - 1;       // fine cell inside the covered neighbour, next to the face
                                const int f_g = side ? c[d]*rd + rd - 1 : c[d]*rd;          // fine ghost cell (inside c), next to the face
                                int t0 = (d + 1) % 3, t1 = (d + 2) % 3;
                                if (t0 > t1) { std::swap(t0, t1); }
                                const int n0 = active[t0] ? rf[t0] : 1, n1 = active[t1] ? rf[t1] : 1;
                                for (int q1 = 0; q1 < n1; ++q1) for (int q0 = 0; q0 < n0; ++q0) {
                                    IntVect fn, fg;
                                    fn[d] = f_n; fg[d] = f_g;
                                    fn[t0] = fg[t0] = c[t0]*rf[t0] + q0;
                                    fn[t1] = fg[t1] = c[t1]*rf[t1] + q1;
                                    const long long fid = id_at(l + 1, fn);
                                    if (fid < 0 || !bdry_stencil(l + 1, fg, d, st)) { bad("reflux face data not available (nesting)"); continue; }
                                    ent.push_back({fid, -coef});
                                    for (auto const& e : st) { ent.push_back({e.id, coef*e.w}); }
                                }
                            }
                        }
                    }
                }
                // finalize: merge duplicates in column order, add the diagonal, drop the pinned column
                std::stable_sort(ent.begin(), ent.end(), [] (Ent const& a, Ent const& b2) { return a.col < b2.col; });
                std::vector<Ent> row;
                bool diag_done = false;
                auto push = [&] (long long col, double v) {
                    if (has_pin && col == pin_id && myid != pin_id) { return; }
                    row.push_back({col, v});
                };
                for (std::size_t m = 0; m < ent.size();) {
                    std::size_t m2 = m; double s = 0.0;
                    while (m2 < ent.size() && ent[m2].col == ent[m].col) { s += ent[m2].val; ++m2; }
                    if (ent[m].col == myid) { s += diag; diag_done = true; }
                    push(ent[m].col, s);
                    m = m2;
                }
                if (!diag_done) {
                    // diagonal entry in column order
                    Ent e{myid, diag};
                    auto it = std::lower_bound(row.begin(), row.end(), e, [] (Ent const& a, Ent const& b2) { return a.col < b2.col; });
                    row.insert(it, e);
                }
                rowids.push_back(HYPRE_BigInt(myid));
                ncols.push_back(HYPRE_Int(row.size()));
                for (auto const& e : row) { cols.push_back(HYPRE_BigInt(e.col)); vals.push_back(e.val); }
            }
        }
    }
    nnz_local = static_cast<long long>(vals.size());
    { Long n = nnz_local; ParallelDescriptor::ReduceLongSum(n); nnz_global = n; }
    if (!ok) { return; }
    MPI_Comm comm = ParallelDescriptor::Communicator();
    int ierr = 0;
    ierr = HYPRE_IJMatrixCreate(comm, HYPRE_BigInt(ilower), HYPRE_BigInt(iupper), HYPRE_BigInt(ilower), HYPRE_BigInt(iupper), &A); hcheck(ierr, "IJMatrixCreate", message, ok);
    ierr = HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR); hcheck(ierr, "IJMatrixSetObjectType", message, ok);
    ierr = HYPRE_IJMatrixSetRowSizes(A, ncols.data()); hcheck(ierr, "IJMatrixSetRowSizes", message, ok);
    ierr = HYPRE_IJMatrixInitialize(A); hcheck(ierr, "IJMatrixInitialize", message, ok);
    if (!rowids.empty()) {
        ierr = HYPRE_IJMatrixSetValues(A, HYPRE_Int(rowids.size()), ncols.data(), rowids.data(), cols.data(), vals.data());
        hcheck(ierr, "IJMatrixSetValues", message, ok);
    }
    ierr = HYPRE_IJMatrixAssemble(A); hcheck(ierr, "IJMatrixAssemble", message, ok);
    ierr = HYPRE_IJMatrixGetObject(A, (void**)&parA); hcheck(ierr, "IJMatrixGetObject", message, ok);

    for (HYPRE_IJVector* v : {&bvec, &xvec}) {
        ierr = HYPRE_IJVectorCreate(comm, HYPRE_BigInt(ilower), HYPRE_BigInt(iupper), v); hcheck(ierr, "IJVectorCreate", message, ok);
        ierr = HYPRE_IJVectorSetObjectType(*v, HYPRE_PARCSR); hcheck(ierr, "IJVectorSetObjectType", message, ok);
        ierr = HYPRE_IJVectorInitialize(*v); hcheck(ierr, "IJVectorInitialize", message, ok);
        ierr = HYPRE_IJVectorAssemble(*v); hcheck(ierr, "IJVectorAssemble", message, ok);
    }
    ierr = HYPRE_IJVectorGetObject(bvec, (void**)&parb); hcheck(ierr, "IJVectorGetObject", message, ok);
    ierr = HYPRE_IJVectorGetObject(xvec, (void**)&parx); hcheck(ierr, "IJVectorGetObject", message, ok);
}

void HypreSystem::Impl::make_solver ()
{
    if (!ok) { return; }
    MPI_Comm comm = ParallelDescriptor::Communicator();
    HypreKrylov k = opt.krylov;
    if (k == HypreKrylov::Auto) { k = (nlev == 1) ? HypreKrylov::PCG : HypreKrylov::GMRES; }
    kind = (k == HypreKrylov::PCG) ? 0 : (k == HypreKrylov::GMRES) ? 1 : 2;
    int ierr = 0;
    ierr = HYPRE_BoomerAMGCreate(&precond); hcheck(ierr, "BoomerAMGCreate", message, ok);
    HYPRE_BoomerAMGSetPrintLevel(precond, 0);
    HYPRE_BoomerAMGSetCoarsenType(precond, opt.coarsen_type);
    HYPRE_BoomerAMGSetRelaxType(precond, opt.relax_type);
    HYPRE_BoomerAMGSetNumSweeps(precond, opt.num_sweeps);
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    if (opt.strong_threshold >= 0.0) { HYPRE_BoomerAMGSetStrongThreshold(precond, opt.strong_threshold); }
    if (opt.interp_type >= 0) { HYPRE_BoomerAMGSetInterpType(precond, opt.interp_type); }
    if (opt.agg_levels > 0) { HYPRE_BoomerAMGSetAggNumLevels(precond, opt.agg_levels); }
    if (kind == 0) {
        ierr = HYPRE_ParCSRPCGCreate(comm, &solver); hcheck(ierr, "PCGCreate", message, ok);
        HYPRE_ParCSRPCGSetMaxIter(solver, 1000);
        HYPRE_ParCSRPCGSetTol(solver, 1.0e-12);
        HYPRE_ParCSRPCGSetTwoNorm(solver, 1);
        HYPRE_ParCSRPCGSetPrintLevel(solver, 0);
        HYPRE_ParCSRPCGSetLogging(solver, 1);
        HYPRE_ParCSRPCGSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
        ierr = HYPRE_ParCSRPCGSetup(solver, parA, parb, parx); hcheck(ierr, "PCGSetup", message, ok);
    } else if (kind == 1) {
        ierr = HYPRE_ParCSRGMRESCreate(comm, &solver); hcheck(ierr, "GMRESCreate", message, ok);
        HYPRE_ParCSRGMRESSetKDim(solver, opt.gmres_kdim);
        HYPRE_ParCSRGMRESSetMaxIter(solver, 1000);
        HYPRE_ParCSRGMRESSetTol(solver, 1.0e-12);
        HYPRE_ParCSRGMRESSetPrintLevel(solver, 0);
        HYPRE_ParCSRGMRESSetLogging(solver, 1);
        HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
        ierr = HYPRE_ParCSRGMRESSetup(solver, parA, parb, parx); hcheck(ierr, "GMRESSetup", message, ok);
    } else {
        ierr = HYPRE_ParCSRBiCGSTABCreate(comm, &solver); hcheck(ierr, "BiCGSTABCreate", message, ok);
        HYPRE_ParCSRBiCGSTABSetMaxIter(solver, 1000);
        HYPRE_ParCSRBiCGSTABSetTol(solver, 1.0e-12);
        HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, 0);
        HYPRE_ParCSRBiCGSTABSetLogging(solver, 1);
        HYPRE_ParCSRBiCGSTABSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
        ierr = HYPRE_ParCSRBiCGSTABSetup(solver, parA, parb, parx); hcheck(ierr, "BiCGSTABSetup", message, ok);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
HypreSystem::HypreSystem (std::vector<HypreLayoutLevel> const& levels, std::array<BC,6> const& ebc, HypreOptions const& opt,
                          bool singular, bool pin_singular)
    : m_impl(std::make_unique<Impl>())
{
    m_impl->lev = levels; m_impl->ebc = ebc; m_impl->opt = opt; m_impl->singular = singular; m_impl->pin_enabled = pin_singular;
    m_impl->build();
}
HypreSystem::~HypreSystem () = default;
bool HypreSystem::ok () const { return m_impl->ok; }
std::string const& HypreSystem::message () const { return m_impl->message; }
long long HypreSystem::rows () const { return m_impl->nglobal; }
long long HypreSystem::nonzeros () const { return m_impl->nnz_global; }
double HypreSystem::setup_seconds () const { return m_impl->setup_s; }
IntVect HypreSystem::pin_cell () const { return m_impl->pin_iv; }
int HypreSystem::pin_level () const { return m_impl->has_pin ? m_impl->pin_lev : -1; }
std::string HypreSystem::method () const { return m_impl->kind == 0 ? "PCG" : m_impl->kind == 1 ? "GMRES" : "BiCGSTAB"; }

bool HypreSystem::matches (std::vector<HypreLayoutLevel> const& levels, std::array<BC,6> const& ebc, HypreOptions const& opt, bool singular) const
{
    Impl const& I = *m_impl;
    if (!(I.ebc == ebc) || !(I.opt == opt) || I.singular != singular || I.lev.size() != levels.size()) { return false; }
    for (std::size_t l = 0; l < levels.size(); ++l) {
        if (!(I.lev[l].ba == levels[l].ba) || !(I.lev[l].dm == levels[l].dm) || !same_geometry(I.lev[l].geom, levels[l].geom)) { return false; }
        if (l > 0 && I.lev[l].ratio != levels[l].ratio) { return false; }
    }
    return true;
}

namespace {
// Copy valid cells of `mf` (scaled) into the rows of an IJ vector.
void put_vector (HypreSystem::Impl const& I, HYPRE_IJVector v, std::vector<MultiFab const*> const& f, std::vector<double> const& scale,
                 bool skip_covered, bool pin_zero, double shift_pin = 0.0)
{
    std::vector<HYPRE_BigInt> idx;
    std::vector<double> val;
    for (int l = 0; l < I.nlev; ++l) {
        if (!f[l]) { continue; }
        for (MFIter mfi(*f[l]); mfi.isValid(); ++mfi) {
            const int b = mfi.index();
            Box const& bx = mfi.validbox();
            auto const& a = f[l]->const_array(mfi);
            std::vector<char> unc(bx.numPts(), 1);
            if (skip_covered && l + 1 < I.nlev) {
                for (auto const& is : I.cfine[l].intersections(bx)) {
                    Box const& ib = is.second;
                    for (int k = ib.smallEnd(2); k <= ib.bigEnd(2); ++k) for (int j = ib.smallEnd(1); j <= ib.bigEnd(1); ++j) for (int i = ib.smallEnd(0); i <= ib.bigEnd(0); ++i) {
                        unc[HypreSystem::Impl::lin(bx, IntVect(i,j,k))] = 0;
                    }
                }
            }
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const Long q = HypreSystem::Impl::lin(bx, IntVect(i,j,k));
                const long long id = I.base[l][b] + q;
                double x = unc[q] ? scale[l]*double(a(i,j,k)) : 0.0;
                if (pin_zero && I.has_pin && id == I.pin_id) { x = 0.0; }
                else if (shift_pin != 0.0 && unc[q]) { x -= shift_pin; }
                idx.push_back(HYPRE_BigInt(id)); val.push_back(x);
            }
        }
    }
    if (!idx.empty()) { HYPRE_IJVectorSetValues(v, HYPRE_Int(idx.size()), idx.data(), val.data()); }
    HYPRE_IJVectorAssemble(v);
}
void get_vector (HypreSystem::Impl const& I, HYPRE_IJVector v, std::vector<MultiFab*> const& f, bool skip_covered)
{
    for (int l = 0; l < I.nlev; ++l) {
        for (MFIter mfi(*f[l]); mfi.isValid(); ++mfi) {
            const int b = mfi.index();
            Box const& bx = mfi.validbox();
            auto const& a = f[l]->array(mfi);
            std::vector<HYPRE_BigInt> idx(bx.numPts());
            std::vector<double> val(bx.numPts());
            for (Long q = 0; q < bx.numPts(); ++q) { idx[q] = HYPRE_BigInt(I.base[l][b] + q); }
            HYPRE_IJVectorGetValues(v, HYPRE_Int(idx.size()), idx.data(), val.data());
            std::vector<char> unc(bx.numPts(), 1);
            if (skip_covered && l + 1 < I.nlev) {
                for (auto const& is : I.cfine[l].intersections(bx)) {
                    Box const& ib = is.second;
                    for (int k = ib.smallEnd(2); k <= ib.bigEnd(2); ++k) for (int j = ib.smallEnd(1); j <= ib.bigEnd(1); ++j) for (int i = ib.smallEnd(0); i <= ib.bigEnd(0); ++i) {
                        unc[HypreSystem::Impl::lin(bx, IntVect(i,j,k))] = 0;
                    }
                }
            }
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const Long q = HypreSystem::Impl::lin(bx, IntVect(i,j,k));
                if (unc[q]) { a(i,j,k) = Real(val[q]); }
            }
        }
    }
}
}

BackendStatus HypreSystem::solve (std::vector<MultiFab*> const& phi, std::vector<MultiFab const*> const& b,
                                  double tol_rel, int max_iter, bool use_guess)
{
    Impl& I = *m_impl;
    BackendStatus s;
    s.hypre_rows = I.nglobal; s.hypre_nnz = nonzeros();
    s.hypre_method = method();
    s.pin_applied = I.has_pin; s.pin_level = I.has_pin ? I.pin_lev : 0; s.pin_cell = I.pin_iv;
    std::vector<double> fscale(I.nlev);
    for (int l = 0; l < I.nlev; ++l) { fscale[l] = -I.vol[l]; }
    put_vector(I, I.bvec, b, fscale, true, true);
    if (use_guess) {
        double pv = 0.0;
        if (I.has_pin) {
            for (MFIter mfi(*phi[I.pin_lev]); mfi.isValid(); ++mfi) {
                if (mfi.validbox().contains(I.pin_iv)) { pv = double(phi[I.pin_lev]->const_array(mfi)(I.pin_iv[0], I.pin_iv[1], I.pin_iv[2])); }
            }
            ParallelDescriptor::ReduceRealSum(pv);
        }
        std::vector<MultiFab const*> pc(phi.begin(), phi.end());
        put_vector(I, I.xvec, pc, std::vector<double>(I.nlev, 1.0), true, true, pv);
    } else {
        std::vector<MultiFab const*> none(I.nlev, nullptr);
        std::vector<HYPRE_BigInt> idx; std::vector<double> val;
        for (long long r = I.ilower; r <= I.iupper; ++r) { idx.push_back(HYPRE_BigInt(r)); val.push_back(0.0); }
        if (!idx.empty()) { HYPRE_IJVectorSetValues(I.xvec, HYPRE_Int(idx.size()), idx.data(), val.data()); }
        HYPRE_IJVectorAssemble(I.xvec);
    }
    const double t0 = amrex::second();
    int ierr = 0;
    HYPRE_Int nit = 0;
    double rel = -1.0;
    // The Krylov recurrence residual is slightly optimistic compared with the independent true residual (pin row, round-off),
    // so the solver is asked for a tenth of the requested tolerance; convergence is still judged against tol_rel.
    const double tol_int = 0.1*tol_rel;
    if (I.kind == 0) {
        HYPRE_ParCSRPCGSetTol(I.solver, tol_int); HYPRE_ParCSRPCGSetMaxIter(I.solver, max_iter);
        ierr = HYPRE_ParCSRPCGSolve(I.solver, I.parA, I.parb, I.parx);
        HYPRE_ParCSRPCGGetNumIterations(I.solver, &nit); HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(I.solver, &rel);
    } else if (I.kind == 1) {
        HYPRE_ParCSRGMRESSetTol(I.solver, tol_int); HYPRE_ParCSRGMRESSetMaxIter(I.solver, max_iter);
        ierr = HYPRE_ParCSRGMRESSolve(I.solver, I.parA, I.parb, I.parx);
        HYPRE_ParCSRGMRESGetNumIterations(I.solver, &nit); HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(I.solver, &rel);
    } else {
        HYPRE_ParCSRBiCGSTABSetTol(I.solver, tol_int); HYPRE_ParCSRBiCGSTABSetMaxIter(I.solver, max_iter);
        ierr = HYPRE_ParCSRBiCGSTABSolve(I.solver, I.parA, I.parb, I.parx);
        HYPRE_ParCSRBiCGSTABGetNumIterations(I.solver, &nit); HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(I.solver, &rel);
    }
    s.hypre_solve_seconds = amrex::second() - t0;
    if (ierr != 0) { HYPRE_ClearAllErrors(); }
    s.iterations = int(nit);
    s.own_residual = rel;
    s.converged = (rel >= 0.0 && rel <= tol_rel*(1.0 + 1.0e-9)) || (ierr == 0 && rel < 0.0);
    get_vector(I, I.xvec, phi, true);
    return s;
}

void HypreSystem::apply (std::vector<MultiFab const*> const& x, std::vector<MultiFab*> const& y) const
{
    Impl const& I = *m_impl;
    MPI_Comm comm = ParallelDescriptor::Communicator();
    HYPRE_IJVector vx, vy;
    HYPRE_IJVectorCreate(comm, HYPRE_BigInt(I.ilower), HYPRE_BigInt(I.iupper), &vx);
    HYPRE_IJVectorSetObjectType(vx, HYPRE_PARCSR); HYPRE_IJVectorInitialize(vx);
    HYPRE_IJVectorCreate(comm, HYPRE_BigInt(I.ilower), HYPRE_BigInt(I.iupper), &vy);
    HYPRE_IJVectorSetObjectType(vy, HYPRE_PARCSR); HYPRE_IJVectorInitialize(vy);
    put_vector(I, vx, x, std::vector<double>(I.nlev, 1.0), false, false);
    std::vector<MultiFab const*> none(I.nlev, nullptr);
    {
        std::vector<HYPRE_BigInt> idx; std::vector<double> val;
        for (long long r = I.ilower; r <= I.iupper; ++r) { idx.push_back(HYPRE_BigInt(r)); val.push_back(0.0); }
        if (!idx.empty()) { HYPRE_IJVectorSetValues(vy, HYPRE_Int(idx.size()), idx.data(), val.data()); }
        HYPRE_IJVectorAssemble(vy);
    }
    HYPRE_ParVector px, py;
    HYPRE_IJVectorGetObject(vx, (void**)&px); HYPRE_IJVectorGetObject(vy, (void**)&py);
    HYPRE_ParCSRMatrixMatvec(1.0, I.parA, px, 0.0, py);
    get_vector(I, vy, y, false);
    HYPRE_IJVectorDestroy(vx); HYPRE_IJVectorDestroy(vy);
}

// ---------------------------------------------------------------------------------------------------------------------
// Single-level PressureBackend adaptor (set-up cached for a repeated layout, as the FFT plan).
namespace {
class HypreBackend final : public PressureBackend {
public:
    const char* name () const override { return "HYPRE"; }
    bool plan_matches (PressureProblem const& p) const override
    {
        if (!m_sys || !m_sys->ok()) { return false; }
        return m_sys->matches(layout(p), effective_bc(p.bc, p.geom.Domain()), m_opt, is_singular(p));
    }
    void prepare (PressureProblem const& p) override { build(p, m_opt); }
    void set_options (HypreOptions const& o) { m_opt = o; }
    BackendStatus solve (PressureProblem const& p, PressureOptions const& o, MultiFab& phi, MultiFab const& rhs) override
    {
        BackendStatus s;
        const bool reused = m_sys && m_sys->ok() && m_sys->matches(layout(p), effective_bc(p.bc, p.geom.Domain()), o.hypre, is_singular(p));
        if (!reused) { build(p, o.hypre); }
        if (!m_sys->ok()) { s.converged = false; s.hypre_method = m_sys->message(); return s; }
        s = m_sys->solve({&phi}, {&rhs}, o.tol_rel, o.max_iter, o.use_initial_guess);
        s.plan_reused = reused;
        s.hypre_setup_seconds = reused ? 0.0 : m_setup;
        return s;
    }
    HypreSystem* system () { return m_sys.get(); }
    double last_setup () const { return m_setup; }
private:
    static std::vector<HypreLayoutLevel> layout (PressureProblem const& p)
    {
        HypreLayoutLevel L; L.ba = p.ba; L.dm = p.dm; L.geom = p.geom;
        return {L};
    }
    static bool is_singular (PressureProblem const& p)
    {
        bool open = false;
        std::array<BC,6> const e = effective_bc(p.bc, p.geom.Domain());
        for (int f = 0; f < 6; ++f) { open = open || (e[f] == BC::Dirichlet); }
        return !open;
    }
    void build (PressureProblem const& p, HypreOptions const& opt)
    {
        m_opt = opt;
        m_sys.reset();
        m_sys = std::make_unique<HypreSystem>(layout(p), effective_bc(p.bc, p.geom.Domain()), opt, is_singular(p));
        m_setup = m_sys->setup_seconds();
    }
    std::unique_ptr<HypreSystem> m_sys;
    HypreOptions m_opt;
    double m_setup = 0.0;
};
}

std::unique_ptr<PressureBackend> make_hypre_backend () { return std::make_unique<HypreBackend>(); }

} // namespace pb

#endif // PB_WITH_HYPRE
