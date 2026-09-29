// selftest_fds.cpp: driver checks that need the FDS objects (run as `fds_amr case.fds --selftest` after the FDS set-up).
//
// Kernel-facing rules (M2a): (a) passive scalars: ZZ/ZZS are bound with ncomp = N_TOTAL_SCALARS as the fourth extent; (b) only uniform
// Cartesian metrics are used.
//
// Checks (per local box, so that 4 ranks each test their own mesh):
//   bounds   the FDS-allocated arrays of MESHES(NM) have exactly the bounds of the field table (fds_window): the table matches init.f90.
//   alias    binding a MultiFab FAB to MESHES(NM)%X (W2 alias, fds_alias.c) makes Fortran see FDS-indexed data: element (I,J,K[,N]) read from
//            Fortran equals the value written at the AMReX index to_amrex(I,J,K); writes from Fortran land in the FAB.
//   pass     the aliased array passed to an explicit-shape dummy (how the kernels use it) hits the right element.
//   strided  INFO: whether a strided window of a larger FAB (RHO/RHOS ng=3 seen with the FDS ng=2 bounds) is honoured by gfortran.
//   side     SideData built from the FDS CELL data (fds_get_cell_walls): counts and the layout-independent hash (compared across
//            runs by run_driver_tests.sh: shunn3_32 on 1 rank vs shunn3_4mesh_32 on 4 ranks must print the same hash).
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "FdsSetup.H"
#include "Fields.H"
#include "SideData.H"
#include "check.H"

extern "C" {
int fds_shim_bind(int nm, const char* name, void* base, const int* lb, const int* ext, const long* stride);
int fds_shim_release(int nm, const char* name);
int fds_shim_bounds(int nm, const char* name, int* lb, int* ub, int* rnk, int* alloc);
int fds_shim_access(int nm, const char* name, int i, int j, int k, int n, double* val, int doset);
int fds_shim_pass_test(int nm, const char* name, double mark, void** addr);
}

namespace {
using namespace fdsamr;

double val(int fid, int a, int b, int c, int n)
{
    std::uint64_t x = (std::uint64_t)(fid + 1) * 1000003ULL ^ ((std::uint64_t)(a + 1000) << 32) ^ ((std::uint64_t)(b + 1000) << 16) ^
                      (std::uint64_t)(c + 1000) ^ ((std::uint64_t)n << 50);
    x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; x ^= x >> 31;
    return static_cast<double>(x & 0xFFFFFFFFFFFFULL);
}

int fid_of(const std::string& n)
{
    int i = 0;
    for (const auto& s : field_table()) { if (n == s.name) return i; ++i; }
    return -1;
}
}  // namespace

int fds_selftest(const Level0& l0)
{
    using namespace fdstest;
    const int ns = l0.dom.n_total;
    Fields F(l0, ns);
    const int me = amrex::ParallelDescriptor::MyProc();
    long n_bounds_checked = 0, n_bounds_skipped = 0;

    for (const auto& s : field_table()) {
        if (!F.has(s.name)) continue;
        amrex::MultiFab& mf = F[s.name];
        const int fid = fid_of(s.name);
        const bool four = s.per_scalar;
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const int nm = mfi.index() + 1;
            if (l0.dm[mfi.index()] != me) continue;
            const amrex::Box cb = l0.ba[mfi.index()];
            amrex::FArrayBox& fab = mf[mfi];

            // ---- bounds: the FDS allocation vs the table ----
            {
                int lb[4], ub[4], rnk, alloc;
                const int ierr = fds_shim_bounds(nm, s.name, lb, ub, &rnk, &alloc);
                CHECK_MSG(ierr == 0, s.name);
                if (alloc) {
                    const FdsBounds w = fds_window(s, cb);
                    CHECK_MSG(rnk == (four ? 4 : 3), s.name);
                    for (int d = 0; d < 3; ++d) CHECK_MSG(lb[d] == w.lb[d] && ub[d] == w.lb[d] + w.ext[d] - 1, std::string(s.name) + " bounds differ from init.f90");
                    if (four) CHECK_MSG(ub[3] - lb[3] + 1 == ns, s.name);
                    ++n_bounds_checked;
                } else {
                    ++n_bounds_skipped;   // FDS does not allocate it in this case (option off)
                }
            }

            // ---- alias: whole native FAB, contiguous ----
            {
                auto a = fab.array();
                for (int n = 0; n < fab.nComp(); ++n)
                    amrex::LoopOnCpu(fab.box(), [&](int i, int j, int k) { a(i, j, k, n) = val(fid, i, j, k, n); });
                const FdsBounds nb = fds_bounds(s, cb);
                int lb[4] = {nb.lb[0], nb.lb[1], nb.lb[2], 1};
                int ext[4] = {nb.ext[0], nb.ext[1], nb.ext[2], fab.nComp()};
                long stride[4] = {1, nb.ext[0], (long)nb.ext[0] * nb.ext[1], (long)nb.ext[0] * nb.ext[1] * nb.ext[2]};
                CHECK_MSG(fds_shim_bind(nm, s.name, fab.dataPtr(), lb, ext, stride) == 0, s.name);
                int rb[4], ru[4], rr, ra;
                fds_shim_bounds(nm, s.name, rb, ru, &rr, &ra);
                CHECK_MSG(ra == 1, s.name);
                for (int d = 0; d < 3; ++d) CHECK_MSG(rb[d] == nb.lb[d] && ru[d] == nb.lb[d] + nb.ext[d] - 1, std::string(s.name) + " aliased bounds");
                bool ok = true;
                for (int n = 1; n <= fab.nComp() && ok; ++n)
                    for (int K = nb.lb[2]; K < nb.lb[2] + nb.ext[2]; ++K)
                        for (int J = nb.lb[1]; J < nb.lb[1] + nb.ext[1]; ++J)
                            for (int I = nb.lb[0]; I < nb.lb[0] + nb.ext[0]; ++I) {
                                double v = 0;
                                fds_shim_access(nm, s.name, I, J, K, n, &v, 0);
                                const int ai = to_amrex(s, 0, cb.smallEnd(0), I), aj = to_amrex(s, 1, cb.smallEnd(1), J), ak = to_amrex(s, 2, cb.smallEnd(2), K);
                                if (v != val(fid, ai, aj, ak, n - 1)) { ok = false; break; }
                            }
                CHECK_MSG(ok, std::string(s.name) + " Fortran read through alias");
                // write from Fortran, read from the FAB (a valid element, a low ghost and a high ghost)
                const int probes[3][3] = {{1, 1, 1}, {nb.lb[0], nb.lb[1], nb.lb[2]}, {nb.lb[0] + nb.ext[0] - 1, nb.lb[1] + nb.ext[1] - 1, nb.lb[2] + nb.ext[2] - 1}};
                for (int p = 0; p < 3; ++p) {
                    double w = -12345.0 - p;
                    fds_shim_access(nm, s.name, probes[p][0], probes[p][1], probes[p][2], 1, &w, 1);
                    CHECK_MSG(fab(amrex::IntVect(to_amrex(s, 0, cb.smallEnd(0), probes[p][0]), to_amrex(s, 1, cb.smallEnd(1), probes[p][1]),
                                                 to_amrex(s, 2, cb.smallEnd(2), probes[p][2])), 0) == w, std::string(s.name) + " Fortran write through alias");
                }
                // ---- pass: explicit-shape dummy ----
                if (!four && (std::strcmp(s.name, "RHO") == 0 || std::strcmp(s.name, "RHOS") == 0 || std::strcmp(s.name, "TMP") == 0 ||
                              std::strcmp(s.name, "U") == 0 || std::strcmp(s.name, "H") == 0)) {
                    void* addr = nullptr;
                    const double mark = 7.25e9 + fid;
                    CHECK_MSG(fds_shim_pass_test(nm, s.name, mark, &addr) == 0, s.name);
                    const amrex::IntVect p(to_amrex(s, 0, cb.smallEnd(0), nb.lb[0] + 1), to_amrex(s, 1, cb.smallEnd(1), nb.lb[1] + 1),
                                           to_amrex(s, 2, cb.smallEnd(2), nb.lb[2] + 1));
                    CHECK_MSG(fab(p, 0) == mark, std::string(s.name) + " explicit-shape dummy hits the aliased element");
                    CHECK_MSG(addr == static_cast<void*>(&fab(p, 0)), std::string(s.name) + " no copy-in/copy-out for the explicit-shape dummy");
                }
                CHECK_MSG(fds_shim_release(nm, s.name) == 0, s.name);
                fds_shim_bounds(nm, s.name, rb, ru, &rr, &ra);
                CHECK_MSG(ra == 0, std::string(s.name) + " released");
            }

            // ---- strided window (INFO, RHO/RHOS only) ----
            if (s.ng != s.ng_fds && (std::strcmp(s.name, "RHO") == 0)) {
                FdsView w = make_fds_view(s, fab, cb, true);
                const FdsBounds fw = fds_window(s, cb);
                int lb[3] = {fw.lb[0], fw.lb[1], fw.lb[2]};
                int ext[3] = {fw.ext[0], fw.ext[1], fw.ext[2]};
                long stride[3] = {w.stride[0], w.stride[1], w.stride[2]};
                fab.setVal(0.0);
                fds_shim_bind(nm, s.name, w.base, lb, ext, stride);
                double v = 0, mark = 4242.0;
                // write element (1,1,1) and (1,2,1) from Fortran; in a strided window (1,2,1) is w.stride[1] elements further
                fds_shim_access(nm, s.name, 1, 1, 1, 1, &mark, 1);
                mark = 4343.0;
                fds_shim_access(nm, s.name, 1, 2, 1, 1, &mark, 1);
                fds_shim_access(nm, s.name, 1, 2, 1, 1, &v, 0);
                const bool honoured = (w(1, 2, 1) == 4343.0) && (w(1, 1, 1) == 4242.0);
                // explicit-shape dummy on the strided window: does gfortran pass the address or make a contiguous temporary?
                void* saddr = nullptr;
                fds_shim_pass_test(nm, s.name, 5151.0, &saddr);
                const bool no_copy = (saddr == static_cast<void*>(&w(0, 0, 0)));
                const bool landed = (w(0, 0, 0) == 5151.0);
                fds_shim_release(nm, s.name);
                if (nm == 1 && amrex::ParallelDescriptor::MyProc() == 0)
                    std::printf("  INFO strided window to explicit-shape dummy: %s, value %s\n", no_copy ? "no copy" : "COPY made by gfortran",
                                landed ? "reaches the FAB" : "LOST (temporary not copied back)");
                if (amrex::ParallelDescriptor::IOProcessor() && amrex::ParallelDescriptor::MyProc() == 0 && nm == 1)
                    std::printf("  INFO strided alias (RHO ng=3 FAB seen with FDS ng=2 bounds): element access %s by gfortran 14\n",
                                honoured ? "HONOURED" : "NOT honoured (use the full-FAB alias with native bounds instead)");
                // Element access through the descriptor is expected to be honoured; passing to an explicit-shape dummy is NOT (checked in S3).
            }
        }
    }

    // ---- side data from the FDS cell data ----
    {
        SideData sd(l0, fds_cell_walls);
        long c[3];
        sd.counts(c);
        const std::uint64_t h = sd.hash();
        CHECK(c[0] == 0);   // no OBST in the M2a cases
        if (amrex::ParallelDescriptor::IOProcessor())
            std::printf("  SIDEDATA_HASH %016llx solid=%ld wallfaces=%ld source=%ld\n", (unsigned long long)h, c[0], c[1], c[2]);
    }
    if (amrex::ParallelDescriptor::IOProcessor())
        std::printf("  bounds cross-check against init.f90 allocations: %ld array(s) checked on this rank, %ld not allocated in this case\n",
                    n_bounds_checked, n_bounds_skipped);
    const long f = report("selftest_fds (FDS bounds table, alias read/write/pass, side data from FDS cells)");
    return static_cast<int>(f);
}
