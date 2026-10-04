// main.cpp: C++ entry point of the FDS-AMReX driver (M2a, S1 skeleton).
//
// Order: MPI_Init_thread(FUNNELED) -> amrex::Initialize -> fds_setup(0) [unchanged FDS set-up, main.f90 patch 0001]
//        -> build the level-0 layout -> (S2 onward: fields, time loop) -> amrex::Finalize -> fds_setup(mode) [FDS end, calls
//        MPI_Finalize and STOP, so it must come last].
//
// Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp (S2); (b) only uniform Cartesian metrics are used.
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <mpi.h>

#include <cstdio>
#include <cstdlib>

#include <cstring>
#include <string>

#include "FdsSetup.H"
#include "PressureBcMap.H"
#include "TimeLoop.H"
#include "TwoLevelRun.H"
#ifdef FDSRT_DRIVER_MODES
#include "DriverModes.H"
#endif
namespace fdsamr { int level_bind_check(TimeLoop& loop, const Level0& l0, double dt0); }   // tests/level_bind_check.cpp

extern "C" void fds_setup(int mode, const char* fname, double* dt_out);
extern "C" void fds_p_params(int* ip, double* rp);
extern "C" void fds_p_mesh_info(int nm, int* mi, double* r);
extern "C" void fds_k_visc(int nm, int est);
#ifdef FDS_FINE_B_DRAFT
extern "C" int fds_fine_b_selftest();
extern "C" void fds_fine_b_abort_test(int kind);
extern "C" int fds_fine_b_shadow_enable(int mode);
extern "C" int fds_fine_b_shadow_report();
#endif

namespace {
// --pressure-bc: one line "PBC ..." with the FDS pressure codes of the level's domain faces and the strings the driver passes per direction (PressureBcMap.H), no time loop.
void print_pressure_bc(const fdsamr::Level0& l0, const char* name)
{
    int ip[32] = {0}; double rp[8] = {0};
    fds_p_params(ip, rp);
    int code[6] = {-1, -1, -1, -1, -1, -1};
    const amrex::Box dom = l0.geom.Domain();
    for (int i = 0; i < static_cast<int>(l0.ba.size()); ++i) {
        if (l0.dm[i] != amrex::ParallelDescriptor::MyProc()) continue;
        int mi[16]; double r[4];
        fds_p_mesh_info(i + 1, mi, r);
        const amrex::Box vb = l0.ba[i];
        for (int d = 0; d < 3; ++d) {
            if (vb.smallEnd(d) == dom.smallEnd(d)) code[d] = mi[3 + d];
            if (vb.bigEnd(d) == dom.bigEnd(d)) code[3 + d] = mi[3 + d];
        }
    }
    for (int f = 0; f < 6; ++f) amrex::ParallelAllReduce::Max(code[f], amrex::ParallelContext::CommunicatorSub());
    std::string line = std::string("PBC ") + name + " twod=" + std::to_string(ip[20]) + " n=" + std::to_string(dom.length(0)) + "x" + std::to_string(dom.length(1)) + "x" + std::to_string(dom.length(2)) +
                       " codes=" + std::to_string(code[0]) + "," + std::to_string(code[1]) + "," + std::to_string(code[2]);
    fdsamr::DirBc dm[3];
    for (int d = 0; d < 3; ++d) {
        const bool ign = (d == 1 && ip[20] && dom.length(1) == 1);
        dm[d] = fdsamr::map_pressure_bc_direction(d, code[d], code[3 + d], dom.length(d), l0.dom.periodic[d] != 0, ign);
    }
    if (dm[1].error.empty() && dm[0].error.empty() && dm[2].error.empty() && ip[20] && dom.length(1) == 1) fdsamr::ignored_direction_follows_open(dm[1], dm[0], dm[2]);
    for (int d = 0; d < 3; ++d) {
        const fdsamr::DirBc& m = dm[d];
        line += std::string(" ") + "xyz"[d] + "=" + (m.error.empty() ? fdsamr::bc_string(m) : std::string("ERROR"));
        if (!m.error.empty()) line += " [" + m.error + "]";
        if (!m.note.empty()) line += " [" + m.note + "]";
    }
    amrex::Print() << line << "\n";
}
}  // namespace
int fds_selftest(const fdsamr::Level0& l0);   // tests/selftest_fds.cpp
int fds_kernelcheck(const fdsamr::Level0& l0, const std::string& dump, bool window, const std::string& ghost);   // tests/kernelcheck.cpp

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: fds_amr case.fds\n");
        return 1;
    }

    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    // Exit code 0 only through the Fortran STOP in fds_setup(mode>0).
    {
        // build_parm_parse = false: argv[1] is the FDS input file, not an AMReX inputs file.
        amrex::Initialize(argc, argv, false, MPI_COMM_WORLD);
        int mode_end = 2;  // 2 = FDS setup-only stop message (S1 has no time loop yet)
        {
            double dt = 0.;
            fds_setup(0, argv[1], &dt);

            fdsamr::Level0 l0 = fdsamr::build_level0();
            fdsamr::print_level0(l0);
            amrex::Print() << "FDS-AMReX: initial dt from FDS set-up = " << dt << "\n";
            if (argc > 2 && std::strcmp(argv[2], "--pressure-bc") == 0) print_pressure_bc(l0, argv[1]);
            if (argc > 2 && std::strcmp(argv[2], "--fine-guard-test") == 0) {
                // D-056: a kernel wrapper called with a fine-level mesh number (above NMESHES) must stop the run with a clear message (non-zero exit); this call does not return
                amrex::Print() << "FINE-GUARD-TEST: calling fds_k_visc with mesh number 1000000\n";
                fds_k_visc(1000000, 0);
                amrex::Print() << "FINE-GUARD-TEST FAIL: returned\n";
            }
#ifdef FDS_FINE_B_DRAFT
            if (argc > 2 && std::strcmp(argv[2], "--fine-b-selftest") == 0) {
                const int nfail = fds_fine_b_selftest();
                amrex::Print() << (nfail == 0 ? "FINE-B SELFTEST PASS" : "FINE-B SELFTEST FAIL") << "\n";
            }
            if (argc > 3 && std::strcmp(argv[2], "--fine-b-abort") == 0) {
                fds_fine_b_abort_test(std::atoi(argv[3]));
                amrex::Print() << "FINE-B ABORT TEST: returned\n";
            }
#endif
            if (argc > 2 && std::strcmp(argv[2], "--selftest") == 0) {
                const int nfail = fds_selftest(l0);
                amrex::Print() << (nfail == 0 ? "SELFTEST PASS" : "SELFTEST FAIL") << "\n";
            }
            if (argc > 2 && std::strcmp(argv[2], "--level-bind-check") == 0) {
                // S12: ratio-1 copy of level 0 bound through TimeLoop::bind_level, stage by stage comparison with level 0 (tests/level_bind_check.cpp)
                fdsamr::RunOptions ro;
                ro.quiet = true;
                fdsamr::TimeLoop loop(l0, dt, ro);
                const int nfail = fdsamr::level_bind_check(loop, l0, dt);
                if (nfail != 0) amrex::Abort("level bind check failed");
            }
#ifdef FDSRT_DRIVER_MODES
            if (argc > 2 && std::strcmp(argv[2], "--rt-e2e") == 0) {
                // Role 3 end-to-end modes (regrid_transport/DriverModes.cpp): level 1 over the level-0 case, stage entry points, FluxStageRunner
                fdsamr::RunOptions ro;
                ro.quiet = true;
                fdsamr::TimeLoop loop(l0, dt, ro);
                const int nfail = fdsrt::driver_mode(argc, argv, loop, l0, dt);
                if (nfail != 0) amrex::Abort("rt-e2e failed");
            }
#endif
            if (argc > 2 && std::strcmp(argv[2], "--two-level-run") == 0) {
                // --two-level-run [--steps N] [--patch i0 i1 k0 k1] [--ratio R] [--maxsize M] [--blocking B] [--projection ON|AUTO|OFF] [--no-overwrite] [--outdir D] [--chid C]: S14
                fdsamr::TwoLevelOptions to;
                fdsamr::RunOptions ro;
                ro.quiet = true;
                for (int i = 3; i < argc; ++i) {
                    if (std::strcmp(argv[i], "--steps") == 0 && i + 1 < argc) to.steps = std::atoi(argv[++i]);
                    else if (std::strcmp(argv[i], "--patch") == 0 && i + 4 < argc) { for (int q = 0; q < 4; ++q) to.patch[q] = std::atoi(argv[++i]); }
                    else if (std::strcmp(argv[i], "--ratio") == 0 && i + 1 < argc) to.ratio = std::atoi(argv[++i]);
                    else if (std::strcmp(argv[i], "--maxsize") == 0 && i + 1 < argc) to.maxsize = std::atoi(argv[++i]);
                    else if (std::strcmp(argv[i], "--blocking") == 0 && i + 1 < argc) to.blocking = std::atoi(argv[++i]);
                    else if (std::strcmp(argv[i], "--projection") == 0 && i + 1 < argc) to.projection = argv[++i];
                    else if (std::strcmp(argv[i], "--no-overwrite") == 0) to.overwrite = false;
                    else if (std::strcmp(argv[i], "--outdir") == 0 && i + 1 < argc) { to.outdir = argv[++i]; ro.outdir = to.outdir; }
                    else if (std::strcmp(argv[i], "--chid") == 0 && i + 1 < argc) { to.chid = argv[++i]; ro.chid = to.chid; }
                    else if (std::strcmp(argv[i], "--log-every") == 0 && i + 1 < argc) to.log_every = std::atoi(argv[++i]);
                }
                fdsamr::TimeLoop loop(l0, dt, ro);
                const int nfail = fdsamr::two_level_run(loop, l0, dt, to);
                amrex::Print() << (nfail == 0 ? "TWO-LEVEL RUN PASS" : "TWO-LEVEL RUN FAIL") << "\n";
                if (nfail != 0) amrex::Abort("two-level run failed");
            }
            if (argc > 2 && std::strcmp(argv[2], "--run") == 0) {
                // --run [--steps N] [--outdir D] [--chid C] [--quiet] [--exact-zone-sums]: the C++ time loop (S5)
                fdsamr::RunOptions ro;
                bool fine_shadow = false;
                int fine_build = 0;
                for (int i = 3; i < argc; ++i) {
                    if (std::strcmp(argv[i], "--steps") == 0 && i + 1 < argc) ro.max_steps = std::atoi(argv[++i]);
                    else if (std::strcmp(argv[i], "--outdir") == 0 && i + 1 < argc) ro.outdir = argv[++i];
                    else if (std::strcmp(argv[i], "--chid") == 0 && i + 1 < argc) ro.chid = argv[++i];
                    else if (std::strcmp(argv[i], "--quiet") == 0) ro.quiet = true;
                    else if (std::strcmp(argv[i], "--exact-zone-sums") == 0) ro.exact_zone_sums = true;
                    else if (std::strcmp(argv[i], "--fine-b-shadow") == 0) fine_shadow = true;
                    else if (std::strcmp(argv[i], "--fine-b-build") == 0) { fine_shadow = true; fine_build = 1; }
                }
#ifdef FDS_FINE_B_DRAFT
                // --fine-b-shadow: every kernel call is repeated on a fine-level copy of the box (draft/fds_fine_box_b.f90) and compared bit by bit
                if (fine_shadow && fds_fine_b_shadow_enable(fine_build) != 0) amrex::Abort("--fine-b-shadow needs a single level-0 mesh");
#else
                if (fine_shadow) amrex::Abort("--fine-b-shadow needs a build with -DFDS_AMR_FINE_B_DRAFT=ON (patches 0007 and 0008)");
                (void)fine_build;   // only the draft build reads it
#endif
                fdsamr::TimeLoop loop(l0, dt, ro);
                int nfail = loop.run();
#ifdef FDS_FINE_B_DRAFT
                if (fine_shadow) nfail += fds_fine_b_shadow_report();
#endif
                amrex::Print() << (nfail == 0 ? "RUN COMPLETE" : "RUN FAILED") << ": steps = " << loop.icyc() << ", T = " << loop.time() << "\n";
                if (nfail == 0) mode_end = 1;
            }
            if (argc > 3 && std::strcmp(argv[2], "--kernelcheck") == 0) {
                bool window = false;
                std::string ghost = "dump";
                for (int i = 4; i < argc; ++i) {
                    if (std::strcmp(argv[i], "--window") == 0) window = true;
                    else if (std::strncmp(argv[i], "--ghost=", 8) == 0) ghost = argv[i] + 8;
                }
                const int nfail = fds_kernelcheck(l0, argv[3], window, ghost);
                amrex::Print() << (nfail == 0 ? "KERNELCHECK PASS" : "KERNELCHECK FAIL") << "\n";
            }
        }
        amrex::Finalize();
        fds_setup(mode_end, argv[1], nullptr);  // does not return
    }
    return 0;
}
