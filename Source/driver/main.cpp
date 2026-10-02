// main.cpp: C++ entry point of the FDS-AMReX driver (M2a, S1 skeleton).
//
// Order: MPI_Init_thread(FUNNELED) -> amrex::Initialize -> fds_setup(0) [unchanged FDS set-up, main.f90 patch 0001]
//        -> build the level-0 layout -> (S2 onward: fields, time loop) -> amrex::Finalize -> fds_setup(mode) [FDS end, calls
//        MPI_Finalize and STOP, so it must come last].
//
// Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp (S2); (b) only uniform Cartesian metrics are used.
#include <AMReX.H>
#include <AMReX_Print.H>

#include <mpi.h>

#include <cstdio>
#include <cstdlib>

#include <cstring>
#include <string>

#include "FdsSetup.H"
#include "TimeLoop.H"

extern "C" void fds_setup(int mode, const char* fname, double* dt_out);
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
            if (argc > 2 && std::strcmp(argv[2], "--selftest") == 0) {
                const int nfail = fds_selftest(l0);
                amrex::Print() << (nfail == 0 ? "SELFTEST PASS" : "SELFTEST FAIL") << "\n";
            }
            if (argc > 2 && std::strcmp(argv[2], "--run") == 0) {
                // --run [--steps N] [--outdir D] [--chid C] [--quiet] [--exact-zone-sums]: the C++ time loop (S5)
                fdsamr::RunOptions ro;
                for (int i = 3; i < argc; ++i) {
                    if (std::strcmp(argv[i], "--steps") == 0 && i + 1 < argc) ro.max_steps = std::atoi(argv[++i]);
                    else if (std::strcmp(argv[i], "--outdir") == 0 && i + 1 < argc) ro.outdir = argv[++i];
                    else if (std::strcmp(argv[i], "--chid") == 0 && i + 1 < argc) ro.chid = argv[++i];
                    else if (std::strcmp(argv[i], "--quiet") == 0) ro.quiet = true;
                    else if (std::strcmp(argv[i], "--exact-zone-sums") == 0) ro.exact_zone_sums = true;
                }
                fdsamr::TimeLoop loop(l0, dt, ro);
                const int nfail = loop.run();
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
