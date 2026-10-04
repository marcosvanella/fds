// Checks of a pressure backend built WITHOUT PB_WITH_HYPRE (the driver-style source list; see README.md, "Build options"): the HYPRE
// selector answers Status::NotBuilt with kHypreNotBuilt for a single level and for a composite request, a workspace never holds a HYPRE
// set-up, and the FFT and MLMG paths still solve. Linked with the sources compiled without the definition (object library pb_nohypre_obj).
#include "PressureIface.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <cmath>
#include <cstdio>
#include <string>

using namespace amrex;

namespace {
int g_fail = 0;
void check (bool ok, const char* what) { std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what); if (!ok) { ++g_fail; } }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        const int n = 16;
        Box dom(IntVect(0), IntVect(n - 1));
        RealBox rb({0., 0., 0.}, {1., 1., 1.});
        Geometry geom(dom, rb, 0, {1, 1, 1});
        BoxArray ba(dom); ba.maxSize(8);
        DistributionMapping dm(ba);
        MultiFab rhs(ba, dm, 1, 0), phi(ba, dm, 1, 1);
        const double tp = 8.0 * std::atan(1.0);
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            auto const& a = rhs.array(mfi);
            Box b = mfi.validbox();
            for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k)
                for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j)
                    for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i)
                        a(i, j, k) = std::sin(tp * (i + 0.5) / n) * std::cos(tp * (j + 0.5) / n) + 0.25 * std::sin(tp * (k + 0.5) / n);
        }
        phi.setVal(0.0);

        pb::PressureProblem p;
        p.ba = ba; p.dm = dm; p.geom = geom; p.rhs = &rhs; p.phi = &phi;
        for (auto& c : p.bc) { c = pb::BC::Periodic; }

        pb::PressureOptions o;
        o.backend = pb::BackendKind::HYPRE;
        pb::PressureResult r = pb::solve_pressure(p, o, nullptr);
        check(r.status == pb::Status::NotBuilt && r.message == pb::kHypreNotBuilt, "single level: HYPRE request is NotBuilt('built without HYPRE')");
        check(std::string(pb::kHypreNotBuilt) == "built without HYPRE", "message text");

        pb::PressureWorkspace ws;
        r = pb::solve_pressure(p, o, &ws);
        check(r.status == pb::Status::NotBuilt && !ws.hypre_built() && !ws.built(), "single level with a workspace: NotBuilt, no set-up cached");

        pb::PressureOptions of; of.backend = pb::BackendKind::FFT;
        r = pb::solve_pressure(p, of, &ws);
        check(r.status == pb::Status::Ok && r.backend == "FFT" && ws.built() && !ws.hypre_built(), "single level FFT still solves, workspace holds the FFT plan only");

        pb::PressureOptions om; om.backend = pb::BackendKind::MLMG;
        phi.setVal(0.0);
        r = pb::solve_pressure(p, om, nullptr);
        check(r.status == pb::Status::Ok && r.backend == "MLMG", "single level MLMG still solves");

        // composite request (one level is enough: the selector answers before the hierarchy is validated)
        pb::PressureProblem pc;
        pb::PressureLevel L; L.ba = ba; L.dm = dm; L.geom = geom; L.rhs = &rhs; L.phi = &phi;
        pc.levels.push_back(L);
        pc.bc = p.bc;
        r = pb::solve_pressure(pc, o, nullptr);
        check(r.status == pb::Status::NotBuilt && r.message == pb::kHypreNotBuilt, "composite: HYPRE request is NotBuilt('built without HYPRE')");
        pb::Selection s = pb::select_backend(pc, pb::BackendKind::HYPRE);
        check(!s.ok && s.message == pb::kHypreNotBuilt, "select_backend(HYPRE) rejects with the same message");
        phi.setVal(0.0);
        r = pb::solve_pressure(pc, om, nullptr);
        check(r.status == pb::Status::Ok && r.backend == "MLMG", "composite (one level) MLMG still solves");
    }
    amrex::Finalize();
    std::printf("SUMMARY failures=%d %s\n", g_fail, g_fail ? "FAIL" : "PASS");
    return g_fail ? 1 : 0;
}
