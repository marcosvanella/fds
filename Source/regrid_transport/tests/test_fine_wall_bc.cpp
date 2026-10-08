// test_fine_wall_bc.cpp: wall BC of the domain faces of a level > 0 (GhostExchange mirror_domain_edges with a FaceBc). Level 0 8^3, level 1 ratio 2 (16^3), periodic x/y,
// z-lo INERT (kFaceInert), z-hi MIRROR (kFaceMirror). Fine box A (0..7, 0..7, 0..15) touches both z walls, fine box B (8..15, 8..15, 4..11) touches none.
// Asserts: W = 0 on the two wall faces of A (nodal k = 0 and k = 16, FDS W(K=0), W(K=KBAR)); U, V z-ghosts of A = -interior at z-lo (no-slip), +interior at z-hi (free-slip);
// a scalar keeps the even copy; every other value of A (x/y ghosts included) and all of B unchanged. Negative control: the default FaceBc{} leaves W on the wall faces nonzero.
// Exit code 0 = pass.
#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <cmath>
#include <string>

#include "FdsAmr.H"
#include "GhostExchange.H"
#include "check.H"

using namespace fdsamr;

namespace {
double v0(int i, int j, int k, int c) { return 1.0 + c + 0.01 * (i + 3) + 0.1 * (j + 3) + 1.0e-3 * (k + 3) * (k + 5); }   // nonzero everywhere, ghosts included

amrex::MultiFab make(const Level& lev, int nodal_dir, int c)
{
    amrex::IntVect t(0);
    if (nodal_dir >= 0) t[nodal_dir] = 1;
    amrex::MultiFab mf(amrex::convert(lev.ba, t), lev.dm, 1, 1);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf.array(mfi);
        amrex::LoopOnCpu(mf[mfi].box(), [&](int i, int j, int k) { a(i, j, k) = v0(i, j, k, c); });
    }
    return mf;
}

// lo, hi: expected z-ghost of box A = lo (hi) * the cell next to the z-lo (z-hi) wall; wall: W of box A is 0 on the wall faces. Everything else keeps its initial value.
void check_field(const amrex::MultiFab& mf, int c, double lo, double hi, bool wall)
{
    const bool w = mf.ixType().nodeCentered(2);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf.const_array(mfi);
        const bool boxA = mfi.validbox().smallEnd(2) == 0;
        amrex::LoopOnCpu(mf[mfi].box(), [&](int i, int j, int k) {
            double e = v0(i, j, k, c);
            if (boxA && w && wall && (k == 0 || k == 16)) e = 0.0;
            if (boxA && !w && k == -1) e = lo * v0(i, j, 0, c);
            if (boxA && !w && k == 16) e = hi * v0(i, j, 15, c);
            CHECK_MSG(a(i, j, k) == e, "c=" + std::to_string(c) + " box " + (boxA ? "A" : "B") + " (" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ") got " +
                                          std::to_string(a(i, j, k)) + " want " + std::to_string(e));
        });
    }
}
}  // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        const amrex::Array<int, 3> per = {1, 1, 0};
        DomainInfo dom{{1, 1, 0}, 0, 1, 1, amrex::ParallelDescriptor::NProcs()};
        const amrex::RealBox rb({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
        const amrex::Geometry g0(amrex::Box(amrex::IntVect(0), amrex::IntVect(7)), rb, 0, per);
        const amrex::Geometry g1(amrex::refine(g0.Domain(), 2), rb, 0, per);
        amrex::BoxList bl;
        bl.push_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 7, 15)));     // A: both z walls
        bl.push_back(amrex::Box(amrex::IntVect(8, 8, 4), amrex::IntVect(15, 15, 11)));   // B: no wall
        const amrex::BoxArray ba(bl);
        const amrex::DistributionMapping dm(ba);
        const Level l1 = make_layout_level(1, g1, ba, dm, amrex::IntVect(2), dom);
        const FaceBc bc = {0, 0, 0, 0, kFaceInert, kFaceMirror};

        amrex::MultiFab U = make(l1, 0, 0), V = make(l1, 1, 1), W = make(l1, 2, 2), R = make(l1, -1, 3), Wd = make(l1, 2, 2);
        mirror_domain_edges(l1, U, 0, bc);
        mirror_domain_edges(l1, V, 1, bc);
        mirror_domain_edges(l1, W, 2, bc);
        mirror_domain_edges(l1, R, -1);    // scalar: default (even copy)
        mirror_domain_edges(l1, Wd, 2);    // negative control: the old default leaves the wall-normal face alone
        check_field(U, 0, -1.0, 1.0, true);
        check_field(V, 1, -1.0, 1.0, true);
        check_field(W, 2, 1.0, 1.0, true);
        check_field(R, 3, 1.0, 1.0, true);
        check_field(Wd, 2, 1.0, 1.0, false);
        long bad = fdstest::report("fine_wall_bc");
        // the negative control must differ from the fixed W on the wall faces (W nonzero there)
        double wmin = 1.0e30;
        for (amrex::MFIter mfi(Wd); mfi.isValid(); ++mfi) {
            if (mfi.validbox().smallEnd(2) != 0) continue;
            auto a = Wd.const_array(mfi);
            for (int k : {0, 16}) amrex::LoopOnCpu(amrex::Box(amrex::IntVect(0, 0, k), amrex::IntVect(7, 7, k)), [&](int i, int j, int kk) { wmin = std::min(wmin, std::abs(a(i, j, kk))); });
        }
        amrex::ParallelDescriptor::ReduceRealMin(wmin);
        const bool neg = wmin > 0.0;
        amrex::Print() << (neg ? "PASS" : "FAIL") << " fine_wall_bc negative control: default FaceBc leaves min|W| on the z wall faces = " << wmin << " (nonzero)\n";
        bad += neg ? 0 : 1;
        amrex::Print() << (bad == 0 ? "TEST PASS" : "TEST FAIL") << " test_fine_wall_bc\n";
        if (bad) { amrex::Finalize(); return 1; }
    }
    amrex::Finalize();
    return 0;
}
