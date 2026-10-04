// GhostShare.cpp: see GhostShare.H.
#include "GhostShare.H"

#include <AMReX_BoxList.H>
#include <AMReX_Loop.H>

#include <sstream>

namespace fdsrt {

std::string GhostShareCount::to_string() const
{
    std::ostringstream os;
    os << "covered coarse cells serving one ghost face: " << serving_one << ", serving two or more (D-059 shared cells): " << serving_several
       << " (of which on two or more different faces: " << serving_two_faces << ")";
    return os.str();
}

GhostShareCount count_shared_ghost_cells(const amrex::BoxArray& fine_ba, const amrex::IntVect& ratio, const amrex::Geometry& cg)
{
    GhostShareCount out;
    if (fine_ba.empty()) return out;
    amrex::BoxArray cfba = fine_ba;
    cfba.coarsen(ratio);
    const amrex::Box dom = cg.Domain();
    auto covered = [&](amrex::IntVect p) {   // 1 covered, 0 not, -1 outside a non-periodic domain edge
        for (int e = 0; e < 3; ++e) {
            if (p[e] < dom.smallEnd(e)) { if (!cg.isPeriodic(e)) return -1; p[e] += dom.length(e); }
            else if (p[e] > dom.bigEnd(e)) { if (!cg.isPeriodic(e)) return -1; p[e] -= dom.length(e); }
        }
        return cfba.contains(p) ? 1 : 0;
    };
    for (int ib = 0; ib < static_cast<int>(cfba.size()); ++ib) {
        const amrex::Box b = cfba[ib];
        amrex::BoxList shell = amrex::boxDiff(b, amrex::grow(b, -2));   // the cells within two of a box face; the others cannot serve a face
        for (const amrex::Box& sb : shell) {
            amrex::LoopOnCpu(sb, [&](int i, int j, int k) {
                const amrex::IntVect x(i, j, k);
                int n = 0;
                bool face_used[3][2] = {{false, false}, {false, false}, {false, false}};
                for (int layer = 1; layer <= 2; ++layer)
                    for (int d = 0; d < 3; ++d)
                        for (int s = 0; s < 2; ++s) {
                            const int nb = (s == 0) ? -1 : 1;
                            amrex::IntVect y1 = x, y2 = x;
                            y1[d] += nb; y2[d] += 2 * nb;
                            const bool hit = (layer == 1) ? (covered(y1) == 0) : (covered(y1) == 1 && covered(y2) == 0);
                            if (hit) { ++n; face_used[d][s] = true; }
                        }
                if (n == 0) return;
                if (n == 1) { ++out.serving_one; return; }
                ++out.serving_several;
                int nf = 0;
                for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) if (face_used[d][s]) ++nf;
                if (nf >= 2) ++out.serving_two_faces;
            });
        }
    }
    return out;
}

}  // namespace fdsrt
