// FdsAmr.cpp: see FdsAmr.H. Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp (S2); (b) only uniform
// Cartesian metrics are used, R(I)/RRN(I) are dropped and CYLINDRICAL/TRN* meshes are rejected. No Fortran symbol is referenced here,
// so the unit tests link this file without the FDS objects (the FDS queries are in FdsSetup.cpp).
#include "FdsAmr.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <algorithm>
#include <cmath>

namespace fdsamr {

Level0 assemble_level0(const amrex::Vector<MeshInfo>& meshes, const DomainInfo& dom)
{
    Level0 l0;
    l0.dom = dom;
    l0.mesh = meshes;
    if (l0.dom.cylindrical) amrex::Abort("M2a: CYLINDRICAL meshes are not supported (IR-002)");
    const int nm = static_cast<int>(l0.mesh.size());
    AMREX_ALWAYS_ASSERT(nm >= 1);

    // Uniform cell size, equal in all meshes (level 0 has one geometry).
    for (int d = 0; d < 3; ++d) {
        l0.dx[d] = (l0.mesh[0].xb[2*d+1] - l0.mesh[0].xb[2*d]) / l0.mesh[0].ijk[d];
        for (int i = 1; i < nm; ++i) {
            const double dxi = (l0.mesh[i].xb[2*d+1] - l0.mesh[i].xb[2*d]) / l0.mesh[i].ijk[d];
            if (std::abs(dxi - l0.dx[d]) > 1.e-9 * l0.dx[d])
                amrex::Abort("M2a: all meshes must have the same cell size (no mesh refinement at level 0)");
        }
    }

    // Global index space: origin at the lowest mesh corner.
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = l0.mesh[0].xb[2*d];
        hi[d] = l0.mesh[0].xb[2*d+1];
        for (int i = 1; i < nm; ++i) {
            lo[d] = std::min(lo[d], l0.mesh[i].xb[2*d]);
            hi[d] = std::max(hi[d], l0.mesh[i].xb[2*d+1]);
        }
    }

    amrex::BoxList bl;
    amrex::Vector<int> pmap(nm);
    amrex::IntVect dom_hi(0);
    for (int i = 0; i < nm; ++i) {
        amrex::IntVect blo, bhi;
        for (int d = 0; d < 3; ++d) {
            const double s = (l0.mesh[i].xb[2*d] - lo[d]) / l0.dx[d];
            const long is = std::lround(s);
            if (std::abs(s - is) > 1.e-6) amrex::Abort("M2a: mesh corner is not on the level-0 cell lattice");
            blo[d] = static_cast<int>(is);
            bhi[d] = blo[d] + l0.mesh[i].ijk[d] - 1;
            dom_hi[d] = std::max(dom_hi[d], bhi[d]);
        }
        bl.push_back(amrex::Box(blo, bhi));
        pmap[i] = l0.mesh[i].rank;
    }
    l0.ba.define(bl);
    if (!l0.ba.isDisjoint()) amrex::Abort("M2a: FDS meshes overlap in index space");
    if (*std::max_element(pmap.begin(), pmap.end()) >= amrex::ParallelDescriptor::NProcs())
        amrex::Abort("FDS PROCESS(NM) exceeds the AMReX rank count");
    l0.dm.define(pmap);

    const amrex::Box domain(amrex::IntVect(0), dom_hi);
    const amrex::RealBox rb({lo[0], lo[1], lo[2]}, {hi[0], hi[1], hi[2]});
    const amrex::Array<int, 3> is_per{l0.dom.periodic[0], l0.dom.periodic[1], l0.dom.periodic[2]};
    l0.geom.define(domain, rb, 0, is_per);
    return l0;
}

void print_level0(const Level0& l0)
{
    amrex::Print() << "FDS-AMReX level 0: " << l0.ba.size() << " box(es), domain " << l0.geom.Domain()
                   << ", dx = (" << l0.dx[0] << ", " << l0.dx[1] << ", " << l0.dx[2] << "), periodic = ("
                   << l0.dom.periodic[0] << l0.dom.periodic[1] << l0.dom.periodic[2] << "), ranks = " << l0.dom.nranks
                   << ", tracked/total scalars = " << l0.dom.n_tracked << "/" << l0.dom.n_total << "\n";
    for (int i = 0; i < static_cast<int>(l0.ba.size()); ++i) {
        amrex::Print() << "  box " << i << " (FDS mesh " << i + 1 << ", rank " << l0.dm[i] << "): " << l0.ba[i]
                       << " cells " << l0.ba[i].numPts() << "\n";
    }
}

}  // namespace fdsamr
