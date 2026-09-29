// FdsSetup.cpp: level-0 layout and cell wall data taken from the FDS set-up (calls the Fortran bind(C) queries of
// fds_mesh_query.f90). Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp; (b) only uniform Cartesian metrics
// are used, nonuniform (TRN*) and cylindrical meshes are rejected.
#include "FdsSetup.H"

#include <AMReX_Print.H>

extern "C" {
int fds_get_nmeshes();
void fds_get_mesh(int nm, int* ijk, double* xb, int* rank, int* nonuniform);
void fds_get_domain(int* periodic, int* cyl, int* n_tracked, int* n_total, int* nranks);
void fds_get_cell_walls(int nm, const int* n, int* flags);
}

namespace fdsamr {

Level0 build_level0()
{
    DomainInfo dom;
    fds_get_domain(dom.periodic, &dom.cylindrical, &dom.n_tracked, &dom.n_total, &dom.nranks);
    const int nm = fds_get_nmeshes();
    amrex::Vector<MeshInfo> meshes(nm);
    for (int i = 0; i < nm; ++i) {
        int nonuniform = 0;
        MeshInfo& m = meshes[i];
        fds_get_mesh(i + 1, m.ijk, m.xb, &m.rank, &nonuniform);
        if (nonuniform) amrex::Abort("M2a: TRNX/TRNY/TRNZ (nonuniform) meshes are not supported (IR-002)");
    }
    return assemble_level0(meshes, dom);
}

// Fills comps 0..6 of the valid cells of the FDS mesh nm (= box index + 1) from MESHES(nm)%CELL (local meshes only).
void fds_cell_walls(int nm, const amrex::Box& vbox, amrex::Array4<int> const& a)
{
    const int n[3] = {vbox.length(0), vbox.length(1), vbox.length(2)};
    std::vector<int> f(static_cast<std::size_t>(7) * n[0] * n[1] * n[2]);
    fds_get_cell_walls(nm, n, f.data());
    const amrex::IntVect lo = vbox.smallEnd();
    for (int k = 0; k < n[2]; ++k)
        for (int j = 0; j < n[1]; ++j)
            for (int i = 0; i < n[0]; ++i)
                for (int c = 0; c < 7; ++c)
                    a(lo[0] + i, lo[1] + j, lo[2] + k, c) = f[c + 7 * (i + n[0] * (j + static_cast<std::size_t>(n[1]) * k))];
}

}  // namespace fdsamr
