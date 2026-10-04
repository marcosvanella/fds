// LevelOps.H: level-to-level operators on cell-centred MultiFabs (R2a). AMReX only, no FDS and no driver types.
//
//  average_down_cells : fine -> coarse volume average on the coarse cells covered by the fine level (conservative: the integral over the
//                       covered region is unchanged up to round-off). For a CONSERVED quantity only (rho, rho*Z, ...): never for mass fractions Z.
//  average_down_species: the restriction to use for rho and Z (rho and rho*Z averaged, Z rebuilt).
//  fill_cf_ghosts_pc  : ghost cells of the fine level that belong to no fine grid (the cells across a coarse-fine face) get the value of the
//                       coarse cell that contains them (piecewise-constant, no slope). With ratio >= 2 and two layers, layer 2 is in the same
//                       coarse cell as layer 1, so layer 2 = layer 1, as the FDS scheme of the numeric-scheme note prescribes (Phase 3).
//                       Ghost cells filled by a neighbouring fine grid (same-level exchange) and ghost cells outside a non-periodic domain
//                       edge are left untouched.
#ifndef FDSRT_LEVELOPS_H_
#define FDSRT_LEVELOPS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <functional>
#include <vector>

namespace fdsrt {

// Components scomp .. scomp+ncomp-1 of `fine` are averaged into the same components of `crse` on every coarse cell that lies under a fine grid.
// Other coarse cells are not touched. Valid cells only.
void average_down_cells(const amrex::MultiFab& fine, amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                        const amrex::IntVect& ratio, int scomp, int ncomp);

// Conservative restriction of the species (mass fractions Z with the density rho): rho and rho*Z_n are volume-averaged, and Z_n of the covered coarse cells is rebuilt
// as (rho*Z_n)/rho. Averaging the mass fractions themselves is NOT conservative (the species mass changes by the covariance of rho and Z inside the parent, ~1e-4 relative
// for a flame front) and must not be used for any restriction. zz_* have ns components (rho*Z_n for all n, passive scalars included), rho_* one. `covered_c` marks the
// coarse cells under the fine level (same BoxArray as the coarse data); only those cells of rho_c and zz_c are written. Fine data: valid cells only.
void average_down_species(const amrex::MultiFab& rho_f, const amrex::MultiFab& zz_f, amrex::MultiFab& rho_c, amrex::MultiFab& zz_c, const amrex::iMultiFab& covered_c,
                          const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom, const amrex::IntVect& ratio);

// nlayers (<= fine.nGrow(), <= ratio) ghost layers of `fine`, components scomp .. scomp+ncomp-1, from `crse`. Run the same-level ghost fill first.
// Returns the number of ghost cells written on this rank.
long fill_cf_ghosts_pc(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                       const amrex::IntVect& ratio, int nlayers, int scomp, int ncomp);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// FDS ghost rules at a coarse-fine face (R2b-prep). Source: wall.f90 ASSIGN_GHOST_VALUE (RHO/RHOS, ZZ/ZZS, RSUM, TMP) and velo.f90 VISCOSITY_BC (MU, KRES, D, DS).
// FDS fills only the face-adjacent ghost cell of every wall cell, and it looks at ONE layer of the other mesh: the cells next to the face (probe points
// in init.f90 CHECK_MESHES), not at the whole r^3 children of a coarse cell. Therefore:
//   fine side   (ghost cell = inside the coarse cell, layer 1 and 2):  RHO = coarse value; ZZ = clip((rho*ZZ)/rho) = coarse ZZ; RSUM from ZZ; TMP = PBAR/(RSUM*RHO);
//                layer 2 of RHO, ZZ, TMP = layer 1 (the 4-equal-cells second-order branch does not apply to unequal cells). MU, KRES, D, DS: layer 1 = coarse value.
//   coarse side (ghost cell = covered coarse cell next to the face):   RHO = area-weighted mean of the r_t^2 fine cells of the layer next to the face (ARO = min(1, area ratio));
//                ZZ = clip(sum(ARO rho ZZ) / sum(ARO rho)) (mass weighted); RSUM, TMP rebuilt as above; layer 2 of RHO, ZZ, TMP = layer 1.
//                MU, KRES, D, DS: layer 1 = arithmetic mean of those fine cells.
// The AMReX average-down (volume mean of all r^3 children) is NOT what FDS reads in layers 1 and 2; layers >= 3 of the covered region are never read and keep it.
// Known limit: a covered coarse cell that is the ghost cell of two faces (edge or corner of the fine patch, or a patch thinner than 4 coarse cells) holds ONE value
// here, FDS has one array per mesh. The first candidate in the order x, y, z, layer 1 before layer 2 wins and `conflicts` counts such cells.
// ---------------------------------------------------------------------------------------------------------------------------------------------------------

// The fields of one stage on one level (nullptr = not present). predictor: rho = RHOS, zz = ZZS; corrector: rho = RHO, zz = ZZ. mean = MU, KRES and D (predictor) or DS (corrector).
struct ScalarStage {
    amrex::MultiFab* rho = nullptr;
    amrex::MultiFab* zz = nullptr;
    amrex::MultiFab* tmp = nullptr;
    amrex::MultiFab* rsum = nullptr;
    std::vector<amrex::MultiFab*> mean;
};

// Equation-of-state callbacks for RSUM and TMP of a ghost cell (FDS: GET_SPECIFIC_GAS_CONSTANT and PBAR_P(KK, zone)). Without them TMP and RSUM of a ghost cell
// are the injected (fine side) or averaged (coarse side) values, a documented fallback that differs from FDS at round-off to O(change of state across the face).
struct Thermo {
    int n_tracked = 0;                                           // number of ZZ components passed to rsum
    std::function<double(const double* zz)> rsum;                // specific gas constant of the mixture ZZ(0..n_tracked-1)
    std::function<double(const amrex::IntVect& cell)> pbar;      // background pressure at the (ghost) cell of the level the call is made for
};

struct CfStats {
    long ghost_cells = 0;   // cells written on this rank
    long conflicts = 0;     // coarse side: cells that are the ghost cell of more than one face (see above)
};

// Fine-level ghost cells (layers 1 and 2; layer 1 only for the fields with one ghost layer) from the coarse level by the FDS rule. Run the same-level fill first.
CfStats fill_fine_ghosts_fds(const ScalarStage& fine, const ScalarStage& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                             const amrex::IntVect& ratio, const Thermo* th_fine);

// Covered coarse cells next to a coarse-fine face (layers 1 and 2) from the fine level by the FDS rule. The caller refreshes the same-level ghosts of the coarse
// level afterwards (a coarse box next to the covered region reads these cells as its ghost cells).
CfStats fill_covered_ghosts_fds(const ScalarStage& crse, const ScalarStage& fine, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                                const amrex::IntVect& ratio, const Thermo* th_crse);

}  // namespace fdsrt

#endif
