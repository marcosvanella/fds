// PbWorkspaceImpl.H: private definition of PressureWorkspace::Impl (composite hierarchy set-up and the single-level FFT
// backend with its cached plan). Included by CompositeSolve.cpp and PressureIface.cpp only.
#ifndef PB_WORKSPACE_IMPL_H_
#define PB_WORKSPACE_IMPL_H_

#include "PressureIface.H"
#include "PressureBackend.H"
#include "CommonLayer.H"
#ifdef PB_WITH_HYPRE
#include "HypreBackend.H"
#endif

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#include <cmath>
#include <memory>

namespace pb {

using namespace amrex;

struct MaskedMlmg;                                    // CompositeSolve.cpp

struct PressureWorkspace::Impl {
    int nlev = 0;
    std::vector<BoxArray> ba;
    std::vector<DistributionMapping> dm;
    std::vector<Geometry> geom;
    std::vector<IntVect> ratio;                       // ratio[l]: level l to level l-1 (ratio[0] = 1)
    std::array<BC,6> bc{};
    int hidden = -1;                                  // one-cell direction (extruded for the MLMG solve), or -1
    // One-cell direction: the solve runs on an extruded copy of the hierarchy (that direction made periodic with
    // ext_n * product(ratios) isotropic cells, same boxes otherwise and same DistributionMapping). AMReX's
    // "hidden dimension" mode was tried first and is not used: with more than one AMR level it converged slowly
    // or diverged (see frozen/composite-notes.md).
    int ext_n = 4;
    std::vector<BoxArray> eba;
    std::vector<Geometry> egeom;
    std::vector<std::unique_ptr<iMultiFab>> unc;      // 1 = uncovered
    std::vector<Real> vol;                            // cell volume per level
    std::vector<Long> nunc;                           // number of uncovered cells per level
    std::unique_ptr<MLCellLinOp> mlp;                 // MLPoisson, or MLABecLaplacian (A = 0, B = -1) with Robin faces for face_neumann
    std::array<std::vector<int>,6> face_neumann;      // p.face_neumann this hierarchy was built for
    std::unique_ptr<MLMG> mlmg;                       // refers to *mlp: declared after it, destroyed first
    // Single level: the FFT backend object, which owns the cached FFT::Poisson plan (null until needed).
    std::unique_ptr<PressureBackend> fft;
    std::unique_ptr<PressureBackend> mlmg1;           // single-level MLMG backend: owns the cached operator (unmasked; the solver study, P2)
    // Masked single level (D1): component labels, singular flags and pins, keyed on the mask content and the layout. The mask is static
    // for a run and changes at a regrid; the flood fill is the expensive part of a masked solve, so it is done once per mask.
    struct MaskCache {
        bool valid = false;
        std::uint64_t hash = 0;
        amrex::Box dom;
        std::array<BC,6> bc{};
        amrex::BoxArray ba;
        amrex::DistributionMapping dm;
        ComponentMap cm;                              // pristine (no removed_mean / gauge_shift)
    } mask;
    // HYPRE backend set-up (matrix + BoomerAMG hierarchy), built on the first HYPRE solve of this layout and dropped with the Impl
    // (rebuild() / discard()). Single level: the backend object (it owns the HypreSystem); composite: the system itself.
    // Masked composite MLMG operator and solver (P3; CompositeSolve.cpp), keyed on the mask content; dropped with the Impl at a regrid.
    std::shared_ptr<MaskedMlmg> mmg;
    // Masked composite component labels (CompositeMask.cpp), valid for the layout of this Impl and the mask stamp lab_stamp (composite_class_hash).
    CompositeComponentMap lab;
    std::uint64_t lab_stamp = 0;
    bool lab_valid = false;
    // Unmasked operator of face_gradient_composite (built once per layout; the gradient is of the unmasked operator, then masked afterwards).
    std::shared_ptr<Impl> grad;
#ifdef PB_WITH_HYPRE
    std::unique_ptr<PressureBackend> hypre;
    std::unique_ptr<HypreSystem> hsys;
#endif
};

} // namespace pb
#endif
