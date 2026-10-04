// PbWorkspaceImpl.H: private definition of PressureWorkspace::Impl (composite hierarchy set-up and the single-level FFT
// backend with its cached plan). Included by CompositeSolve.cpp and PressureIface.cpp only.
#ifndef PB_WORKSPACE_IMPL_H_
#define PB_WORKSPACE_IMPL_H_

#include "PressureIface.H"
#include "PressureBackend.H"
#include "HypreBackend.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#include <cmath>

namespace pb {

using namespace amrex;

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
    std::unique_ptr<MLPoisson> mlp;
    std::unique_ptr<MLMG> mlmg;                       // refers to *mlp: declared after it, destroyed first
    // Single level: the FFT backend object, which owns the cached FFT::Poisson plan (null until needed).
    std::unique_ptr<PressureBackend> fft;
    // HYPRE backend set-up (matrix + BoomerAMG hierarchy), built on the first HYPRE solve of this layout and dropped with the Impl
    // (rebuild() / discard()). Single level: the backend object (it owns the HypreSystem); composite: the system itself.
    std::unique_ptr<PressureBackend> hypre;
    std::unique_ptr<HypreSystem> hsys;
};

} // namespace pb
#endif
