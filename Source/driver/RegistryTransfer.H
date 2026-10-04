// RegistryTransfer.H: the adapter from the LevelRegistry fields to Role 3's regrid transfer (fdsrt::LevelDataTransfer, RegridInterface.H), S11.2.
//
// It maps the FDS field set of the registry (RHO, RHOS, ZZ, ZZS: mass fractions, TMP, face velocities U V W and US VS WS in the AMReX face index) to the operators of
// CellTransfer.H / FaceTransfer.H / LevelOps.H, which know nothing of FDS:
//   fill_initial_level   t = 0: LevelRegistry::fill_initial_level, the initial conditions evaluated directly on the new level (not interpolated);
//   fill_new_level       conserved densities rho*Z_n (n = 1..N_TOTAL_SCALARS) by prolong_conserved (limited slopes, positivity floor 0, parent-sum clip, clip counter), then
//                        rho = sum of the tracked rho*Z_n, ZZ_n = rho*Z_n / rho (so the mass fractions sum to 1 to round-off), RHOS = RHO, ZZS = ZZ; TMP by injection of the
//                        parent value unless `derive` is given; face velocities by prolong_faces_normal_linear (interface faces = coarse value, every fine cell has the
//                        divergence of its parent), US VS WS = U V W;
//   fill_remade_level    as fill_new_level, then the valid data of the previous level objects (LevelRegistry::retired_fields, kept until end_regrid) are copied over the overlap
//                        BITWISE (all cell fields and faces), so only cells without fine data are interpolated;
//   hierarchy_done       one average-down finest to coarsest (cell scalars by average_down_cells, faces by AMReX average_down_faces: coarse face = mean of the fine faces),
//                        then the same-level and periodic ghost refill of every transferred field. Physical-boundary ghost cells are the boundary-condition step's business.
// Preconditions: the coarse level's valid data are final; for a non-periodic domain edge its ghost layer holds boundary values (the slope limiter reads one layer).
// Not transferred (rebuilt by the kernels in the next stage): the stage arrays H, HS, KRES, D, DS, DDDT, MU, RSUM, FV*, and everything per-box scratch. `derive(level)`, if set, is
// called after the transfer of a level (the place for the equation-of-state recomputation of TMP/RSUM, run_divergence_part1 ...).
#ifndef FDS_REGISTRYTRANSFER_H_
#define FDS_REGISTRYTRANSFER_H_

#include <functional>

#include "LevelRegistry.H"
#include "RegridInterface.H"

namespace fdsamr {

struct RegistryTransferStats {
    long parents = 0, limited = 0, clips = 0, unfixable = 0;   // sums over the prolong_conserved calls
    long copied_old_cells = 0;                                 // valid cells of retired fields copied over the overlap (remade levels)
    int n_initial = 0, n_new = 0, n_remade = 0, n_done = 0;
};

class RegistryTransfer : public fdsrt::LevelDataTransfer {
public:
    // n_tracked: number of ZZ components that sum to 1 (rho is their rho*Z sum); the other components up to N_TOTAL_SCALARS are passive scalars.
    RegistryTransfer(LevelRegistry& reg, int n_tracked) : m_reg(reg), m_ntracked(n_tracked) {}

    LevelRegistry::InitialFill initial;                    // t = 0 evaluation (required for fill_initial_level)
    std::function<void(int level)> derive;                 // optional, called after a level was filled by fill_new_level / fill_remade_level / fill_initial_level
    std::function<void(bool initial)> after_hierarchy;     // optional, called at the end of hierarchy_done

    void fill_initial_level(const fdsrt::LevelLayout& l) override;
    void fill_new_level(const fdsrt::LevelLayout& l) override;
    void fill_remade_level(const fdsrt::LevelLayout& l) override;
    void level_cleared(int) override {}
    void hierarchy_done(bool initial) override;

    const RegistryTransferStats& stats() const { return m_st; }

private:
    void prolong(int level, bool remade);
    LevelRegistry& m_reg;
    int m_ntracked;
    RegistryTransferStats m_st;
};

}  // namespace fdsamr

#endif
