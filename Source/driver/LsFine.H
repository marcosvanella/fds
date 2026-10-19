// LsFine.H: level set on one global 2-D grid (fds_ls_grid.f90): the C++ side. Every rank holds the whole LS grid (flat terrain, ratio R to the level-0 cells, R = the cumulative refinement
// ratio of the finest level that can exist). The wind at the reference height of every LS cell is gathered from the LEAF cells of all levels (cells not covered by a finer level), one
// MPI_Allreduce(SUM) of the column samples per stage, so every rank advances the same LS grid with the same input.
#ifndef FDS_LSFINE_H_
#define FDS_LSFINE_H_

#include "LevelRegistry.H"

namespace fdsamr {

struct LsFineConfig {
    int mode = 0;    // 0 off, 2 shadow beside the original with the gathered wind (ratio 1, one level-0 mesh), 3 the LS grid drives the fire (--ls-fine)
    int ratio = 1;   // LS cell = level-0 cell / ratio
};
LsFineConfig& ls_fine_config();

// Gathers the wind at the reference height of every spreading LS cell from the leaf cells of `reg` and hands it to the Fortran LS grid (fds_ls_wind_from_slots). Collective.
// Returns the number of spreading columns without a sample (0 expected).
int ls_fine_gather_wind(const LevelRegistry& reg);

}  // namespace fdsamr

#endif
