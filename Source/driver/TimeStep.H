#ifndef FDS_TIMESTEP_H_
#define FDS_TIMESTEP_H_
// TimeStep.H: replay hooks for T and DT of the FDS main loop (main.f90:763-767, 927-940, 1007), header only.
//
// Kernel-facing rules (M2a): (a) passive scalars are not involved; (b) only uniform Cartesian metrics are used.
//
// The C++ driver owns the time loop (S4). The rules below are exactly those of FDS, expressed on plain doubles so that they can be checked
// bitwise against the STEP and PASS records of the reference dump (tests/README section "Reference dump"):
//   step start        DT = MINVAL(DT_NEW) if every mesh asked for an increase (CHANGE_TIME_STEP_INDEX==1 for all), else DT is kept
//   final time step   DT = MAX(T_END-T+TWENTY_EPSILON_EB, DT_END_MINIMUM) if T+DT+DT_END_FILL > T_END
//   pass of the loop  if any mesh returned index -1 (VELOCITY_PREDICTOR -> CHECK_STABILITY) DT = MINVAL(DT_NEW), FIRST_PASS = false, repeat
//   advance           T = T + DT at the start of the corrector
#include <algorithm>
#include <limits>
#include <vector>

namespace fdsamr {

struct TimeConsts {
    double t_end = 0.0;
    double dt_end_fill = 1.e-6;                                      // DT_END_FILL (cons.f90)
    double dt_end_minimum = 20.0 * std::numeric_limits<double>::epsilon();   // DT_END_MINIMUM = TWENTY_EPSILON_EB
};

// DT for the first pass of a time step, from the DT of the previous step and its final DT_NEW / change index of every mesh.
inline double dt_at_step_start(double t, double dt_prev, const std::vector<double>& dt_new, const std::vector<int>& change_index, const TimeConsts& c)
{
    bool all_increase = !change_index.empty();
    for (int i : change_index) all_increase = all_increase && (i == 1);
    double dt = dt_prev;
    if (all_increase) dt = *std::min_element(dt_new.begin(), dt_new.end());
    if ((t + dt + c.dt_end_fill) > c.t_end) dt = std::max(c.t_end - t + 20.0 * std::numeric_limits<double>::epsilon(), c.dt_end_minimum);
    return dt;
}

// After a pass: true if the pass must be repeated (some mesh asked for a smaller step); dt is then updated to MINVAL(DT_NEW).
inline bool dt_after_pass(double& dt, const std::vector<double>& dt_new, const std::vector<int>& change_index)
{
    bool any_decrease = false;
    for (int i : change_index) any_decrease = any_decrease || (i == -1);
    if (any_decrease) dt = *std::min_element(dt_new.begin(), dt_new.end());
    return any_decrease;
}

inline double time_advance(double t, double dt) { return t + dt; }

}  // namespace fdsamr

#endif
