// TimeLoop.H: the C++ time loop of the AMR driver (M2a, S5): MAIN_LOOP of main.f90 (predictor / corrector, CHANGE_TIME_STEP loop, PRESSURE_ITERATION_SCHEME)
// on the level-0 MultiFabs, with the unmodified FDS kernels called box by box (S3), the ghost fill and boundary-condition step of S4, global DT control by MIN/MAX
// reductions, exact fixed-point divergence integrals (D-028) and the pressure solve through pressure_backend/PressureIface.H.
//
// Kernel-facing rules (M2a):
//  (a) Passive scalars: N_TOTAL_SCALARS beyond the tracked species are handled by Fields.cpp: ZZ/ZZS carry ncomp = N_TOTAL_SCALARS and every exchange, copy and clip
//      here loops over all components. The time loop itself only special-cases the tracked species where the D-031 clip does.
//  (b) Only uniform Cartesian metrics are used; R(I)/RRN(I) are dropped, CYLINDRICAL and TRN* meshes are rejected at level-0 assembly (IR-002); the pressure
//      call-through assumes IPS=0 (no stretched direction) and unit-coefficient Cartesian Laplacian.
//
// Scope of this loop (refused with a message otherwise): gas phase, FFT pressure flag, one mesh level, no complex geometry, HVAC, radiation, reactions, particles,
// level set, tunnel preconditioner, obstruction creation, no inhomogeneous pressure boundary data (periodic or closed faces; closed faces only where the boundary
// data is zero or the direction is not solved: TWO_D). S6 replaces these limits.
#ifndef FDS_TIMELOOP_H_
#define FDS_TIMELOOP_H_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Fields.H"
#include "GhostExchange.H"
#include "SideData.H"
#include "TimeStep.H"

namespace fdsamr {

struct RunOptions {
    std::string chid = "run";      // prefix of the files written
    std::string outdir = ".";
    int max_steps = -1;            // <0: until T_END
    std::string ref_dir;           // baseline directory of the case (optional): MMS and restart comparison data are read by tests/compare_run.py
    bool quiet = false;
};

struct StepRecord {
    int step = 0;
    double t0 = 0.0, dt = 0.0;     // time at the start of the step and DT of the step (final pass)
    int passes = 0;                // passes of the CHANGE_TIME_STEP loop
    int it_pred = 0, it_corr = 0;  // pressure iterations of the last pass (predictor) and of the corrector
    double perr_pred = 0.0, verr_pred = 0.0, perr_corr = 0.0, verr_corr = 0.0;
    double zone_rel = 0.0;         // max relative difference FDS accumulation vs exact fixed-point zone sums in this step (diagnostic)
};

class TimeLoop {
public:
    TimeLoop(const Level0& l0, double dt_setup, const RunOptions& o);
    ~TimeLoop();
    // Hook used by the kernel check (tests/kernelcheck.cpp, --chain): called instead of the pressure iteration scheme of a stage; returns true if it took over
    // (it must leave H/HS, FVX.. as FDS has them after PRESSURE_ITERATION_SCHEME).
    using PressureHook = std::function<bool(bool predictor, int pass)>;
    void set_pressure_hook(PressureHook h) { m_hook = std::move(h); }
    // State set by the kernel check before a chained step: time, DT of the step, cycle counter; dt_fixed: the step start uses this DT unchanged.
    void set_state(double t, double dt, int icyc);
    int run();                 // the whole run (writes the outputs); returns the number of failures (0 = completed)
    bool advance();            // one MAIN_LOOP iteration; false if the run must stop (instability)
    // Bind the level-0 fields to the FDS arrays of every local box. Copies the FDS set-up state into the FABs first (the binding replaces the FDS arrays).
    void bind_fields(bool copy_setup_state);
    Fields& fields() { return *m_F; }
    BcStep& bcstep() { return *m_bc; }
    double time() const { return m_t; }
    double dt() const { return m_dt; }
    int icyc() const { return m_icyc; }
    const std::vector<StepRecord>& steps() const { return m_steps; }
    void set_rho_limits(double rmin, double rmax) { m_rmin = rmin; m_rmax = rmax; }

private:
    struct Impl;
    std::unique_ptr<Impl> m;
    const Level0& m_l0;
    RunOptions m_o;
    std::unique_ptr<Fields> m_F;
    std::unique_ptr<SideData> m_sd;
    std::unique_ptr<BcStep> m_bc;
    PressureHook m_hook;
    double m_t = 0.0, m_dt = 0.0, m_rmin = 0.0, m_rmax = 0.0;
    int m_icyc = 0;
    std::vector<StepRecord> m_steps;
};

}  // namespace fdsamr

#endif
