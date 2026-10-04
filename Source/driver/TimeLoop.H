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
#include "LevelRegistry.H"
#include "SideData.H"
#include "TimeStep.H"

namespace fdsamr {

struct RunOptions {
    std::string chid = "run";      // prefix of the files written
    std::string outdir = ".";
    int max_steps = -1;            // <0: until T_END
    std::string ref_dir;           // baseline directory of the case (optional): MMS and restart comparison data are read by tests/compare_run.py
    bool quiet = false;
    // D-053: the pressure-zone sums USUM, DSUM, PSUM keep FDS's own summation order (per rank in box order, then MPI_ALLREDUCE) by default. true (--exact-zone-sums or
    // env FDSTL_EXACT_ZONES=1) replaces them by the decomposition-independent fixed-point sums of D-028 (ExactSum); a GPU option, off by default, tested separately
    // (tests/run_zone_sum_check.sh). The set-up volumes and areas (zone_setup) are D-028 sums in both modes.
    bool exact_zone_sums = false;
};

struct StepRecord {
    int step = 0;
    double t0 = 0.0, dt = 0.0;     // time at the start of the step and DT of the step (final pass)
    int passes = 0;                // passes of the CHANGE_TIME_STEP loop
    int it_pred = 0, it_corr = 0;  // pressure iterations of the last pass (predictor) and of the corrector
    double perr_pred = 0.0, verr_pred = 0.0, perr_corr = 0.0, verr_corr = 0.0;
    double zone_rel = 0.0;         // max relative difference FDS accumulation vs exact fixed-point zone sums in this step (diagnostic; computed only with exact_zone_sums or FDSTL_ZONE_DIAG, else 0)
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
    Fields& fields() { return *m_F; }      // level 0 (same as fields(0))
    BcStep& bcstep() { return *m_bc; }     // level 0

    // ---- S9: per-level interface (notes/level-interface.md; Role 3 plan section 4 items 1-5). Stable header: these names and signatures are what Role 3 codes against.
    // registry(): the per-level Level/Fields/SideData/covered mask and the fdsrt::LevelListener for Role 3's make/remake/clear_level. Valid after construction of the loop.
    LevelRegistry& registry() { return *m_reg; }
    // Levels that advance() steps: those with an FDS binding. 1 until fine-level FDS MESH_TYPE objects exist; the registry may hold more levels (layout, Fields, SideData).
    int num_levels() const;
    // ---- S12 level binding. bind_level(lev): the registry level `lev` (made by LevelRegistry::make_level, data filled by the transfer) gets FDS fine-box mesh objects (BUILD_FINE_BOX of
    // fds_fine_level.f90, patch 0007 DRAFT), views onto its Fields (FDS_FINE_B_SET_VIEW), a BcStep (ext_ghost on) and the clip scratch; afterwards every per-level stage entry point and the
    // flux hooks accept `lev` (FDS mesh number of box i is fds0 + i + 1 > NMESHES). Levels are bound in order; only the top bound level can be rebound (call it again after a regrid that remade it,
    // or unbind_level first when it was cleared). A tree without patch 0007 aborts with a message. The registry level object must stay alive while bound. advance() with a bound level > 0
    // runs the composite step (S14, below). The single-level path (no bind_level call) is unchanged.
    // The global FDS state the kernels read (PREDICTOR/CORRECTOR, FIRST_PASS, ICYC, density limits) as advance() sets it before a stage, and the reset of the divergence integrals
    // (INIT_DIV) that advance() does before WALL_BC of the predictor. A driver that runs the per-level stage entry points itself calls them at the same places.
    void stage_state(bool predictor, bool first_pass);
    void stage_init_divergence();
    void bind_level(int lev);
    void unbind_level(int lev);
    bool level_bound(int lev) const;
    // Per-level stage entry points. Each is the exact body that advance() runs for that level at that position of MAIN_LOOP (same calls, same order), on the level's boxes
    // of this rank. `predictor` selects the stage (true = predictor). t and dt are the loop's current time and step. A level without an FDS binding aborts with a
    // message naming the missing piece (no fine-level physics is faked). The pressure solve, the zone sums and the outputs are NOT per level: they stay global
    // (the composite pressure solve across levels is Role 2/3 work, notes/level-interface.md section 4).
    void stage_viscosity(int lev, bool predictor);          // COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES
    void stage_density(int lev, bool predictor);            // DENSITY incl. the D-031 clip over this level (the pre/clip/post split form)
    void stage_exchange(int lev, int code, bool predictor); // MESH_EXCHANGE code on the level: same-level fill, then the coarse-fine ghost hook (BcStep::exchange)
    void stage_boundary(int lev, int code);                 // the boundary routines FDS runs after exchange `code` (BcStep::after_exchange)
    void stage_velocity_flux(int lev, bool predictor);      // VELOCITY_FLUX (+ the interface walls as regular faces)
    void stage_wall_bc(int lev, bool predictor);            // WALL_BC
    void stage_divergence1(int lev);                        // DIVERGENCE_PART_1
    void stage_divergence2(int lev);                        // DIVERGENCE_PART_2
    // VELOCITY_PREDICTOR of the level: smallest stable dt of the level's boxes of this rank (a large number if the rank has none) and CHANGE_TIME_STEP flag (OR).
    void stage_velocity_update(int lev, bool predictor, double* local_dt, int* change_time_step);
    void stage_velocity_correct(int lev);                   // VELOCITY_CORRECTOR
    // ---- S10.2 interface flux hooks (notes/flux-hooks-design.md). kind 0 = ADV (stage product face value x face velocity, kg/(m2 s) of rho*Z_n), kind 1 = DIF (RHO_D_DZDX/Y/Z as the divergence reads them).
    // flux_array(lev, kind, dir): nodal MultiFab in `dir`, ncomp N_TOTAL_SCALARS (component n-1 = scalar n), no ghost cells, the level's BoxArray and DistributionMapping; AMReX face a is FDS face I = a-lo.
    // flux_readout_adv(lev): after stage_viscosity; fills the ADV arrays, no cell update. flux_readout_dif(lev): runs DIVERGENCE_PART_1 (= stage_divergence1) and fills the DIF arrays.
    // flux_set_override(lev, kind, lists): replaces the lists (one vector of FluxOverride per local box in MFIter order; an empty vector clears) and switches the kind to override mode: the next
    // stage_density (ADV) and flux_apply_dif (DIF) use them. flux_apply_dif(lev): DIVERGENCE_PART_1 again with the DIF list applied (zone sums restored first), a no-op without a list.
    void flux_readout_adv(int lev);
    void flux_readout_dif(int lev);
    void flux_apply_dif(int lev);
    void flux_set_override(int lev, int kind, const std::vector<std::vector<fdsrt::FluxOverride>>& per_local_box);
    const amrex::MultiFab& flux_array(int lev, int kind, int dir);
    // Global dt: MIN of the per-level local dts over levels and ranks. Same value for both stages of a step (D-050: dt is decided once per pass, both stages use it).
    double global_dt(const std::vector<double>& local_dt_per_level) const;
    // Coarse-fine ghost hook of level `lev` (Role 3 plan item 4): called by stage_exchange after the same-level fill of every exchange code. See CfGhostRequest.
    void set_cf_ghost_hook(int lev, CfGhostHook h);
    // ---- S14 composite step. With a second bound level advance() runs the hierarchy: per-level stages, ONE composite pressure solve (pressure_backend, PressureProblem::levels, one
    // PressureWorkspace kept per layout and rebuilt after bind_level), fine faces overwrite the coarse faces under them (velocity, FVX/FVY/FVZ), fine ghost faces from the coarse level,
    // and the interface flux overwrite of the scalar fluxes (Role 3 FluxStageRunner; FDSTL_OVERWRITE=0 or set_flux_overwrite(false) is the negative control). One dt for all levels.
    // Not built for several levels: the pressure iteration (PRESSURE_ITERATIONS), BAROCLINIC, zone sums over the fine level, outputs of the fine level (the files hold level 0).
    // prime_levels(): one corrector-form pass without the density update on every level (D, DS, KRES, MU of a newly bound level); call it after bind_level and set_state.
    struct CompositeReport {
        long solves = 0, iterations = 0, h_ghost_cells = 0;
        double removed_mean = 0.0, removed_rel = 0.0;   // largest over the solves (backend: mean actually removed from the RHS, relative to its rms)
        std::string backend;
    };
    const CompositeReport& composite_report() const;
    void prime_levels();
    void set_flux_overwrite(bool on);
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
    std::unique_ptr<LevelRegistry> m_reg;
    PressureHook m_hook;
    double m_t = 0.0, m_dt = 0.0, m_rmin = 0.0, m_rmax = 0.0;
    int m_icyc = 0;
    std::vector<double> m_last_local_dt;   //!< per level: result of the last stage_velocity_update (FluxStages::local_dt)
    friend class FluxStages;
    std::vector<StepRecord> m_steps;
};

// Adapter that implements the fdsrt::FluxAccess protocol (RegridInterface.H) on a TimeLoop. FDS orders the work of one stage as viscosity, DENSITY (advective fluxes, cell update),
// exchanges / boundary steps / velocity flux / WALL_BC, DIVERGENCE_PART_1 (diffusive fluxes); the adapter maps the phases as follows (D-061, notes/flux-stage-wiring.md, S12 order):
//   compute_stage_fluxes(l, p)  flux_readout_adv: ADV face fluxes of the stage, no cell update (stage_viscosity must have run); stage_flux(l, Adv, d) valid afterwards;
//   apply_flux_divergence(l, p) stage_density with the ADV override list (cell update only; NO divergence-1 here);
//   readout_dif(l)              after the exchange / boundary / velocity-flux / WALL_BC positions: DIVERGENCE_PART_1 with the DIF read-out; stage_flux(l, Dif, d) valid afterwards;
//   run_divergence_part1(l)     after readout_dif: the DIF override list is applied (DIVERGENCE_PART_1 re-run, zone sums restored; a no-op without a list). Without a preceding
//                               readout_dif on the level it is a plain stage_divergence1 (new levels after a regrid).
// The caller therefore reads ADV and sets the ADV overrides before apply_flux_divergence, and reads DIF and sets the DIF overrides between readout_dif and run_divergence_part1, every
// phase over all levels (finest first for the apply phases) before the next. With no override set the sequence is bitwise the plain stage (tests/run_flux_hook_check.sh).
// max_divergence_error (D-058) is max |div u - D| over the valid cells of the level on this rank (uniform Cartesian differences of the level's U, V, W and D).
class FluxStages : public fdsrt::FluxAccess {
public:
    explicit FluxStages(TimeLoop& tl) : m_tl(tl) {}
    void compute_stage_fluxes(int level, bool predictor) override;
    void apply_flux_divergence(int level, bool predictor) override;
    // D-061 order (notes/flux-stage-wiring.md): the DIF read-out is its own phase, after the exchange / boundary / velocity-flux / WALL_BC positions of the stage. readout_dif(level) runs it
    // (DIVERGENCE_PART_1 with the DIF read-out); stage_flux(level, Dif, d) is valid after it; run_divergence_part1(level) then applies the DIF override (re-run) or, without a read-out, is a plain divergence-1.
    void readout_dif(int level);
    const amrex::MultiFab& stage_flux(int level, fdsrt::FluxKind kind, int dir) const override;
    void set_flux_override(int level, fdsrt::FluxKind kind, const std::vector<std::vector<fdsrt::FluxOverride>>& per_local_box) override;
    double local_dt(int level) const override;
    void run_divergence_part1(int level) override;
    double max_divergence_error(int level) const override;
private:
    TimeLoop& m_tl;
    std::vector<char> m_div1_pending;   //!< per level: DIVERGENCE_PART_1 ran inside apply_flux_divergence and its DIF override is still open
};

}  // namespace fdsamr

#endif
