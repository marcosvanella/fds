# M2 notes: solve triggers (FR-039), FFT plan cache, mixed faces, boundary data, D-057 mapping

This file is appended to as the M2 items land. Sections are independent.

## 1. Solve-trigger points (FR-039)

`PressureOptions::trigger` is a bit set of `SolveTrigger` values describing why this solve is being made:

| flag | value | meaning |
|---|---|---|
| `SolveRoutine` | 0 | ordinary time-step solve |
| `SolveFirst` | 1 | first solve with a workspace (detected automatically from `PressureWorkspace::num_solves() == 0`) |
| `SolveFirstAfterRegrid` | 2 | first solve after a regrid (set automatically by `PressureWorkspace::rebuild()` and when a solve sees a layout that differs from the workspace; the caller may also set it) |
| `SolveDebug` | 4 | debug run, check everything |

`PressureOptions::full_checks_on` (default `SolveFirst | SolveFirstAfterRegrid | SolveDebug`) selects which flags switch the full checks on.
`trigger` defaults to `SolveDebug`, so a caller that never touches the option gets the full checks on every solve, as before.
The effective flags (caller flags plus the ones the workspace detected) are returned in `PressureResult::triggers`;
`PressureResult::full_checks` says which path ran.

Full checks: structural layout validation, pin scan, rms / removed-mean compatibility diagnostic and its warning,
true residual (`residual_checked`).
Cheap path: selector, rhs/phi layout checks, mean removal, backend solve, gauge, convergence status.
`NotBuilt`, `InvalidInput` and `NotConverged` are reported on both paths.
The solution, the removed mean and the gauge constants are bitwise identical on both paths (tests below).

Tests: `pb_trigger_single_NNNNNN`, `pb_trigger_single_PPPPPP` (FFT and MLMG, np 1 and 2), `pb_comp_trigger` (hierarchy, np 1 and 2).

## 2. FFT plan caching

`PressureWorkspace` owns an `FFTBackend` that caches the `amrex::FFT::Poisson` object, keyed by BoxArray, DistributionMapping,
Geometry and the six boundary types. A repeated solve with the same key reuses the plan (`BackendStatus::plan_reused`,
`PressureWorkspace::fft_plan_reuses()`); any key difference, `rebuild()`, `discard()` or a composite rebuild builds a new plan
(`fft_plan_builds()`). A workspace result is bitwise identical to a no-workspace result (test `pb_fftcache_*`).
Single-level MLMG has no cache (the MLMG object is cheap to build; not requested).

Timing, 20 solves each, 1 core machine with oversubscribed ranks (so noisy; ratios are what matter), Neumann FFT, Routine trigger:

| grid | ranks | fresh plan per solve | cached plan per solve | ratio |
|---|---|---|---|---|
| 32^3 | 1 | 5.5 ms | 4.4 ms | 1.26 |
| 32^3 | 4 | 5.2 ms | 4.9 ms | 1.07 |
| 64^3 | 1 | 45.5 ms | 36.0 ms | 1.26 |
| 64^3 | 4 | 14.9 ms | 12.2 ms | 1.22 |
| 128^3 | 1 | 340 ms | 324 ms | 1.05 |
| 128^3 | 4 | 104 ms | 97 ms | 1.07 |

The saving is the plan construction, a few percent to about a quarter of a solve; it is not large because the transforms
dominate. `plan_build_ms` in the harness line is the time of an explicit `rebuild()`.
