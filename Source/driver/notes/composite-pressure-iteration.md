# Composite pressure iteration (not built)

Pressure Lead finding: with refinement the multi-level branch (`TimeLoop::pressure_scheme_composite`, around the composite solve) never computes the residual diagnostic of loop L1208 (`pres.f90`), so `VELOCITY_TOLERANCE` and `PRESSURE_TOLERANCE` have no effect on refined runs. The host functions for that check are ready when the iteration is built. Do not implement the iteration here; record only.
