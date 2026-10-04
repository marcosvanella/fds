// Composite.H: internal entry points of the composite (multi-level) pressure path. Not part of the public
// interface (PressureIface.H); used by PressureIface.cpp.
#ifndef PB_COMPOSITE_H_
#define PB_COMPOSITE_H_

#include "PressureIface.H"

namespace pb {

// Selector for a problem with a non-empty `levels`. NotBuilt reasons only (structure and flags); layout errors are
// InvalidInput from validate_composite().
Selection select_composite (PressureProblem const& p, BackendKind requested);
// Layout checks (BoxArrays, domains, nesting, ratios, matching rhs/phi). Empty string = valid. structural = false skips
// the BoxArray/domain/nesting checks (cheap path with a workspace that was built, hence validated, for this layout).
std::string validate_composite (PressureProblem const& p, bool structural = true);
// Full composite solve. The selector and the validation have already passed.
PressureResult solve_composite (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws, bool full = true);
// Single-level rebuild of a PressureWorkspace (defined in PressureIface.cpp).
Status rebuild_single_level (PressureWorkspace& ws, PressureProblem const& p, std::string* message);

} // namespace pb
#endif
