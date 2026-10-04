// Composite.H: internal entry points of the composite (multi-level) pressure path. Not part of the public
// interface (PressureIface.H); used by PressureIface.cpp.
#ifndef PB_COMPOSITE_H_
#define PB_COMPOSITE_H_

#include "PressureIface.H"

namespace pb {

// Selector for a problem with a non-empty `levels`. NotBuilt reasons only (structure and flags); layout errors are
// InvalidInput from validate_composite().
Selection select_composite (PressureProblem const& p, BackendKind requested);
// Layout checks (BoxArrays, domains, nesting, ratios, matching rhs/phi). Empty string = valid.
std::string validate_composite (PressureProblem const& p);
// Full composite solve. The selector and the validation have already passed.
PressureResult solve_composite (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws);

} // namespace pb
#endif
