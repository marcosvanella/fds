# Generator-side decisions awaiting entry in the decision log

The decision log (`docs/README.md`, "Decision log", D-xxx rows) is owned by the Spec & Program Lead and the numbering by the
Architect's log; no D-number is assigned here. Each entry below is a plain-language text ready to be added as a row, with the
documents it already changed.

## Decision A: the lower-bound argument change is deferred

**Text.** The change that would add box lower/upper bound arguments (`ILO, IHI, JLO, JHI, KLO, KHI`) to every generated kernel is
deferred. Nobody starts it. The question goes back to the CEO only when a specific loop truly needs it.

**Consequences.**
- The kernels keep the FDS index space of one mesh: arrays are explicit-shape dummies whose bounds are the FDS bounds, and the
  first interior cell has index 1. A kernel cannot be told where its box sits in the mesh.
- Running a mesh as several boxes therefore requires the driver to hand each box to a kernel as its own array with local origin 1
  (box-local arrays with the same ghost widths as the FDS allocation, box-local wall tables, and every 1-D metric array such as
  `R`, `RDX..RDZN`, `X`, `XC`, `Z`, `ZC` sliced to the box). A box whose lower corner is not 1 is not handled by the kernels.
- Functions of position called inside a kernel see box-local coordinates.
- Single-box (whole mesh) runs are unaffected. The existing bitwise tests are unchanged.
- The cost that disappears from the plan: 4 to 5 working days for the shared change and 0.5 day of rebasing per other engineer.

**Changed documents.** `generator-projection.md` (assumption becomes "bounds change not done"), `stage1-generator-answers.md` Q2.

## Decision B: summation order for zone sums and other sums

**Text.**
- Default: keep the FDS order. For `USUM`, `DSUM`, `PSUM` (the zone sums) the GPU computes the individual terms (they are
  elementwise and bitwise equal to FDS) and the additions are done serially in the FDS order, either on the host or in a
  single-thread device pass. The result is bitwise equal to FDS for one process.
- Option: a compile-time or macro switch selects exact fixed-point sums on the GPU instead. It is off by default and tested
  separately. Its result does not depend on the box layout or thread count and is the correctly rounded exact sum; it is
  documented as not bitwise equal to FDS.
- The GPU Generator Engineer implements the switch and its separate tests.

**Changed documents.** `blocked-loop-families.md` P1 and `stage1-generator-answers.md` Q3 (and the matching exact-sum rows of
`generator-projection.md`).
