# Standing coding rules

Owner: project owner (recorded as D-084). Applies to kernels, host paths, the driver, and generator output.

## Simplicity (D-084)

Write code as simple as possible. Prefer the shorter clear form over clever or copy-heavy alternatives, and keep it readable and maintainable.

Minimise:

- lines
- FLOPs
- conditionals
- copy operations

The best FLOP is the one we don't compute, the best conditional is the one we don't test, and the best byte is the one we don't copy.

### Data layout

Always look for the data layout that maximizes usage and memory-management optimizations — prefetching, single-direction messaging, and similar. Prefer layouts that enable those over copy-heavy or round-trip patterns when correctness allows.

On conflict, correctness and existing gates win: bitwise and class tolerances, and ADR-001 K2 rules, are not overridden by this page. D-084 only guides how those requirements are implemented.

See also: decision log D-084 in `README.md`; implementer pointer in `next-work.md`; ADR-001 for K2.
