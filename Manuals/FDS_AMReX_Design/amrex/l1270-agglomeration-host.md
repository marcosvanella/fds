# L1270 CALC_AGGLOMERATION: host path (short term)

Decision: L1270 stays on the FDS host path for now. Do not translate
`CALC_AGGLOMERATION` through the generator. A later device path, if any, is a
hand-written K2-style kernel with a fixed maximum number of sub-steps (or a
host-chosen counted `N_SUB`), not an s5-gen translation of the adaptive
`DO WHILE`.

## Why

`work-species/handoff-s4/next/L1270/ANALYSIS.md` confirms the nest at
`soot.f90` SUBROUTINE `CALC_AGGLOMERATION` **265–391** (geometry nest **282–388**,
`STEPLOOP: DO WHILE` **356–383**) on FireX `8ad685f435` and on this repository
HEAD. Outer `SPEC_LOOP` over agglomeration species; per gas cell Brownian/grav/shear
`PHI`, then per-cell adaptive agglomeration; write-back to `ZZ`.

A private generator probe of lines 282–388 hit, in order:

1. Intrinsic `ALL` in the loop body (`IF (ALL(N0 < MIN_AGGLOMERATION)) CYCLE` at
   line 289) — refused as a reduction/inquiry/array intrinsic.
2. After replacing `ALL` with a scalar compare: `DO WHILE` is not a supported loop
   form (`TypeError` in `leaf_events` unpacking a `DO WHILE` bound).
3. After replacing `STEPLOOP: DO WHILE` with a fixed counted DO: intrinsic
   `MAXVAL` in the loop body; further strips hit `SUM` and other unresolved names.

Exact blockers recorded in the analysis: array intrinsics `ALL` / `MAXVAL` /
`SUM` (and likely others on bin vectors); adaptive `DO WHILE (DT_SUM < DT)`;
private arrays sized by `MAXVAL(N_PARTICLE_BINS)`; rank-2 `PHI` and bin tables;
callees `GET_VISCOSITY` / `CUNNINGHAM`. Outer `SPEC_LOOP` can remain on the host.

Claim status in `amrex/loop_claims.csv`: **blocked** (host or hand-written; no
s5-gen land from this analysis).

## Preferred short term / device later

- **Short term:** host path, with a private bitwise unit test of the agglomeration
  step against upstream on fixed `PHI` / `N0` fixtures (see
  `work-species/handoff-s4/next/L1270/` when present).
- **Device later:** hand-written K2-style kernel for the per-cell sub-step with a
  fixed max sub-steps (or host-computed `N_SUB` counted DO), plus either generator
  admission of `ALL`/`MAXVAL`/`SUM` on private vectors of known extent or
  precomputed scalar predicates. Do not commit a generator change for this nest
  until that design is settled.

## Scope of this note

Docs only. No change to `CALC_AGGLOMERATION`, markers, or generated kernels. Line
numbers are those cited in the analysis; this note does not invent new citations.

## Related

- `work-species/handoff-s4/next/L1270/ANALYSIS.md` (probe log `s5gen_probe.log`)
- `amrex/loop_claims.csv` row L1270 (status blocked)
- L1273 (`SOOT_SURFACE_OXIDATION`) is a separate nest; still blocked on its own
  grounds (`next/L1273/ANALYSIS.md`)
