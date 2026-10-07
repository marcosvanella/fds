# L0596 COMBUSTION_BC: host path (same-mesh empty under D-055)

Decision: L0596 stays on the FDS host path. Not planned for the generator. Under
D-055 the same-mesh body is empty (`NOM==0` does `CYCLE`); neighbour-mesh
(`OMESH`) work is excluded from the device route.

## Why

`work-species/handoff-s4/next/L0596/ANALYSIS.md` confirms
`COMBUSTION_BC` at `fire.f90:1902-1918` on FireX `8ad685f435` and on
this repository HEAD. The wall loop body runs only when
`EWC%NOM /= 0`; `NOM == 0` does `CYCLE WALL_LOOP`. The body sums
`OMESH(NOM)%Q` over the neighbour box and stores the average into
`Q(BC%II,BC%JJ,BC%KK)`.

Ruling D-055: no GPU kernels for neighbour-mesh (`NOM>0`) branches. The
same-mesh part is exactly the `NOM == 0` path, which is empty — nothing to
translate for same-mesh. A private generator probe of the OMESH designator
refused (`derived-type designator OMESH(NOM) % Q(...)`), matching
`handoff-s4/next/blockers_L0596_L1273.md`.

Claim status in `amrex/loop_claims.csv`: **not planned (host)**.

## Scope of this note

Docs only. No flat per-wall Q table, no s5-gen commit, no bitwise test (same-mesh
has nothing to compare). Any later host-fill + wall-store rewrite for `NOM>0`
after exchange is outside D-055 device scope and is not planned here.

## Related

- `work-species/handoff-s4/next/L0596/ANALYSIS.md`
- `work-species/handoff-s4/next/blockers_L0596_L1273.md`
- `amrex/o2-edge-exchange-scope.md` (D-055 neighbour-mesh host rulings; different family)
- `amrex/loop_claims.csv` row L0596 (status not planned (host))
