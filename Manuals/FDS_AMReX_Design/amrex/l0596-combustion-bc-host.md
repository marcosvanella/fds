# L0596 COMBUSTION_BC: host path (same-mesh empty under D-055)

Decision: L0596 stays on the FDS host path. Not planned for the generator. Under
D-055 the same-mesh body is empty (`NOM==0` does `CYCLE`); neighbour-mesh
(`OMESH`) work is excluded from the device route.

## Why

`work-species/handoff-s4/next/L0596/ANALYSIS.md` confirms
`COMBUSTION_BC` at `fire.f90:1902-1918` on FireX `8ad685f435` and on
this repository HEAD. Re-read on src tip (WALL_LOOP 1902–1918): the
wall loop body runs only when `EWC%NOM /= 0`; `NOM == 0` does `CYCLE WALL_LOOP`.
The body sums `OMESH(NOM)%Q` over the neighbour box and stores the average into
`Q(BC%II,BC%JJ,BC%KK)`. No same-mesh body remains when `NOM==0`.

Ruling D-055: no GPU kernels for neighbour-mesh (`NOM>0`) branches. The
same-mesh part is exactly the `NOM == 0` path, which is empty — nothing to
translate for same-mesh. A private generator probe of the OMESH designator
refused (`derived-type designator OMESH(NOM) % Q(...)`), matching
`handoff-s4/next/blockers_L0596_L1273.md`.

Claim status in `amrex/loop_claims.csv`: **not planned (host)**.

## Structural / bitwise host fixture (private)

`work-species/handoff-s4/next/L0596/host_test/` encodes the WALL_LOOP filter
shape and proves:

1. Fixture with only `NOM==0` walls: claimed same-mesh path executes **zero
   body iterations**; sentinel `Q` is unchanged bit for bit (`-O0`, 1 case).
2. Mixed fixture (`NOM` in `{0,2,0,3}`): claimed body count stays 0; the two
   `NOM>0` walls are counted as **excluded** from the D-055 claim (OMESH body
   not run). Sentinel `Q` unchanged.

Build: `gfortran -O0 -g -o comb_bc_same_mesh_empty comb_bc_same_mesh_empty.F90`
→ exit 0, `PASS: L0596 same-mesh filter empty (2 cases)`. Not landed on s5-gen
(nothing to generate). Upstream snippet:
`combustion_bc_upstream_snippet.f90.txt`.

## Scope of this note

Docs + private structural fixture only. No flat per-wall Q table, no s5-gen
commit. Any later host-fill + wall-store rewrite for `NOM>0` after exchange is
outside D-055 device scope and is not planned here.

## Related

- `work-species/handoff-s4/next/L0596/ANALYSIS.md`
- `work-species/handoff-s4/next/L0596/host_test/` (structural fixture)
- `work-species/handoff-s4/next/blockers_L0596_L1273.md`
- `amrex/o2-edge-exchange-scope.md` (D-055 neighbour-mesh host rulings; different family)
- `amrex/loop_claims.csv` row L0596 (status not planned (host))
