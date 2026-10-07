# L0630 BLOCK_CELL: host path (D-083 / next-work r6)

Decision: L0630 stays on the FDS host path. Do not translate the new single-cell
consumable-obstruction branch. MESH_EXCHANGE codes 66 and 67 stay on the host
path too.

## Why

next-work r6 (Chief Architect refresh through D-083) lists L0630 as a Mesh Data
item: decide translation or host path for the new single-cell branch after the
FireX merge, with the recommendation host path because the routine is
event-driven (it runs only on create or remove of an obstruction). The Architect
and CEO recommend the same host path.

The merge impact note (`upstream-1fb5c300-loop-impact.md` section 2) holds L0630
because the body grew a single-cell consumable-obstruction branch that reads
`WALL_INDEX`, neighbour `TMP`, `PBAR_S`, `PRESSURE_ZONE`, `RSUM0` and
`OBSTRUCTION%CONSUMABLE`, and writes `TMP`, `RHO` and `RHOS`. A translation
would need those as kernel and shim arguments and should be paired with L0817
(`CREATE_OR_REMOVE_OBST`). Keeping the loop on the host avoids that surface
change while the event-driven call rate stays low.

MESH_EXCHANGE codes 66 (full `PACK_WALL` at set-up) and 67 (`PACK_WALL_LITE`
every step) carry the back-wall exchange that left code 6 in the same merge.
next-work r6 and the merge notes say the driver's code 6 fills only `U`, `V`,
`W`, `H` ghosts; codes 66 and 67 must remain on the FDS host path and must not
be routed through the driver ghost fill.

## Scope of this note

Docs only. No change to `BLOCK_CELL`, `MESH_EXCHANGE`, markers, or generated
kernels. Line numbers for the nest are those in next-work r6 and the impact
report (after the merge); this note does not restate them as new citations.

## Related

- D-083 (GPU-only two-pass `CHECK_MASS_DENSITY`; separate from this hold)
- next-work r6, Mesh Data Loops Engineer item 2 and Role 3 item on exchange codes
- `amrex/upstream-1fb5c300-loop-impact.md` (L0630 hold; codes 66/67 outside the work list)
- `amrex/o2-edge-exchange-scope.md` (O2 neighbour-mesh host rulings under D-055; different family)

## Mapper notes (further reasons for the host path)

If L0630 were translated later it would need L0817 (`CREATE_OR_REMOVE_OBST`,
`init.f90` after the merge) in the same order: that routine resets `RHO`/`RHOS`
of the same cells after `BLOCK_CELL`, so a kernel alone would not match the
create/remove path. Also, in the new single-cell branch the average
`TMP_AVG/REAL(COUNT,EB)` (cited in `upstream-1fb5c300-loop-impact.md` section 2
as the statement at `func.f90:5616` on the merged tip) divides by zero when every
`WALL_INDEX(±1..±3)` entry is non-zero. Keeping the nest on the host avoids both
the paired L0817 surface and that divide-by-zero edge until a joint design exists.
