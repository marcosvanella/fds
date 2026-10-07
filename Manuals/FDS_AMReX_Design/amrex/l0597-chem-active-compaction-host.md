# L0597 chemically active cell compaction: host path

Decision: L0597 stays on the FDS host path. It is not planned for the generator.
Do not translate the chemically active cell mark-and-append nest.

## Why

`work-species/handoff-s4/next/L0597/ANALYSIS.md` confirms the nest at
`fire.f90` **123–137** inside `COMBUSTION_GENERAL_LOAD_BALANCED` (routine 89–338)
on FireX `8ad685f435` and on this repository HEAD. For each gas cell the
body builds `ZZ_GET` / `PRES`, calls `CHECK_CHEMICALLY_ACTIVE_STATE`, may
**RETURN** on `STOP_STATUS/=NO_STOP`, and if active appends `[I,J,K]` (and
optional ignition zone) to `CHEM_ACTIVE_CELLS` in source order (stream
compaction).

A private generator probe of lines 123–137 refused:

`callee CHECK_CHEMICALLY_ACTIVE_STATE uses module variable CELL ... no rule for it`

Blockers recorded in the analysis:

1. Callee reads module `CELL` / other mesh state — not an admissible device
   callee without a large flatten.
2. Stream compaction with a live counter into
   `CHEM_ACTIVE_CELLS(NCHEM_ACTIVE_CELLS,:)` — race under a parallel K,J,I loop;
   needs a serial host pass, a two-pass compact, or a generator append-in-source-
   order feature.
3. `RETURN` on `STOP_STATUS` (fire.f90:129) — early exit of the whole routine,
   not a cell `CYCLE`.
4. Outer `DO NM=...` mesh loop stays on the host.

Claim status in `amrex/loop_claims.csv`: **not planned (host)**.

## Optional later (not this note)

Host builds a compact index list; a separate reaction kernel runs over that list.
That is a host rewrite plus a new kernel contract, not a translation of the
current nest. Not landed; no s5-gen commit.

## Scope of this note

Docs only. No change to `COMBUSTION_GENERAL_LOAD_BALANCED`, markers, or generated
kernels. Line numbers are those cited in the analysis; this note does not invent
new citations. The cut-cell companion nest at fire.f90:153+ is a separate nest
and is not claimed here.

## Related

- `work-species/handoff-s4/next/L0597/ANALYSIS.md` (probe log `s5gen_probe.log`)
- `amrex/loop_claims.csv` row L0597 (status not planned (host))
