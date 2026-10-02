# Upstream patches (decision D-051)

Proposed changes to upstream FDS source (FireX or master) that come out of the AMR/GPU work, for the owner to review and commit upstream himself. Nobody else commits or pushes them upstream.

## Rules for a patch
1. One short-scope change per patch (one marker set, one loop cleanup, one routine). No mixed patches.
2. Delivered as a file in this folder: `git format-patch` output or a unified diff against FireX, applicable to the upstream source file as it is, without our `WITH_AMREX` code.
3. A short plain-language rationale in the patch description: what changes, why the GPU or AMR work needs it, what it does not change.
4. A behavior-unchanged check: the exact command or test and its result (for example the V&V baseline comparison against the reference binary, bitwise or within the stated tolerance), run on the patched upstream file.
5. A target line: **FireX**, **master**, or **both**. Master is merged into FireX periodically, so write the patch to apply cleanly to both and say which one it was tested against.
6. Name: `NNNN-short-name.patch` (four-digit sequence), and add a row to the index below in the same change.

## Index
| No. | File | Target | Source file(s) | Summary | Behavior-unchanged check | Author | State (proposed / reviewed / committed upstream by owner) |
|---|---|---|---|---|---|---|---|
| (none yet) | | | | | | | |

## Blocked-loop family reviews (decision D-051, item 1)
Each blocked-loop family is signed off by a domain lead before its rewrite becomes a patch.

| Family | Reviewer | Status |
|---|---|---|
| Zone sums with a pressure-zone index (`USUM(IPZ)`, `DSUM`, `PSUM`) | Pressure Solver Lead | requested 2026-10-02 |
| `DELTA_RHO_ZZ` scatter, species and combustion loops | Species & Combustion Lead | requested 2026-10-02 |
| Solid-phase counters and wall loops (`CELL_COUNTER` and similar) | Solid Phase Lead | requested 2026-10-02 |
| Radiation loops | Radiation Lead | requested 2026-10-02 |
| `CHECK_STABILITY` reductions and other cross-cutting reductions | V&V Lead | requested 2026-10-02 |
