# Gap: FDS H dump with nonzero wall data for fold sign-off (pressure/07 §11.4)

**Status: open / blocked on FDS run capacity — Role 2 did not invent a dump.**

## What §11.4 asks for

Repeat the A-56 fold sign check (`docs/pressure/07` §11.2) against **true FDS H**, not only the dense FDS-style matrix. That needs an FDS run with:

- Nonzero wall data (`BXS..BZF` from nonzero normal velocity or `DUNDT` on a wall), and
- An H (or HS) dump of that solve.

A-09b H-dump request was recorded as not yet requested in the V&V area. A-38 runs are periodic (no wall data).

## What already exists (not sufficient for §11.4)

| Asset | Why it does not close §11.4 |
|---|---|
| `frozen/fds_csmag32_periodic/` + UP-0011 `FDS_PDUMP_STEPS` | Real FDS `rhs.bin` / `phi.bin`, but **all faces periodic**, `bc_data_maxabs_* = 0` — no inhomogeneous fold to check |
| Dense-matrix + MMS fold tests (`pb_bcdata_*`, A-56 §11.2–11.3) | Confirm FISHPAK signs without FDS H |
| UP-0011 patch (`upstream-patches/UP-0011-pressure-rhs-dump.patch`) | **Mechanism** to dump `phi.bin` (= H/HS) + `bc.bin` (= BXS..BZF) once a suitable INPUT is run |

## How to close (for Mapper / V&V / Role 1 with run capacity)

1. Scratch FireX (or tip) + apply UP-0011; build.
2. Choose a **single-mesh** input with at least one OPEN or moving-wall face so `bc_data_maxabs_* > 0` in `meta.txt` (avoid multi-mesh until `pb_fds_frozen` multi-mesh is agreed).
3. `FDS_PDUMP_STEPS=<icy>` run; keep `*_meta.txt`, `*_rhs.bin`, `*_phi.bin`, `*_bc.bin`.
4. Land under `frozen/` (or vv-runs) with a README modelled on `fds_csmag32_periodic/README.md`.
5. Optional: extend `pb_fds_frozen` or a small fold harness to: fold `bc.bin` into stored rhs (or rebuild rhs), solve with `solve_pressure`, compare to `phi.bin` within eps_H.

Until that dump exists, fold signs remain accepted on dense/MMS evidence only (§11.5). Do not treat periodic pdumps as fold validation.
