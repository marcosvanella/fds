# A-14 / D-019 — GNU HDF5 VTKHDF smoke (V&V)

Owner: AMR V&V Lead. Joint with GNU Build Chief. **Not** the owner show-and-tell viz gate
(that remains native AMReX plotfile; see `show_and_tell_notes.md`).

## Status

**PASS** on tip `047adca89bc614fec143c702aed73139ed2cbc5b` (short `047adca89bc6`).

Hub (GNU): `vv-runs/a14-hdf5-gnu/RESULT.md` (+ `PROVENANCE.md`). Owner-facing environment
pointer already landed in docs commit `02f7197` (`vv/environment.md`); this note is the
V&V verify record.

## Spot-check (V&V)

| Check | Result |
|---|---|
| Hub RESULT / PROVENANCE | Present; tip, HDF5 1.14.5, CMake scratch patch, binary path, smoke case agree |
| Binary | `vv-runs/a14-hdf5-gnu/builds/b_hdf5_gnu/fds`; sha256 `c99737e643abb03b5ac49fcf14bddcbef43ac59ae5677ef9a0a9ae9344c538a2` |
| HDF5 install | `(local GNU third-party library tree)/libs/hdf5/v1.14.5` — `libhdf5.a`, `libhdf5_fortran.a`, `libhdf5_f90cstub.a` present |
| Enable method | Scratch CMake patch `patches/cmake-with-hdf5.diff` only; Makefile fallback not used |
| Smoke log | Banner `HDF5 library version: v1.14.5`; `STOP: FDS completed successfully` |
| GEOM `.vtkhdf` | `runs/tiny_vtk_cold/tiny_vtk_cold_GEOM.vtkhdf` — 17876 bytes; `file` = HDF5; magic `\x89HDF\r\n\x1a\n` |
| BNDF `.vtkhdf` | `..._BNDF_0000000000.vtkhdf` and `..._0010.vtkhdf` — 17276 bytes each; `file` = HDF5 |
| `paraview.py` | `tiny_vtk_cold_paraview.py` — 18183 bytes; loads GEOM/BNDF via `VTKHDFReader` |
| A-45 licence | `firemodels-gnu/hdf5/COPYING` (and install `share/COPYING`) on disk; BSD-like HDF Group text; no SPDX asserted here |

## Scope relative to owner show

A-14 proves FR-074 / FireX VTKHDF capability on a GNU Release, `USE_AMREX=OFF` smoke.
Owner show-and-tell critical-path viz stays **AMReX multilevel plotfile** (`WriteMultiLevelPlotfile`).
A-14 does not block the pack.

## Follow-ons (not this verify)

- Formal Spec close of D-019 / Q5 if still open beyond the GNU evidence.
- Optional `USE_AMREX=ON` + same HDF5 patch (called out in the hub RESULT).
- Intel HDF5 variant not in this hub.
