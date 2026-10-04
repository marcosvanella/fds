# Owner summary: D-074, the 23 FDS-only pressure inputs

Status: v0.4.38 draft, one page, for the project owner. Source: decision D-074 (`README.md`), requirements FR-003, FR-006 and FR-030, `vv/scope_case_list.csv` and `vv/case_inventory.csv`.

## What was decided
23 verification inputs have a pressure code of 0 on a non-periodic direction: the 21 `soborot_*` inputs plus `bound_test_1` and `bound_test_2`. In these inputs FDS skips the pressure solve. Under D-074 they are FDS-only in every mode, not only in AMR mode.

- The AMReX driver aborts on them with a message that names the FDS executable, so nobody gets a silent wrong answer.
- No no-pressure-solve mode is added to the driver. The recorded rationale is that such a mode would cost effort and have no test value, since no refinement or GPU feature needs it.
- FR-003 (uniform-mode non-regression) now covers "the verification set minus these 23" and the other FDS-only inputs, and the V&V Lead has moved them out of every denominator.

## What it changes in the numbers
- Verification scope list: 941 inputs in total. 682 are in scope and run in both modes (uniform and AMR). The rest are deferred (175), out of scope (49, which includes the 23) or unclear (35). The 23 are marked "FDS-only in every mode, no pressure solve" in the list.
- Inventory rows: 152 in scope of 187 in AMR mode.
- The exact split is in `vv/scope_case_list.csv` (columns `cls` and `amr_mode_status`). The 175 deferred inputs are features postponed by owner decisions (for example geometry cut cells, HT3D and cylindrical meshes), not failures.

## What you are asked
Nothing is blocked on you for D-074; it is recorded as accepted by the Chief Architect. Two things to know:
1. If you want these 23 to run in the AMReX code later, the work is a no-pressure-solve mode in the driver plus its tests. It is not in any phase now.
2. Related decisions that do need an answer from you are unchanged: the budget and date question (Q7), the FFTW licence notice (A-45) and the GPU parity tolerance class.

## Also carried in v0.4.38 (no action by you)
- D-076: an input with level-0 meshes of different resolution is FDS-only as written. A converter turns the finer meshes into refinement levels, and an abort guard (driver patch 0010) stops anything that slips through.
- D-077: the host/device flip-budget gate for the solid-phase solve is ratified.
