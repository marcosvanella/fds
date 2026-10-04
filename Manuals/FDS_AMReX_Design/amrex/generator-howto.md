# Generator how-to: add a loop, regenerate, check it bitwise

Owner: FDS Legacy Mapper. Audience: every owner on `docs/amrex/loop-work-list.md` who adds a loop to the S5 kernel generator. Paths are relative to the root of the FDS-AMReX checkout; commands run from the generator directory unless a path says otherwise. The commands and options below were checked against the script sources and `--help` output; no test suite was run to write this page.

Rules of the shared worktree (commit own files by path, `shared:` commits, the generator lock) are in section 4 of `loop-work-list.md`. The short form: take the lock around every full run, and never run a generator or test that writes into `generated/` or `build/` while another owner holds it.

## 1. Where the generator lives

Generator directory: `amrex/s4_mass/s5_gen/`.

| File or directory | What it is |
|---|---|
| `s5gen.py` | The generator. Reads the upstream FDS source of the pinned commit, parses the marked loops with `fparser`, writes the kernels. |
| `s5_fsrc.py` | Fortran source helpers (statement normalisation, anchors) used by `s5gen.py`. |
| `s5_markers.toml` | The sidecar: `[policy.*]` tables and one `[[kernel]]` entry per translated loop. The upstream source is never edited. |
| `markers/pvf_kernels.toml`, `markers/gsfv_field_kernels.toml` | Sidecars of the PATCH_VELOCITY_FLUX kernels and of the whole-field GET_SCALAR_FACE_VALUE kernels. Same entry format. |
| `markers/golden_signatures.json` | The golden file: one pinned signature per kernel (see section 6). |
| `red_argmax.py`, `red_extrema.py`, `edge_tables.py`, `s5_pvf.py`, `s5_gsfv.py`, `s5_gsfv_field.py` | Builders for special kernel kinds (argmax reductions, extrema reductions, edge tables, patch velocity, one-face and whole-field scalar face values). A kernel picks one with a flag in its entry (`argmax`, `extrema`, `edge`, `pvf`, `gsfv_field`). |
| `generated/` | Output: `s5gen_k2.F90`, `s5gen_k2.H`, `s5gen_report.md`, `s5gen_args.json`. Committed, so a regeneration shows up as a diff. |
| `port_merge_check.py`, `merge_check/` | The upstream-merge classifier and its replay driver (section 6). |
| `build_and_test.sh`, `test/` | Bitwise tests and mutation checks (section 5). |
| `README.md`, `README_pvf.md`, `README_gsfv.md` | The generator's own documentation. Read the README before adding a new kernel kind. |

Upstream line numbers in the sidecar refer to the pinned upstream commit (`policy.upstream_commit` in `s5_markers.toml`, the same FireX source as the survey `docs/inventory/gpu_generator_loop_classes.csv`).

## 2. Set-up

The generator needs Python 3 and `fparser` (0.2.5). A local wheel of it is kept in the repository under `src/`.

```
cd amrex/s4_mass/s5_gen
python3 -m venv venv
venv/bin/pip install fparser                 # or: venv/bin/pip install <path to the local fparser-0.2.5 wheel>
venv/bin/python s5gen.py --help
venv/bin/python port_merge_check.py --help
venv/bin/python merge_check/replay.py --help
```

The bitwise tests need `gfortran` (variable `FC`, default `gfortran`) and, for the make-scripts, a Python with `fparser` (variable `PY`, default `python3`).

## 3. Add a loop and regenerate

1. **Claim the loop** (protocol in `loop-work-list.md`). Read the loop in the pinned upstream source and the blocker family in `blocked-loop-families.md`.
2. **Add a `[[kernel]]` entry** to `s5_markers.toml` (or to your own sidecar file if you add a new kernel kind; see the READMEs). Fields:

   | Field | Meaning |
   |---|---|
   | `name` | Kernel name; becomes the C symbol suffix after `policy.prefix` (`s5gen_`). Unique. |
   | `file`, `routine` | Upstream source file and routine. |
   | `lines = [first, last]` | Statement range in the pinned commit: one K,J,I nest, whole-array assignments, or a wall or index loop. |
   | `anchor_line`, `anchor` | The first statement, whitespace-normalised, and its line. If the anchor is not at `anchor_line`, the generator relocates by a unique anchor match and warns; otherwise it fails. This is the drift check. |
   | `args = [...]` | The C-visible argument order the generator must derive. The generator fails (exit code 2) when its derived argument set differs. Omit it only for kernel kinds whose builder derives the order (see the READMEs). |
   | `layout` | `"exact"` for arrays whose allocated bounds are read from the `ALLOCATE` text (wall and exact-shape kernels); default is the FDS index box. |
   | `private = [...]` | Scalars assigned before they are read in every pass and not read after the loop. |
   | `unique = true` | Every write goes to a target no other iteration writes (wall-subscript stores). The generator checks what it can; the proof is yours. |
   | `idempotent = true` | Several iterations may write the same target, but always the same bits (the stored value depends only on the target). Justify it from the source in the entry comment. |
   | builder flags | `argmax`, `extrema`, `edge`, `pvf`, `gsfv_field`: pick the builder. |

   New module tables, component tables, wall tables or pointer targets go in the matching `[policy.*]` table (for example `[policy.tables]`, `[policy.wall]`, `[policy.assoc]`), each with the `ALLOCATE` line that backs it.
3. **Regenerate.** From `amrex/s4_mass/s5_gen`:

   ```
   venv/bin/python s5gen.py                       # uses the pinned commit via git (default --repo is the FDS-AMReX checkout)
   venv/bin/python s5gen.py --src-dir <dir>       # or read the .f90 files of a directory instead of a git commit
   venv/bin/python s5gen.py --callee-directive bind   # default switch for regions that call a routine: dpd (default) or bind
   ```

   It writes `generated/s5gen_k2.F90`, `generated/s5gen_k2.H`, `generated/s5gen_report.md` and `generated/s5gen_args.json`. Exit codes: `0` success; `2` a loop cannot be translated, or the argument set differs from the pinned `args`; `3` the golden signature differs; `4` a parse error. Every failure names the kernel and the reason.
4. **Pin the new signature** only after you have reviewed the kernel: `venv/bin/python s5gen.py --update-golden`. Never use it to silence an exit code 3 you do not understand (section 6).
5. **Regenerate one family only.** The families with their own builder and test have their own sidecar and test driver, so you can work on one without touching the others:

   | Family | Sidecar | Generate and test |
   |---|---|---|
   | main K2 kernels | `s5_markers.toml` | `venv/bin/python s5gen.py`, then `build_and_test.sh` |
   | PATCH_VELOCITY_FLUX | `markers/pvf_kernels.toml` | `test/make_pvf_tests.py --out <dir>`, `test/run_pvf.sh`, `test/test_mesh_pvf.py` |
   | one-face scalar face value | (builder `s5_gsfv.py`) | `test/make_gsfv_tests.py --out <dir>`, `test/run_gsfv_bitwise.sh`, `test/test_mesh_gsfv.py` |
   | whole-field scalar face value | `markers/gsfv_field_kernels.toml` | `test/make_gsfv_field_tests.py --out <dir>`, `test/run_gsfv_field.sh`, `test/test_mesh_gsfv_field.py` |

   `--help` of the make-scripts lists their options; `test/make_r2_tests.py [--src-dir DIR]` writes `s5_r2_ref.F90` and `s5_r2_bitwise.F90` from verbatim upstream line ranges.

## 4. What "bitwise vs upstream" means

The reference of every test is the **verbatim upstream loop text**, copied by line range from the pinned source by the `make_*` script (never retyped). The generated kernel and the reference run on the same test inputs (the cases are built to reach the `CYCLE` branches, for example internal walls on the mesh boundary); the outputs must be equal bit for bit. The build uses `-ffp-contract=off`, no `-ffast-math`, and gfortran's default evaluation order, so the comparison does not depend on fused multiply-add.

A loop is accepted when all of these hold:
- all six flag sets: `O0`, `O2`, `O0omp`, `O2omp`, `O2omp_off` (host fallback of the offload branch), `O2omp_dpd` (forced `distribute parallel do`);
- both callee switches: `dpd` and `bind` (variable `CALLEES`, default `"dpd bind"`);
- 1, 4 and 8 threads for the OpenMP sets (the `run_*` drivers do this inside the program);
- the mutants are caught (section 5) and the negative checks still reject what they reject;
- reductions: sums keep the FDS order unless the exact fixed-point switch is on (Decision B in `generator-decisions.md`, GPU Generator Engineer).

## 5. Test entry points

Run from `amrex/s4_mass/s5_gen`. Take the generator lock around full runs; all of them write below `build/` or `generated/`:

```
flock ../../../.s5gen.lock ./build_and_test.sh all
```

(`../../../.s5gen.lock` is `.s5gen.lock` at the repository root; the `run_*.sh` headers show the same form as `flock <repo>/.s5gen.lock test/run_pvf.sh all`.)

| Command | What it does |
|---|---|
| `./build_and_test.sh [flag-set ...]` | Builds the S4d reference, the hand-made K2, the generated K2 and the round 1 and 2 bitwise drivers with gfortran for each flag set (default `all` = the six sets) and each callee switch, runs the driver. Variables: `FC`, `CALLEES`, `OUT`. Argument `O0` alone is the quick one. |
| `test/run_pvf.sh [flag-set ...]` | PATCH_VELOCITY_FLUX bitwise test at 1, 4, 8 threads. Variables: `FC`, `OUT`, `PY`, `SKIP_MAKE=1`, `MUT=<file>` (replacement kernel file for a mutant), `XDEF` (extra `-D` flags, e.g. host-rule mutants). |
| `test/run_gsfv_bitwise.sh [flag-set ...]` | One-face scalar-face-value bitwise test. Variables: `FC`, `OUT`, `PY`. |
| `test/run_gsfv_field.sh [flag-set ...]` | Whole-field scalar-face-value bitwise test. Variables: `FC`, `OUT`, `PY`, `SKIP_MAKE`, `MUT`. |
| `test/mutation_check.sh` | Each mutant edits one statement of a copy of the generated module; the driver must report failing cases. Run after `./build_and_test.sh O0`; override the build directory with `B=<dir>`. It writes a scratch file under `generated/`, so take the lock. |
| `python3 test/test_generator_checks.py` | The negative checks of the generator (loops it must refuse: loop-carried scalars, reads before writes, possible aliasing, stencil reach beyond the box, derived-type designators and so on). It takes no arguments and runs when started, so do not start it just to see its options. |
| `python3 test/test_merge_check.py` | Tests of `port_merge_check.py` on edited scratch copies of the pinned source (cases the text heuristic misses, body-only controls, CLI exit codes, golden signatures). |
| `python3 test/test_edge_checks.py`, `test_argmax_checks.py`, `test_extrema_checks.py` | Negative tests of the edge-table, argmax and extrema builders: each edits a scratch copy of the pinned source and the generator must stop with the stated exit code. |
| `python3 test/test_mesh_pvf.py [--build\|--mutants\|--all\|--write-golden]` | No option: fast checks (generation, kernel text against the family golden file, argument list and shapes, stencils, signatures, refusals on edited copies of the upstream text). `--build` adds the Fortran bitwise test in all six flag sets (1, 4, 8 threads, under the generator lock); `--mutants` adds mutation of the generated kernels (the bitwise test must fail for every mutant); `--all` is everything; `--write-golden` rewrites the family golden after a reviewed, intended change. Same options in `test_mesh_gsfv.py` and `test_mesh_gsfv_field.py`. |

Run the single test of your own kernel first, then the full set under the lock before you commit.

## 6. What golden and port_merge_check do

**Golden** (`markers/golden_signatures.json`, plus the `test/*.golden` files of the builders). The generator writes, for every kernel, its full C-visible signature: each argument with name, shape and bounds, type and intent (`IBAR:scalar:integer:in`, `ADV_FX:array(0:IBAR+1,...):real:inout`, ...). On every run `s5gen.py` compares what it derived with the pinned signature and stops with exit code 3 when they differ. That catches an unintended change of a kernel's interface (a new read, an intent that changed, a bound that moved) that would break the C++ shim. After a reviewed, intended change run `s5gen.py --update-golden` and commit the new golden in the same commit as the sidecar entry. The `args` list of a `[[kernel]]` entry pins the argument order the same way, with exit code 2.

**port_merge_check** (`port_merge_check.py OLD NEW`). When upstream FDS moves from one commit to another, this tool tells which upstream hunks touch the translated loops. It builds every marked kernel in OLD and in NEW (reusing the generator, so every generator rule is also a classifier rule), diffs the upstream text, and classifies each hunk against the manifest into one of four classes: `body-only` (same manifest, only the body or a callee body differs: paste through the rename map, regenerate, run the bitwise gate), `surface-changing` (the manifest differs or the loop is no longer translatable: new argument, new private scalar, bound, offset, reordered access, new guard, changed callee interface; the sidecar and the golden need review), `outside-kernel` (hunk in a marked routine but outside every marked range, manifest unchanged) and `unmarked-text-heuristic` (hunk in a candidate routine with no marker). Exit codes: `0` nothing needs action in a marked kernel, `1` a surface-changing hunk touches a marked kernel, `2` tool error.

```
venv/bin/python port_merge_check.py <old ref or dir> <new ref or dir> [--repo R] [--markers FILE] [--inline FILE ...] \
    [--candidates port_kernel_map.csv] [--golden markers/golden_signatures.json] [--all-kernels] [--format json|csv] [-o OUT]
venv/bin/python merge_check/replay.py [--repo R] [--since S] [--limit N] [--jobs J] [--out DIR]
```

`merge_check/replay.py` replays the upstream history through the same classifier to measure how often a merge touches a marked kernel; its results are kept in `merge_check/results/`. Run `port_merge_check.py` before every upstream merge into the working branch and after adding kernels in a hot region (DIVERGENCE_PART_1, VELOCITY_FLUX); a `1` means stop and review.

## 7. Refresh the work list after you commit

```
python3 tools/inventory/loop_work_list.py            # regenerate docs/amrex/loop_work_list.csv and the generated tables
python3 tools/inventory/loop_work_list.py --check    # exit 1 when a committed output is stale
```

The script reads the committed sidecars at the generator branch head with `git show`, so a loop turns to `translated` only once your kernel entry is committed. Uncommitted edits in the shared worktree are never counted.
