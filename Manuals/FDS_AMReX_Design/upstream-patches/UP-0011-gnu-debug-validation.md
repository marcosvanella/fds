# GNU Build Chief validation note: UP-0011 stock Debug target

Patch: `UP-0011-pressure-rhs-dump.patch`. Companion: `UP-0011-pressure-rhs-dump.signoff.md` (Debug-build revision section).

## Finding (makefile Debug, master `ce1f659cd4`)

On the GNU stock Debug target `ompi_gnu_linux_db` (`-Wall -Werror -fcheck=all` and the usual FP traps), an earlier apply of this patch failed to compile:

```
pres.f90:1098: Error: Unused variable 'fn' declared at (1) [-Werror=unused-variable]
```

The declaration was `CHARACTER(256) :: BASE,FN` in `PRESSURE_SOLVER_DUMP`. `FN` was never used. Log retained under `vv-runs/upstream-patch-validation/logs/build_pat11_ompi_gnu_linux_db.log`. Release (`ompi_gnu_linux`) built and passed the switch-off / switch-on checks (`vv-runs/upstream-patch-validation/RESULT.md`).

A scratch makefile that demoted that one warning (`-Wno-error=unused-variable`) built cleanly and the Debug runs matched the Release behaviour (bitwise with the switch off; dump files only with the switch on).

## Current on-disk patch

The patch file in this folder now declares only `CHARACTER(256) :: BASE` (no `FN`), matching the sign-off "Debug-build revision". That revision is the intended fix. This note does not rewrite the patch.

## Ask

Re-run stock `ompi_gnu_linux_db` on master (or FireX) with the current patch file and close this finding if the build is clean. Do not treat the demoted-warning scratch makefile as the official Debug result.

## Related

Dry-run on tip `8ad685f435`: applies (main.f90 offset +76; `pres.f90` clean). See `line-check-firex-1fb5c300.md` and `vv-runs/firex-merge-8ad685f435/logs/dry_apply_0010_up0011.log`.

## Re-run (stock `ompi_gnu_linux_db`) — CLOSED

Tree: master `ce1f659cd4` archive under `build/scratch-up0011-dbg`, patch applied with offsets −131 (main.f90) / −4 (pres.f90 hunk 2). Declaration in patched `pres.f90` is `CHARACTER(256) :: BASE` only. Docs note commit `05d78e2`.

| Check | Result |
|---|---|
| Stock `ompi_gnu_linux_db` (`-Wall -Werror -fcheck=all`, no demoted warnings) | **PASS** — builds and links; no unused-variable diagnostic |
| Smoke `ns2d_16` 1 rank, switch off | **PASS** — exit 0, `STOP: FDS completed successfully`, no Fortran runtime check |

Evidence: `vv-runs/up0011-stock-debug/logs/build_ompi_gnu_linux_db.log`, `vv-runs/up0011-stock-debug/runs/ns2d/`. Finding CLOSED.
