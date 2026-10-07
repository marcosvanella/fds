# D-072 (c) co-sign: 17 `rad_*` kernels (synthetic hash meshes)

Owner: V&V Lead (reviewer of record under D-077). Date: 2026-10-06.
Full write-up: `docs/vv/ported-review.md` section 10. Proposed rows: `docs/tools/patches/ported-rad-proposed.toml` (docs commit `e3d241d`, co-signed wording in a follow-up docs commit). Tip: s5-gen `25aad430b9`.

## Verdict

**ACCEPT-WITH-LIMITS** — co-sign for `[[ported]]` land with the scope-limited note below. Integration may land `[[device_run]]` immediately. D-085 scheme-3 is separate (already ratified).

## Checklist

| # | Condition | Result |
|---|---|---|
| 1 | Device run on record within class tolerance (bitwise) | **PASS** on SYNTHETIC hash meshes (smoke 8×6×4 + big 24×20×16, dflt+dpd, 17/17 EQUAL, IN match, `check_device_logs` PASS). **Not** real-field / FDS-truth overlay. |
| 2 | Kernel owner reviewed | **PASS** — AMR Radiation Lead |
| 3 | `ci_checks.sh --strict` | **PASS** on `src-s5gen` `25aad430b9` (ZS-05 waived approved; prior `rad_wall_qin_zero` / `W_SURF_INDEX` FAIL cleared) |

## Evidence verified by V&V

- Tarball `scratch/devrad_logs_25aad.tgz` (extract: `scratch/devrad_logs_extract/logs/`). Summary: `scratch/rad_device_summary_25aad.txt`.
- Spot-check `rad_gray_kfst4` big REF/GEN `2F3F57E4408AB2D5 87CE5BCCF856BF5F`.
- 17 names match both proposed blocks and all four bitwise files.
- Stale pre-fix extract (`dcada34a2e` DIFF run) was replaced; do not use the old summary/run_*.out content.

## Recommended `[[ported]]` fields (Integration owns the live file)

```
reviewed_by = "AMR Radiation Lead (device hash green at s5-gen 25aad430b9); V&V Lead co-sign D-072 (c) ACCEPT-WITH-LIMITS (docs/vv/ported-review.md section 10; docs/vv/rad_device_d072c.md)"
set_by = "AMReX Integration Lead"
note = "Bitwise class. Device hash driver gpu/mkdev_rad.py (s5-gen 25aad430b9): device GEN == gfortran REF and IN equal (dflt and dpd) on SYNTHETIC smoke 8x6x4 and big 24x20x16 (RTX 4070, nvfortran 26.9, cc89,nofma,mem:managed); check_device_logs --logs-only PASS; rad_gray_kfst4 integer TMP**4 bitwise EQUAL (big REF/GEN 2F3F57E4408AB2D5 87CE5BCCF856BF5F). Host smoke REF=GEN 17/17; ci_checks --strict production PASS (ZS-05 s5rad_wall1 approved waiver; prior rad_wall_qin_zero W_SURF_INDEX FAIL cleared). SCOPE LIMIT: SYNTHETIC hash meshes only — same evidence class as gate Form A / rounds 4-7 hash coverage; NOT claiming real-field / FDS-truth overlay (open follow-up). D-085 scheme-3 is separate."
```

## Open follow-up (not a blocker for this co-sign)

Real FDS case overlay (device = host on dumped radiation fields) would close the literal D-072 (c) "real fields" reading; until then the ported note must keep the SYNTHETIC scope limit.
