# Waiver list: FDS read-before-set sites (trap-free NFR, spec v0.4.6)

Owner: AMR V&V Lead. Base: FireX 36975d765f (unmodified).
Scope: the Intel checked-debug build with `-init=snan,arrays` (plus integer init where noted). New AMR code must not trap. Only the existing FDS sites listed here are waived.
How the check runs: use `build/env/build-intel-trapfree-diag.sh --smoke <src-copy> <out> <patch>`, or by hand: apply the cumulative patch covering only the waived sites (`build/firex-36975d7/diag_intel/init_arrays_patched/provenance/cumulative.patch`) to the milestone source in a scratch copy, build with the full -init flags, and run the FR-016 gating case and the FR-005 (v) cases. Any trap outside this list fails the check.

| # | Site (file:line) | Variable | Trigger | Affects results? | Found | Evidence |
|---|---|---|---|---|---|---|
| W1 | read.f90:1846 (alloc), read.f90:17450 (read) | H_V_H2O(0) | snan,arrays | Probably not (element 0 is read only for temperatures below 1 K) | 2026-09-25 | impi_intel_db_chk/provenance/acceptance/SNAN_TRAP_FINDING.md |
| W2 | cons.f90:640, read.f90:8024 | N_MATL, used to size CHILD_LAYER/CHILD_SURF before READ_MATL sets it | integer huge | No (only an oversized allocation) | 2026-09-25 | impi_intel_db_chk_noarr/provenance/acceptance/HUGE_INT_SEGV_FINDING.md |
| W3 | radi.f90:2779 / 3999 / 5181 | N_RADCAL_ARRAY_SIZE (with huge, the `>0` test uses an unallocated KAPPA_COND) | integer huge | No | 2026-09-25 | impi_intel_db_chk/provenance/acceptance/SNAN_TRAP_FINDING.md |

All three are candidates for upstream issues under A-20. Only 2 smoke cases have been swept so far, so there may be more sites.

Patch sha256: 54789a51eae569819a4a1add3ca8fef723b37d03aed9e025327afabb5bca6c5b (its 3 hunks map to W1 to W3, see build/env/README-intel.md).
Requirement: NFR-046 (requirements v0.4.7).
