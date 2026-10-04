# Waiver list: FDS read-before-set sites (trap-free NFR, spec v0.4.6)

Owner: AMR V&V Lead. Base: FireX 36975d765f (unmodified).
Scope: the Intel checked-debug build with `-init=snan,arrays` (plus integer init where noted). New AMR code must not trap. Only the existing FDS sites listed here are waived.
How the check runs: use `vv-runs/refbin/impi_intel_firex-36975d7/diag/build-intel-trapfree-diag.sh --smoke <src-copy> <out> <patch>`, or by hand: apply the cumulative patch covering only the waived sites (`vv-runs/refbin/impi_intel_firex-36975d7/diag/cumulative.patch`) to the milestone source in a scratch copy, build with the full -init flags, and run the FR-016 gating case and the FR-005 (v) cases. Any trap outside this list fails the check.

| # | Site (file:line) | Variable | Trigger | Affects results? | Found | Evidence |
|---|---|---|---|---|---|---|
| W1 | read.f90:1846 (alloc), read.f90:17450 (read) | H_V_H2O(0) | snan,arrays | Probably not (element 0 is read only for temperatures below 1 K) | 2026-09-25 | impi_intel_db_chk/provenance/acceptance/SNAN_TRAP_FINDING.md |
| W2 | cons.f90:640, read.f90:8024 | N_MATL, used to size CHILD_LAYER/CHILD_SURF before READ_MATL sets it | integer huge | No (only an oversized allocation) | 2026-09-25 | impi_intel_db_chk_noarr/provenance/acceptance/HUGE_INT_SEGV_FINDING.md |
| W3 | radi.f90:2779 / 3999 / 5190 | N_RADCAL_ARRAY_SIZE (with huge, the `>0` test uses an unallocated KAPPA_COND; the loop bound at 5190 reads it too) | integer huge | No | 2026-09-25 | impi_intel_db_chk/provenance/acceptance/SNAN_TRAP_FINDING.md |

All three are candidates for upstream issues under A-20. Only 2 smoke cases have been swept so far, so there may be more sites.

Patch sha256: `7d75cc2307c0af1498b547b3b6715a27d8204d34e64222e553a19f3a94cd6a7a` (its 3 hunks map to W1 to W3). The patch and the diagnostic build were **reconstructed** by the Intel Build Chief after an earlier environment reset; they are not byte-identical to the lost original (`54789a51…ca6c5b`, diagnostic fds `4ffc71c3…13c4`). Location: `vv-runs/refbin/impi_intel_firex-36975d7/diag/` (`cumulative.patch`, `build-intel-trapfree-diag.sh`, `impi_intel_db_trapfree/fds`, sha256 `0191d22e1d3cf07fa81acd98911a3c08925b68dc3b63a54c7e930d30149ff8c6`, `-init=snan,arrays -init=huge`, version string `FDS-6.11.1-1244-g36975d765f-trapfree-diag`). The evidence files named in the Evidence column (`SNAN_TRAP_FINDING.md`, `HUGE_INT_SEGV_FINDING.md`) were lost in the reset; the W3 lines were re-checked against the 36975d7 source (radi.f90:5190 is the `DO N = 1, N_RADCAL_ARRAY_SIZE` loop; the `Z_IN = Z_IN * KAPPA_COND` line is 5181).
Requirement: NFR-046 (requirements v0.4.7).
