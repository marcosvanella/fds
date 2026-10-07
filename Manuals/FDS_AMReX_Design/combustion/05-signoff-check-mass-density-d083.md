# 05 — D-051 device sign-off: CHECK_MASS_DENSITY two-pass K2 (D-083)

Owner: AMR Species & Combustion Lead · Status: ACCEPT-with-conditions · Kernel packet: `s5_clip_gather` (L0876 / L0877).

## Verdict

**ACCEPT-with-conditions** for the GPU-only two-pass hand-written K2 draft of `CHECK_MASS_DENSITY` under D-083. CPU `CHECK_MASS_DENSITY` stays unchanged; no upstream FDS edit (`UP-0003` remains tabled).

Conditions:
1. Kernel `nclip` OUT is a **0/1 any-clip flag** (`reduction(max:)`), not exact `NCLIP` / `NLO` / `NHI`. Callers must treat it as presence only. Exact counts stay on the CPU CLIP flags and the P1-oracle counters (as in the packet notes). Acceptable for the device gating path; not a reject.
2. **Single-box** CPU driver in this packet. Multi-box `FillBoundary` of clipped `RHOP` between density apply and species gather is **out of scope here** (P1 driver / Species Lead cmd2p). Must land with the multi-box driver before claiming multi-box equivalence.
3. Do **not** mark ported until a device run on real fields is on record and the Integration Lead adds the `ported.toml` entry (see Recommendation).

## Checklist vs `combustion/02` (clip gather prototype conditions)

| # | Check (from `02` clip-gather / UP-0003 record) | Result |
|---|---|---|
| 1 | Scatter / gather order: sources `(k-1)`, `(j-1)`, `(i-1)`, self, `(i+1)`, `(j+1)`, `(k+1)`; `OD` = 3,2,1,0,−1,−2,−3 (FDS K,J,I delivery) | PASS — `OI`/`OJ`/`OK`/`OD` match `clip_gather.f90` and `blocked-loop-families` S1 |
| 2 | Rounding: explicit left-to-right `SUM(MASS_N)` over −3..3 with parentheses; FMA off in build | PASS — `sum_mass_n = (((((mass_n(-3)+…)+mass_n(3)`; `-ffp-contract=off` in `run_clip.sh` |
| 3 | Flags / threading: integer / max reductions for clip flags; no `reduction(+)` of domain flags; host OR between passes | PASS — local `reduction(max:)` then host `pack_or` / `or_into`; apply/renorm gated |
| 4 | `CLIP_TERMS` + face mask: species SOLID-first / density range-first; `MASK(:,:,:,0:7)` = SOLID, faces −1,+1,−2,+2,−3,+3, clip-source | PASS — matches `clip_gather.f90` and FireX `WALL_INDEX` / `SOLID` order |

`combustion/03` covered/outside scatter targets: this packet’s `MASK(:,:,:,7)` enforces single-mesh clip-source membership (layout-independent gather). Per-level “outside/covered target counted, not redistributed” counters from `03` are **not** in this packet; they remain a driver / AMR accounting item, not a kernel reject.

## Evidence

| Item | Result |
|------|--------|
| Kernel commit | `8edd573d331f47df660cc05361b4295d74e3ec58` (`s5-cmd-clip` from s5-gen `7062ec2dce`) |
| FireX reference | `36975d765f` `mass.f90:775-963` (sha256 prefix `8bbacbbe…`); density L0876, species L0877 |
| Prototype | `prototypes/p1_mass_shim/clip_gather.f90` (`CLIP_TERMS` + gather) |
| Host bitwise (gfortran, FMA off) | Flag sets `O0 O2 O0omp O2omp O2omp_off O2omp_dpd` @ 1/4/8 threads: P1-oracle **150** cases and CPU **540** cases PASS (author build; O2omp re-run confirmed) |
| Mutants | **13/13** killed (re-run) |
| `k2_ci_check.py` | PASS |
| `kernel_lint.py --clauses-only` | PASS |
| Upstream FDS edits | None in the commit (packet under `s5_clip_gather/` only) |
| `ported.toml` | Not entered (correct) |

### Mutant list (13)

| Mutant | Role |
|--------|------|
| `gather_order_self_first` | wrong gather order |
| `gather_order_reversed` | wrong gather order |
| `sum_mass_n_regroup` | wrong `SUM(MASS_N)` association |
| `species_apply_ungated` | apply without host flag |
| `renorm_ungated` | renorm without host flag |
| `species_one_ungated` | one-species path without density flag |
| `flag_from_ghost` | density clip flag from non-valid (ghost) |
| `species_flag_from_ghost` | species clip flag from non-valid |
| `species_qmax_self` | species `QMAX` from target instead of source |
| `source_mask_ignored` | dropped neighbour / ignored clip-source mask |
| `face_mask_ignored_xm` | dropped neighbour / ignored face mask |
| `species_apply_no_ceiling` | apply without `min(RHOP,…)` ceiling |
| `pass_order_species_before_apply` | wrong pass order (driver) |

Wrong sum order and dropped-neighbour mutants are present (`sum_mass_n_regroup`, `source_mask_ignored`, `face_mask_ignored_xm`). Exact counter-race mutants are **not** present: the device path does not expose exact `NCLIP`/`NLO`/`NHI` (0/1 `max` only), so a `reduction(+)` race on counts cannot fire there; ghost-flag mutants cover wrong flag sources. Acceptable under condition 1.

## Body vs FireX / P1 (spot checks)

- Density gather/apply and species gather/apply expressions match `mass.f90:799-854` / `:858-937` via the gather rewrite (same `CONST*MASS_N/VC` terms, gated apply).
- Renorm matches `mass.f90:947-961` with hand-written left-associated species sum and first-occurrence max (`>`), gated by reduced density OR any species flag.
- Host `s5clip_build_mask_box` maps `SOLID` / `WALL_INDEX(-3:3)` / source → `MASK(:,:,:,0:7)` as in ADR-001 / D-031.

## Open conditions (summary)

1. `nclip` is 0/1 presence only; exact counts via CPU / P1-oracle.
2. Multi-box `FillBoundary` of clipped `RHOP` deferred to the driver packet.
3. `ported.toml` entry: Species Lead clears for Integration Lead to set (device evidence below). Conditions 1–2 still in force.

## Device evidence

**PASS** — evidence commit `b582e79027cd8d7bae3e5c5656a202ccbf122c8e` on branch `s5-cmd-clip` (`DEVICE.md`; tip at evidence `f5cbeec5e8`). Kernel commit unchanged: `8edd573d331f47df660cc05361b4295d74e3ec58`.

| Item | Value |
|------|--------|
| Flags | `nvfortran -O2 -Minline -mp=gpu -gpu=cc89,nofma,mem:managed -Mnofma -Mnovect -Mextend -Mnoflushz -Mnodaz -Mpreprocess`; `OMP_TARGET_OFFLOAD=MANDATORY` |
| Threads | 1 and 8 |
| P1-oracle | **150** PASS (rand / negzero / ring) at both thread counts |
| CPU FireX path | **540** PASS (rand / dyad / open / ns1 / nsN / noclip) at both thread counts |
| Coverage | `density_pass_runs=487`, `renorm_runs=440` |
| Offload | `-Minfo=mp`: GPU kernels generated for `s5clip_density_gather/apply`, `s5clip_species_gather/apply/one`, `s5clip_renorm` |
| Conditions 1–2 | Still in force (`nclip` 0/1 presence; multi-box FillBoundary not claimed) |

Spot-check: `device_run.log` and `minfo_k2.log` under the Mesh Data scratch match `DEVICE.md` (no invented numbers; GPU kernels present).

## Recommendation for `ported.toml`

**Species Lead clears** for Integration Lead to add the `[[ported]]` entry (Species Lead = `reviewed_by`; Integration Lead = `set_by`). Conditions 1–2 remain in force and must stay in the entry note. Also record the device run in `tools/device_runs.toml` (PM-08) if not already present. Exact entry text: `work-species/handoff-s4/cmd_ported_clearance.md`. Do **not** claim multi-box FillBoundary.

## Notify

Integration Lead: set `ported.toml` (+ `device_runs.toml` if needed) from the clearance handoff. Architect: this file remains the record. Mesh Data Loops Engineer: see `work-species/handoff-s4/cmd_ported_clearance.md`.

## Caller-contract mirror

The ACCEPT-with-conditions caller rules are also stated in `amrex/s4_mass/s5_gen/s5_clip_gather/README.md` (Caller contract) and `REVIEW.md` on branch `s5-cmd-clip`, so drivers do not have to re-read this sign-off alone.
