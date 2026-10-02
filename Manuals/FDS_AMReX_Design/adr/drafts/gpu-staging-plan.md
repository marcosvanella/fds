# Staged plan: moving the FDS time step onto the GPU

Draft for the project owner, plain language. Time shares come from `docs/inventory/gpu_generator_coverage.md` (FireX `36975d7`) and are the survey's **modelled** shares of CPU run time, not measurements. Every speed-up below is an **estimate** (method in section 4).

**Terms.** *Device / host*: the GPU and its memory / the CPU and its memory; copying between them is slow, so fields should stay on the device for the whole step ("resident"). *Kernel*: one loop nest run on the GPU. *Generated kernel*: a kernel produced by our generator from the upstream FDS loop text, so it follows upstream changes. *FillBoundary*: AMReX's routine that fills ghost (halo) cells from neighbouring boxes and ranks; it replaces FDS's `MESH_EXCHANGE`. *Wall table*: one flat array per wall quantity, replacing the `WALL(IW)%...` records.

## 1. Where we are
34 generated kernels (about 1.2% of modelled run time) are bitwise-tested against the upstream loop text on the CPU only. The hand-ported mass and divergence prototypes ran on the GPU test machine and matched the CPU bit for bit. The driver (M2a) runs the full step for periodic cases on the CPU, with pressure from the AMReX FFT solver. Nothing runs a whole step on the GPU yet.

## 2. What goes to the device
| Modelled share | What | Plan |
|---|---|---|
| 8.2% | 147 plain cell-array loops (classes A and B: mass, density, divergence part 1, viscosity) | GPU, Stage 1 |
| 23.9% | 32 wall loops that need only flat wall tables (tier T1; includes the large `DIVERGENCE_PART_1`, `MASS_FINITE_DIFFERENCES`, `SPECIES_ADVECTION` nests) | GPU, Stage 2 |
| 8.2% | 14 wall loops that call other routines, mainly `WALL_BC` (tier T3) | GPU, Stage 3 (hand port, generated arithmetic) |
| 8.9% | 37 loops that read neighbouring meshes or per-wall layer data (tiers T4, T2) | GPU, Stage 4, after the exchange-buffer layout exists |
| 16.5% + 1.8% | `MESH_EXCHANGE` (its timer includes MPI waiting, so it overstates compute) and the `pois.f90` pencil loops | Not ported: replaced by FillBoundary and the AMReX FFT/multigrid solvers, which run on the device |
| 23.9% | Cut-cell, cut-face and CFACE loops, `pres.f90` matrix assembly, `MATCH_VELOCITY` (geometry) | Deferred to the geometry refactor (ADR-003); host |
| 1.9% | Particles, control flow, file I/O | Host |
| about 6.7% | Time the 829-loop model assigns to no loop | Assumed host |

## 3. Stages
**Stage 1: whole step on the device for wall-free periodic cases.** Cell kernels, FillBoundary on the device (CUDA-aware MPI), FFT pressure on the device, `dt` and the exact sums as device reductions. *Resident:* every cell field and work array, exact-sum buffers, FFT plans. *Transfers:* output and diagnostics only, plus a few scalars per step. *Exit test:* bitwise equal to the CPU kernels on frozen input, and the M2a cases unchanged.

**Stage 2: simple walls.** Per-mesh wall tables with per-box index lists (ADR-001 ruling W1), the 32 T1 wall loops. *Resident:* the wall tables the kernels write. *Transfers:* `UVW_SAVE`, `U_GHOST`, `V_GHOST`, `W_GHOST` go host-to-device once per step because their producers stay on the host (ruling W2; one real per external wall cell each); tables are re-gathered after obstruction events; **`WALL_BC` still runs on the host, so the wall tables it reads and writes round-trip every step**. That round trip is the main extra cost of this stage, and it scales with the number of wall cells, not cells.

**Stage 3: `WALL_BC` on the device.** Removes the Stage 2 round trip. Port after the upstream interface settles (`SURFACE_HEAT_TRANSFER` changes often).

**Stage 4: neighbour-mesh and per-wall layer data.** `VELOCITY_BC`, `NO_FLUX`, `HT3D` exchange and the other T4 loops need the exchange-buffer layout (what `OMESH(NOM)%...` becomes) and flat tables for the 1-D solid conduction layers.

**Beyond this plan (ADR-003).** Cut cells move to flat form, `UVW_SAVE` and `U_GHOST` are produced on the device and the Stage 2 upload is deleted. The regrid-time rebuild runs on the host until Phase 11 (D-047).

**Transfers that remain by design:** output and restart, host-only loops, geometry routines still on the host, and MPI between nodes when it is not CUDA-aware.

## 4. First real speed-up, and how it is estimated
New time as a percentage of the old CPU time = host-only share + device share / S + 2, where the retired 18.3% (`MESH_EXCHANGE`, `pois`) is assumed to cost 2 points after FillBoundary and the FFT solver replace it. S is the kernel speed-up against one FDS CPU process. **S is not measured**; 10 and 30 are assumed. Against a multi-core CPU node the real figure is lower.

| After stage | Device share | Estimate, S = 10 | S = 30 |
|---|---|---|---|
| 1 | 8.2% | 1.3x | 1.3x |
| 2 | 32.1% | 1.8x | 1.9x |
| 3 | 40.3% | 2.1x | 2.2x |
| 4 | 49.2% | 2.5x | 2.8x |

The ceiling with geometry on the host is about 2.9x. The survey mix is not any one case: a periodic wall-free case has almost no wall or geometry time, so its Stage 1 share is higher than 8.2% [VERIFY with the per-case timers in `gpu_timing_cases.csv`].

- **First end-to-end speed-up:** expected after Stage 1 on a large periodic case (128^3 or more). At `shunn3_32` size (32^3) launch latency dominates and no speed-up is expected. For typical wall-bounded cases the first estimate above 1.5x is Stage 2.
- **To replace the estimate with measurements** (not done): time the generated kernels against the CPU loops on one core and one GPU to get S per kernel group; run the per-case timers on a periodic 128^3 case; measure the Stage 2 round-trip volume on a wall case.
