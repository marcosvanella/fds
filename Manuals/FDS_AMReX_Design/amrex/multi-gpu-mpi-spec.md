# Multi-CPU / multi-GPU MPI execution: spec input (owner requirement Q6)

| | |
|---|---|
| Author | AMReX Integration Lead |
| Status | **DRAFT** (input to spec v0.4.11; requirement and action numbers are placeholders `FR-0xx`/`NFR-0xx`/`A-0xx` for the Spec Lead) |
| Date | 2026-09-25 (rev. 2: §5 added for A-43, several CPU cores per GPU; §3, §6 and §7 tables adjusted and §8 extended; old §5-§7 renumbered to §6-§8) |
| AMReX base | `(local AMReX checkout)` @ `99ddfda` (26.09), read-only |
| Companion docs | `driver-options.md` §5 (GPU implications), `docs/adr/ADR-001-driver-architecture.md` (IR-007 kernel rules), `docs/requirements.md` §2.2 (tolerance classes T0-T3), `docs/risks.md` R-31, `(local NVIDIA HPC SDK directory)/README.md` |

AMReX paths are relative to `(local AMReX checkout)/`; `GPU.rst`, `RuntimeParameters.rst` and `Faq.rst` are in `Docs/sphinx_documentation/source/`.
Nothing here was built or run. Claims about AMReX come from reading source at `99ddfda`. Claims about MPS, Open MPI and the HPC SDK come from the vendor pages cited. **Assumption** marks anything not verified from source or docs.

**Requirement (Q6, the project owner).** The code runs on multiple CPUs and multiple NVIDIA GPUs with MPI, across nodes. It uses at least one MPI rank per GPU, and possibly several ranks sharing one GPU. The only GPU hardware coming online is an owner-provided NVIDIA test machine with one GPU. No cluster has been named.

---

## 1. How AMReX maps MPI ranks to GPUs

**Initialization order** (`Src/Base/AMReX.cpp`):
1. MPI start: `StartParallel` `:403`.
2. Device selection: `Gpu::Device::Initialize(minimal, a_device_id)` `:562-565`.
3. `ParallelDescriptor::Initialize()` `:693`, which fixes `use_gpu_aware_mpi` (§4).
4. `Arena::Initialize` `:696`.

As a result, both the device choice and the count of ranks sharing a device are known before the arenas are sized.

**Node-local rank** (`Src/Base/AMReX_ParallelDescriptor.cpp:346-356`). When NProcs > 1, AMReX calls `MPI_Comm_split_type`. With Open MPI the split type is `OMPI_COMM_TYPE_NODE`; with other MPIs it is `MPI_COMM_TYPE_SHARED`. The result gives `NProcsPerNode()`/`MyRankInNode()`. A second count, `NProcsPerProcessor()`/`MyRankInProcessor()`, groups ranks by `MPI_Get_processor_name` (`:358-381`). Both are declared in `AMReX_ParallelDescriptor.H:226-243`.

**Device selection rule** (`Src/Base/AMReX_GpuDevice.cpp:283-303`). `gpu_device_count` comes from `cudaGetDeviceCount` (`:269-273`), so it counts only the devices that `CUDA_VISIBLE_DEVICES` leaves visible. If it is 0, AMReX aborts with "No GPU device found" (`:271-272`). The rules below are checked in order:

| Condition | `device_id` | Line |
|---|---|---|
| `Init_minimal` path | whatever device is current (`cudaGetDevice`) | `:283-286` |
| caller passed `a_device_id >= 0` to `amrex::Initialize` | `a_device_id` | `:287-288` (API: `AMReX.H:106,113`, default `-1`) |
| 1 MPI rank | 0 | `:289-290` |
| 1 visible device | 0 | `:292-293` |
| ranks per node (`MPI_Comm_split_type`) == visible devices | `MyRankInNode()` | `:296-297` |
| ranks per processor name == visible devices | `MyRankInProcessor()` | `:298-299` |
| otherwise | `MyProc() % gpu_device_count` (**global** rank modulo) | `:300-301` |

After the table applies:
- With more than one visible device and a non-minimal init, AMReX applies site heuristics for `nersc.perlmutter` (reversed order) and `olcf.frontier` (fixed order for 8 GCDs). On any other machine it only prints a warning, and only when verbose: "Multiple GPUs are visible to each MPI rank … may lead to incorrect or suboptimal rank-to-GPU mapping" (`:305-327`).
- Then `cudaSetDevice(device_id)` (`:330-333`).

**Runtime knobs.**
- There is **no** `amrex.the_device_id` or other ParmParse parameter for device selection. `Device::Initialize` reads only `amrex.max_gpu_streams` and `device.verbose` (`:240-248`), and `grep` finds no `the_device_id` or `CUDA_VISIBLE_DEVICES` anywhere under `Src/`.
- The only override is the C++ argument `a_device_id`. The FDS driver could compute it from `OMPI_COMM_WORLD_LOCAL_RANK`/`SLURM_LOCALID` and pass it in.

**Ranks per node > devices per node.**
- No abort. The modulo rule assigns devices round-robin over the global rank (`:301`).
- With the same rank count on every node and ranks placed in blocks, this is balanced when ranks-per-node is a multiple of devices-per-node. Example: 4 ranks and 2 GPUs per node gives 0,1,0,1 on every node.
- It is unbalanced otherwise. Example: 3 ranks and 2 GPUs per node, over 2 nodes, gives 0,1,0 and 1,0,1. This case is derived from the rule, not run.

**Detecting shared devices** (`:342-378`).
- Every rank all-gathers its device UUID. From these AMReX sets `num_devices_used` (unique UUIDs) and `num_device_partners` (ranks on my device).
- In verbose mode, if devices < ranks, the IO rank prints the warning "There are more MPI processes than the number of unique GPU devices. This is not necessarily a problem." (`:428-442`).
- There is **no abort** when ranks share a device.
- `num_device_partners` is later used to divide the default arena size (§2). It is based on UUIDs, so it stays correct when a wrapper hides devices through `CUDA_VISIBLE_DEVICES`.

**AMReX guidance.** "Each MPI rank offloads its work to a single GPU. Multiple ranks can share the same device, but for best performance we usually recommend (MPI ranks == Number of GPUs)" (`GPU.rst:41-42`). Several ranks per GPU "can make sense … when you have some portion of your code that is not GPU accelerated" (`GPU.rst:2157-2168`). Round robin "would not be aware of locality benefits" (`GPU.rst:2170-2173`).

## 2. Several ranks per GPU: NVIDIA MPS and AMReX memory

**What MPS does** (NVIDIA MPS docs, Release 615, 2026-09-09: https://docs.nvidia.com/deploy/mps/index.html).
- A control daemon (`nvidia-cuda-mps-control`) and a per-user server (`nvidia-cuda-mps-server`) give the client processes a shared connection to the GPU.
- With MPS, kernels from different processes can run concurrently, which "remove[s] an unnecessary point of serialization" (https://docs.nvidia.com/deploy/mps/when-to-use-mps.html, "Identifying Candidate Applications").
- Without MPS, several processes can still share a GPU in the DEFAULT compute mode (same page, "GPU Compute Modes"). NVIDIA's pages do not say how their kernels are scheduled then. **Assumption:** the kernels are time-sliced rather than run concurrently.
- Clients have fully isolated GPU address spaces. A fatal GPU fault in one client is reported to every client on that GPU (same page, "Memory Protection and Error Containment").

**Start and stop, single user** (https://docs.nvidia.com/deploy/mps/common-tasks.html, "On a Single-User System"; https://docs.nvidia.com/deploy/mps/quick-start.html):
```
export CUDA_VISIBLE_DEVICES=0                     # daemon side only
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps    # default /tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log     # default /var/log/nvidia-mps
nvidia-cuda-mps-control -d                        # start (legacy v2; v3: add -p 3)
echo ps | nvidia-cuda-mps-control                 # list clients
echo quit | nvidia-cuda-mps-control               # stop (v3: nvidia-cuda-mps-control -q)
```
- Clients must set the same `CUDA_MPS_PIPE_DIRECTORY`/`CUDA_MPS_LOG_DIRECTORY`.
- **"CUDA_VISIBLE_DEVICES should not be set in the client's environment"** (common-tasks). A daemon started with `CUDA_VISIBLE_DEVICES` remaps device indices for its clients; NVIDIA recommends using UUIDs (https://docs.nvidia.com/deploy/mps/appendix-environment-variables.html, `CUDA_VISIBLE_DEVICES`).
- This rules out a per-rank `CUDA_VISIBLE_DEVICES` wrapper under MPS. Use AMReX's node-local rule or `a_device_id` instead.
- `CUDA_MPS_ACTIVE_THREAD_PERCENTAGE` caps the share of threads (SMs) a client may use. It does not reserve them. NVIDIA suggests 100%/n or 100%/(0.5 n) for n clients (when-to-use, "Dynamic Execution Resource Provisioning").

**Limits.**

| Limit | Value / behaviour | Source |
|---|---|---|
| OS | Linux and QNX only | when-to-use, "Limitations" |
| Clients | up to 60 client CUDA contexts per device with the default `CUDA_DEVICE_MAX_CONNECTIONS=2` (fewer if raised) | when-to-use, "Client-Server Connection Limits" |
| Users | one user per system may have an active MPS server; same UID as the server | when-to-use, "Limitations" |
| Device memory | **not partitioned unless configured.** The optional limit is set by `set_default_device_pinned_mem_limit` / `set_device_pinned_mem_limit` (control daemon) or `CUDA_MPS_PINNED_DEVICE_MEM_LIMIT` (client, e.g. `0=4G`); over the limit, allocations return out-of-memory | when-to-use, "MPS Device Memory Limit"; appendix-environment-variables |
| Per-client overhead | each client allocates its own context storage, which scales with the threads it is allowed | when-to-use, "Dynamic Execution Resource Provisioning" |
| Abnormal exit | killing a client without synchronizing "can leave the MPS server and other MPS clients in an undefined state" | when-to-use, "Limitations" |
| Debugging | cuda-gdb only without MPS | when-to-use, "CUDA-GDB" |

**Consumer and mobile GPUs.**
- The current MPS documentation (Release 615) names **no product-line restriction**. It neither lists GeForce or mobile GPUs as supported nor excludes them. The hardware conditions it states are Linux, 64-bit and same UID. Static SM partitioning needs Ampere or newer.
- **MPS on the GPU of the owner-provided test machine is therefore unverified** and is a test item (§7, L-3).
- A developer-forum answer reports no observed difference between consumer and data-centre GPUs, but says this "isn't published anywhere" (https://forums.developer.nvidia.com/t/cannot-use-stream-ordered-async-memory-allocator-with-cuda-mps/230229). That is anecdotal, not a spec.

**AMReX memory under sharing** (`Src/Base/AMReX_Arena.cpp`).
- Defaults: `the_arena_init_size = 8 MiB` (`:40`) and `the_arena_is_managed = false` (`:59`).
- On a non-minimal GPU init: **`the_arena_init_size = totalGlobalMem / numDevicePartners / 4 * 3`** (`:427`). Each rank takes 3/4 of the **total**, not the free, device memory, divided by the number of ranks on that GPU. The partner count comes from the UUID gather in §1.
- The pinned-arena release threshold is `totalGlobalMem / numDevicePartners / 2` (`:435`).
- The environment variable `AMREX_THE_ARENA_INIT_SIZE` overrides the computed default (`:439-443`), and the ParmParse `amrex.the_arena_init_size` overrides both (`:446`).
- The chunk is allocated at startup with `cudaMalloc` (`:476-481`, `:323-334`). If that fails, `out_of_memory_abort` stops the run (`:341-350`).
- On CPU builds `The_Arena()` is `The_BArena()` and the init size is unused (`:483`).
- Docs discrepancy: `GPU.rst:631-636` and `RuntimeParameters.rst:1044-1050` say "3/4 of the total device memory" and do not mention the division by partners that the source does.

**What to set when N ranks share one GPU.**
- AMReX already divides the default by N, so the ranks do **not** each grab 3/4 of the GPU.
- The default still assumes the ranks own the whole GPU. It ignores the CUDA context of each rank (larger per client under MPS, see above), a display on the test machine's GPU, and other processes. **Assumption:** these can make the startup `cudaMalloc` fail.
- Spec proposal (FR-0xx):
  - with N > 1 ranks per GPU, set `amrex.the_arena_init_size` explicitly to roughly (free memory − headroom)/N, and record the value in the run log;
  - keep `amrex.the_arena_is_managed=0`: managed memory is rejected for production and allowed only for debugging (ADR-001, "Kernel interface and device data for K2");
  - optionally set `amrex.abort_on_out_of_gpu_memory=1` (`Arena.cpp:60,258-262`) in tests. It only acts on managed-memory arenas (`:258-259`). A device arena already aborts on `cudaMalloc` failure.
- The headroom value is not specified by AMReX or NVIDIA; measure it (L-2).

## 3. CPU-only MPI + OpenMP fallback

- **Build.** Same sources, with `AMReX_GPU_BACKEND=NONE`, `AMReX_MPI=ON` and `AMReX_OMP=ON` (`Tools/CMake/AMReXOptions.cmake:124-125, 265, 276`). This matches both existing installs, `(local AMReX install)` and `(local AMReX install built against HYPRE 2.32)`: their `AMReX_Config.H` defines `AMREX_USE_MPI` and `AMREX_USE_OMP` and leaves `AMREX_USE_GPU`/`AMREX_USE_CUDA` undefined.
- **Tiling.**
  - `TilingIfNotGPU()` is `constexpr true` in CPU builds and `Gpu::notInLaunchRegion()` in GPU builds (`Src/Base/AMReX_MFIter.H:11-15`; `GPU.rst:1286-1288`).
  - Default tile size: on CPU in 3-D, tiling is on with size 8 in y and z; on GPU, tiling is off (`RuntimeParameters.rst:1247-1254`).
  - OpenMP threads share the tiles of an `MFIter` loop.
- **Threads.**
  - `amrex.omp_threads` is `system` (default: use `OMP_NUM_THREADS`), `nosmt` (physical cores), or an integer. An integer overrides `OMP_NUM_THREADS` (`Src/Base/AMReX_OpenMP.cpp:156-194`; `RuntimeParameters.rst:660-671`).
  - In verbose mode AMReX warns when threads × ranks-per-node exceeds the core count (`AMReX.cpp:541-556`).
  - Binding: `mpirun --map-by ppr:R:node --bind-to core` plus `OMP_NUM_THREADS=cores/R`. Pin the thread count and use `OMP_DYNAMIC=false` for reproducibility (FR-005 (iv)). Open MPI 5 options: https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/mpirun.1.html.
  - The HPC SDK's HPC-X 2.50 (Open MPI 5, PRRTE) maps by core by default (https://docs.nvidia.com/hpc-sdk/release-notes/index.html, "Communication libraries"). Set the mapping explicitly.
- **Kernel rule (IR-007, ADR-001 v0.3 wording).**
  - No host-threading OpenMP (`parallel do`) inside device kernels, and `target` directives only in kernel files. AMReX does host threading and tiling from the driver.
  - This replaces the earlier "no `!$omp` inside device kernels" wording.
  - On CPU, the threading level is the `MFIter` tile. Kernels must accept `lo/hi` tile bounds (`driver-options.md` §2). The P1 shim runs whole boxes with tiling off (`driver-options.md` §4).
- **GPU and OpenMP are *not* mutually exclusive in this AMReX CMake. This corrects the brief.**
  - `AMReX_OMP` is a plain option (`AMReXOptions.cmake:276`), and `AMReXParallelBackends.cmake:48-58` links OpenMP without checking the GPU backend.
  - AMReX recommends host OpenMP in GPU builds only for host work, via `#pragma omp parallel if (Gpu::notInLaunchRegion())` (`GPU.rst:1631-1648`).
  - Proposal (revised for A-43): `AMReX_OMP` in the GPU build is a **measured choice** (§5.2, §5.4). Start with OFF, and switch to ON only if threaded host sections pay off. **Open:** how nvfortran `-mp=gpu` (K2 offload) interacts with host OpenMP in the same executable (see §8).

## 4. GPU-aware MPI

**Default and detection** (`Src/Base/AMReX_ParallelDescriptor.cpp`).
- The static default is `false` (`:67-71`).
- In `ParallelDescriptor::Initialize` (`:1550-1590`), a CUDA build then **auto-detects** support:
  - Open MPI: if `mpi-ext.h` defines `OMPI_HAVE_MPI_EXT_CUDA`, AMReX uses `MPIX_Query_cuda_support()` (`:14-15, 1554-1555`);
  - MPICH: `MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA)` (`:1556-1560`).
- `amrex.use_gpu_aware_mpi` overrides the detected value (`:1586-1587`).
- `GPU.rst:2189-2193` agrees ("MPI-dependent"). `RuntimeParameters.rst:682-693` still says `false` and "does not enable GPU-aware MPI by itself".
- AMReX does not check a forced `1` against the library's capability. **Assumption:** forcing it on a non-CUDA-aware MPI would hand device pointers to MPI and crash or corrupt data. Spec proposal: never force it; let detection decide, and log the value.

**What changes.**
- **On:** communication buffers come from device memory. `The_Comms_Arena()` is the device arena or a separate device `CArena` (`Arena.cpp:530-543`). `FillBoundary`/`ParallelCopy` pack and send device buffers directly (`AMReX_FabArrayCommI.H:884, 982`).
- **Off:** `The_Comms_Arena()` is `The_Pinned_Arena()` (`Arena.cpp:544-545`), which is host memory from `cudaHostAlloc(..., cudaHostAllocMapped)` (`:221`). Device kernels pack into mapped pinned host buffers, and MPI sees host pointers. Collectives on device data stage through pinned host memory:
  - `ParallelAllReduce` on `Gpu::DeviceVector` (`AMReX_GpuParallelReduce.H:22-24, 41-47`);
  - `FabArray` broadcast (`AMReX_FabArrayCommI.H:792-796`).
- Under IR-007, FDS domain-wide reductions run on the host between kernel passes, so they do not depend on this flag.

**What the MPI library needs.**
- Open MPI has to be built with CUDA support, preferably through a CUDA-enabled UCX. Check it with (https://docs.open-mpi.org/en/v5.0.x/tuning-apps/networking/cuda.html, §11.2.6.1-11.2.6.2):
  ```
  ompi_info --parsable --all | grep mpi_built_with_cuda_support:value   # expect ...:value:true
  ompi_info | grep "MPI extensions"                                      # expect "cuda" in the list
  ucx_info -v                                                            # configure line must show --with-cuda
  ```
- Intra-node transfers use CUDA IPC:
  - the smcuda BTL has `btl_smcuda_use_cuda_ipc` and `btl_smcuda_use_cuda_ipc_same_gpu` (on by default), plus `btl_smcuda_cuda_ipc_verbose 100` for diagnostics;
  - with UCX, the documented example is `mpirun --mca pml ucx -x UCX_TLS=rc,sm,cuda_copy,gdr_copy,cuda_ipc` (same page, §11.2.6.8 and §11.2.6.13).
- Non-blocking reduction collectives and one-sided operations do not accept CUDA buffers with CUDA-aware UCX (§11.2.6.15). AMReX's halo exchange is point-to-point; **not verified:** whether AMReX calls any of those operations on device buffers.
- Slurm sites: if `cgroup.conf` has `ConstrainDevices=yes`, IPC can break. Use `--gpu-bind=none` rather than `closest` (`Faq.rst:180-182`).

**The development machine's Open MPI 5.0.7 (apt) is not CUDA-aware at compile time.** Its header `/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/mpiext/mpiext_cuda_c.h` has `#define MPIX_CUDA_AWARE_SUPPORT 0`. The CUDA extension is present (`mpi-ext.h`: `OMPI_HAVE_MPI_EXT_CUDA 1`), so AMReX will call `MPIX_Query_cuda_support()`. **Assumption to verify (A-0xx):** it returns 0 and AMReX picks `use_gpu_aware_mpi=false`. Confirm with `ompi_info` once a GPU build exists. Nothing was run.

**NVIDIA HPC SDK 26.9** (install planned for `(local NVIDIA HPC SDK directory)`, 2026-09-26).
- Bundles CUDA 12.9U1 and 13.3U1 and HPC-X 2.50 (Open MPI 5, UCX). HPC-X 2.20 is the fallback for CUDA 12.0-12.3 drivers.
- The MPI wrappers pick the MPI matching the installed driver (https://docs.nvidia.com/hpc-sdk/release-notes/index.html, §2 and "Communication libraries").
- HPC-X documents GPU transports and GPUDirect RDMA through UCX (https://docs.nvidia.com/networking/display/hpcxv226/HPC-X-General-Support; "Known Issues" lists GPU-buffer workarounds such as `UCX_TLS=^gdr_copy`).
- The HPC SDK release notes do not literally say "CUDA-aware". **Treat HPC-X as CUDA-aware only after the `ompi_info` check above** (A-0xx).
- Release-note caveats that matter here:
  - `MPI_Send`/`MPI_Recv` throughput "can have significantly degraded throughput" on `cudaMallocAsync` memory. AMReX's main arenas use plain `cudaMalloc`/`cudaHostAlloc` (`Arena.cpp:221, 323`).
  - To cut HPC-X init time when GPU communication is not used: `OMPI_MCA_coll_ucc_enable=0` and `UCX_MODULES=^cuda`.

## 5. Several CPU cores per GPU (A-43)

**Why.** The owner (A-43) finds one CPU core per GPU too slow for FDS because of host-side work:
- Smokeview and VTK output (host-side per D-027);
- setup;
- anything not yet ported.

In FDS today that work is serial within a rank:
- `dump.f90` (13,276 lines), `vtkf.f90`, `read.f90` and `init.f90` contain **no** `!$OMP` directives. All 356 `!$OMP` lines are in the compute files: `divg`, `mass`, `velo`, `pres`, `wall`, `radi`, `turb`, `func`, `pois`, `fire` (`grep` count, `Source`).
- Output runs mesh by mesh: `main.f90:623-627` calls `DUMP_MESH_OUTPUTS` per mesh. That routine starts with `CALL POINT_TO_MESH(NM)` (`dump.f90:89`), which sets process-global module pointers (`mesh.f90:505-916`; R-26).

The two routes compared below can be combined: async output works in both.

### 5.1 What AMReX offers for asynchronous output

| Item | Behaviour (AMReX source) |
|---|---|
| Switches | `amrex.async_out` (default `false`) and `amrex.async_out_nfiles` (default 64, clamped to [1, NProcs]), read in `AsyncOut::Initialize` (`Src/Base/AMReX_AsyncOut.cpp:13-14, 27-32`). It is called from `amrex::Initialize` except under `Init_minimal` (`AMReX.cpp:716-718`). |
| Thread | One `BackgroundThread` per rank: a single `std::thread` with a FIFO job queue (`AMReX_AsyncOut.cpp:51-53`; `AMReX_BackgroundThread.cpp:5-8, 19-50`). It is **not** OpenMP, so it works with `AMReX_OMP=OFF`. It is joined at finalize (`AMReX_BackgroundThread.cpp:10-17`). |
| API | `AsyncOut::Submit(std::function<void()>)` runs any callable on that thread. `Finish()` blocks until the queue is drained. `Wait()`/`Notify()` order the ranks that share a file (`AMReX_AsyncOut.H:50-83`; `.cpp:95-140`). There is no function named `WriteAsyncPlotfile`: async is a mode of the normal writers. |
| MPI thread level | If `async_out_nfiles < NProcs`, AMReX requires `MPI_THREAD_MULTIPLE` at runtime and aborts otherwise (`AMReX_AsyncOut.cpp:35-44`). The reason: `Wait`/`Notify` post `Abarrier`s on a split communicator from the background thread (`:47, 112-139`). With `nfiles ≥ NProcs` each rank has its own file (`ispot = 0`) and the thread makes no MPI calls. AMReX requests `MPI_THREAD_MULTIPLE` only when built with `AMReX_MPI_THREAD_MULTIPLE=ON` (default OFF, `AMReXOptions.cmake:268-271`; `AMReX_ParallelDescriptor.cpp:307-315`), or when the driver initialises MPI itself. See also `IO.rst:104-120`. |
| Writers covered | `WriteMultiLevelPlotfile`/`WriteSingleLevelPlotfile`: header job submitted (`AMReX_PlotFileUtil.cpp:198-216`), level data through `VisMF::AsyncWrite` (`:219-225`; `AMReX_VisMF.cpp:2334-2364`). Particles (`AMReX_ParticleIO.H:433-437`). `Amr`/`AmrLevel`/`StateData`/`FabSet` (`IO.rst:124-137`). Native FAB format only; `fab.format` is ignored (`IO.rst:122-123`). |
| What is copied | `VisMF::AsyncWriteDoit` (`AMReX_VisMF.cpp:2367-2598`), on the **calling** thread: <br>1. per-FAB min/max, on the device if the data is there (`:2397-2398, 2425-2429`); <br>2. `MPI_Gatherv` of header data to the IO rank (`:2462-2464`); <br>3. a snapshot of every local FAB: device/managed data goes into a new FAB in `The_Pinned_Arena()` by device copy or `dtoh_memcpy_async` (`:2471-2481`); host data is copied into the CPU arena, or moved for rvalues (`:2482-2490`). <br>The `MFIter` destructor synchronises the GPU streams (`AMReX_MFIter.cpp:246-252`), so the thread receives only complete host memory at `Submit` (`:2494`). The thread then writes the header and the FABs (`:2574-2596`). |
| Cost | The main thread still pays the device-to-host copy and the gather. Every pending write holds a full host (pinned) copy of the output fields. `IO.rst:140-145` warns that the extra thread can oversubscribe cores when OpenMP uses all of them. |
| Not offered | No Smokeview or FDS-VTK writer: the only VTK-related code in `Src/` is EB→PVD (`Src/EB/AMReX_EBToPVD.*`) and SENSEI adaptors. `grep` found no facility for dedicated I/O ranks. |

**Does it help FDS's Smokeview/VTK writers?**
- **AMReX-native plotfiles and checkpoints:** yes, with `amrex.async_out=1` and no code change.
- **FDS's Fortran writers** (`DUMP_SLCF`, `DUMP_SMOKE3D`, `DUMP_BNDF`, the VTKHDF variants, per-mesh files): only the mechanism carries over. `AsyncOut::Submit` could run them, but three things are needed first:
  1. a host-side snapshot of the fields each writer reads, taken the way `VisMF` does it;
  2. writers refactored to take that snapshot as arguments, without `POINT_TO_MESH` or module state, so they cannot race the main thread;
  3. confirmation that the Fortran runtime's I/O is safe on a non-main thread. **Unverified** for gfortran and nvfortran.
- The alternative is dedicated I/O ranks that receive snapshots over MPI and run the writers in their own process. AMReX has no support for this; it would be project work.

### 5.2 Host OpenMP in a GPU build
- `AMReX_OMP=ON` together with CUDA is allowed (`AMReXOptions.cmake:276`; `AMReXParallelBackends.cmake:48-58`).
- In a GPU build `Gpu::inLaunchRegion()` is `true` by default (`AMReX_GpuControl.cpp:7`; `AMReX_GpuControl.H:86-89`). As a result:
  - `TilingIfNotGPU()` returns `false` (`AMReX_MFIter.H:11-12`);
  - AMReX's own loops are written `#pragma omp parallel if (Gpu::notInLaunchRegion())`, e.g. `AMReX_FabArray.H:3040-3043`, `AMReX_Geometry.cpp:244`, `AMReX_DistributionMapping.cpp:1718`; 167 occurrences in `Src/Base` + `Src/LinearSolvers`. On the host these loops run **single-threaded**.
- So `AMReX_OMP=ON` speeds up only host code that we thread explicitly: `!$OMP` in Fortran host routines, or `MFIter` loops over host (pinned) data with tiling requested and an unconditional `omp parallel`.
- `Gpu::setLaunchRegion(false)` sets a process-wide variable (`AMReX_GpuControl.cpp:7`; `AMReX_GpuControl.H:110-114`) and does not affect `ParallelFor` (`:94-97`). It is not a per-thread tool.
- **Reconciliation with §3:** `AMReX_OMP` in the GPU build is now a **measured choice**:
  - OFF for route (a) and for (b) with async output only;
  - ON only if §5.4 shows that threaded host sections pay off (FR-0xx).

  The nvfortran `-mp=gpu` interaction (§8) has to be settled before ON is chosen.

### 5.3 Route comparison

| | (a) N ranks per GPU under MPS | (b) 1 rank per GPU + host threads and/or async I/O |
|---|---|---|
| Host work (output, setup, unported) | N-way parallel with **no code change**: each rank handles only its own boxes, as FDS does per mesh today (`main.f90:623-627`) | async: AMReX-native output overlaps compute for free (§5.1). FDS writers need snapshot + refactor. OpenMP: output/setup have no `!$OMP` today, so threading must be written |
| Thread-safety risk | none added (separate processes) | **risk to verify:** unported FDS host code relies on process-global module pointers (`POINT_TO_MESH`, `dump.f90:89`; R-26), so threading its mesh loop or running a writer concurrently with host code is unsafe until refactored. Fortran runtime I/O thread safety is unverified |
| MPI / ghost exchange | cross-rank ghost exchanges become MPI messages, staged through pinned host memory unless GPU-aware (§4). With 1 rank they are local device copies (`FB_local_copy_gpu`, `AMReX_FBI.H:546`) | no extra MPI |
| Box count | needs ≥ N boxes per GPU at the 16-cell box floor (guidance, `requirements.md:356`). Example: D-024 level 0 (32×32×80) gives at most 20 boxes of 16³ (derived). Smaller boxes per rank give fewer threads per kernel launch. **Assumption:** that means lower GPU utilisation | largest boxes per GPU |
| GPU memory | one CUDA context (and MPS client storage) per rank; arena split N ways (`AMReX_Arena.cpp:427`); limit of 60 MPS clients per GPU | one context |
| Dependencies | MPS, **unverified** on consumer GPUs (§2). **Assumption:** without MPS the kernels are time-sliced | `MPI_THREAD_MULTIPLE` only if `async_out_nfiles < ranks` |
| Reproducibility | rank-count dependence already covered by FR-005 (iii) | thread count fixed per FR-005 (iv) |

**Default hypothesis (to confirm or reject in §5.4).**
- **Transition period, while host code is unported and FDS writers are unchanged:** route (a) with a small N (2-4) under MPS, plus `amrex.async_out=1` for AMReX-native output. (a) is the only route that parallelises unmodified serial Fortran host code without a thread-safety refactor.
- **Later, as kernels move to the device and output moves to snapshot-based writers:** route (b) with async output becomes the target (one context, no MPS, no extra MPI).
- **If MPS does not work on the test machine:** measure (a) without MPS. If it still helps on output steps, keep it; otherwise use (b) with async output.

### 5.4 Measurement plan (A-43)
- **Owner:** AMReX Integration Lead. The method is reviewed by the V&V Lead.
- **When and where:** on the owner-provided test machine, once the test machine and a GPU build are ready. Numbers land in spec v0.4.12.
- **Case:** `<anchor case TBD>`. Same inputs, same output interval (`<Δt_out>`, at least `<K>` output steps) and the same `max_grid_size` in every configuration. Route (a) additionally needs at least N boxes.
- **Configurations**, for N ∈ {1, 2, 4, …, `<P>` = test-machine physical cores}:
  - **C1:** 1 rank, 1 thread, synchronous output (baseline).
  - **C2a:** N ranks under MPS. **C2b:** N ranks without MPS (control).
  - **C3a:** 1 rank, N threads (`AMReX_OMP=ON` build, threaded host sections).
  - **C3b:** 1 rank + `amrex.async_out=1` (AMReX-native plotfiles; the FDS-writer path only if a snapshot prototype exists).
  - **C3c:** C3a + C3b.
- **Timing:**
  - TinyProfiler with `tiny_profiler.device_synchronize_around_region=1` (`AMReX_TinyProfiler.cpp:353`) so that device time is attributed to regions;
  - `BL_PROFILE_REGION` regions `<Step::compute>` and `<Step::output>` (`AMReX_BLProfiler.H:422`);
  - wall time per step, reported separately for compute steps and output steps, plus the average over the run;
  - for async runs, also the time spent in `AsyncOut::Finish()` at the end. **Assumption:** background-thread write time does not appear in TinyProfiler regions;
  - peak device memory per rank, host RSS and pinned memory;
  - `<R>` repeats, median reported.
- **Correctness, against C1:**
  - explicit stages bitwise (FR-005 (i));
  - pressure within eps_H (FR-005 (iii));
  - async plotfile data bitwise identical to the synchronous write;
  - Smokeview files byte-identical where the writer is unchanged;
  - no aborts, and no MPS client hangs at exit.
- **Decision:**
  - The default route is the one with the lowest output-step time and average step time, among the configurations that pass every correctness check.
  - A route counts only if it improves the output-step time by at least `<X %>` over C1 and does not worsen the average step time by more than `<Y %>`.
  - If (a) and (b) are within `<Z %>`, prefer (b): no MPS dependency and no extra MPI.

## 6. Recommended configuration (spec proposal)

| Item | Default | Notes |
|---|---|---|
| Ranks per GPU | **at least 1; several CPU cores per GPU required (A-43).** Default route is chosen by §5.4. Hypothesis: 2-4 ranks per GPU under MPS during the transition (§5.3) | AMReX recommends 1 rank per GPU for GPU-bound work (`GPU.rst:41-42`). The owner's host-side I/O and unported code make one core per GPU too slow |
| Rank→GPU binding, multi-GPU node, no MPS | launcher places R = GPUs-per-node ranks per node in blocks (`mpirun --map-by ppr:R:node --bind-to core`, or Slurm `--ntasks-per-node=R --gpus-per-task=1`); AMReX then applies `MyRankInNode()` (`GpuDevice.cpp:296-297`) | Alternative: a wrapper sets `CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK` (or `$SLURM_LOCALID`), so each rank sees one device and AMReX picks 0 (`:292-293`). Variables are documented at https://docs.open-mpi.org/en/v5.0.x/tuning-apps/environment-var.html and https://slurm.schedmd.com/srun.html. **Not under MPS** (§2). For NUMA-correct pairing, prefer the wrapper, a Slurm GPU binding, or `a_device_id` from the driver, because AMReX's rule ignores topology (`GPU.rst:2170-2173`). |
| GPU-aware MPI | auto-detected; on only when `ompi_info` shows CUDA support | log `ParallelDescriptor::UseGpuAwareMpi()` at startup (FR-0xx) |
| Several ranks per GPU (route a) | with MPS running where available, `amrex.the_arena_init_size` set explicitly (§2), N ≤ 60 per device and N ≤ boxes per GPU at the 16-cell floor | Parallelises unmodified host Fortran (output, setup, unported code). Without MPS it runs only as the §5.4 control or as a measured fallback. Gain unmeasured (§5.4; `driver-options.md` §5, R-31) |
| One rank per GPU + host threads/async (route b) | `amrex.async_out=1`; `amrex.async_out_nfiles ≥` ranks, unless MPI is initialised with `MPI_THREAD_MULTIPLE` (`AMReX_AsyncOut.cpp:35-44`); host OpenMP only if §5.4 selects it | Async output helps AMReX-native plotfiles/checkpoints now. FDS Smokeview/VTK writers need a snapshot and a refactor first (§5.1) |
| Async output (both routes) | `amrex.async_out=1` for AMReX-native output | leave at least one core per rank for the writer thread (`IO.rst:140-145`) |
| CPU-only | `GPU_BACKEND=NONE`, MPI+OMP, `--bind-to core`, `OMP_NUM_THREADS=cores/ranks`, `OMP_DYNAMIC=false` | same inputs, same drivers |
| GPU build | `AMReX_GPU_BACKEND=CUDA`, explicit `CMAKE_CUDA_ARCHITECTURES`, HPC SDK 26.9 (CUDA 12.9, fallback 13.3), `AMReX_OMP` OFF by default, ON only if §5.4 selects threaded host sections | §3, §5.2 |
| Verbosity for runs under test | `amrex.v=1` (prints "CUDA initialized with N devices" and the sharing warning, `GpuDevice.cpp:428-442`) | |

## 7. What needs testing

Tolerance classes are those of `requirements.md` §2.2 (T0 bitwise, T2 verification, T3 physical) and §2.1a (eps_H for pressure). GPU vs CPU agreement cannot be expected at T0:
- AMReX device sum reductions combine partial sums with atomics (`AMReX_GpuReduce.H:126, 202`);
- **Assumption:** compilers may contract to FMA differently on device and host.

The class for kernel-level GPU vs CPU parity is an open question (§8).

### 7.1 Testable on one GPU (the owner-provided test machine)

| ID | Test | Checks | Pass criterion | Hardware |
|---|---|---|---|---|
| L-1 | 1 rank, 1 GPU, anchor cases | GPU build runs end to end; device selection | runs to completion; log shows "CUDA initialized with 1 device"; results vs CPU build per L-6 | test-machine GPU, HPC SDK |
| L-2 | 2 and 4 ranks sharing the GPU, **without** MPS | sharing works; arena sizing; no aborts | no abort; verbose warning "more MPI processes than … unique GPU devices" appears (expected); per-rank arena = `totalGlobalMem/N·3/4` unless set; explicit stages bitwise vs L-1 at the same BoxArray (FR-005 (i)); pressure within eps_H (FR-005 (iii)); record time and peak memory | same |
| L-3 | MPS availability on this GPU | `nvidia-cuda-mps-control -d` starts; server appears in `nvidia-smi` as M+C | daemon and server start, `echo ps` lists clients; failure recorded in `server.log`/`control.log` | same |
| L-4 | L-2 **with** MPS (if L-3 passes) | correctness and concurrency under MPS; clients do **not** set `CUDA_VISIBLE_DEVICES` | results identical to L-2 at the same rank count and layout, except device-atomic reductions (then eps_H/T2); no hang on normal exit; timing vs L-2 recorded, no pass threshold | same |
| L-5 | GPU-aware MPI on vs off, 2 ranks on the one GPU (CUDA IPC, same GPU) | detection and both comm paths | with HPC-X: `ompi_info` shows CUDA support; `UseGpuAwareMpi()` reports detected value; forced `amrex.use_gpu_aware_mpi=0` and detected `1` give identical fields (buffers change, arithmetic does not; **assumption**); with apt Open MPI detection gives 0 | same, HPC-X |
| L-6 | CPU-only MPI+OMP build vs GPU build, same inputs and layout | the same code works on both backends | whole runs at T2 (T3 where the case uses it); pressure within eps_H; kernel-level class per §7 | test-machine CPU + GPU (CPU half also on the development machine) |
| L-7 | Device selection with `CUDA_VISIBLE_DEVICES` | rule rows of §1 | `CUDA_VISIBLE_DEVICES=""` → abort "No GPU device found"; `=0` → device 0; wrapper per rank → device 0 each, `num_device_partners` still N (arena default divides by N) | same |
| L-8 | Memory headroom | arena init vs free memory with display on | find largest `amrex.the_arena_init_size` per rank for N=1,2,4 that starts reliably; record it as the test machine default | same |
| L-10 | Cores-per-GPU matrix (A-43) | C1, C2a/C2b, C3a/C3b/C3c over N = 1, 2, 4 … `<P>` (§5.4) | correctness per §5.4 against C1; timing table complete; default route picked by the §5.4 decision rule | test-machine GPU + all physical cores |
| L-11 | Async output | `amrex.async_out=1` with `async_out_nfiles ≥` ranks, and `<` ranks with and without `MPI_THREAD_MULTIPLE` | async plotfile bitwise identical to the synchronous one; `nfiles < ranks` without `THREAD_MULTIPLE` aborts with the `AsyncOut` message (`AMReX_AsyncOut.cpp:39-43`); snapshot pinned-memory peak recorded | same |
| L-12 | `AMReX_OMP=ON` + CUDA build | the build links; threaded host sections run | results bitwise identical to `AMReX_OMP=OFF` at fixed thread count (FR-005 (i), (iv)); no oversubscription warning with the async thread counted | same (build also possible on the development machine, compile only) |
| L-9 | Abnormal termination under MPS | one rank aborts (`amrex::Abort`) | the other ranks and the MPS server recover or are restarted cleanly; procedure documented | same |

### 7.2 Needs multi-GPU or multi-node hardware

| ID | Test | Checks | Pass criterion | Hardware |
|---|---|---|---|---|
| M-1 | Rank→GPU round robin, R ranks = G GPUs per node | `MyRankInNode()` row | each rank on a distinct GPU (UUIDs logged per rank), `num_devices_used = G·nodes` | ≥1 node with ≥2 GPUs |
| M-2 | R = 2G and R = 3G per node, and a non-multiple (e.g. 3 ranks, 2 GPUs) | modulo row and imbalance | per-GPU rank counts match the rule in §1; imbalance documented | ≥2 GPUs/node, 2 nodes for the non-multiple case |
| M-3 | Inter-node GPU-aware MPI (GPUDirect RDMA via UCX) | detection and fields on vs off | identical fields on vs off (assumption as in L-5); no UCX errors; transport confirmed with UCX diagnostics | ≥2 nodes with IB/RoCE and GPUDirect |
| M-4 | Rank-count invariance across nodes | halo and reductions over the network | explicit stages bitwise vs 1 GPU at the same BoxArray; pressure within eps_H; whole runs T2 | ≥2 nodes |
| M-5 | Strong and weak scaling | performance (NFR-0xx) | recorded efficiencies; thresholds set by the Spec Lead | ≥2 nodes × ≥2 GPUs |
| M-6 | NUMA/GPU affinity | wrapper or `a_device_id` vs AMReX default | each rank's cores are on the GPU's NUMA node (`nvidia-smi topo -m`, `--report-bindings`); timing difference recorded | multi-socket, multi-GPU node |
| M-7 | Load balancing across GPUs with AMR | `DistributionMapping` over GPUs as levels change | per-GPU work and memory within a stated imbalance after regrid | ≥4 GPUs |
| M-9 | Chosen cores-per-GPU route on a multi-GPU node | the route selected in L-10 with G GPUs and N ranks or threads per GPU; async output to a shared filesystem | L-10 correctness criteria; output-step time per GPU no worse than on the test machine at the same N (`<tolerance>`) | ≥1 node with ≥2 GPUs |
| M-8 | Several ranks per GPU with MPS on a data-centre GPU | MPS at scale, 60-client limit | as L-4, plus the per-node MPS start procedure under the scheduler | cluster with MPS allowed |

## 8. Open questions and risks

1. **MPS on a consumer/mobile GPU.** NVIDIA's MPS docs do not say whether it is supported (§2). L-3 decides. If MPS is unavailable, several ranks per GPU on the test machine run without MPS. **Assumption:** the kernels are then time-sliced (§2), so expect little or no speed-up from sharing.
2. **Test-machine GPU memory is unknown** (model not yet named), and the display takes part of it. The 3/4-of-total default may not fit with N ranks plus contexts. L-8 sets the test machine values. Large cases may not fit on the test machine's GPU at all.
3. **HYPRE is CPU-only** in the project's builds. On the GPU path, MLMG's own bottom solver is the default (ruling). A HYPRE bottom solver in a GPU build would need a CUDA HYPRE matching `AMReX_GPU_BACKEND` (`driver-options.md` §6, R-30). Not planned.
4. **The Fortran shim prototype (P1) is CPU-only.** It forecloses tiling and GPU (`driver-options.md` §4). Multi-GPU tests need ported K2 kernels.
5. **HPC-X CUDA-awareness** is assumed from HPC-X docs, not from an HPC SDK statement. Confirm with `ompi_info` after the 2026-09-26 install. The HPC-X 2.50 default needs a CUDA ≥ 12.4 driver (`(local NVIDIA HPC SDK directory)/README.md`) on the test machine.
6. **nvfortran OpenMP offload + host OpenMP.** How `-mp=gpu` Fortran kernels coexist with an AMReX build that has `AMReX_OMP` ON or OFF has not been checked (R-39, A-31).
7. **Kernel-level GPU vs CPU parity class.** T0 (D-022) is defined against CPU FDS. Whether GPU kernels can meet T0 against the CPU build, or only T1/T2, needs a Spec Lead ruling once L-6 has data.
8. **Global-rank modulo mapping** (`GpuDevice.cpp:301`) can be unbalanced when ranks-per-node is not a multiple of GPUs-per-node. The driver should warn about or reject such layouts (FR-0xx), or set `a_device_id` itself.
9. **No cluster is named.** Every item in §7.2 is untestable until one is available (charter Q6, R-12).
10. **FDS host code thread safety (A-43, route b).** FDS output and setup rely on process-global module pointers (`POINT_TO_MESH`, `dump.f90:89`; `mesh.f90:505-916`). Threading them, or running a writer beside other host code, is a **risk to verify**, not a known fact. The runtime I/O thread safety of gfortran and nvfortran on a non-main thread is also unverified.
11. **FDS writers under async output** need a host snapshot and explicit-argument refactor of `DUMP_SLCF`/`DUMP_SMOKE3D`/`DUMP_BNDF`/VTKHDF (§5.1), or dedicated I/O ranks. AMReX offers neither. The size of that work is not estimated.
12. **Host memory on the test machine.** Every pending async write holds a pinned host copy of the output fields (`AMReX_VisMF.cpp:2471-2481`). Route (a) multiplies per-rank host state by N. Test-machine RAM and physical core count (`<P>`) are unknown.
