# Proposed ADR-001 entry: data-movement kernels in the C++ driver layer (for the Architect to paste; ADR-001 itself is not edited here)

Status: accepted as D-066 and placed in ADR-001 as the subsection "Data-movement kernels in the C++ driver layer" (with the streams-per-box constraint of D-068). The text below is the original proposal, kept for the record. Original status: proposal from the stage-1 spike (plan `stage1-gpu-spike-plan.md`, sections 7.12 and 4, WP9, WP10, WP12). It needs the Architect's ruling and the owner's sign-off; nothing below is in force until then.

## Proposed text (place under "Kernel style", after the K2 default / K1 per-kernel fallback ruling)

**Ruling (data-movement kernels).** A *data-movement kernel* is a device loop in the C++ driver layer that only copies, gathers, scatters, packs, unpacks, fills or checksums values, with no physics: no FDS formula, no property lookup, no branch on a physical state, and no arithmetic other than the index arithmetic and the bit operations of a checksum. Data-movement kernels in the C++ driver layer (examples: the flux read-out pack and the interface-face gather and scatter of the flux hooks, the wall-state upload staging and the device checksum of the wall seam, the SOLID mask copy, ghost and periodic-partner copies) **may be written as restricted C++ `ParallelFor` (K1)** and need no per-kernel K1 exception. **All physics kernels stay K2** (generated Fortran with OpenMP `target`), unchanged by this entry; a loop that contains any physics term, however small, is a physics kernel.

**Conditions.**
1. The kernel is bitwise testable against a host reference: values are moved or bit-combined, never re-computed. (A checksum is an integer XOR/rotate fold; it is compared bit for bit with the host fold.)
2. The GPU build flags rule (K2 `nofma`, nvcc `--fmad=false`, fast math off) applies to these translation units as well, so that a neighbouring physics expression is never fused across the hand-off.
3. Launch and ordering follow the launch model below.

**Reason.** The driver layer is already C++ and AMReX owns its streams, arenas and `ParallelFor`; a K2 kernel for a pure copy would add an nvfortran translation unit and a blocking launch (see below) for no maintenance gain, and the single-source argument for K2 (one FDS loop, one kernel) does not apply because there is no FDS loop behind a copy.

## Facts the ruling relies on (run, cc 8.9 test-machine GPU, nvfortran 26.9, nvcc, AMReX CUDA build; plan section 7.12)
- A K2 target launch (`target teams distribute parallel do`, no `nowait`) is **blocking**: a kernel of about 285 ms made the call return after about 285 ms, and the following device synchronisation took about 3.5 us. A K1 launch returns in about 1.6 to 3 us.
- A K2 launch does **not** wait for work on other streams: launched right after a kernel on another stream with no sync, it read stale data in 10 of 10 runs. The host must synchronise that stream before a K2 launch that reads what an AMReX stream wrote.
- A `cudaMemcpy` placed after a K2 launch is safe (0 of 10 stale; the launch already completed); after a kernel on a non-blocking stream it was stale in 10 of 10 runs.
- `nowait` K2 launches cost about 2 ms each to issue and are unordered against blocking K2 launches (10 of 10 runs saw an unfinished writer without a task wait). They are not usable.
- Blocking K2 launches issued from different host threads overlap on the device (20 small launches: 25.3 ms on 1 thread, 26.5 ms on 2 threads, 41 ms on 3 threads). "One stream per box" for K2 therefore means one host thread per box; K2 cannot be placed on an AMReX stream.

## Consequences for the driver layer
- Data-movement kernels written as K1 can be enqueued on the development machine stream and overlap with other boxes; K2 kernels are host-synchronous points. The ordering rule for a mixed sequence is: stream sync before each K2 launch that consumes a stream result; no extra sync after a K2 launch; copy-back (`cudaMemcpy`) after a K2 launch needs none.
- The checksum copy-back of the wall-state test (8 bytes) is the only per-stage sync it adds, and only in the first-run/CI check mode.
