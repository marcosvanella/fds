# FDS-AMReX design documents: overview

This directory is a snapshot of the design and specification documents of the FDS-AMReX effort, written by the project's planning team (roles in `charter.md` §5). This page was added for the snapshot, and `MANIFEST.md` lists the files that were left out.

## What the effort is

FDS today runs on a fixed set of user-defined rectilinear meshes whose resolution cannot change during a run. The goal is to give FDS block-structured adaptive mesh refinement through the AMReX library: a hierarchy of properly nested, uniform Cartesian levels that are refined and coarsened during the run according to solution criteria, while FDS physics and its verification record stay intact (`charter.md` §1-2). The code base is FireX, the FDS development line that can, among other things, offload the HYPRE pressure solve to GPUs. Existing FDS input files are meant to run unchanged, with AMR enabled by a new namelist. CPU runs use MPI and OpenMP. With an NVIDIA GPU build the full time step is to run on the device, with I/O as the only exception, and the kernels must stay readable for FDS Fortran developers (`charter.md` §3, objectives O6-O8; D-027).

## Current phase

The project is in its specification phase (roadmap Phase 0, discovery, and Phase 1, data-structure and interface specification). No AMR source code has been added to this repository; the documents specify that code lands on the `FDS-AMReX` branch only after the specification phase (D-005, D-034, D-037). D-044 allows documentation snapshots such as this one during the specification phase. The small prototypes used to answer specific questions were built outside the repository. Milestone M1 ("Architecture decided") requires the three ADRs to be accepted; in this snapshot they are still Proposed, with some sub-decisions accepted. Unless a document says otherwise, `file:line` citations refer to FireX commit `36975d765f`.

## How the documents are organized

Read these first, in this order:

1. `README.md`: document index, decision log, action log and changelog.
2. `charter.md`: problem, goal, objectives, scope, constraints and open questions for the project owner.
3. `roadmap.md`: principles, phases and milestones with their gates.
4. `requirements.md`: numbered requirements, each with a verification method and acceptance criterion.
5. `adr/README.md`, then ADR-001 (driver architecture), ADR-002 (time stepping across levels) and ADR-003 (solid geometry on refined levels).

`risks.md` is the risk register, and `spec-responses.md` answers questions raised in the ADRs. The subdirectories each cover one topic: `amrex/` (AMReX integration, driver options, prototype findings, multi-GPU and MPI), `pressure/`, `combustion/`, `radiation/` and `solid/` (how each subsystem maps onto the AMR hierarchy), `inventory/` (a map of FDS mesh data structures and where they are used, with generated CSV tables) and `vv/` (verification and validation plan, case inventory and baseline status; `vv/archive/` keeps earlier plan versions).

## Status labels

Each document header names an owner and a status: *draft* (written, not reviewed), *in progress*, *reviewed* or *approved* (by the project owner). ADRs are Proposed, Accepted, Superseded or Rejected. Decision-log entries are marked proposed, accepted (with who accepted them) or superseded. Within the requirements, *proposed* marks a number nobody has confirmed, *TBD(owner)* a value the named lead still has to supply, and *ASSUMPTION* a scope assumption awaiting confirmation; none of these is a commitment until the project owner approves it.

## Conventions

- **FR-**, **NFR-** and **IR-** numbers are functional, non-functional and interface requirements (`requirements.md`). IDs are never reused.
- **D-** numbers are decisions and **A-** numbers are actions, both logged in `README.md`.
- **R-** numbers are risks (`risks.md`), rated by likelihood and severity.
- **Q1-Q12** are the open questions for the project owner (`charter.md` §9).
- **ADR-** numbers are architecture decision records (`adr/`).
- Also used: milestones M0-M11 (`roadmap.md`), tolerance classes T0-T3 (`requirements.md` §2.2) and prototypes P1-P3 (`amrex/driver-options.md`).

Paths written as `docs/...` inside the documents refer to this directory. Paths on the development machine were replaced by descriptions in parentheses, such as `(local build directory)`; where `(repo root)` remains inside a command, it stands for the top of this repository. Events are dated by day only: times of day, names of individuals and details of the development machine that do not affect results were removed from this snapshot.
