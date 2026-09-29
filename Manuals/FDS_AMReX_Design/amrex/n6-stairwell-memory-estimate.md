# N6: memory of tall stairwell cases under the non-box level-0 ruling

Answers ruling N6 (`adr/drafts/ruling-nonbox-level0.md:15`). Estimate only: arithmetic (`prototypes/n6_mem/n6_estimate.py`,
output `n6_estimate.out`) plus one tiny instrumented AMReX run (§5). Source pins: FireX `36975d765f` (`src/Source`), AMReX `99ddfda`.
Labels: **[src]** read from source, **[meas]** measured in the tiny run, **[der]** arithmetic, **[A#]** assumption.

## 0. Verdict
- **The covering level 0 does not fit `stairwell` on the 16 GB development machine.** Level 0 alone is 15.1–19.7 GB (NS = 1, 16³ boxes) against NFR-031's
  12 GB cap (`requirements.md:442`). About 90 % of that is gap cells. **Omitting all-gap boxes brings it to 2.4–3.2 GB**, so the N6 fallback is needed.
- **Both N6 conditions can be met.** (i) `stairwell` has no gap-face `OPEN` vent: its `OPEN` vent (`stairwell.fds:451`, y = −7.0) is on the
  bounding box's low-y face, and the `Extract` vent (`:455`, z = 31.0) is on the top face. Both are domain faces **[der, from coordinates; snapping not run]**.
  (ii) AmrCore supports a level 0 that does not cover the domain through a documented hook (§3). MLMG handles it too, but with one gap:
  keeping the overset gap mask in the partly-gap boxes is outside what MLMG asserts (§3).
- **Terms and assumptions.** "All-gap box" [A1]: a box of the level-0 BoxArray (`MakeBaseGrids`, run's `max_grid_size`) with no cell inside any
  `&MESH`. "Fallback" [A2]: exactly the ruling's wording. The kept boxes still contain gap cells, which stay masked as in N1/N3. It is not a BoxArray
  that follows the mesh union (that variant is noted in §3).

## 1. Per-cell byte model
AMReX storage [src]: a FAB holds `nvar × numPts(grown box) × sizeof(T)` (`AMReX_BaseFab.H:1353,1379`), with the grown box
`grow(ba[K], n_grow)` (`AMReX_FabArrayBase.cpp:217-219`). Real = 8 B, int (`iMultiFab`) = 4 B, and ghost cells cost the same as valid cells.
Only boxes in the BoxArray are allocated. NS = `N_TOTAL_SCALARS` (`read.f90:3032`); `stairwell` has no `&SPEC`/`&REAC`, so NS = 1.

| Group (ng) | Fields, source | Reals/cell |
|---|---|---|
| Density, ng=3 | `RHO`, `RHOS` (`init.f90:525-526`); ng=3 native layout (`adr/ADR-001-driver-architecture.md:308`) | 2 |
| Cell state, ng=2 | `TMP` (`init.f90:524`); `ZZ`, `ZZS` (`:527-528`), species ng=2 (same line) | 1 + 2NS |
| Faces, ng=1, nodal | `U,V,W,US,VS,WS` (`:533-538`; ng=1 per `p1-findings.md:71`) | 6 |
| Cell, ng=1 (persistent) | `H,HS,KRES,DDDT,D,DS,MU,MU_DNS` (`:556-563`), `STRAIN_RATE` (:569), `Q,MIX_TIME` (:581-582), `RSUM` (:599), `D_SOURCE` (:670), `CHI_R,QR,KAPPA_GAS,UII` (:677-680), `LES_FILTER_WIDTH` (:799), `FVX..FVZ, FVX_B..FVZ_B` (:543-548), `PRHS` (:2445); `DEL_RHO_D_DEL_Z` (:590), `M_DOT_PPP` (:671) | 25 + 2NS |
| Radiation, ng=1 | `UIID` (`init.f90:854`), UIIDIM = min(5, 100/15) = 5 by default (`read.f90:10225,10301-10302`) | 5 |
| Scratch (ng 1–2) | `WORK_PAD` (:529), `WORK_U/V/W` (:539-541), `WORK1-9` (:691-699), `SWORK1-4` (:702-705), `FX/FY/FZ` (:611-613, per-box scratch, `p1-findings.md:72`) | 13 + 7NS |
| Old state, ng=0 [A3] | ADR-002 option C keeps per-level old/new state (`adr/ADR-002-time-stepping.md:118-119`); assumed `RHO`, `ZZ`, U/V/W | 4 + NS |
| Integers (4 B) | `PRESSURE_ZONE` (:586), `CELL_INDEX` (`read.f90:11622`), `IWORK1` (:701), ng=1; `CHEM_ACTIVE_CELLS` ×4 (:578); face mask, 8 comps, ng=2 (`p1-findings.md:381`), plus 1 gap mask [A4] | 16 ints |
| MLMG, masked single level | `rhs`, `res`, `rescor`, `cor` (ng=1), `cor_hold` (ng=1) (`AMReX_MLMG.H:1522,1561-1562,1585,1601`; `sol` aliases `H` when ng matches, :1499-1503); `a` (`AMReX_MLABecLaplacian.H:356`), 3 face `b` (:363); overset mask per MG level (`AMReX_MLCellABecLap.H:183,247-256`); × 8/7 for the MG hierarchy | 9 × 8/7, + 1 int |
| MLMG, measured extra | boundary registers and masks (`AMReX_MLCellLinOp.H:436-496`); the tiny run measured 129–150 B/cell against 92–98 modelled | +40 B |

Species add 12 reals per cell per species. Not modelled: optional arrays (`PARTICLE_DRAG`, `CHECK_VN`, droplet and `STORE_*` arrays), wall/OBST records,
and the P1 shim's per-cell side data. The shim's `CELL_TYPE` is 384 B/cell and its walls 540 B/entry, which gives **+760–1000 B/cell** (`p1-findings.md:124-126`).

**Result [der; the field set is confirmed to the byte by the run, §5].** B = bytes per allocated cell; lean = scratch not persistent (allocated per box/tile):

| NS | valid only | 16³ boxes: full / lean | 32³ boxes: full / lean | ghost overhead 16³ / 32³ |
|---|---|---|---|---|
| 1 | 735 | **1042** / 808 | **877** / 682 | 1.42 / 1.19 |
| 3 | 927 | 1326 / 932 | 1111 / 782 | 1.43 / 1.20 |
| 6 | 1215 | 1751 / 1118 | 1462 / 932 | 1.44 / 1.20 |

Breakdown at NS = 1 and 16³ (B/cell): persistent ng=1 cell fields 308, scratch 234, MLMG 138, integers 103, density and cell state 89, faces 72,
radiation 57, old state 42.
- **Flux registers.** `YAFluxRegister` allocates `nvar` components plus an int flag on the **whole coarse BoxArray** (`AMReX_YAFluxRegister.H:218,220`).
  A covering level 0 therefore pays (1+NS)·8 + ~6 B per gap cell once a level 1 exists. The classic `FluxRegister` is a `BndryRegister` on the fine
  boundary only (`AMReX_FluxRegister.H:22-24`), so it is negligible. MLMG's composite solve always adds one `YAFluxRegister` (`AMReX_MLCellLinOp.H:496`).
- **Metadata.** A `CopyComTag` is two Boxes and two ints (`AMReX_FabArrayBase.H:210-215`, about 64 B). The FB/CPC caches (`:565-605`, `:628-670`)
  hold about 26 tags per box per key, which is about 8 MB per key for `stairwell`'s 4,620 boxes. BoxArray/DistributionMapping are about 30 B per box.
  Both are below 1 % **[der]**.

## 2. Representative cases (level-0 MultiFab bytes)
Case selection. `stairwell` (`Verification/Pressure_Solver/stairwell.fds:5-15`) is the ruling's case. `hallways` (`hallways.fds:5-9`) is the other real
low-fill input, and its gap-face `OPEN` vent (`:15`) excludes it from the fallback. NRCC Smoke Tower `BK-R` (`Validation/NRCC_Smoke_Tower/.../BK-R.fds:6-10`)
is a real 10-storey stair tower whose `MULT` meshes union to a box (fill 1), so it is the control. S1–S3 are synthetic parametric cases, with meshes
deliberately not aligned to 16-cell box edges.

`stairwell` meshes (dx = 0.1 m; cell offsets from `vv-runs/inputs/A-46/stairwell_np7_check.txt`), cell counts in thousands:
59×27×178 (284), 70×37×35 (91), 54×37×22 (44), 3 × 61×26×95 (151 each), 17×81×30 (41), 49×173×30 (254), 76×51×34 (132), 82×78×29 (185), 67×77×30 (155).
That is 1,637,828 covered cells. The bounding box is 189×173×549 = 17,950,653 cells, so **fill = 0.091**.

GB = 10⁹ B. Full model, lean model in brackets. "kept" = cells in kept boxes as a fraction of the bounding box. Gap share = the part of the covering
level 0 spent on gap cells (≈ covering × (1 − fill)).

| Case (NS, box) | Level-0 domain | N_bb | Covered | Fill | Kept (16³) | Covering L0 GB | Fallback L0 GB | Gap share of covering |
|---|---|---|---|---|---|---|---|---|
| `stairwell` exact bb (1, 16³) | 189×173×549 | 17.95 M | 1.64 M | 0.091 | 0.162 | **19.7 (15.5)** | **3.19 (2.50)** | ≈ 17.9 GB |
| `stairwell` padded to bf 8 [A5] (1, 16³) | 192×176×552 | 18.65 M | 1.64 M | 0.088 | 0.162 | 19.5 (15.1) | 3.15 (2.44) | ≈ 17.8 GB |
| same (1, 32³) | 192×176×552 | 18.65 M | | | 0.223 | 16.5 (12.8) | 3.67 (2.86) | ≈ 15.0 GB |
| same (3, 16³) | 192×176×552 | 18.65 M | | | 0.162 | 24.8 (17.4) | 4.00 (2.82) | ≈ 22.6 GB |
| `hallways` (1, 16³); fallback not allowed | 128×64×64 | 0.52 M | 73,728 | 0.141 | 0.141 | 0.55 (0.42) | (0.08) | 0.47 GB |
| NRCC `BK-R` (3, 32³), box domain | 168×112×306 | 5.76 M | 5.76 M | 1.000 | 1.000 | 6.52 (4.59) | same | 0 |
| S1: 18×18 shaft H=512 beside a 96×64×16 floor (3, 16³) | 96×64×512 | 3.15 M | 0.26 M | 0.082 | 0.193 | 4.34 (3.10) | 0.84 (0.60) | ≈ 4.0 GB |
| S2: 36×28 shaft H=512 over a 120×120×30 lobby (3, 16³) | 120×120×512 | 7.37 M | 0.92 M | 0.124 | 0.180 | 10.3 (7.38) | 1.85 (1.33) | ≈ 9.0 GB |
| S3: C-shaped stacked stair, corridor length 64/128/256 (3, 16³) | L×30×510 | 0.98/1.96/3.92 M | | 0.37/0.24/0.18 | 0.56/0.34/0.23 | 1.31/2.62/5.24 | 0.73/0.89/1.21 | 0.8/2.0/4.3 GB |

**Refinement scenarios** (ratio 2; level 1 lies inside the refinable region, which excludes gaps (IR-008, N1), so it costs the same in both designs):

| Case | Level 1 | Fine cells | L1 GB | Total covering | Total fallback |
|---|---|---|---|---|---|
| `stairwell` padded (1, 16³) | upper stair shaft (mesh 1) | 2.27 M | 2.03 (1.58) | 22.2 (17.4) | 5.29 (4.13) |
| S1 H=512 | fire in the floor (32×32×16) + whole shaft | 1.42 M | 1.71 | 6.05 | 2.55 |
| S2 H=512 | lobby core (40×40×30) + whole shaft | 4.27 M | 4.89 | 15.2 (10.8) | 6.74 (4.76) |

Observations **[der]**
- Kept-box fraction ≈ 1.3–2.3 × fill with 16³ boxes. Thin shafts (18 cells wide, not box-aligned) straddle box edges. Larger boxes keep more gap
  (stairwell: 0.162 at 16³, 0.22–0.25 at 32³). `stairwell`'s kept boxes still hold 1.27–2.8 M gap cells.
- A covering level 0 pays about 1 KB per gap cell. `stairwell` covering costs about 11.9 KB per gas cell; the fallback costs 1.9 KB.
- Refinement dilutes the saving (S2: 8.3× on level 0, 2.3× in total), but it never removes it.
- **NFR-031's 1.3× baseline rule is also at stake.** FDS itself holds roughly 1 KB per gas cell (about 63 reals plus `CELL_TYPE` 384 B, `p1-findings.md:125`),
  so the `stairwell` baseline is ≈ 1.7–2 GB **[estimate, not measured]**. Covering is ≈ 8–11× baseline; the fallback is ≈ 1.2–1.9×. In general a
  covering level 0 stays within 1.3× baseline only when fill ≳ 0.75. Many of the 40 masked cases are below that, even though they are small in
  absolute terms (≤ 0.52 M bounding-box cells each except `stairwell`, `pressure/04-masked-domain-pressure.md:52,222`).

## 3. Can a covering level 0 drop all-gap boxes? (the key point)
**The ruling** allows it only as the N6 fallback, and only after AmrCore support is confirmed and for a case with no gap-face `OPEN` vent
(`ruling-nonbox-level0.md:15`). §3 of the ruling lists the non-covering level 0 as "kept only as the N6 fallback" (`:25`).

**AmrCore/AmrMesh: supported [src]**
- `AmrMesh::PostProcessBaseGrids(BoxArray&)` is a virtual hook documented "for example … to remove covered grids on the coarsest refinement level"
  (`AMReX_AmrMesh.H:412-421`). It runs on the level-0 BoxArray after chopping and before `MakeNewLevelFromScratch` (`AMReX_AmrMesh.cpp:574`;
  `AMReX_AmrMeshGridding.cpp:300`). Removing all-gap boxes there is the sanctioned route.
- Proper nesting for finer levels is built from the complement of `grids[lbase]`, not from the domain (`AMReX_AmrMesh.cpp:631-644`). Level 1 therefore
  stays inside the kept level-0 boxes. The domain-divisibility check still applies to the Geometry domain (`:1252-1259`).
- What the driver must add: ghost cells of kept boxes that fall in dropped boxes are inside the domain. `FillBoundary` has no source for them, and
  `PhysBCFunct` fills only ghost boxes that are not contained in the grown domain (`AMReX_PhysBCFunct.H:208-233`). So after every FillPatch the driver
  has to set them to the gap-wall ghost values itself. In the covering design they are valid solid cells that the N2 gap-wall code already maintains.

**MLMG: works, but mask + non-covering is outside the asserted contract [src + meas]**
- `MLLinOp` detects a non-covering level 0 (`m_domain_covered[0]`, `AMReX_MLLinOp.H:1238-1240`) and sets `m_needs_coarse_data_for_bc` (:1514).
  `setCoarseFineBC(nullptr, r, Neumann)` (:251-266) then gives homogeneous Neumann on the uncovered boundary. That is exactly the gas–gap no-flux
  of N2/N3, but it is one BC type for the whole uncovered boundary, which is why gap-face `OPEN` vents are excluded. Agglomeration is off unless the
  BoxArray fills its minimal box (`:1261-1264`).
- The kept boxes still contain gap cells, which need the overset mask. In the non-covering + Neumann C/F branch, both `MLABecLaplacian` and
  `MLPoisson` do `AMREX_ASSERT(m_overset_mask[0][0] == nullptr)` (`AMReX_MLABecLaplacian.H:847-850`, `AMReX_MLPoisson.H:233-236`). The assert is
  active in debug builds only. `(local AMReX install)` is built with `AMReX_ASSERTIONS OFF` (`AMReXConfig.cmake:179`).
- Measured in release (§5): the combination runs, but MLMG flags the problem singular (`isSingular(0)=1`) even though a pinned cell exists. It
  converges more slowly (13–18 iterations against 3), and it agrees with the covering solution **only after gas-mean removal** (max difference
  2.7–3.5e-9 at a relative tolerance of 1e-10; raw difference 0.11–0.14, a constant). The D-032 gauge (per-zone mean removal) makes this acceptable
  for a single sealed component. With several components or zones, MLMG's one global singular fix is wrong in principle **[der]**.
- Options, in order:
  - (a) keep the mask and patch or relax the assert, with a test (an upstream question);
  - (b) as an unruled variant, build the level-0 BoxArray from the `&MESH` union chopped by `max_grid_size`, so there are no gap cells, no mask and
    no pin (option a1 in `pressure/04-masked-domain-pressure.md:155-161`). Memory then ≈ fill × covering ≈ 1.7–2.0 GB for `stairwell`, but boxes
    follow the mesh extents (27-, 37- and 26-cell sides give 13–14-cell boxes, below the 16-cell guidance in `requirements.md:443`).
- A mask-free trick tried in the same run (a = 1 in gap cells) aborted with "MLMG failed". It was not investigated.

## 4. When the fallback is worth it
Memory model **[der]**: M_cov ≈ b·N_bb + M_L≥1 and M_fb ≈ b·f_box·N_bb + M_L≥1. Here b ≈ 1.0 KB per allocated level-0 cell (NS = 1, 16³, full;
0.68–1.75 KB across NS, box size and lean/full), and f_box ≈ 1.3–2.3 × fill at 16³ boxes.
Usable budget U: 12 GB on the 16 GB workstation (NFR-031). On a GPU, the arena's 3/4 of device memory (`AMReX_Arena.cpp:427`), divided by the
ranks sharing the GPU.

| Budget | U | Largest N_bb for covering at 0.75·U (b = 1.04 KB) | `stairwell` covering (15–20 GB) | `stairwell` fallback (2.4–3.2 GB) |
|---|---|---|---|---|
| 16 GB workstation | 12 GB | 8.6 M | no | yes (20–27 % of U) |
| 8 GB GPU | 6 GB | 4.3 M | no | yes (40–53 %) |
| 16 GB GPU | 12 GB | 8.6 M | no | yes |
| 24 GB GPU | 18 GB | 13 M | lean only, above 75 % of U | yes |
| 80 GB GPU | 60 GB | 43 M | yes | yes |

**Rule.** Enable the fallback for a case when **covered cells are below 30 % of the bounding box** (so f_box ≲ 0.5 and level 0 at least halves) **and the
covering estimate exceeds 75 % of U**. The 25 % headroom is for memory the FAB model leaves out (§5). Above a fill of about 0.5 the fallback saves ≲ 30 %,
which does not justify a second code path. On the 16 GB development machine only `stairwell` qualifies among the 40 cases (the next largest is ≤ 0.52 M bounding-box
cells, about 0.55 GB). On an 8 GB GPU the rule triggers from about 4.3 M bounding-box cells.

## 5. Caveats, the tiny run, and what would firm this up
**Tiny run [meas]** (`prototypes/n6_mem/src/n6_mem.cpp`). Domain 64×64×128 with gas = a 64×64×16 lobby plus an 18×18×112 shaft (fill 0.194).
The field set of §1 at NS = 1 (MultiFabs with the modelled ng and nodality) and a masked `MLABecLaplacian` (β = 0 on gas–gap faces, one pinned cell,
Neumann domain). Covering: 128 boxes at 16³, 16 at 32³. All-gap boxes dropped: 44 and 7.
```
cmake -S . -B build -DCMAKE_PREFIX_PATH=(local AMReX install) -DHYPRE_ROOT=(local GNU third-party library tree)/libs/hypre/v3.0.0 \
      -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_BUILD_TYPE=Release && cmake --build build -j4
OMP_NUM_THREADS=1 mpirun -np 4 --bind-to none build/n6_mem max_grid_size=16|32 [bottom=hypre]   # each run < 5 s
```
- `TotalBytesAllocatedInFabs` for the field set: 903.95 B per allocated cell at 16³ and 744.25 at 32³. Both layouts give the same value, and the
  script predicts 904 and 745.
- MLMG FAB high-water above the pre-solve level (`TotalBytesAllocatedInFabsHWM`, summed over ranks): 129–150 B/cell with the default bottom
  solver and 166–205 B/cell with the HYPRE bottom solver. HYPRE's own matrix is not in FABs and was not measured.
- Solver finding for the Pressure Solver Lead (O3): with β = 0 gas–gap faces and the default BiCGStab bottom solver, **the covering masked solve
  returned NaN** (`run_mgs16.log`, `run_mgs32_bicg.log`). The HYPRE bottom solver converged in 3 iterations. With β = 1 everywhere, BiCGStab converged.
  The covering path therefore needs its bottom solver chosen or fixed, independent of N6.

**Caveats**
- The model counts FAB data only. Not included: wall/OBST records and `CELL_TYPE` side data (+760–1000 B/cell while the P1 shim is used, which
  would push `stairwell` covering to about 30–38 GB and the fallback to about 5–6 GB), particles, output snapshot buffers (ADR-004), MPI buffers, arena
  fragmentation, and HYPRE matrices. HYPRE is ≈ 250–500 B/cell when its level is large (`pressure/02-performance-expectations.md:128-129`).
- Which FDS arrays stay persistent (full vs lean, a 20–35 % spread) and the option-C old-state set [A3] are not decided.
- [A5] Domain divisibility. `stairwell`'s bounding box (189×173×549) is divisible by no `blocking_factor` ≥ 2 (A-46 check file), so AMReX needs
  bf = 1 (`AMReX_AmrMesh.cpp:1252-1259`). With every dimension odd, MLMG cannot coarsen at all, so the bottom solver runs on 18 M cells (+0.9 GB for
  BiCGStab's 6 MultiFabs, `AMReX_MLCGSolver.H:194-209`; far more if HYPRE). Padding the high sides to 192×176×552 fixes this and keeps the `OPEN`
  vent on a domain face, but it makes the top `Extract` vent a gap-face vent (allowed by N2 since it is not `OPEN`) and departs from N1's
  "bounding box". This needs a ruling.
- The kept-box fraction depends on `max_grid_size` and on how the meshes align with box edges.

**To firm up (all cheap)**
1. The FDS baseline peak RSS for `stairwell` (needs about 10 ranks, so it was not run under the 4-core limit). This fixes the NFR-031 ratio.
2. Instrument the P1 shim on a 1–2 M-cell masked case for real side-data bytes per cell.
3. Decide full vs lean scratch and the option-C field list.
4. Ask upstream whether the assert in `MLABecLaplacian.H:850` can be relaxed for overset + Neumann C/F, or adopt variant (b).
