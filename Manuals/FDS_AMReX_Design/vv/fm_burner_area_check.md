# FM_Burner cylinder area check (FDS_AREA / AREA_ADJUST vs mesh split)

Owner: AMR V&V Lead · 2026-09-25 · Build: GNU FireX Release `36975d765f` (`(local build directory)/firex-36975d7/ompi_gnu_rel/fds`) · Related: requirements FR-005 (v) / D-028, upstream issue candidate #1 (`docs/inventory/upstream_issue_candidates.md`, not filed) · Work folder: `(local V&V run directory)/inputs/fm_burner_area/`

## Summary

- **The committed family-A FM_Burner cases (12 of the 21 inputs) are split-dependent.** Their burning area, and so the fuel flow and HRR, depend on how meshes map to MPI processes.
  - **1 mesh per process**, the layout upstream uses (`Run_All.sh`: `-p 12/96/208`): the burner top delivers exactly the specified fuel, **1.00000 × MASS_FLUX × πR²**.
  - **All meshes on 1 process:** it delivers **0.52083 × (= 25/48)**, the same factor at every resolution (2 cm, 1 cm, 5 mm) and for every fuel. For CH4 that is 8.46e-5 kg/s instead of 1.624e-4 kg/s, i.e. a nominal fuel-supply HRR of **4.23 kW instead of 8.13 kW**.
  - A merged **single-mesh** copy gives 1.00000.
- **Where they can be used as references:** only in the upstream 1-mesh-per-process layout, or as a single-mesh or MPI_PROCESS copy that puts each burner-top mesh on its own process. On this 8-core development machine the committed 12/96/208-mesh inputs can run only as np=1 (FDS ERROR(115): np must be 1 or ≥ number of meshes), which is the wrong (0.52×) layout. **Do not use them as FDS baselines here as committed.**
- **They are not layout-independent references for the AMR code** (FR-005 (v), D-028): an AMR run that sums areas exactly over the domain must match the 1.00000 single-mesh value, not any FDS multi-mesh value.
- **Family B (9 inputs, `C2H4_{15p2,16p8,20p9}_{2cm,1cm,5mm}`) is not subject to this defect.** The burner is a `&GEOM` (the OBST cylinder is commented out). It can't be set up here, because the `../../../../cad/…bingeom` files are missing (ERROR(705)), and GEOM/CC_IBM is out of scope anyway (FR-044).
- No FM_Burner inputs exist under `Verification/`. All 21 are under `Validation/FM_Burner/FDS_Input_Files/`.

## Inputs

| Family | Inputs | Burner | Meshes (upstream np) |
|---|---|---|---|
| A (12) | `FM_15cm_Burner_{CH4,C2H4,C3H6,C3H8}_{2cm,1cm,5mm}.fds` | MULT-voxel `&OBST SHAPE='CYLINDER'`, R=0.0685 m, H=0.254 m, XYZ 0,0,-0.254, `SURF_IDS='burner','wall','wall'` (only the top burns; `burner` is MASS_FLUX, e.g. CH4 0.011019529971700 kg/m²/s, C3H8 0.012138440573284) | 2cm: 12 (2×2×3 of 30³, dx 0.02, 324k cells, `-p 12`); 1cm: 96 (2.59M, `-p 96`); 5mm: 208 (80 at 1 cm + 128 at 5 mm, ratio 2; 5.6M; `-p 208`) |
| B (9) | `FM_15cm_Burner_C2H4_{15p2,16p8,20p9}_{2cm,1cm,5mm}.fds` | `&GEOM` burner and hood from `cad/` bingeom files | `-p 96` (2cm, 1cm), `-p 180` (5mm) |

Side observations, not area issues:
- The family-A C2H4 inputs have `&SPEC ID='EHTYLENE'`, which produces a warning.
- The input comment "D = 0.137 m (10 kW), HRRPUA=678.4 kW/m2" doesn't match the input itself. MASS_FLUX × FDS ΔHc (50,027 kJ/kg for the 2-step CH4 scheme) = 551.3 kW/m², i.e. 8.13 kW over πR².

## Method

1. **Setup-check all 21 inputs as committed** (`setup_check/`, 1 rank, T_END=0).
   - Family A: all 12 give "STOP: Set-up only".
   - Family B: ERROR(705), bingeom files not found.
   - A setup-only run can't show the adjusted area. T_END=0 sets SETUP_ONLY (read.f90:1510), which skips INITIALIZE_WALL_ARRAY and stops before INITIALIZE_MESH_VARIABLES_1 (main.f90:207, 238-248). So `ADJUST_OBST_SHAPE_AREA` (init.f90:339-432, called at init.f90:917) never runs, and neither the .out nor the .smv reports FDS_AREA or AREA_ADJUST.
2. **Face-count model from the setup .smv** (`burner_top_faces_from_smv.py`). For each mesh, it counts the voxel top faces (IOR +3, burner SURF) and the target area. The target is circle ∩ mesh for the top (read.f90:10962) and the whole side 2πRH (read.f90:10966). It then applies FDS's accumulation: `MULTIPLIER%FDS_AREA` accumulates per mesh in MULT order, K outer / J / I inner (read.f90:637-639), and is never reset or reduced across processes (init.f90:366-424).
   - In all 12 family-A cases the top faces sit in exactly 4 meshes (the quadrants), each with the same voxel area: 0.0032 m² at 2 cm, 0.0037 m² at 1 cm and 5 mm. Each has a target of πR²/4 = 0.003685 m².
   - The burner-top meshes are #1-4 (2cm), #22, 23, 26, 27 (1cm) and #118, 119, 122, 123 (5mm).
   - **All on 1 process:** the k-th top mesh in order sees a cumulative FDS_AREA of k quadrants, so AREA_ADJUST_k = 1/k. The effective top area is (1 + 1/2 + 1/3 + 1/4)/4 · πR² = **25/48 · πR² = 0.52083 πR²** (0.0076777 m²), independent of resolution.
   - **1 mesh per process:** each process sees only its own quadrant, so the result is exactly 1.
3. **Very short runs**, few time steps (`runs/`, T_END 0.03–0.1 s, 1 rank unless noted, machine load 3–5 of 8 cores). The HRR in `_hrr.csv` is still in the ignition transient (0.27–0.50 kW at 0.1 s), so it can't be compared with nominal. The **fuel MLR** in `_hrr.csv` is compared instead. It is the direct measure of the burning area: ratio = MLR_fuel / (MASS_FLUX · πR² · tanh(t/1 s)), where tanh is FDS's default mass-flux ramp (TAU = 1 s), which fits exactly. The pilot (`PART` 'ignitor', separate SURF) serves as a control and has identical MLR_PILOT FUEL in every layout.

## Results

| Run (input in the work folder) | Layout | Fuel MLR at t | Ratio to MASS_FLUX·πR²·tanh(t) |
|---|---|---|---|
| `CH4_2cm_committed` (committed input, CHID changed) | 12 meshes, all on 1 process | MLR_METHANE 8.4323491E-06 kg/s at 0.1 s | **0.520833** (the same at every output time: 0.0238, 0.0476, 0.0714, 0.0952 and 0.1 s) |
| `CH4_2cm_1mesh` (same domain and cells merged into one 60×60×90 mesh) | 1 mesh, 1 process | 1.6190110E-05 at 0.1 s | **1.000000** |
| `CH4_2cm_mpiproc4` (MPI_PROCESS copy: top meshes 1-4 on ranks 0-3, meshes 5-12 on rank 3) | 4 processes, each top mesh first on its process (equivalent to 1 mesh per process for this area) | 1.6190110E-05 at 0.1 s (identical to 1mesh at all times) | **1.000000** |
| `C3H8_2cm_committed` | 12 meshes, 1 process | MLR_PROPANE 4.6558715E-06 at 0.05 s | **0.520833** |
| `CH4_1cm_region8` (8 burner-region 1 cm meshes, MULT order kept) | 8 meshes, 1 process | 4.2266974E-06 at 0.05 s | **0.520833** |
| `CH4_5mm_region12` (12 burner-region 5 mm meshes, MULT order kept) | 12 meshes, 1 process | 2.5373703E-06 at 0.03 s | **0.520833** |

HRR at 0.1 s (CH4 2cm): 0.274 kW for the committed layout on 1 process, 0.501 kW for 1mesh, 0.487 kW for mpiproc4. The fire is only starting, but the 1-process run already lags by about 45%.

Nominal vs delivered (fuel supply × FDS ΔHc, steady state, i.e. after the 1 s ramp):

| Fuel | Nominal fuel flow MASS_FLUX·πR² | 1 mesh/process or single mesh | All meshes on 1 process |
|---|---|---|---|
| CH4 | 1.62440E-4 kg/s → 8.13 kW | 8.13 kW (×1.00000) | 4.23 kW (×0.52083) |
| C3H8 | 1.78934E-4 kg/s | ×1.00000 (predicted; 1-process run measured ×0.52083) | ×0.52083 |
| C2H4, C3H6 | — | ×1.00000 (face model) | ×0.52083 (face model, same meshing) |

Whether the rest of the fuel burns doesn't change the conclusion: the fuel supply itself is 48% short.

### Comparison with the Legacy Mapper reproducer

The upstream issue candidate #1 cases (R=H=0.1 m cylinder) showed a top area ×0.75 for 2 meshes on 1 process and ×1.0 for 2 meshes on 2 processes. The FM_Burner result follows the same mechanism with 4 top meshes: (1+1/2+1/3+1/4)/4 = 0.52083, where 2 meshes give (1+1/2)/2 = 0.75. The 0.75x–2x range in #1 is therefore not the worst case: splitting a burner top 4 ways on one process gives 0.52x.

### Side and bottom faces (inert 'wall'; predicted, not measured)

The side target is the whole cylinder side, 2πRH, not the in-mesh part. With 1 mesh per process each mesh scales its own side faces up to the whole side, so the effective side area is multiplied by the number of meshes holding side faces: ×4 (2cm), ×8 (1cm), ×12 (5mm). On 1 process the accumulation gives a different, still wrong, factor. This affects only heat exchange with the inert burner walls, not the fuel flow. Bottom faces (IOR -3) behave like the top. Both are predicted from source and the Legacy Mapper's reproducer (variant c), not measured here.

## Conclusion

1. **Committed upstream FM_Burner validation results** (run 1 mesh per process) have the correct burner-top area. Their HRR and fuel flow are not affected by this defect. Only the inert side-wall area (predicted ×4/×8/×12) is affected.
2. **Here**, the committed family-A inputs can run only as np=1, which gives 52.1% of the specified fuel and HRR (CH4: 4.23 vs 8.13 kW). They are **not safe as FDS baselines as committed**. If we need them, use one of these:
   - a single-mesh copy (feasible at 2 cm: 324k cells);
   - an MPI_PROCESS copy with each burner-top mesh first on its own process (2cm: np=4 is enough, as shown);
   - 1 mesh per process at 12 ranks. That exceeds the development machine's 8 cores and is not allowed under the current load rules.
3. **For the AMR code**, FM_Burner (family A) is a natural FR-005 (v) test input, a cylinder spanning boxes. The reference is the single-mesh value (ratio 1.00000). FDS multi-mesh runs are not a reference for setup areas. The AMR code must give 1.00000 for every box split, rank count and thread count.
4. **Family B** is unaffected by this defect, but it can't run here (missing `cad` repo) and is excluded as GEOM/CC_IBM (FR-044).
5. None of the FM_Burner inputs are in our Tier 1/Tier 2/anchor lists (`docs/vv/case_inventory.csv`), so no existing gate changes. D-030 already records that FM_Burner is unaffected by the {2,4} ratio limit (the 5mm inputs are 2:1). This note adds that it is affected by the area defect whenever it is run with more than one burner-top mesh on a process.

## Files (`(local V&V run directory)/inputs/fm_burner_area/`)

- `make_fm_burner_variants.py`: writes the committed, 1mesh and mpiprocN copies of the 2cm inputs. MPI_PROCESS must be contiguous and monotonic, otherwise ERROR(117).
- `make_reduced_region.py`: burner-region-only copies for 1cm/5mm, with MULT order kept. For the area check only; the 5mm copy warns "DEVC XO2 not within any mesh".
- `burner_top_faces_from_smv.py`: the face-count model.
- `setup_check/`: setup-only copies of all 21 inputs with .out/.smv/.log.
- `runs/`: short-run outputs (`*_hrr.csv`, `.out`, `.log`).
