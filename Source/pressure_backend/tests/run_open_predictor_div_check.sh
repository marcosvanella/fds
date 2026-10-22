#!/bin/bash
# Predictor divergence on a 2-level OPEN-face case (C064 cut-down c10 with the grass bed, wind, OPEN XMIN/XMAX/YMIN/YMAX/ZMAX, level 1 over the
# grass bed, fine boxes on YMIN/YMAX). The fine OPEN-face H ghost must use the Dirichlet data folded into the same composite solve; the
# level-1 max|div(US)-DS| after the predictor (FDSTL_GHOSTDBG DIVERR) must be at solver level, like the corrector.
# usage: run_open_predictor_div_check.sh <fds_amr binary> [work-dir] [steps=4]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); BIN=$(readlink -f "$1"); W=${2:-$PWD/open_pred_div}; N=${3:-4}; rm -rf "$W"; mkdir -p "$W"; cd "$W"
cp "$HERE/../../../Validation/CSIRO_Grassland_Fires/FDS_Input_Files/vegetation_model.txt" . || exit 1
cat > o.fds <<'EOC'
&HEAD CHID='o', TITLE='Australian Grassland Fire Experiment, Case C064' /
&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /
&MESH IJK=20,20,20, XB=-10.0,10.0,-20.0,0.0,0.0,20.0, MULT_ID='mesh' /
&MULT ID='mesh', DX=20., DY=20., DZ=20., I_LOWER=0, I_UPPER=1, J_LOWER=0, J_UPPER=1, K_UPPER=0 /
&TIME T_END=3. /
&MISC SIMULATION_MODE='SVLES', TMPA=32. /
&WIND SPEED=4.6, Z_REF=2., DIRECTION=270., L=-500., Z_0=0.03 /
&SURF ID='IGNITION', HRRPUA=1500., COLOR='RED', TMP_FRONT=100., RAMP_Q='ignitor' /
&RAMP ID='ignitor', T= 0., F=0. /
&RAMP ID='ignitor', T= 1., F=1. /
&RAMP ID='ignitor', T= 5., F=1. /
&RAMP ID='ignitor', T= 6., F=0. /
&VENT XB=0.0,1.6,-10.,10.,0.,0., SURF_ID='IGNITION', XYZ=0,0,0, SPREAD_RATE=1. /
&VENT MB='XMIN', SURF_ID='OPEN' /
&VENT MB='XMAX', SURF_ID='OPEN' /
&VENT MB='YMIN', SURF_ID='OPEN' /
&VENT MB='YMAX', SURF_ID='OPEN' /
&VENT MB='ZMAX', SURF_ID='OPEN' /
&SURF ID='ROAD'
      COLOR='BEIGE'
      MATL_ID='DIRT'
      THICKNESS=0.2 /
&VENT XB= -3,  0,-20, 20,0,0, SURF_ID='ROAD' /
&MATL ID='DIRT'
      CONDUCTIVITY = 0.25
      SPECIFIC_HEAT = 2.
      DENSITY = 1300. /
&SURF ID='GROUND'
      DEFAULT=T
      MATL_ID='DIRT'
      THICKNESS=0.2 /
&SURF ID                        = 'wet vegetation'
      MATL_ID(1,1:1)            = 'GENERIC VEGETATION'
      MOISTURE_CONTENT         = 0.063
      SURFACE_VOLUME_RATIO      = 9770.
      LENGTH                    = 0.20
      GEOMETRY                  = 'CYLINDRICAL' /
&CATF OTHER_FILES='vegetation_model.txt' /
&PART ID='blade of grass'
      SAMPLING_FACTOR=5
      SURF_ID='wet vegetation'
      PROP_ID='grass image'
      DRAG_COEFFICIENT=2.8
      QUANTITIES='PARTICLE TEMPERATURE','PARTICLE MASS'
      STATIC=T
      COLOR='BROWN' /
&INIT PART_ID='blade of grass', XB=0.,25.,-15.,15.,0.00,0.21, N_PARTICLES_PER_CELL=1, CELL_CENTERED=T, MASS_PER_VOLUME=1.33, DRY=T /
&PROP ID='grass image', SMOKEVIEW_ID='TUBE', SMOKEVIEW_PARAMETERS='L=0.20','D=0.0004','DIRZ=1' /
&DUMP DT_HRR=0.1, DT_DEVC=0.1, MASS_FILE=T, DT_MASS=0.1, SIG_FIGS=15, SIG_FIGS_EXP=3 /
&DEVC ID='x', QUANTITY='TEMPERATURE', SPATIAL_STATISTIC='MAXLOC X', XB=-10,30,-1,0,0,1 /
&DEVC ID='D', QUANTITY='FIRE DEPTH', XB=-10,30,0,0,0.1,0.1 /
&DEVC ID='T1', QUANTITY='TEMPERATURE', XYZ=2.0,0.,0.5 /
&DEVC ID='T2', QUANTITY='TEMPERATURE', XYZ=5.0,0.,1.0 /
&DEVC ID='U1', QUANTITY='U-VELOCITY', XYZ=5.0,0.,2.0 /
&DEVC ID='W1', QUANTITY='W-VELOCITY', XYZ=2.0,0.,2.0 /
&TAIL /
EOC
FDSTL_GHOSTDBG=1 OMP_NUM_THREADS=1 timeout 900 "$BIN" o.fds --two-level-run --levels 2 --blocking 4 --outdir . --steps $N --chid o > run.out 2> run.err || { echo "FAIL: run failed (see $W/run.err)"; exit 1; }
python3 - "$W/run.out" <<'EOP'
import re, sys
p = []; c = []
for l in open(sys.argv[1], errors='replace'):
    m = re.match(r'DIVERR after (predictor|corrector) \S+ \(L0 \S+, L1 (\S+)\)', l)
    if m: (p if m.group(1) == 'predictor' else c).append(float(m.group(2)))
if not p or not c: print("FAIL: no DIVERR lines"); sys.exit(1)
lim = 1e-8
ok = max(p) < lim and max(c) < lim
print(("PASS" if ok else "FAIL") + ": L1 max|div(US)-DS| after predictor %.3e, after corrector %.3e (limit %g, %d steps)" % (max(p), max(c), lim, len(p)))
sys.exit(0 if ok else 1)
EOP
