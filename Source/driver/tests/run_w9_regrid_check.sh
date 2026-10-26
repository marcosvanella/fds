#!/bin/bash
# run_w9_regrid_check.sh <fds_amr> <c10 cut-down C064 input (4 meshes at 1 m, &AMR line, vegetation_model.txt beside it)> [work-dir] [ranks] [steps]
# P0 W9: regrids that change the hierarchy over a burning particle bed, without the --bed level-1 pin: 2 levels tagged on temperature (--tag-temp 320), regrid
# every 5 steps, bed particles from the input. PASS: at least 2 changed regrids; on the 2 steps after each changed regrid every level's max|div u - D| <= 1e-9
# (two_level_multi_level.csv) and the mass-balance residual change <= 1e-6 x (initial mass + through-flow so far) per step; W8 count and hash SAME at every
# regrid; TWO-LEVEL gate PASS (fine div, mass and species balance, species check 7 over all steps).
set -u
BIN=$1; IN=$2; W=${3:-$PWD/w9_regrid}; NP=${4:-1}; N=${5:-40}
rm -rf "$W"; mkdir -p "$W"; cp "$IN" "$W/c10.fds"; cp "$(dirname "$IN")/vegetation_model.txt" "$W/"
cd "$W" && FDSTL_W8_CHECK=1 OMP_NUM_THREADS=1 timeout 1800 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" c10.fds --two-level-run --levels 2 --tag-temp 320 \
    --regrid-every 5 --steps "$N" > log 2>&1; rc=$?
rm -f *.s3d* *.sf *.bf *.prt5* *.restart
python3 - "$rc" "$NP" <<'P'
import csv,re,sys
rc,np=int(sys.argv[1]),sys.argv[2]; t=open('log',errors='ignore').read(); why=[]
if rc: why.append('rc %d'%rc)
if 'RESULT gate PASS' not in t: why.append('gate not PASS')
ch=[int(m.group(1)) for m in re.finditer(r'mid-regrid after step (\d+) t \S+ changed 1',t)]
if len(ch)<2: why.append('changed regrids %d < 2'%len(ch))
w8=re.findall(r'W8 CHECK .*',t)
if not w8 or any('SAME' not in x for x in w8): why.append('W8 checks %d not all SAME'%len(w8))
m=re.search(r'composite mass (\S+) ->',t); m0=float(m.group(1)) if m else 0
r=list(csv.DictReader(open('two_level_multi_level.csv'))); by={int(x['step']):x for x in r}
dmax=0.0; jmax=0.0
for s in ch:
    for k in (s+1,s+2):
        if k not in by or k-1 not in by: continue
        x=by[k]; dmax=max([dmax]+[abs(float(v)) for c,v in x.items() if c.startswith('div_l')])
        res=lambda y: float(y['mass'])-m0-float(y['inflow'])
        jmax=max(jmax,abs(res(x)-res(by[k-1])))
lim=1e-6*m0   # the mass part of the gate limit (through-flow adds to it)
if dmax>1e-9: why.append('div after regrid %.3e'%dmax)
if jmax>lim: why.append('mass residual jump after regrid %.3e kg > %.3e'%(jmax,lim))
s='np%s: %d changed regrids (after steps %s), max div on the 2 steps after %.3e, max residual jump %.3e kg (limit %.3e), W8 checks %d SAME'%(np,len(ch),ch,dmax,jmax,lim,len(w8))
print(('PASS w9_regrid ' if not why else 'FAIL w9_regrid ')+s+('' if not why else ': '+'; '.join(why))); sys.exit(1 if why else 0)
P
