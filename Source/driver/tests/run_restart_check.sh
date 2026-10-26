#!/bin/bash
# run_restart_check.sh <fds_amr> <c10 input as run_w12_smoke.sh> [work-dir] [ranks] [checkpoint step=30] [last step=40]
# Restart M3: c10 2 levels, bed particles, regrid every 5. A straight run to the last step against a run that writes the checkpoint at the given step and stops there,
# then continues with --restart-from (RESTART=T added to &MISC, same directory). PASS when the HRR, DEVC, mass and multi-level CSVs agree to 1e-11 (relative, floor 1e-300).
set -u
BIN=$1; IN=$2; W=${3:-$PWD/restart}; NP=${4:-1}; S=${5:-30}; E=${6:-40}
run() { (cd "$1" && OMP_NUM_THREADS=1 timeout 3600 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" c10.fds --two-level-run --bed -10 30 -20 20 0 20 \
    --regrid-every 5 --tag-temp 330 $2 >> log 2>&1; echo $? >> rc); }
clean() { rm -f "$1"/*.s3d* "$1"/*.sf "$1"/*.bf "$1"/*.prt5*; }   # stock restart appends to these: kept until the restart has run
rm -rf "$W"; for d in s c; do mkdir -p "$W/$d"; cp "$IN" "$W/$d/c10.fds"; cp "$(dirname "$IN")/vegetation_model.txt" "$W/$d/"; done
run "$W/s" "--steps $E"; clean "$W/s"
run "$W/c" "--steps $S --checkpoint-at $S"
sed -i 's/&MISC /\&MISC RESTART=.TRUE., /' "$W/c/c10.fds"
run "$W/c" "--steps $E --restart-from $(cd "$W/c" && ls -d chk_*_s$S)"
clean "$W/c"; rm -rf "$W"/c/chk_* "$W"/c/*.restart
r=$(python3 - "$W" <<'PY'
import sys,glob,os
w=sys.argv[1]; worst=0.0; out=[]
for f in sorted(sum([glob.glob(w+'/s/*'+k+'.csv') for k in ('_hrr','_devc','_mass','multi_level')],[])):
    g=w+'/c/'+os.path.basename(f)
    if not os.path.exists(g): out.append(os.path.basename(f)+' missing'); worst=1; continue
    a=[l.split(',') for l in open(f)]; b=[l.split(',') for l in open(g)]
    m=0.0
    if len(a)!=len(b): out.append('%s rows %d/%d'%(os.path.basename(f),len(a),len(b))); worst=1
    for x,y in zip(a,b):
        for p,q in zip(x,y):
            try: p=float(p); q=float(q)
            except ValueError: continue
            m=max(m,abs(p-q)/max(abs(p),abs(q),1e-300))
    worst=max(worst,m); out.append('%s %.2e'%(os.path.basename(f),m))
print(('PASS' if worst<=1e-11 else 'FAIL')+' max rel %.2e: '%worst+', '.join(out))
PY
)
echo "$r" | grep -q ^PASS && [ "$(tr -d '\n' < "$W/s/rc")$(tr -d '\n' < "$W/c/rc")" = 000 ] && { echo "PASS restart np$NP ($S->$E): ${r#PASS }"; rm -rf "$W"; exit 0; }
echo "FAIL restart np$NP rc $(tr -d '\n' < "$W/s/rc")/$(tr -d '\n' < "$W/c/rc"): $r"; exit 1
