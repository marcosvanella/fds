#!/bin/bash
# Fine ground stress on a seam: the level-set ground case with one static level-1 patch run as several fine boxes (MAX_GRID_SIZE 8: seams across the
# ground) and as one box (32), radiation and initial noise off. A ground face shared by two boxes gets the same stress from both, computed from ghost-filled data, so
# the two layouts must agree: every device (c_devc.csv) within 1e-11 relative over STEPS steps. The case needs devices; the default cases/ground_seam_c064.fds (trimmed C064, the gate's ref_L2_nf0 input) has them.
# usage: run_ground_stress_seam_check.sh <fds_amr binary> [work-dir] [steps=6] [mpirun prefix, e.g. "mpirun -np 1"] [case=cases/ground_seam_c064.fds]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); BIN=$(readlink -f "$1"); W=${2:-$PWD/ground_seam}; N=${3:-6}; MPI=${4:-}; CASE=$(readlink -f "${5:-$HERE/cases/ground_seam_c064.fds}")
rm -rf "$W"; mkdir -p "$W"
for m in 8 32; do
  d=$W/m$m; mkdir -p $d; cp "$CASE" $d/c.fds; cp "$(dirname "$CASE")"/vegetation_model.txt $d/ 2>/dev/null
  if grep -q '^&MISC' $d/c.fds; then grep -q 'NOISE' $d/c.fds || sed -i '0,/^&MISC/s//\&MISC NOISE=F,/' $d/c.fds; else echo "&MISC NOISE=F /" >> $d/c.fds; fi   # INITIAL_NOISE is drawn per fine box (e2c7779431): layout-dependent by design
  grep -q '^&RADI' $d/c.fds || echo "&RADI RADIATION=F /" >> $d/c.fds   # the radiation solve depends on the box layout by itself (Q_RADI)
  grep -q '^&AMR' $d/c.fds || sed -i '1a &AMR MAX_LEVEL=1, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /' $d/c.fds
  (cd $d && $MPI "$BIN" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize $m --regrid-every 1000000 --tag-phi --ls-fine --nofine 0 \
     --steps $N --outdir . --chid c > run.log 2>&1) || { echo "FAIL run m$m (see $d/run.log)"; exit 1; }
  grep -a "level 1 = " $d/run.log | head -1
done
python3 - "$W/m8" "$W/m32" <<'EOP'
import sys, csv, glob, os
a, b = sys.argv[1], sys.argv[2]; worst = 0.0; bad = 0; nd = 0
for f in sorted(glob.glob(a + "/c_*.csv")):
    n = os.path.basename(f)
    if "perf" in n or "steps" in n or "multi_level" in n: continue   # timings; per-box diagnostics
    g = os.path.join(b, n)
    if not os.path.exists(g): print("MISSING", n); bad += 1; continue
    ra, rb = list(csv.reader(open(f))), list(csv.reader(open(g)))
    if len(ra) != len(rb): print("ROWS", n, len(ra), len(rb)); bad += 1; continue
    m = 0.0
    for x, y in zip(ra[2:], rb[2:]):
        for p, q in zip(x, y):
            try: p, q = float(p), float(q)
            except ValueError: continue
            s = max(abs(p), abs(q))
            if s > 0: m = max(m, abs(p - q) / s)
    gated = n == "c_devc.csv"   # field probes; c_hrr integrals (Q_COND, Q_ENTH near 0, HRR 3e-9) also differ by box layout outside the ground stress: info only
    print(f"{n}: max rel diff {m:.3e}" + ("" if gated else " (info)"))
    if gated: worst = max(worst, m); nd += 1
ok = bad == 0 and nd > 0 and worst <= 1e-11
print(("PASS" if ok else "FAIL") + f" ground stress seam: worst rel diff {worst:.3e} (<= 1e-11)"); sys.exit(0 if ok else 1)
EOP
