#!/bin/bash
# Smoke check of the WP1b edge-table dump (notes/wp1b-edge-tables.md): FDSTL_STAGES writes scr_c_vfluxE_{CE_EDGE_INDEX,ED_OMEGA,ED_TAU}_b0.bin
# with the flat layout of section 7.13b. usage: run_edge_tables_dump_check.sh <build-dir> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/edge_dump}; rm -rf "$W"; mkdir -p "$W"
cp "$HERE/cases/ns2d_16_l0.fds" "$W/"; cd "$W"
FDSTL_STAGES=1 FDSTL_STAGEG=1 timeout 120 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" ns2d_16_l0.fds --run --steps 1 --outdir . --chid edge > run.log 2>&1 || { echo "FAIL run"; tail -20 run.log; exit 1; }
python3 - <<'PY' || exit 1
import os, struct, sys
need = ["CE_EDGE_INDEX", "ED_OMEGA", "ED_TAU", "WORK1", "DS", "FVX", "FVY", "FVZ"]
for n in need:
    fn = f"scr_c_vfluxE_{n}_b0.bin"
    if not os.path.isfile(fn): print("FAIL missing", fn); sys.exit(1)
def hdr(fn):
    with open(fn, "rb") as f: return struct.unpack("<16i", f.read(64))
# CE_EDGE_INDEX: rank 4, lb=(0,0,0,1), ub=(IBP1,JBP1,KBP1,12) with IBP1=IBAR+1=17 for 16 cells
h = hdr("scr_c_vfluxE_CE_EDGE_INDEX_b0.bin")
assert h[0] == 4 and h[1:5] == (0, 0, 0, 1) and h[5] >= 2 and h[6] >= 2 and h[7] >= 2 and h[8] == 12, h  # 0:IBP1,0:JBP1,0:KBP1,1:12
# ED_OMEGA / ED_TAU: rank 2, lb=(0,-2), ub=(NEDGE,2)
for n in ("ED_OMEGA", "ED_TAU"):
    h = hdr(f"scr_c_vfluxE_{n}_b0.bin")
    assert h[0] == 2 and h[1:3] == (0, -2) and h[6] == 2 and h[5] >= 0, (n, h)
    nedge = h[5]
print(f"PASS edge dump: CE_EDGE_INDEX (0:{h[5]},0:{h[6]},0:{h[7]},12), ED_OMEGA/ED_TAU (0:{nedge},-2:2), WORK/DS/FV* present")
# companion raw inputs the harness needs
for n in ("US", "VS", "WS", "MU", "RHOS"):
    assert os.path.isfile(f"raw_c_prevflux_{n}_b0.bin"), n
print("PASS raw_c_prevflux companion fields present")
PY
