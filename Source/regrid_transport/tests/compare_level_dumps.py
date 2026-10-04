#!/usr/bin/env python3
"""Compare two dumps of DriverModes (RHO and ZZ on the valid cells of one level).
usage: compare_level_dumps.py <prefixA> <levelA> <prefixB> <levelB> [tol]
Files <prefix>.L<level>.r<rank>.bin: int32 nboxes; per box 6 int32 (lo, hi) and the doubles RHO, ZZ_1..ZZ_n, i fastest. The number n is not stored: it is solved
from the file length. Prints the largest differences (relative to the largest value); exit 1 when above tol (default 1e-12)."""
import glob, struct, sys
import numpy as np


def parse(fn):
    b = open(fn, "rb").read()
    nb = struct.unpack_from("<i", b, 0)[0]
    if nb == 0:
        return []
    for ns in range(1, 9):
        o, ok, boxes = 4, True, []
        for _ in range(nb):
            if o + 24 > len(b):
                ok = False
                break
            lo = struct.unpack_from("<3i", b, o)
            hi = struct.unpack_from("<3i", b, o + 12)
            o += 24
            n = (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1)
            boxes.append((lo, hi, o, n))
            o += 8 * n * (1 + ns)
        if ok and o == len(b):
            return [(lo, hi, np.frombuffer(b, dtype="<f8", count=n * (1 + ns), offset=off).reshape((1 + ns, hi[2] - lo[2] + 1, hi[1] - lo[1] + 1, hi[0] - lo[0] + 1)))
                    for lo, hi, off, n in boxes]
    raise SystemExit("cannot parse " + fn)


def load(prefix, lev):
    """dict component -> {(i,j,k): value}"""
    out = {}
    for fn in sorted(glob.glob(f"{prefix}.L{lev}.r*.bin")):
        for lo, hi, arr in parse(fn):
            for c in range(arr.shape[0]):
                d = out.setdefault(c, {})
                for k in range(arr.shape[1]):
                    for j in range(arr.shape[2]):
                        for i in range(arr.shape[3]):
                            d[(lo[0] + i, lo[1] + j, lo[2] + k)] = arr[c, k, j, i]
    return out


a, b = load(sys.argv[1], int(sys.argv[2])), load(sys.argv[3], int(sys.argv[4]))
tol = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-12
if not a or not b or set(a[0]) != set(b[0]):
    print("RTE2E-CMP cell sets differ or missing dumps:", len(a.get(0, {})), "vs", len(b.get(0, {})))
    sys.exit(2)
names = ["RHO"] + [f"ZZ{c}" for c in range(1, len(a))]
bad = 0
for c in sorted(a):
    keys = list(a[c])
    va = np.array([a[c][k] for k in keys])
    vb = np.array([b[c][k] for k in keys])
    d = np.abs(va - vb)
    ref = max(1e-300, np.max(np.abs(va)))
    print(f"RTE2E-CMP {names[c]}: cells {len(keys)}, max abs diff {d.max():.3e}, relative to the largest value {d.max() / ref:.3e}, cells above tol {int((d > tol * ref).sum())}")
    bad += int(d.max() > tol * ref)
sys.exit(1 if bad else 0)
