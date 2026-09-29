#!/usr/bin/env python3
"""Cross-check the registered field table against the field inventory (mesh_fields.csv): declared/allocated bounds must match.
usage: check_inventory.py <driver_unit_tests binary> <mesh_fields.csv>
Passive scalars: ZZ/ZZS are checked with the fourth extent N_TOTAL_SCALARS. Uniform Cartesian metrics only: no metric arrays are registered."""
import csv, re, subprocess, sys
exe, csvf = sys.argv[1], sys.argv[2]
rows = {r['field']: r for r in csv.DictReader(open(csvf)) if r['parent_type'] == 'MESH_TYPE'}
out = subprocess.run([exe, '--dump-table'], capture_output=True, text=True).stdout
bad = n = 0
for l in out.splitlines():
    if not l.startswith('TABLE'):
        continue
    p = l.split(None, 6)
    name, bounds = p[1], p[6]
    r = rows.get(name)
    m = re.search(r'allocated as (\(.*\))', r['declared_bounds']) if r else None
    n += 1
    if not m or m.group(1).replace(' ', '') != bounds.replace(' ', ''):
        bad += 1
        print('MISMATCH', name, bounds, '|', r['declared_bounds'] if r else 'not in inventory')
print(('PASS' if bad == 0 else 'FAIL'), 'inventory cross-check: %d fields, %d mismatches' % (n, bad))
sys.exit(1 if bad else 0)
