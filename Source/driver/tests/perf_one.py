#!/usr/bin/env python3
"""One timed run for NFR-030 / NFR-031 (measured, not required): usage perf_one.py <workdir> <np> <cmd...>
Runs <cmd> under mpirun in <workdir>, prints one CSV line: wall_s, max_rss_kb_of_one_rank (RUSAGE_CHILDREN maxrss = largest descendant), exit code."""
import os, resource, subprocess, sys, time
wd, np_ = sys.argv[1], sys.argv[2]
cmd = ['mpirun', '--bind-to', 'none', '--oversubscribe', '-np', np_] + sys.argv[3:]
t0 = time.time()
with open(os.path.join(wd, 'stdout.txt'), 'w') as so, open(os.path.join(wd, 'stderr.txt'), 'w') as se:
    rc = subprocess.call(cmd, cwd=wd, stdout=so, stderr=se)
wall = time.time() - t0
print('%.3f,%d,%d' % (wall, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, rc))
