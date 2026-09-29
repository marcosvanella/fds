# frozen/

M1 uses synthetic frozen cases generated on the fly: `pb_harness mode=gen n_cell="..." bc=... out=<prefix>`
writes the RHS (`<prefix>_rhs.bin`) and the FFT reference H (`<prefix>_H.bin`) as raw float64, x fastest, and
the `frozen` CTest solves the stored RHS with MLMG and FFT and compares against the stored H. The RHS is the P2
smooth non-separable field (`rhs_func` in `harness/main.cpp`), made compatible once with the common-layer mean
removal. The committed store (manifest, small stored cases, reviewed generator script) is milestone M3;
Dirichlet and FDS-derived cases are M2 and M4.
