# Where the fine-level excess error sits for mixed Neumann/Dirichlet directions on a hierarchy

Note for the Pressure Lead. Measured with `pb_harness mode=err_map` (harness/composite_modes.cpp, `run_err_map`); nothing in the
library was changed for it. The numbers below are from the MLMG composite backend; the HYPRE backend gives the same solution to
round-off (see hypre-notes.md), so the picture is the same for both.

## Question

On a two-level hierarchy with a mixed pair of faces in one direction (N at the low face and D at the high face, "ND", or "DN"),
the composite error on the fine level exceeds the error of the uniform coarse solve at one resolution (compconv therefore checks
the convergence order instead for mixed faces). Is the excess located at the coarse/fine (C/F) faces or at the domain faces?

## Experiment

Manufactured solution of `composite_modes.cpp` (half-integer cosine/sine in the mixed direction, so it satisfies N at one face and D
at the other, cosine/cosine in the others, wave numbers (2,1,2)). n coarse cells per direction, ratio r, one refined patch.
Every uncovered cell of every level is classed by its nearest feature, taking the distance in cells of its own level: `CF<m>` next to
a C/F face (m = 0 adjacent, on either side of the interface), `domN<m>`/`domD<m>` next to a closed/open domain face, and so on;
`<m>` is capped at 3, so `domN3`/`domD3` are in practice "bulk, three or more cells from any feature". For the fine level three
errors are reported: `comp` (composite minus exact), `u` (uniform grid at the fine resolution, FFT, minus exact) and `d` (composite
minus uniform, the part of the composite error that is not discretisation error of the fine grid). All values are rms over the class.

```
mpirun -np 2 pb_harness mode=err_map n=32 nlev=2 ratio=2 mgs=16 bcfaces=ND,NN,NN layout=0     # patch in the middle
mpirun -np 2 pb_harness mode=err_map n=32 nlev=2 ratio=2 mgs=16 bcfaces=DN,NN,NN layout=1     # patch at the low corner
mpirun -np 2 pb_harness mode=err_map n=32 nlev=2 ratio=2 mgs=16 bcfaces=ND,NN,NN full=1       # fine level over the whole domain
```

## Results (rms, fine level unless said otherwise)

| case | all fine: comp / u / d | fine at C/F (m=0): comp / u / d | fine at domain faces, m=0: comp / u / d | coarse level: C/F m=0, domN m=0, domD m=0 |
|---|---|---|---|---|
| n=32 ND, patch in the middle | 3.5e-4 / 2.1e-4 / 1.9e-4 | 4.6e-4 / 2.4e-4 / 2.9e-4 | (patch does not reach a domain face) | 1.0e-3, 2.0e-3, 2.5e-4 |
| n=32 ND, low corner (touches the N face) | 3.3e-4 / 3.8e-4 / 4.1e-4 | 4.9e-4 / 3.8e-4 / 6.7e-4 | N face: 4.9e-4 / 5.1e-4 / 4.2e-4 | 1.4e-3, 2.0e-3, 2.5e-4 |
| n=32 DN, low corner (touches the D face) | 2.9e-4 / 3.3e-4 / 1.6e-4 | 3.6e-4 / 3.3e-4 / 2.8e-4 | D face: 2.7e-5 / 3.1e-5 / 1.1e-5; N face: 4.5e-4 / 4.7e-4 / 1.9e-4 | 1.8e-3, 2.0e-3, 2.5e-4 |
| n=32 ND and DN, fine level over the whole domain (no C/F face) | 3.5e-4 / 3.5e-4 / 3e-15 | none | N and D faces: d between 6e-17 and 1e-14 | no uncovered coarse cell |
| n=64 ND, patch in the middle | 9.0e-5 / 5.3e-5 / 4.9e-5 | 1.2e-4 / 5.9e-5 / 8.0e-5 | none | 2.5e-4, 5.0e-4, 3.1e-5 |
| n=32 NN, patch in the middle (control) | 2.8e-4 / 1.6e-4 / 1.5e-4 | 3.5e-4 / 1.5e-4 / 2.3e-4 | none | 6.2e-4, 1.4e-3 |
| n=32 DD, patch in the middle (control) | 2.3e-4 / 1.6e-4 / 1.0e-4 | 3.2e-4 / 1.9e-4 / 1.7e-4 | none | 8.9e-4, 1.5e-3, 1.5e-4 |
| ratio 4, n=16 ND, patch in the middle | 9.9e-4 / 2.1e-4 / 8.6e-4 | 1.5e-3 / 2.4e-4 / 1.3e-3 | none | 4.3e-3, 7.8e-3, 2.0e-3 |

## Reading

1. The excess is not produced at the domain faces. With the fine level over the whole domain (no C/F face) the composite and the
   uniform solution agree to round-off (`d` about 1e-15) for both ND and DN, so the closed and the open faces are discretised
   identically on the fine level. With the patch touching a domain face, `d` at that face (4.2e-4 for N, 1.1e-5 for D) is of the same
   size as `d` at the C/F faces of the same patch; it is the interface error arriving there, not a boundary defect.
2. The excess is a smooth error field that starts at the C/F faces and is inherited from the coarse level. It is largest in the
   first fine cells at the interface (d = 2.9e-4 at m=0 against 1.9e-4 at m=2 and 1.2e-4 in the bulk for n=32 ND) and decays into
   the patch. Its size follows the coarse-level error, which is the dominant error of the whole hierarchy: the coarse error is
   largest next to the closed (N) domain face (2.0e-3 against 2.5e-4 next to the D face; the manufactured solution has non-zero
   normal curvature there), then next to the C/F faces.
3. The mixed direction is not special. The pure NN and DD controls show the same pattern with the same ratio of `d` to the uniform
   error (NN 0.94, DD 0.64, ND 0.9 for all fine cells, n=32). What differs for mixed faces is only that the half-integer
   manufactured solution does not vanish at the patch edges, so the fine-level error at one resolution is not guaranteed to be below
   the coarse uniform error, which is why compconv checks the order.
4. Convergence is second order for every class: ND, patch in the middle, n=32 to 64: all fine cells 3.5e-4 to 9.0e-5 (3.9x), `d` 1.9e-4 to
   4.9e-5 (3.9x); C/F m=0: 4.6e-4 to 1.2e-4. Ratio 4 has a larger excess (d/u about 4 against 0.9) because the interface
   interpolation then spans four fine cells per coarse cell on a coarse grid that is still 16 cells wide.

Conclusion for the Lead: the fine-level excess for mixed N/D directions sits at the C/F faces (and decays over a few fine cells into
the patch), with its magnitude set by the coarse-level truncation error near the closed domain face. No domain-face treatment of
the fine level needs to change. Limits of the experiment: one manufactured solution, one patch shape, ND/DN in x only, MLMG backend.
