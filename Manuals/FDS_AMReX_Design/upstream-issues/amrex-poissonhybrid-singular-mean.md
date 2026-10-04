# Draft upstream issue for AMReX: `FFT::PoissonHybrid` returns a wrong solution for a fully singular problem unless the right-hand side is already compatible

**State: draft for the project owner to review. Not sent, not filed, nothing posted anywhere.** The text below is written so it can be pasted into a GitHub issue on AMReX-Codes/amrex as it stands; the proposed fix could be attached as a pull request.

Tested against AMReX 26.09 at commit 99ddfda (`SpMatrix: 32-bit row offsets in the split blocks (#5968)`), CPU build (GNU 14.2, Open MPI 5.0.7, OpenMP on, `AMReX_FFT=ON`), install `(local AMReX install built against HYPRE 2.32)`. An earlier study of the same problem on a CUDA build (test-machine GPU) gives the same numbers; the fix below changes both the GPU and the CPU branch of the function.

---

## Title

`FFT::PoissonHybrid`: fully singular problem (periodic or Neumann in x and y, Neumann in z) silently solves a perturbed system when the right-hand side is not compatible (`FFT::Poisson` handles the same case)

## Summary

For a Poisson problem whose operator is singular (periodic or Neumann at both ends in x and in y, and Neumann at both ends in z, so the constant is in the null space), `FFT::PoissonHybrid::solve` pins the zero Fourier mode by doubling the last diagonal entry of the tridiagonal system in z. This gives the right answer only if the right-hand side already satisfies the discrete compatibility condition `sum_k dz_k f_k = 0` for the zero mode. If it does not, the solver does not fail, warn, or drop the incompatible part (as `FFT::Poisson` does, which skips the zero mode); it solves a perturbed system whose residual sits in the last z layer and is about `nz` times the right-hand-side mean. With a mean of 2e-2 (relative to max |f|) on a 32^3 grid the residual of the discrete equation is 0.75 of max |f| instead of 1e-14.

Affected: all boundary combinations with `is_singular` true in `PoissonHybrid::solve_z`, i.e. (periodic or Neumann/Neumann) in x and y, and `Boundary::even` at both z ends. Not affected: any problem with a Dirichlet side or a half-integer offset (`offset != 0`); those are non-singular.

## Minimal reproducer

One file, links against an installed AMReX (3D, FFT enabled). It builds a 32^3 grid, dx = dy = 1/32, periodic in x and y, Neumann (`Boundary::even`) at both z ends, with a uniform z spacing (`dz=u`) or a geometrically stretched one (`dz=geo`, dz_k = dx * 1.01^k). The right-hand side is a hash-noise plus smooth field with a non-zero mean ("raw"), and the same field with the dz-weighted mean removed ("compatible"). It solves with `FFT::Poisson` (uniform z only) and `FFT::PoissonHybrid`, then evaluates the discrete operator

    (A p)_{ijk} = (p_{i+1} - 2 p_i + p_{i-1})/dx^2 + (same in y) + [2/(dz_k (dz_k+dz_{k+1})) (p_{k+1}-p_k) - 2/(dz_k (dz_k+dz_{k-1})) (p_k-p_{k-1})],  zero flux at both z ends

and prints `max|A p - (f - <f>)| / max|f|`, where `<f>` is the dz-weighted mean (for a singular problem the usual convention is that the incompatible constant part of f is dropped, which is what `FFT::Poisson` does), and `max|A p - f| / max|f|` for the raw right-hand side.

Run: `mpirun -np 1 ./hybrid_singular_repro dz=u` and `... dz=geo`.

```cpp
// Minimal reproducer: amrex::FFT::PoissonHybrid on a fully singular problem (periodic x and y, Neumann z at both ends).
// The discrete operator is A p = (p_{i+1}-2p_i+p_{i-1})/dx^2 + (same in y) + z finite-volume term with cell widths dz_k:
//   (A p)_k = 2/(dz_k (dz_k+dz_{k+1})) (p_{k+1}-p_k) - 2/(dz_k (dz_k+dz_{k-1})) (p_k-p_{k-1}),  zero flux at both z ends.
// A p = f is solvable only if sum_k dz_k f_k = 0 (summed over x,y too). For a right-hand side that does not satisfy this, the usual
// convention (and what amrex::FFT::Poisson does) is to solve A p = f - <f>_dz, with <f>_dz the dz-weighted mean.
// Usage: hybrid_singular_repro [dz=u|geo]    prints one line per (solver, rhs) with the residuals relative to max|f|.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT_Poisson.H>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
using namespace amrex;

static double hashf (int i, int j, int k) {
    unsigned h = (unsigned)i*73856093u ^ (unsigned)j*19349663u ^ (unsigned)k*83492791u;
    h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
    return double(h & 0xFFFFFFu)/16777216.0 - 0.5;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        std::string dzmode = "u"; { ParmParse pp; pp.query("dz", dzmode); }
        const int n = 32; const double dx = 1.0/n;
        std::vector<double> dz(n); double zl = 0.0;
        for (int k = 0; k < n; ++k) { dz[k] = (dzmode == "geo") ? dx*std::pow(1.01, k) : dx; zl += dz[k]; }
        Box dom(IntVect(0,0,0), IntVect(n-1,n-1,n-1));
        RealBox rb({0.,0.,0.}, {1.0, 1.0, zl});
        Geometry geom(dom, rb, 0, Array<int,3>{1,1,0});
        BoxArray ba(dom); DistributionMapping dm(ba);
        Array<std::pair<FFT::Boundary,FFT::Boundary>,3> bc{
            std::make_pair(FFT::Boundary::periodic, FFT::Boundary::periodic),
            std::make_pair(FFT::Boundary::periodic, FFT::Boundary::periodic),
            std::make_pair(FFT::Boundary::even, FFT::Boundary::even)};   // even = Neumann

        for (int compatible = 0; compatible < 2; ++compatible)
        {
            MultiFab rhs(ba, dm, 1, 0), soln(ba, dm, 1, 0);
            std::vector<double> f(size_t(n)*n*n);
            auto at = [&](int i, int j, int k) -> size_t { return size_t(i) + size_t(n)*(j + size_t(n)*k); };
            double sw = 0.0, swf = 0.0;
            for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
                f[at(i,j,k)] = hashf(i,j,k) + 0.3*std::sin(2.0*M_PI*(i+0.5)/n)*std::cos(2.0*M_PI*(k+0.5)/n) + 0.02;   // nonzero mean
                sw += dz[k]; swf += dz[k]*f[at(i,j,k)];
            }
            const double fmean = swf/sw;
            if (compatible) for (auto& v : f) { v -= fmean; }
            double fmax = 0.0; for (double v : f) fmax = std::max(fmax, std::abs(v));
            {
                auto const& a = rhs.array(0);
                for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) a(i,j,k) = f[at(i,j,k)];
            }
            for (int solver = 0; solver < 2; ++solver)
            {
                if (solver == 0 && dzmode != "u") continue;      // FFT::Poisson has no stretched-z option
                soln.setVal(0.0);
                if (solver == 0) { FFT::Poisson<MultiFab> s(geom, bc); s.solve(soln, rhs); }
                else             { FFT::PoissonHybrid<MultiFab> s(geom, bc); s.solve(soln, rhs, Vector<double>(dz.begin(), dz.end())); }
                auto const& p = soln.const_array(0);
                double r_comp = 0.0, r_raw = 0.0;
                for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
                    const int ip = (i+1)%n, im = (i+n-1)%n, jp = (j+1)%n, jm = (j+n-1)%n;
                    double Ap = (p(ip,j,k) - 2*p(i,j,k) + p(im,j,k))/(dx*dx) + (p(i,jp,k) - 2*p(i,j,k) + p(i,jm,k))/(dx*dx);
                    const double c = (k < n-1) ? 2.0/(dz[k]*(dz[k]+dz[k+1]))*(p(i,j,k+1)-p(i,j,k)) : 0.0;
                    const double a = (k > 0)   ? 2.0/(dz[k]*(dz[k]+dz[k-1]))*(p(i,j,k)-p(i,j,k-1)) : 0.0;
                    Ap += c - a;
                    const double fv = f[at(i,j,k)];
                    r_raw  = std::max(r_raw,  std::abs(Ap - fv));
                    r_comp = std::max(r_comp, std::abs(Ap - (fv - (compatible ? 0.0 : fmean))));
                }
                std::printf("SINGULAR dz=%-3s rhs=%-10s solver=%-12s dz-weighted mean(f)=% .3e  max|A p - (f - mean)|/max|f| = %.3e   max|A p - f|/max|f| = %.3e\n",
                            dzmode.c_str(), compatible ? "compatible" : "raw", solver == 0 ? "Poisson" : "PoissonHybrid",
                            compatible ? 0.0 : fmean, r_comp/fmax, r_raw/fmax);
            }
        }
    }
    amrex::Finalize();
    return 0;
}
```

## Expected and observed

Residuals relative to max |f|, 32^3, one rank, one thread, Release build.

Shipped header (`(local AMReX install built against HYPRE 2.32)/include/AMReX_FFT_Poisson.H`):

| z spacing | right-hand side | solver | dz-weighted mean of f | max abs(A p - (f - mean)) / max abs(f) | max abs(A p - f) / max abs(f) |
|---|---|---|---|---|---|
| uniform | raw | `Poisson` | 1.944e-02 | 8.9e-14 | 2.4e-02 (the dropped mean, as intended) |
| uniform | raw | **`PoissonHybrid`** | 1.944e-02 | **7.5e-01** | 7.7e-01 |
| uniform | compatible | `Poisson` | 0 | 1.6e-14 | 1.6e-14 |
| uniform | compatible | `PoissonHybrid` | 0 | 3.3e-14 | 3.3e-14 |
| geometric 1.01 | raw | **`PoissonHybrid`** | 1.946e-02 | **6.4e-01** | 6.7e-01 |
| geometric 1.01 | compatible | `PoissonHybrid` | 0 | 3.2e-14 | 3.2e-14 |

With the proposed change (`PoissonHybrid` rows only; `Poisson` is untouched and identical):

| z spacing | right-hand side | solver | max abs(A p - (f - mean)) / max abs(f) |
|---|---|---|---|
| uniform | raw | `PoissonHybrid` | 1.4e-14 |
| uniform | compatible | `PoissonHybrid` | 1.7e-14 |
| geometric 1.01 | raw | `PoissonHybrid` | 1.7e-14 |
| geometric 1.01 | compatible | `PoissonHybrid` | 1.7e-14 |

A wider study with an independent reference (separable eigen-decomposition in numpy, full discrete operator including the stretched z term), 158 cases: grids 32^3 and 128^3, boundary types DD/DN/ND/NN/PP in each direction, `Poisson` and `PoissonHybrid`, z spacing uniform / symmetric stretch / geometric stretch, right-hand sides raw (non-zero mean), dz-weighted-mean-removed, and "raw plus 1": as shipped, 24 of 158 fail on the CPU build (and the same 24 on the GPU build) and all 24 are singular `PoissonHybrid` cases with a non-compatible right-hand side (PP,PP,NN / NN,NN,NN / PP,NN,NN at both sizes, e.g. 32^3 PP,PP,NN raw: residual 1.7e-2, relative error against the reference 0.17; "raw plus 1": residual 31, error 303); with the change, 158 of 158 pass on both builds (residual 1.5e-15, error 1e-14 level). Compatible right-hand sides and every non-singular boundary combination give identical results with and without the change.

## Root cause

`Src/FFT/AMReX_FFT_Poisson.H`, `PoissonHybrid<MF>::solve_z`, branch for non-trivial z coefficients (the `else` of the `Tri_Zero` test):

* `is_singular` is set at lines 701-702: `(offset[0] == 0) && (offset[1] == 0) && zlo_neumann && zhi_neumann`.
* GPU branch, lines 748-752, and CPU branch, lines 811-815 (same text with `bd[k]` instead of `bd(i,j,k)`):

  ```cpp
  } else if (k == nz-1) {
      ald(i,j,k) = tria(i,j,k);
      cud(i,j,k) = T(0.);
      if (zhi_neumann) {
          bd(i,j,k) = k2 - ald(i,j,k);
          if (i == 0 && j == 0 && is_singular) {
              bd(i,j,k) *= T(2.0);           // <- pins the zero mode
          }
  ```

For the zero Fourier mode (k2 = 0) the tridiagonal matrix T (rows `c_k p_{k+1} - (a_k + c_k) p_k + a_k p_{k-1}`, with `a_0 = c_{nz-1} = 0`) is singular. With weights `w_k` that satisfy `w_{k+1} a_{k+1} = w_k c_k` (for the finite-volume form, `w_k = dz_k`), `w^T T = 0`, so `T p = f` is solvable only if `w^T f = 0`. Doubling the last diagonal replaces T by `T' = T - a_{nz-1} e e^T` (e = last unit vector). Multiplying `T' p = f` by `w^T` gives `p_{nz-1} = - (w^T f) / (w_{nz-1} a_{nz-1})`:

* if `w^T f = 0`, then `p_{nz-1} = 0`, `T p = f` holds exactly, and the doubled entry is just a pin (this is the case the code was written for, and why every test with a compatible right-hand side passes);
* if `w^T f != 0`, then `p_{nz-1} != 0` and `T p = f + a_{nz-1} p_{nz-1} e`: the whole zero-mode mean of f ends up as a residual in the top cell, of size `(w^T f)/w_{nz-1}`, i.e. about `nz` times the mean for uniform z (0.019 * 32 = 0.62 here, observed 0.75 relative to max |f|).

The zero mode of `FFT::Poisson` is handled by skipping the division when `k2 == 0` (`PoissonHybrid`'s own 2-D branch does the same, lines 692-694; `Poisson::solve`, lines 353-355), which is the implicit removal of the mean; `PoissonHybrid` has no equivalent for the singular mode, and the pin does not provide one. The unweighted mean is not enough for a stretched z: the compatibility weights are the `dz_k`, which is why the sum is built from the matrix coefficients (`w_{k+1}/w_k = c_k/a_{k+1}`) and not from a plain mean. The existing test (`Tests/FFT/Poisson/main.cpp`, `make_rhs`) subtracts the plain mean of the right-hand side whenever there is no Dirichlet side, and calls `PoissonHybrid` only with uniform `dz`, where the plain and weighted means coincide, so it cannot see the problem.

## Proposed fix

Remove the weighted mean of the zero-mode column of the spectral right-hand side before the tridiagonal solve, only when `is_singular`, in both branches (hunks below; `patch -p1` against the AMReX source tree, applies to 99ddfda). The weights are built from `tria` and `tric`, so it works for any spacing the solver accepts; the cost is one extra pass over `nz` entries for the single (i, j) = (0, 0) column. Nothing changes when `is_singular` is false or the right-hand side is already compatible (the mean found is round-off and the subtraction changes the column by round-off only). The doubled last diagonal stays as the pin.

```diff
--- a/Src/FFT/AMReX_FFT_Poisson.H
+++ b/Src/FFT/AMReX_FFT_Poisson.H
@@ -758,6 +758,21 @@
                     }
                 }

+                if (i == 0 && j == 0 && is_singular) {
+                    // Fully singular problem: the zero mode of the right-hand side must satisfy the discrete
+                    // compatibility condition sum_k w_k f_k = 0, w_{k+1}/w_k = c_k/a_{k+1} (w_k = dz_k).
+                    // The doubled last diagonal above pins the solution but silently solves a perturbed
+                    // problem unless that condition holds, so remove the weighted mean here.
+                    T wk = T(1), sw = T(0);
+                    auto swf = spectral(i,j,0) * T(0);
+                    for (int k = 0; k < nz; ++k) {
+                        sw += wk; swf += spectral(i,j,k) * wk;
+                        if (k < nz-1) { wk *= cud(i,j,k)/ald(i,j,k+1); }
+                    }
+                    auto fmean = swf * (T(1)/sw);
+                    for (int k = 0; k < nz; ++k) { spectral(i,j,k) -= fmean; }
+                }
+
                 scratch(i,j,0) = cud(i,j,0)/bd(i,j,0);
                 spectral(i,j,0) = spectral(i,j,0)/bd(i,j,0);

@@ -821,6 +836,21 @@
                     }
                 }

+                if (i == 0 && j == 0 && is_singular) {
+                    // Fully singular problem: the zero mode of the right-hand side must satisfy the discrete
+                    // compatibility condition sum_k w_k f_k = 0, w_{k+1}/w_k = c_k/a_{k+1} (w_k = dz_k).
+                    // The doubled last diagonal above pins the solution but silently solves a perturbed
+                    // problem unless that condition holds, so remove the weighted mean here.
+                    T wk = T(1), sw = T(0);
+                    auto swf = spectral(i,j,0) * T(0);
+                    for (int k = 0; k < nz; ++k) {
+                        sw += wk; swf += spectral(i,j,k) * wk;
+                        if (k < nz-1) { wk *= cud[k]/ald[k+1]; }
+                    }
+                    auto fmean = swf * (T(1)/sw);
+                    for (int k = 0; k < nz; ++k) { spectral(i,j,k) -= fmean; }
+                }
+
                 scratch[0] = cud[0]/bd[0];
                 spectral(i,j,0) = spectral(i,j,0)/bd[0];

```

Alternative if the maintainers prefer a different place: zero the (0, 0) mode of the weighted-mean part in the caller before the transform, or drop the zero mode as `FFT::Poisson` does (set the solution of that column to its pinned particular solution) and document the convention. The above keeps the existing pin and the existing answer for compatible input.

## Test suggestion

Extend `Tests/FFT/Poisson/main.cpp` (the `PoissonHybrid` loop):

1. Add a variant that does not apply the plain-mean shift of `make_rhs` for the hybrid solver when the case is singular (all of x, y non-Dirichlet and z `even/even`), and check the residual against `rhs - <rhs>_w` with `check_convergence` modified to subtract the dz-weighted mean of the right-hand side for singular cases (`rnorm < eps * bnorm`, as the other cases). Without the fix this case fails at about 0.5 of `bnorm` (this reproducer); with it it passes at the usual 2e-10.
2. Add a stretched `dz` (for example `dz_k = dz_0 * 1.01^k`, with `geom` z extent set to the sum) next to the existing uniform `dz`, for all `PoissonHybrid` cases, with `check_convergence` using the matching variable-spacing z term (the finite-volume form given above). The existing uniform `dz` cannot distinguish the plain and the weighted mean.
3. Optionally assert that a compatible right-hand side gives the same answer to round-off with and without the weighted-mean step (guards against changing the existing behaviour).

## Notes for the reviewer (not part of the issue text)

* The reproducer and numbers above were produced with the shipped header and with the patched header shadowing it through `-I`; the installed AMReX and the AMReX source tree were not modified.
* Earlier measurements of the 158-case study are in `prototypes/s4_cuda_mass/stage1/fft/results/` (`sing_compare_host.log`, `sing_compare_host_fix.log`, `sing_compare_gpu.log`, `sing_compare_gpu_fix.log`, `sing_table.md`); the fix as a file is `prototypes/s4_cuda_mass/stage1/fft/amrex_fix/hybrid_singular_mean.patch`.
* Why it matters for the FDS-to-AMReX work: the hybrid solver is a candidate for a stretched-z periodic or closed box, where the right-hand side of the pressure equation is compatible only up to round-off. The pressure backend does not rely on `PoissonHybrid` for singular components (it uses `FFT::Poisson` or MLMG and removes the mean itself in its common layer), so nothing here is blocking; the issue is about correctness of the library.
* If the fix is filed as a pull request, the project owner may want to check whether the documentation (`Docs/sphinx_documentation/source/FFT.rst`, the `PoissonHybrid` paragraph) should state the convention (the incompatible constant part of the right-hand side is dropped for a singular problem, as for `FFT::Poisson`).
