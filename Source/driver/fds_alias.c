/* fds_alias.c: the ONLY place where an unallocated Fortran ALLOCATABLE array (a MESH_TYPE component such as MESHES(NM)%RHO) is made to
 * describe externally owned memory (an AMReX FAB) with FDS lower bounds, by writing the C descriptor directly (decision 1, "W2 alias").
 *
 * Kernel-facing rules (M2a): (a) passive scalars: a 4-D array (ZZ, ZZS) is aliased with rank 4 and the scalar count as the last extent;
 * nothing here depends on N_TOTAL_SCALARS. (b) only uniform Cartesian metrics are used; this file handles no metric data.
 *
 * Non-standard: F2018 18.5.5 allows CFI_allocate/CFI_deallocate only. It works with gfortran 14 because the CFI descriptor is converted
 * back into the gfortran descriptor on return (prototype P1, p1_alias.c). The array must never be DEALLOCATEd or reallocated by Fortran
 * while aliased; fds_alias_release() nulls it first.
 *
 * Extension over P1: element strides are passed explicitly, so a window of a larger FAB (RHO/RHOS with ng=3 seen as FDS ng=2) can be
 * described. CAUTION (tests/test_driver.cpp, section "strided alias"): the Fortran standard makes an allocated allocatable array
 * contiguous, so a compiler may pass a strided alias to an explicit-shape dummy without copy-in/copy-out. Use strides other than the
 * contiguous ones only after that test passes for the routines concerned. */
#include <ISO_Fortran_binding.h>
#include <stddef.h>

void fds_alias_alloc(CFI_cdesc_t *a, void *base, const int *lb, const int *ext, const long *stride)
{
    for (int i = 0; i < a->rank; ++i) {
        a->dim[i].lower_bound = lb[i];
        a->dim[i].extent      = ext[i];
        a->dim[i].sm          = (CFI_index_t)stride[i] * (CFI_index_t)a->elem_len;
    }
    a->base_addr = base;
}

void fds_alias_release(CFI_cdesc_t *a)
{
    a->base_addr = NULL;
}
