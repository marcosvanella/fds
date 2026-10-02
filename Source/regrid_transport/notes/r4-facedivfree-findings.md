# R4: what the installed AMReX FaceDivFree interpolator does (probe results)

Test: `tests/probe_facedivfree.cpp` (ctest `regrid_transport_facedivfree_probe`, `regrid_transport_facedivfree_ratio1_aborts`). It calls the interpolator directly on a coarse face field and one fine patch of 6x6x6 coarse cells. The coarse field is either exactly divergence-free (discrete curl of a random potential) or random with a nonzero divergence. "Divergence of a child" means the divergence of one fine cell computed from its six fine faces. AMReX 26.09 as installed, built for 3-D only (a 2-D run of FDS is a 3-D run with one cell in the hidden direction).

| Case | Result |
|---|---|
| Ratio 2 in all directions, 3-D | Works. Divergence-free coarse gives divergence-free fine (max divergence 6e-16 of the scale u/dx). Every child has the divergence of its parent (difference 7e-16) also for a coarse field with divergence. |
| Ratio 4 in all directions | Works, same properties (difference 3e-14). |
| Mixed ratio (2,4,2) | Works, same properties (difference 8e-15). |
| Ratio 1 in a direction (single-cell direction of a 2-D run) | **Not supported.** AMReX aborts: "Only refinement ratio of 2 or 4 is supported" (AMReX_Interpolater.cpp:1378, always-on assertion, not catchable). It also needs a ratio of 2 or 4 in every direction, so a one-cell direction cannot be refined by 2 either without changing FDS 2-D semantics. |
| Faces on the patch boundary | AMReX does not copy the coarse value there: it puts the coarse value plus the transverse slope, so the fine values differ from the coarse face value by a sizeable fraction (about a quarter of the largest face value in the random-field case; their mean over the coarse face equals the coarse value to 1e-15). To get the ruled behaviour (interface faces take the coarse value) pre-fill those fine faces with the coarse value and pass a mask that is 0 on them: AMReX then skips them and computes the interior faces from the pre-filled values. Tested: interface faces equal the coarse value exactly, children keep the parent divergence to round-off, face means unchanged. |
| Fallback `prolong_faces_normal_linear` (`FaceTransfer.H`) | Each fine face takes the coarse value of the coarse cell that contains it in the transverse directions, and is linear between the two coarse faces in the normal direction. Works for any ratio including 1. Every child has exactly the divergence of its parent (difference 4e-16 to 1e-15, ratios 2, 4, (2,4,2), hidden direction (2,1,2) and (4,1,4)). Interface faces equal the coarse value, face means are unchanged. The price is a velocity that is constant across a coarse face in the transverse directions (less smooth than FaceDivFree). |

Consequence for the design (proposal, Architect to confirm):
1. 3-D and any run with no single-cell direction: `FaceDivFree` with the interface faces pre-filled and masked.
2. 2-D (single-cell direction, ratio 1): `prolong_faces_normal_linear`. It is not "injection plus projection": the children keep the parent divergence exactly, so no extra projection is needed for the face transfer itself. What remains is the variation of the target divergence D inside a parent cell; the number to report after each regrid is max|div u - D| on the new levels, as ruled.
3. The post-regrid projection stays off by default in both cases.

Not verified: behaviour with non-uniform cell sizes (FDS levels are uniform), with a GPU build (`FaceDivFree` needs the device for ratio other than 2 on GPU), and the effect on the later pressure solve (Phase 4).
