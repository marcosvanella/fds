// TagOps.cpp: see TagOps.H.
#include "TagOps.H"

#include <AMReX_ParallelDescriptor.H>

// rt_tag_kernels.F90 writes the literal 2 (RT_TAG_SET): AMReX's buffer() only grows from tags equal to TagBox::SET
static_assert(amrex::TagBox::SET == 2 && amrex::TagBox::CLEAR == 0, "TagBox values assumed by rt_tag_kernels.F90");

namespace {
struct Lay {
    int lo[3], hi[3];   // bounds of the whole array (valid + ghosts)
};
template <class A>
Lay bounds(const A& a)
{
    return Lay{{a.begin[0], a.begin[1], a.begin[2]}, {a.end[0] - 1, a.end[1] - 1, a.end[2] - 1}};
}
}  // namespace

extern "C" {
void rt_tag_cells(const int* lo, const int* hi, const double* q, const int* qlo, const int* qhi, const double* den, const int* use_den, const int* cov,
                  const int* clo, const int* chi, const int* use_cov, signed char* tag, const int* tlo, const int* thi, const int* mode, double thr,
                  double base, double keepfac, const int* rel, const int* dirs);
void rt_tag_box(const int* lo, const int* hi, signed char* tag, const int* tlo, const int* thi, const int* blo, const int* bhi);
void rt_tag_copy_box(const int* lo, const int* hi, const signed char* src, const int* slo, const int* shi, signed char* dst, const int* dlo, const int* dhi,
                     const int* blo, const int* bhi);
void rt_tag_count(const int* lo, const int* hi, const signed char* tag, const int* tlo, const int* thi, long* n);
}

namespace fdsrt {

namespace {
signed char* tag_ptr(amrex::TagBoxArray& tags, const amrex::MFIter& mfi, Lay& l)
{
    auto a = tags.array(mfi);
    l = bounds(a);
    return reinterpret_cast<signed char*>(a.p);
}
// Same grids and same box-to-rank map (equal, not necessarily the same reference-counted object: a field of the driver's level-0 registry is built from its own BoxArray).
bool same_layout(const amrex::FabArrayBase& x, const amrex::FabArrayBase& y)
{
    return x.DistributionMap() == y.DistributionMap() && x.boxArray().CellEqual(y.boxArray());
}
}  // namespace

void tag_cells(amrex::TagBoxArray& tags, const amrex::MultiFab& q, const amrex::MultiFab* den, int den_comp, const amrex::iMultiFab* covered,
               const TagCriterion& c)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(same_layout(tags, q), "tag_cells: field and tags must share BoxArray and DistributionMapping");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(c.keepfac > 0.0 && c.keepfac <= 1.0, "tag_cells: TAG_KEEP must be in (0,1]");
    for (int d = 0; d < 3; ++d)
        if (c.dirs[d] && c.mode == TagMode::Diff)
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(q.nGrowVect()[d] >= 1, "tag_cells: the field needs one ghost layer in every active direction");
    if (den) AMREX_ALWAYS_ASSERT_WITH_MESSAGE(same_layout(tags, *den) && den->nGrowVect().allGE(q.nGrowVect()), "tag_cells: den layout");
    if (covered) AMREX_ALWAYS_ASSERT_WITH_MESSAGE(same_layout(tags, *covered), "tag_cells: mask layout");
    const int mode = c.mode == TagMode::Above ? 0 : 1, rel = c.relative ? 1 : 0, use_den = den ? 1 : 0, use_cov = covered ? 1 : 0;
    const int dirs[3] = {c.dirs[0], c.dirs[1], c.dirs[2]};
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const int lo[3] = {vb.smallEnd(0), vb.smallEnd(1), vb.smallEnd(2)}, hi[3] = {vb.bigEnd(0), vb.bigEnd(1), vb.bigEnd(2)};
        Lay tl;
        signed char* tp = tag_ptr(tags, mfi, tl);
        auto qa = q.const_array(mfi, c.comp);
        const Lay ql = bounds(qa);
        const double* dp = qa.p;
        if (den) {
            auto da = den->const_array(mfi, den_comp);
            const Lay dl = bounds(da);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dl.lo[0] == ql.lo[0] && dl.hi[0] == ql.hi[0] && dl.lo[1] == ql.lo[1] && dl.hi[1] == ql.hi[1] &&
                                                 dl.lo[2] == ql.lo[2] && dl.hi[2] == ql.hi[2],
                                             "tag_cells: den needs the same ghost width as q");
            dp = da.p;
        }
        Lay cl = tl;
        const int* cp = nullptr;
        if (covered) {
            auto ca = covered->const_array(mfi);
            cl = bounds(ca);
            cp = ca.p;
        } else {
            static const int dummy = 0;
            cp = &dummy;
            cl = Lay{{lo[0], lo[1], lo[2]}, {lo[0], lo[1], lo[2]}};   // never read (use_cov = 0)
        }
        rt_tag_cells(lo, hi, qa.p, ql.lo, ql.hi, dp, &use_den, cp, cl.lo, cl.hi, &use_cov, tp, tl.lo, tl.hi, &mode, c.thr, c.base, c.keepfac, &rel, dirs);
    }
}

void tag_boxes(amrex::TagBoxArray& tags, const std::vector<amrex::Box>& boxes)
{
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const int lo[3] = {vb.smallEnd(0), vb.smallEnd(1), vb.smallEnd(2)}, hi[3] = {vb.bigEnd(0), vb.bigEnd(1), vb.bigEnd(2)};
        Lay tl;
        signed char* tp = tag_ptr(tags, mfi, tl);
        for (const amrex::Box& b : boxes) {
            if (!b.intersects(vb)) continue;
            const int blo[3] = {b.smallEnd(0), b.smallEnd(1), b.smallEnd(2)}, bhi[3] = {b.bigEnd(0), b.bigEnd(1), b.bigEnd(2)};
            rt_tag_box(lo, hi, tp, tl.lo, tl.hi, blo, bhi);
        }
    }
}

long count_tags(const amrex::TagBoxArray& tags)
{
    long total = 0;
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const int lo[3] = {vb.smallEnd(0), vb.smallEnd(1), vb.smallEnd(2)}, hi[3] = {vb.bigEnd(0), vb.bigEnd(1), vb.bigEnd(2)};
        auto a = tags.const_array(mfi);
        const Lay tl = bounds(a);
        long n = 0;
        rt_tag_count(lo, hi, reinterpret_cast<const signed char*>(a.p), tl.lo, tl.hi, &n);
        total += n;
    }
    amrex::ParallelAllReduce::Sum(total, amrex::ParallelContext::CommunicatorSub());
    return total;
}

long clip_tags_to_boxes(amrex::TagBoxArray& tags, const std::vector<amrex::Box>& allowed)
{
    const long before = count_tags(tags);
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const int lo[3] = {vb.smallEnd(0), vb.smallEnd(1), vb.smallEnd(2)}, hi[3] = {vb.bigEnd(0), vb.bigEnd(1), vb.bigEnd(2)};
        auto a = tags.array(mfi);
        const Lay tl = bounds(a);
        signed char* tp = reinterpret_cast<signed char*>(a.p);
        // copy of the tags of this box (valid cells), clear, copy back the cells inside the allowed boxes
        amrex::BaseFab<char> saved(tags[mfi].box(), 1);
        saved.copy<amrex::RunOn::Host>(tags[mfi]);
        auto sa = saved.array();
        const Lay sl = bounds(sa);
        tags[mfi].setVal<amrex::RunOn::Host>(amrex::TagBox::CLEAR);
        for (const amrex::Box& b : allowed) {
            if (!b.intersects(vb)) continue;
            const int blo[3] = {b.smallEnd(0), b.smallEnd(1), b.smallEnd(2)}, bhi[3] = {b.bigEnd(0), b.bigEnd(1), b.bigEnd(2)};
            rt_tag_copy_box(lo, hi, reinterpret_cast<const signed char*>(sa.p), sl.lo, sl.hi, tp, tl.lo, tl.hi, blo, bhi);
        }
    }
    return before - count_tags(tags);
}

amrex::iMultiFab make_covered_mask(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::BoxArray& fine_ba, const amrex::IntVect& ratio)
{
    amrex::iMultiFab m(ba, dm, 1, 0);
    m.setVal(0);
    if (fine_ba.empty()) return m;
    amrex::BoxArray cov = fine_ba;
    cov.coarsen(ratio);
    for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
        for (const auto& is : cov.intersections(mfi.validbox())) m[mfi].setVal<amrex::RunOn::Host>(1, is.second, 0, 1);
    }
    return m;
}

}  // namespace fdsrt
