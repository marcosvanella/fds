// FluxStageRunner.cpp: see FluxStageRunner.H.
#include "FluxStageRunner.H"

#include <algorithm>
#include <cmath>

namespace fdsrt {

std::vector<std::vector<FluxOverride>> FluxStageRunner::lists_for_level(const FluxAccess& fa, int l, FluxKind kind, OverrideStats* st) const
{
    const StageLevel& c = lv_[l];
    const StageLevel& f = lv_[l + 1];
    const amrex::MultiFab* fl[3] = {&fa.stage_flux(l + 1, kind, 0), &fa.stage_flux(l + 1, kind, 1), &fa.stage_flux(l + 1, kind, 2)};
    return build_flux_overrides(c.ba, c.dm, c.geom, f.ba, f.dm, f.ratio_from_parent, fl, st);
}

void FluxStageRunner::match(const FluxAccess& fa, int l, FluxKind kind, const std::vector<std::vector<FluxOverride>>& lists)
{
    const int k = kind == FluxKind::Adv ? 0 : 1;
    for (int lev : {l, l + 1})
        for (int d = 0; d < 3; ++d)
            if (amrex::convert(fa.stage_flux(lev, kind, d).boxArray(), amrex::IntVect::TheCellVector()) != lv_[lev].ba) ++stats_.stale[k];
    const amrex::MultiFab* cf[3] = {&fa.stage_flux(l, kind, 0), &fa.stage_flux(l, kind, 1), &fa.stage_flux(l, kind, 2)};
    std::vector<double>& pre = stats_.pre[k];
    std::size_t b = 0;
    for (amrex::MFIter mfi(*cf[0]); mfi.isValid(); ++mfi, ++b) {
        if (b >= lists.size()) break;
        for (const FluxOverride& fo : lists[b]) {
            if (pre.size() < static_cast<std::size_t>(fo.nscal) + 1) pre.resize(fo.nscal + 1, 0.0);
            auto a = cf[fo.dir]->const_array(mfi);
            for (std::size_t e = 0; e < fo.face.size(); ++e) {
                double rho = 0;
                for (int n = 0; n < fo.nscal; ++n) {
                    const double dlt = a(fo.face[e][0], fo.face[e][1], fo.face[e][2], n) - fo.value[e * fo.nscal + n];
                    pre[n] = std::max(pre[n], std::abs(dlt));
                    rho += dlt;
                }
                pre[fo.nscal] = std::max(pre[fo.nscal], std::abs(rho));
            }
        }
    }
}

void FluxStageRunner::set_overrides(FluxAccess& fa, FluxKind kind)
{
    ++stats_.calls;
    const int top = num_levels() - 1;
    for (int l = top; l >= 0; --l) {
        if (!overwrite_ || l == top) {
            fa.set_flux_override(l, kind, {});
            continue;
        }
        OverrideStats st;
        auto lists = lists_for_level(fa, l, kind, &st);
        match(fa, l, kind, lists);
        (kind == FluxKind::Adv ? stats_.adv_entries : stats_.dif_entries) += st.entries;
        fa.set_flux_override(l, kind, lists);
    }
}

}  // namespace fdsrt
