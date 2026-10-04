// FluxStageRunner.cpp: see FluxStageRunner.H.
#include "FluxStageRunner.H"

namespace fdsrt {

std::vector<std::vector<FluxOverride>> FluxStageRunner::lists_for_level(const FluxAccess& fa, int l, FluxKind kind, OverrideStats* st) const
{
    const StageLevel& c = lv_[l];
    const StageLevel& f = lv_[l + 1];
    const amrex::MultiFab* fl[3] = {&fa.stage_flux(l + 1, kind, 0), &fa.stage_flux(l + 1, kind, 1), &fa.stage_flux(l + 1, kind, 2)};
    return build_flux_overrides(c.ba, c.dm, c.geom, f.ba, f.dm, f.ratio_from_parent, fl, st);
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
        (kind == FluxKind::Adv ? stats_.adv_entries : stats_.dif_entries) += st.entries;
        fa.set_flux_override(l, kind, lists);
    }
}

}  // namespace fdsrt
