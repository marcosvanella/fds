// Fields.cpp: see Fields.H. Kernel-facing rules (M2a): (a) passive scalars are handled here, as extra components of ZZ/ZZS (ncomp =
// N_TOTAL_SCALARS); (b) only uniform Cartesian metrics are used, R(I)/RRN(I) are dropped.
#include "Fields.H"

#include <AMReX_BLassert.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>

namespace fdsamr {

namespace {
// name, staggering, ng, ng_fds, per_scalar, FDS bounds (init.f90 line in inventory/mesh_fields.csv)
const std::vector<FieldSpec> kTable = {
    {"RHO",    Stag::Cell,  3, 2, false, "(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)"},   // init.f90:525
    {"RHOS",   Stag::Cell,  3, 2, false, "(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)"},   // :526
    {"TMP",    Stag::Cell,  2, 2, false, "(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)"},   // :524
    {"ZZ",     Stag::Cell,  2, 2, true,  "(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1,N_TOTAL_SCALARS)"},   // :527
    {"ZZS",    Stag::Cell,  2, 2, true,  "(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1,N_TOTAL_SCALARS)"},   // :528
    {"U",      Stag::FaceX, 1, 1, false, "(-1:IBP1,0:JBP1,0:KBP1)"},           // :533
    {"V",      Stag::FaceY, 1, 1, false, "(0:IBP1,-1:JBP1,0:KBP1)"},           // :534
    {"W",      Stag::FaceZ, 1, 1, false, "(0:IBP1,0:JBP1,-1:KBP1)"},           // :535
    {"US",     Stag::FaceX, 1, 1, false, "(-1:IBP1,0:JBP1,0:KBP1)"},           // :536
    {"VS",     Stag::FaceY, 1, 1, false, "(0:IBP1,-1:JBP1,0:KBP1)"},           // :537
    {"WS",     Stag::FaceZ, 1, 1, false, "(0:IBP1,0:JBP1,-1:KBP1)"},           // :538
    {"H",      Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :556
    {"HS",     Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :557
    {"KRES",   Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :558
    {"DDDT",   Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :559
    {"D",      Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :560
    {"DS",     Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :561
    {"MU",     Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :562
    {"MU_DNS", Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :563
    {"RSUM",   Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :599
    {"Q",      Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :581
    {"CSD2",   Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :572
    {"FVX",    Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :543 (cell-shaped momentum flux terms)
    {"FVY",    Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :544
    {"FVZ",    Stag::Cell,  1, 1, false, "(0:IBP1,0:JBP1,0:KBP1)"},            // :545
};

const std::vector<std::string> kDefault = {"RHO", "RHOS", "TMP", "ZZ", "ZZS", "U", "V", "W", "US", "VS", "WS",
                                           "H", "HS", "KRES", "D", "DS", "DDDT", "MU", "RSUM", "FVX", "FVY", "FVZ"};

// Per-box scratch: stays in Fortran/side data of the shim (S3), never a MultiFab (p1-findings section 3).
const char* const kScratch[] = {"FX", "FY", "FZ", "ADV_FX", "ADV_FY", "ADV_FZ", "DIF_FX", "DIF_FY", "DIF_FZ",
                                "DIF_FXS", "DIF_FYS", "DIF_FZS", "SWORK1", "SWORK2", "SWORK3", "SWORK4"};
}  // namespace

const std::vector<FieldSpec>& field_table() { return kTable; }
const std::vector<std::string>& default_field_names() { return kDefault; }

const FieldSpec* find_field(const std::string& name)
{
    for (const auto& s : kTable) if (name == s.name) return &s;
    return nullptr;
}

bool is_per_box_scratch(const std::string& name)
{
    for (const char* s : kScratch) if (name == s) return true;
    return false;
}

FdsView make_fds_view(const FieldSpec& s, amrex::FArrayBox& fab, const amrex::Box& cell_box, bool window)
{
    // The FAB box must be exactly the grown box the maps assume.
    amrex::Box fb = cell_box;
    if (s.stag != Stag::Cell) fb.surroundingNodes(static_cast<int>(s.stag) - 1);
    fb.grow(s.ng);
    AMREX_ALWAYS_ASSERT(fab.box() == fb);

    const FdsBounds nat = fds_bounds(s, cell_box);
    const FdsBounds win = window ? fds_window(s, cell_box) : nat;
    FdsView v;
    v.ncomp = fab.nComp();
    v.stride[0] = 1;
    v.stride[1] = nat.ext[0];
    v.stride[2] = static_cast<long>(nat.ext[0]) * nat.ext[1];
    v.cstride = static_cast<long>(nat.ext[0]) * nat.ext[1] * nat.ext[2];
    long off = 0;
    for (int d = 0; d < 3; ++d) {
        v.lb[d] = win.lb[d];
        v.ext[d] = win.ext[d];
        off += static_cast<long>(win.lb[d] - nat.lb[d]) * v.stride[d];
    }
    v.base = fab.dataPtr() + off;
    v.contiguous = (off == 0) && (win.ext[0] == nat.ext[0]) && (win.ext[1] == nat.ext[1]) && (win.ext[2] == nat.ext[2]);
    return v;
}

Fields::Fields(const Level0& l0, int nscalars, const std::vector<std::string>& names)
    : m_l0(l0), m_ns(nscalars)
{
    AMREX_ALWAYS_ASSERT(nscalars >= 1);
    const std::vector<std::string>& list = names.empty() ? default_field_names() : names;
    for (const auto& n : list) {
        const FieldSpec* s = find_field(n);
        if (s == nullptr) amrex::Abort("Fields: '" + n + "' is not a registered field" +
                                       (is_per_box_scratch(n) ? " (per-box scratch, not a MultiFab)" : ""));
        amrex::BoxArray ba = l0.ba;
        if (s->stag != Stag::Cell) ba.surroundingNodes(static_cast<int>(s->stag) - 1);
        auto mf = std::make_unique<amrex::MultiFab>(ba, l0.dm, s->ncomp(m_ns), amrex::IntVect(s->ng));
        mf->setVal(0.0);
        m_mf[n] = std::move(mf);
    }
}

amrex::MultiFab& Fields::operator[](const std::string& name)
{
    auto it = m_mf.find(name);
    if (it == m_mf.end()) amrex::Abort("Fields: field '" + name + "' not allocated");
    return *it->second;
}

const amrex::MultiFab& Fields::operator[](const std::string& name) const
{
    auto it = m_mf.find(name);
    if (it == m_mf.end()) amrex::Abort("Fields: field '" + name + "' not allocated");
    return *it->second;
}

const FieldSpec& Fields::spec(const std::string& name) const { return *find_field(name); }

void Fields::fill_ghosts(const std::string& name) const
{
    amrex::MultiFab& mf = const_cast<amrex::MultiFab&>((*this)[name]);
    mf.FillBoundary(0, mf.nComp(), mf.nGrowVect(), m_l0.geom.periodicity(), /*cross=*/false);
}

long Fields::bytes() const
{
    long b = 0;
    for (const auto& kv : m_mf)
        for (amrex::MFIter mfi(*kv.second); mfi.isValid(); ++mfi)
            b += static_cast<long>((*kv.second)[mfi].box().numPts()) * (*kv.second).nComp() * static_cast<long>(sizeof(double));
    return b;
}

std::vector<std::string> Fields::names() const
{
    std::vector<std::string> v;
    for (const auto& kv : m_mf) v.push_back(kv.first);
    return v;
}

}  // namespace fdsamr
