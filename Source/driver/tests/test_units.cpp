// test_units.cpp: driver unit tests that need AMReX only (no FDS objects). Run by run_driver_tests.sh on 1 and 4 ranks, threads=1.
//
// Kernel-facing rules (M2a): (a) passive scalars: the species arrays are tested with 2 (tracked only) and 5 (with passive scalars)
// components, see test_registry; (b) only uniform Cartesian metrics are used: the layouts have one cell size per direction.
//
// Tests:
//   index_map   IR-005 round trip FDS index -> AMReX index -> back for every registered array (cell and staggered), lower-bound
//               remap, the +1 face offset, valid/window ranges, on several layouts and periodicity settings.
//   ghosts      MultiFab ghost widths and staggering per field (D-031: RHO/RHOS ng=3, ZZ/ZZS/TMP ng=2, U/V/W nodal ng=1).
//   fill        ghost fill (faces, edges, corners, periodic images; nothing beyond a non-periodic domain edge).
//   sidedata    per-box side data: valid+2 mask, interface faces open, domain faces closed, layout-independent hash.
//   registry    default field set, scratch arrays refused, passive-scalar component count.
//   tile_race   IR-007 skeleton: a reduction-free kernel gives bitwise equal results for tile sizes, box layouts and (optionally) threads.
#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "FdsAmr.H"
#include "Fields.H"
#include "SideData.H"
#include "check.H"

using namespace fdsamr;

namespace {

// ---------- synthetic layouts ----------
struct Layout {
    std::string name;
    int n[3];       // cells of the whole domain
    int cuts[3];    // boxes per direction
    int per[3];     // periodicity
};

std::vector<Layout> layouts()
{
    return {
        {"single 16x12x8 per(1,1,1)", {16, 12, 8}, {1, 1, 1}, {1, 1, 1}},
        {"2x2x1 24x18x12 per(1,0,1)", {24, 18, 12}, {2, 2, 1}, {1, 0, 1}},
        {"3x1x2 30x8x16 nonperiodic", {30, 8, 16}, {3, 1, 2}, {0, 0, 0}},
        {"2x1x2 thin-y 16x1x16 per(1,0,1)", {16, 1, 16}, {2, 1, 2}, {1, 0, 1}},
    };
}

Level0 make_level0(const Layout& L, int nranks_used)
{
    amrex::Vector<MeshInfo> meshes;
    DomainInfo dom{};
    for (int d = 0; d < 3; ++d) dom.periodic[d] = L.per[d];
    dom.cylindrical = 0;
    dom.n_tracked = 2;
    dom.n_total = 2;
    dom.nranks = nranks_used;
    const double h[3] = {0.5, 0.25, 0.125};   // one cell size per direction, all different (catches swapped directions)
    int nb = 0;
    for (int k = 0; k < L.cuts[2]; ++k)
        for (int j = 0; j < L.cuts[1]; ++j)
            for (int i = 0; i < L.cuts[0]; ++i) {
                const int ci[3] = {i, j, k};
                MeshInfo m{};
                for (int d = 0; d < 3; ++d) {
                    const int lo = L.n[d] * ci[d] / L.cuts[d];
                    const int hi = L.n[d] * (ci[d] + 1) / L.cuts[d];
                    m.ijk[d] = hi - lo;
                    m.xb[2 * d] = -1.0 + lo * h[d];
                    m.xb[2 * d + 1] = -1.0 + hi * h[d];
                }
                m.rank = nb % nranks_used;
                meshes.push_back(m);
                ++nb;
            }
    return assemble_level0(meshes, dom);
}

std::uint64_t mix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Exactly representable value that depends on the field, the AMReX (global) index and the component
double val(int fid, int a, int b, int c, int n)
{
    const std::uint64_t h = mix64((std::uint64_t)(fid + 1) * 1000003ULL ^ ((std::uint64_t)(a + 1000) << 32) ^
                                  ((std::uint64_t)(b + 1000) << 16) ^ (std::uint64_t)(c + 1000) ^ ((std::uint64_t)n << 50));
    return static_cast<double>(h & 0xFFFFFFFFFFFFULL);
}

int field_id(const std::string& name)
{
    int i = 0;
    for (const auto& s : field_table()) {
        if (name == s.name) return i;
        ++i;
    }
    return -1;
}

// ---------- index_map ----------
long test_index_map(int nranks)
{
    for (const auto& L : layouts()) {
        Level0 l0 = make_level0(L, nranks);
        Fields F(l0, 5);   // 5 = 2 tracked + 3 passive components in ZZ, ZZS
        for (const auto& s : field_table()) {
            if (!F.has(s.name)) continue;
            amrex::MultiFab& mf = F[s.name];
            const int fid = field_id(s.name);
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const amrex::Box cell_box = l0.ba[mfi.index()];
                amrex::FArrayBox& fab = mf[mfi];
                // fill every element of the grown FAB from its AMReX global index
                {
                    auto a = fab.array();
                    const amrex::Box fb = fab.box();
                    for (int n = 0; n < fab.nComp(); ++n)
                        amrex::LoopOnCpu(fb, [&](int i, int j, int k) { a(i, j, k, n) = val(fid, i, j, k, n); });
                }
                for (int window = 0; window < 2; ++window) {
                    FdsView v = make_fds_view(s, fab, cell_box, window == 1);
                    const FdsBounds b = window ? fds_window(s, cell_box) : fds_bounds(s, cell_box);
                    // lower-bound remap: the FDS lower bound maps to the AMReX index of the FAB (or window) start
                    const amrex::Box fabbox = fab.box();
                    for (int d = 0; d < 3; ++d) {
                        const int ng_used = window ? s.ng_fds : s.ng;
                        CHECK(to_amrex(s, d, cell_box.smallEnd(d), b.lb[d]) == cell_box.smallEnd(d) - ng_used);
                        CHECK(b.ext[d] == cell_box.length(d) + 2 * ng_used + (s.nodal(d) ? 1 : 0));
                        if (!window) CHECK(to_amrex(s, d, cell_box.smallEnd(d), b.lb[d]) == fabbox.smallEnd(d));
                        if (!window) CHECK(b.ext[d] == fabbox.length(d));
                        CHECK(v.lb[d] == b.lb[d] && v.ext[d] == b.ext[d]);
                    }
                    // round trip over the whole (window of the) array, every component
                    for (int n = 1; n <= fab.nComp(); ++n)
                        for (int K = b.lb[2]; K < b.lb[2] + b.ext[2]; ++K)
                            for (int J = b.lb[1]; J < b.lb[1] + b.ext[1]; ++J)
                                for (int I = b.lb[0]; I < b.lb[0] + b.ext[0]; ++I) {
                                    const int ai = to_amrex(s, 0, cell_box.smallEnd(0), I);
                                    const int aj = to_amrex(s, 1, cell_box.smallEnd(1), J);
                                    const int ak = to_amrex(s, 2, cell_box.smallEnd(2), K);
                                    if (to_fds(s, 0, cell_box.smallEnd(0), ai) != I || to_fds(s, 1, cell_box.smallEnd(1), aj) != J ||
                                        to_fds(s, 2, cell_box.smallEnd(2), ak) != K) { CHECK(false); continue; }
                                    if (v(I, J, K, n) != val(fid, ai, aj, ak, n - 1)) { CHECK(false); continue; }
                                }
                    CHECK(true);
                    // window is contiguous unless it is a strict sub-block of the FAB (RHO, RHOS)
                    CHECK(v.contiguous == (s.ng == s.ng_fds || !window));
                }
                // valid range and the +1 face offset
                int vlo[3], vhi[3];
                fds_valid_range(s, cell_box, vlo, vhi);
                for (int d = 0; d < 3; ++d) {
                    const int lo = cell_box.smallEnd(d), hi = cell_box.bigEnd(d);
                    if (s.nodal(d)) {
                        CHECK(to_amrex(s, d, lo, 0) == lo);            // FDS U(0): low face of the box
                        CHECK(to_amrex(s, d, lo, vhi[d]) == hi + 1);   // FDS U(IBAR): high face of the last cell
                        CHECK(to_amrex(s, d, lo, 1) == lo + 1);        // FDS U(1) is the high face of cell 1 = AMReX face lo+1
                    } else {
                        CHECK(to_amrex(s, d, lo, 1) == lo);
                        CHECK(to_amrex(s, d, lo, vhi[d]) == hi);
                    }
                }
            }
        }
    }
    return fdstest::report("index_map (IR-005 round trip, all registered arrays, 4 layouts)");
}

// ---------- ghost widths and staggering ----------
long test_ghosts(int nranks)
{
    Layout L = layouts()[1];
    Level0 l0 = make_level0(L, nranks);
    struct Exp { const char* name; int ng; int nodal_dir; };   // nodal_dir -1 = cell-centred
    const Exp exp[] = {{"RHO", 3, -1}, {"RHOS", 3, -1}, {"TMP", 2, -1}, {"ZZ", 2, -1}, {"ZZS", 2, -1},
                       {"U", 1, 0}, {"V", 1, 1}, {"W", 1, 2}, {"US", 1, 0}, {"VS", 1, 1}, {"WS", 1, 2},
                       {"H", 1, -1}, {"HS", 1, -1}, {"KRES", 1, -1}, {"D", 1, -1}, {"DS", 1, -1}, {"MU", 1, -1}};
    for (int ns : {2, 5}) {
        Fields F(l0, ns);
        for (const auto& e : exp) {
            CHECK_MSG(F.has(e.name), e.name);
            if (!F.has(e.name)) continue;
            const amrex::MultiFab& mf = F[e.name];
            CHECK_MSG(mf.nGrowVect() == amrex::IntVect(e.ng), e.name);
            for (int d = 0; d < 3; ++d) CHECK_MSG(mf.ixType().nodeCentered(d) == (d == e.nodal_dir), e.name);
            const bool sp = (std::strcmp(e.name, "ZZ") == 0 || std::strcmp(e.name, "ZZS") == 0);
            CHECK_MSG(mf.nComp() == (sp ? ns : 1), e.name);
            CHECK_MSG(F.spec(e.name).ng == e.ng, e.name);
            // native ghost width is never smaller than the FDS allocation
            CHECK_MSG(F.spec(e.name).ng >= F.spec(e.name).ng_fds, e.name);
        }
        // only RHO and RHOS carry the extra native layer
        for (const auto& s : field_table()) CHECK_MSG((s.ng != s.ng_fds) == (std::strcmp(s.name, "RHO") == 0 || std::strcmp(s.name, "RHOS") == 0), s.name);
    }
    return fdstest::report("ghosts (D-031 ghost widths, staggering, passive-scalar components)");
}

// ---------- ghost fill ----------
long test_fill(int nranks)
{
    for (const auto& L : layouts()) {
        Level0 l0 = make_level0(L, nranks);
        Fields F(l0, 3);
        const amrex::Box dom = l0.geom.Domain();
        for (const char* nm : {"RHO", "ZZ", "U", "V", "W", "H"}) {
            {
                const FieldSpec& s = F.spec(nm);
                const int fid = field_id(nm);
                amrex::MultiFab& mf = F[nm];
                // domain extent of the index space of this field, and periodic wrap
                int dlo[3], dn[3];   // lowest index and count of the field's index space
                for (int d = 0; d < 3; ++d) { dlo[d] = dom.smallEnd(d); dn[d] = dom.length(d) + (s.nodal(d) ? 1 : 0); }
                auto wrap = [&](int d, int a) {   // periodic image; nodal: the period is the cell count (face n = face 0)
                    if (!L.per[d]) return a;
                    const int per = dom.length(d);
                    return dlo[d] + ((a - dlo[d]) % per + per) % per;
                };
                const double sentinel = -7.0;
                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    auto a = mf.array(mfi);
                    const amrex::Box vb = mfi.validbox();
                    amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                        for (int n = 0; n < mf.nComp(); ++n)
                            a(i, j, k, n) = vb.contains(amrex::IntVect(i, j, k)) ? val(fid, wrap(0, i), wrap(1, j), wrap(2, k), n) : sentinel;
                    });
                }
                F.fill_ghosts(nm);
                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    auto a = mf.const_array(mfi);
                    const amrex::Box vb = mfi.validbox();
                    amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                        const int p[3] = {i, j, k};
                        int n_out = 0;
                        bool src_ok = true;
                        for (int d = 0; d < 3; ++d) {
                            if (p[d] < vb.smallEnd(d) || p[d] > vb.bigEnd(d)) ++n_out;
                            const int w = wrap(d, p[d]);
                            if (w < dlo[d] || w >= dlo[d] + dn[d]) src_ok = false;
                        }
                        for (int n = 0; n < mf.nComp(); ++n) {
                            double expect;
                            if (n_out == 0) expect = val(fid, wrap(0, i), wrap(1, j), wrap(2, k), n);
                            else if (!src_ok) expect = sentinel;                                   // outside a non-periodic domain
                            else expect = val(fid, wrap(0, i), wrap(1, j), wrap(2, k), n);
                            if (a(i, j, k, n) != expect) { CHECK_MSG(false, std::string(nm) + " " + L.name); return; }
                        }
                    });
                    CHECK(true);
                }
            }
        }
    }
    return fdstest::report("fill (ghost fill faces+edges+corners, periodic images, nothing outside a closed domain)");
}

// ---------- side data ----------
struct GlobalWalls {   // synthetic single-mesh wall data on the global domain
    amrex::Box dom;
    bool solid(int i, int j, int k) const { return (mix64((std::uint64_t)(i + 50) * 73856093ULL ^ (std::uint64_t)(j + 50) * 19349663ULL ^
                                                          (std::uint64_t)(k + 50) * 83492791ULL) % 7) == 0; }
};

std::uint64_t run_sidedata(const Layout& L, int nranks, bool verify)
{
    Level0 l0 = make_level0(L, nranks);
    GlobalWalls g{l0.geom.Domain()};
    const amrex::Box dom = g.dom;
    CellWallProvider prov = [&](int, const amrex::Box& vb, amrex::Array4<int> const& a) {
        amrex::LoopOnCpu(vb, [&](int i, int j, int k) {
            a(i, j, k, 0) = g.solid(i, j, k) ? 1 : 0;
            for (int f = 0; f < 6; ++f) {
                amrex::IntVect nb(i, j, k);
                nb[f / 2] += (f % 2 == 0) ? -1 : 1;
                a(i, j, k, 1 + f) = !dom.contains(nb) ? kWall : (!vb.contains(nb) ? kInterface : kNoWall);   // like FDS: wall at the domain edge, interface between meshes
            }
        });
    };
    SideData sd(l0, prov);
    if (verify) {
        auto wrapc = [&](amrex::IntVect p) {
            for (int d = 0; d < 3; ++d)
                if (L.per[d]) { const int n = dom.length(d); p[d] = dom.smallEnd(d) + ((p[d] - dom.smallEnd(d)) % n + n) % n; }
            return p;
        };
        for (amrex::MFIter mfi(sd.mask()); mfi.isValid(); ++mfi) {
            auto m = sd.mask().const_array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                const amrex::IntVect p(i, j, k);
                const bool inside = dom.contains(p);
                CHECK(m(i, j, k, 7) == (inside ? 1 : 0));
                const amrex::IntVect w = wrapc(p);
                if (!dom.contains(w)) {   // beyond a non-periodic edge: never filled
                    for (int c = 0; c < 7; ++c) CHECK(m(i, j, k, c) == 0);
                    return;
                }
                CHECK(m(i, j, k, 0) == (g.solid(w[0], w[1], w[2]) ? 1 : 0));
                for (int f = 0; f < 6; ++f) {
                    amrex::IntVect nb = w;
                    nb[f / 2] += (f % 2 == 0) ? -1 : 1;
                    CHECK(m(i, j, k, 1 + f) == (dom.contains(nb) ? 0 : 1));   // closed at the domain edge (periodic included), open elsewhere
                }
            });
        }
    }
    return sd.hash();
}

long test_sidedata(int nranks)
{
    for (const auto& L : layouts()) {
        const std::uint64_t h1 = run_sidedata(L, nranks, true);
        Layout single = L;
        single.cuts[0] = single.cuts[1] = single.cuts[2] = 1;
        const std::uint64_t h2 = run_sidedata(single, nranks, true);
        Layout fine = L;
        for (int d = 0; d < 3; ++d) fine.cuts[d] = std::max(1, std::min(L.n[d] / 4, 2 * L.cuts[d]));
        const std::uint64_t h3 = run_sidedata(fine, nranks, true);
        CHECK_MSG(h1 == h2, L.name + ": hash differs between the layout and one box");
        CHECK_MSG(h1 == h3, L.name + ": hash differs between the layout and a finer split");
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  sidedata hash [%s] = %016llx\n", L.name.c_str(), (unsigned long long)h1);
    }
    return fdstest::report("sidedata (valid+2 mask, open interfaces, closed domain faces, layout-independent hash)");
}

// ---------- registry ----------
long test_registry(int nranks)
{
    Level0 l0 = make_level0(layouts()[1], nranks);
    for (const char* n : {"FX", "FY", "FZ", "ADV_FX", "DIF_FZS", "SWORK4"}) {
        CHECK_MSG(is_per_box_scratch(n), n);
        CHECK_MSG(find_field(n) == nullptr, n);   // per-box scratch is never a registered MultiFab
    }
    for (const auto& s : field_table()) CHECK_MSG(!is_per_box_scratch(s.name), s.name);
    for (int ns : {1, 2, 5, 9}) {
        Fields F(l0, ns);
        CHECK(F["ZZ"].nComp() == ns && F["ZZS"].nComp() == ns);
        CHECK(F["RHO"].nComp() == 1 && F["U"].nComp() == 1);
        for (const auto& n : default_field_names()) CHECK_MSG(F.has(n), n);
    }
    {
        Fields F(l0, 2, {"RHO", "U"});
        CHECK(F.has("RHO") && F.has("U") && !F.has("V"));
        CHECK(F.bytes() >= 0);
    }
    return fdstest::report("registry (default set, scratch refused, passive-scalar component count)");
}

// ---------- IR-007 tile race skeleton ----------
// Frozen reduction-free kernels (no floating-point sum across the box): (1) a cell kernel with a 7-point stencil of RHO; (2) a face
// kernel U(i) = 0.5*(rho(i-1)+rho(i)) written on nodaltilebox. Compared bitwise for different tile sizes, box layouts and threads,
// and a coverage check (every valid element written exactly once) that would expose overlapping tiles.
struct Tiling { const char* name; bool tiling; amrex::IntVect tile; };

double rho_of(int i, int j, int k) { return 1.0 + 1.0e-3 * (i * 7 + j * 13 + k * 17) + 1.0e-6 * ((i * i + 3 * j * k) % 11); }

// Returns the result on a single-box MultiFab (gathered) so different layouts compare cell by cell.
void run_kernels(const Layout& L, int nranks, const Tiling& T, int nthreads, amrex::MultiFab& out_div, amrex::MultiFab& out_u, long& coverage_bad)
{
    Level0 l0 = make_level0(L, nranks);
    Fields F(l0, 2, {"RHO", "U"});
    amrex::MultiFab& rho = F["RHO"];
    amrex::MultiFab& U = F["U"];
    for (amrex::MFIter mfi(rho); mfi.isValid(); ++mfi) {
        auto a = rho.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) { a(i, j, k) = rho_of(i, j, k); });   // ghosts filled analytically (frozen input)
    }
    amrex::MultiFab div(l0.ba, l0.dm, 1, 0);
    amrex::MultiFab cnt_c(l0.ba, l0.dm, 1, 0);
    amrex::MultiFab cnt_f(U.boxArray(), l0.dm, 1, 0);
    div.setVal(0.0); cnt_c.setVal(0.0); cnt_f.setVal(0.0);
    amrex::MFItInfo info;
    if (T.tiling) info.EnableTiling(T.tile);
#ifdef AMREX_USE_OMP
    const int save = omp_get_max_threads();
    omp_set_num_threads(nthreads);
#pragma omp parallel
#else
    (void)nthreads;
#endif
    for (amrex::MFIter mfi(rho, info); mfi.isValid(); ++mfi) {
        const amrex::Box tb = mfi.tilebox();
        auto r = rho.const_array(mfi);
        auto d = div.array(mfi);
        auto cc = cnt_c.array(mfi);
        amrex::LoopOnCpu(tb, [&](int i, int j, int k) {
            d(i, j, k) = (r(i + 1, j, k) - 2.0 * r(i, j, k) + r(i - 1, j, k)) * 4.0 + (r(i, j + 1, k) - 2.0 * r(i, j, k) + r(i, j - 1, k)) * 16.0 +
                         (r(i, j, k + 1) - 2.0 * r(i, j, k) + r(i, j, k - 1)) * 64.0 + (r(i + 1, j + 1, k) - r(i - 1, j - 1, k)) * 0.25;
            cc(i, j, k) += 1.0;
        });
    }
#ifdef AMREX_USE_OMP
    omp_set_num_threads(save);
#endif
    // face kernel: iterate on U's own MFIter (nodal boxes) with nodaltilebox(0)
#ifdef AMREX_USE_OMP
    omp_set_num_threads(nthreads);
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(U, info); mfi.isValid(); ++mfi) {
        const amrex::Box tb = mfi.tilebox();   // for a nodal-in-x MultiFab tilebox is already the x-nodal tile
        auto r = rho.const_array(mfi);
        auto u = U.array(mfi);
        auto cf = cnt_f.array(mfi);
        amrex::LoopOnCpu(tb, [&](int i, int j, int k) {
            u(i, j, k) = 0.5 * (r(i - 1, j, k) + r(i, j, k));
            cf(i, j, k) += 1.0;
        });
    }
#ifdef AMREX_USE_OMP
    omp_set_num_threads(save);
#endif
    // coverage: every valid element written exactly once
    coverage_bad = 0;
    for (amrex::MFIter mfi(cnt_c); mfi.isValid(); ++mfi) {
        auto a = cnt_c.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (a(i, j, k) != 1.0) ++coverage_bad; });
    }
    for (amrex::MFIter mfi(cnt_f); mfi.isValid(); ++mfi) {
        auto a = cnt_f.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (a(i, j, k) != 1.0) ++coverage_bad; });
    }
    // gather to one box per rank-independent layout for bitwise comparison
    const amrex::Box dom = l0.geom.Domain();
    amrex::BoxArray one(dom);
    amrex::DistributionMapping dm1(amrex::Vector<int>{0});   // the single gather box always lives on rank 0
    out_div.define(one, dm1, 1, 0);
    out_div.ParallelCopy(div, 0, 0, 1);
    amrex::BoxArray one_f(amrex::convert(dom, amrex::IntVect(1, 0, 0)));
    out_u.define(one_f, dm1, 1, 0);
    // valid faces of U live on the nodal boxes; the shared face is written identically by both boxes (same expression, same inputs)
    out_u.ParallelCopy(U, 0, 0, 1);
}

bool same_bits(const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::Box& region)
{
    bool same = true;
    for (amrex::MFIter mfi(a); mfi.isValid(); ++mfi) {
        auto x = a.const_array(mfi);
        auto y = b.const_array(mfi);
        const amrex::Box bx = mfi.validbox() & region;
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) { if (std::memcmp(&x(i, j, k), &y(i, j, k), sizeof(double)) != 0) same = false; });
    }
    return same;
}

long test_tile_race(int nranks, bool thread_sweep)
{
    Layout L = layouts()[1];
    Level0 l0 = make_level0(L, nranks);
    const amrex::Box interior = amrex::grow(l0.geom.Domain(), -1);   // stencil needs valid ghosts: compare where inputs are analytic everywhere (all cells)
    (void)interior;
    amrex::MultiFab ref_div, ref_u;
    long bad = 0;
    run_kernels(L, nranks, {"no tiling", false, amrex::IntVect(0)}, 1, ref_div, ref_u, bad);
    CHECK_MSG(bad == 0, "coverage, no tiling");
    const std::vector<Tiling> tilings = {
        {"tile 8x8x8", true, amrex::IntVect(8, 8, 8)},
        {"tile 4x4x4", true, amrex::IntVect(4, 4, 4)},
        {"tile 1024000x4x4 (AMReX default)", true, amrex::IntVect(1024000, 4, 4)},
        {"tile 5x3x7 (uneven)", true, amrex::IntVect(5, 3, 7)},
        {"tile 1x1x1 (extreme)", true, amrex::IntVect(1, 1, 1)},
    };
    std::vector<int> threads = {1};
    if (thread_sweep) threads = {1, 2, 4, 8};
    for (int nt : threads)
        for (const auto& T : tilings) {
            amrex::MultiFab d, u;
            run_kernels(L, nranks, T, nt, d, u, bad);
            CHECK_MSG(bad == 0, std::string("coverage (each valid element written exactly once): ") + T.name + " threads " + std::to_string(nt));
            CHECK_MSG(same_bits(ref_div, d, l0.geom.Domain()), std::string("cell kernel bitwise: ") + T.name + " threads " + std::to_string(nt));
            CHECK_MSG(same_bits(ref_u, u, amrex::convert(l0.geom.Domain(), amrex::IntVect(1, 0, 0))), std::string("face kernel bitwise: ") + T.name + " threads " + std::to_string(nt));
        }
    // box layouts (FR-005 (i)): same physical problem, different splits, bitwise equal cell by cell
    for (const int cuts : {1, 2, 3}) {
        Layout M = L;
        M.cuts[0] = cuts; M.cuts[1] = 1; M.cuts[2] = cuts;
        amrex::MultiFab d, u;
        run_kernels(M, nranks, {"tile 8x8x8", true, amrex::IntVect(8, 8, 8)}, 1, d, u, bad);
        CHECK_MSG(bad == 0, "coverage, layout");
        CHECK_MSG(same_bits(ref_div, d, l0.geom.Domain()), "cell kernel bitwise across box layouts");
        CHECK_MSG(same_bits(ref_u, u, amrex::convert(l0.geom.Domain(), amrex::IntVect(1, 0, 0))), "face kernel bitwise across box layouts");
    }
    return fdstest::report(thread_sweep ? "tile_race (IR-007 skeleton: tile sizes, box layouts, threads 1/2/4/8)" : "tile_race (IR-007 skeleton: tile sizes, box layouts; threads=1)");
}

}  // namespace

int main(int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);
    long fails = 0;
    {
        bool thread_sweep = false;
        bool dump = false;
        for (int i = 1; i < argc; ++i) {
            if (std::strcmp(argv[i], "--thread-sweep") == 0) thread_sweep = true;
            if (std::strcmp(argv[i], "--dump-table") == 0) dump = true;
        }
        if (dump) {   // machine-readable table for the inventory cross-check (tests/check_inventory.py)
            if (amrex::ParallelDescriptor::IOProcessor())
                for (const auto& s : field_table())
                    std::printf("TABLE %s %d %d %d %d %s\n", s.name, static_cast<int>(s.stag), s.ng, s.ng_fds, s.per_scalar ? 1 : 0, s.fds_bounds);
        } else {
            const int nranks = amrex::ParallelDescriptor::NProcs();
            if (amrex::ParallelDescriptor::IOProcessor()) std::printf("driver_unit_tests: %d rank(s), OMP threads %s\n", nranks, std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset");
            fails += test_index_map(nranks);
            fails += test_ghosts(nranks);
            fails += test_fill(nranks);
            fails += test_sidedata(nranks);
            fails += test_registry(nranks);
            fails += test_tile_race(nranks, thread_sweep);
            if (amrex::ParallelDescriptor::IOProcessor()) std::printf("%s: %ld failing checks\n", fails == 0 ? "ALL PASS" : "SOME FAILED", fails);
        }
    }
    amrex::Finalize();
    return fails == 0 ? 0 : 1;
}
