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
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
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
#include "ExactSum.H"
#include "LevelRegistry.H"
#include "PressureBcMap.H"

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


// D-028 / FR-005 (v): exact fixed-point sums are independent of the box split and of the rank count. A 32 x 4 x 32 periodic level is cut into 1, 2x2, 4x4 and 8x8 boxes
// (max_grid_size 32/16/8/4 of the 32 cells per side) and run on the ranks given; every split must return the SAME bits, equal to a serial reference that does not
// use the box layout. Terms are values of wide dynamic range and mixed sign, so the plain double sum is order dependent (printed, not asserted).
long test_exact_sum(int nranks)
{
    const int n[3] = {32, 4, 32};
    auto term = [&](int i, int j, int k, int c) {
        const std::uint64_t h = mix64((std::uint64_t)(i + 7) * 73856093ULL ^ (std::uint64_t)(j + 11) * 19349663ULL ^ (std::uint64_t)(k + 13) * 83492791ULL ^ ((std::uint64_t)c << 40));
        const double m = static_cast<double>(h & 0xFFFFFFFFFFFFFULL) / 4503599627370496.0 - 0.5;   // [-0.5,0.5)
        const int e = static_cast<int>((h >> 52) % 40) - 20;                                   // 40 binades
        return std::ldexp(m, e);
    };
    // serial reference in __int128 at the documented scale, plus a long double sum for a magnitude check
    double mx = 0.0; long double ld = 0.0L;
    for (int k = 0; k < n[2]; ++k) for (int j = 0; j < n[1]; ++j) for (int i = 0; i < n[0]; ++i) { const double t = term(i, j, k, 0) * 0.125; mx = std::max(mx, std::abs(t)); ld += t; }
    int ex = 0; std::frexp(mx, &ex);
    const int sc = 62 - ex;
    __int128 acc = 0;
    for (int k = 0; k < n[2]; ++k) for (int j = 0; j < n[1]; ++j) for (int i = 0; i < n[0]; ++i) acc += static_cast<__int128>(std::llround(std::ldexp(term(i, j, k, 0) * 0.125, sc)));
    const double ref = std::ldexp(static_cast<double>(acc), -sc);
    CHECK_MSG(std::abs(ref - static_cast<double>(ld)) <= 1e-12 * std::abs(static_cast<double>(ld)) + 1e-300, "exact sum reference is close to the long double sum");
    double naive_first = 0.0; bool naive_differs = false;
    std::vector<double> got;
    for (const int cuts : {1, 2, 4, 8}) {
        Layout L{"exact", {n[0], n[1], n[2]}, {cuts, 1, cuts}, {1, 1, 1}};
        Level0 l0 = make_level0(L, nranks);
        amrex::MultiFab mf(l0.ba, l0.dm, 2, 0);
        double naive = 0.0;
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto a = mf.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, 0) = term(i, j, k, 0); a(i, j, k, 1) = term(i, j, k, 1); naive += term(i, j, k, 0) * 0.125; });
        }
        amrex::ParallelAllReduce::Sum(naive, amrex::ParallelContext::CommunicatorSub());
        if (got.empty()) naive_first = naive; else if (naive != naive_first) naive_differs = true;
        const double ex0 = exact_sum(mf, 0, 0.125);
        got.push_back(ex0);
        CHECK_MSG(std::memcmp(&ex0, &ref, sizeof(double)) == 0, std::string("exact_sum equals the serial reference bitwise, cuts ") + std::to_string(cuts));
        // group form: two groups, the second is the plain sum of component 1 with the other weight; the result of group 0 must not depend on group 1
        std::vector<int> g; std::vector<double> v;
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto a = mf.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { g.push_back(0); v.push_back(a(i, j, k, 0) * 0.125); g.push_back(1); v.push_back(a(i, j, k, 1) * 1024.0); });
        }
        const std::vector<double> gs = exact_group_sums(2, g, v);
        CHECK_MSG(std::memcmp(&gs[0], &ref, sizeof(double)) == 0, std::string("exact_group_sums group 0 equals the reference, cuts ") + std::to_string(cuts));
        // cancellation: adjacent pairs (x, -x) inside and across box cuts sum to exactly 0
        amrex::MultiFab zero(l0.ba, l0.dm, 1, 0);
        for (amrex::MFIter mfi(zero); mfi.isValid(); ++mfi) {
            auto a = zero.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { const double t = term(i, j, k, 1); a(i, j, k) = (i % 2 == 0) ? t : -term(i - 1, j, k, 1); });   // pairs (i even, i+1 odd) cancel exactly
        }
        const double zs = exact_sum(zero, 0, 1.0);
        CHECK_MSG(zs == 0.0, std::string("exact cancellation (x + (-x) over box cuts) is exactly 0, cuts ") + std::to_string(cuts) + " got " + std::to_string(zs));
        // mean of a known field: constant 3.0 over 32*4*32 cells of volume 2^-7 is exactly 3*4096/128 = 96
        amrex::MultiFab cst(l0.ba, l0.dm, 1, 0);
        cst.setVal(3.0);
        CHECK_MSG(exact_sum(cst, 0, 0.0078125) == 96.0, std::string("sum of a constant field is exact, cuts ") + std::to_string(cuts));
    }
    for (std::size_t q = 1; q < got.size(); ++q) CHECK_MSG(std::memcmp(&got[0], &got[q], sizeof(double)) == 0, "exact sum bitwise equal across box splits");
    if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  exact_sum: %.17g; plain double sum %s across splits (ranks %d)\n", got[0], naive_differs ? "DIFFERS" : "happens to agree", nranks);
    return fdstest::report("exact_sum (D-028: bitwise across box splits 32/16/8/4 cells; ranks as launched)");
}


// S9: per-level registry (Role 3 plan section 4 items 1, 2, 5; RegridInterface.H LevelListener). Level 0 is a 16x8x16 domain in 2x1x2 boxes (periodic in y), level 1 a
// refined (ratio 2) region over its lower-left part cut into two boxes: make_level, covered mask, layout-built SideData at coarse-fine and box-box faces, exact sums over
// uncovered cells (bitwise independent of the fine box split), remake_level with the old objects readable, clear_level.
long test_levels(int nranks)
{
    Layout L0{"lv", {16, 8, 16}, {2, 1, 2}, {0, 1, 0}};
    Level0 l0 = make_level0(L0, nranks);
    Fields F0(l0, 2);
    SideData sd0(l0, layout_cell_walls(l0));
    LevelRegistry reg(l0.dom, 2);
    reg.adopt_level0(l0, F0, sd0);
    CHECK_MSG(reg.num_levels() == 1 && reg.has_level(0) && !reg.has_level(1), "registry starts with level 0 only");
    CHECK_MSG(reg.fds_bound(0), "level 0 is FDS-bound (mesh list present)");
    CHECK_MSG(reg.covered_mask(0) == nullptr, "no covered mask without a finer level");

    const amrex::IntVect rr(2, 2, 2);
    auto fine_layout = [&](int chop) {
        fdsrt::LevelLayout fl;
        fl.level = 1; fl.ref_ratio_from_parent = rr;
        amrex::Box fdom = amrex::refine(l0.geom.Domain(), rr);
        amrex::RealBox rb(l0.geom.ProbLo(), l0.geom.ProbHi());
        amrex::Array<int, 3> per{l0.dom.periodic[0], l0.dom.periodic[1], l0.dom.periodic[2]};
        fl.geom.define(fdom, rb, 0, per);
        amrex::BoxList bl;
        // refined part: fine cells (0..31, 0..7, 0..15) = level-0 cells (0..15, 0..3, 0..7)
        for (int q = 0; q < chop; ++q) bl.push_back(amrex::Box(amrex::IntVect(32 * q / chop, 0, 0), amrex::IntVect(32 * (q + 1) / chop - 1, 7, 15)));
        fl.ba.define(bl);
        fl.dm.define(fl.ba);
        return fl;
    };
    fdsrt::LevelLayout f2 = fine_layout(2);
    reg.make_level(f2);
    CHECK_MSG(reg.num_levels() == 2 && reg.has_level(1), "make_level(1)");
    CHECK_MSG(!reg.fds_bound(1), "level 1 has no FDS binding (layout only)");
    CHECK_MSG(reg.level(1).level == 1 && reg.level(1).ref_ratio_from_parent == rr, "level index and ratio recorded");
    CHECK_MSG(std::abs(reg.level(1).dx[0] - 0.5 * l0.dx[0]) < 1e-15 && std::abs(reg.level(1).dx[2] - 0.5 * l0.dx[2]) < 1e-15, "fine cell size is half the coarse one");
    CHECK_MSG(reg.fields(1).has("RHO") && reg.fields(1)["RHO"].boxArray() == f2.ba && reg.fields(1)["RHO"].nGrow() == 3, "Fields of level 1 on the level's BoxArray, D-031 ghost width 3 for RHO");
    CHECK_MSG(reg.fields(1)["ZZ"].nComp() == 2, "ZZ carries N_TOTAL_SCALARS components on level 1");
    CHECK_MSG(reg.n_side_rebuilds == 1, "make_level rebuilt the level-1 SideData");

    // SideData of level 1: wall at the domain edge (periodic edge too), open at box-box and coarse-fine faces, source flag 1 inside the domain
    {
        const amrex::iMultiFab& m = reg.side_data(1).mask();
        for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
            const auto a = m.const_array(mfi);
            const amrex::Box vb = mfi.validbox();
            auto at = [&](int i, int j, int k, int c) { return vb.contains(amrex::IntVect(i, j, k)) ? a(i, j, k, c) : -99; };
            const int w1 = at(0, 0, 0, 1);    if (w1 != -99) CHECK_MSG(w1 == 1, "level 1: -x face at the domain edge is a wall");
            const int wy = at(0, 0, 0, 3);    if (wy != -99) CHECK_MSG(wy == 1, "level 1: periodic y edge keeps the wall code (as FDS)");
            const int bb = at(15, 0, 0, 2);   if (bb != -99) CHECK_MSG(bb == 0, "level 1: box-box +x face is open");
            const int cf = at(0, 7, 0, 4);    if (cf != -99) CHECK_MSG(cf == 0, "level 1: coarse-fine +y face is open");
            const int cz = at(0, 0, 15, 6);   if (cz != -99) CHECK_MSG(cz == 0, "level 1: coarse-fine +z face is open");
            const int sf = at(0, 0, 0, 0);    if (sf != -99) CHECK_MSG(sf == 0, "level 1: no solid cells");
            const int src = a(vb.smallEnd(0), vb.smallEnd(1), vb.smallEnd(2), 7); CHECK_MSG(src == 1, "level 1: source flag set inside the domain");
        }
    }
    // covered mask of level 0: the refined part, 16 x 4 x 8 = 512 cells
    {
        const amrex::iMultiFab* cov = reg.covered_mask(0);
        CHECK_MSG(cov != nullptr && reg.covered_mask(1) == nullptr, "level 0 covered mask exists, the finest level has none");
        long n = 0;
        if (cov) for (amrex::MFIter mfi(*cov); mfi.isValid(); ++mfi) { const auto a = cov->const_array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { n += a(i, j, k); }); }
        amrex::ParallelAllReduce::Sum(n, amrex::ParallelContext::CommunicatorSub());
        CHECK_MSG(n == 512, "covered cells of level 0 = 512, got " + std::to_string(n));
    }
    // exact sums over uncovered cells
    auto term = [&](int i, int j, int k, int c) {
        const std::uint64_t h = mix64((std::uint64_t)(i + 5) * 73856093ULL ^ (std::uint64_t)(j + 3) * 19349663ULL ^ (std::uint64_t)(k + 9) * 83492791ULL ^ ((std::uint64_t)c << 40));
        return std::ldexp(static_cast<double>(h & 0xFFFFFFFFFFFFFULL) / 4503599627370496.0 - 0.5, static_cast<int>((h >> 52) % 24) - 12);
    };
    amrex::MultiFab& R0 = F0["RHO"];
    R0.setVal(0.0);
    for (amrex::MFIter mfi(R0); mfi.isValid(); ++mfi) { auto a = R0.array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = term(i, j, k, 0); }); }
    const double vol0 = l0.dx[0] * l0.dx[1] * l0.dx[2];
    {
        const double all_plain = exact_sum(R0, 0, vol0);
        const double nomask = exact_sum_uncovered(R0, 0, vol0, nullptr);
        CHECK_MSG(std::memcmp(&all_plain, &nomask, sizeof(double)) == 0, "exact_sum_uncovered with no mask is bitwise exact_sum");
        // serial reference over the uncovered cells (fixed point at the documented scale)
        double mx = 0.0; long double ld = 0.0L;
        auto unc = [&](int i, int j, int k) { return !(i <= 15 && j <= 3 && k <= 7); };
        for (int k = 0; k < 16; ++k) for (int j = 0; j < 8; ++j) for (int i = 0; i < 16; ++i) if (unc(i, j, k)) { const double t = term(i, j, k, 0) * vol0; mx = std::max(mx, std::abs(t)); ld += t; }
        int ex = 0; std::frexp(mx, &ex); const int sc = 62 - ex;
        __int128 acc = 0;
        for (int k = 0; k < 16; ++k) for (int j = 0; j < 8; ++j) for (int i = 0; i < 16; ++i) if (unc(i, j, k)) acc += static_cast<__int128>(std::llround(std::ldexp(term(i, j, k, 0) * vol0, sc)));
        const double ref = std::ldexp(static_cast<double>(acc), -sc);
        const double got = exact_sum_uncovered(R0, 0, vol0, reg.covered_mask(0));
        CHECK_MSG(std::memcmp(&got, &ref, sizeof(double)) == 0, "exact_sum_uncovered equals the serial reference over the uncovered cells (bitwise)");
        CHECK_MSG(std::abs(ref - static_cast<double>(ld)) <= 1e-12 * std::abs(static_cast<double>(ld)), "reference close to the long double sum");
        const double pr = exact_sum_product_uncovered(R0, 0, R0, 0, 1.0, reg.covered_mask(0)); (void)pr;
    }
    // hierarchy sum: level 1 holds the piecewise-constant prolongation of level 0; the sum over the uncovered coarse cells plus all fine cells (weight vol/8) is the level-0 sum
    auto fill_fine = [&]() {
        amrex::MultiFab& R1 = reg.fields(1)["RHO"];
        R1.setVal(0.0);
        for (amrex::MFIter mfi(R1); mfi.isValid(); ++mfi) { auto a = R1.array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = term(i / 2, j / 2, k / 2, 0); }); }
    };
    fill_fine();
    double h2;
    {
        std::vector<const amrex::MultiFab*> mf{&R0, &reg.fields(1)["RHO"]};
        std::vector<double> w{vol0, vol0 / 8.0};
        std::vector<const amrex::iMultiFab*> cv{reg.covered_mask(0), nullptr};
        h2 = exact_sum_hierarchy(mf, 0, w, cv);
        const double lvl0 = exact_sum(R0, 0, vol0);
        CHECK_MSG(std::abs(h2 - lvl0) <= 1e-12 * std::abs(lvl0) + 1e-300, "hierarchy exact sum equals the level-0 sum of the same field (rounding only)");
    }
    // remake_level with another fine box split: the old objects stay readable, the hierarchy sum is bitwise the same (decomposition independent)
    {
        fdsrt::LevelLayout f4 = fine_layout(4);
        reg.remake_level(f4);
        CHECK_MSG(reg.retired_fields(1) != nullptr && reg.retired_level(1) != nullptr && reg.retired_level(1)->ba == f2.ba, "remake_level keeps the previous level readable");
        CHECK_MSG(reg.level(1).ba == f4.ba && reg.level(1).ba.size() == 4, "remake_level installed the new layout");
        reg.release_retired(1);
        CHECK_MSG(reg.retired_fields(1) == nullptr, "release_retired frees the previous level");
        fill_fine();
        std::vector<const amrex::MultiFab*> mf{&R0, &reg.fields(1)["RHO"]};
        std::vector<double> w{vol0, vol0 / 8.0};
        std::vector<const amrex::iMultiFab*> cv{reg.covered_mask(0), nullptr};
        const double h4 = exact_sum_hierarchy(mf, 0, w, cv);
        CHECK_MSG(std::memcmp(&h2, &h4, sizeof(double)) == 0, "hierarchy exact sum is bitwise independent of the fine box split (2 vs 4 boxes)");
        CHECK_MSG(reg.n_remake == 1 && reg.n_side_rebuilds == 2, "one remake, SideData rebuilt per make/remake");
    }
    // level 0 cannot be remade differently; the same layout is accepted
    {
        fdsrt::LevelLayout l0l; l0l.level = 0; l0l.geom = l0.geom; l0l.ba = l0.ba; l0l.dm = l0.dm;
        reg.make_level(l0l);
        CHECK_MSG(reg.num_levels() == 2, "make_level(0) with the adopted layout is a no-op");
    }
    // D-058 (R4): regrid bracket, old objects of remade and of cleared levels stay readable until end_regrid
    {
        CHECK_MSG(!reg.in_regrid(), "not in a regrid");
        reg.begin_regrid();
        CHECK_MSG(reg.in_regrid() && reg.n_begin_regrid == 1, "begin_regrid sets the bracket");
        const amrex::BoxArray ba_before = reg.level(1).ba;
        reg.clear_level(1);
        CHECK_MSG(reg.num_levels() == 1 && !reg.has_level(1), "clear_level inside the bracket removes the level");
        CHECK_MSG(reg.retired_level(1) != nullptr && reg.retired_fields(1) != nullptr && reg.retired_level(1)->ba == ba_before, "the cleared level stays readable inside the bracket");
        CHECK_MSG(reg.retired_fields(1)->has("RHO") && reg.retired_fields(1)->operator[]("RHO").boxArray() == ba_before, "the cleared level's fields are intact");
        reg.make_level(f2);
        CHECK_MSG(reg.has_level(1) && reg.level(1).ba == f2.ba, "make_level of the same level number inside the bracket");
        reg.end_regrid();
        CHECK_MSG(!reg.in_regrid() && reg.retired_level(1) == nullptr && reg.retired_fields(1) == nullptr && reg.n_end_regrid == 1, "end_regrid frees everything retired");
        reg.begin_regrid();
        fdsrt::LevelLayout f4c = fine_layout(4);
        reg.remake_level(f4c);
        CHECK_MSG(reg.retired_level(1) != nullptr && reg.retired_level(1)->ba == f2.ba && reg.level(1).ba == f4c.ba, "remake inside the bracket keeps the previous level");
        reg.end_regrid();
        CHECK_MSG(reg.retired_level(1) == nullptr, "end_regrid frees the remade level's previous objects (no release_retired call needed)");
        // runtime build on a different DistributionMapping (round robin, reversed) of the same BoxArray
        fdsrt::LevelLayout f4r = fine_layout(4);
        amrex::Vector<int> pm(f4r.ba.size());
        for (int i = 0; i < static_cast<int>(f4r.ba.size()); ++i) pm[i] = (static_cast<int>(f4r.ba.size()) - 1 - i) % nranks;
        f4r.dm = amrex::DistributionMapping(pm);
        reg.begin_regrid(); reg.remake_level(f4r); reg.end_regrid();
        CHECK_MSG(reg.level(1).dm.ProcessorMap() == f4r.dm.ProcessorMap() && reg.fields(1)["RHO"].DistributionMap().ProcessorMap() == f4r.dm.ProcessorMap(), "Fields are built on the level's own DistributionMapping at run time");
    }
    // D-058 (R4): initial fields on level 1 by direct evaluation at cell and face centres
    {
        LevelRegistry::InitialFill ini;
        ini.cell = [](double x, double y, double z, double* rho, double* tmp, double* zz) { *rho = 1.0 + x + 2.0 * y + 3.0 * z; *tmp = 300.0 + x; zz[0] = 0.25; zz[1] = 0.75; };
        ini.velocity = [](int d, double x, double y, double z) { return (d + 1) * (x + 10.0 * y + 100.0 * z); };
        reg.fill_initial_level(1, ini);
        const Level& lv = reg.level(1);
        const double* plo = lv.geom.ProbLo();
        long bad = 0, n = 0;
        for (amrex::MFIter mfi(reg.fields(1)["RHO"]); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            const auto r = reg.fields(1)["RHO"].const_array(mfi); const auto rs = reg.fields(1)["RHOS"].const_array(mfi);
            const auto t = reg.fields(1)["TMP"].const_array(mfi); const auto z = reg.fields(1)["ZZ"].const_array(mfi); const auto zs = reg.fields(1)["ZZS"].const_array(mfi);
            amrex::LoopOnCpu(vb, [&](int i, int j, int k) {
                const double x = plo[0] + (i + 0.5) * lv.dx[0], y = plo[1] + (j + 0.5) * lv.dx[1], zc = plo[2] + (k + 0.5) * lv.dx[2];
                ++n;
                if (r(i, j, k) != 1.0 + x + 2.0 * y + 3.0 * zc || rs(i, j, k) != r(i, j, k) || t(i, j, k) != 300.0 + x || z(i, j, k, 0) != 0.25 || z(i, j, k, 1) != 0.75 || zs(i, j, k, 1) != 0.75) ++bad;
            });
        }
        const char* fn[3] = {"U", "V", "W"}; const char* fs[3] = {"US", "VS", "WS"};
        for (int d = 0; d < 3; ++d)
            for (amrex::MFIter mfi(reg.fields(1)[fn[d]]); mfi.isValid(); ++mfi) {
                const amrex::Box fb = mfi.validbox();
                const auto a = reg.fields(1)[fn[d]].const_array(mfi); const auto as = reg.fields(1)[fs[d]].const_array(mfi);
                amrex::LoopOnCpu(fb, [&](int i, int j, int k) {
                    const int ix[3] = {i, j, k}; double x[3];
                    for (int e = 0; e < 3; ++e) x[e] = plo[e] + (ix[e] + (e == d ? 0.0 : 0.5)) * lv.dx[e];
                    ++n;
                    if (a(i, j, k) != (d + 1) * (x[0] + 10.0 * x[1] + 100.0 * x[2]) || as(i, j, k) != a(i, j, k)) ++bad;
                });
            }
        amrex::ParallelAllReduce::Sum(bad, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Sum(n, amrex::ParallelContext::CommunicatorSub());
        CHECK_MSG(n > 0 && bad == 0, "fill_initial_level: " + std::to_string(bad) + " of " + std::to_string(n) + " values differ from the direct evaluation at cell and face centres");
    }
    reg.clear_level(1);
    CHECK_MSG(reg.num_levels() == 1 && !reg.has_level(1) && reg.covered_mask(0) == nullptr && reg.n_clear == 2, "clear_level(1) removes the level and the covered mask below it");
    return fdstest::report("levels (S9: registry make/remake/clear, covered mask, layout SideData, uncovered exact sums; ranks as launched)");
}

// S9 thin-direction check: FDS pressure codes -> pb::BC strings (PressureBcMap.H). amrex::FFT::Poisson ignores a one-cell direction, so a Dirichlet face there must not be passed on.
long test_pressure_bc_map(int)
{
    auto str = [](int d, int lo, int hi, int n, bool per, bool ign) { const DirBc m = map_pressure_bc_direction(d, lo, hi, n, per, ign); return m.error.empty() ? bc_string(m) : std::string("ERR"); };
    // thick directions: code -> string
    CHECK_MSG(str(0, 3, 3, 32, false, false) == "NN", "closed N-N");
    CHECK_MSG(str(0, 1, 1, 32, false, false) == "DD", "open D-D");
    CHECK_MSG(str(2, 2, 2, 32, false, false) == "DN", "code 2 = D low, N high");
    CHECK_MSG(str(2, 4, 4, 32, false, false) == "ND", "code 4 = N low, D high");
    CHECK_MSG(str(0, 0, 0, 32, true, false) == "PP", "level periodic");
    CHECK_MSG(str(1, 0, 0, 32, false, false) == "ERR", "code 0 on a thick non-periodic direction is refused");
    CHECK_MSG(str(0, 5, 5, 32, false, false) == "ERR" && str(0, -1, -1, 32, false, false) == "ERR", "codes 5/6 and missing faces are refused");
    // one-cell y of a TWO_D case: FDS has no y operator, every code maps to Neumann (no coupling), periodic stays periodic
    for (int c = 1; c <= 4; ++c) CHECK_MSG(str(1, c, c, 1, false, true) == "NN", std::string("TWO_D y, code ") + std::to_string(c) + " -> NN");
    CHECK_MSG(str(1, 0, 0, 1, true, true) == "PP" && str(1, 0, 0, 1, false, true) == "NN", "TWO_D y periodic -> PP, code 0 -> NN");
    CHECK_MSG(!map_pressure_bc_direction(1, 1, 1, 1, false, true).note.empty(), "the Dirichlet-to-Neumann change in the ignored direction is reported");
    // code 0 on a non-periodic thick direction (soborot_*, bound_test_*: no pressure solve): refused, and the message names the case type and the AMR-mode limit
    {
        const DirBc m = map_pressure_bc_direction(0, 0, 0, 32, false, false);
        CHECK_MSG(m.error.find("no pressure solve") != std::string::npos && m.error.find("AMR mode") != std::string::npos && m.error.find("soborot_") != std::string::npos && m.error.find("PERIODIC_TEST=13") != std::string::npos,
                  "code 0 on a non-periodic direction: refusal message names the case type and says that no pressure solve is supported in AMR mode");
    }
    // one-cell x or z (FDS solves it with the Dirichlet term): Neumann/periodic exact, Dirichlet refused
    CHECK_MSG(str(0, 3, 3, 1, false, false) == "NN" && str(2, 0, 0, 1, false, false) == "NN" && str(0, 0, 0, 1, true, false) == "PP", "one-cell x/z: N and P are exact");
    for (int c : {1, 2, 4}) CHECK_MSG(str(2, c, c, 1, false, false) == "ERR" && str(0, c, c, 1, false, false) == "ERR", std::string("one-cell x/z with code ") + std::to_string(c) + " (Dirichlet face) is refused");
    // ignored (TWO_D y) direction follows an all-open pair of solved directions, otherwise stays Neumann; periodic is never changed
    {
        auto dd = [](int c) { return map_pressure_bc_direction(0, c, c, 32, false, false); };
        DirBc y = map_pressure_bc_direction(1, 3, 3, 1, false, true);
        CHECK_MSG(ignored_direction_follows_open(y, dd(1), dd(1)) && bc_string(y) == "DD" && !y.note.empty(), "TWO_D, x and z both D-D: y takes D (six uniform open faces for the selector)");
        DirBc y2 = map_pressure_bc_direction(1, 3, 3, 1, false, true);
        CHECK_MSG(!ignored_direction_follows_open(y2, dd(1), dd(4)) && bc_string(y2) == "NN", "TWO_D, z only half open: y stays N");
        DirBc y3 = map_pressure_bc_direction(1, 3, 3, 1, false, true);
        CHECK_MSG(!ignored_direction_follows_open(y3, dd(3), dd(3)) && bc_string(y3) == "NN", "TWO_D, closed box: y stays N");
        DirBc y4 = map_pressure_bc_direction(1, 0, 0, 1, true, true);
        CHECK_MSG(!ignored_direction_follows_open(y4, dd(1), dd(1)) && bc_string(y4) == "PP", "periodic ignored direction is never changed");
    }
    return fdstest::report("pressure_bc_map (S9 thin-direction rules: FDS codes -> FFT boundary strings)");
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
            fails += test_exact_sum(nranks);
            fails += test_levels(nranks);
            fails += test_pressure_bc_map(nranks);
            if (amrex::ParallelDescriptor::IOProcessor()) std::printf("%s: %ld failing checks\n", fails == 0 ? "ALL PASS" : "SOME FAILED", fails);
        }
    }
    amrex::Finalize();
    return fails == 0 ? 0 : 1;
}
