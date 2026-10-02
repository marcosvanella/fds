// Hierarchy.cpp: IR-002 grouping, input checks (FR-010, FR-013), static hierarchy and dump (R1). See Hierarchy.H.
#include "Hierarchy.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ostream>

namespace fdsrt {

namespace {

const double kTol = 1.0e-3;  // alignment tolerance in cells, same value as FDS (ALIGNMENT_TOLERANCE)
const char* kDirName[3] = {"x", "y", "z"};

std::string str(int v) { return std::to_string(v); }
std::string mesh_name(int m) { return "mesh " + str(m + 1); }

struct NMesh {  // normalised copy of a mesh
    double lo[3], hi[3], dx[3];
    int ijk[3];
};

NMesh normalise(const MeshInput& m)
{
    NMesh n;
    for (int d = 0; d < 3; ++d) {
        n.lo[d] = std::min(m.xb[2 * d], m.xb[2 * d + 1]);
        n.hi[d] = std::max(m.xb[2 * d], m.xb[2 * d + 1]);
        n.ijk[d] = m.ijk[d];
        n.dx[d] = (n.hi[d] - n.lo[d]) / m.ijk[d];
    }
    return n;
}

// Do two meshes touch at a face (zero overlap in one direction, positive overlap in the other two)? Periodic shifts are tried.
bool touch(const NMesh& a, const NMesh& b, const std::array<bool, 3>& per, const double len[3])
{
    int ns[3];
    for (int d = 0; d < 3; ++d) ns[d] = per[d] ? 3 : 1;
    for (int i0 = 0; i0 < ns[0]; ++i0)
        for (int i1 = 0; i1 < ns[1]; ++i1)
            for (int i2 = 0; i2 < ns[2]; ++i2) {
                const int s[3] = {i0 - 1 + (per[0] ? 0 : 1), i1 - 1 + (per[1] ? 0 : 1), i2 - 1 + (per[2] ? 0 : 1)};
                int ztouch = 0, npos = 0;
                for (int d = 0; d < 3; ++d) {
                    double sh = s[d] * len[d];
                    double ov = std::min(a.hi[d], b.hi[d] + sh) - std::max(a.lo[d], b.lo[d] + sh);
                    double tol = kTol * std::min(a.dx[d], b.dx[d]);
                    if (std::fabs(ov) <= tol) ++ztouch;
                    else if (ov > tol) ++npos;
                }
                if (ztouch == 1 && npos == 2) return true;
            }
    return false;
}

bool overlap(const NMesh& a, const NMesh& b)
{
    for (int d = 0; d < 3; ++d) {
        double ov = std::min(a.hi[d], b.hi[d]) - std::max(a.lo[d], b.lo[d]);
        if (ov <= kTol * std::min(a.dx[d], b.dx[d])) return false;
    }
    return true;
}

std::string ratio_list(const double r[3], const std::array<bool, 3>& hidden)
{
    std::string s;
    for (int d = 0; d < 3; ++d) {
        if (hidden[d]) continue;
        char b[32];
        std::snprintf(b, sizeof b, "%s%s=%.4g", s.empty() ? "" : ", ", kDirName[d], r[d]);
        s += b;
    }
    return s;
}

}  // namespace

Grouping group_meshes(const std::vector<MeshInput>& meshes, const AmrParams& p, const std::array<bool, 3>& periodic, Report& rep)
{
    Grouping g;
    g.periodic = periodic;
    const int n = static_cast<int>(meshes.size());
    const size_t nerr0 = rep.errors.size();
    if (n == 0) {
        rep.error("no &MESH entries");
        return g;
    }
    std::vector<NMesh> nm(n);
    for (int m = 0; m < n; ++m) {
        for (int d = 0; d < 3; ++d)
            if (meshes[m].ijk[d] < 1 || !(std::max(meshes[m].xb[2 * d], meshes[m].xb[2 * d + 1]) > std::min(meshes[m].xb[2 * d], meshes[m].xb[2 * d + 1]))) {
                rep.error(mesh_name(m) + ": IJK must be >= 1 and XB must have positive extent in direction " + kDirName[d]);
                return g;
            }
        nm[m] = normalise(meshes[m]);
    }
    for (int d = 0; d < 3; ++d) {
        g.hidden[d] = true;
        for (int m = 0; m < n; ++m) g.hidden[d] = g.hidden[d] && nm[m].ijk[d] == 1;
    }
    // Level-0 cell size: the coarsest cell size per direction. A hidden direction keeps its extent (ratio 1, one cell).
    for (int d = 0; d < 3; ++d) {
        g.dx0[d] = 0;
        for (int m = 0; m < n; ++m) g.dx0[d] = std::max(g.dx0[d], nm[m].dx[d]);
        if (g.hidden[d])
            for (int m = 0; m < n; ++m)
                if (std::fabs(nm[m].dx[d] - g.dx0[d]) > kTol * g.dx0[d]) {
                    rep.error(mesh_name(m) + " has a different thickness in the single-cell direction " + kDirName[d] + " than the other meshes");
                    return g;
                }
    }
    double len[3];  // bounding box of all meshes, for the periodic shifts of the pair check
    for (int d = 0; d < 3; ++d) {
        double lo = nm[0].lo[d], hi = nm[0].hi[d];
        for (int m = 1; m < n; ++m) {
            lo = std::min(lo, nm[m].lo[d]);
            hi = std::max(hi, nm[m].hi[d]);
        }
        len[d] = hi - lo;
    }

    // --- mesh pair check (FR-010): touching meshes of different cell size ---
    bool any_finer = false;
    int first_pair_a = -1, first_pair_b = -1;
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b) {
            double q[3];
            bool differ = false;
            for (int d = 0; d < 3; ++d) {
                q[d] = g.hidden[d] ? 1.0 : nm[a].dx[d] / nm[b].dx[d];
                if (std::fabs(q[d] - 1.0) > kTol) differ = true;
            }
            if (!differ) continue;
            any_finer = true;
            if (!touch(nm[a], nm[b], periodic, len)) continue;
            if (first_pair_a < 0) { first_pair_a = a; first_pair_b = b; }
            // ratio per direction as coarse/fine, with a sign for which mesh is finer
            int sign = 0;
            bool mixed = false, dir_dep = false, nonint = false;
            double rr[3] = {1, 1, 1};
            int rref = 0;
            for (int d = 0; d < 3; ++d) {
                if (g.hidden[d]) continue;
                int s = (q[d] > 1.0 + kTol) ? 1 : (q[d] < 1.0 - kTol ? -1 : 0);
                rr[d] = (q[d] >= 1.0) ? q[d] : 1.0 / q[d];
                if (s != 0) {
                    if (sign != 0 && s != sign) mixed = true;
                    sign = s;
                }
                if (s == 0 && sign != 0) dir_dep = true;
            }
            for (int d = 0; d < 3; ++d) {
                if (g.hidden[d]) continue;
                int s = (q[d] > 1.0 + kTol) ? 1 : (q[d] < 1.0 - kTol ? -1 : 0);
                if (s == 0 && sign != 0) dir_dep = true;
                if (std::fabs(rr[d] - std::round(rr[d])) > kTol * rr[d]) nonint = true;
                int ri = static_cast<int>(std::round(rr[d]));
                if (rref == 0) rref = ri;
                else if (ri != rref) dir_dep = true;
            }
            const std::string pair = mesh_name(a) + " and " + mesh_name(b);
            if (mixed || dir_dep)
                rep.error("refinement ratio between " + pair + " depends on direction (" + ratio_list(rr, g.hidden) +
                          "); AMR mode needs the same ratio, 2 or 4, in every direction (FR-010)");
            else if (nonint)
                rep.error("refinement ratio between " + pair + " is not an integer (" + ratio_list(rr, g.hidden) + "); AMR mode supports ratios 2 and 4 (FR-010)");
            else if (rref != 2 && rref != 4)
                rep.error("refinement ratio " + str(rref) + " between " + pair + " is not supported; AMR mode supports ratios 2 and 4 (FR-010)");
        }
    if (any_finer && !p.present) {
        std::string pr = first_pair_a >= 0 ? (" (for example " + mesh_name(first_pair_a) + " and " + mesh_name(first_pair_b) + ")") : std::string();
        rep.error("&MESH entries have unequal cell sizes" + pr + " but the input has no &AMR line; AMR mode is never inferred: add an &AMR line (for example &AMR MAX_LEVEL=1 /)");
        return g;
    }
    if (rep.errors.size() > nerr0) return g;

    // --- level of each mesh ---
    const int ml = p.present ? p.max_level : 0;
    g.cum_ratio.assign(ml + 1, IVec{1, 1, 1});
    g.ref_ratio.assign(ml + 1, IVec{1, 1, 1});
    for (int l = 1; l <= ml; ++l)
        for (int d = 0; d < 3; ++d) {
            g.ref_ratio[l][d] = g.hidden[d] ? 1 : p.ref_ratio_at(l - 1);
            g.cum_ratio[l][d] = g.cum_ratio[l - 1][d] * g.ref_ratio[l][d];
        }
    g.level_of_mesh.assign(n, -1);
    g.box_of_mesh.assign(n, IBox());
    for (int m = 0; m < n; ++m) {
        double r[3];
        int ri = 0;
        bool bad = false;
        for (int d = 0; d < 3; ++d) {
            r[d] = g.hidden[d] ? 1.0 : g.dx0[d] / nm[m].dx[d];
            if (g.hidden[d]) continue;
            int x = static_cast<int>(std::round(r[d]));
            if (std::fabs(r[d] - x) > kTol * x) bad = true;
            if (ri == 0) ri = x;
            else if (x != ri) bad = true;
        }
        if (bad) {
            rep.error(mesh_name(m) + " has cell sizes that are not the same integer fraction of the coarsest cell size in every direction (" +
                      ratio_list(r, g.hidden) + ")");
            continue;
        }
        int lev = -1;
        const int rd = g.hidden[0] ? (g.hidden[1] ? 2 : 1) : 0;  // a direction that is not hidden (all three if the input is 3-D)
        for (int l = 0; l <= ml; ++l)
            if (ri == g.cum_ratio[l][rd] || (ri == 0 && l == 0)) { lev = l; break; }
        if (lev < 0) {
            rep.error(mesh_name(m) + " has cell size 1/" + str(ri) + " of the coarsest meshes, which is not a level of the hierarchy (MAX_LEVEL=" + str(ml) +
                      ", level ratios 1" + [&] { std::string s; for (int l = 1; l <= ml; ++l) s += "," + str(g.cum_ratio[l][rd]); return s; }() +
                      "); raise MAX_LEVEL or change REF_RATIO");
            continue;
        }
        g.level_of_mesh[m] = lev;
    }
    if (rep.errors.size() > nerr0) return g;

    // --- overlaps ---
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            if (overlap(nm[a], nm[b]))
                rep.error(mesh_name(a) + " and " + mesh_name(b) + " overlap; overlapping or embedded meshes are not supported in AMR mode (IR-002)");
    if (rep.errors.size() > nerr0) return g;

    // --- level-0 domain: bounding box of the level-0 meshes ---
    bool have0 = false;
    for (int d = 0; d < 3; ++d) { g.dom_lo[d] = 1e300; g.dom_hi[d] = -1e300; }
    for (int m = 0; m < n; ++m) {
        if (g.level_of_mesh[m] != 0) continue;
        have0 = true;
        for (int d = 0; d < 3; ++d) {
            g.dom_lo[d] = std::min(g.dom_lo[d], nm[m].lo[d]);
            g.dom_hi[d] = std::max(g.dom_hi[d], nm[m].hi[d]);
        }
    }
    if (!have0) {
        rep.error("no &MESH entry has the coarsest cell size; level 0 would be empty");
        return g;
    }
    for (int d = 0; d < 3; ++d) {
        double cells = (g.dom_hi[d] - g.dom_lo[d]) / g.dx0[d];
        g.n0[d] = g.hidden[d] ? 1 : static_cast<int>(std::round(cells));
        if (!g.hidden[d] && std::fabs(cells - g.n0[d]) > kTol)
            rep.error(std::string("level-0 meshes do not end on a common cell lattice in direction ") + kDirName[d]);
    }
    if (rep.errors.size() > nerr0) return g;

    // --- index boxes ---
    int top = 0;
    for (int m = 0; m < n; ++m) {
        const int l = g.level_of_mesh[m];
        top = std::max(top, l);
        IBox b;
        for (int d = 0; d < 3; ++d) {
            if (g.hidden[d]) { b.lo[d] = 0; b.hi[d] = 0; continue; }
            double dxl = g.dx0[d] / g.cum_ratio[l][d];
            double x = (nm[m].lo[d] - g.dom_lo[d]) / dxl;
            int xi = static_cast<int>(std::round(x));
            if (std::fabs(x - xi) > kTol) {
                rep.error(mesh_name(m) + " (level " + str(l) + ") has its " + kDirName[d] + " edge off the cell lattice of its level (IR-002)");
                continue;
            }
            b.lo[d] = xi;
            b.hi[d] = xi + nm[m].ijk[d] - 1;
            int nl = g.n0[d] * g.cum_ratio[l][d];
            if (b.lo[d] < 0 || b.hi[d] >= nl)
                rep.error(mesh_name(m) + " (level " + str(l) + ") lies outside the level-0 domain in direction " + kDirName[d] + " (IR-002)");
            if (l > 0) {
                int rr = g.ref_ratio[l][d];
                if (b.lo[d] % rr != 0 || (b.hi[d] + 1) % rr != 0)
                    rep.error(mesh_name(m) + " (level " + str(l) + ") has " + kDirName[d] + " edges that are not on level " + str(l - 1) +
                              " cell faces (cells " + str(b.lo[d]) + ".." + str(b.hi[d]) + ", ratio " + str(rr) + "; D-030)");
            }
        }
        g.box_of_mesh[m] = b;
    }
    if (rep.errors.size() > nerr0) return g;
    g.meshes_of_level.assign(top + 1, std::vector<int>());
    for (int m = 0; m < n; ++m) g.meshes_of_level[g.level_of_mesh[m]].push_back(m);

    // --- gaps: level-0 cells under no mesh at any level (non-box level 0, waits for M2) ---
    IBox dom0;
    for (int d = 0; d < 3; ++d) { dom0.lo[d] = 0; dom0.hi[d] = g.n0[d] - 1; }
    std::vector<IBox> foot;
    for (int m = 0; m < n; ++m) foot.push_back(coarsen(g.box_of_mesh[m], g.cum_ratio[g.level_of_mesh[m]]));
    g.gap_boxes = subtract(dom0, foot);
    if (!g.gap_boxes.empty()) {
        long long cells = 0;
        for (const auto& b : g.gap_boxes) cells += b.numPts();
        rep.error("level 0 is not a box: " + std::to_string(cells) + " level-0 cells lie under no &MESH; non-box level 0 (gap cells, domain padding of ruling N1a) is not supported before M2");
        return g;
    }
    g.ok = true;
    return g;
}

int mlmg_coarsening(const Grouping& g, IVec& coarsest, Report& rep)
{
    int k = 99;
    for (int d = 0; d < 3; ++d) {
        if (g.hidden[d]) continue;
        int c = 0, v = g.n0[d];
        while (v % 2 == 0 && v > 1) { v /= 2; ++c; }
        k = std::min(k, c);
    }
    if (k == 99) k = 0;
    for (int d = 0; d < 3; ++d) coarsest[d] = g.hidden[d] ? 1 : g.n0[d] >> k;
    if (k < 3)
        rep.warn("WARNING: the level-0 domain " + str(g.n0[0]) + "x" + str(g.n0[1]) + "x" + str(g.n0[2]) + " can be halved only " + str(k) +
                 " times in every direction; the coarsest MLMG grid is " + str(coarsest[0]) + "x" + str(coarsest[1]) + "x" + str(coarsest[2]) +
                 " (FR-010 coarsening warning)");
    return k;
}

namespace {

std::vector<IBox> chop(const IBox& b, const IVec& mgs)
{
    std::vector<IBox> cur{b};
    for (int d = 0; d < 3; ++d) {
        std::vector<IBox> next;
        for (const IBox& c : cur) {
            for (int lo = c.lo[d]; lo <= c.hi[d]; lo += mgs[d]) {
                IBox s = c;
                s.lo[d] = lo;
                s.hi[d] = std::min(c.hi[d], lo + mgs[d] - 1);
                next.push_back(s);
            }
        }
        cur.swap(next);
    }
    return cur;
}

IBox snap(const IBox& b, const IVec& bf)
{
    IBox o;
    for (int d = 0; d < 3; ++d) {
        o.lo[d] = floor_div(b.lo[d], bf[d]) * bf[d];
        o.hi[d] = (floor_div(b.hi[d], bf[d]) + 1) * bf[d] - 1;
    }
    return o;
}

}  // namespace

bool build_static_hierarchy(const Grouping& g, const AmrParams& p, Hierarchy& h, Report& rep)
{
    const size_t nerr0 = rep.errors.size();
    if (!g.ok) return false;
    h = Hierarchy();
    h.max_level = p.present ? p.max_level : 0;
    h.top = g.top_level();
    h.periodic = g.periodic;
    h.hidden = g.hidden;
    h.n_proper = p.n_proper;
    for (int d = 0; d < 3; ++d) { h.dom_lo[d] = g.dom_lo[d]; h.dom_hi[d] = g.dom_hi[d]; }
    std::array<bool, 3> act{!g.hidden[0], !g.hidden[1], !g.hidden[2]};
    h.levels.resize(h.max_level + 1);
    for (int l = 0; l <= h.max_level; ++l) {
        Level& L = h.levels[l];
        L.level = l;
        L.ref_from_parent = g.ref_ratio[l];
        L.cum_ratio = g.cum_ratio[l];
        for (int d = 0; d < 3; ++d) {
            L.dx[d] = g.hidden[d] ? g.dx0[d] : g.dx0[d] / g.cum_ratio[l][d];
            L.domain.lo[d] = 0;
            L.domain.hi[d] = g.n0[d] * g.cum_ratio[l][d] - 1;
            L.blocking_factor[d] = g.hidden[d] ? 1 : p.blocking_factor_at(l);
            L.max_grid_size[d] = p.max_grid_size_at(l);
        }
    }
    // Blocking-factor checks (FR-010). Level 0: the domain must divide by the blocking factor (AMReX aborts otherwise, AMReX_AmrMesh.cpp
    // checkInput). Level 0 meshes are kept as given (no re-boxing), so their edges are not checked. Level l >= 1: mesh edges aligned.
    for (int d = 0; d < 3; ++d) {
        if (g.hidden[d]) continue;
        const Level& L = h.levels[0];
        if (L.blocking_factor[d] <= L.max_grid_size[d] && g.n0[d] % L.blocking_factor[d] != 0)
            rep.error(std::string("the level-0 domain is ") + str(g.n0[d]) + " cells in " + kDirName[d] + ", not divisible by BLOCKING_FACTOR " +
                      str(L.blocking_factor[d]) + " on level 0; choose a blocking factor that divides it (for example BLOCKING_FACTOR=2)");
    }
    for (int l = 1; l <= g.top_level(); ++l)
        for (int m : g.meshes_of_level[l])
            for (int d = 0; d < 3; ++d) {
                if (g.hidden[d]) continue;
                int bf = h.levels[l].blocking_factor[d];
                const IBox& b = g.box_of_mesh[m];
                if (b.lo[d] % bf != 0 || (b.hi[d] + 1) % bf != 0)
                    rep.error(mesh_name(m) + " on level " + str(l) + " spans cells " + str(b.lo[d]) + ".." + str(b.hi[d]) + " in " + kDirName[d] +
                              ", not aligned with BLOCKING_FACTOR " + str(bf) + " on level " + str(l));
            }
    if (rep.errors.size() > nerr0) return false;

    // Grids, finest level first: meshes of the level plus the cover needed by the next finer level.
    std::vector<std::vector<GridBox>> grids(h.top + 1);
    for (int l = h.top; l >= 0; --l) {
        for (int m : g.meshes_of_level[l]) grids[l].push_back({g.box_of_mesh[m], m});
        if (l == h.top) continue;
        std::vector<IBox> have;
        for (const auto& gb : grids[l]) have.push_back(gb.box);
        std::vector<GridBox> added;
        for (const auto& fine : grids[l + 1]) {
            IBox c = coarsen(fine.box, g.ref_ratio[l + 1]);
            std::vector<IBox> req;
            if (l == 0) req.push_back(c);  // level 0 is the whole domain: only the part under finer meshes has to be added
            else {
                IBox gr = snap(grow(c, p.n_proper, act), h.levels[l].blocking_factor);
                req = wrap_clip(gr, h.levels[l].domain, h.periodic);
            }
            for (const IBox& r : req)
                for (const IBox& piece : subtract(r, have)) {
                    added.push_back({piece, -1});
                    have.push_back(piece);
                }
        }
        std::sort(added.begin(), added.end(), [](const GridBox& a, const GridBox& b) { return a.box < b.box; });
        long long cells = 0;
        for (const auto& a : added) cells += a.box.numPts();
        if (l > 0 && cells > 0)
            rep.warn("WARNING: level " + str(l) + ": " + std::to_string(cells) + " cells in " + std::to_string(added.size()) +
                     " box(es) are added to this level because the finer level must be properly nested (N_PROPER=" + str(p.n_proper) + ")");
        grids[l].insert(grids[l].end(), added.begin(), added.end());
    }
    for (int l = 0; l <= h.top; ++l) {
        Level& L = h.levels[l];
        if (l == 0) {
            L.grids = grids[0];  // level-0 meshes in input order, then the cover boxes
            continue;
        }
        std::sort(grids[l].begin(), grids[l].end(), [](const GridBox& a, const GridBox& b) { return a.box < b.box; });
        for (const auto& gb : grids[l])
            for (const IBox& piece : chop(gb.box, L.max_grid_size)) L.grids.push_back({piece, gb.mesh});
    }

    // Refinable region (IR-008): finer mesh footprints plus declared boxes, per level, for cells that may be tagged.
    for (int l = 0; l < h.max_level; ++l) {
        Level& L = h.levels[l];
        for (int f = l + 1; f <= h.top; ++f) {
            IVec r;
            for (int d = 0; d < 3; ++d) r[d] = g.cum_ratio[f][d] / g.cum_ratio[l][d];
            for (int m : g.meshes_of_level[f]) L.taggable.push_back({coarsen(g.box_of_mesh[m], r), m, -1});
        }
        for (size_t k = 0; k < p.regions.size(); ++k) {
            const auto& rg = p.regions[k];
            int cap = rg.level < 0 ? p.max_level : rg.level;
            if (cap < l + 1) continue;
            IBox b;
            for (int d = 0; d < 3; ++d) {
                if (g.hidden[d]) { b.lo[d] = 0; b.hi[d] = 0; continue; }
                b.lo[d] = static_cast<int>(std::floor((rg.xb[2 * d] - g.dom_lo[d]) / L.dx[d] + kTol));
                b.hi[d] = static_cast<int>(std::ceil((rg.xb[2 * d + 1] - g.dom_lo[d]) / L.dx[d] - kTol)) - 1;
            }
            // Snap outward to the blocking factor of level l+1 (measured in level l+1 cells), then back to level l cells.
            const IVec rr = g.ref_ratio[l + 1];
            IBox s = coarsen(snap(refine(b, rr), h.levels[l + 1].blocking_factor), rr);
            s = intersect(s, L.domain);
            if (s.empty()) continue;
            if (!(s == intersect(b, L.domain)))
                rep.warn("WARNING: &AMR_REGION " + str(static_cast<int>(k) + 1) + " snapped outward to the blocking factor on level " + str(l + 1) +
                         ": level-" + str(l) + " cells (" + str(b.lo[0]) + "," + str(b.lo[1]) + "," + str(b.lo[2]) + ")-(" + str(b.hi[0]) + "," + str(b.hi[1]) + "," +
                         str(b.hi[2]) + ") became (" + str(s.lo[0]) + "," + str(s.lo[1]) + "," + str(s.lo[2]) + ")-(" + str(s.hi[0]) + "," + str(s.hi[1]) + "," +
                         str(s.hi[2]) + "); the snapped box is the refinable region");
            L.taggable.push_back({s, -1, static_cast<int>(k)});
        }
    }
    if (h.max_level > 0 && h.levels[0].taggable.empty())
        rep.warn("WARNING: the refinable region is empty (no finer &MESH and no &AMR_REGION); only level 0 will exist");
    check_nesting(h, rep);
    return rep.errors.size() == nerr0;
}

void check_nesting(const Hierarchy& h, Report& rep)
{
    std::array<bool, 3> act{!h.hidden[0], !h.hidden[1], !h.hidden[2]};
    for (int l = 0; l <= h.top; ++l) {
        const Level& L = h.levels[l];
        const std::string lv = "hierarchy check, level " + str(l) + ": ";
        int shown = 0;
        auto fail = [&](const std::string& s) { if (shown++ < 3) rep.error(lv + s); };
        long long pts = 0;
        for (size_t i = 0; i < L.grids.size(); ++i) {
            const IBox& b = L.grids[i].box;
            pts += b.numPts();
            if (b.empty() || !contains(L.domain, b)) fail("box " + str(static_cast<int>(i)) + " is empty or outside the domain");
            for (size_t j = i + 1; j < L.grids.size(); ++j)
                if (intersects(b, L.grids[j].box)) fail("boxes " + str(static_cast<int>(i)) + " and " + str(static_cast<int>(j)) + " overlap");
            if (l >= 1)
                for (int d = 0; d < 3; ++d) {
                    if (!act[d]) continue;
                    if (b.lo[d] % L.blocking_factor[d] != 0 || (b.hi[d] + 1) % L.blocking_factor[d] != 0)
                        fail("box " + str(static_cast<int>(i)) + " is not aligned with the blocking factor in " + kDirName[d]);
                    if (b.length(d) > L.max_grid_size[d]) fail("box " + str(static_cast<int>(i)) + " is larger than MAX_GRID_SIZE in " + kDirName[d]);
                }
            if (l >= 1) {
                IBox c = coarsen(b, L.ref_from_parent);
                for (const IBox& piece : wrap_clip(grow(c, h.n_proper, act), h.levels[l - 1].domain, h.periodic)) {
                    std::vector<IBox> have;
                    for (const auto& pg : h.levels[l - 1].grids) have.push_back(pg.box);
                    if (!subtract(piece, have).empty()) {
                        fail("box " + str(static_cast<int>(i)) + " is not properly nested in level " + str(l - 1) + " (N_PROPER=" + str(h.n_proper) + ")");
                        break;
                    }
                }
            }
        }
        if (l == 0 && pts != L.domain.numPts()) fail("level-0 boxes do not tile the domain");
    }
}

bool build_hierarchy_from_meshes(const std::vector<MeshInput>& meshes, const AmrParams& p, const std::array<bool, 3>& periodic, Hierarchy& out,
                                 Report& rep)
{
    Grouping g = group_meshes(meshes, p, periodic, rep);
    if (!g.ok) return false;
    IVec coarsest;
    mlmg_coarsening(g, coarsest, rep);
    return build_static_hierarchy(g, p, out, rep);
}

void dump_hierarchy(std::ostream& os, const Hierarchy& h)
{
    char buf[256];
    std::snprintf(buf, sizeof buf, "HIERARCHY max_level=%d top=%d hidden=%d%d%d periodic=%d%d%d n_proper=%d\n", h.max_level, h.top, h.hidden[0], h.hidden[1],
                  h.hidden[2], h.periodic[0], h.periodic[1], h.periodic[2], h.n_proper);
    os << buf;
    for (int l = 0; l <= h.top; ++l) {
        const Level& L = h.levels[l];
        std::snprintf(buf, sizeof buf, "LEVEL %d dx=%.9g,%.9g,%.9g ratio=%d,%d,%d domain=(%d,%d,%d)-(%d,%d,%d) bf=%d,%d,%d mgs=%d,%d,%d boxes=%d\n", l, L.dx[0], L.dx[1],
                      L.dx[2], L.ref_from_parent[0], L.ref_from_parent[1], L.ref_from_parent[2], L.domain.lo[0], L.domain.lo[1], L.domain.lo[2], L.domain.hi[0],
                      L.domain.hi[1], L.domain.hi[2], L.blocking_factor[0], L.blocking_factor[1], L.blocking_factor[2], L.max_grid_size[0], L.max_grid_size[1],
                      L.max_grid_size[2], static_cast<int>(L.grids.size()));
        os << buf;
        for (size_t i = 0; i < L.grids.size(); ++i) {
            const GridBox& gb = L.grids[i];
            std::snprintf(buf, sizeof buf, "  BOX %d (%d,%d,%d)-(%d,%d,%d) ", static_cast<int>(i), gb.box.lo[0], gb.box.lo[1], gb.box.lo[2], gb.box.hi[0], gb.box.hi[1],
                          gb.box.hi[2]);
            os << buf;
            if (gb.mesh >= 0) os << "mesh " << gb.mesh + 1 << "\n";
            else os << "added\n";
        }
    }
    for (int l = 0; l < h.max_level; ++l) {  // effective refinable region (IR-008), after snapping of declared boxes
        const Level& L = h.levels[l];
        std::snprintf(buf, sizeof buf, "REGION level %d boxes=%d\n", l, static_cast<int>(L.taggable.size()));
        os << buf;
        for (const RegionBox& rb : L.taggable) {
            std::snprintf(buf, sizeof buf, "  (%d,%d,%d)-(%d,%d,%d) ", rb.box.lo[0], rb.box.lo[1], rb.box.lo[2], rb.box.hi[0], rb.box.hi[1], rb.box.hi[2]);
            os << buf;
            if (rb.mesh >= 0) os << "finer mesh " << rb.mesh + 1 << "\n";
            else os << "region " << rb.region + 1 << " (snapped)\n";
        }
    }
}

long long TagClipper::clip_tags(int lev, std::vector<IVec>& tags)
{
    const auto& ok = h_.levels[lev].taggable;
    long long removed = 0;
    std::vector<IVec> keep;
    for (const IVec& t : tags) {
        bool in = false;
        for (const RegionBox& rb : ok) {
            const IBox& b = rb.box;
            if (t[0] >= b.lo[0] && t[0] <= b.hi[0] && t[1] >= b.lo[1] && t[1] <= b.hi[1] && t[2] >= b.lo[2] && t[2] <= b.hi[2]) { in = true; break; }
        }
        if (in) keep.push_back(t);
        else ++removed;
    }
    tags.swap(keep);
    count_[lev] += removed;
    return removed;
}

bool TagClipper::first_discard(int lev)
{
    if (count_[lev] == 0 || reported_[lev]) return false;
    reported_[lev] = true;
    return true;
}

}  // namespace fdsrt
