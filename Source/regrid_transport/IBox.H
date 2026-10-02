// IBox.H: integer index box and the few box operations the Role 3 hierarchy code needs (no AMReX dependency).
//
// An IBox is a set of cells lo..hi inclusive in the index space of one AMR level (the same convention as amrex::Box without
// index type). Keeping this class AMReX-free lets the input checks, the IR-002 grouping and the nesting checks be compiled and
// tested without the AMReX library; RegridAmrCore.cpp converts to amrex::Box at the boundary.
#ifndef FDSRT_IBOX_H_
#define FDSRT_IBOX_H_

#include <algorithm>
#include <array>
#include <vector>

namespace fdsrt {

using IVec = std::array<int, 3>;

inline int floor_div(int a, int b) { return (a >= 0) ? a / b : -((-a + b - 1) / b); }

struct IBox {
    IVec lo{0, 0, 0};
    IVec hi{-1, -1, -1};  // default-constructed box is empty

    bool empty() const { return hi[0] < lo[0] || hi[1] < lo[1] || hi[2] < lo[2]; }
    int length(int d) const { return hi[d] - lo[d] + 1; }
    long long numPts() const
    {
        if (empty()) return 0;
        return 1LL * length(0) * length(1) * length(2);
    }
    bool operator==(const IBox& o) const { return lo == o.lo && hi == o.hi; }
    bool operator<(const IBox& o) const { return (lo != o.lo) ? lo < o.lo : hi < o.hi; }
};

inline IBox intersect(const IBox& a, const IBox& b)
{
    IBox r;
    for (int d = 0; d < 3; ++d) {
        r.lo[d] = std::max(a.lo[d], b.lo[d]);
        r.hi[d] = std::min(a.hi[d], b.hi[d]);
    }
    return r;
}
inline bool intersects(const IBox& a, const IBox& b) { return !intersect(a, b).empty(); }

inline bool contains(const IBox& outer, const IBox& inner)
{
    if (inner.empty()) return true;
    for (int d = 0; d < 3; ++d)
        if (inner.lo[d] < outer.lo[d] || inner.hi[d] > outer.hi[d]) return false;
    return true;
}

// Cell-by-cell refinement: each coarse cell becomes r[d] fine cells in direction d.
inline IBox refine(const IBox& b, const IVec& r)
{
    IBox o;
    for (int d = 0; d < 3; ++d) {
        o.lo[d] = b.lo[d] * r[d];
        o.hi[d] = (b.hi[d] + 1) * r[d] - 1;
    }
    return o;
}
// Smallest coarse box that contains the fine box.
inline IBox coarsen(const IBox& b, const IVec& r)
{
    IBox o;
    for (int d = 0; d < 3; ++d) {
        o.lo[d] = floor_div(b.lo[d], r[d]);
        o.hi[d] = floor_div(b.hi[d], r[d]);
    }
    return o;
}
// Grow by n cells in the directions flagged in `dirs` (a direction with a single cell is never grown).
inline IBox grow(const IBox& b, int n, const std::array<bool, 3>& dirs)
{
    IBox o = b;
    for (int d = 0; d < 3; ++d)
        if (dirs[d]) {
            o.lo[d] -= n;
            o.hi[d] += n;
        }
    return o;
}

// a minus b, as at most six disjoint boxes.
inline std::vector<IBox> subtract(const IBox& a, const IBox& b)
{
    std::vector<IBox> out;
    if (a.empty()) return out;
    IBox c = intersect(a, b);
    if (c.empty()) {
        out.push_back(a);
        return out;
    }
    IBox rest = a;
    for (int d = 0; d < 3; ++d) {
        if (rest.lo[d] < c.lo[d]) {
            IBox s = rest;
            s.hi[d] = c.lo[d] - 1;
            out.push_back(s);
            rest.lo[d] = c.lo[d];
        }
        if (rest.hi[d] > c.hi[d]) {
            IBox s = rest;
            s.lo[d] = c.hi[d] + 1;
            out.push_back(s);
            rest.hi[d] = c.hi[d];
        }
    }
    return out;
}

// a minus the union of `others`: the part of `a` not covered, as disjoint boxes.
inline std::vector<IBox> subtract(const IBox& a, const std::vector<IBox>& others)
{
    std::vector<IBox> cur;
    if (!a.empty()) cur.push_back(a);
    for (const IBox& o : others) {
        std::vector<IBox> next;
        for (const IBox& c : cur) {
            std::vector<IBox> p = subtract(c, o);
            next.insert(next.end(), p.begin(), p.end());
        }
        cur.swap(next);
        if (cur.empty()) break;
    }
    return cur;
}

// Brings a (possibly grown) box back into the domain: pieces that stick out of a periodic direction are shifted to the other side,
// pieces that stick out of a non-periodic direction are cut off. Requires the box to stick out by less than one domain length.
inline std::vector<IBox> wrap_clip(const IBox& b, const IBox& domain, const std::array<bool, 3>& periodic)
{
    std::vector<IBox> cur{b};
    for (int d = 0; d < 3; ++d) {
        const int n = domain.length(d);
        std::vector<IBox> next;
        for (const IBox& p : cur) {
            IBox in = p;
            in.lo[d] = std::max(p.lo[d], domain.lo[d]);
            in.hi[d] = std::min(p.hi[d], domain.hi[d]);
            if (!in.empty()) next.push_back(in);
            if (!periodic[d]) continue;
            if (p.lo[d] < domain.lo[d]) {
                IBox s = p;
                s.hi[d] = std::min(p.hi[d], domain.lo[d] - 1);
                s.lo[d] += n;
                s.hi[d] += n;
                if (!s.empty()) next.push_back(s);
            }
            if (p.hi[d] > domain.hi[d]) {
                IBox s = p;
                s.lo[d] = std::max(p.lo[d], domain.hi[d] + 1);
                s.lo[d] -= n;
                s.hi[d] -= n;
                if (!s.empty()) next.push_back(s);
            }
        }
        cur.swap(next);
    }
    return cur;
}

}  // namespace fdsrt

#endif
