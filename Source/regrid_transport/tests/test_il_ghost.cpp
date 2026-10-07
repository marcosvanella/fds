// test_il_ghost.cpp: L1251 ghost table (IlGhost.H) against a coordinate-based reimplementation of the FDS IIO_R search (init.f90:3139-3176,
// main.f90:2191-2254). Per case: every box-face ghost cell is a wall; the table's sender cells must equal the FDS cells in the FDS order; values are
// packed as s5rad_ils does, copied ILS -> ILR (the lag), and the s5rad_wall1 mean must equal the FDS mean of IL_R bitwise. Plus pair contiguity,
// NIC counts (coarse <- fine: tangential ratio product, fine <- coarse: 1), non-negativity, and the level search order on full (overlapping) levels.
#include "IlGhost.H"

#include <AMReX.H>
#include <AMReX_Loop.H>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace fdsrt;

namespace {
int fails = 0, checks = 0;
void check(bool ok, const char* what, int w = -1)
{
    ++checks;
    if (!ok && ++fails <= 20) std::printf("FAIL %s (wall %d)\n", what, w);
}

struct Case {
    const char* name;
    std::vector<IlBox> boxes;
    std::vector<amrex::IntVect> ratio;
};

amrex::IntVect rr(const Case& c, int l)
{
    amrex::IntVect r(1);
    for (int m = 0; m < l; ++m) r *= c.ratio[m];
    return r;
}

// every box-face ghost cell, faces in the order IOR = 1,-1,2,-2,3,-3
std::vector<IlWall> face_walls(const Case& c)
{
    std::vector<IlWall> w;
    for (int b = 0; b < static_cast<int>(c.boxes.size()); ++b)
        for (int ior : {1, -1, 2, -2, 3, -3}) {
            const int d = std::abs(ior) - 1;
            amrex::Box f = c.boxes[b].box;
            f.setRange(d, ior > 0 ? f.smallEnd(d) - 1 : f.bigEnd(d) + 1);
            amrex::LoopOnCpu(f, [&](int i, int j, int k) { w.push_back({b, ior, amrex::IntVect(i, j, k)}); });
        }
    return w;
}

double il_value(int box, int i, int j, int k, int n, int ib)   // sender box and sender-local cell (1-based); positive, not exactly representable sums
{
    return std::exp(-0.1 * box - 0.013 * i - 0.0071 * j - 0.0029 * k) * (1.0 + 0.11 * n + 0.05 * ib) / 3.0;
}

// FDS: two points on the wall face at +-0.475 of the receiver cell, MESH_SEPARATION_DISTANCE outside; the sender is the box containing them;
// the cells are the min..max of the two hits, listed K, J, I with I fastest. Level-0 cell size 1. Returns the cells (sender box, 1-based local).
std::vector<std::array<int, 4>> fds_cells(const Case& c, const IlWall& w)
{
    const int lr = c.boxes[w.box].level, d = std::abs(w.ior) - 1;
    const amrex::IntVect hr = rr(c, lr);
    int s = -1, lo[3] = {1 << 30, 1 << 30, 1 << 30}, hi[3] = {-(1 << 30), -(1 << 30), -(1 << 30)};
    for (int iter = 0; iter <= 1; ++iter) {
        double x[3];
        for (int e = 0; e < 3; ++e) x[e] = (w.g[e] + 0.5 + (iter * 0.95 - 0.475)) / hr[e];
        x[d] = (w.ior > 0 ? w.g[d] + 1 : w.g[d]) / double(hr[d]) - (w.ior > 0 ? 1e-6 : -1e-6) / hr[d];
        for (int b = 0; b < static_cast<int>(c.boxes.size()); ++b) {
            if (b == w.box) continue;
            const amrex::IntVect hs = rr(c, c.boxes[b].level);
            const amrex::Box& B = c.boxes[b].box;
            int cell[3];
            bool in = true;
            for (int e = 0; e < 3; ++e) {
                cell[e] = static_cast<int>(std::floor(x[e] * hs[e]));
                in = in && cell[e] >= B.smallEnd(e) && cell[e] <= B.bigEnd(e);
            }
            if (!in) continue;
            check(s < 0 || s == b, "fds: both points in one box", -1);
            s = b;
            for (int e = 0; e < 3; ++e) {
                lo[e] = std::min(lo[e], cell[e] - B.smallEnd(e) + 1);
                hi[e] = std::max(hi[e], cell[e] - B.smallEnd(e) + 1);
            }
        }
    }
    std::vector<std::array<int, 4>> v;
    if (s >= 0)
        for (int k = lo[2]; k <= hi[2]; ++k)
            for (int j = lo[1]; j <= hi[1]; ++j)
                for (int i = lo[0]; i <= hi[0]; ++i) v.push_back({s + 1, i, j, k});
    return v;
}

void run(const Case& c, int expect_pairs)
{
    const std::vector<IlWall> walls = face_walls(c);
    const IlGhostTable t = build_il_ghost_table(c.boxes, c.ratio, walls);
    const int nicr = t.nicr(), nra = 3, nsb = 2, nw = static_cast<int>(walls.size());

    // pack (s5rad_ils, radi.f90:4831) for every angle and band, then the lagged copy ILS -> ILR
    std::vector<double> ils(static_cast<size_t>(nicr) * nra * nsb);
    for (int ib = 0; ib < nsb; ++ib)
        for (int n = 0; n < nra; ++n)
            for (int ll = 0; ll < nicr; ++ll) ils[ll + nicr * (n + nra * ib)] = il_value(t.ils_box[ll], t.iio[ll], t.jjo[ll], t.kko[ll], n, ib);
    const std::vector<double> ilr = ils;

    int nint = 0;
    for (int w = 0; w < nw; ++w) {
        const auto ref = fds_cells(c, walls[w]);
        const int nic = t.nic_max[w] - t.nic_min[w] + 1;
        check(nic == static_cast<int>(ref.size()), "NIC equals the FDS count", w);
        if (nic != static_cast<int>(ref.size()) || nic == 0) continue;
        ++nint;
        for (int q = 0; q < nic; ++q) {
            const int ll = t.nic_min[w] - 1 + q;
            check(t.ils_box[ll] == ref[q][0] && t.iio[ll] == ref[q][1] && t.jjo[ll] == ref[q][2] && t.kko[ll] == ref[q][3], "sender cell and order equal FDS IIO_R", w);
        }
        const int ls = c.boxes[ref[0][0] - 1].level, lr = c.boxes[walls[w].box].level, d = std::abs(walls[w].ior) - 1;
        int expect = 1;
        if (ls > lr) for (int e = 0; e < 3; ++e) if (e != d) expect *= rr(c, ls)[e] / rr(c, lr)[e];
        check(nic == expect, "NIC: 1 from equal or coarser level, tangential ratio product from finer", w);
        for (int ib = 0; ib < nsb; ++ib)
            for (int n = 0; n < nra; ++n) {
                double g = 0.0, r = 0.0;   // s5rad_wall1 (radi.f90:4357-4361) on ILR, and the FDS mean of IL_R in IIO_R order
                for (int ll = t.nic_min[w]; ll <= t.nic_max[w]; ++ll) g = g + ilr[(ll - 1) + nicr * (n + nra * ib)];
                g = g / double(nic);
                for (const auto& e : ref) r = r + il_value(e[0], e[1], e[2], e[3], n, ib);
                r = r / double(ref.size());
                check(g == r && g >= 0.0, "sweep mean equals FDS mean bitwise, non-negative", w);
                if (ls == lr) check(g == il_value(ref[0][0], ref[0][1], ref[0][2], ref[0][3], n, ib), "equal level: raw neighbour value", w);
            }
    }
    // pairs: contiguous, one sender box each, ascending (receiver, sender); walls of a pair consecutive in caller order
    const int np = static_cast<int>(t.pair_lo.size()) - 1;
    check(np == expect_pairs && t.pair_lo.front() == 0 && t.pair_lo.back() == nicr, "pair count and cover");
    std::vector<int> recv(nicr, -1);
    for (int w = 0; w < nw; ++w)
        for (int ll = t.nic_min[w] - 1; ll < t.nic_max[w]; ++ll) recv[ll] = walls[w].box;
    std::pair<int, int> prev(-1, -1);
    for (int p = 0; p < np; ++p) {
        const std::pair<int, int> key(recv[t.pair_lo[p]], t.ils_box[t.pair_lo[p]]);
        check(prev < key, "pairs ascending (receiver, sender)");
        for (int ll = t.pair_lo[p]; ll < t.pair_lo[p + 1]; ++ll) check(recv[ll] == key.first && t.ils_box[ll] == key.second, "one pair per range");
        prev = key;
    }
    int last = 0, lastw = -1;
    for (int w = 0; w < nw; ++w) {   // within a (receiver, sender) pair the walls follow the caller's order
        if (t.nic_max[w] < t.nic_min[w]) continue;
        if (lastw >= 0 && walls[lastw].box == walls[w].box && t.ils_box[t.nic_min[lastw] - 1] == t.ils_box[t.nic_min[w] - 1])
            check(t.nic_min[w] == last + 1, "walls of a pair consecutive", w);
        lastw = w;
        last = t.nic_max[w];
    }
    std::printf("%-28s walls %4d interpolated %4d NICR %4d pairs %2d\n", c.name, nw, nint, nicr, np);
}
}  // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    {
        using B = amrex::Box;
        using I = amrex::IntVect;
        // equal level, three boxes; the low y face of box 2 spans boxes 0 and 1 (interleaved senders in caller order)
        run({"equal level 3D", {{0, B(I(0, 0, 0), I(3, 2, 1))}, {0, B(I(4, 0, 0), I(7, 2, 1))}, {0, B(I(0, 3, 0), I(7, 4, 1))}}, {}}, 6);
        // coarse/fine ratio 2: coarse | fine | coarse along x
        run({"coarse/fine 3D ratio 2", {{0, B(I(0, 0, 0), I(3, 3, 3))}, {1, B(I(8, 0, 0), I(11, 7, 7))}, {0, B(I(6, 0, 0), I(7, 3, 3))}}, {I(2)}}, 4);
        // 2D (one cell in y), ratio (4,1,4): coarse | fine along x, a coarse box above the fine one in z
        run({"coarse/fine 2D ratio 4", {{0, B(I(0, 0, 0), I(3, 0, 3))}, {1, B(I(16, 0, 0), I(23, 0, 15))}, {0, B(I(4, 0, 4), I(5, 0, 5))}}, {I(4, 1, 4)}}, 4);
        // two levels, ratio 2 then 3 (level 2 inside level 1 tiling): coarse | level-1 | level-2 along x
        run({"three levels ratio 2,3", {{0, B(I(0, 0, 0), I(1, 1, 1))}, {1, B(I(4, 0, 0), I(5, 3, 3))}, {2, B(I(18, 0, 0), I(23, 11, 11))}}, {I(2), I(3)}}, 4);

        // full levels (coarse under fine): a box takes its ghost from its own level first, then the coarser one, never from the finer one
        const std::vector<IlBox> full = {{0, B(I(0, 0, 0), I(3, 3, 3))}, {0, B(I(4, 0, 0), I(7, 3, 3))}, {1, B(I(6, 0, 0), I(9, 7, 7))}};
        const std::vector<IlWall> fw = {{0, -1, I(4, 1, 1)}, {2, 1, I(5, 3, 2)}, {2, -1, I(10, 3, 2)}, {2, -2, I(7, 8, 2)}};
        const IlGhostTable ft = build_il_ghost_table(full, {I(2)}, fw);
        check(ft.nic_min[0] == ft.nic_max[0] && ft.ils_box[ft.nic_min[0] - 1] == 2 && ft.iio[ft.nic_min[0] - 1] == 1, "full levels: coarse wall from the coarse neighbour");
        check(ft.nic_min[1] == ft.nic_max[1] && ft.ils_box[ft.nic_min[1] - 1] == 1 && ft.iio[ft.nic_min[1] - 1] == 3 && ft.jjo[ft.nic_min[1] - 1] == 2,
              "full levels: fine wall from the containing coarse cell");
        check(ft.nic_min[2] == ft.nic_max[2] && ft.ils_box[ft.nic_min[2] - 1] == 2 && ft.iio[ft.nic_min[2] - 1] == 2, "full levels: fine wall, other coarse box");
        check(ft.nic_max[3] == ft.nic_min[3] - 1, "full levels: outside the domain, no neighbour");
    }
    amrex::Finalize();
    std::printf("test_il_ghost: %d checks, %d failed\n", checks, fails);
    return fails == 0 ? 0 : 1;
}
