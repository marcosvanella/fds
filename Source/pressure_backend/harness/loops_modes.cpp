// loops_modes.cpp: pb_fds_loops, bitwise check of the host translations in FdsPressureLoops.cpp against the Fortran reference
// (tests/fds_loops/ref_loops.f90.in with the verbatim upstream loops, run by tests/fds_loops_test.py).
//   pb_fds_loops dir=<directory of the case files written by ref_loops> [verbose=1]
// Each case file holds the inputs and the outputs of the verbatim loop. This program rebuilds the inputs, runs the C++ function and compares
// EVERY element of every output array bitwise (the untouched parts of an output array carry random or sentinel values, so a store that
// is too large or too small is a difference). Prints max abs and max rel difference per loop and a SUMMARY line; exit 0 only if all cases are bitwise equal.
// No AMReX, no MPI.
#include "FdsPressureLoops.H"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace pressure_backend::fdsloops;

namespace {
struct Rec { int kind = 0; std::vector<int> dims; std::vector<double> d; std::vector<int> i; };
using Case = std::map<std::string, Rec>;

bool read_case (std::string const& fn, Case& c)
{
    std::ifstream f(fn, std::ios::binary);
    if (!f) return false;
    while (f.peek() != EOF) {
        char nm[16]; std::int32_t kind, rank;
        f.read(nm, 16); f.read(reinterpret_cast<char*>(&kind), 4); f.read(reinterpret_cast<char*>(&rank), 4);
        if (!f) return false;
        Rec r; r.kind = kind; std::size_t n = 1;
        for (int k = 0; k < rank; ++k) { std::int32_t d; f.read(reinterpret_cast<char*>(&d), 4); r.dims.push_back(d); n *= static_cast<std::size_t>(d); }
        if (kind == 2) { r.d.resize(n); f.read(reinterpret_cast<char*>(r.d.data()), static_cast<std::streamsize>(8 * n)); }
        else { std::vector<std::int32_t> t(n); f.read(reinterpret_cast<char*>(t.data()), static_cast<std::streamsize>(4 * n)); r.i.assign(t.begin(), t.end()); }
        if (!f) return false;
        std::string s(nm, 16); while (!s.empty() && s.back() == ' ') s.pop_back();
        c[s] = std::move(r);
    }
    return true;
}

static bool g_debug = false;
static bool g_uu_lb0 = false;   // negative control: read the FDS velocity arrays with lower bound 0 (they start at -1)
struct Stat { long cases = 0, bad_cases = 0; long n = 0, nbad = 0; double max_abs = 0, max_rel = 0; };

void compare (Stat& st, std::string const& what, std::vector<double> const& got, std::vector<double> const& ref, bool verbose, std::string const& file)
{
    bool bad = got.size() != ref.size();
    if (bad) { std::printf("  %s %s: size %zu vs %zu\n", file.c_str(), what.c_str(), got.size(), ref.size()); st.nbad++; return; }
    long nb = 0;
    for (std::size_t q = 0; q < got.size(); ++q) {
        std::uint64_t a, b; std::memcpy(&a, &got[q], 8); std::memcpy(&b, &ref[q], 8);
        st.n++;
        if (a != b) {
            ++nb; st.nbad++;
            if (g_debug && nb <= 6) std::printf("    [%zu] got %.17g ref %.17g\n", q, got[q], ref[q]);
            const double ad = std::abs(got[q] - ref[q]), m = std::max(std::abs(got[q]), std::abs(ref[q]));
            st.max_abs = std::max(st.max_abs, ad); st.max_rel = std::max(st.max_rel, m > 0 ? ad / m : 0.0);
        }
    }
    if (nb && verbose) std::printf("  %s %s: %ld of %zu elements differ\n", file.c_str(), what.c_str(), nb, got.size());
    if (nb) st.bad_cases += 0;
}

F3 f3 (Rec& r, int l0, int l1, int l2) { F3 f; f.p = r.d.data(); f.lo[0] = l0; f.lo[1] = l1; f.lo[2] = l2; f.n[0] = r.dims[0]; f.n[1] = r.dims[1]; f.n[2] = r.dims[2]; return f; }
F2 f2 (Rec& r) { F2 f; f.p = r.d.data(); f.lo[0] = 1; f.lo[1] = 1; f.n[0] = r.dims[0]; f.n[1] = r.dims[1]; return f; }
F1 f1 (Rec& r, int lo) { F1 f; f.p = r.d.data(); f.lo = lo; f.n = r.dims[0]; return f; }
int ival (Case& c, const char* k) { return c.at(k).i.at(0); }
double dval (Case& c, const char* k) { return c.at(k).d.at(0); }

void run_l1211 (Case& c, Stat& st, bool vb, std::string const& fn)
{
    const int ib = ival(c, "IBAR"), jb = ival(c, "JBAR"), kb = ival(c, "KBAR");
    Rec prhs = c.at("PRHS_IN");
    pres_compute_rhs_div(fds_interior(ib, jb, kb), f3(c.at("FVX"), 0, 0, 0), f3(c.at("FVY"), 0, 0, 0), f3(c.at("FVZ"), 0, 0, 0), f3(c.at("DDDT"), 0, 0, 0),
                         f1(c.at("RDX"), 0), f1(c.at("RDY"), 0), f1(c.at("RDZ"), 0), f3(prhs, 1, 1, 1));
    compare(st, "PRHS", prhs.d, c.at("PRHS_OUT").d, vb, fn);
}

void run_l1207 (Case& c, Stat& st, bool vb, std::string const& fn)
{
    const int ib = ival(c, "IBAR"), jb = ival(c, "JBAR"), kb = ival(c, "KBAR");
    Rec p = c.at("H"); std::fill(p.d.begin(), p.d.end(), -7.77E+77);
    pres_p_from_h(fds_with_ghosts(ib, jb, kb), f3(c.at("RHO"), -1, -1, -1), f3(c.at("H"), 0, 0, 0), f3(c.at("KRES"), 0, 0, 0), f3(p, 0, 0, 0));
    compare(st, "P", p.d, c.at("P_OUT").d, vb, fn);
}

void run_hbc (Case& c, Stat& st, bool vb, std::string const& fn)
{
    const int ib = ival(c, "IBAR"), jb = ival(c, "JBAR"), kb = ival(c, "KBAR");
    Rec h = c.at("H_IN");
    const Box3 b = fds_interior(ib, jb, kb);
    F3 hp = f3(h, 0, 0, 0);
    pres_h_bc_x(b, ival(c, "LBC"), dval(c, "DXI"), f2(c.at("BXS")), f2(c.at("BXF")), hp);
    pres_h_bc_y(b, ival(c, "MBC"), dval(c, "DETA"), f2(c.at("BYS")), f2(c.at("BYF")), hp);
    pres_h_bc_z(b, ival(c, "NBC"), dval(c, "DZETA"), f2(c.at("BZS")), f2(c.at("BZF")), hp);
    compare(st, "H", h.d, c.at("H_OUT").d, vb, fn);
}

void run_l1209 (Case& c, Stat& st, bool vb, std::string const& fn)
{
    const int ib = ival(c, "IBAR"), jb = ival(c, "JBAR"), kb = ival(c, "KBAR"), nw = ival(c, "NWALL");
    PoissonBcContext x;
    x.ibar = ib; x.jbar = jb; x.kbar = kb;
    // FDS allocates U/US in x as -1:IBP1, V/VS in y as -1:JBP1 and W/WS in z as -1:KBP1 (init.f90), and the dump hook writes the
    // arrays as shape-only records, so the dumped velocity arrays start at index -1 in their own direction. The synthetic cases
    // (ref_loops.f90.in) use 0:IBP1 for all three. Reading the FDS files with lower bound 0 shifts every velocity by one cell (the
    // cause of the former "open wall differs from FDS" finding).
    const bool fds_file = c.count("WV_RAMP") != 0 && !g_uu_lb0;   // uu_lb0=1: negative control, the former (wrong) lower bound 0
    x.hp = f3(c.at("HP"), 0, 0, 0); x.kres = f3(c.at("KRES"), 0, 0, 0);
    x.uu = f3(c.at("UU"), fds_file ? -1 : 0, 0, 0); x.vv = f3(c.at("VV"), 0, fds_file ? -1 : 0, 0); x.ww = f3(c.at("WW"), 0, 0, fds_file ? -1 : 0);
    x.fvx = f3(c.at("FVX"), 0, 0, 0); x.fvy = f3(c.at("FVY"), 0, 0, 0); x.fvz = f3(c.at("FVZ"), 0, 0, 0);
    x.hx = f1(c.at("HX"), 0); x.hy = f1(c.at("HY"), 0); x.hz = f1(c.at("HZ"), 0);
    x.dx = f1(c.at("DX"), 1); x.dy = f1(c.at("DY"), 1); x.dz = f1(c.at("DZ"), 1);
    x.rdxn = f1(c.at("RDXN"), 0); x.rdyn = f1(c.at("RDYN"), 0); x.rdzn = f1(c.at("RDZN"), 0);
    x.u_wind = f1(c.at("U_WIND"), 0); x.v_wind = f1(c.at("V_WIND"), 0); x.w_wind = f1(c.at("W_WIND"), 0);
    x.u0 = dval(c, "U0"); x.v0 = dval(c, "V0"); x.w0 = dval(c, "W0"); x.open_wind_boundary = ival(c, "OPEN_WIND") != 0;
    x.t = dval(c, "T"); x.dt = dval(c, "DT"); x.t_begin = dval(c, "T_BEGIN");
    // stand-in EVALUATE_RAMP of ref_loops.f90.in (the real one is the caller's)
    x.evaluate_ramp = [](double X, int ri) -> double { return ri < 1 ? 1.0 + 0.01 * X : (X * X + 0.75 * static_cast<double>(ri)) / (1.0 + X * X); };
    Rec bxs = c.at("BXS_IN"), bxf = c.at("BXF_IN"), bys = c.at("BYS_IN"), byf = c.at("BYF_IN"), bzs = c.at("BZS_IN"), bzf = c.at("BZF_IN");
    x.bxs = f2(bxs); x.bxf = f2(bxf); x.bys = f2(bys); x.byf = f2(byf); x.bzs = f2(bzs); x.bzf = f2(bzf);

    const bool fds = c.count("WV_RAMP") != 0;   // file written by the FDS dump hook (per-wall vent fields, FDS-evaluated ramp values)
    std::vector<PoissonVent> vents;
    std::vector<PoissonWall> walls(nw);
    if (!fds) {
        const int nv = static_cast<int>(c.at("V_RAMP").i.size());
        vents.resize(nv);
        Rec &ue = c.at("V_UE"), &ve = c.at("V_VE"), &we = c.at("V_WE");
        for (int v = 0; v < nv; ++v) {
            PoissonVent& t = vents[v];
            t.pressure_ramp_index = c.at("V_RAMP").i[v]; t.n_eddy = c.at("V_NEDDY").i[v]; t.ior = c.at("V_IOR").i[v]; t.dynamic_pressure = c.at("V_DYNP").d[v];
            auto slab = [&](Rec& r, int n0, int n1) { F2 f; f.p = r.d.data() + static_cast<std::size_t>(n0) * n1 * v; f.lo[0] = 1; f.lo[1] = 1; f.n[0] = n0; f.n[1] = n1; return f; };
            t.u_eddy = slab(ue, ue.dims[0], ue.dims[1]); t.v_eddy = slab(ve, ve.dims[0], ve.dims[1]); t.w_eddy = slab(we, we.dims[0], we.dims[1]);
        }
    } else {
        vents.resize(nw);
        for (int w = 0; w < nw; ++w) {
            PoissonVent& t = vents[w];
            t.pressure_ramp_index = c.at("WV_RAMP").i[w]; t.n_eddy = c.at("WV_NEDDY").i[w]; t.ior = c.at("WV_IOR").i[w]; t.dynamic_pressure = c.at("WV_DYNP").d[w];
        }
    }
    for (int w = 0; w < nw; ++w) {
        PoissonWall& q = walls[w];
        q.i = c.at("W_I").i[w]; q.j = c.at("W_J").i[w]; q.k = c.at("W_K").i[w]; q.ior = c.at("W_IOR").i[w];
        q.pressure_bc_type = c.at("W_PBC").i[w]; q.boundary_type = c.at("W_BT").i[w];
        q.dundt = c.at("W_DUNDT").d[w]; q.wall_work1 = c.at("WALL_WORK1").d[w];
        if (fds) {
            q.other_d = c.at("W_OTHERD").d[w];
            q.vent = &vents[w];
        } else {
            Rec &mdx = c.at("M_DX"), &mdy = c.at("M_DY"), &mdz = c.at("M_DZ");
            const int nom = c.at("W_NOM").i[w] - 1, a = std::abs(q.ior);
            const int nd = mdx.dims[0];
            q.other_d = a == 1 ? mdx.d[(c.at("W_IIO").i[w] - 1) + static_cast<std::size_t>(nd) * nom]
                      : a == 2 ? mdy.d[(c.at("W_JJO").i[w] - 1) + static_cast<std::size_t>(nd) * nom]
                               : mdz.d[(c.at("W_KKO").i[w] - 1) + static_cast<std::size_t>(nd) * nom];
            q.vent = &vents[c.at("W_VENT").i[w] - 1];
        }
        q.t_ign = c.at("W_TIGN").d[w]; q.rho_f = c.at("W_RHOF").d[w];
    }
    // FDS-derived file: EVALUATE_RAMP is FDS's own; the caller (here) returns the value FDS computed for the same call, in call order,
    // and checks that the time argument (TSI) the loop forms is bitwise the one FDS formed.
    std::size_t ncall = 0; long tsi_bad = 0;
    if (fds) {
        const std::vector<double>& rf = c.at("RAMPF_SEQ").d; const std::vector<double>& ts = c.at("TSI_SEQ").d;
        x.evaluate_ramp = [&, rf, ts](double tsi, int) -> double {
            if (ncall >= rf.size()) { ++tsi_bad; return 0.0; }
            std::uint64_t a, b; std::memcpy(&a, &tsi, 8); std::memcpy(&b, &ts[ncall], 8);
            if (a != b) ++tsi_bad;
            return rf[ncall++];
        };
    }
    pres_poisson_boundary_arrays(x, walls.data(), nw);
    if (fds) {
        if (tsi_bad || static_cast<int>(ncall) != ival(c, "N_RAMP")) {
            std::printf("  %s: ramp calls %zu (FDS made %d), TSI mismatches %ld\n", fn.c_str(), ncall, ival(c, "N_RAMP"), tsi_bad);
            st.nbad += tsi_bad ? tsi_bad : 1;
        }
    }
    compare(st, "BXS", bxs.d, c.at("BXS_OUT").d, vb, fn); compare(st, "BXF", bxf.d, c.at("BXF_OUT").d, vb, fn);
    compare(st, "BYS", bys.d, c.at("BYS_OUT").d, vb, fn); compare(st, "BYF", byf.d, c.at("BYF_OUT").d, vb, fn);
    compare(st, "BZS", bzs.d, c.at("BZS_OUT").d, vb, fn); compare(st, "BZF", bzf.d, c.at("BZF_OUT").d, vb, fn);
}
}   // namespace

int main (int argc, char** argv)
{
    std::string dir; bool vb = false;
    for (int a = 1; a < argc; ++a) {
        std::string s = argv[a];
        if (s.rfind("dir=", 0) == 0) dir = s.substr(4);
        else if (s == "uu_lb0=1") g_uu_lb0 = true;
        else if (s.rfind("verbose=", 0) == 0) { vb = s.substr(8) != "0"; g_debug = s.substr(8) == "2"; }
    }
    if (dir.empty()) { std::fprintf(stderr, "usage: pb_fds_loops dir=<case directory> [verbose=1]\n"); return 2; }
    std::vector<std::string> files;
    if (DIR* d = opendir(dir.c_str())) { while (dirent* e = readdir(d)) { std::string n = e->d_name; if (n.size() > 4 && n.substr(n.size() - 4) == ".bin") files.push_back(n); } closedir(d); }
    std::sort(files.begin(), files.end());
    if (files.empty()) { std::fprintf(stderr, "no case files in %s\n", dir.c_str()); return 2; }
    std::map<std::string, Stat> stats;
    int rc = 0;
    for (auto const& fn : files) {
        Case c;
        if (!read_case(dir + "/" + fn, c)) { std::printf("CASE %s: unreadable\n", fn.c_str()); rc = 1; continue; }
        const std::string kind = fn.substr(0, fn.find('_'));
        Stat& st = stats[kind]; st.cases++;
        const long before = st.nbad;
        try {
            if (kind == "L1211") run_l1211(c, st, vb, fn);
            else if (kind == "L1207") run_l1207(c, st, vb, fn);
            else if (kind == "HBC") run_hbc(c, st, vb, fn);
            else if (kind == "L1209") run_l1209(c, st, vb, fn);
            else { std::printf("CASE %s: unknown kind\n", fn.c_str()); rc = 1; }
        } catch (std::exception const& e) { std::printf("CASE %s: exception %s\n", fn.c_str(), e.what()); st.nbad++; }
        if (st.nbad != before) st.bad_cases++;
    }
    long tc = 0, tb = 0;
    for (auto const& kv : stats) {
        const char* loops = kv.first == "HBC" ? "L1220+L1221+L1222" : kv.first.c_str();
        std::printf("LOOP %-18s cases %3ld  elements %9ld  bitwise-different %8ld  max_abs %.3e  max_rel %.3e  %s\n", loops, kv.second.cases, kv.second.n, kv.second.nbad,
                    kv.second.max_abs, kv.second.max_rel, kv.second.nbad == 0 ? "BITWISE EQUAL" : "DIFFERENT");
        tc += kv.second.cases; tb += kv.second.nbad;
    }
    std::printf("SUMMARY cases=%ld bitwise_different_elements=%ld %s\n", tc, tb, (tb == 0 && rc == 0) ? "PASS" : "FAIL");
    return (tb == 0 && rc == 0) ? 0 : 1;
}
