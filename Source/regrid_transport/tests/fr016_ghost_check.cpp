// fr016_ghost_check.cpp: FR-016 baseline comparison for the scalar ghost cells at the coarse-fine interface of ns2d_16_int_1to2_refinement.
// Reads the reference dumps of the instrumented, unmodified FDS (A-09b style; Source/driver/tests/refdump), loads the valid cells of every mesh into
// two AMReX levels (12 coarse meshes + the hole under the fine mesh = level 0, mesh 13 = level 1), runs fill_covered_ghosts_fds and fill_fine_ghosts_fds, and compares
// the ghost cells FDS itself wrote (layers 1 and 2, face-adjacent) with ours: coarse-mesh ghost cells inside the hole, fine-mesh ghost cells outside the fine box.
// Usage: fr016_ghost_check <cases dir> <dump prefix> [step]      (dump files <prefix>, <prefix>.2 ... <prefix>.13, from `FDSREF_FILE=<prefix> mpirun -np 13 fds case.fds`)
// Not part of ctest (needs the dump files); exit code 0 = every non-conflict ghost cell agrees to rounding.
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include "LevelOps.H"
#include "RegridAmrCore.H"
#include "dump_reader.H"
#include "mesh_text.H"

using namespace fdsrt;

struct Stat { long n = 0, bad = 0; double maxrel = 0; };

int main(int argc, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    int rc = 0;
    {
        if (argc < 3) { std::fprintf(stderr, "usage: fr016_ghost_check <cases dir> <dump prefix> [step]\n"); return 2; }
        const std::string dir = argv[1], prefix = argv[2];
        const int want_step = argc > 3 ? std::atoi(argv[3]) : 3;
        const std::string text = fdsrt_test::read_file(dir + "/ns2d_16_int_1to2.fds");
        auto meshes = fdsrt_test::meshes_from_text(text);
        Report rep;
        AmrParams p = parse_amr_params(text, rep);
        Hierarchy h;
        if (!build_hierarchy_from_meshes(meshes, p, {true, false, true}, h, rep)) { std::fprintf(stderr, "hierarchy failed\n"); return 2; }
        RegridAmrCore core(h, p);
        struct Null : LevelListener { void make_level(const LevelLayout&) override {} void remake_level(const LevelLayout&) override {} void clear_level(int) override {} } nl;
        amrex::BoxList bl0;
        for (const GridBox& g : h.levels[0].grids) bl0.push_back(to_amrex(g.box));
        amrex::DistributionMapping dm0(amrex::BoxArray{bl0});
        core.init_static(nl, false, RegridAmrCore::DmFn(), &dm0);
        const int nm_total = static_cast<int>(meshes.size());          // 13: meshes 1..12 coarse, 13 fine
        const amrex::Box hole = to_amrex(h.levels[0].grids.back().box);
        const amrex::Box fbox = to_amrex(h.levels[1].grids[0].box);
        const amrex::IntVect ratio = core.refRatio(0);

        struct Case { const char* label; const char* rec; bool pred; std::vector<std::string> rho_zz; std::vector<std::string> extra; };
        // corrector-state ghosts (set by the last WALL_BC of the corrector): record DENS_P; predictor-state ghosts: DENS_C; viscosity ghosts: VISC_P
        const std::vector<Case> cases = {
            {"corrector RHO/ZZ/TMP/RSUM", "DENS_P", false, {"RHO", "ZZ"}, {"TMP", "RSUM"}},
            {"predictor RHOS/ZZS", "DENS_C", true, {"RHOS", "ZZS"}, {}},
            {"viscosity MU/KRES", "VISC_P", false, {}, {"MU", "KRES"}},
        };
        std::map<int, std::vector<fdsrt_test::DumpRec>> dumps;
        for (int m = 1; m <= nm_total; ++m) dumps[m] = fdsrt_test::read_dump(m == 1 ? prefix : prefix + "." + std::to_string(m));
        for (const Case& cs : cases) {
            std::map<int, const fdsrt_test::DumpRec*> rec;
            for (int m = 1; m <= nm_total; ++m)
                for (const auto& r : dumps[m])
                    if (r.name == cs.rec && r.icyc == want_step && (std::string(cs.rec) != "DENS_C" || true)) { rec[m] = &r; break; }
            if (static_cast<int>(rec.size()) != nm_total) { std::printf("SKIP %s: record %s step %d not found on all meshes\n", cs.label, cs.rec, want_step); continue; }
            const bool have_rz = !cs.rho_zz.empty();
            const std::string rn = have_rz ? cs.rho_zz[0] : "RHO", zn = have_rz ? cs.rho_zz[1] : "ZZ";
            // arrays on both levels, FDS-sized ghost widths (ZZ and TMP 2, RHO 3, RSUM/MU/KRES 1)
            struct Lev { amrex::MultiFab rho, zz, tmp, rsum, mu, kres; };
            auto mk = [&](int l) {
                const amrex::BoxArray& ba = core.boxArray(l);
                const amrex::DistributionMapping& dm = core.DistributionMap(l);
                return std::make_unique<Lev>(Lev{amrex::MultiFab(ba, dm, 1, 3), amrex::MultiFab(ba, dm, 1, 2), amrex::MultiFab(ba, dm, 1, 2), amrex::MultiFab(ba, dm, 1, 1),
                                                 amrex::MultiFab(ba, dm, 1, 1), amrex::MultiFab(ba, dm, 1, 1)});
            };
            auto L0 = mk(0), L1 = mk(1);
            for (Lev* L : {L0.get(), L1.get()}) for (auto* m : {&L->rho, &L->zz, &L->tmp, &L->rsum, &L->mu, &L->kres}) m->setVal(-1.0e30);   // hole interior and every ghost cell: poison
            // valid cells of mesh m -> AMReX box m-1 of its level
            auto load = [&](amrex::MultiFab& mf, int lev, int m, const char* name) {
                auto it = rec[m]->bef.find(name);
                if (it == rec[m]->bef.end()) return false;
                const auto& arr = it->second;
                const int bi = (m == nm_total) ? 0 : m - 1;
                const GridBox& gb = h.levels[lev].grids[bi];
                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    if (mfi.index() != bi) continue;
                    auto a = mf.array(mfi);
                    const amrex::Box vb = mfi.validbox();
                    for (int k = vb.smallEnd(2); k <= vb.bigEnd(2); ++k)
                        for (int j = vb.smallEnd(1); j <= vb.bigEnd(1); ++j)
                            for (int i = vb.smallEnd(0); i <= vb.bigEnd(0); ++i) a(i, j, k) = arr.at(i - gb.box.lo[0] + 1, j - gb.box.lo[1] + 1, k - gb.box.lo[2] + 1);
                }
                return true;
            };
            // the ghost values come from the neighbouring level's cells; mesh fields of a mesh that lives on another rank are not loaded here (this harness runs on one rank)
            for (int m = 1; m <= nm_total; ++m) {
                const int lev = (m == nm_total) ? 1 : 0;
                Lev& L = lev ? *L1 : *L0;
                if (have_rz) { load(L.rho, lev, m, rn.c_str()); load(L.zz, lev, m, zn.c_str()); }
                for (const auto& e : cs.extra) {
                    if (e == "TMP") load(L.tmp, lev, m, "TMP");
                    if (e == "RSUM") load(L.rsum, lev, m, "RSUM");
                    if (e == "MU") load(L.mu, lev, m, "MU");
                    if (e == "KRES") load(L.kres, lev, m, "KRES");
                }
            }
            // equation of state of the single-species case: RSUM and PBAR are constants of the run (valid-cell values of the dump)
            const double rsum0 = rec[nm_total]->bef.count("RSUM") ? rec[nm_total]->bef.at("RSUM").at(1, 1, 1) : 0.0;
            const double pbar0 = rec[nm_total]->bef.count("PBAR") ? rec[nm_total]->bef.at("PBAR").at(1, 1, 1) : 0.0;
            Thermo th;
            th.n_tracked = 1; th.rsum = [rsum0](const double*) { return rsum0; }; th.pbar = [pbar0](const amrex::IntVect&) { return pbar0; };
            const bool use_th = have_rz && rsum0 > 0 && pbar0 > 0 && !cs.pred;
            // the fine level's own same-level exchange is not needed (one fine box); the coarse ghosts across the hole are compared directly in the hole's cells
            ScalarStage cs0, fs0;
            if (have_rz) { cs0.rho = &L0->rho; cs0.zz = &L0->zz; fs0.rho = &L1->rho; fs0.zz = &L1->zz; }
            else { cs0.rho = &L0->mu; fs0.rho = &L1->mu; }   // viscosity-only case: a dummy "rho" slot carries MU; KRES goes through mean
            ScalarStage cst, fst;
            if (have_rz) {
                cst = cs0; fst = fs0;
                for (const auto& e : cs.extra) {
                    if (e == "TMP") { cst.tmp = &L0->tmp; fst.tmp = &L1->tmp; }
                    if (e == "RSUM") { cst.rsum = &L0->rsum; fst.rsum = &L1->rsum; }
                }
            } else {
                // MU and KRES are plain mean-fields; ScalarStage needs a rho slot: use a unit-density helper
                static amrex::MultiFab dummy0, dummy1;
                dummy0.define(core.boxArray(0), core.DistributionMap(0), 1, 3); dummy0.setVal(1.0);
                dummy1.define(core.boxArray(1), core.DistributionMap(1), 1, 3); dummy1.setVal(1.0);
                cst.rho = &dummy0; fst.rho = &dummy1;
                cst.mean = {&L0->mu, &L0->kres}; fst.mean = {&L1->mu, &L1->kres};
            }
            const CfStats sc = fill_covered_ghosts_fds(cst, fst, core.Geom(1), core.Geom(0), ratio, use_th ? &th : nullptr);
            const CfStats sf = fill_fine_ghosts_fds(fst, cst, core.Geom(1), core.Geom(0), ratio, use_th ? &th : nullptr);
            std::printf("%s (record %s, step %d): covered cells written %ld, fine ghost cells written %ld, conflict cells %ld%s\n", cs.label, cs.rec, want_step, sc.ghost_cells,
                        sf.ghost_cells, sc.conflicts, use_th ? "" : " (no EOS: TMP/RSUM not rebuilt)");

            // ---- compare ----
            std::map<std::string, Stat> stat;
            auto cmp = [&](const char* tag, double ours, double fds, bool conflict_zone) {
                const double rel = std::abs(ours - fds) / std::max(1e-300, std::abs(fds));
                Stat& s = stat[std::string(tag) + (conflict_zone ? " [hole edge/corner]" : "")];
                ++s.n; s.maxrel = std::max(s.maxrel, rel);
                if (rel > 1e-12) ++s.bad;
            };
            auto in_hole_corner_zone = [&](const amrex::IntVect& c) {   // within 2 layers of two faces of the hole: one covered cell serves two coarse meshes
                int near = 0;
                for (int d = 0; d < 3; ++d) { if (hole.length(d) <= 1) continue; if (c[d] - hole.smallEnd(d) <= 1 || hole.bigEnd(d) - c[d] <= 1) ++near; }
                return near >= 2;
            };
            for (int m = 1; m <= nm_total; ++m) {
                const int lev = (m == nm_total) ? 1 : 0;
                const int bi = (m == nm_total) ? 0 : m - 1;
                const amrex::Box vb = to_amrex(h.levels[lev].grids[bi].box);
                const IBox& gbx = h.levels[lev].grids[bi].box;
                Lev& L = lev ? *L1 : *L0;
                const amrex::Box lo_bnd = amrex::grow(vb, 2);
                const auto& r = *rec[m];
                for (int k = lo_bnd.smallEnd(2); k <= lo_bnd.bigEnd(2); ++k)
                    for (int j = lo_bnd.smallEnd(1); j <= lo_bnd.bigEnd(1); ++j)
                        for (int i = lo_bnd.smallEnd(0); i <= lo_bnd.bigEnd(0); ++i) {
                            const amrex::IntVect iv(i, j, k);
                            if (vb.contains(iv)) continue;
                            int nout = 0, layer = 0;
                            for (int d = 0; d < 3; ++d) { const int dd = std::max({vb.smallEnd(d) - iv[d], iv[d] - vb.bigEnd(d), 0}); if (dd > 0) { ++nout; layer = dd; } }
                            if (nout != 1) continue;
                            // coarse mesh: only ghost cells inside the hole; fine mesh: every ghost cell inside the domain (the fine box is interior)
                            if (lev == 0 && !hole.contains(iv)) continue;
                            if (lev == 1 && !core.Geom(1).Domain().contains(iv)) continue;
                            const int fi = i - gbx.lo[0] + 1, fj = j - gbx.lo[1] + 1, fk = k - gbx.lo[2] + 1;
                            const bool cz = (lev == 0) && in_hole_corner_zone(iv);
                            const std::string side = lev ? "fine" : "coarse";
                            const std::string ly = " layer " + std::to_string(layer);
                            // coarse side: the ghost cell of mesh m is a valid cell of the hole box; fine side: a ghost cell of the (only) fine box
                            auto get = [&](const amrex::MultiFab& mf) {
                                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
                                    if (lev == 0 ? mfi.validbox().contains(iv) : (mfi.index() == bi && mf.fabbox(mfi.index()).contains(iv))) return mf.const_array(mfi)(i, j, k);
                                return std::nan("");
                            };
                            auto one = [&](const char* nm, const amrex::MultiFab& mf, int maxlayer) {
                                if (layer > maxlayer || !r.bef.count(nm)) return;
                                const auto& arr = r.bef.at(nm);
                                if (!arr.has(fi, fj, fk)) return;
                                cmp((side + " " + nm + ly).c_str(), get(mf), arr.at(fi, fj, fk), cz);
                            };
                            if (have_rz) { one(rn.c_str(), L.rho, 2); one(zn.c_str(), L.zz, 2); }
                            for (const auto& e : cs.extra) {
                                if (e == "TMP") one("TMP", L.tmp, 2);
                                if (e == "RSUM") one("RSUM", L.rsum, 1);
                                if (e == "MU") one("MU", L.mu, 1);
                                if (e == "KRES") one("KRES", L.kres, 1);
                            }
                        }
            }
            for (auto& kv : stat) {
                std::printf("  %-44s cells %5ld  mismatches %5ld  max rel diff %.3g\n", kv.first.c_str(), kv.second.n, kv.second.bad, kv.second.maxrel);
                if (kv.first.find("[hole edge/corner]") == std::string::npos && kv.second.bad > 0) rc = 1;
            }
        }
    }
    amrex::Finalize();
    return rc;
}
