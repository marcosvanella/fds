// test_hierarchy.cpp: unit tests for the IR-002 grouping, input checks, static hierarchy, nesting checks and dump (R1). No AMReX.
// Run from ctest with the directory of the case files as argument. Exit code 0 = pass.
#include <cstdio>
#include <sstream>
#include <string>

#include "Hierarchy.H"
#include "mesh_text.H"

using namespace fdsrt;

static int g_checks = 0, g_fail = 0;
#define CHECK(c)                                                              \
    do {                                                                      \
        ++g_checks;                                                           \
        if (!(c)) {                                                           \
            ++g_fail;                                                         \
            std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); \
        }                                                                     \
    } while (0)

static void show(const Report& r)
{
    for (const auto& e : r.errors) std::fprintf(stderr, "    error: %s\n", e.c_str());
}

static MeshInput mesh(int i, int j, int k, double x0, double x1, double y0, double y1, double z0, double z1)
{
    MeshInput m;
    m.ijk[0] = i; m.ijk[1] = j; m.ijk[2] = k;
    m.xb[0] = x0; m.xb[1] = x1; m.xb[2] = y0; m.xb[3] = y1; m.xb[4] = z0; m.xb[5] = z1;
    return m;
}

static AmrParams amr(const std::string& line)
{
    Report r;
    AmrParams p = parse_amr_params(line, r);
    if (!r.ok()) { std::fprintf(stderr, "bad test input '%s'\n", line.c_str()); show(r); }
    return p;
}

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : "cases";
    const std::array<bool, 3> noper{false, false, false};

    // ---- IBox helpers ----
    {
        IBox a; a.lo = {0, 0, 0}; a.hi = {3, 3, 3};
        IBox b; b.lo = {1, 1, 1}; b.hi = {2, 2, 2};
        auto rest = subtract(a, b);
        long long n = 0;
        for (auto& x : rest) n += x.numPts();
        CHECK(n == 64 - 8);
        IBox g; g.lo = {-1, 0, 0}; g.hi = {1, 0, 0};
        IBox dom; dom.lo = {0, 0, 0}; dom.hi = {7, 0, 0};
        auto w = wrap_clip(g, dom, {true, false, false});  // sticks out by one cell on the low side: wraps to the high side
        CHECK(w.size() == 2);
        auto c = wrap_clip(g, dom, noper);                 // not periodic: cut off
        CHECK(c.size() == 1 && c[0].lo[0] == 0 && c[0].hi[0] == 1);
        CHECK(coarsen(IBox{{-3, 0, 0}, {-1, 0, 0}}, {2, 1, 1}).lo[0] == -2);
    }

    // ---- ns2d_16_int_1to2_refinement: 13 meshes, hidden y, default blocking factor 8 must pass ----
    const std::string ns2d_text = fdsrt_test::read_file(dir + "/ns2d_16_int_1to2.fds");
    auto ns2d = fdsrt_test::meshes_from_text(ns2d_text);
    CHECK(ns2d.size() == 13);
    {
        Report r;
        AmrParams p = parse_amr_params(ns2d_text, r);
        CHECK(r.ok() && p.present && p.max_level == 1);
        Hierarchy h;
        bool ok = build_hierarchy_from_meshes(ns2d, p, {true, false, true}, h, r);
        show(r);
        CHECK(ok && r.ok());
        CHECK(h.hidden[1] && !h.hidden[0]);
        CHECK(h.top_level() == 1);
        CHECK(h.levels[0].grids.size() == 13);  // 12 meshes + the 8x8 hole under the fine patch
        CHECK(h.levels[1].grids.size() == 1);
        std::ostringstream os;
        dump_hierarchy(os, h);
        // Level 0 lists the 12 meshes (MULT order) and then the 8x8 cover; the checks below pick the stable lines.
        CHECK(os.str().find("LEVEL 0 dx=0.392699082,0.1,0.392699082 ratio=1,1,1 domain=(0,0,0)-(15,0,15) bf=8,1,8 mgs=32,32,32 boxes=13") != std::string::npos);
        CHECK(os.str().find("LEVEL 1 dx=0.196349541,0.1,0.196349541 ratio=2,1,2 domain=(0,0,0)-(31,0,31) bf=8,1,8 mgs=32,32,32 boxes=1") != std::string::npos);
        CHECK(os.str().find("BOX 0 (8,0,8)-(23,0,23) mesh 13") != std::string::npos);
        CHECK(os.str().find("(4,0,4)-(11,0,11) added") != std::string::npos);
        // negative controls for the assertions of FR-013
        Report rc;
        check_nesting(h, rc);
        CHECK(rc.ok());
        Hierarchy bad = h;
        bad.levels[0].grids.pop_back();  // remove the cover: level 0 no longer tiles the domain
        Report r1;
        check_nesting(bad, r1);
        CHECK(r1.has_error_containing("do not tile the domain"));
        bad = h;
        bad.levels[1].grids.push_back(bad.levels[1].grids[0]);  // duplicate box
        Report r2;
        check_nesting(bad, r2);
        CHECK(r2.has_error_containing("overlap"));
        // refinable region: finer mesh footprint, tags outside discarded and counted once per level (IR-008)
        const IBox& tg = h.levels[0].taggable[0];
        CHECK(h.levels[0].taggable.size() == 1 && tg.lo[0] == 4 && tg.hi[0] == 11 && tg.lo[2] == 4 && tg.hi[2] == 11 && tg.hi[1] == 0);
        TagClipper tc(h);
        std::vector<IVec> tags{{5, 0, 5}, {0, 0, 0}, {11, 0, 11}, {12, 0, 4}};
        CHECK(tc.clip_tags(0, tags) == 2 && tags.size() == 2 && tc.discarded(0) == 2);
        CHECK(tc.first_discard(0) && !tc.first_discard(0) && !tc.first_discard(1));
    }
    {   // declared region adds to the refinable region
        Report r;
        AmrParams p = amr("&AMR MAX_LEVEL=1 /\n&AMR_REGION XB=0,3,-0.05,0.05,0,3 /");
        Hierarchy h;
        CHECK(build_hierarchy_from_meshes(ns2d, p, {true, false, true}, h, r));
        CHECK(h.levels[0].taggable.size() == 2 && h.levels[0].taggable[1].lo[0] == 0 && h.levels[0].taggable[1].hi[0] == 7);
    }
    {   // finer meshes without an &AMR line: error naming the missing line and the mesh pair; never inferred
        Report r;
        Hierarchy h;
        AmrParams none;
        CHECK(!build_hierarchy_from_meshes(ns2d, none, {true, false, true}, h, r));
        CHECK(r.has_error_containing("no &AMR line") && r.has_error_containing(" and mesh 13"));
        // same-size meshes without &AMR are fine (uniform mode)
        Report r2;
        std::vector<MeshInput> two{mesh(8, 8, 8, 0, 1, 0, 1, 0, 1), mesh(8, 8, 8, 1, 2, 0, 1, 0, 1)};
        CHECK(build_hierarchy_from_meshes(two, none, noper, h, r2) && h.top_level() == 0 && h.levels[0].grids.size() == 2);
    }
    {   // blocking factor 16: level-0 domain 16 divides, but the fine patch (cells 8..23) is not aligned
        Report r;
        Hierarchy h;
        AmrParams p = amr("&AMR MAX_LEVEL=1, BLOCKING_FACTOR=16 /");
        CHECK(!build_hierarchy_from_meshes(ns2d, p, {true, false, true}, h, r));
        CHECK(r.has_error_containing("mesh 13 on level 1 spans cells 8..23 in x, not aligned with BLOCKING_FACTOR 16 on level 1"));
    }

    // ---- race_test_1 remesh (A-35): level 0 = 34x18x32, five 4:1 interfaces ----
    auto rt1 = fdsrt_test::meshes_from_text(fdsrt_test::read_file(dir + "/race_test_1_r4.fds"));
    CHECK(rt1.size() == 6);
    {
        Report r;
        Hierarchy h;
        AmrParams p = amr("&AMR MAX_LEVEL=1, REF_RATIO=4 /");  // default blocking factor 8
        CHECK(!build_hierarchy_from_meshes(rt1, p, noper, h, r));
        CHECK(r.has_error_containing("level-0 domain is 34 cells in x, not divisible by BLOCKING_FACTOR 8"));
        CHECK(r.has_error_containing("level-0 domain is 18 cells in y, not divisible by BLOCKING_FACTOR 8"));
        CHECK(!r.has_error_containing("32 cells in z"));
        CHECK(r.warnings.size() == 1 && r.warnings[0].find("17x9x16") != std::string::npos);  // FR-010 coarsening warning
    }
    {
        Report r;
        Hierarchy h;
        AmrParams p = amr("&AMR MAX_LEVEL=1, REF_RATIO=4, BLOCKING_FACTOR=2,8 /");
        bool ok = build_hierarchy_from_meshes(rt1, p, noper, h, r);
        show(r);
        CHECK(ok);
        CHECK(h.levels[0].grids.size() == 6 && h.levels[1].grids.size() == 1);
        const IBox& b = h.levels[1].grids[0].box;
        CHECK(b.lo[0] == 56 && b.hi[0] == 79 && b.lo[1] == 24 && b.hi[1] == 47 && b.lo[2] == 0 && b.hi[2] == 15);
        CHECK(h.levels[0].domain.hi[0] == 33 && h.levels[0].domain.hi[1] == 17 && h.levels[0].domain.hi[2] == 31);
        CHECK(!r.warnings.empty() && r.warnings[0].find("17x9x16") != std::string::npos);
    }
    {   // ratio 2 hierarchy cannot hold a 4:1 mesh; and '2,8' with ratio 2 is rejected by the parser
        Report r;
        Hierarchy h;
        AmrParams p = amr("&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2,4 /");
        CHECK(!build_hierarchy_from_meshes(rt1, p, noper, h, r));
        CHECK(r.has_error_containing("not a level of the hierarchy"));
        Report rp;
        parse_amr_params("&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2,8 /", rp);
        CHECK(!rp.ok());
    }
    {   // three levels with ratio 2: the 4:1 mesh is a level-2 mesh and level 1 must be added around it for proper nesting
        Report r;
        Hierarchy h;
        AmrParams p = amr("&AMR MAX_LEVEL=2, REF_RATIO=2,2, BLOCKING_FACTOR=2,4,8 /");
        bool ok = build_hierarchy_from_meshes(rt1, p, noper, h, r);
        show(r);
        CHECK(ok && h.top_level() == 2);
        CHECK(h.levels[2].grids.size() == 1 && h.levels[1].grids.size() >= 1);
        bool added = false;
        for (auto& gb : h.levels[1].grids) added = added || gb.mesh < 0;
        CHECK(added);
        Report rc;
        check_nesting(h, rc);
        CHECK(rc.ok());
        Hierarchy bad = h;
        bad.levels[1].grids.pop_back();  // drop part of the nesting cover
        Report rb;
        check_nesting(bad, rb);
        CHECK(rb.has_error_containing("not properly nested in level 1"));
    }

    // ---- ratio checks on synthetic meshes (FR-010) ----
    {
        AmrParams p = amr("&AMR MAX_LEVEL=2, REF_RATIO=2,2 /");
        Hierarchy h;
        {   // 3:1 in all directions
            Report r;
            std::vector<MeshInput> m{mesh(3, 3, 3, 0, 3, 0, 3, 0, 3), mesh(3, 3, 3, 3, 4, 0, 1, 0, 1), mesh(3, 3, 3, 3, 4, 1, 3, 0, 3)};
            CHECK(!build_hierarchy_from_meshes(m, p, noper, h, r));
            CHECK(r.has_error_containing("refinement ratio 3 between mesh 1 and mesh 2 is not supported"));
        }
        {   // 5:1
            Report r;
            std::vector<MeshInput> m{mesh(5, 5, 5, 0, 5, 0, 5, 0, 5), mesh(5, 5, 5, 5, 6, 0, 1, 0, 1)};
            CHECK(!group_meshes(m, p, noper, r).ok);
            CHECK(r.has_error_containing("refinement ratio 5 between mesh 1 and mesh 2"));
        }
        {   // direction-dependent: 2:1 in x only
            Report r;
            std::vector<MeshInput> m{mesh(4, 4, 4, 0, 4, 0, 4, 0, 4), mesh(8, 4, 4, 4, 8, 0, 4, 0, 4)};
            CHECK(!group_meshes(m, p, noper, r).ok);
            CHECK(r.has_error_containing("between mesh 1 and mesh 2 depends on direction"));
        }
        {   // 8:1 neighbours: pair error even though 8 = 2*2*2 would need three levels
            Report r;
            std::vector<MeshInput> m{mesh(2, 2, 2, 0, 2, 0, 2, 0, 2), mesh(8, 8, 8, 2, 3, 0, 1, 0, 1), mesh(2, 2, 2, 2, 4, 0, 2, 0, 2)};
            CHECK(!group_meshes(m, p, noper, r).ok);
            CHECK(r.has_error_containing("refinement ratio 8 between mesh 1 and mesh 2"));
        }
        {   // fine mesh edge not on a level-0 face
            Report r;
            std::vector<MeshInput> m{mesh(4, 4, 4, 0, 4, 0, 4, 0, 4), mesh(3, 8, 8, 4, 5.5, 0, 4, 0, 4), mesh(2, 4, 4, 6, 8, 0, 4, 0, 4)};
            CHECK(!group_meshes(m, p, noper, r).ok);
            CHECK(r.has_error_containing("not on level 0 cell faces"));
        }
        {   // overlapping meshes
            Report r;
            std::vector<MeshInput> m{mesh(4, 4, 4, 0, 4, 0, 4, 0, 4), mesh(4, 4, 4, 3, 7, 0, 4, 0, 4)};
            CHECK(!group_meshes(m, p, noper, r).ok);
            CHECK(r.has_error_containing("mesh 1 and mesh 2 overlap"));
        }
        {   // periodic pair: a fine mesh at the low x end touches the coarse mesh at the high x end across the boundary (3:1)
            Report r;
            std::vector<MeshInput> m{mesh(3, 3, 3, 0, 1, 0, 1, 0, 1), mesh(1, 1, 1, 1, 2, 0, 1, 0, 1)};
            CHECK(!group_meshes(m, p, {true, false, false}, r).ok);
            CHECK(r.has_error_containing("refinement ratio 3 between mesh 1 and mesh 2"));
        }
        {   // gap cells (non-box level 0) are reported, not silently accepted
            Report r;
            std::vector<MeshInput> m{mesh(4, 4, 4, 0, 4, 0, 4, 0, 4), mesh(4, 4, 4, 8, 12, 0, 4, 0, 4)};
            CHECK(!group_meshes(m, p, noper, r).ok);
            CHECK(r.has_error_containing("not a box"));
        }
    }

    std::printf("test_hierarchy: %d checks, %d failures\n", g_checks, g_fail);
    return g_fail == 0 ? 0 : 1;
}
