// test_input_converter.cpp: the D-076 input converter (InputConverter.H): classification, level-0-only text, hierarchy, round trip, error cases with negative controls.
// No AMReX, no FDS. Argument 1: directory of the case files. Exit code 0 = pass.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

#include "InputConverter.H"
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

static bool has(const std::string& s, const std::string& sub) { return s.find(sub) != std::string::npos; }
static int count_char(const std::string& s, char c) { int n = 0; for (char x : s) n += x == c; return n; }
static void show(const Report& r) { for (const auto& e : r.errors) std::fprintf(stderr, "    error: %s\n", e.c_str()); }

// Small 2-D example: domain [0,2] x [0,2], level-0 cells 0.25 (8 x 8), a ring of four level-0 meshes around the hole x,z in [0.5,1.5] that a ratio-2 mesh fills.
// `fine_x0` shifts the fine mesh (0.5 = aligned), `amr` is the &AMR line, `ranks` adds MPI_PROCESS, `extra` is appended before &TAIL.
static std::string example(double fine_x0 = 0.5, const std::string& amr = "&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /", bool ranks = false,
                           const std::string& extra = "", const std::string& fine_id = "", const std::string& fine_ijk = "8,1,8", double fine_w = 1.0)
{
    std::ostringstream s;
    s << "&HEAD CHID='conv_small' /\n"
      << "! &MESH IJK=8,1,8, XB=0,2,-0.05,0.05,0,2 /   (single-mesh form, commented out)\n"
      << "&MESH IJK=2,1,4, XB=0.0,0.5,-0.05,0.05,0.5,1.5" << (ranks ? ", MPI_PROCESS=0" : "") << " /  left\n"
      << "&MESH IJK=" << fine_ijk << ", XB=" << fine_x0 << "," << fine_x0 + fine_w << ",-0.05,0.05,0.5,1.5" << (ranks ? ", MPI_PROCESS=1" : "") << (fine_id.empty() ? "" : ", ID='" + fine_id + "'")
      << " /  fine, ratio 2\n"
      << "&MESH IJK=2,1,4, XB=1.5,2.0,-0.05,0.05,0.5,1.5" << (ranks ? ", MPI_PROCESS=2" : "") << " /  right\n"
      << "&MESH IJK=8,1,2, XB=0.0,2.0,-0.05,0.05,0.0,0.5" << (ranks ? ", MPI_PROCESS=3" : "") << " /  bottom\n"
      << "&MESH IJK=8,1,2, XB=0.0,2.0,-0.05,0.05,1.5,2.0" << (ranks ? ", MPI_PROCESS=4" : "") << " /  top\n"
      << amr << "\n"
      << "&TIME T_END=1.0 /\n"
      << "&VENT DB='XMIN', SURF_ID='PERIODIC' /\n&VENT DB='XMAX', SURF_ID='PERIODIC' /\n&VENT DB='ZMIN', SURF_ID='PERIODIC' /\n&VENT DB='ZMAX', SURF_ID='PERIODIC' /\n"
      << extra << "&TAIL /\n";
    return s.str();
}

static bool same_box(const IBox& a, const IBox& b)
{
    for (int d = 0; d < 3; ++d) if (a.lo[d] != b.lo[d] || a.hi[d] != b.hi[d]) return false;
    return true;
}

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : "cases";

    // ---- namelist spans and mesh expansion
    {
        Report r;
        auto sp = find_group_spans("&HEAD CHID='a/b' /\n! &MESH IJK=1,1,1 /\n  &MESH IJK=4,4,4 /  tail text\n&MISC X=1\n/\n&TAIL /\n&MESH IJK=9,9,9 /\n", r);
        CHECK(r.ok() && sp.size() == 3);
        CHECK(sp.size() == 3 && sp[0].name == "HEAD" && sp[1].name == "MESH" && sp[1].line == 3 && sp[2].name == "MISC" && sp[2].line == 4);   // quote with '/', commented group, &TAIL stop
        Report r2;
        find_group_spans("&MESH IJK=4,4,4\n", r2);
        CHECK(r2.has_error_containing("no closing '/'"));   // unterminated group
    }
    {
        // the 13-mesh input of the FR-016 gate: 16 MULT copies minus the 2 x 2 skip box = 12, plus the fine mesh; mesh order as READ_MESH (I fastest)
        const std::string txt = fdsrt_test::read_file(dir + "/ns2d_16_int_1to2.fds");
        Report r;
        auto sp = find_group_spans(txt, r);
        std::vector<MeshLine> lines;
        std::vector<MeshInput> m;
        CHECK(parse_meshes(txt, sp, lines, m, r) && m.size() == 13 && lines.size() == 2 && lines[0].count == 12 && lines[1].first == 12);
        auto ref = fdsrt_test::meshes_from_text(txt);   // the test helper expands the same way
        bool same = ref.size() == m.size();
        for (size_t i = 0; same && i < m.size(); ++i)
            for (int q = 0; q < 6; ++q) same = same && std::fabs(m[i].xb[q] - ref[i].xb[q]) < 1e-14;
        CHECK(same);
        CHECK(m.size() == 13 && std::fabs(m[1].xb[0] - 1.5707963267949) < 1e-12 && std::fabs(m[4].xb[4] - 1.5707963267949) < 1e-12);   // copy 1: I+1; copy 4: K+1 after the I row of 4 ... (I fastest)
        // MULT_ID that does not exist
        Report r3;
        std::vector<MeshLine> l3; std::vector<MeshInput> m3;
        std::string bad = "&MESH IJK=4,1,4, XB=0,1,-0.05,0.05,0,1, MULT_ID='nope' /\n";
        auto sp3 = find_group_spans(bad, r3);
        CHECK(!parse_meshes(bad, sp3, l3, m3, r3) && r3.has_error_containing("MULT_ID 'nope' not found"));
    }

    // ---- the gate input: ns2d_16_int_1to2_refinement (mesh lines + &AMR line), periodic x and z
    {
        const std::string txt = fdsrt_test::read_file(dir + "/ns2d_16_int_1to2.fds") + "&VENT DB='XMIN', SURF_ID='PERIODIC' /\n&VENT DB='XMAX', SURF_ID='PERIODIC' /\n&VENT DB='ZMIN', SURF_ID='PERIODIC' /\n&VENT DB='ZMAX', SURF_ID='PERIODIC' /\n&TAIL /\n";
        Report r;
        ConvertResult c;
        const bool ok = convert_input(txt, c, r);
        show(r);
        CHECK(ok && r.ok() && c.ok);
        CHECK(c.periodic[0] && !c.periodic[1] && c.periodic[2]);
        CHECK(c.meshes.size() == 13 && c.level0_meshes.size() == 12 && c.removed_meshes.size() == 1 && c.removed_meshes[0] == 12);
        CHECK(c.n_removed_lines == 1 && c.cover_boxes.size() == 1);
        CHECK(c.hierarchy.top == 1 && c.hierarchy.levels[1].grids.size() == 1 && c.hierarchy.levels[0].grids.size() == 13);
        CHECK(!has(c.level0_text, "&AMR") && !has(c.level0_text, "IJK=16,1,16") && has(c.level0_text, "PERIODIC") && has(c.level0_text, "MULT_ID='m1'"));
        CHECK(has(c.level0_text, "IJK=8,1,8"));   // the cover: 8 x 8 level-0 cells under the fine mesh
        // the converted input read again: 13 level-0 meshes, one level, same level-0 boxes as the hierarchy, cover XB = XB of the removed fine mesh
        Report r2;
        ConvertResult c2;
        const bool ok2 = convert_input(c.level0_text, c2, r2);
        show(r2);
        CHECK(ok2 && c2.meshes.size() == 13 && c2.removed_meshes.empty() && c2.hierarchy.top == 0 && c2.cover_boxes.empty());
        bool boxes = ok2 && c2.hierarchy.levels[0].grids.size() == c.hierarchy.levels[0].grids.size();
        for (size_t i = 0; boxes && i < c2.hierarchy.levels[0].grids.size(); ++i) boxes = same_box(c2.hierarchy.levels[0].grids[i].box, c.hierarchy.levels[0].grids[i].box);
        CHECK(boxes);
        bool xb = ok2 && c2.meshes.size() == 13;
        for (int q = 0; xb && q < 6; ++q) xb = std::fabs(c2.meshes[12].xb[q] - c.meshes[12].xb[q]) < 1e-12;   // converted mesh 13 = cover = footprint of the removed mesh
        CHECK(xb);
        CHECK(ok2 && c2.level0_text == c.level0_text);   // idempotent: a level-0-only input converts to itself
        // everything that is not a mesh or AMR line is unchanged: same number of lines plus one cover line
        CHECK(count_char(c.level0_text, '\n') == count_char(txt, '\n') + 1);
    }

    // ---- small multi-mesh example: round trip, ranks, hierarchy
    {
        Report r;
        ConvertResult c;
        const bool ok = convert_input(example(), c, r);
        show(r);
        CHECK(ok && c.meshes.size() == 5 && c.level0_meshes.size() == 4 && c.removed_meshes.size() == 1 && c.removed_meshes[0] == 1);   // the commented-out single mesh is not read
        CHECK(c.hierarchy.top == 1 && c.hierarchy.levels[1].grids.size() == 1 && c.cover_boxes.size() == 1);
        CHECK(c.grouping.level_of_mesh[1] == 1 && c.grouping.level_of_mesh[0] == 0 && c.grouping.level_of_mesh[4] == 0);
        CHECK(same_box(c.cover_boxes[0], IBox{{2, 0, 2}, {5, 0, 5}}));   // coarse cells 2..5 in x and z
        CHECK(!has(c.level0_text, "IJK=8,1,8, XB=0.5") && has(c.level0_text, "! converter: finer MESH line removed (level 1, input mesh 2)"));
        CHECK(!has(c.level0_text, "fine, ratio 2"));                    // the trailing comment of the removed line goes with it
        CHECK(has(c.level0_text, "left") && has(c.level0_text, "bottom") && has(c.level0_text, "single-mesh form"));
        Report r2;
        ConvertResult c2;
        CHECK(convert_input(c.level0_text, c2, r2) && c2.meshes.size() == 5 && c2.hierarchy.top == 0 && c2.level0_text == c.level0_text);
        // meshes of the converted input: level-0 meshes in original order, then the cover
        bool order = c2.meshes.size() == 5;
        const int idx[4] = {0, 2, 3, 4};
        for (int k = 0; order && k < 4; ++k)
            for (int q = 0; q < 6; ++q) order = order && std::fabs(c2.meshes[k].xb[q] - c.meshes[idx[k]].xb[q]) < 1e-14;
        for (int q = 0; order && q < 6; ++q) order = order && std::fabs(c2.meshes[4].xb[q] - c.meshes[1].xb[q]) < 1e-14;
        CHECK(order);
        // with explicit ranks that stay continuous after the removal (fine mesh last in rank order): ok; the cover takes the rank of the last level-0 mesh
        Report r3;
        ConvertResult c3;
        std::string ex = example(0.5, "&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /", true);
        CHECK(!convert_input(ex, c3, r3) && r3.has_error_containing("MPI_PROCESS is not continuous"));   // fine mesh had rank 1 between 0 and 2
    }

    // ---- error cases, each with its positive control above (aligned, &AMR line present, ratio 2, ids not referenced, ranks not given)
    {   // misaligned fine mesh: its high x edge at 1.375 lies in the middle of a level-0 cell (7 fine cells wide); the aligned control is the base example
        Report r;
        ConvertResult c;
        CHECK(!convert_input(example(0.5, "&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /", false, "", "", "7,1,8", 0.875), c, r));
        CHECK(r.has_error_containing("not on level 0 cell faces") && r.has_error_containing("mesh 2"));
    }
    {   // finer mesh without an &AMR line: error names the line to add and a mesh pair
        Report r;
        ConvertResult c;
        CHECK(!convert_input(example(0.5, ""), c, r));
        CHECK(r.has_error_containing("no &AMR line") && r.has_error_containing("mesh 1 and mesh 2") && r.has_error_containing("never inferred"));
        CHECK(c.level0_text.empty());
    }
    {   // ratio 3
        Report r;
        ConvertResult c;
        CHECK(!convert_input(example(0.5, "&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /", false, "", "", "12,1,12"), c, r));
        CHECK(r.has_error_containing("ratio 3") && r.has_error_containing("not supported"));
    }
    {   // overlapping meshes (a level-0 mesh over the fine mesh)
        Report r;
        ConvertResult c;
        CHECK(!convert_input(example(0.5, "&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /", false, "&MESH IJK=2,1,2, XB=0.75,1.25,-0.05,0.05,0.75,1.25 /\n"), c, r));
        CHECK(!r.ok());
    }
    {   // a DEVC-like group that names the removed fine mesh by MESH_ID is an error; control: the same fine mesh with an ID that nothing refers to converts
        const std::string amr = "&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /";
        Report r;
        ConvertResult c;
        CHECK(!convert_input(example(0.5, amr, false, "&DEVC ID='d', QUANTITY='DENSITY', XYZ=1,0,1, MESH_ID='FINE' /\n", "FINE"), c, r));
        CHECK(r.has_error_containing("MESH_ID='FINE'") && r.has_error_containing("removes"));
        Report r2;
        ConvertResult c2;
        CHECK(convert_input(example(0.5, amr, false, "", "FINE"), c2, r2) && c2.removed_meshes.size() == 1);
    }
    {   // no finer mesh at all, no &AMR line: converts to itself (control for the "no &AMR line" error)
        Report r;
        ConvertResult c;
        const std::string t = "&HEAD CHID='one' /\n&MESH IJK=8,8,8, XB=0,1,0,1,0,1 /\n&TAIL /\n";
        CHECK(convert_input(t, c, r) && c.level0_text == t && c.hierarchy.top == 0 && c.removed_meshes.empty());
    }
    {   // &AMR group is removed from the text even when no mesh is finer; an unknown &AMR name is still an error with its line
        Report r;
        ConvertResult c;
        CHECK(convert_input("&MESH IJK=8,8,8, XB=0,1,0,1,0,1 /\n&AMR MAX_LEVEL=0 /\n&TAIL /\n", c, r) && !has(c.level0_text, "&AMR"));
        Report r2;
        ConvertResult c2;
        CHECK(!convert_input("&MESH IJK=8,8,8, XB=0,1,0,1,0,1 /\n\n&AMR MAX_LVL=1 /\n&TAIL /\n", c2, r2) && r2.has_error_containing("line 3"));
    }

    {   // race_test_1_r4 (thread-check input with mesh 3 remeshed to a 4:1 pair): must convert to a true equal-level-0 input. The case file has no &AMR line (AMR mode is never
        // inferred), so the test adds one. The FDS-side guard (upstream patch 0010, ERROR 9001 on any level-0 NIC > 1 face) is NOT committed or validated yet: its check is
        // PENDING (see notes/input-converter.md); the converter-side assertion below is the check in force.
        const std::string full = fdsrt_test::read_file(dir + "/race_test_1_r4_full.fds");
        const std::string amr_line = "&AMR MAX_LEVEL=1, REF_RATIO=4, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /\n";
        const size_t at = full.find("&TIME");
        CHECK(!full.empty() && at != std::string::npos);
        const std::string txt = full.substr(0, at) + amr_line + full.substr(at);
        Report r;
        ConvertResult c;
        const bool ok = convert_input(txt, c, r);
        if (!ok) show(r);
        CHECK(ok && r.ok() && c.meshes.size() == 6 && c.removed_meshes.size() == 1 && c.removed_meshes[0] == 2 && c.level0_meshes.size() == 5);
        CHECK(ok && c.cover_boxes.size() == 1 && c.hierarchy.top == 1);
        CHECK(ok && c.grouping.n0[0] == 34 && c.grouping.n0[1] == 18 && c.grouping.n0[2] == 32);
        // re-read the converted text as FDS would: six meshes, all at the level-0 cell size, tiling 34 x 18 x 32 cells, no ratio other than 1 anywhere
        Report r2;
        std::vector<GroupSpan> sp = find_group_spans(c.level0_text, r2);
        std::vector<MeshLine> ml;
        std::vector<MeshInput> m2;
        CHECK(parse_meshes(c.level0_text, sp, ml, m2, r2) && m2.size() == 6);
        long long cells = 0;
        double worst = 0;
        for (const MeshInput& mi : m2) {
            cells += static_cast<long long>(mi.ijk[0]) * mi.ijk[1] * mi.ijk[2];
            for (int d = 0; d < 3; ++d) worst = std::max(worst, std::fabs(std::fabs(mi.xb[2 * d + 1] - mi.xb[2 * d]) / mi.ijk[d] - 0.05));
        }
        CHECK(cells == 34LL * 18 * 32 && worst < 1e-12);
        CHECK(verify_equal_level0(m2, c.grouping, r2) && r2.ok());
        // the added cover is the footprint of mesh 3: x,y in [-0.15,0.15], z in [0,0.2], 6 x 6 x 4 level-0 cells
        const IBox& cb = c.cover_boxes[0];
        CHECK(cb.hi[0] - cb.lo[0] + 1 == 6 && cb.hi[1] - cb.lo[1] + 1 == 6 && cb.hi[2] - cb.lo[2] + 1 == 4 && m2.size() == 6 &&
              std::fabs(m2[5].xb[0] + 0.15) < 1e-12 && std::fabs(m2[5].xb[1] - 0.15) < 1e-12 && std::fabs(m2[5].xb[5] - 0.20) < 1e-12);
        // hierarchy level 0 = the five kept meshes and the cover, the finer mesh is level 1 (24 x 24 x 16 cells at ratio 4)
        CHECK(c.hierarchy.levels[0].grids.size() == 6 && c.hierarchy.levels[1].grids.size() == 1);
        // the five 4:1 interfaces of the input: mesh 3 against meshes 1, 2, 4, 5, 6 (each pair touches one face), none left in the converted text
        int fine_faces = 0;
        for (int m : {0, 1, 3, 4, 5}) {
            int touch = 0;
            for (int d = 0; d < 3; ++d) {
                const double a0 = std::min(c.meshes[2].xb[2 * d], c.meshes[2].xb[2 * d + 1]), a1 = std::max(c.meshes[2].xb[2 * d], c.meshes[2].xb[2 * d + 1]);
                const double b0 = std::min(c.meshes[m].xb[2 * d], c.meshes[m].xb[2 * d + 1]), b1 = std::max(c.meshes[m].xb[2 * d], c.meshes[m].xb[2 * d + 1]);
                if (std::fabs(a1 - b0) < 1e-9 || std::fabs(b1 - a0) < 1e-9) ++touch;
            }
            fine_faces += touch == 1;
        }
        CHECK(fine_faces == 5);
        // negative controls: the unconverted meshes are not equal-level-0 (mesh 3 has cell size 0.0125); removing the cover leaves a hole; a duplicated mesh overlaps
        Report n1;
        CHECK(!verify_equal_level0(c.meshes, c.grouping, n1) && n1.has_error_containing("mesh 3 has cell size") && n1.has_error_containing("NIC > 1"));
        Report n2;
        std::vector<MeshInput> hole(m2.begin(), m2.begin() + 5);
        CHECK(!verify_equal_level0(hole, c.grouping, n2) && n2.has_error_containing("do not tile"));
        Report n3;
        std::vector<MeshInput> dup = m2;
        dup.push_back(m2[0]);
        CHECK(!verify_equal_level0(dup, c.grouping, n3) && n3.has_error_containing("overlap"));
        Report n4;
        ConvertResult c4;
        CHECK(!convert_input(full, c4, n4) && n4.has_error_containing("&AMR"));   // without the &AMR line the converter refuses (not inferred)
        // a blocking factor that does not fit mesh 3 (24 x 24 x 16 cells) is refused, so the input is never converted half way
        Report n5;
        ConvertResult c5;
        CHECK(!convert_input(full.substr(0, at) + "&AMR MAX_LEVEL=1, REF_RATIO=4, BLOCKING_FACTOR=32, MAX_GRID_SIZE=64 /\n" + full.substr(at), c5, n5));
        std::printf("PENDING: patch 0010 guard check (ERROR 9001 on race_test_1_r4 unconverted): patch not committed or validated; converter-side equal-level-0 assertion in force\n");
    }

    std::printf("test_input_converter: %d checks, %d failures\n", g_checks, g_fail);
    return g_fail == 0 ? 0 : 1;
}
