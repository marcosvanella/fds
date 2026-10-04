// test_amr_input.cpp: unit tests for the &AMR parser (R0, IR-003, FR-010). Pure C++, no AMReX. Exit code 0 = pass.
#include <cstdio>
#include <string>

#include "AmrInput.H"

using namespace fdsrt;

static int g_checks = 0, g_fail = 0;
#define CHECK(c)                                                                  \
    do {                                                                          \
        ++g_checks;                                                               \
        if (!(c)) {                                                               \
            ++g_fail;                                                             \
            std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c);     \
        }                                                                         \
    } while (0)

static AmrParams parse(const std::string& t, Report& r) { return parse_amr_params(t, r); }

int main()
{
    {  // no &AMR: defaults, not present, no errors
        Report r;
        AmrParams p = parse("&HEAD CHID='x'/\n&MESH IJK=4,4,4, XB=0,1,0,1,0,1/\n", r);
        CHECK(r.ok());
        CHECK(!p.present);
        CHECK(p.max_level == 0);
    }
    {  // full line, FDS style (comments outside groups, &TAIL stops scanning)
        Report r;
        AmrParams p = parse(
            "this is a comment & with an ampersand\n"
            "&AMR MAX_LEVEL=2, REF_RATIO=2,4, REGRID_INTERVAL=8, BLOCKING_FACTOR=4,8,16, MAX_GRID_SIZE=32,64,64,\n"
            "     N_ERROR_BUF=2, N_PROPER=1, GRID_EFF=0.8, OUTPUT_LEVEL_CAP=1, VELOCITY_TRANSFER='FACE_LINEAR' /\n"
            "&AMR_REGION XB=0,1,0,1,0,2, LEVEL=2 /\n&AMR_REGION XB=0.5,1,0,1,0,1 /\n&TAIL /\n&AMR MAX_LEVEL=9 /\n", r);
        CHECK(r.ok());
        CHECK(p.present && p.max_level == 2);
        CHECK(p.ref_ratio_at(0) == 2 && p.ref_ratio_at(1) == 4);
        CHECK(p.regrid_interval == 8 && p.n_error_buf == 2 && p.grid_eff == 0.8 && p.output_level() == 1);
        CHECK(p.blocking_factor_at(2) == 16 && p.max_grid_size_at(0) == 32);
        CHECK(p.velocity_transfer == VelocityTransfer::FaceLinear);
        CHECK(p.regions.size() == 2 && p.regions[0].level == 2 && p.regions[1].level == -1 && p.regions[1].xb[0] == 0.5);
    }
    {  // defaults when only MAX_LEVEL is given: ratio 2, FaceDivFree (provisional default, FR-012)
        Report r;
        AmrParams p = parse("&AMR MAX_LEVEL=1 /", r);
        CHECK(r.ok() && p.ref_ratio_at(0) == 2 && p.velocity_transfer == VelocityTransfer::FaceDivFree && p.output_level() == 1);
    }
    {  // refinement ratios other than 2 and 4 are rejected (FR-010)
        for (int rr : {1, 3, 5, 8}) {
            Report r;
            parse("&AMR MAX_LEVEL=1, REF_RATIO=" + std::to_string(rr) + " /", r);
            CHECK(r.has_error_containing("supports refinement ratios 2 and 4"));
        }
        Report r;
        parse("&AMR MAX_LEVEL=2, REF_RATIO=2,3 /", r);
        CHECK(r.has_error_containing("REF_RATIO entry 2 is 3"));
    }
    {  // blocking factor rules (A-35 caveat): '2,8' needs ratio 4; ratio 2 needs '2,4'
        Report r1, r2, r3;
        parse("&AMR MAX_LEVEL=1, REF_RATIO=4, BLOCKING_FACTOR=2,8 /", r1);
        CHECK(r1.ok());
        parse("&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2,8 /", r2);
        CHECK(r2.has_error_containing("does not allow blocking factor 8 on level 1"));
        parse("&AMR MAX_LEVEL=1, REF_RATIO=2, BLOCKING_FACTOR=2,4 /", r3);
        CHECK(r3.ok());
        Report r4;
        parse("&AMR MAX_LEVEL=1, BLOCKING_FACTOR=6 /", r4);
        CHECK(r4.has_error_containing("power of 2"));
    }
    {  // list lengths, ranges, unknown names, duplicates
        Report r;
        parse("&AMR MAX_LEVEL=2, REF_RATIO=2,2,2 /", r);
        CHECK(r.has_error_containing("REF_RATIO needs one value or MAX_LEVEL"));
        Report r2;
        parse("&AMR MAX_LEVEL=1, BLOCKING_FACTOR=8,8,8 /", r2);
        CHECK(r2.has_error_containing("BLOCKING_FACTOR needs one value or MAX_LEVEL+1"));
        Report r3;
        parse("&AMR MAX_LEVEL=1, MAX_GRID_SIZE=12 /", r3);
        CHECK(r3.has_error_containing("not a multiple of the blocking factor"));
        Report r4;
        parse("&AMR MAX_LEVEL=1, FOO=3 /", r4);
        CHECK(r4.has_error_containing("unknown parameter FOO"));
        Report r5;
        parse("&AMR MAX_LEVEL=1 /\n&AMR MAX_LEVEL=2 /", r5);
        CHECK(r5.has_error_containing("more than one &AMR"));
        Report r6;
        parse("&AMR MAX_LEVEL=-1 /", r6);
        CHECK(r6.has_error_containing("MAX_LEVEL must be >= 0"));
        Report r7;
        parse("&AMR MAX_LEVEL=1, GRID_EFF=1.5, N_PROPER=0, REGRID_INTERVAL=-1, N_ERROR_BUF=-1, OUTPUT_LEVEL_CAP=3 /", r7);
        CHECK(r7.errors.size() == 5);
        Report r8;
        parse("&AMR MAX_LEVEL=1, VELOCITY_TRANSFER='NOPE' /", r8);
        CHECK(r8.has_error_containing("VELOCITY_TRANSFER"));
        {   // POST_REGRID_PROJECTION (D-063): default AUTO, ON/OFF/logicals accepted, anything else is an error
            Report ra, rb, rc, rd, re;
            CHECK(parse("&AMR MAX_LEVEL=1 /", ra).post_regrid_projection == PostRegridProjection::Auto);
            CHECK(parse("&AMR MAX_LEVEL=1, POST_REGRID_PROJECTION='ON' /", rb).post_regrid_projection == PostRegridProjection::On);
            CHECK(parse("&AMR MAX_LEVEL=1, POST_REGRID_PROJECTION='off' /", rc).post_regrid_projection == PostRegridProjection::Off);
            CHECK(parse("&AMR MAX_LEVEL=1, POST_REGRID_PROJECTION=.FALSE. /", rd).post_regrid_projection == PostRegridProjection::Off);
            parse("&AMR MAX_LEVEL=1, POST_REGRID_PROJECTION='SOMETIMES' /", re);
            CHECK(re.has_error_containing("POST_REGRID_PROJECTION"));
            CHECK(!ra.has_error_containing("POST_REGRID") && !rb.has_error_containing("POST_REGRID"));
        }
        Report r9;
        parse("&AMR MAX_LEVEL=1, MAX_LEVEL=abc /", r9);
        CHECK(r9.has_error_containing("MAX_LEVEL must be one integer"));
    }
    {  // regions
        Report r;
        parse("&AMR MAX_LEVEL=1 /\n&AMR_REGION XB=0,1,0,1,1,0 /", r);
        CHECK(r.has_error_containing("zero or negative extent in direction 3"));
        Report r2;
        parse("&AMR MAX_LEVEL=1 /\n&AMR_REGION XB=0,1,0,1 /", r2);
        CHECK(r2.has_error_containing("XB must be six numbers"));
        Report r3;
        parse("&AMR MAX_LEVEL=1 /\n&AMR_REGION XB=0,1,0,1,0,1, LEVEL=2 /", r3);
        CHECK(r3.has_error_containing("LEVEL must be between 1 and MAX_LEVEL"));
        Report r4;
        parse("&AMR_REGION XB=0,1,0,1,0,1 /", r4);
        CHECK(r4.has_error_containing("without an &AMR line"));
    }
    std::printf("test_amr_input: %d checks, %d failures\n", g_checks, g_fail);
    return g_fail == 0 ? 0 : 1;
}
