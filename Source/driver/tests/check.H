// check.H: test macros only, no kernel code. Kernel-facing rules (M2a): passive scalars are handled by Fields.cpp; only uniform Cartesian metrics are used.
// check.H: minimal check counter shared by the driver tests (no framework dependency). Failures are counted per rank and reduced.
#ifndef FDS_TEST_CHECK_H_
#define FDS_TEST_CHECK_H_

#include <AMReX_ParallelDescriptor.H>

#include <cstdio>
#include <string>

namespace fdstest {

struct Counter {
    long checks = 0;
    long failures = 0;
    int shown = 0;
};
inline Counter& counter() { static Counter c; return c; }

inline void check(bool ok, const char* expr, const char* file, int line, const std::string& detail = std::string())
{
    Counter& c = counter();
    ++c.checks;
    if (!ok) {
        ++c.failures;
        if (c.shown++ < 20)
            std::fprintf(stderr, "  FAIL [rank %d] %s:%d: %s %s\n", amrex::ParallelDescriptor::MyProc(), file, line, expr, detail.c_str());
    }
}

// Sums over ranks; returns total failures. Prints one summary line on the I/O rank.
inline long report(const char* name)
{
    long v[2] = {counter().checks, counter().failures};
    amrex::ParallelDescriptor::ReduceLongSum(v, 2);
    if (amrex::ParallelDescriptor::IOProcessor())
        std::printf("%s %s: %ld checks, %ld failures\n", v[1] == 0 ? "PASS" : "FAIL", name, v[0], v[1]);
    counter() = Counter();
    return v[1];
}

}  // namespace fdstest

#define CHECK(cond) fdstest::check((cond), #cond, __FILE__, __LINE__)
#define CHECK_MSG(cond, msg) fdstest::check((cond), #cond, __FILE__, __LINE__, (msg))

#endif
