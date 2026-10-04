// ConvertInputTool.cpp: command line front end of the D-076 input converter (InputConverter.H), for checking an input by hand and for the driver pre-pass design.
//   fds_amr_convert_input <input.fds> [-o <level0.fds>] [--dump]
// Prints errors and warnings (exit 1 on error), the hierarchy dump with --dump, and writes the level-0-only input with -o. No FDS, no AMReX.
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "InputConverter.H"

int main(int argc, char** argv)
{
    const char* in = nullptr;
    const char* outpath = nullptr;
    bool dump = false;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "-o") && i + 1 < argc) outpath = argv[++i];
        else if (!std::strcmp(argv[i], "--dump")) dump = true;
        else if (!in) in = argv[i];
        else { std::fprintf(stderr, "usage: %s <input.fds> [-o <level0.fds>] [--dump]\n", argv[0]); return 2; }
    }
    if (!in) { std::fprintf(stderr, "usage: %s <input.fds> [-o <level0.fds>] [--dump]\n", argv[0]); return 2; }
    std::ifstream f(in);
    if (!f) { std::fprintf(stderr, "cannot read %s\n", in); return 2; }
    std::stringstream ss;
    ss << f.rdbuf();
    fdsrt::Report rep;
    fdsrt::ConvertResult r;
    const bool ok = fdsrt::convert_input(ss.str(), r, rep);
    for (const auto& w : rep.warnings) std::fprintf(stderr, "warning: %s\n", w.c_str());
    for (const auto& e : rep.errors) std::fprintf(stderr, "error: %s\n", e.c_str());
    if (!ok || !rep.ok()) return 1;
    std::printf("meshes %zu, level 0 meshes %zu, finer meshes removed %zu (%d lines), cover meshes added %zu, levels %d\n", r.meshes.size(), r.level0_meshes.size(),
                r.removed_meshes.size(), r.n_removed_lines, r.cover_boxes.size(), r.hierarchy.top + 1);
    if (dump) fdsrt::dump_hierarchy(std::cout, r.hierarchy);
    if (outpath) {
        std::ofstream o(outpath);
        o << r.level0_text;
        if (!o) { std::fprintf(stderr, "cannot write %s\n", outpath); return 2; }
    }
    return 0;
}
