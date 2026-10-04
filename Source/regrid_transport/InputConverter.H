// InputConverter.H: the input converter of D-076 as pure functions (no AMReX, no FDS). The driver pre-pass reads a multi-mesh FDS input whose &MESH lines have
// different resolutions, classifies the meshes into levels (IR-002, Hierarchy.H), checks nesting and alignment with the existing rules, and writes a level-0-only
// input text for the unchanged FDS set-up: the finer &MESH lines (and the &AMR / &AMR_REGION groups, which FDS does not know) are removed and everything else is
// copied byte for byte. The level-0 cells that lay under a removed finer mesh get one added &MESH line each (the cover), so that the level 0 FDS sets up tiles the domain. The hierarchy built from the removed finer meshes is returned with the parsed &AMR parameters, so that the driver can hand the level-0
// meshes to assemble_level0 and the hierarchy to RegridAmrCore. Nothing here writes a file or calls FDS.
//
// Reading rules (the ones FDS uses): a namelist group starts at an '&' that is the first non-blank character of a line and ends at the first '/' outside quotes
// (so a commented-out "! &MESH ..." line is not a group); scanning stops at &TAIL. &MESH lines are expanded with their &MULT exactly as READ_MESH does (copies with
// I fastest, then J, then K; the skip box of the I/J/K_LOWER_SKIP and _UPPER_SKIP; N_LOWER/N_UPPER sequential form), so mesh i of the result is FDS mesh i+1 of the
// ORIGINAL input. In the emitted input the level-0 meshes keep their order, so level-0 mesh i is FDS mesh i+1 of the converted input (Grouping::level0_meshes()).
#ifndef FDSRT_INPUTCONVERTER_H_
#define FDSRT_INPUTCONVERTER_H_

#include <array>
#include <string>
#include <vector>

#include "AmrInput.H"
#include "Hierarchy.H"

namespace fdsrt {

// A namelist group found in the text: name (upper case), the text range [begin, end) from the '&' through the closing '/', and its 1-based first line.
struct GroupSpan {
    std::string name;
    size_t begin = 0;
    size_t end = 0;
    int line = 0;
};

// Finds the groups of `text` with the FDS line-start rule. An unterminated group is an error.
std::vector<GroupSpan> find_group_spans(const std::string& text, Report& rep);

// One &MESH line of the input and the meshes it expands to (`first` is the 0-based index of its first mesh in the expanded list, `count` the number of copies).
struct MeshLine {
    GroupSpan span;
    int first = 0;
    int count = 0;
    std::string id;           // ID= of the line ("" if none); with MULT_ID the copies share it
    bool explicit_rank = false;
};

// Expands the &MESH lines (with &MULT) into meshes in FDS order. `lines` receives one entry per &MESH group. Errors: unknown MULT_ID, missing &MESH.
bool parse_meshes(const std::string& text, const std::vector<GroupSpan>& spans, std::vector<MeshLine>& lines, std::vector<MeshInput>& meshes, Report& rep);

// Directions (0 x, 1 y, 2 z) in which the input has a &VENT with SURF_ID='PERIODIC' (DB='XMIN' ..., PBX/PBY/PBZ, or a planar XB).
std::array<bool, 3> detect_periodic(const std::string& text, const std::vector<GroupSpan>& spans);

struct ConvertResult {
    bool ok = false;
    AmrParams params;                      // parsed &AMR / &AMR_REGION (present == false without an &AMR line)
    std::vector<MeshLine> mesh_lines;      // all &MESH lines of the original input
    std::vector<MeshInput> meshes;         // all meshes of the original input, FDS order
    std::vector<int> level_of_mesh;        // level of each original mesh (0-based mesh index)
    std::vector<int> level0_meshes;        // original 0-based mesh indices of the level-0 meshes, in order = the first meshes of the converted input
    std::vector<IBox> cover_boxes;         // level-0 cells under the finer meshes: the converted input holds one added &MESH line per box, after the level-0 meshes (converted
                                           // mesh number n_level0_meshes + i + 1 for cover box i), so the converted level 0 tiles the domain and equals hierarchy.levels[0].grids
    std::vector<int> removed_meshes;       // original 0-based mesh indices of the finer meshes
    std::array<bool, 3> periodic{false, false, false};
    Grouping grouping;
    Hierarchy hierarchy;                   // built from the original meshes; level-0 boxes are the level-0 meshes in the converted order, then the cover under finer meshes
    std::string level0_text;               // the converted input, for fds_setup
    int n_removed_lines = 0;
};

// The whole pre-pass on the text of an FDS input. Returns false with errors in `rep` when the input must be rejected:
//  - finer meshes but no &AMR line: the error names the &AMR line and a mesh pair (AMR mode is never inferred),
//  - a finer mesh not aligned with its parent level's cell faces, ratios other than 2 and 4, overlapping meshes, nesting and blocking-factor violations (Hierarchy.H),
//  - a &MESH line whose copies fall on different levels, a MESH_ID= of another group that names a removed mesh, MPI_PROCESS numbers that are no longer continuous
//    after the removal (all level-0 lines must then give MPI_PROCESS).
// An input without finer meshes is converted too (nothing removed but the &AMR groups) so that the driver has one path.
bool convert_input(const std::string& text, ConvertResult& out, Report& rep);

// Text of the level-0-only input: the groups of `spans` copied unchanged, except the &AMR and &AMR_REGION groups and the &MESH lines of level > 0, which become
// one comment line each (a '!' line; the rest of the removed line is dropped). `level_of_mesh` has one entry per expanded mesh. `extra_mesh_lines` (the cover
// meshes, one line each, no final newline) is inserted after the last &MESH line. Pure; used by convert_input.
std::string emit_level0_text(const std::string& text, const std::vector<GroupSpan>& spans, const std::vector<MeshLine>& lines, const std::vector<int>& level_of_mesh,
                             const std::string& extra_mesh_lines, int& n_removed_lines);

}  // namespace fdsrt

#endif
