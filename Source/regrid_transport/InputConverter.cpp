// InputConverter.cpp: the D-076 input converter (pure functions); see InputConverter.H.
#include "InputConverter.H"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <set>
#include <sstream>

namespace fdsrt {

namespace {

std::string up(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// The single group parsed from the text of a span (the span starts at its '&').
NamelistGroup parse_span(const std::string& text, const GroupSpan& s)
{
    auto v = scan_namelists(text.substr(s.begin, s.end - s.begin));
    NamelistGroup g;
    if (!v.empty()) g = v.front();
    g.name = s.name;
    g.line = s.line;
    return g;
}

std::string where(const GroupSpan& s) { return "&" + s.name + " (input line " + std::to_string(s.line) + ")"; }

bool num(const NamelistGroup& g, const char* key, size_t i, double& v)
{
    auto it = g.values.find(key);
    if (it == g.values.end() || it->second.size() <= i) return false;
    std::string t = it->second[i];
    for (char& c : t)
        if (c == 'd' || c == 'D') c = 'e';
    char* e = nullptr;
    v = std::strtod(t.c_str(), &e);
    return e != t.c_str() && *e == '\0';
}
double numd(const NamelistGroup& g, const char* key, double dflt)
{
    double v = dflt;
    return num(g, key, 0, v) ? v : dflt;
}
int numi(const NamelistGroup& g, const char* key, int dflt) { return static_cast<int>(std::lround(numd(g, key, dflt))); }
std::string str1(const NamelistGroup& g, const char* key)
{
    auto it = g.values.find(key);
    return it == g.values.end() || it->second.empty() ? std::string() : it->second[0];
}

// MULT of READ_MULT: offsets, index ranges and the skip box, as FDS builds them.
struct Mult {
    double dx0[3] = {0, 0, 0};
    double dxb[6] = {0, 0, 0, 0, 0, 0};
    int lo[3] = {0, 0, 0}, up[3] = {0, 0, 0};
    bool sequential = false;
    int skip_lo[3] = {-999, -999, -999}, skip_hi[3] = {999, 999, 999};
    bool has_skip = false;
    bool skipped(int i, int j, int k) const
    {
        if (!has_skip) return false;
        const int ijk[3] = {i, j, k};
        for (int d = 0; d < 3; ++d) {
            const int a = std::max(lo[d], skip_lo[d]), b = std::min(up[d], skip_hi[d]);
            if (ijk[d] < a || ijk[d] > b) return false;
        }
        return true;
    }
};

Mult make_mult(const NamelistGroup& g)
{
    Mult m;
    m.dx0[0] = numd(g, "DX0", 0); m.dx0[1] = numd(g, "DY0", 0); m.dx0[2] = numd(g, "DZ0", 0);
    for (int q = 0; q < 6; ++q) { double v = 0; if (num(g, "DXB", q, v)) m.dxb[q] = v; }
    const double d[3] = {numd(g, "DX", 0), numd(g, "DY", 0), numd(g, "DZ", 0)};
    for (int a = 0; a < 3; ++a)
        if (std::fabs(d[a]) > 1e-18) m.dxb[2 * a] = m.dxb[2 * a + 1] = d[a];
    m.lo[0] = numi(g, "I_LOWER", 0); m.up[0] = numi(g, "I_UPPER", 0);
    m.lo[1] = numi(g, "J_LOWER", 0); m.up[1] = numi(g, "J_UPPER", 0);
    m.lo[2] = numi(g, "K_LOWER", 0); m.up[2] = numi(g, "K_UPPER", 0);
    m.skip_lo[0] = numi(g, "I_LOWER_SKIP", -999); m.skip_hi[0] = numi(g, "I_UPPER_SKIP", 999);
    m.skip_lo[1] = numi(g, "J_LOWER_SKIP", -999); m.skip_hi[1] = numi(g, "J_UPPER_SKIP", 999);
    m.skip_lo[2] = numi(g, "K_LOWER_SKIP", -999); m.skip_hi[2] = numi(g, "K_UPPER_SKIP", 999);
    if (g.values.count("N_LOWER") || g.values.count("N_UPPER")) {
        m.sequential = true;
        m.lo[0] = numi(g, "N_LOWER", 0); m.up[0] = numi(g, "N_UPPER", 0);
        m.lo[1] = m.up[1] = m.lo[2] = m.up[2] = 0;
        m.skip_lo[0] = numi(g, "N_LOWER_SKIP", -999); m.skip_hi[0] = numi(g, "N_UPPER_SKIP", 999);
        m.skip_lo[1] = m.skip_lo[2] = -999; m.skip_hi[1] = m.skip_hi[2] = 999;
    }
    for (int d = 0; d < 3; ++d) m.has_skip = m.has_skip || m.skip_lo[d] >= m.lo[d] || m.skip_hi[d] <= m.up[d];
    return m;
}

}  // namespace

std::vector<GroupSpan> find_group_spans(const std::string& text, Report& rep)
{
    std::vector<GroupSpan> out;
    size_t i = 0;
    int line = 1;
    while (i < text.size()) {
        size_t j = i;
        while (j < text.size() && (text[j] == ' ' || text[j] == '\t')) ++j;
        if (j + 1 < text.size() && text[j] == '&' && (std::isalpha(static_cast<unsigned char>(text[j + 1])))) {
            size_t k = j + 1;
            std::string name;
            while (k < text.size() && (std::isalnum(static_cast<unsigned char>(text[k])) || text[k] == '_')) name += text[k++];
            name = up(name);
            if (name == "TAIL") break;
            // closing '/' outside quotes
            char q = 0;
            size_t e = k;
            int nl = 0;
            bool closed = false;
            for (; e < text.size(); ++e) {
                const char c = text[e];
                if (c == '\n') ++nl;
                if (q) { if (c == q) q = 0; continue; }
                if (c == '\'' || c == '"') { q = c; continue; }
                if (c == '/') { closed = true; break; }
            }
            if (!closed) {
                rep.error("&" + name + " (input line " + std::to_string(line) + ") has no closing '/'");
                return out;
            }
            GroupSpan s;
            s.name = name;
            s.begin = j;
            s.end = e + 1;
            s.line = line;
            out.push_back(s);
            line += nl;
            i = e + 1;  // the rest of this line (a comment) is skipped by the line loop below
        }
        // advance to the start of the next line
        while (i < text.size() && text[i] != '\n') ++i;
        if (i < text.size()) { ++i; ++line; }
    }
    return out;
}

bool parse_meshes(const std::string& text, const std::vector<GroupSpan>& spans, std::vector<MeshLine>& lines, std::vector<MeshInput>& meshes, Report& rep)
{
    const size_t nerr0 = rep.errors.size();
    lines.clear();
    meshes.clear();
    std::vector<std::pair<std::string, Mult>> mults;
    for (const GroupSpan& s : spans)
        if (s.name == "MULT") mults.emplace_back(str1(parse_span(text, s), "ID"), make_mult(parse_span(text, s)));
    for (const GroupSpan& s : spans) {
        if (s.name != "MESH") continue;
        const NamelistGroup g = parse_span(text, s);
        MeshInput m;
        m.ijk[0] = m.ijk[1] = m.ijk[2] = 10;
        m.xb[0] = m.xb[2] = m.xb[4] = 0.0;
        m.xb[1] = m.xb[3] = m.xb[5] = 1.0;
        for (int d = 0; d < 3; ++d) { double v = 0; if (num(g, "IJK", d, v)) m.ijk[d] = static_cast<int>(std::lround(v)); }
        for (int q = 0; q < 6; ++q) { double v = 0; if (num(g, "XB", q, v)) m.xb[q] = v; }
        auto ri = g.values.find("MPI_PROCESS");
        MeshLine ml;
        ml.span = s;
        ml.id = str1(g, "ID");
        ml.explicit_rank = ri != g.values.end() && !ri->second.empty();
        m.rank = ml.explicit_rank ? numi(g, "MPI_PROCESS", 0) : 0;
        ml.first = static_cast<int>(meshes.size());
        const std::string mid = str1(g, "MULT_ID");
        if (mid.empty() || up(mid) == "NULL") {
            meshes.push_back(m);
        } else {
            const Mult* mr = nullptr;
            for (const auto& pr : mults)
                if (pr.first == mid) mr = &pr.second;
            if (!mr) {
                rep.error(where(s) + ": MULT_ID '" + mid + "' not found");
                return false;
            }
            for (int k = mr->lo[2]; k <= mr->up[2]; ++k)
                for (int j = mr->lo[1]; j <= mr->up[1]; ++j)
                    for (int i = mr->lo[0]; i <= mr->up[0]; ++i) {
                        if (mr->skipped(i, j, k)) continue;
                        MeshInput c = m;
                        const int ii[3] = {i, mr->sequential ? i : j, mr->sequential ? i : k};
                        for (int d = 0; d < 3; ++d) {
                            c.xb[2 * d] += mr->dx0[d] + ii[d] * mr->dxb[2 * d];
                            c.xb[2 * d + 1] += mr->dx0[d] + ii[d] * mr->dxb[2 * d + 1];
                        }
                        meshes.push_back(c);
                    }
        }
        ml.count = static_cast<int>(meshes.size()) - ml.first;
        lines.push_back(ml);
    }
    if (lines.empty()) rep.error("no &MESH line in the input");
    return rep.errors.size() == nerr0;
}

std::array<bool, 3> detect_periodic(const std::string& text, const std::vector<GroupSpan>& spans)
{
    std::array<bool, 3> per{false, false, false};
    for (const GroupSpan& s : spans) {
        if (s.name != "VENT") continue;
        const NamelistGroup g = parse_span(text, s);
        if (up(str1(g, "SURF_ID")) != "PERIODIC") continue;
        const std::string db = up(str1(g, "DB"));
        if (!db.empty()) {
            if (db[0] == 'X') per[0] = true;
            else if (db[0] == 'Y') per[1] = true;
            else if (db[0] == 'Z') per[2] = true;
        } else if (g.values.count("PBX")) per[0] = true;
        else if (g.values.count("PBY")) per[1] = true;
        else if (g.values.count("PBZ")) per[2] = true;
        else {
            for (int d = 0; d < 3; ++d) {
                double a = 0, b = 0;
                if (num(g, "XB", 2 * d, a) && num(g, "XB", 2 * d + 1, b) && a == b) per[d] = true;
            }
        }
    }
    return per;
}

std::string emit_level0_text(const std::string& text, const std::vector<GroupSpan>& spans, const std::vector<MeshLine>& lines, const std::vector<int>& level_of_mesh,
                             const std::string& extra_mesh_lines, int& n_removed_lines)
{
    n_removed_lines = 0;
    std::string out;
    size_t pos = 0;
    size_t mesh_line = 0;
    for (const GroupSpan& s : spans) {
        std::string note;
        bool last_mesh = false;
        if (s.name == "AMR" || s.name == "AMR_REGION") {
            note = "! the " + s.name + " group is read by the AMR driver pre-pass and removed from this input";
        } else if (s.name == "MESH") {
            const MeshLine& ml = lines[mesh_line];
            last_mesh = ++mesh_line == lines.size();
            const int lev = level_of_mesh[ml.first];
            if (lev > 0) {
                note = "! converter: finer MESH line removed (level " + std::to_string(lev) + ", input mesh " + std::to_string(ml.first + 1) +
                       (ml.count > 1 ? " to " + std::to_string(ml.first + ml.count) : std::string()) + ")";
                ++n_removed_lines;
            }
        }
        const bool insert = last_mesh && !extra_mesh_lines.empty();
        if (note.empty() && !insert) continue;
        out.append(text, pos, s.begin - pos);
        size_t e = s.end;
        while (e < text.size() && text[e] != '\n') ++e;     // end of the line of the group
        if (!note.empty()) {
            out += note;
            pos = (text.find('&', s.end) >= e) ? e : s.end;  // drop the rest of the line (a trailing comment) unless it holds another group
        } else {
            out.append(text, s.begin, e - s.begin);
            pos = e;
        }
        if (insert) out += "\n" + extra_mesh_lines;
    }
    out.append(text, pos, std::string::npos);
    return out;
}

bool convert_input(const std::string& text, ConvertResult& out, Report& rep)
{
    out = ConvertResult();
    const size_t nerr0 = rep.errors.size();
    const std::vector<GroupSpan> spans = find_group_spans(text, rep);
    if (rep.errors.size() > nerr0) return false;

    // &AMR, &AMR_REGION, &MISC (EXACT_SUMS) with their original line numbers: blank lines stand in for the other groups
    {
        std::string sub;
        int cur = 1;
        for (const GroupSpan& s : spans) {
            if (s.name != "AMR" && s.name != "AMR_REGION" && s.name != "MISC") continue;
            if (s.line > cur) sub.append(static_cast<size_t>(s.line - cur), '\n');
            cur = std::max(cur, s.line);
            sub.append(text, s.begin, s.end - s.begin);
            cur += static_cast<int>(std::count(text.begin() + s.begin, text.begin() + s.end, '\n'));
            sub += ' ';
        }
        out.params = parse_amr_params(sub, rep);
        if (rep.errors.size() > nerr0) return false;
    }

    if (!parse_meshes(text, spans, out.mesh_lines, out.meshes, rep)) return false;
    out.periodic = detect_periodic(text, spans);

    out.grouping = group_meshes(out.meshes, out.params, out.periodic, rep);
    if (!out.grouping.ok || rep.errors.size() > nerr0) return false;
    {
        IVec coarsest;
        mlmg_coarsening(out.grouping, coarsest, rep);
    }
    if (!build_static_hierarchy(out.grouping, out.params, out.hierarchy, rep) || rep.errors.size() > nerr0) return false;

    out.level_of_mesh = out.grouping.level_of_mesh;
    const int nm = static_cast<int>(out.meshes.size());
    for (int m = 0; m < nm; ++m) (out.level_of_mesh[m] == 0 ? out.level0_meshes : out.removed_meshes).push_back(m);

    // a &MESH line (MULT copies) must not straddle levels
    for (const MeshLine& ml : out.mesh_lines)
        for (int c = 1; c < ml.count; ++c)
            if (out.level_of_mesh[ml.first + c] != out.level_of_mesh[ml.first]) {
                rep.error(where(ml.span) + ": its MULT copies fall on different levels (" + std::to_string(out.level_of_mesh[ml.first]) + " and " +
                          std::to_string(out.level_of_mesh[ml.first + c]) + "); a line must hold one resolution");
                return false;
            }

    // other groups must not name a removed mesh by MESH_ID
    {
        std::set<std::string> gone;
        for (const MeshLine& ml : out.mesh_lines)
            if (!ml.id.empty() && ml.count > 0 && out.level_of_mesh[ml.first] > 0) gone.insert(ml.id);
        if (!gone.empty())
            for (const GroupSpan& s : spans) {
                if (s.name == "MESH") continue;
                if (text.find("MESH_ID", s.begin) >= s.end && text.find("mesh_id", s.begin) >= s.end) continue;
                const NamelistGroup g = parse_span(text, s);
                const std::string v = str1(g, "MESH_ID");
                if (gone.count(v)) {
                    rep.error(where(s) + ": MESH_ID='" + v + "' names a finer mesh that the converter removes from the level-0 input (it exists only on level " +
                              "1 or above); attach this line to a level-0 mesh or to the AMR level");
                    return false;
                }
            }
    }

    // MPI_PROCESS of the remaining meshes must stay continuous (FDS ERROR(117)/(118))
    if (!out.removed_meshes.empty()) {
        bool all_explicit = true;
        for (const MeshLine& ml : out.mesh_lines)
            if (out.level_of_mesh[ml.first] == 0 && !ml.explicit_rank) all_explicit = false;
        if (all_explicit) {
            int prev = -1;
            for (int m : out.level0_meshes) {
                const int r = out.meshes[m].rank;
                if ((prev < 0 && r != 0) || (prev >= 0 && (r < prev || r - prev > 1))) {
                    rep.error("MPI_PROCESS is not continuous after the finer meshes are removed (level-0 mesh " + std::to_string(m + 1) + " of the input has MPI_PROCESS " +
                              std::to_string(r) + " after " + std::to_string(prev) + "); renumber the ranks of the level-0 meshes");
                    return false;
                }
                prev = r;
            }
        }
    }

    // cover meshes: the added level-0 grids of the hierarchy (level-0 cells under finer meshes), in grid order
    std::ostringstream cover;
    cover.precision(17);
    int last_rank = 0;
    bool all_explicit = !out.mesh_lines.empty();
    for (const MeshLine& ml : out.mesh_lines)
        if (out.level_of_mesh[ml.first] == 0) { all_explicit = all_explicit && ml.explicit_rank; last_rank = out.meshes[ml.first].rank; }
    for (const GridBox& gb : out.hierarchy.levels[0].grids) {
        if (gb.mesh >= 0) continue;
        out.cover_boxes.push_back(gb.box);
        const IBox& b = gb.box;
        cover << (out.cover_boxes.size() > 1 ? "\n" : "") << "&MESH IJK=" << b.hi[0] - b.lo[0] + 1 << "," << b.hi[1] - b.lo[1] + 1 << "," << b.hi[2] - b.lo[2] + 1 << ", XB=";
        for (int d = 0; d < 3; ++d)
            cover << out.grouping.dom_lo[d] + b.lo[d] * out.grouping.dx0[d] << "," << out.grouping.dom_lo[d] + (b.hi[d] + 1) * out.grouping.dx0[d] << (d < 2 ? "," : "");
        if (all_explicit) cover << ", MPI_PROCESS=" << last_rank;
        cover << " / ! converter: level-0 cover under finer meshes";
    }
    out.level0_text = emit_level0_text(text, spans, out.mesh_lines, out.level_of_mesh, cover.str(), out.n_removed_lines);

    // the converted text, read back as FDS reads it, must be an equal-level-0 input (no mesh pair with a ratio other than 1)
    {
        Report r2;
        std::vector<GroupSpan> spans2 = find_group_spans(out.level0_text, r2);
        std::vector<MeshLine> lines2;
        std::vector<MeshInput> meshes2;
        if (!r2.ok() || !parse_meshes(out.level0_text, spans2, lines2, meshes2, r2)) {
            rep.error("converter self-check: the converted input cannot be read back: " + (r2.errors.empty() ? std::string("?") : r2.errors[0]));
            return false;
        }
        if (meshes2.size() != out.level0_meshes.size() + out.cover_boxes.size()) {
            rep.error("converter self-check: the converted input holds " + std::to_string(meshes2.size()) + " meshes, expected " +
                      std::to_string(out.level0_meshes.size() + out.cover_boxes.size()));
            return false;
        }
        if (!verify_equal_level0(meshes2, out.grouping, rep)) return false;
    }
    out.ok = true;
    return true;
}

bool verify_equal_level0(const std::vector<MeshInput>& meshes, const Grouping& g, Report& rep)
{
    const size_t nerr0 = rep.errors.size();
    std::vector<IBox> boxes(meshes.size());
    for (size_t m = 0; m < meshes.size(); ++m) {
        const MeshInput& mi = meshes[m];
        for (int d = 0; d < 3; ++d) {
            boxes[m].lo[d] = 0;
            boxes[m].hi[d] = 0;
            if (g.hidden[d]) continue;
            const double lo = std::min(mi.xb[2 * d], mi.xb[2 * d + 1]), hi = std::max(mi.xb[2 * d], mi.xb[2 * d + 1]);
            const double dx = (hi - lo) / mi.ijk[d];
            if (std::fabs(dx - g.dx0[d]) > 1e-9 * g.dx0[d]) {
                rep.error("mesh " + std::to_string(m + 1) + " has cell size " + std::to_string(dx) + " in direction " + "xyz"[d] + ", not the level-0 size " +
                          std::to_string(g.dx0[d]) + ": the input is not an equal-level-0 input (an interface with NIC > 1)");
                continue;
            }
            const double q = (lo - g.dom_lo[d]) / g.dx0[d];
            const long k = std::lround(q);
            if (std::fabs(q - k) > 1e-6) {
                rep.error("mesh " + std::to_string(m + 1) + " starts off the level-0 lattice in direction " + "xyz"[d] + " (" + std::to_string(q) + " cells from the domain corner)");
                continue;
            }
            boxes[m].lo[d] = static_cast<int>(k);
            boxes[m].hi[d] = static_cast<int>(k) + mi.ijk[d] - 1;
        }
    }
    if (rep.errors.size() > nerr0) return false;
    long long cells = 0;
    for (size_t a = 0; a < boxes.size(); ++a) {
        long long v = 1;
        for (int d = 0; d < 3; ++d) v *= boxes[a].hi[d] - boxes[a].lo[d] + 1;
        cells += v;
        for (size_t b = a + 1; b < boxes.size(); ++b) {
            bool overlap = true;
            for (int d = 0; d < 3; ++d) overlap = overlap && boxes[a].lo[d] <= boxes[b].hi[d] && boxes[b].lo[d] <= boxes[a].hi[d];
            if (overlap) rep.error("meshes " + std::to_string(a + 1) + " and " + std::to_string(b + 1) + " overlap on the level-0 lattice");
        }
    }
    long long total = 1;
    for (int d = 0; d < 3; ++d) total *= g.hidden[d] ? 1 : g.n0[d];
    if (rep.errors.size() == nerr0 && cells != total)
        rep.error("the meshes hold " + std::to_string(cells) + " level-0 cells, the level-0 domain " + std::to_string(total) + ": they do not tile the domain");
    return rep.errors.size() == nerr0;
}

}  // namespace fdsrt
