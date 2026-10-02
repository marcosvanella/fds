// AmrInput.cpp: namelist scanner and &AMR / &AMR_REGION parser (R0, IR-003, FR-010). See AmrInput.H.
#include "AmrInput.H"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <sstream>

namespace fdsrt {

bool Report::has_error_containing(const std::string& s) const
{
    for (const auto& e : errors)
        if (e.find(s) != std::string::npos) return true;
    return false;
}

namespace {

std::string upper(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

// Token reader inside one group. Returns false at the closing '/'.
struct Cursor {
    const std::string& s;
    size_t i;
    int line;
    explicit Cursor(const std::string& str, size_t pos, int ln) : s(str), i(pos), line(ln) {}
    void skip_ws()
    {
        while (i < s.size() && (is_space(s[i]) || s[i] == ',')) {
            if (s[i] == '\n') ++line;
            ++i;
        }
    }
};

bool read_token(Cursor& c, std::string& tok, bool& quoted)
{
    c.skip_ws();
    if (c.i >= c.s.size() || c.s[c.i] == '/') return false;
    quoted = false;
    tok.clear();
    char q = c.s[c.i];
    if (q == '\'' || q == '"') {
        quoted = true;
        ++c.i;
        while (c.i < c.s.size() && c.s[c.i] != q) tok += c.s[c.i++];
        if (c.i < c.s.size()) ++c.i;
        return true;
    }
    if (q == '=') {  // stray '='
        tok = "=";
        ++c.i;
        return true;
    }
    while (c.i < c.s.size() && !is_space(c.s[c.i]) && c.s[c.i] != ',' && c.s[c.i] != '=' && c.s[c.i] != '/') tok += c.s[c.i++];
    return true;
}

// True when the next non-blank character is '=' (the token just read was a key).
bool next_is_equals(const Cursor& c)
{
    size_t j = c.i;
    while (j < c.s.size() && is_space(c.s[j])) ++j;
    return j < c.s.size() && c.s[j] == '=';
}

bool to_int(const std::string& t, int& v)
{
    char* e = nullptr;
    errno = 0;
    long x = std::strtol(t.c_str(), &e, 10);
    if (e == t.c_str() || *e != '\0' || errno != 0) return false;
    v = static_cast<int>(x);
    return true;
}
bool to_double(const std::string& t, double& v)
{
    char* e = nullptr;
    errno = 0;
    double x = std::strtod(t.c_str(), &e);
    if (e == t.c_str() || *e != '\0' || errno != 0) return false;
    v = x;
    return true;
}

std::string at_line(const NamelistGroup& g) { return "&" + g.name + " (input line " + std::to_string(g.line) + "): "; }

bool get_int(const NamelistGroup& g, const std::string& key, int& out, Report& rep)
{
    auto it = g.values.find(key);
    if (it == g.values.end()) return false;
    if (it->second.size() != 1 || !to_int(it->second[0], out)) {
        rep.error(at_line(g) + key + " must be one integer");
        return false;
    }
    return true;
}
bool get_double(const NamelistGroup& g, const std::string& key, double& out, Report& rep)
{
    auto it = g.values.find(key);
    if (it == g.values.end()) return false;
    if (it->second.size() != 1 || !to_double(it->second[0], out)) {
        rep.error(at_line(g) + key + " must be one number");
        return false;
    }
    return true;
}
bool get_int_list(const NamelistGroup& g, const std::string& key, std::vector<int>& out, Report& rep)
{
    auto it = g.values.find(key);
    if (it == g.values.end()) return false;
    std::vector<int> v;
    for (const auto& t : it->second) {
        int x;
        if (!to_int(t, x)) {
            rep.error(at_line(g) + key + " must be a list of integers, found '" + t + "'");
            return false;
        }
        v.push_back(x);
    }
    if (v.empty()) {
        rep.error(at_line(g) + key + " has no value");
        return false;
    }
    out = v;
    return true;
}

void check_keys(const NamelistGroup& g, const std::vector<std::string>& allowed, Report& rep)
{
    for (const auto& kv : g.values) {
        bool known = false;
        for (const auto& a : allowed) known = known || (a == kv.first);
        if (!known) rep.error(at_line(g) + "unknown parameter " + kv.first);
    }
}

}  // namespace

std::vector<NamelistGroup> scan_namelists(const std::string& text)
{
    std::vector<NamelistGroup> out;
    size_t i = 0;
    int line = 1;
    while (i < text.size()) {
        if (text[i] == '\n') ++line;
        if (text[i] != '&') {
            ++i;
            continue;
        }
        size_t j = i + 1;
        std::string name;
        while (j < text.size() && (std::isalnum(static_cast<unsigned char>(text[j])) || text[j] == '_')) name += text[j++];
        name = upper(name);
        if (name.empty()) {  // an ampersand in a comment is not a namelist
            i = j;
            continue;
        }
        if (name == "TAIL") break;
        NamelistGroup g;
        g.name = name;
        g.line = line;
        Cursor c(text, j, line);
        std::string key;
        std::string tok;
        bool quoted = false;
        while (read_token(c, tok, quoted)) {
            if (!quoted && next_is_equals(c)) {
                key = upper(tok);
                while (c.i < text.size() && is_space(text[c.i])) ++c.i;
                ++c.i;  // '='
                g.values[key];
            } else if (!key.empty() && tok != "=") {
                g.values[key].push_back(tok);
            }
        }
        out.push_back(g);
        i = c.i < text.size() ? c.i + 1 : c.i;
        line = c.line;
    }
    return out;
}

AmrParams parse_amr_params(const std::string& text, Report& rep)
{
    AmrParams p;
    const std::vector<std::string> keys = {"MAX_LEVEL", "REF_RATIO", "REGRID_INTERVAL", "BLOCKING_FACTOR", "MAX_GRID_SIZE", "N_ERROR_BUF",
                                           "N_PROPER", "GRID_EFF", "OUTPUT_LEVEL_CAP", "VELOCITY_TRANSFER"};
    int n_amr = 0;
    for (const auto& g : scan_namelists(text)) {
        if (g.name == "AMR") {
            if (++n_amr > 1) {
                rep.error(at_line(g) + "more than one &AMR line");
                continue;
            }
            p.present = true;
            check_keys(g, keys, rep);
            get_int(g, "MAX_LEVEL", p.max_level, rep);
            get_int_list(g, "REF_RATIO", p.ref_ratio, rep);
            get_int(g, "REGRID_INTERVAL", p.regrid_interval, rep);
            get_int_list(g, "BLOCKING_FACTOR", p.blocking_factor, rep);
            get_int_list(g, "MAX_GRID_SIZE", p.max_grid_size, rep);
            get_int(g, "N_ERROR_BUF", p.n_error_buf, rep);
            get_int(g, "N_PROPER", p.n_proper, rep);
            get_double(g, "GRID_EFF", p.grid_eff, rep);
            get_int(g, "OUTPUT_LEVEL_CAP", p.output_level_cap, rep);
            auto vt = g.values.find("VELOCITY_TRANSFER");
            if (vt != g.values.end()) {
                std::string v = vt->second.size() == 1 ? upper(vt->second[0]) : std::string();
                if (v == "FACE_DIV_FREE") p.velocity_transfer = VelocityTransfer::FaceDivFree;
                else if (v == "FACE_LINEAR") p.velocity_transfer = VelocityTransfer::FaceLinear;
                else if (v == "FACE_CONSERVATIVE") p.velocity_transfer = VelocityTransfer::FaceConservative;
                else rep.error(at_line(g) + "VELOCITY_TRANSFER must be 'FACE_DIV_FREE', 'FACE_LINEAR' or 'FACE_CONSERVATIVE'");
            }
        }
    }
    // Regions are read after &AMR so that LEVEL can default to MAX_LEVEL.
    for (const auto& g : scan_namelists(text)) {
        if (g.name != "AMR_REGION") continue;
        check_keys(g, {"XB", "LEVEL"}, rep);
        RegionSpec r;
        auto xb = g.values.find("XB");
        bool ok = xb != g.values.end() && xb->second.size() == 6;
        for (size_t k = 0; ok && k < 6; ++k) ok = to_double(xb->second[k], r.xb[k]);
        if (!ok) {
            rep.error(at_line(g) + "XB must be six numbers");
            continue;
        }
        for (int d = 0; d < 3; ++d)
            if (r.xb[2 * d] >= r.xb[2 * d + 1]) rep.error(at_line(g) + "XB has zero or negative extent in direction " + std::to_string(d + 1));
        int lv = -1;
        get_int(g, "LEVEL", lv, rep);
        r.level = lv;
        p.regions.push_back(r);
    }
    if (!p.regions.empty() && !p.present) rep.error("&AMR_REGION found without an &AMR line");
    validate_amr_params(p, rep);
    return p;
}

void validate_amr_params(const AmrParams& p, Report& rep)
{
    if (p.max_level < 0) rep.error("&AMR: MAX_LEVEL must be >= 0");
    if (p.max_level > 8) rep.error("&AMR: MAX_LEVEL = " + std::to_string(p.max_level) + " is above the supported limit of 8");
    const int ml = std::max(p.max_level, 0);
    auto len_ok = [&](const std::vector<int>& v, size_t n) { return v.size() == 1 || v.size() == n; };
    if (!len_ok(p.ref_ratio, static_cast<size_t>(ml)) && !(ml == 0))
        rep.error("&AMR: REF_RATIO needs one value or MAX_LEVEL = " + std::to_string(ml) + " values, found " + std::to_string(p.ref_ratio.size()));
    if (!len_ok(p.blocking_factor, static_cast<size_t>(ml + 1)))
        rep.error("&AMR: BLOCKING_FACTOR needs one value or MAX_LEVEL+1 = " + std::to_string(ml + 1) + " values, found " + std::to_string(p.blocking_factor.size()));
    if (!len_ok(p.max_grid_size, static_cast<size_t>(ml + 1)))
        rep.error("&AMR: MAX_GRID_SIZE needs one value or MAX_LEVEL+1 = " + std::to_string(ml + 1) + " values, found " + std::to_string(p.max_grid_size.size()));
    // Ratios: FR-010 supports 2 and 4 only (MLMG coarse-fine interpolation asserts ratio <= 4; ratio 3 is untested).
    for (size_t k = 0; k < p.ref_ratio.size(); ++k)
        if (p.ref_ratio[k] != 2 && p.ref_ratio[k] != 4)
            rep.error("&AMR: REF_RATIO entry " + std::to_string(k + 1) + " is " + std::to_string(p.ref_ratio[k]) +
                      "; AMR mode supports refinement ratios 2 and 4 only (FR-010)");
    for (size_t k = 0; k < p.blocking_factor.size(); ++k) {
        int b = p.blocking_factor[k];
        if (b < 1 || (b & (b - 1)) != 0)
            rep.error("&AMR: BLOCKING_FACTOR entry " + std::to_string(k + 1) + " is " + std::to_string(b) + "; it must be a power of 2");
    }
    for (size_t k = 0; k < p.max_grid_size.size(); ++k)
        if (p.max_grid_size[k] < 1) rep.error("&AMR: MAX_GRID_SIZE entry " + std::to_string(k + 1) + " must be >= 1");
    const bool lens_ok0 = (p.blocking_factor.size() == 1 || p.blocking_factor.size() == static_cast<size_t>(ml + 1)) &&
                          (p.max_grid_size.size() == 1 || p.max_grid_size.size() == static_cast<size_t>(ml + 1));
    for (int l = 0; lens_ok0 && l <= ml; ++l) {
        int b = AmrParams::at(p.blocking_factor, l), m = AmrParams::at(p.max_grid_size, l);
        if (b >= 1 && m >= 1 && m % b != 0)
            rep.error("&AMR: MAX_GRID_SIZE " + std::to_string(m) + " on level " + std::to_string(l) + " is not a multiple of the blocking factor " +
                      std::to_string(b));
    }
    // AMReX: the blocking factor may not grow faster than blocking_factor * ref_ratio between levels (AMReX_AmrMesh.cpp checkInput).
    const bool rr_ok = p.ref_ratio.size() == 1 || p.ref_ratio.size() == static_cast<size_t>(ml);
    for (int l = 0; lens_ok0 && rr_ok && l < ml; ++l) {
        int rr = AmrParams::at(p.ref_ratio, l), b0 = AmrParams::at(p.blocking_factor, l), b1 = AmrParams::at(p.blocking_factor, l + 1);
        if (b0 * rr < b1)
            rep.error("&AMR: BLOCKING_FACTOR " + std::to_string(b0) + " on level " + std::to_string(l) + " with ratio " + std::to_string(rr) +
                      " does not allow blocking factor " + std::to_string(b1) + " on level " + std::to_string(l + 1) +
                      " (blocking factors may grow at most by the refinement ratio; for example '2,8' needs ratio 4, ratio 2 needs '2,4')");
    }
    if (p.regrid_interval < 0) rep.error("&AMR: REGRID_INTERVAL must be >= 0 (0 = no regrid)");
    if (p.n_error_buf < 0) rep.error("&AMR: N_ERROR_BUF must be >= 0");
    if (p.n_proper < 1) rep.error("&AMR: N_PROPER must be >= 1");
    if (!(p.grid_eff > 0.0 && p.grid_eff <= 1.0)) rep.error("&AMR: GRID_EFF must be in (0, 1]");
    if (p.output_level_cap < -1 || p.output_level_cap > p.max_level)
        rep.error("&AMR: OUTPUT_LEVEL_CAP must be between 0 and MAX_LEVEL");
    for (size_t k = 0; k < p.regions.size(); ++k) {
        int lv = p.regions[k].level;
        if (lv != -1 && (lv < 1 || lv > p.max_level))
            rep.error("&AMR_REGION " + std::to_string(k + 1) + ": LEVEL must be between 1 and MAX_LEVEL");
    }
    if (p.max_level == 0 && p.present && !p.regions.empty())
        rep.warn("&AMR_REGION lines are ignored with MAX_LEVEL = 0");
}

}  // namespace fdsrt
