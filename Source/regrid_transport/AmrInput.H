// AmrInput.H: the &AMR input of the AMR code (IR-003) as a stand-alone parser. No AMReX dependency and no edit of read.f90 needed:
// the parser takes the text of the FDS input file (or any string holding &AMR / &AMR_REGION lines) and returns a struct plus a list
// of errors and warnings. Wiring into the FDS set-up (so that the same lines are accepted by the FDS namelist reader) is the patch
// list in notes/readf90-cmake-patch-list.md.
//
// Parameter names are working names (IR-003: final once they appear in the User Guide draft).
#ifndef FDSRT_AMRINPUT_H_
#define FDSRT_AMRINPUT_H_

#include <cmath>
#include <limits>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace fdsrt {

// Collects input-time messages. Errors stop the run during set-up; warnings are printed once.
struct Report {
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    void error(const std::string& s) { errors.push_back(s); }
    void warn(const std::string& s) { warnings.push_back(s); }
    bool ok() const { return errors.empty(); }
    bool has_error_containing(const std::string& s) const;
};

// A refinable-region box declared by the user (IR-008 (a)): physical corners and the finest level that may exist inside it.
struct RegionSpec {
    double xb[6] = {0, 0, 0, 0, 0, 0};
    int level = -1;  // finest level allowed in this box; -1 in the input means MAX_LEVEL
};


// POST_REGRID_PROJECTION (ruling D-063): Auto = run the composite projection after a regrid that creates new fine faces next to retained fine faces
// (some level >= 1 keeps part of its old grids and also gets new cells); On = after every regrid that changed the grids; Off = never (the Phase 3 behaviour before D-063).
enum class PostRegridProjection { Auto, On, Off };

// &AMR. Per-level vectors hold the values as entered (one value for all levels, or one per level; see expand()).
struct AmrParams {
    bool present = false;        // an &AMR line was found
    int max_level = 0;           // MAX_LEVEL: finest level allowed (0 = uniform)
    std::vector<int> ref_ratio{2};          // REF_RATIO: ratio from level l to l+1, {2,4} (one value or max_level values)
    int regrid_interval = 4;     // REGRID_INTERVAL: steps between regrids; 0 = hierarchy fixed after set-up
    std::vector<int> blocking_factor{8};    // BLOCKING_FACTOR: per level (one value or max_level+1 values)
    std::vector<int> max_grid_size{32};     // MAX_GRID_SIZE: per level (one value or max_level+1 values)
    int n_error_buf = 1;         // N_ERROR_BUF: cells added around tagged cells
    int n_proper = 1;            // N_PROPER: proper-nesting width in coarse cells
    double grid_eff = 0.7;       // GRID_EFF: grid efficiency of the clustering
    int output_level_cap = -1;   // OUTPUT_LEVEL_CAP (L_out, IR-008 (b)); -1 = MAX_LEVEL
    PostRegridProjection post_regrid_projection = PostRegridProjection::Auto;  // POST_REGRID_PROJECTION: 'AUTO' (default, D-063), 'ON' or 'OFF'
    bool exact_sums = false;     // &MISC EXACT_SUMS: read so that the line is accepted; the AMR code ignores it (its sums follow the FDS-order default of D-028/D-053)
    std::vector<RegionSpec> regions;        // &AMR_REGION lines

    // Tagging (FR-011 / IR-003 working names). Unset = off (NaN / empty). TAG_DIST deferred.
    double tag_tmp_diff = std::numeric_limits<double>::quiet_NaN();   // TAG_TMP_DIFF (K), undivided
    double tag_tmp_above = std::numeric_limits<double>::quiet_NaN();  // TAG_TMP_ABOVE (K above ambient)
    double tag_rho_diff = std::numeric_limits<double>::quiet_NaN();   // TAG_RHO_DIFF (relative)
    std::string tag_spec_id;                                         // TAG_SPEC_ID (character); needs TAG_SPEC_DIFF
    double tag_spec_diff = std::numeric_limits<double>::quiet_NaN();  // TAG_SPEC_DIFF (undivided on Y)
    double tag_hrrpuv = std::numeric_limits<double>::quiet_NaN();     // TAG_HRRPUV (kW/m3)
    double tag_phi_diff = std::numeric_limits<double>::quiet_NaN();   // TAG_PHI_DIFF (undivided PHI jump: the level-set front)
    double tag_vort = std::numeric_limits<double>::quiet_NaN();       // TAG_VORT (|ω| undivided threshold)
    double tag_keep = 0.8;                                           // TAG_KEEP hysteresis factor in (0,1]

    static bool tag_set(double v) { return !std::isnan(v); }
    bool has_tag_tmp_diff() const { return tag_set(tag_tmp_diff); }
    bool has_tag_tmp_above() const { return tag_set(tag_tmp_above); }
    bool has_tag_rho_diff() const { return tag_set(tag_rho_diff); }
    bool has_tag_spec_diff() const { return tag_set(tag_spec_diff); }
    bool has_tag_hrrpuv() const { return tag_set(tag_hrrpuv); }
    bool has_tag_phi_diff() const { return tag_set(tag_phi_diff); }
    bool has_tag_vort() const { return tag_set(tag_vort); }
    bool has_any_tag_criterion() const
    {
        return has_tag_tmp_diff() || has_tag_tmp_above() || has_tag_rho_diff() || has_tag_hrrpuv() ||
               has_tag_vort() || has_tag_phi_diff() || (has_tag_spec_diff() && !tag_spec_id.empty());
    }

    // Value for level `lev` of a per-level vector (a single entry applies to every level).
    static int at(const std::vector<int>& v, int lev) { return v.size() == 1 ? v[0] : v[lev]; }
    int ref_ratio_at(int lev) const { return at(ref_ratio, lev); }          // ratio between lev and lev+1
    int blocking_factor_at(int lev) const { return at(blocking_factor, lev); }
    int max_grid_size_at(int lev) const { return at(max_grid_size, lev); }
    int output_level() const { return output_level_cap < 0 ? max_level : output_level_cap; }
};

// One parsed namelist group: name plus key -> list of value tokens (strings keep their content without quotes).
struct NamelistGroup {
    std::string name;
    std::map<std::string, std::vector<std::string>> values;
    int line = 0;
};

// Finds every &NAME ... / group in `text`. Text outside groups is ignored (as FDS does); scanning stops at &TAIL.
std::vector<NamelistGroup> scan_namelists(const std::string& text);

// Parses the &AMR and &AMR_REGION groups of `text` (and reads EXACT_SUMS from &MISC, logical, default F, without any other &MISC check) and validates them (IR-003, FR-010: ratios 2 and 4 only; ranges; list lengths;
// unknown names are errors). `present` is false and defaults are returned when there is no &AMR line.
AmrParams parse_amr_params(const std::string& text, Report& rep);

// Range and consistency checks that do not need the meshes (also called by parse_amr_params).
void validate_amr_params(const AmrParams& p, Report& rep);

}  // namespace fdsrt

#endif
