// mesh_text.H: test helper that reads &MESH and &MULT lines from FDS input text into MeshInput (mesh order: MULT copies with I fastest).
// Production code does not use this: the driver gets the meshes from the FDS set-up.
#ifndef FDSRT_TEST_MESH_TEXT_H_
#define FDSRT_TEST_MESH_TEXT_H_

#include <cstdlib>
#include <fstream>
#include <sstream>

#include "Hierarchy.H"

namespace fdsrt_test {

inline std::string read_file(const std::string& path)
{
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

inline std::vector<fdsrt::MeshInput> meshes_from_text(const std::string& text)
{
    using namespace fdsrt;
    auto groups = scan_namelists(text);
    auto num = [](const NamelistGroup& g, const char* k, double dflt) {
        auto it = g.values.find(k);
        return it == g.values.end() || it->second.empty() ? dflt : std::atof(it->second[0].c_str());
    };
    std::vector<MeshInput> out;
    for (const auto& g : groups) {
        if (g.name != "MESH") continue;
        MeshInput m;
        for (int d = 0; d < 3; ++d) m.ijk[d] = std::atoi(g.values.at("IJK")[d].c_str());
        for (int k = 0; k < 6; ++k) m.xb[k] = std::atof(g.values.at("XB")[k].c_str());
        m.rank = static_cast<int>(num(g, "MPI_PROCESS", 0));
        auto mi = g.values.find("MULT_ID");
        if (mi == g.values.end()) { out.push_back(m); continue; }
        for (const auto& mg : groups) {
            if (mg.name != "MULT" || mg.values.at("ID")[0] != mi->second[0]) continue;
            const double dd[3] = {num(mg, "DX", 0), num(mg, "DY", 0), num(mg, "DZ", 0)};
            const int lo[3] = {(int)num(mg, "I_LOWER", 0), (int)num(mg, "J_LOWER", 0), (int)num(mg, "K_LOWER", 0)};
            const int up[3] = {(int)num(mg, "I_UPPER", 0), (int)num(mg, "J_UPPER", 0), (int)num(mg, "K_UPPER", 0)};
            const char* sk[3][2] = {{"I_LOWER_SKIP", "I_UPPER_SKIP"}, {"J_LOWER_SKIP", "J_UPPER_SKIP"}, {"K_LOWER_SKIP", "K_UPPER_SKIP"}};
            bool has_skip = false;
            int sl[3], su[3];
            for (int d = 0; d < 3; ++d) {
                has_skip = has_skip || mg.values.count(sk[d][0]) || mg.values.count(sk[d][1]);
                sl[d] = (int)num(mg, sk[d][0], -1000000);
                su[d] = (int)num(mg, sk[d][1], 1000000);
            }
            for (int k = lo[2]; k <= up[2]; ++k)
                for (int j = lo[1]; j <= up[1]; ++j)
                    for (int i = lo[0]; i <= up[0]; ++i) {
                        const int ijk[3] = {i, j, k};
                        if (has_skip && ijk[0] >= sl[0] && ijk[0] <= su[0] && ijk[1] >= sl[1] && ijk[1] <= su[1] && ijk[2] >= sl[2] && ijk[2] <= su[2]) continue;
                        MeshInput c = m;
                        for (int d = 0; d < 3; ++d) { c.xb[2 * d] += ijk[d] * dd[d]; c.xb[2 * d + 1] += ijk[d] * dd[d]; }
                        out.push_back(c);
                    }
        }
    }
    return out;
}

}  // namespace fdsrt_test

#endif
