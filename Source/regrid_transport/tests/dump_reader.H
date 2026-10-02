// dump_reader.H: reader for the reference dumps of the instrumented FDS build (format: Source/driver/README.md, "Reference dump"). Test helper only.
#ifndef FDSRT_TEST_DUMP_READER_H_
#define FDSRT_TEST_DUMP_READER_H_

#include <cstdint>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fdsrt_test {

struct DumpArr {
    int rank = 0, lb[4] = {0, 0, 0, 0}, ub[4] = {0, 0, 0, 0};
    std::vector<double> d;
    bool has(int i, int j, int k, int n = 1) const { return i >= lb[0] && i <= ub[0] && j >= lb[1] && j <= ub[1] && k >= lb[2] && k <= ub[2] && n >= lb[3] && n <= ub[3]; }
    double at(int i, int j, int k, int n = 1) const
    {
        const long e0 = ub[0] - lb[0] + 1, e1 = ub[1] - lb[1] + 1, e2 = ub[2] - lb[2] + 1;
        return d[(i - lb[0]) + e0 * ((j - lb[1]) + e1 * ((k - lb[2]) + e2 * (long)(n - lb[3])))];
    }
};
struct DumpRec {
    std::string name;
    int icyc = 0, pred = 0, first = 0;
    std::map<std::string, DumpArr> bef;
};

inline std::string trim16(const char* p) { std::string s(p, 16); while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back(); return s; }

// All records of one dump file (one mesh). Throws on a malformed file.
inline std::vector<DumpRec> read_dump(const std::string& path)
{
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[8];
    int32_t nrec = 0;
    if (std::fread(magic, 1, 8, f) != 8 || std::string(magic, 8) != "FDSDUMP1" || std::fread(&nrec, 4, 1, f) != 1) throw std::runtime_error("bad dump header " + path);
    std::vector<DumpRec> out;
    for (;;) {
        char nm[16];
        if (std::fread(nm, 1, 16, f) != 16) break;
        int32_t h4[4];
        double td[8];
        int32_t nbk[2];
        if (std::fread(h4, 4, 4, f) != 4 || std::fread(td, 8, 8, f) != 8 || std::fread(nbk, 4, 2, f) != 2) throw std::runtime_error("truncated record " + path);
        DumpRec r;
        r.name = trim16(nm); r.icyc = h4[0]; r.pred = h4[1]; r.first = h4[2];
        for (int part = 0; part < 2; ++part) {
            for (int a = 0; a < nbk[part]; ++a) {
                char an[16];
                int32_t rank, lb[4], ub[4];
                if (std::fread(an, 1, 16, f) != 16 || std::fread(&rank, 4, 1, f) != 1 || std::fread(lb, 4, 4, f) != 4 || std::fread(ub, 4, 4, f) != 4) throw std::runtime_error("truncated array " + path);
                DumpArr arr;
                arr.rank = rank;
                long n = 1;
                for (int i = 0; i < 4; ++i) { arr.lb[i] = lb[i]; arr.ub[i] = ub[i]; if (i < rank) n *= (ub[i] - lb[i] + 1); }
                // arrays of lower rank carry lb=ub=1 in the unused slots (see PUT): treat them as extent 1
                for (int i = rank; i < 4; ++i) { arr.lb[i] = lb[i]; arr.ub[i] = ub[i]; }
                arr.d.resize(n);
                if (std::fread(arr.d.data(), 8, n, f) != static_cast<size_t>(n)) throw std::runtime_error("truncated data " + path);
                if (part == 0) r.bef[trim16(an)] = std::move(arr);
            }
        }
        out.push_back(std::move(r));
        if (nrec > 0 && static_cast<int>(out.size()) >= nrec) break;
    }
    std::fclose(f);
    return out;
}

}  // namespace fdsrt_test

#endif
