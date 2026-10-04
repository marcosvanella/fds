// test_species_avgdown.cpp: restriction of species is mass weighted (rho and rho*Z averaged, Z rebuilt), never a linear average of Z.
//  A  LevelOps::average_down_species on bare MultiFabs: the species mass of the covered coarse cells equals the fine species mass to round-off, sum Z = 1, 0 <= Z <= 1;
//     negative control: a linear average of Z (what average_down_cells on ZZ does) misses the mass by far more than the budget.
//  B  through the registry (LevelRegistry, average_down_registry with default names and the coarse-fine ghost hook step): same checks on the real field set.
//  C  agreement with the driver's RegistryTransfer::hierarchy_done (Role 1, S11.1) on an identical copy: RHO and ZZ of the covered cells, bitwise.
// Prints the measured numbers. Runs on 1 and 4 ranks.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

#include "DriverAdapter.H"
#include "LevelOps.H"
#include "RegistryTransfer.H"
#include "check.H"

using namespace fdsrt;

namespace {

constexpr int NS = 3;
// a "flame-like" state: density and the first mass fraction strongly correlated, sharp variation on the fine scale
void state(int i, int j, int k, double* rho, double* z)
{
    const double s = 0.5 + 0.5 * std::tanh(1.5 * std::sin(0.9 * i + 0.4 * j) + 1.2 * std::cos(0.7 * k - 0.3 * i));   // 0..1
    z[0] = 0.02 + 0.9 * s;
    z[1] = 0.5 * (1.0 - z[0]) * (0.6 + 0.3 * std::sin(1.1 * j + 0.5 * k));
    z[2] = 1.0 - z[0] - z[1];
    *rho = 1.2 / (1.0 + 4.0 * s);   // hot (light) where Z0 is large
}

void fill_fine(amrex::MultiFab& rho, amrex::MultiFab& zz)
{
    for (amrex::MFIter mfi(rho); mfi.isValid(); ++mfi) {
        auto r = rho.array(mfi); auto z = zz.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { double zv[NS], rv; state(i, j, k, &rv, zv); r(i, j, k) = rv; for (int n = 0; n < NS; ++n) z(i, j, k, n) = zv[n]; });
    }
}
void fill_coarse(amrex::MultiFab& rho, amrex::MultiFab& zz)
{
    for (amrex::MFIter mfi(rho); mfi.isValid(); ++mfi) {
        auto r = rho.array(mfi); auto z = zz.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { r(i, j, k) = 0.7 + 0.01 * i; for (int n = 0; n < NS; ++n) z(i, j, k, n) = (n == 0) ? 0.5 : (n == 1 ? 0.25 : 0.25); });
    }
}

struct Mass { long double fine[NS], crse_cov[NS], sumz_dev = 0, zmin = 1e30, zmax = -1e30; };

// species mass of the fine level and of the covered coarse cells (weights: fine cell volume = coarse volume / nref)
Mass measure(const amrex::MultiFab& rf, const amrex::MultiFab& zf, const amrex::MultiFab& rc, const amrex::MultiFab& zc, const amrex::iMultiFab& cov, double nref)
{
    Mass m{};
    for (int n = 0; n < NS; ++n) m.fine[n] = m.crse_cov[n] = 0;
    for (amrex::MFIter mfi(rf); mfi.isValid(); ++mfi) {
        auto r = rf.const_array(mfi); auto z = zf.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < NS; ++n) m.fine[n] += (long double)r(i, j, k) * z(i, j, k, n) / nref; });
    }
    for (amrex::MFIter mfi(rc); mfi.isValid(); ++mfi) {
        auto r = rc.const_array(mfi); auto z = zc.const_array(mfi); auto c = cov.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            if (!c(i, j, k)) return;
            long double s = 0;
            for (int n = 0; n < NS; ++n) { m.crse_cov[n] += (long double)r(i, j, k) * z(i, j, k, n); s += z(i, j, k, n); m.zmin = std::min<long double>(m.zmin, z(i, j, k, n)); m.zmax = std::max<long double>(m.zmax, z(i, j, k, n)); }
            m.sumz_dev = std::max<long double>(m.sumz_dev, std::abs(s - 1.0L));
        });
    }
    for (int n = 0; n < NS; ++n) {
        double v[2] = {(double)m.fine[n], (double)m.crse_cov[n]};   // sums are reduced in double: the relative budget is 1e-12
        amrex::ParallelAllReduce::Sum(v, 2, amrex::ParallelContext::CommunicatorSub());
        m.fine[n] = v[0]; m.crse_cov[n] = v[1];
    }
    double x[3] = {(double)m.sumz_dev, (double)m.zmax, (double)(-m.zmin)};
    amrex::ParallelDescriptor::ReduceRealMax(x, 3);
    m.sumz_dev = x[0]; m.zmax = x[1]; m.zmin = -x[2];
    return m;
}
double worst_rel(const Mass& m)
{
    double w = 0;
    for (int n = 0; n < NS; ++n) w = std::max(w, (double)(std::abs(m.crse_cov[n] - m.fine[n]) / m.fine[n]));
    return w;
}

void test_bare()
{
    amrex::Box cdom(amrex::IntVect(0), amrex::IntVect(15));
    amrex::RealBox rb({0, 0, 0}, {1, 1, 1});
    amrex::Array<int, 3> per{1, 1, 1};
    amrex::Geometry cg(cdom, rb, 0, per);
    const amrex::IntVect R(2);
    amrex::Geometry fg = amrex::refine(cg, R);
    amrex::BoxArray cba(cdom); cba.maxSize(8);
    amrex::DistributionMapping cdm(cba);
    amrex::BoxList bl; bl.push_back(amrex::refine(amrex::Box(amrex::IntVect(2, 3, 4), amrex::IntVect(8, 10, 11)), R));
    amrex::BoxArray fba(bl); fba.maxSize(8);
    amrex::DistributionMapping fdm(fba);
    amrex::iMultiFab cov = fdsamr::build_covered_mask(cba, cdm, fba, R);

    for (int mode = 0; mode < 2; ++mode) {
        amrex::MultiFab rf(fba, fdm, 1, 0), zf(fba, fdm, NS, 0), rc(cba, cdm, 1, 0), zc(cba, cdm, NS, 0);
        fill_fine(rf, zf); fill_coarse(rc, zc);
        if (mode == 0) average_down_species(rf, zf, rc, zc, cov, fg, cg, R);
        else {   // negative control: linear average of the mass fractions (and of rho)
            average_down_cells(zf, zc, fg, cg, R, 0, NS);
            average_down_cells(rf, rc, fg, cg, R, 0, 1);
        }
        const Mass m = measure(rf, zf, rc, zc, cov, 8.0);
        const double w = worst_rel(m);
        if (amrex::ParallelDescriptor::IOProcessor())
            std::printf("  A %-26s worst relative species-mass error of the covered cells %.2e, max|sum Z - 1| %.1e, Z in [%.3f, %.3f]\n", mode == 0 ? "rho*Z average (new)" : "linear Z average (control)", w, (double)m.sumz_dev,
                        (double)m.zmin, (double)m.zmax);
        if (mode == 0) {
            CHECK_MSG(w < 1e-14, "A: mass-weighted restriction conserves each species to round-off, got " + std::to_string(w));
            CHECK_MSG(m.sumz_dev < 1e-14 && m.zmin >= 0.0 && m.zmax <= 1.0, "A: Z sums to 1 and stays in [0,1]");
        } else {
            CHECK_MSG(w > 1e-5, "A: negative control, linear averaging of Z breaks the species mass (" + std::to_string(w) + ")");
            CHECK_MSG(m.sumz_dev < 1e-14, "A: negative control keeps sum Z = 1 (so a sum check alone would not catch it)");
        }
    }
}

// ---- registry level ----
struct Reg {
    fdsamr::Level0 l0;
    std::unique_ptr<fdsamr::Fields> F0;
    std::unique_ptr<fdsamr::SideData> sd0;
    std::unique_ptr<fdsamr::LevelRegistry> reg;
};

std::unique_ptr<Reg> make_reg()
{
    auto r = std::make_unique<Reg>();
    amrex::Vector<fdsamr::MeshInfo> meshes;
    fdsamr::DomainInfo dom{};
    dom.periodic[0] = dom.periodic[1] = dom.periodic[2] = 1;
    dom.n_tracked = NS; dom.n_total = NS; dom.nranks = amrex::ParallelDescriptor::NProcs();
    for (int q = 0; q < 2; ++q) {   // two level-0 boxes
        fdsamr::MeshInfo m{};
        m.ijk[0] = 8; m.ijk[1] = 16; m.ijk[2] = 16;
        m.xb[0] = q * 0.5; m.xb[1] = (q + 1) * 0.5; m.xb[2] = 0; m.xb[3] = 1; m.xb[4] = 0; m.xb[5] = 1;
        m.rank = q % dom.nranks;
        meshes.push_back(m);
    }
    r->l0 = fdsamr::assemble_level0(meshes, dom);
    r->F0 = std::make_unique<fdsamr::Fields>(r->l0, NS);
    r->sd0 = std::make_unique<fdsamr::SideData>(r->l0, fdsamr::layout_cell_walls(r->l0));
    r->reg = std::make_unique<fdsamr::LevelRegistry>(r->l0.dom, NS);
    r->reg->adopt_level0(r->l0, *r->F0, *r->sd0);
    fdsrt::LevelLayout fl;
    fl.level = 1; fl.ref_ratio_from_parent = amrex::IntVect(2);
    fl.geom = amrex::refine(r->l0.geom, fl.ref_ratio_from_parent);
    amrex::BoxList bl;
    bl.push_back(amrex::refine(amrex::Box(amrex::IntVect(2, 3, 4), amrex::IntVect(5, 10, 11)), fl.ref_ratio_from_parent));   // left of the level-0 box boundary
    bl.push_back(amrex::refine(amrex::Box(amrex::IntVect(6, 3, 4), amrex::IntVect(9, 10, 11)), fl.ref_ratio_from_parent));   // crosses it
    fl.ba.define(bl);
    fl.dm.define(fl.ba);
    r->reg->make_level(fl);
    return r;
}

void fill_reg(Reg& r)
{
    fill_coarse(r.reg->fields(0)["RHO"], r.reg->fields(0)["ZZ"]);
    fill_fine(r.reg->fields(1)["RHO"], r.reg->fields(1)["ZZ"]);
    for (int l = 0; l < 2; ++l) {
        fdsamr::Fields& F = r.reg->fields(l);
        amrex::MultiFab::Copy(F["RHOS"], F["RHO"], 0, 0, 1, 0);
        amrex::MultiFab::Copy(F["ZZS"], F["ZZ"], 0, 0, NS, 0);
        F["TMP"].setVal(300.0);
    }
}

void test_registry()
{
    auto a = make_reg();
    auto b = make_reg();
    fill_reg(*a); fill_reg(*b);
    const amrex::IntVect R(2);
    // B: through average_down_registry (default names: every transfer scalar)
    average_down_registry(*a->reg);
    {
        const Mass m = measure(a->reg->fields(1)["RHO"], a->reg->fields(1)["ZZ"], a->reg->fields(0)["RHO"], a->reg->fields(0)["ZZ"], *a->reg->covered_mask(0), 8.0);
        const Mass ms = measure(a->reg->fields(1)["RHOS"], a->reg->fields(1)["ZZS"], a->reg->fields(0)["RHOS"], a->reg->fields(0)["ZZS"], *a->reg->covered_mask(0), 8.0);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  B average_down_registry: species-mass error RHO/ZZ %.2e, RHOS/ZZS %.2e, max|sum Z - 1| %.1e\n", worst_rel(m), worst_rel(ms), (double)m.sumz_dev);
        CHECK_MSG(worst_rel(m) < 1e-14 && worst_rel(ms) < 1e-14, "B: registry restriction conserves the species mass of RHO/ZZ and RHOS/ZZS");
        CHECK_MSG(m.sumz_dev < 1e-14, "B: sum Z = 1 on the covered cells");
    }
    // B2: the hook's step 1 (it restricts before rebuilding the covered ghost layers)
    {
        auto c = make_reg(); fill_reg(*c);
        CfHookStats st;
        auto hook = make_cf_ghost_hook(*c->reg, ThermoProvider(), &st);
        fdsamr::CfGhostRequest rq; rq.level = 1; rq.fields = {"RHO", "ZZ"};
        hook(rq);
        // layers 1 and 2 of the covered cells next to a face are rebuilt by the FDS rule; deeper covered cells keep the step-1 restriction: check those only
        const amrex::MultiFab& r0 = c->reg->fields(0)["RHO"]; const amrex::MultiFab& z0 = c->reg->fields(0)["ZZ"];
        const amrex::MultiFab& r1 = c->reg->fields(1)["RHO"]; const amrex::MultiFab& z1 = c->reg->fields(1)["ZZ"];
        long n_deep = 0;
        CHECK_MSG(st.calls == 1, "B2: hook ran");
        // direct check of the whole-covered-region species mass restricted by step 1: layers >= 3 of the covered region are untouched by the FDS rule, so the integral over the
        // cells deeper than layer 2 must equal the fine integral over their children
        double mc[NS] = {0, 0, 0}, mf[NS] = {0, 0, 0};
        const amrex::Box deep(amrex::IntVect(5, 6, 7), amrex::IntVect(6, 7, 8));   // coarse cells at least 3 cells from the patch edge (patches span 2..9, 3..10, 4..11)
        for (amrex::MFIter mfi(r0); mfi.isValid(); ++mfi) {
            auto ra = r0.const_array(mfi); auto za = z0.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox() & deep, [&](int i, int j, int k) { ++n_deep; for (int n = 0; n < NS; ++n) mc[n] += ra(i, j, k) * za(i, j, k, n); });
        }
        for (amrex::MFIter mfi(r1); mfi.isValid(); ++mfi) {
            auto ra = r1.const_array(mfi); auto za = z1.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox() & amrex::refine(deep, R), [&](int i, int j, int k) { for (int n = 0; n < NS; ++n) mf[n] += ra(i, j, k) * za(i, j, k, n) / 8.0; });
        }
        amrex::ParallelAllReduce::Sum(mc, NS, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Sum(mf, NS, amrex::ParallelContext::CommunicatorSub());
        double w = 0; for (int n = 0; n < NS; ++n) w = std::max(w, std::abs(mc[n] - mf[n]) / mf[n]);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  B2 hook step 1, deep covered cells (%ld on this rank): species-mass error %.2e\n", n_deep, w);
        CHECK_MSG(w < 1e-14, "B2: the cf ghost hook restricts the species mass-weighted (deep covered cells), got " + std::to_string(w));
    }
    // C: Role 1's hierarchy_done on the identical copy
    {
        fdsamr::RegistryTransfer rt(*b->reg, NS);
        rt.hierarchy_done(true);
        double dr = 0, dz = 0;
        const amrex::iMultiFab& cov = *a->reg->covered_mask(0);
        for (amrex::MFIter mfi(a->reg->fields(0)["RHO"]); mfi.isValid(); ++mfi) {
            auto ra = a->reg->fields(0)["RHO"].const_array(mfi); auto rb = b->reg->fields(0)["RHO"].const_array(mfi);
            auto za = a->reg->fields(0)["ZZ"].const_array(mfi); auto zb = b->reg->fields(0)["ZZ"].const_array(mfi);
            auto zsa = a->reg->fields(0)["ZZS"].const_array(mfi); auto zsb = b->reg->fields(0)["ZZS"].const_array(mfi);
            auto c = cov.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (!c(i, j, k)) return;
                dr = std::max(dr, std::abs(ra(i, j, k) - rb(i, j, k)));
                for (int n = 0; n < NS; ++n) { dz = std::max(dz, std::abs(za(i, j, k, n) - zb(i, j, k, n))); dz = std::max(dz, std::abs(zsa(i, j, k, n) - zsb(i, j, k, n))); }
            });
        }
        double x[2] = {dr, dz}; amrex::ParallelDescriptor::ReduceRealMax(x, 2);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  C covered cells, Role 3 average_down_registry vs Role 1 RegistryTransfer::hierarchy_done: max|d RHO| %.1e, max|d ZZ/ZZS| %.1e\n", x[0], x[1]);
        CHECK_MSG(x[0] == 0.0 && x[1] == 0.0, "C: Role 3 and Role 1 restrictions agree bitwise on RHO, ZZ, ZZS of the covered cells");
        const Mass m = measure(b->reg->fields(1)["RHO"], b->reg->fields(1)["ZZ"], b->reg->fields(0)["RHO"], b->reg->fields(0)["ZZ"], cov, 8.0);
        CHECK_MSG(worst_rel(m) < 1e-14, "C: Role 1's restriction conserves the species mass too");
    }
}

}  // namespace

int main(int, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    {
        test_bare();
        test_registry();
    }
    const long nfail = fdstest::report("test_species_avgdown");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
