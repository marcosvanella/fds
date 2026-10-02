// RegridAmrCore.H: the AmrCore subclass of Role 3 (R2a). AmrCore is the AMReX class that owns the grids of all levels (BoxArray and
// DistributionMapping per level) and calls back to make, remake and clear a level and to tag cells for refinement.
//
// R2a builds the STATIC hierarchy: the grids come from fdsrt::Hierarchy (R1: finer &MESH boxes plus nesting boxes), not from tagging.
// Every created, changed or removed level is reported to a fdsrt::LevelListener (RegridInterface.H; the driver's LevelRegistry implements it).
// Dynamic regrid (tagging, transfer of data) is R4; ErrorEst tags nothing yet.
#ifndef FDSRT_REGRIDAMRCORE_H_
#define FDSRT_REGRIDAMRCORE_H_

#include <AMReX_AmrCore.H>

#include <functional>

#include "Hierarchy.H"
#include "RegridInterface.H"

namespace fdsrt {

// amrex::AmrInfo (ratios, blocking factors, max grid size, nesting width, per direction) from the hierarchy and the &AMR parameters.
amrex::AmrInfo make_amr_info(const Hierarchy& h, const AmrParams& p);
// Geometry of level 0 (physical box = bounding box of the level-0 meshes, periodicity from the hierarchy).
amrex::Geometry make_level0_geometry(const Hierarchy& h);
amrex::Box to_amrex(const IBox& b);

class RegridAmrCore : public amrex::AmrCore {
public:
    // Distribution of the grids of `level` over ranks. Default (empty function): AMReX space-filling-curve mapping, deterministic.
    using DmFn = std::function<amrex::DistributionMapping(int level, const amrex::BoxArray&)>;

    RegridAmrCore(const Hierarchy& h, const AmrParams& p);

    // Installs the static hierarchy: SetBoxArray / SetDistributionMap / SetFinestLevel and calls listener.make_level for levels 1..top, and
    // for level 0 when notify_level0 is true (the driver owns level 0 and adopts it, so the driver passes false and level0_dm = its mapping).
    void init_static(LevelListener& listener, bool notify_level0 = false, const DmFn& dm_fn = DmFn(),
                     const amrex::DistributionMapping* level0_dm = nullptr);

    LevelLayout layout(int level) const;     // what the registry gets for `level`
    const Hierarchy& hierarchy() const { return h_; }

protected:
    // AmrCore callbacks. R2a forwards to the listener; the data transfer of a real regrid is R4.
    void MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) override;
    void MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) override;
    void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) override;
    void ClearLevel(int lev) override;
    void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;   // R2a: tags nothing

private:
    Hierarchy h_;
    LevelListener* listener_ = nullptr;
};

}  // namespace fdsrt

#endif
