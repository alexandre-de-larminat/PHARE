#ifndef PHARE_SOLVER_PIC_HPP
#define PHARE_SOLVER_PIC_HPP


#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/Patch.h>
#include <cstddef>
#include <iomanip>
#include <sstream>


#include "amr/messengers/pic_messenger.hpp"
#include "amr/messengers/pic_messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "amr/physical_models/physical_model.hpp"
#include "amr/physical_models/pic_model.hpp"
#include "amr/solvers/solver.hpp"

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/hybrid/hybrid_quantities.hpp"


#include "core/numerics/fermion_updater/fermion_updater.hpp"
#include "core/numerics/ampere/maxwell_ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

namespace PHARE::solver
{

template<typename PICModel, typename AMR_Types>
class SolverPIC : public ISolver<AMR_Types>
{
private:
    using Electromag       = typename PICModel::electromag_type;
    using Ions             = typename PICModel::ions_type;
    using PICElectrons     = typename PICModel::pic_electrons_type;
    using ParticleArray    = typename Ions::particle_array_type;
    using VecFieldT        = typename PICModel::vecfield_type;
    using GridLayout       = typename PICModel::gridlayout_type;
    using ResourcesManager = typename PICModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using PICMessenger     = amr::PICMessenger<PICModel>;

    Electromag electromagNew_{"EMNew"};
    Electromag electromagAvg_{"EMAvg"};

    // are these ever called elsewhere? maybe delete
    static constexpr auto dimension    = PICModel::dimension;
    static constexpr auto interp_order = PICModel::gridlayout_type::interp_order;

    PHARE::core::Faraday<GridLayout> faraday_;
    PHARE::core::MaxwellAmpere<GridLayout> maxwellampere_;
    PHARE::core::FermionUpdater<Ions, PICElectrons, Electromag, VecFieldT, GridLayout> fermionUpdater_;


public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

    SolverPIC(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"PICSolver"}
        , fermionUpdater_{dict["ion_updater"]}
    {
    }


    virtual ~SolverPIC() = default;

    virtual std::string modelName() const override { return PICModel::model_name; }


    virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    virtual void registerResources(IPhysicalModel_t& model) override;


    virtual void allocate(IPhysicalModel_t& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override;
                        

    virtual void advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy, int const levelNumber,
                              IPhysicalModel_t& model, IMessenger& fromCoarserMessenger,
                              double const currentTime, double const newTime) override;

private:
    using Messenger = PICMessenger;

    void restartJ_(level_t& level, PICModel& model, Messenger& fromCoarser, double const currentTime);
                   
    void MAMF_(level_t& level, PICModel& model,  Electromag& electromag, Messenger& fromCoarser, double const currentTime, double const newTime);

    void MF_(level_t& level, PICModel& model, Messenger& fromCoarser, double const currentTime, double const newTime);

    void moveFermions_(level_t& level, PICModel& model, Electromag& electromag,
                       ResourcesManager& rm, Messenger& fromCoarser, double const currentTime,
                       double const newTime);

    void average_(level_t& level, PICModel& model, Messenger& fromCoarser);

    void saveState_(level_t& level, PICModel& model, ResourcesManager& rm);

    void restoreState_(level_t& level, PICModel& model, ResourcesManager& rm);

    // extend lifespan
    std::unordered_map<std::string, ParticleArray> tmpDomain;
    std::unordered_map<std::string, ParticleArray> patchGhost;


}; // end SolverPIC


template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::registerResources(IPhysicalModel_t& model)
{
    auto& pmodel = dynamic_cast<PICModel&>(model);
    pmodel.resourcesManager->registerResources(electromagNew_);
    pmodel.resourcesManager->registerResources(electromagAvg_);
}



template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::allocate(IPhysicalModel_t& model,
                                                 SAMRAI::hier::Patch& patch,
                                                 double const allocateTime) const
{
    auto& pmodel = dynamic_cast<PICModel&>(model);
    pmodel.resourcesManager->allocate(electromagNew_, patch, allocateTime);
    pmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);
}


// TODO: this stays for now, but should probably be edited
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    printf("SolverPIC::fillMessengerInfo\n");
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    auto const& Eavg  = electromagAvg_.E;
    auto const& Bnew = electromagNew_.B;

    //hybridInfo.ghostElectric.emplace_back(core::VecFieldNames{Eavg});
    //hybridInfo.initMagnetic.emplace_back(core::VecFieldNames{Bnew});
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                                                  int const levelNumber, IPhysicalModel_t& model,
                                                  IMessenger& fromCoarserMessenger,
                                                  double const currentTime, double const newTime)
{
    //printf("advanceLevel\n");
    PHARE_LOG_SCOPE("SolverPIC::advanceLevel");

    auto& PICmodel         = dynamic_cast<PICModel&>(model);
    auto& fromCoarser      = dynamic_cast<PICMessenger&>(fromCoarserMessenger);
    auto& resourcesManager = *PICmodel.resourcesManager;
    auto level             = hierarchy->getPatchLevel(levelNumber);
    auto& electromag       = PICmodel.state.electromag;

    //printf("restartJ\n");
    restartJ_(*level, PICmodel, fromCoarser, currentTime);

    //MF_(*level, PICmodel, fromCoarser, currentTime, newTime);
    //printf("saveState\n");
    saveState_(*level, PICmodel, resourcesManager);

    //printf("moveFermions\n");
    // Push particles, project current onto the grid
    moveFermions_(*level, PICmodel, electromag, resourcesManager, fromCoarser, currentTime,
                  newTime);

    //printf("MAMF\n");
   // Solve Faraday, MaxwellAmpere
    MAMF_(*level, PICmodel, electromag, fromCoarser, currentTime, newTime);

    //printf("restoreState\n");
    restoreState_(*level, PICmodel, resourcesManager);

    //printf("restart2\n");
    restartJ_(*level, PICmodel, fromCoarser, currentTime);

    //printf("move2\n");
    moveFermions_(*level, PICmodel, electromagAvg_, resourcesManager, fromCoarser, currentTime,
                  newTime);

    //printf("MAMF2\n");
    MAMF_(*level, PICmodel, electromagAvg_, fromCoarser, currentTime, newTime);

    //printf("average\n");
    // Set Bnew to B and Enew to E
    average_(*level, PICmodel, fromCoarser);

}



// Restart current density (there's probably a better way to do this)
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::restartJ_(level_t& level, PICModel& model,
                                              Messenger& fromCoarser, double const currentTime)
{
    auto& PICState         = model.state;
    auto& resourcesManager = model.resourcesManager;

    auto& J = PICState.J;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, J);
        J.zero();
    }
    // TODO: check whether this is necessary
    // fromCoarser.fillCurrentGhosts(J, level.getLevelNumber(), currentTime);
}


template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::MF_(level_t& level, PICModel& model, Messenger& fromCoarser, 
                                           double const currentTime, double const newTime)
{
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;

    {
        PHARE_LOG_SCOPE("SolverPIC::MF_");

        auto& Bavg = electromagAvg_.B;
        auto& E    = model.state.electromag.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bavg, E);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(Bavg, E, Bavg, dt);
            resourcesManager->setTime(Bavg, *patch, newTime);
        }
        //fromCoarser.fillCurrentGhosts(E, level.getLevelNumber(), newTime);
    }
}    


template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::saveState_(level_t& level, PICModel& model, ResourcesManager& rm)
{
    auto& ions = model.state.ions;
    auto& electrons = model.state.pic_electrons;

    for (auto& patch : level)
    {
        std::stringstream ss;
        ss << patch->getGlobalId();

        auto _ = rm.setOnPatch(*patch, ions, electrons);
        for (auto& pop : ions)
        {
            auto key = ss.str() + "_" + pop.name();
            if (!tmpDomain.count(key))
                tmpDomain.emplace(key, pop.domainParticles());
            else
                tmpDomain.at(key) = pop.domainParticles();
            if (!patchGhost.count(key))
                patchGhost.emplace(key, pop.patchGhostParticles());
            else
                patchGhost.at(key) = pop.patchGhostParticles();
        }
        for (auto& pop : electrons)
        {
            auto key = ss.str() + "_" + pop.name();
            if (!tmpDomain.count(key))
                tmpDomain.emplace(key, pop.domainParticles());
            else
                tmpDomain.at(key) = pop.domainParticles();
            if (!patchGhost.count(key))
                patchGhost.emplace(key, pop.patchGhostParticles());
            else
                patchGhost.at(key) = pop.patchGhostParticles();
        }
    
    }
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::restoreState_(level_t& level, PICModel& model,
                                                      ResourcesManager& rm)
{
    auto& ions = model.state.ions;
    auto& electrons = model.state.pic_electrons;

    for (auto& patch : level)
    {
        std::stringstream ss;
        ss << patch->getGlobalId();

        auto _ = rm.setOnPatch(*patch, ions, electrons);
        for (auto& pop : ions)
        {
            pop.domainParticles()     = std::move(tmpDomain.at(ss.str() + "_" + pop.name()));
            pop.patchGhostParticles() = std::move(patchGhost.at(ss.str() + "_" + pop.name()));
        }
        for (auto& pop : electrons)
        {
            pop.domainParticles()     = std::move(tmpDomain.at(ss.str() + "_" + pop.name()));
            pop.patchGhostParticles() = std::move(patchGhost.at(ss.str() + "_" + pop.name()));
        }
    }
}

// Solve Maxell-Faraday and Maxwell-Ampere
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::MAMF_(level_t& level, PICModel& model, Electromag& electromag_, Messenger& fromCoarser, 
                                           double const currentTime, double const newTime)
{
    auto& PICState         = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto& electromag       = PICState.electromag;

    {
        PHARE_LOG_SCOPE("SolverPIC::MAMF_.maxwellampere");

        auto& B_    = electromag_.B;
        auto& E     = electromag.E;
        auto& Enew  = electromagNew_.E;
        auto& Eavg  = electromagAvg_.E;
        auto& J     = PICState.J;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, B_, E, Enew, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, maxwellampere_);
            maxwellampere_(B_, E, J, Enew, dt);

            resourcesManager->setTime(Enew, *patch, newTime);
        }
        //fromCoarser.fillCurrentGhosts(E, level.getLevelNumber(), newTime);

        for (auto& patch : level)
        {
            auto _ = resourcesManager->setOnPatch(*patch, E, Enew, Eavg);
            PHARE::core::average(E, Enew, Eavg);
        }
    }

    {
        PHARE_LOG_SCOPE("SolverPIC::MAMF_.faraday");

        auto& Bnew  = electromagNew_.B;
        auto& B     = electromag.B;
        auto& Bavg  = electromagAvg_.B;
        auto& Eavg  = electromagAvg_.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bnew, B, Eavg);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, Eavg, Bnew, dt);
            resourcesManager->setTime(Bnew, *patch, newTime);
        }
        
        for (auto& patch : level)
        {
            auto _ = resourcesManager->setOnPatch(*patch, B, Bnew, Bavg);
            PHARE::core::average(B, Bnew, Bavg);
        }
    }
}

/* 
TODO: electrons should be dynamically 
created with the level and deleted with it.
*/
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::moveFermions_(level_t& level, PICModel& model, 
                                                  Electromag& electromag, ResourcesManager& rm,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime)
{
    PHARE_LOG_SCOPE("SolverPIC::moveFermions_");
    auto& ions = model.state.ions;
    auto& electrons = model.state.pic_electrons;

    PHARE_DEBUG_DO(std::size_t nbrDomainParticles = 0; std::size_t nbrPatchGhostParticles = 0;
                   std::size_t nbrLevelGhostNewParticles                                  = 0;
                   std::size_t nbrLevelGhostOldParticles                                  = 0;
                   std::size_t nbrLevelGhostParticles = 0; for (auto& patch
                                                                : level) {
                       auto _ = rm.setOnPatch(*patch, ions, electrons);

                       for (auto& pop : ions)
                        {
                           nbrDomainParticles += pop.domainParticles().size();
                           nbrPatchGhostParticles += pop.patchGhostParticles().size();
                           nbrLevelGhostNewParticles += pop.levelGhostParticlesNew().size();
                           nbrLevelGhostOldParticles += pop.levelGhostParticlesOld().size();
                           nbrLevelGhostParticles += pop.levelGhostParticles().size();
                           nbrPatchGhostParticles += pop.patchGhostParticles().size();

                           if (nbrLevelGhostOldParticles < nbrLevelGhostParticles
                               and nbrLevelGhostOldParticles > 0)
                               throw std::runtime_error(
                                   "Error - there are less old level ghost particles ("
                                   + std::to_string(nbrLevelGhostOldParticles) + ") than pushable ("
                                   + std::to_string(nbrLevelGhostParticles) + ")");
                       }
                        for (auto& pop : electrons)
                        {
                           nbrDomainParticles += pop.domainParticles().size();
                           nbrPatchGhostParticles += pop.patchGhostParticles().size();
                           nbrLevelGhostNewParticles += pop.levelGhostParticlesNew().size();
                           nbrLevelGhostOldParticles += pop.levelGhostParticlesOld().size();
                           nbrLevelGhostParticles += pop.levelGhostParticles().size();
                           nbrPatchGhostParticles += pop.patchGhostParticles().size();

                           if (nbrLevelGhostOldParticles < nbrLevelGhostParticles
                               and nbrLevelGhostOldParticles > 0)
                               throw std::runtime_error(
                                   "Error - there are less old level ghost particles ("
                                   + std::to_string(nbrLevelGhostOldParticles) + ") than pushable ("
                                   + std::to_string(nbrLevelGhostParticles) + ")");
                       }
                       
                   })


    auto dt        = newTime - currentTime;
    auto& PICState = model.state;
    auto& J        = PICState.J;

    for (auto& patch : level)
    {
        auto _ = rm.setOnPatch(*patch, electromag, ions, electrons, J);

        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        fermionUpdater_.updatePopulations(ions, electrons, electromag, J, layout, dt);

        // this needs to be done before calling the messenger
        rm.setTime(ions, *patch, newTime);
        rm.setTime(electrons, *patch, newTime);
        rm.setTime(J, *patch, newTime);
    }

    // Only ion ghosts are filled since electrons don't exist on the coarser level
    //fromCoarser.fillIonGhostParticles(ions, level, newTime);
    //fromCoarser.fillIonPopMomentGhosts(ions, level, newTime);

    for (auto& patch : level)
    {
        auto _      = rm.setOnPatch(*patch, electromag, ions, electrons);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        fermionUpdater_.updateFermions(ions, electrons, layout);

        // no need to update time, since it has been done before
    }

    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    //fromCoarser.fillIonMomentGhosts(ions, level, newTime);
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::average_(level_t& level, PICModel& model,
                                                 Messenger& fromCoarser)
{
    PHARE_LOG_SCOPE("SolverPIC::average_");

    auto& electromag      = model.state.electromag;
    auto& resourcesManager = model.resourcesManager;

    auto& Bnew  = electromagNew_.B;
    auto& B     = electromag.B;
    auto& Enew  = electromagNew_.E;    
    auto& E     = electromag.E;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, Bnew, B, E, Enew);
        PHARE::core::average(Bnew, Bnew, B);
        PHARE::core::average(Enew, Enew, E);
    }
}

} // namespace PHARE::solver

#endif