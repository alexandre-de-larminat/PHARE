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
        , maxwellampere_{dict["maxwell_ampere"]}
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
                   
    void MF_(level_t& level, PICModel& model,  Messenger& fromCoarser, double const currentTime, double const newTime);

    void MA_(level_t& level, PICModel& model,  Messenger& fromCoarser, double const currentTime, double const newTime);

    void moveFermions_(level_t& level, PICModel& model, Electromag& electromag,
                       ResourcesManager& rm, Messenger& fromCoarser, double const currentTime,
                       double const newTime);

    void averageAndSet_(level_t& level, PICModel& model, Messenger& fromCoarser);


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


// TODO: is this needed?
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& PICInfo = dynamic_cast<amr::PICMessengerInfo&>(*info);

    auto const& Bavg  = electromagAvg_.B;
    //auto const& E     = electromag.E;

    //PICInfo.ghostElectric.emplace_back(core::VecFieldNames{E});
    PICInfo.initMagnetic.emplace_back(core::VecFieldNames{Bavg});
}


template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                                                  int const levelNumber, IPhysicalModel_t& model,
                                                  IMessenger& fromCoarserMessenger,
                                                  double const currentTime, double const newTime)
{
    printf("advanceLevel\n");
    PHARE_LOG_SCOPE("SolverPIC::advanceLevel");

    auto& PICmodel         = dynamic_cast<PICModel&>(model);
    auto& fromCoarser      = dynamic_cast<PICMessenger&>(fromCoarserMessenger);
    auto& resourcesManager = *PICmodel.resourcesManager;
    auto level             = hierarchy->getPatchLevel(levelNumber);

    printf("restartJ_\n");
    restartJ_(*level, PICmodel, fromCoarser, currentTime);

    // Solve Faraday 
    printf("MF_\n");
    MF_(*level, PICmodel, fromCoarser, currentTime, newTime);

    printf("averageAndSet_\n");
    // Set Eavg to E, Bnew to B and calculate Bavg
    averageAndSet_(*level, PICmodel, fromCoarser);
    printf("done\n");

    printf("moveFermions_\n");
    // Push particles, project current onto the grid
    moveFermions_(*level, PICmodel, electromagAvg_, resourcesManager, fromCoarser, currentTime,
                  newTime);

    printf("MA_\n");
   // Solve MaxwellAmpere
    MA_(*level, PICmodel, fromCoarser, currentTime, newTime);


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
}
  


// Solve Maxell-Faraday and Maxwell-Ampere
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::MA_(level_t& level, PICModel& model, Messenger& fromCoarser, 
                                           double const currentTime, double const newTime)
{
    auto& PICState         = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto& electromag       = PICState.electromag;
    auto levelNumber       = level.getLevelNumber();

    {
        PHARE_LOG_SCOPE("SolverPIC::MA_");

        auto& B    = electromag.B;
        auto& E     = electromag.E;
        auto& J     = PICState.J;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, B, E, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, maxwellampere_);
            maxwellampere_(B, E, J, E, dt);

            resourcesManager->setTime(E, *patch, newTime);
        }
        printf("fillElectricGhosts\n");
        fromCoarser.fillElectricGhosts(E, levelNumber, newTime);

    }
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::MF_(level_t& level, PICModel& model, Messenger& fromCoarser, 
                                           double const currentTime, double const newTime)
{
    auto& PICState         = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto& electromag       = PICState.electromag;

    {
        PHARE_LOG_SCOPE("SolverPIC::MF_");

        auto& B    = electromag.B;
        auto& E    = electromag.E;
        auto& Bnew = electromagNew_.B;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, B, Bnew, E);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, E, Bnew, dt);

            resourcesManager->setTime(Bnew, *patch, newTime);
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

    PHARE_DEBUG_DO(std::size_t nbrDomainParticles        = 0;
                   std::size_t nbrPatchGhostParticles    = 0;
                   std::size_t nbrLevelGhostNewParticles = 0;
                   std::size_t nbrLevelGhostOldParticles = 0;
                   std::size_t nbrLevelGhostParticles    = 0;

                    for (auto& patch : level) {
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

   
    fromCoarser.fillGhostParticles(ions, electrons, level, newTime);
    fromCoarser.fillPopMomentGhosts(ions, electrons, level, newTime);

    for (auto& patch : level)
    {
        auto _      = rm.setOnPatch(*patch, electromag, ions, electrons);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        fermionUpdater_.updateFermions(ions, electrons, layout);

        // no need to update time, since it has been done before
    }

    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    fromCoarser.fillParticleMomentGhosts(ions, electrons, level, newTime);
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::averageAndSet_(level_t& level, PICModel& model,
                                                 Messenger& fromCoarser)
{
    PHARE_LOG_SCOPE("SolverPIC::averageAndSet_");

    auto& resourcesManager = model.resourcesManager;

    auto& Bnew  = electromagNew_.B;
    auto& B     = model.state.electromag.B;
    auto& Bavg  = electromagAvg_.B;
    auto& E     = model.state.electromag.E;
    auto& Eavg  = electromagAvg_.E;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, Bnew, B, Bavg, E, Eavg);
        PHARE::core::average(B, Bnew, Bavg);
        PHARE::core::average(Bnew, Bnew, B);
        PHARE::core::average(E, E, Eavg); // Unrelated to Faraday, we are just using the loop for convenience, to deal with E and B as a unit.
    }
}

} // namespace PHARE::solver

#endif