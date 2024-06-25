#ifndef PHARE_PIC_LEVEL_INITIALIZER_HPP
#define PHARE_PIC_LEVEL_INITIALIZER_HPP

#include "amr/level_initializer/level_initializer.hpp"
#include "amr/messengers/pic_messenger.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/physical_models/pic_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/moments/moments.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::solver
{
    template<typename PICModel>
    class PICLevelInitializer : public LevelInitializer<typename PICModel::amr_types>
    {
        using amr_types                    = typename PICModel::amr_types;
        using hierarchy_t                  = typename amr_types::hierarchy_t;
        using level_t                      = typename amr_types::level_t;
        using patch_t                      = typename amr_types::patch_t;
        using IPhysicalModelT              = IPhysicalModel<amr_types>;
        using IMessengerT                  = amr::IMessenger<IPhysicalModelT>;
        using PICMessenger                 = amr::PICMessenger<PICModel>;
        using GridLayoutT                  = typename PICModel::gridlayout_type;
        static constexpr auto dimension    = GridLayoutT::dimension;
        static constexpr auto interp_order = GridLayoutT::interp_order;

        PHARE::core::Ohm<GridLayoutT> ohm_;

        inline bool isRootLevel(int levelNumber) const { return levelNumber == 0; }

    public:
        explicit PICLevelInitializer(PHARE::initializer::PHAREDict const& dict)
        : ohm_{dict["algo"]["ohm"]}
        {
        }
        virtual void initialize(std::shared_ptr<hierarchy_t> const& hierarchy, int levelNumber,
                                std::shared_ptr<level_t> const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) override
        {
            core::Interpolator<dimension, interp_order> interpolate_;
            auto& picModel = static_cast<PICModel&>(model);
            auto& level       = amr_types::getLevel(*hierarchy, levelNumber);

            auto& picMessenger = dynamic_cast<PICMessenger&>(messenger);

            if (isRootLevel(levelNumber))
            {
                PHARE_LOG_START("picLevelInitializer::initialize : root level init");
                model.initialize(level);
                messenger.fillRootGhosts(model, level, initDataTime);
                PHARE_LOG_STOP("picLevelInitializer::initialize : root level init");
            }

            else
            {
                if (isRegridding)
                {
                    std::cout << "regriding level " << levelNumber << "\n";
                    PHARE_LOG_START("picLevelInitializer::initialize : regriding block");
                    messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
                    PHARE_LOG_STOP("picLevelInitializer::initialize : regriding block");
                }
                else
                {
                    PHARE_LOG_START("picLevelInitializer::initialize : initlevel");
                    messenger.initLevel(model, level, initDataTime);
                    PHARE_LOG_STOP("picLevelInitializer::initialize : initlevel");
                }
            }

            // now all particles are here
            // we must compute moments.
            

            auto& ions             = picModel.state.ions;
            auto& electrons        = picModel.state.pic_electrons;

            for (auto& patch : level)
            {
                auto& resourcesManager = picModel.resourcesManager;
                auto dataOnPatch       = resourcesManager->setOnPatch(*patch, ions, electrons);
                auto layout            = amr::layoutFromPatch<GridLayoutT>(*patch);

                core::resetMoments(ions, electrons);
                core::depositParticles(ions, electrons, layout, interpolate_, core::DomainDeposit{});
                core::depositParticles(ions, electrons, layout, interpolate_, core::PatchGhostDeposit{});

                if (!isRootLevel(levelNumber))
                {
                    core::depositParticles(ions, electrons, layout, interpolate_, core::LevelGhostDeposit{});
                }

                ions.computeDensity();
                ions.computeBulkVelocity();
                electrons.computeDensity();
                electrons.computeBulkVelocity();
            }

            // on level i>0, this relies on 'prepareStep' having been called on when
            // level i-1 was initialized (at the end of this function)
            // it seems SAMRAI does not call timeInterpolate() at this point although
            // both moments and J need time interpolation. It probably knows that
            // we are at a sync time across levels and that the time interpolation
            // is not needed. But is still seems to use the messenger temporaries like
            // NiOld etc. so prepareStep() must be called, see end of the function.
            
            picMessenger.fillParticleMomentGhosts(picModel.state.ions, picModel.state.pic_electrons, level, initDataTime);


            // now moments are known everywhere, compute J and E
            // via Ampere and Ohm
            // this only needs to be done for the root level
            // since otherwise initLevel has done it already

            if (isRootLevel(levelNumber))
            {
                auto& E = picModel.state.electromag.E;

                for (auto& patch : level)
                {
                    auto _      = picModel.resourcesManager->setOnPatch(*patch, E);
                    E.zero();

                    picModel.resourcesManager->setTime(E, *patch, 0.);
                }

                auto& B = picModel.state.electromag.B;


                for (auto& patch : level)
                {
                    printf("PIC level initializer: setOnPatch\n");
                    auto _      = picModel.resourcesManager->setOnPatch(*patch, E, B);
                    auto layout = PHARE::amr::layoutFromPatch<GridLayoutT>(*patch);
                    auto __     = core::SetLayout(&layout, ohm_);
                    
                    auto& Ex = E(core::Component::X);
                    auto& Ey = E(core::Component::Y);
                    auto& Ez = E(core::Component::Z);

                    auto& Bx = B(core::Component::X);
                    auto& By = B(core::Component::Y);
                    auto& Bz = B(core::Component::Z);

                    average(Bx, Bx, Ey);
                    average(By, By, Ez);
                    average(Bz, Bz, Ex);

                    picModel.resourcesManager->setTime(E, *patch, 0.);
                }

                if (electrons.isSettable() and 0)
                {
                    for (auto& patch : level)
                    {
                        auto layout = PHARE::amr::layoutFromPatch<GridLayoutT>(*patch);
                        auto _      = picModel.resourcesManager->setOnPatch(*patch, B, E, electrons);
                        auto& Ve    = electrons.velocity();
                        auto __     = core::SetLayout(&layout, ohm_);
                        ohm_(E, Ve, B);
                        picModel.resourcesManager->setTime(E, *patch, 0.);
                    }
                }

                picMessenger.fillElectricGhosts(E, levelNumber, 0.);
            }

            // quantities have been computed on the level, like the moments and J
            // that we later in the code need to get on level ghost nodes via
            // space and TIME interpolation. We thus need to save current values
            // in "old" messenger temporaries.
            // NOTE :  this may probably be skipped for finest level since, TBC at some point
            printf("PIC level initializer: prepareStep\n");
            picMessenger.prepareStep(picModel, level, initDataTime);
        }
    };
} // namespace PHARE::solver

#endif
