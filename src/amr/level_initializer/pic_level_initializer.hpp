#ifndef PHARE_PIC_LEVEL_INITIALIZER_HPP
#define PHARE_PIC_LEVEL_INITIALIZER_HPP

#include "amr/level_initializer/level_initializer.hpp"
#include "amr/messengers/pic_messenger.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/physical_models/pic_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/fermions/fermions.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
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

        inline bool isRootLevel(int levelNumber) const { return levelNumber == 0; }

    public:
        explicit PICLevelInitializer(PHARE::initializer::PHAREDict const& dict)
        {
        }
        virtual void initialize(std::shared_ptr<hierarchy_t> const& hierarchy, int levelNumber,
                                std::shared_ptr<level_t> const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) override
        {
            printf("PICLevelInitializer::initialize\n");
            core::Interpolator<dimension, interp_order> interpolate_;
            auto& picModel = static_cast<PICModel&>(model);
            auto& level       = amr_types::getLevel(*hierarchy, levelNumber);

            auto& picMessenger = dynamic_cast<PICMessenger&>(messenger);

            if (isRootLevel(levelNumber))
            {
                printf("Is root level\n");
                PHARE_LOG_START("picLevelInitializer::initialize : root level init");
                printf("Logging started; initializing\n");
                model.initialize(level);
                printf("Root level initialized\n");
                messenger.fillRootGhosts(model, level, initDataTime);
                printf("Root ghosts filled\n");
                PHARE_LOG_STOP("picLevelInitializer::initialize : root level init");
                std::cout << "Initialized level " << levelNumber << "\n";
            }

            else
            {
                printf("Not root level\n");
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
            printf("Decalring particles\n");
            auto& ions             = picModel.state.ions;
            auto& electrons        = picModel.state.pic_electrons;

            for (auto& patch : level)
            {
                printf("Computing moments\n");
                auto& resourcesManager = picModel.resourcesManager;
                auto dataOnPatch       = resourcesManager->setOnPatch(*patch, ions, electrons);
                auto layout            = amr::layoutFromPatch<GridLayoutT>(*patch);

                core::resetMoments(ions); //TODO add electrons
                core::depositParticles(ions, layout, interpolate_, core::DomainDeposit{});
                 core::depositParticles(ions, layout, interpolate_, core::PatchGhostDeposit{});

                if (!isRootLevel(levelNumber))
                {
                    core::depositParticles(ions, layout, interpolate_, core::LevelGhostDeposit{});
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
            // is not needed. But is still seems to use the messenger tempoeraries like
            // NiOld etc. so prepareStep() must be called, see end of the function.
            
            //picMessenger.fillIonMomentGhosts(ions, level, initDataTime);


            // now moments are known everywhere, compute J and E
            // via Ampere and Ohm
            // this only needs to be done for the root level
            // since otherwise initLevel has done it already

            if (isRootLevel(levelNumber))
            {
                printf("Computing E\n");
                auto& E = picModel.state.electromag.E;

                for (auto& patch : level)
                {
                    auto _      = picModel.resourcesManager->setOnPatch(*patch, E);
                    E.zero();

                    picModel.resourcesManager->setTime(E, *patch, 0.);
                }

                //picMessenger.fillElectricGhosts(E, levelNumber, 0.);
            }

            // quantities have been computed on the level, like the moments and J
            // that we later in the code need to get on level ghost nodes via
            // space and TIME interpolation. We thus need to save current values
            // in "old" messenger temporaries.
            // NOTE :  this may probably be skipped for finest level since, TBC at some point
            //picMessenger.prepareStep(picModel, level, initDataTime);
        }
    };
} // namespace PHARE::solver

#endif
