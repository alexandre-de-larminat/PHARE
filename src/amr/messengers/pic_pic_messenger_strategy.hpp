#ifndef PHARE_PIC_PIC_MESSENGER_STRATEGY_HPP
#define PHARE_PIC_PIC_MESSENGER_STRATEGY_HPP

#include "core/def/phare_mpi.hpp"

#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/IntVector.h"
#include "refiner_pool.hpp"
#include "synchronizer_pool.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/field/coarsening/magnetic_field_coarsener.hpp"
#include "amr/data/field/refine/field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_refiner.hpp"
#include "amr/data/field/refine/electric_field_refiner.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/pic_messenger_info.hpp"
#include "amr/messengers/pic_messenger_strategy.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/utilities/point/point.hpp"



#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>


#include <iterator>
#include <optional>
#include <utility>
#include <iomanip>
#include <iostream>
#include <string>


namespace PHARE
{
namespace amr
{

    /** \brief An PICMessenger is the specialization of a PICMessengerStrategy for PIC to
     * PIC data communications.
     */
    template<typename PICModel, typename RefinementParams>
    class PICPICMessengerStrategy : public PICMessengerStrategy<PICModel>
    {
        using IonsT                              = typename PICModel::ions_type;
        using PICElectronsT                      = typename PICModel::pic_electrons_type;
        using ElectromagT                        = typename PICModel::electromag_type;
        using VecFieldT                          = typename PICModel::vecfield_type;
        using GridLayoutT                        = typename PICModel::gridlayout_type;
        using FieldT                             = typename VecFieldT::field_type;
        using ResourcesManagerT                  = typename PICModel::resources_manager_type;
        static constexpr std::size_t dimension   = GridLayoutT::dimension;
        static constexpr std::size_t interpOrder = GridLayoutT::interp_order;
        using IPhysicalModel                     = typename PICModel::Interface;

        using InteriorParticleRefineOp = typename RefinementParams::InteriorParticleRefineOp;
        using CoarseToFineRefineOpOld  = typename RefinementParams::CoarseToFineRefineOpOld;
        using CoarseToFineRefineOpNew  = typename RefinementParams::CoarseToFineRefineOpNew;

        template<typename Policy>
        using BaseRefineOp          = FieldRefineOperator<GridLayoutT, FieldT, Policy>;
        using DefaultFieldRefineOp  = BaseRefineOp<DefaultFieldRefiner<dimension>>;
        using MagneticFieldRefineOp = BaseRefineOp<MagneticFieldRefiner<dimension>>;
        using ElectricFieldRefineOp = BaseRefineOp<ElectricFieldRefiner<dimension>>;
        using FieldTimeInterp       = FieldLinearTimeInterpolate<GridLayoutT, FieldT>;

        template<typename Policy>
        using BaseCoarsenOp     = FieldCoarsenOperator<GridLayoutT, FieldT, Policy>;
        using MagneticCoarsenOp = BaseCoarsenOp<MagneticFieldCoarsener<dimension>>;
        using DefaultCoarsenOp  = BaseCoarsenOp<DefaultFieldCoarsener<dimension>>;

    public:
        static const std::string stratName;
        static constexpr std::size_t rootLevelNumber = 0;


        PICPICMessengerStrategy(std::shared_ptr<ResourcesManagerT> manager,
                                      int const firstLevel)
            : PICMessengerStrategy<PICModel>{stratName}
            , resourcesManager_{std::move(manager)}
            , firstLevel_{firstLevel}
        {
            resourcesManager_->registerResources(NiOldUser_);
            resourcesManager_->registerResources(ViOld_);
            resourcesManager_->registerResources(NeOldUser_);
            resourcesManager_->registerResources(VeOld_);
        }

        virtual ~PICPICMessengerStrategy() = default;



        /* ------------------------------------------------------------------------
                        methods used for the IMessenger interface
           ------------------------------------------------------------------------ */


        /**
         * @brief allocate the messenger strategy internal variables to the model
         * resourceManager
         */
        void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            resourcesManager_->allocate(NiOldUser_, patch, allocateTime);
            resourcesManager_->allocate(ViOld_, patch, allocateTime);
            resourcesManager_->allocate(NeOldUser_, patch, allocateTime);
            resourcesManager_->allocate(VeOld_, patch, allocateTime);
        }



        /**
         * @brief setup creates all SAMRAI algorithms to communicate data involved in a
         * messenger between the coarse and fine levels.
         *
         * This method creates the SAMRAI algorithms for communications associated between pairs
         * of variables. The function does not create the SAMRAI schedules since they depend on
         * the levels
         */
        void registerQuantities([[maybe_unused]] std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<PICMessengerInfo> PICInfo{
                dynamic_cast<PICMessengerInfo*>(fromFinerInfo.release())};

            registerGhostComms_(PICInfo);
            registerInitComms(PICInfo);
            registerSyncComms(PICInfo);
        }



        /**
         * @brief all RefinerPool must be notified the level levelNumber now exist.
         * not doing so will result in communication to/from that level being impossible
         */
        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magSharedNodesRefiners_.registerLevel(hierarchy, level);
            elecSharedNodesRefiners_.registerLevel(hierarchy, level);

            magPatchGhostsRefiners_.registerLevel(hierarchy, level);
            magGhostsRefiners_.registerLevel(hierarchy, level);
            elecGhostsRefiners_.registerLevel(hierarchy, level);

            rhoGhostsRefiners_.registerLevel(hierarchy, level);
            velGhostsRefiners_.registerLevel(hierarchy, level);
            rhoEGhostsRefiners_.registerLevel(hierarchy, level);
            velEGhostsRefiners_.registerLevel(hierarchy, level);

            patchGhostPartRefiners_.registerLevel(hierarchy, level);

            // root level is not initialized with a schedule using coarser level data
            // so we don't create these schedules if root level
            // TODO this 'if' may not be OK if L0 is regrided
            if (levelNumber != rootLevelNumber)
            {
                // those are for refinement
                magneticInitRefiners_.registerLevel(hierarchy, level);
                electricInitRefiners_.registerLevel(hierarchy, level);
                domainParticlesRefiners_.registerLevel(hierarchy, level);
                lvlGhostPartOldRefiners_.registerLevel(hierarchy, level);
                lvlGhostPartNewRefiners_.registerLevel(hierarchy, level);

                // and these for coarsening
                magnetoSynchronizers_.registerLevel(hierarchy, level);
                electroSynchronizers_.registerLevel(hierarchy, level);
                densitySynchronizers_.registerLevel(hierarchy, level);
                ionBulkVelSynchronizers_.registerLevel(hierarchy, level);
                electronDensitySynchronizers_.registerLevel(hierarchy, level);
                electronBulkVelSynchronizers_.registerLevel(hierarchy, level);
            }
        }



        /**
         * @brief regrid performs the regriding communications for PIC to PIC messengers
         , all quantities that are in initialization refiners need to be regridded
         */
        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                    const int levelNumber,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                    IPhysicalModel& model, double const initDataTime) override
        {
            auto& picModel = dynamic_cast<PICModel&>(model);
            auto level        = hierarchy->getPatchLevel(levelNumber);
            magneticInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            electricInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            domainParticlesRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            patchGhostPartRefiners_.fill(levelNumber, initDataTime);

            // regriding will fill the new level wherever it has points that overlap
            // old level. This will include its level border points.
            // These new level border points will thus take values that where previous
            // domain values. Magnetic flux is thus not necessarily consistent with
            // the Loring et al. method to sync the induction between coarse and fine faces.
            // Specifically, we need all fine faces to have equal magnetic field and also
            // equal to that of the shared coarse face.
            // This means that we now need to fill ghosts and border included
            auto& B = picModel.state.electromag.B;
            auto& E = picModel.state.electromag.E;
            // magSharedNodesRefiners_.fill(B, levelNumber, initDataTime);
            magGhostsRefiners_.fill(B, levelNumber, initDataTime);
            // elecSharedNodesRefiners_.fill(E, levelNumber, initDataTime);
            elecGhostsRefiners_.fill(E, levelNumber, initDataTime);

            fix_magnetic_divergence_(*hierarchy, levelNumber, B);

            // we now call only levelGhostParticlesOld.fill() and not .regrid()
            // regrid() would refine from next coarser in regions of level not overlaping
            // oldLevel, but copy from domain particles of oldLevel where there is an
            // overlap while we do not a priori see why this could be wrong,but this led to
            // occasional failures of the SAMRAI MPI module. See
            // https://github.com/PHAREHUB/PHARE/issues/604 calling .fill() ensures that
            // levelGhostParticlesOld particles are filled exclusively from spliting next
            // coarser domain ones like when a new finest level is created.
            lvlGhostPartOldRefiners_.fill(levelNumber, initDataTime);
            copyLevelGhostOldToPushable_(*level, model);

            // computeIonMoments_(*level, model);
            // levelGhostNew will be refined in next firstStep
        }

        std::string fineModelName() const override { return PICModel::model_name; }

        std::string coarseModelName() const override { return PICModel::model_name; }


        std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<PICMessengerInfo>();
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<PICMessengerInfo>();
        }

        /**
         * @brief initLevel is used to initialize PIC data on the level levelNumer at
         * time initDataTime from PIC coarser data.
         */
        void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                       double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();

            magneticInitRefiners_.fill(levelNumber, initDataTime);
            electricInitRefiners_.fill(levelNumber, initDataTime);

            // no need to call these :
            // magGhostsRefiners_.fill(levelNumber, initDataTime);
            // elecGhostsRefiners_.fill(levelNumber, initDataTime);
            // because the SAMRAI schedules in the 'init' communicators
            // already fill the patch ghost box from the neighbor interior box.
            // so ghost nodes are already filled .

            PHARE_LOG_START("picpicmessengerStrat::initLevel : interior part fill schedule");
            domainParticlesRefiners_.fill(levelNumber, initDataTime);
            PHARE_LOG_STOP("picpicmessengerStrat::initLevel : interior part fill schedule");
            // however we need to call the ghost communicator for patch ghost particles
            // since the interior schedules have a restriction to the interior of the patch.
            PHARE_LOG_START("picpicmessengerStrat::initLevel : patch ghost part fill schedule");
            patchGhostPartRefiners_.fill(levelNumber, initDataTime);
            PHARE_LOG_STOP("picpicmessengerStrat::initLevel : patch ghost part fill schedule");


            lvlGhostPartOldRefiners_.fill(levelNumber, initDataTime);


            // levelGhostParticles will be pushed during the advance phase
            // they need to be identical to levelGhostParticlesOld before advance
            copyLevelGhostOldToPushable_(level, model);

            // computeIonMoments_(level, model);
        }



        /* ------------------------------------------------------------------------
                     methods used for the PICMessenger interface
           ------------------------------------------------------------------------ */



        void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE("PICPICMessengerStrategy::fillElectricGhosts");
            elecSharedNodesRefiners_.fill(E, levelNumber, fillTime);
            elecGhostsRefiners_.fill(E, levelNumber, fillTime);
        }




        /**
         * @brief fillParticleGhostParticles will fill the interior ghost particle array from
         * neighbor patches of the same level. Before doing that, it empties the array for
         * all populations
         */
        void fillGhostParticles(IonsT& ions, PICElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                   double const fillTime) override
        {
            PHARE_LOG_SCOPE("PICPICMessengerStrategy::fillGhostParticles");

            for (auto patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions, electrons);
                for (auto& pop : ions)
                {
                    empty(pop.patchGhostParticles());
                }
                for (auto& pop : electrons)
                {
                    empty(pop.patchGhostParticles());
                }
            }
            patchGhostPartRefiners_.fill(level.getLevelNumber(), fillTime);
        }




        /**
         * @brief fillPopMomentGhosts works on moment ghost nodes
         *
         * patch border node moments are completed by the deposition of patch ghost
         * particles for all populations level border nodes are completed by the deposition
         * of level ghost [old,new] particles for all populations, linear time interpolation
         * is used to get the contribution of old/new particles
         */
        void fillPopMomentGhosts(IonsT& ions, PICElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                    double const afterPushTime) override
        {
            PHARE_LOG_SCOPE("PICPICMessengerStrategy::fillPopMomentGhosts");

            auto alpha = timeInterpCoef_(afterPushTime, level.getLevelNumber());
            if (level.getLevelNumber() > 0 and (alpha < 0 or alpha > 1))
            {
                std::cout << std::setprecision(12) << alpha << "\n";
                throw std::runtime_error("ion moment ghost time interp coef invalid : alpha: "
                                         + std::to_string(alpha) + " beforePushTime "
                                         + std::to_string(afterPushTime) + " on level "
                                         + std::to_string(level.getLevelNumber()));
            }
            for (auto patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions, electrons);
                auto layout      = layoutFromPatch<GridLayoutT>(*patch);

                for (auto& pop : ions)
                {
                    // first thing to do is to project patchGhostParticles moments
                    auto& patchGhosts = pop.patchGhostParticles();
                    auto& density     = pop.density();
                    auto& flux        = pop.flux();

                    interpolate_(makeRange(patchGhosts), density, flux, layout);

                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // then grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        auto& levelGhostOld = pop.levelGhostParticlesOld();
                        interpolate_(makeRange(levelGhostOld), density, flux, layout, 1. - alpha);

                        auto& levelGhostNew = pop.levelGhostParticlesNew();
                        interpolate_(makeRange(levelGhostNew), density, flux, layout, alpha);
                    }
                }
                for (auto& pop : electrons)
                {
                    auto& patchGhosts = pop.patchGhostParticles();
                    auto& density     = pop.density();
                    auto& flux        = pop.flux();

                    interpolate_(makeRange(patchGhosts), density, flux, layout);

                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // then grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        auto& levelGhostOld = pop.levelGhostParticlesOld();
                        interpolate_(makeRange(levelGhostOld), density, flux, layout, 1. - alpha);

                        auto& levelGhostNew = pop.levelGhostParticlesNew();
                        interpolate_(makeRange(levelGhostNew), density, flux, layout, alpha);
                    }
                }
            }
        }


        /* pure (patch and level) ghost nodes are filled by applying a regular ghost
         * schedule i.e. that does not overwrite the border patch node previously well
         * calculated from particles Note : the ghost schedule only fills the total density
         * and bulk velocity and NOT population densities and fluxes. These partial
         * densities and fluxes are thus not available on ANY ghost node.*/
        virtual void fillParticleMomentGhosts(IonsT& ions, PICElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                         double const afterPushTime) override
        {
            rhoGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
            velGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
            rhoEGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
            velEGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
        }

        /**
         * @brief firstStep : in the PICPICMessengerStrategy, the firstStep method is
         * used to get level border ghost particles from the next coarser level. These
         * particles are defined in the future at the time the method is called because the
         * coarser level is ahead in time. These particles are communicated only at first
         * step of a substepping cycle. They will be used with the levelGhostParticlesOld
         * particles to get the moments on level border nodes. The method is does nothing if
         * the level is the root level because the root level cannot get levelGhost from
         * next coarser (it has none).
         */
        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& level,
                       std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                       double const currentTime, double const prevCoarserTime,
                       double const newCoarserTime) override
        {
            PHARE_LOG_SCOPE("PICPICMessengerStrategy::firstStep");

            auto levelNumber = level.getLevelNumber();
            if (newCoarserTime < prevCoarserTime)
                throw std::runtime_error(
                    "Error : prevCoarserTime (" + std::to_string(prevCoarserTime)
                    + ") should be < newCoarserTime (" + std::to_string(prevCoarserTime) + ")");

            // root level has no levelghost particles
            if (levelNumber != 0)
            {
                PHARE_LOG_START("PICPICMessengerStrategy::firstStep.fill");
                lvlGhostPartNewRefiners_.fill(levelNumber, currentTime);
                PHARE_LOG_STOP("PICPICMessengerStrategy::firstStep.fill");

                // during firstStep() coarser level and current level are at the same time
                // so 'time' is also the beforePushCoarseTime_
                beforePushCoarseTime_[levelNumber] = prevCoarserTime;
                afterPushCoarseTime_[levelNumber]  = newCoarserTime;
            }
        }


        /**
         * @brief lastStep is used to perform operations at the last step of a substepping
         * cycle. It is called after the level is advanced. Here for PIC-PIC messages,
         * the method moves levelGhostParticlesNew particles into levelGhostParticlesOld
         * ones. Then levelGhostParticlesNew are emptied since it will be filled again at
         * firstStep of the next substepping cycle. the new CoarseToFineOld content is then
         * copied to levelGhostParticles so that they can be pushed during the next subcycle
         */
        void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) override
        {
            if (level.getLevelNumber() > 0)
            {
                PHARE_LOG_SCOPE("PICPICMessengerStrategy::lastStep");

                auto& picModel = static_cast<PICModel&>(model);
                for (auto& patch : level)
                {
                    auto& ions       = picModel.state.ions;
                    auto& electrons  = picModel.state.pic_electrons;
                    auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                    for (auto& pop : ions)
                    {
                        auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                        auto& levelGhostParticlesNew = pop.levelGhostParticlesNew();
                        auto& levelGhostParticles    = pop.levelGhostParticles();

                        core::swap(levelGhostParticlesNew, levelGhostParticlesOld);
                        core::empty(levelGhostParticlesNew);
                        core::empty(levelGhostParticles);
                        std::copy(std::begin(levelGhostParticlesOld),
                                  std::end(levelGhostParticlesOld),
                                  std::back_inserter(levelGhostParticles));

                        if (level.getLevelNumber() == 0)
                        {
                            if (levelGhostParticlesNew.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesNew detected in root level : "
                                    + std::to_string(levelGhostParticlesNew.size()));
                            if (levelGhostParticles.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticles detected in root level : "
                                    + std::to_string(levelGhostParticles.size()));
                            if (levelGhostParticlesOld.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesOld detected in root level : "
                                    + std::to_string(levelGhostParticlesOld.size()));
                        }
                    }
                    for (auto& pop : electrons)
                    {
                        auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                        auto& levelGhostParticlesNew = pop.levelGhostParticlesNew();
                        auto& levelGhostParticles    = pop.levelGhostParticles();

                        core::swap(levelGhostParticlesNew, levelGhostParticlesOld);
                        core::empty(levelGhostParticlesNew);
                        core::empty(levelGhostParticles);
                        std::copy(std::begin(levelGhostParticlesOld),
                                  std::end(levelGhostParticlesOld),
                                  std::back_inserter(levelGhostParticles));

                        if (level.getLevelNumber() == 0)
                        {
                            if (levelGhostParticlesNew.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesNew detected in root level : "
                                    + std::to_string(levelGhostParticlesNew.size()));
                            if (levelGhostParticles.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticles detected in root level : "
                                    + std::to_string(levelGhostParticles.size()));
                            if (levelGhostParticlesOld.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesOld detected in root level : "
                                    + std::to_string(levelGhostParticlesOld.size()));
                        }
                    }
                }
            }
        }



        /**
         * @brief prepareStep is the concrete implementation of the
         * PICMessengerStrategy::prepareStep method For PIC-PIC communications.
         * This method copies the density and bulk velocity, defined at t=n.
         * Since prepareStep() is called just before advancing the level, this operation
         * actually saves the t=n versions of Ni, Vi into the messenger. When the time comes
         * that the next finer level needs to time interpolate the electromagnetic field
         * at its ghost nodes, this level will be able to interpolate at required time
         * because the t=n Vi,Ni fields of previous next coarser step will be in the
         * messenger.
         */
        void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                         double currentTime) override
        {
            PHARE_LOG_SCOPE("PICPICMessengerStrategy::prepareStep");

            auto& picModel = static_cast<PICModel&>(model);
            auto& picState = picModel.state;
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, picState.electromag, picState.ions, 
                     picState.pic_electrons, NiOldUser_, ViOld_, NeOldUser_, VeOld_);

                resourcesManager_->setTime(NiOldUser_, *patch, currentTime);
                resourcesManager_->setTime(ViOld_, *patch, currentTime);
                resourcesManager_->setTime(NeOldUser_, *patch, currentTime);
                resourcesManager_->setTime(VeOld_, *patch, currentTime);

                auto& Vi = picState.ions.velocity();
                auto& Ni = picState.ions.density();
                auto& Ve = picState.pic_electrons.velocity();
                auto& Ne = picState.pic_electrons.density();

                printf("prepareStep : copying ViOld\n");
                ViOld_.copyData(Vi);
                printf("prepareStep : copying NiOld\n");
                NiOldUser_.copyData(Ni);
                printf("prepareStep : copying VeOld\n");
                VeOld_.copyData(Ve);
                printf("prepareStep : copying NeOld\n");
                NeOldUser_.copyData(Ne);
            }
            printf("End prepareStep\n");
        }




        void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                            double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();
            assert(levelNumber == 0);

            auto& picModel = static_cast<PICModel&>(model);

            elecSharedNodesRefiners_.fill(picModel.state.electromag.E, levelNumber,
                                          initDataTime);

            elecGhostsRefiners_.fill(picModel.state.electromag.E, levelNumber, initDataTime);
            patchGhostPartRefiners_.fill(levelNumber, initDataTime);

            // at some point in the future levelGhostParticles could be filled with injected
            // particles depending on the domain boundary condition.
            //
            // do we need moment ghosts filled here?
            // a priori no because those are at this time only needed for coarsening, which
            // will not happen before the first advance
        }



        void synchronize(SAMRAI::hier::PatchLevel& level) override
        {
            PHARE_LOG_SCOPE("PICPICMessengerStrategy::synchronize");

            auto levelNumber = level.getLevelNumber();
            std::cout << "synchronizing level " << levelNumber << "\n";

            // call coarsening schedules...
            magnetoSynchronizers_.sync(levelNumber);
            electroSynchronizers_.sync(levelNumber);
            densitySynchronizers_.sync(levelNumber);
            ionBulkVelSynchronizers_.sync(levelNumber);
            electronDensitySynchronizers_.sync(levelNumber);
            electronBulkVelSynchronizers_.sync(levelNumber);
        }

        // after coarsening, domain nodes have been updated and therefore patch ghost nodes
        // will probably stop having the exact same value as their overlapped neighbor
        // domain node we thus fill ghost nodes. note that we first fill shared border nodes
        // with the sharedNode refiners so that these shared nodes agree on their value at
        // MPI process boundaries. then regular refiner fill are called, which fill only
        // pure ghost nodes. note also that moments are not filled on border nodes since
        // already OK from particle deposition
        void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                             double const time) override
        {
            auto levelNumber  = level.getLevelNumber();
            auto& picModel = static_cast<PICModel&>(model);

            std::cout << "postSynchronize level " << levelNumber << "\n";

            magSharedNodesRefiners_.fill(picModel.state.electromag.B, levelNumber, time);
            elecSharedNodesRefiners_.fill(picModel.state.electromag.E, levelNumber, time);

            // we fill magnetic field ghosts only on patch ghost nodes and not on level
            // ghosts the reason is that 1/ filling ghosts is necessary to prevent mismatch
            // between ghost and overlaped neighboring patch domain nodes resulting from
            // former coarsening which does not occur for level ghosts and 2/ overwriting
            // level border with next coarser model B would invalidate divB on the first
            // fine domain cell since its border face only received a fraction of the
            // induction that has occured on the shared coarse face.
            magPatchGhostsRefiners_.fill(picModel.state.electromag.B, levelNumber, time);
            elecGhostsRefiners_.fill(picModel.state.electromag.E, levelNumber, time);
            rhoGhostsRefiners_.fill(levelNumber, time);
            velGhostsRefiners_.fill(picModel.state.ions.velocity(), levelNumber, time);
            rhoEGhostsRefiners_.fill(levelNumber, time);
            velEGhostsRefiners_.fill(picModel.state.pic_electrons.velocity(), levelNumber, time);
        }

    private:
        void registerGhostComms_(std::unique_ptr<PICMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& vecFieldNames) {
                std::vector<std::string> keys;
                std::transform(std::begin(vecFieldNames), std::end(vecFieldNames),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };
            magSharedNodesRefiners_.addStaticRefiners(info->ghostMagnetic, BfieldNodeRefineOp_,
                                                      makeKeys(info->ghostMagnetic));

            magGhostsRefiners_.addStaticRefiners(info->ghostMagnetic, BfieldRefineOp_,
                                                 makeKeys(info->ghostMagnetic));

            magPatchGhostsRefiners_.addStaticRefiner(info->modelMagnetic, BfieldRefineOp_,
                                                     info->modelMagnetic.vecName);

            elecSharedNodesRefiners_.addStaticRefiners(info->ghostElectric, EfieldNodeRefineOp_,
                                                       makeKeys(info->ghostElectric));

            elecGhostsRefiners_.addStaticRefiners(info->ghostElectric, EfieldRefineOp_,
                                                  makeKeys(info->ghostElectric));

            rhoGhostsRefiners_.addTimeRefiner(info->modelIonDensity, info->modelIonDensity,
                                              NiOldUser_.name, fieldRefineOp_, fieldTimeOp_,
                                              info->modelIonDensity);

            rhoEGhostsRefiners_.addTimeRefiner(info->modelElectronDensity, info->modelElectronDensity,
                                               NeOldUser_.name, fieldRefineOp_, fieldTimeOp_,
                                               info->modelElectronDensity);

            velGhostsRefiners_.addTimeRefiners(info->ghostBulkVelocity, info->modelIonBulkVelocity,
                                               core::VecFieldNames{ViOld_}, fieldRefineOp_,
                                               fieldTimeOp_);

            velEGhostsRefiners_.addTimeRefiners(info->ghostBulkVelocity, info->modelElectronBulkVelocity,
                                                core::VecFieldNames{VeOld_}, fieldRefineOp_,
                                                fieldTimeOp_);
        }




        void registerInitComms(std::unique_ptr<PICMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& descriptor) {
                std::vector<std::string> keys;
                std::transform(std::begin(descriptor), std::end(descriptor),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            magneticInitRefiners_.addStaticRefiners(info->initMagnetic, BfieldRefineOp_,
                                                    makeKeys(info->initMagnetic));

            electricInitRefiners_.addStaticRefiners(info->initElectric, EfieldRefineOp_,
                                                    makeKeys(info->initElectric));


            domainParticlesRefiners_.addStaticRefiners(
                info->interiorParticles, interiorParticleRefineOp_, info->interiorParticles);


            lvlGhostPartOldRefiners_.addStaticRefiners(info->levelGhostParticlesOld,
                                                       levelGhostParticlesOldOp_,
                                                       info->levelGhostParticlesOld);


            lvlGhostPartNewRefiners_.addStaticRefiners(info->levelGhostParticlesNew,
                                                       levelGhostParticlesNewOp_,
                                                       info->levelGhostParticlesNew);


            patchGhostPartRefiners_.addStaticRefiners(info->patchGhostParticles, nullptr,
                                                      info->patchGhostParticles);
        }




        void registerSyncComms(std::unique_ptr<PICMessengerInfo> const& info)
        {
            magnetoSynchronizers_.add(info->modelMagnetic, magneticCoarseningOp_,
                                      info->modelMagnetic.vecName);

            electroSynchronizers_.add(info->modelElectric, fieldCoarseningOp_,
                                      info->modelElectric.vecName);

            ionBulkVelSynchronizers_.add(info->modelIonBulkVelocity, fieldCoarseningOp_,
                                         info->modelIonBulkVelocity.vecName);

            electronBulkVelSynchronizers_.add(info->modelElectronBulkVelocity, fieldCoarseningOp_,
                                              info->modelElectronBulkVelocity.vecName);

            densitySynchronizers_.add(info->modelIonDensity, fieldCoarseningOp_,
                                      info->modelIonDensity);
            
            electronDensitySynchronizers_.add(info->modelElectronDensity, fieldCoarseningOp_,
                                        info->modelElectronDensity);
        }




        void copyLevelGhostOldToPushable_(SAMRAI::hier::PatchLevel& level, IPhysicalModel& model)
        {
            auto& picModel = static_cast<PICModel&>(model);
            for (auto& patch : level)
            {
                auto& ions       = picModel.state.ions;
                auto& electrons  = picModel.state.pic_electrons;
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions, electrons);
                for (auto& pop : ions)
                {
                    auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                    auto& levelGhostParticles    = pop.levelGhostParticles();

                    core::empty(levelGhostParticles);
                    std::copy(std::begin(levelGhostParticlesOld), std::end(levelGhostParticlesOld),
                              std::back_inserter(levelGhostParticles));
                }
                for (auto& pop : electrons)
                {
                    auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                    auto& levelGhostParticles    = pop.levelGhostParticles();

                    core::empty(levelGhostParticles);
                    std::copy(std::begin(levelGhostParticlesOld), std::end(levelGhostParticlesOld),
                              std::back_inserter(levelGhostParticles));
                }
            }
        }




        double timeInterpCoef_(double const afterPushTime, std::size_t levelNumber)
        {
            return (afterPushTime - beforePushCoarseTime_[levelNumber])
                   / (afterPushCoarseTime_[levelNumber] - beforePushCoarseTime_[levelNumber]);
        }


        void debug_print(VecFieldT const& B, GridLayoutT const& layout, int loc, int ix, int iy,
                         std::string const& aftbef)
        {
            auto& Bx       = B(core::Component::X);
            auto& By       = B(core::Component::Y);
            auto& Bz       = B(core::Component::Z);
            auto const& dx = layout.meshSize()[0];
            auto const& dy = layout.meshSize()[1];

            if (loc == 3) // w hi, y hi
            {
                std::cout << aftbef << "\n";
                std::cout << "cell 4 : "
                          << (Bx(ix, iy) - Bx(ix - 1, iy)) / dx
                                 + (By(ix - 1, iy + 1) - By(ix - 1, iy)) / dy
                          << "\n";

                std::cout << "cell 2 : "
                          << (Bx(ix + 1, iy - 1) - Bx(ix, iy - 1)) / dx
                                 + (By(ix, iy) - By(ix, iy - 1)) / dy
                          << "\n";

                std::cout << "cell 1 : "
                          << (Bx(ix + 1, iy) - Bx(ix, iy)) / dx + (By(ix, iy + 1) - By(ix, iy)) / dy
                          << "\n";

                std::cout << "cell 3 : "
                          << (Bx(ix, iy - 1) - Bx(ix - 1, iy - 1)) / dx
                                 + (By(ix - 1, iy) - By(ix - 1, iy - 1)) / dy
                          << "\n";
            }
        }

        /*
         *
         * */
        void fix_magnetic_divergence_(SAMRAI::hier::PatchHierarchy const& hierarchy,
                                      int levelNumber, VecFieldT& B)
        {
            auto lvlBoundary = SAMRAI::hier::CoarseFineBoundary(
                hierarchy, levelNumber,
                SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension}, 0});

            for (auto& patch : *hierarchy.getPatchLevel(levelNumber))
            {
                if constexpr (dimension == 2)
                {
                    auto _         = resourcesManager_->setOnPatch(*patch, B);
                    auto layout    = layoutFromPatch<GridLayoutT>(*patch);
                    auto& Bx       = B(core::Component::X);
                    auto& By       = B(core::Component::Y);
                    auto& Bz       = B(core::Component::Z);
                    auto const& dx = layout.meshSize()[0];
                    auto const& dy = layout.meshSize()[1];

                    auto boundaries = lvlBoundary.getEdgeBoundaries(patch->getGlobalId());
                    for (auto& boundary : boundaries)
                    {
                        int loc         = boundary.getLocationIndex();
                        auto const& box = boundary.getBox();

                        if (loc == 0) // x_lo
                        {
                            // we're on the left side border outside the domain
                            // we need to get to the first domain cell
                            auto ixAMR = box.lower()[0] + 1;

                            // this is a X edge we need to fix By
                            for (int iyAMR = box.lower()[1]; iyAMR <= box.upper()[1]; ++iyAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(iyAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    By(ix, iy)
                                        = By(ix, iy + 1) + dy / dx * (Bx(ix + 1, iy) - Bx(ix, iy));
                                }
                            }
                        }
                        else if (loc == 1) // x_hi
                        {
                            // we're on the right side border outside the domain
                            // we need to get to the last domain cell
                            auto ixAMR = box.upper()[0] - 1;

                            // this is a X edge we need to fix By
                            for (int iyAMR = box.lower()[1]; iyAMR <= box.upper()[1]; ++iyAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(iyAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    By(ix, iy)
                                        = By(ix, iy + 1) + dy / dx * (Bx(ix + 1, iy) - Bx(ix, iy));
                                }
                            }
                        }
                        else if (loc == 2) // y_lo
                        {
                            // we're on the bottom edge, we need the first domain cell
                            auto iyAMR = box.lower()[1] + 1;

                            // this is a Y edge we need to fix Bx
                            for (int ixAMR = box.lower()[0]; ixAMR <= box.upper()[0]; ++ixAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(ixAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    Bx(ix, iy)
                                        = Bx(ix + 1, iy) + dx / dy * (By(ix, iy + 1) - By(ix, iy));
                                }
                            }
                        }
                        else if (loc == 3) // y_hi
                        {
                            // we're on the top edge, we need the last domain cell
                            auto iyAMR = box.upper()[1] - 1;

                            // this is a Y edge we need to fix Bx
                            for (int ixAMR = box.lower()[0]; ixAMR <= box.upper()[0]; ++ixAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(ixAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    Bx(ix, iy)
                                        = Bx(ix + 1, iy) + dx / dy * (By(ix, iy + 1) - By(ix, iy));
                                }
                            }
                        }
                    } // end boundary loop

                    // above we have treated boundaries as 1D lines
                    // this is not treating corners well and we need to deal with them separatly
                    //
                    //__________________
                    //      |    |     |
                    //      | 4  BX 1  |
                    //      |    |     |
                    //      |_By_|_BY___
                    //      |    |     |
                    //      | 3 Bx  2  |
                    //      |    |     |
                    //      |____|_____
                    //                 |
                    //                 |
                    //                 |
                    //
                    // above are the 4 fine top right corner cells and
                    // above code have changed BX and BY but not Bx and By
                    // faces. But these cells are not independant.
                    // To fix divB in these 4 cells we will re-assigne the four
                    // fine faces.
                    // The idea is to give them the same values they would have had
                    // if refined. By and BY would have the same value, equal to the
                    // average of the coarse By on top and bottom faces of the coarse cell
                    // we do not have access to the coarse cell here but we know that because
                    // of the previous coarsening, these coarse faces have the same flux
                    // as the average of the 2 shared fine faces
                    // So By and BY will be 1/4 * sum of top and bottom fine By
                    // Similarly, BX and Bx will take 1/4 the four fine Bx faces share with the
                    // coarse ones
                    //
                    auto corners = lvlBoundary.getNodeBoundaries(patch->getGlobalId());
                    for (auto& corner : corners)
                    {
                        int loc = corner.getLocationIndex();

                        // the box we get should consist of just 1 cell
                        // i.e. with lower==upper so it should not matter which
                        // one we take in the following.
                        auto const& box = corner.getBox();

                        //* x_lo, y_lo: 0
                        //* x_hi, y_lo: 1
                        //* x_lo, y_hi: 2
                        // * x_hi, y_hi: 3

                        if (loc == 3) // x_hi, y_hi
                        {
                            // we're on the top right corner
                            // and we want the domain cell
                            // that is labeled cell 1 in above drawing
                            // we fix BX and BY
                            auto ixAMR    = box.lower()[0] - 1;
                            auto iyAMR    = box.lower()[1] - 1;
                            auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                            auto ix       = localIdx[0];
                            auto iy       = localIdx[1];

                            // maybe we should keep these for some time
                            // as comments in case they are useful again
                            PHARE_DEBUG_DO(const std::string before = "BEFORE";
                                           debug_print(B, layout, loc, ix, iy, before);)
                            Bx(ix, iy - 1)
                                = Bx(ix + 1, iy - 1) + dx / dy * (By(ix, iy) - By(ix, iy - 1));

                            By(ix - 1, iy)
                                = By(ix - 1, iy + 1) + dy / dx * (Bx(ix, iy) - Bx(ix - 1, iy));

                            PHARE_DEBUG_DO(const std::string after = "AFTER";
                                           debug_print(B, layout, loc, ix, iy, after);)
                        }
                    } // end corner loops
                }     // end if 2D
            }         // end patch loop
        }


        VecFieldT ViOld_{stratName + "_VBulkOld", core::HybridQuantity::Vector::V};
        FieldT* NiOld_{nullptr};
        FieldUser<FieldT> NiOldUser_{stratName + "_NiOld", NiOld_,
                                     core::HybridQuantity::Scalar::rho};
        VecFieldT VeOld_{stratName + "_VeOld", core::HybridQuantity::Vector::Ve};
        FieldT* NeOld_{nullptr};
        FieldUser<FieldT> NeOldUser_{stratName + "_NeOld", NeOld_,
                                     core::HybridQuantity::Scalar::rhoE};


        //! ResourceManager shared with other objects (like the PICModel)
        std::shared_ptr<ResourcesManagerT> resourcesManager_;


        int const firstLevel_;
        std::unordered_map<std::size_t, double> beforePushCoarseTime_;
        std::unordered_map<std::size_t, double> afterPushCoarseTime_;

        core::Interpolator<dimension, interpOrder> interpolate_;

        using rm_t                    = ResourcesManagerT;
        using RefineOperator          = SAMRAI::hier::RefineOperator;
        using TimeInterpolateOperator = SAMRAI::hier::TimeInterpolateOperator;

        // these refiners are used to initialize electromagnetic fields when creating
        // a new level (initLevel) or regridding (regrid)
        using InitRefinerPool           = RefinerPool<rm_t, RefinerType::InitField>;
        using SharedNodeRefinerPool     = RefinerPool<rm_t, RefinerType::SharedBorder>;
        using GhostRefinerPool          = RefinerPool<rm_t, RefinerType::GhostField>;
        using PatchGhostRefinerPool     = RefinerPool<rm_t, RefinerType::PatchGhostField>;
        using InitDomPartRefinerPool    = RefinerPool<rm_t, RefinerType::InitInteriorPart>;
        using PatchGhostPartRefinerPool = RefinerPool<rm_t, RefinerType::InteriorGhostParticles>;

        InitRefinerPool magneticInitRefiners_{resourcesManager_};
        InitRefinerPool electricInitRefiners_{resourcesManager_};

        //! store communicators for magnetic fields that need ghosts to be filled
        SharedNodeRefinerPool magSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool magGhostsRefiners_{resourcesManager_};
        PatchGhostRefinerPool magPatchGhostsRefiners_{resourcesManager_};


        //! store refiners for electric fields that need ghosts to be filled
        SharedNodeRefinerPool elecSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};


        // moment ghosts
        // these do not need sharedNode refiners. The reason is that
        // the border node is already complete by the deposit of ghost particles
        // these refiners are used to fill ghost nodes, and therefore, owing to
        // the GhostField tag, will only assign pur ghost nodes. Border nodes will
        // be overwritten only on level borders, which does not seem to be an issue.
        GhostRefinerPool rhoGhostsRefiners_{resourcesManager_};
        GhostRefinerPool velGhostsRefiners_{resourcesManager_};
        GhostRefinerPool rhoEGhostsRefiners_{resourcesManager_};
        GhostRefinerPool velEGhostsRefiners_{resourcesManager_};

        // pool of refiners for interior particles of each population
        // and the associated refinement operator
        InitDomPartRefinerPool domainParticlesRefiners_{resourcesManager_};

        using RefOp_ptr = std::shared_ptr<RefineOperator>;

        RefOp_ptr interiorParticleRefineOp_{std::make_shared<InteriorParticleRefineOp>()};

        //! store communicators for coarse to fine particles old
        // pools of refiners to fill level ghost particles, old and new ones
        // and their associated refinement operator
        static auto constexpr LGRefT = RefinerType::LevelBorderParticles;
        RefinerPool<rm_t, LGRefT> lvlGhostPartOldRefiners_{resourcesManager_};
        RefinerPool<rm_t, LGRefT> lvlGhostPartNewRefiners_{resourcesManager_};
        RefOp_ptr levelGhostParticlesOldOp_{std::make_shared<CoarseToFineRefineOpOld>()};
        RefOp_ptr levelGhostParticlesNewOp_{std::make_shared<CoarseToFineRefineOpNew>()};


        // this contains refiners for each population to exchange patch ghost particles
        PatchGhostPartRefinerPool patchGhostPartRefiners_{resourcesManager_};

        SynchronizerPool<rm_t> densitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> ionBulkVelSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> electronDensitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> electronBulkVelSynchronizers_{resourcesManager_};     
        SynchronizerPool<rm_t> electroSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> magnetoSynchronizers_{resourcesManager_};


        RefOp_ptr fieldRefineOp_{std::make_shared<DefaultFieldRefineOp>()};
        // see field_variable_fill_pattern.hpp for explanation about this "node_only" flag
        // note that refinement operator, via the boolean argument, serve as a relay for the
        // the refinealgorithm to get the correct variablefillpattern
        RefOp_ptr BfieldNodeRefineOp_{std::make_shared<MagneticFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr BfieldRefineOp_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr EfieldNodeRefineOp_{std::make_shared<ElectricFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr EfieldRefineOp_{std::make_shared<ElectricFieldRefineOp>()};

        std::shared_ptr<TimeInterpolateOperator> fieldTimeOp_{std::make_shared<FieldTimeInterp>()};

        using CoarsenOperator_ptr = std::shared_ptr<SAMRAI::hier::CoarsenOperator>;
        CoarsenOperator_ptr fieldCoarseningOp_{std::make_shared<DefaultCoarsenOp>()};
        CoarsenOperator_ptr magneticCoarseningOp_{std::make_shared<MagneticCoarsenOp>()};
    };

    template<typename PICModel, typename RefinementParams>
    const std::string PICPICMessengerStrategy<PICModel, RefinementParams>::stratName
        = "PICModel-PICModel";

} // namespace amr

} // namespace PHARE

#endif
