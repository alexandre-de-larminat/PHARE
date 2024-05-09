#ifndef PHARE_PIC_MESSENGER_STRATEGY_HPP
#define PHARE_PIC_MESSENGER_STRATEGY_HPP

#include "amr/messengers/messenger_info.hpp"

#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>


#include <utility>


namespace PHARE
{
namespace amr
{
    template<typename PICModel>
    class PICMessengerStrategy
    {
        using IonsT          = decltype(std::declval<PICModel>().state.ions);
        using ElectronsT     = decltype(std::declval<PICModel>().state.pic_electrons);
        using VecFieldT      = decltype(std::declval<PICModel>().state.electromag.E);
        using IPhysicalModel = typename PICModel::Interface;

    public:
        /**
         * @brief allocate PICMessengerStrategy internal resources on a given patch for a given
         * allocation time. The method is virtual and is implemented by a concrete
         * PICMessengerStrategy
         */
        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const = 0;


        /**
         * @brief setup prepares the PICMessengerStrategy to communicate specific data between
         * two levels. The method is called by a PICMessenger::allocate() and its concrete
         * implementation is found in concrete strategies
         */
        virtual void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo)
            = 0;



        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() = 0;
        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner()   = 0;


        virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   int const levelNumber)
            = 0;

        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            const int levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            IPhysicalModel& model, double const initDataTime)
            = 0;

        // domain filling methods
        // used during level creation
        virtual void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               double const initDataTime)
            = 0;

        // ghost filling

        virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime)
            = 0;


        virtual void fillCurrentGhosts(VecFieldT& J, int const levelNumber, double const fillTime)
            = 0;


        virtual void fillIonGhostParticles(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                           double const fillTime)
            = 0;


        virtual void fillIonPopMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                            double const afterPushTime)
            = 0;


        virtual void fillIonMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                         double const afterPushTime)
            = 0;


        virtual void fillElectronGhostParticles(ElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                           double const fillTime)
            = 0;


        virtual void fillElectronPopMomentGhosts(ElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                            double const afterPushTime)
            = 0;


        virtual void fillElectronMomentGhosts(ElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                         double const afterPushTime)
            = 0;

        virtual std::string fineModelName() const = 0;

        virtual std::string coarseModelName() const = 0;

        virtual void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               double const currentTime, double const prevCoarserTime,
                               double const newCoarserTime)
            = 0;

        virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) = 0;



        virtual void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                 double currentTime)
            = 0;


        virtual void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                    double const initDataTime)
            = 0;


        virtual void synchronize(SAMRAI::hier::PatchLevel& level) = 0;

        virtual void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                     double const time)
            = 0;


        std::string name() const { return stratname_; }

        virtual ~PICMessengerStrategy() = default;


    protected:
        explicit PICMessengerStrategy(std::string stratName)
            : stratname_{stratName}
        {
        }

        std::string stratname_;
    };
} // namespace amr
} // namespace PHARE

#endif
