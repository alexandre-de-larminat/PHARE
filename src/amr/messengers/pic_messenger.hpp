
#ifndef PHARE_PIC_MESSENGER_HPP
#define PHARE_PIC_MESSENGER_HPP

#include <memory>
#include <string>
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/RefineOperator.h>

#include "core/hybrid/hybrid_quantities.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/pic_messenger_info.hpp"

namespace PHARE
{
namespace amr
{
    template<typename PICModel>
    class PICMessenger : public IMessenger<typename PICModel::Interface>
    {
    public:
        using IPhysicalModel = typename PICModel::Interface;
        PICMessenger(std::shared_ptr<typename PICModel::resources_manager_type> resourcesManager,
                     int const firstLevel)
            : resourcesManager_{std::move(resourcesManager)}
            , firstLevel_{firstLevel}
        {
        }

        void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            printf("PICMessenger::registerQuantities\n");
            std::unique_ptr<HybridMessengerInfo> xInfo{
                dynamic_cast<HybridMessengerInfo*>(fromCoarserInfo.release())};
        }



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                           int const /*levelNumber*/) override
        {
        }


        static const std::string stratName;

        std::string fineModelName() const override { return PICModel::model_name; }

        std::string coarseModelName() const override { return PICModel::model_name; }

        void allocate(SAMRAI::hier::Patch& /*patch*/, double const /*allocateTime*/) const override
        {
        }

        void initLevel(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       double const /*initDataTime*/) override
        {
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<HybridMessengerInfo>();
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<HybridMessengerInfo>();
        }




        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                    const int /*levelNumber*/,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*oldLevel*/,
                    IPhysicalModel& /*model*/, double const /*initDataTime*/) override
        {
        }


        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& /*hierarchy*/,
                       double const /*currentTime*/, double const /*prevCoarserTIme*/,
                       double const /*newCoarserTime*/) final
        {
        }


        void lastStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/) final {}


        void prepareStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                         double /*currentTime*/) final
        {
        }

        void fillRootGhosts(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                            double const /*initDataTime*/) final
        {
        }



        void synchronize(SAMRAI::hier::PatchLevel& /*level*/) final
        {
            // call coarsning schedules...
        }

        void postSynchronize(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                             double const /*time*/) override
        {
        }


        std::string name() override { return stratName; }

        virtual ~PICMessenger() = default;


    private:
        std::shared_ptr<typename PICModel::resources_manager_type> resourcesManager_;
        int const firstLevel_;
    };


    template<typename PICModel>
    const std::string PICMessenger<PICModel>::stratName = "PICModel-PICModel";
} // namespace amr
} // namespace PHARE
#endif
