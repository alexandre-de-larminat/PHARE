#ifndef PHARE_PIC_MESSENGER_STRATEGY_HPP
#define PHARE_PIC_MESSENGER_STRATEGY_HPP

#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/messenger_strategy.hpp"
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>


#include <utility>


namespace PHARE
{
namespace amr
{
    template<typename PICModel>
    class PICMessengerStrategy : public MessengerStrategy<PICModel>
    {
        using IonsT          = decltype(std::declval<PICModel>().state.ions);
        using PICElectronsT    = decltype(std::declval<PICModel>().state.pic_electrons);

    public:


        virtual void fillGhostParticles(IonsT& ions, PICElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                           double const fillTime)
            = 0;


        virtual void fillPopMomentGhosts(IonsT& ions, PICElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                            double const afterPushTime)
            = 0;


        virtual void fillParticleMomentGhosts(IonsT& ions, PICElectronsT& electrons, SAMRAI::hier::PatchLevel& level,
                                         double const afterPushTime)
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
