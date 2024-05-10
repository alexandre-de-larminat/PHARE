#ifndef PHARE_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_MESSENGER_STRATEGY_HPP

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
    template<typename HybridModel>
    class HybridMessengerStrategy : public MessengerStrategy<HybridModel>
    {
        using VecFieldT      = decltype(std::declval<HybridModel>().state.electromag.E);
        using IonsT          = decltype(std::declval<HybridModel>().state.ions);

    public:

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



        std::string name() const { return stratname_; }

        virtual ~HybridMessengerStrategy() = default;


    protected:
        explicit HybridMessengerStrategy(std::string stratName)
            : stratname_{stratName}
        {
        }

        std::string stratname_;
    };
} // namespace amr
} // namespace PHARE

#endif
