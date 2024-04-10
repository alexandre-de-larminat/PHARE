#ifndef PHARE_PIC_STATE_HPP
#define PHARE_PIC_STATE_HPP

#include "core/hybrid/hybrid_quantities.hpp"
#include "core/models/physical_state.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"
#include "core/utilities/algorithm.hpp"

#include <cstddef>
#include <sstream>
#include <string>
#include <utility>
#include <functional>

namespace PHARE::core
{

/**
 * @brief The PICState class is a concrete implementation of a IPhysicalState.
 * It holds an Electromag and an Fermion object manipulated by a  PIC concrete type of
 * ISolver
 */
template<typename Electromag, typename Fermions>
class PICState : public IPhysicalState
{

    using VecField = typename Electromag::vecfield_type;
    
public:
    static constexpr auto dimension = Fermions::dimension;

    PICState(PHARE::initializer::PHAREDict const& dict)
        : electromag{dict["electromag"]}
        , fermions{dict}
        , J{"J", HybridQuantity::Vector::J}
    {
    }

    Electromag electromag;
    Fermions fermions;
    VecField J;

/*
    NO_DISCARD std::string to_str() // CHECK
    {
        std::stringstream ss;
        ss << "PIC State\n";
        ss << "------------------------------------\n";
        for (auto& particletype : fermions)
            ss << core::to_str(particletype);
        return ss.str();
    }
*/

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const
        {
            return electromag.isUsable() && fermions.isUsable() && J.isUsable();
        }



        NO_DISCARD bool isSettable() const
        {
            return electromag.isSettable() && fermions.isSettable() && J.isSettable();
        }


        NO_DISCARD auto getCompileTimeResourcesUserList() const
        {
            return std::forward_as_tuple(electromag, fermions);
        }

        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(electromag, fermions);
        }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

};// end PICState 

} // namespace PHARE::core


#endif // PHARE_PIC_STATE_HPP
