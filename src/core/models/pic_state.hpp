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
 * It holds an Electromag, Ions and PICElectrons objects manipulated by a  PIC concrete type of
 * ISolver
 */
template<typename Electromag, typename Ions, typename PICElectrons>
class PICState : public IPhysicalState
{

    using VecField = typename Electromag::vecfield_type;
    
public:
    static constexpr auto dimension = Ions::dimension;

    PICState(PHARE::initializer::PHAREDict const& dict)
        : electromag{dict["electromag"]}
        , ions{dict["ions"]}
        , pic_electrons{dict["pic_electrons"]}
        , J{"J", HybridQuantity::Vector::J}
    {
    }

    Electromag electromag;
    PICElectrons pic_electrons;
    Ions ions;
    VecField J;


    NO_DISCARD std::string to_str()
    {
        std::stringstream ss;
        ss << "PIC State\n";
        ss << "------------------------------------\n";
        ss << core::to_str(ions);
        ss << core::to_str(pic_electrons);
        return ss.str();
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const
        {
            return electromag.isUsable() && ions.isUsable() && pic_electrons.isUsable() && J.isUsable();
        }



        NO_DISCARD bool isSettable() const
        {
            return electromag.isSettable() && ions.isSettable() && pic_electrons.isSettable() && J.isSettable();
        }


        NO_DISCARD auto getCompileTimeResourcesUserList() const
        {
            return std::forward_as_tuple(electromag, ions, pic_electrons, J);
        }

        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(electromag, ions, pic_electrons, J);
        }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

};// end PICState 

} // namespace PHARE::core


#endif // PHARE_PIC_STATE_HPP
