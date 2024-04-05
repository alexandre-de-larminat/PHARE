#ifndef PHARE_FERMIONS_HPP
#define PHARE_FERMIONS_HPP

#include "core/data/ions/ions.hpp"
#include "core/data/fermions/pic_electrons.hpp"
#include "core/data/electrons/electrons.hpp"

namespace PHARE::core
{
    template<typename Ions, typename PICElectrons>
    class Fermions
    {
    public:
        using vecfield_type       = typename Ions::vecfield_type;
        using particle_array_type = typename Ions::particle_array_type;
        using ions_type           = Ions;
        using fluid_electrons_type      = Electrons<ions_type>;
        static constexpr auto dimension = Ions::dimension;

        explicit Fermions(PHARE::initializer::PHAREDict const& dict)
        : ions{dict["ions"]}
        , electrons{dict["pic_electrons"]}
        {
        }

        Ions ions;
        PICElectrons electrons;

        // ResourcesUser interface??
        
    };
} // namespace PHARE::core

#endif
