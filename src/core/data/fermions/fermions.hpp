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
        using vecfield_type             = typename Ions::vecfield_type;
        using particle_array_type       = typename Ions::particle_array_type;
        using ions_type                 = Ions;
        using pic_electrons_type         = PICElectrons;
        static constexpr auto dimension = Ions::dimension;

    
/*
        NO_DISCARD std::string to_str() // CHECK
        {
            std::stringstream ss;
            ss << "Fermions\n";
            ss << "------------------------------------\n";
            ss << "number of particle types  : " << 2 << "\n";
            for (auto& particletype : (ions, electrons))
                ss << core::to_str(particletype);
            return ss.str();
        }
*/
        // ResourcesUser interface??
        
    };
} // namespace PHARE::core

#endif
