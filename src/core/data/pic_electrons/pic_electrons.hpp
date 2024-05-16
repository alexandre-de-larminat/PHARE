#ifndef PHARE_PIC_ELECTRONS_HPP
#define PHARE_PIC_ELECTRONS_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>


#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/data/pic_electrons/electron_population.hpp"
#include "core/data/pic_electrons/plasma_components.hpp"

namespace PHARE::core
{
    template<typename ElectronPopulation, typename GridLayout>
    class PICElectrons : public PlasmaComponents<ElectronPopulation, GridLayout>
    {
    public:
        using field_type          = typename ElectronPopulation::field_type;
        using vecfield_type       = typename ElectronPopulation::vecfield_type;
        using Float               = typename field_type::type;
        using particle_array_type = typename ElectronPopulation::particle_array_type;
        using gridlayout_type     = GridLayout;

        using Super = PlasmaComponents<ElectronPopulation, GridLayout>;
        using Super::bulkVelocity_;
        using Super::populations_;

        explicit PICElectrons(PHARE::initializer::PHAREDict const& dict)
            : Super{dict}
        {
            if (!Super::sameMasses_)
                throw std::runtime_error("Error - all electron populations must have the same mass");
        }



        NO_DISCARD std::string densityName() const override { return "rhoE"; }

        NO_DISCARD std::string bulkVelocityName() const override { return "bulkVelE"; }

        NO_DISCARD field_type* rhoPtr() override { return rho_; }


        NO_DISCARD bool isUsable() const override
        {
            bool usable
                = rho_ != nullptr && bulkVelocity_.isUsable();

            for (auto const& pop : populations_)
            {
                usable = usable && pop.isUsable();
            }
            return usable;
        }



        NO_DISCARD bool isSettable() const override
        {
            bool settable
                = rho_ == nullptr && bulkVelocity_.isSettable();

            for (auto const& pop : populations_)
            {
                settable = settable && pop.isSettable();
            }
            return settable;
        }


        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------


        void setBuffer(std::string const& bufferName, field_type* field)
        {
            if (bufferName == densityName())
                rho_ = field;
            else
                throw std::runtime_error("Error - invalid density buffer name : " + bufferName);
        }


        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(bulkVelocity_);
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------


        NO_DISCARD std::string plasmaComponentName() const override { return "Electrons"; }  


        field_type* rho_{nullptr};
    };
} // namespace PHARE::core

#endif
