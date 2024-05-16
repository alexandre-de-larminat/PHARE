#ifndef PHARE_IONS_HPP
#define PHARE_IONS_HPP

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
#include "particle_initializers/particle_initializer_factory.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/data/pic_electrons/plasma_components.hpp"

namespace PHARE::core
{
    template<typename IonPopulation, typename GridLayout>
    class Ions : public PlasmaComponents<IonPopulation, GridLayout>
    {
    public:
        using value_type          = IonPopulation;
        using field_type          = typename IonPopulation::field_type;
        using vecfield_type       = typename IonPopulation::vecfield_type;
        using Float               = typename field_type::type;
        using tensorfield_type    = typename IonPopulation::tensorfield_type;
        using particle_array_type = typename IonPopulation::particle_array_type;
        using ParticleInitializerFactoryT
            = ParticleInitializerFactory<particle_array_type, GridLayout>;
        using gridlayout_type           = GridLayout;

        using Super = PlasmaComponents<IonPopulation, GridLayout>;
        using Super::bulkVelocity_;
        using Super::sameMasses_;
        using Super::populations_;
        using Super::massDensity_;


        explicit Ions(PHARE::initializer::PHAREDict const& dict)
            : Super{dict}
            , momentumTensor_{"momentumTensor", HybridQuantity::Tensor::M}
        {
        }


        NO_DISCARD field_type* rhoPtr() override { return rho_; }

        NO_DISCARD std::string densityName() const override { return "rho"; }

        NO_DISCARD std::string bulkVelocityName() const override { return "bulkVel"; }

        tensorfield_type const& momentumTensor() const { return momentumTensor_; }

        tensorfield_type& momentumTensor() { return momentumTensor_; }



        void computeFullMomentumTensor()
        {
            momentumTensor_.zero();
            auto& mom = momentumTensor_;

            for (auto& pop : populations_)
            {
                auto& p_mom = pop.momentumTensor();
                for (auto p_mij = p_mom.begin(), mij = mom.begin(); p_mij != p_mom.end();
                     ++p_mij, ++mij)
                {
                    std::transform(std::begin(**mij), std::end(**mij), std::begin(**p_mij),
                                   std::begin(**mij), std::plus<typename field_type::type>{});
                }
            }
        }


        // in the following isUsable and isSettable the massDensity_ is not checked
        // because it is for internal use only so no object will ever need to access it.
        NO_DISCARD bool isUsable() const override
        {
            bool usable
                = rho_ != nullptr and bulkVelocity_.isUsable() and momentumTensor_.isUsable();

            // if all populations have the same mass, we don't need the massDensity_
            usable &= (sameMasses_) ? true : massDensity_ != nullptr;

            for (auto const& pop : populations_)
            {
                usable = usable and pop.isUsable();
            }
            return usable;
        }



        NO_DISCARD bool isSettable() const override
        {
            bool settable
                = rho_ == nullptr and bulkVelocity_.isSettable() and momentumTensor_.isSettable();

            // if all populations have the same mass, we don't need the massDensity_
            settable &= (sameMasses_) ? true : massDensity_ == nullptr;

            for (auto const& pop : populations_)
            {
                settable = settable and pop.isSettable();
            }
            return settable;
        }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------


        void setBuffer(std::string const& bufferName, field_type* field)
        {
            if (bufferName == densityName())
            {
                rho_ = field;
            }
            else if (bufferName == Super::massDensityName())
            {
                assert(sameMasses_ == false);
                massDensity_ = field;
            }
            else
            {
                throw std::runtime_error("Error - invalid density buffer name : " + bufferName);
            }
        }



        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(bulkVelocity_, momentumTensor_);
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------



        NO_DISCARD std::string plasmaComponentName() const override { return "Ions"; }

        NO_DISCARD bool sameMasses() const { return sameMasses_; }


    private:

        field_type* rho_{nullptr};
        tensorfield_type momentumTensor_;
    };
} // namespace PHARE::core

#endif
