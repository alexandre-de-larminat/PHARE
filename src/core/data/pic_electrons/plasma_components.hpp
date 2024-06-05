#ifndef PHARE_PLASMA_COMPONENTS_HPP
#define PHARE_PLASMA_COMPONENTS_HPP

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
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"
#include "core/utilities/algorithm.hpp"

namespace PHARE::core
{
    template<typename ParticlePopulation, typename GridLayout>
    class PlasmaComponents
    {
    public:
        using field_type          = typename ParticlePopulation::field_type;
        using vecfield_type       = typename ParticlePopulation::vecfield_type;
        using Float               = typename field_type::type;

        static constexpr auto dimension = GridLayout::dimension;


        explicit PlasmaComponents(PHARE::initializer::PHAREDict const& dict)
            : populations_{generate(
                  [&dict](auto ipop) { //
                      return ParticlePopulation{dict["pop" + std::to_string(ipop)]};
                  },
                  dict["nbrPopulations"].template to<std::size_t>())}
            , bulkVelocity_{bulkVelocityName(), HybridQuantity::Vector::V}
            , sameMasses_{allSameMass_()}
        {
        }


        NO_DISCARD std::size_t nbrPopulations() const { return populations_.size(); }

        NO_DISCARD vecfield_type const& velocity() const { return bulkVelocity_; }

        NO_DISCARD vecfield_type& velocity() { return bulkVelocity_; }

        NO_DISCARD std::string massDensityName() const { return "massDensity"; }

        virtual std::string densityName() const { return ""; }

        virtual std::string bulkVelocityName() const { return ""; }
                
        virtual field_type* rhoPtr() { return rho_; }


        virtual bool isUsable() const{ return false; }

        virtual bool isSettable() const { return false; }


        NO_DISCARD field_type const& density() const
        {
            if (isUsable())
                return *rhoPtr();
            else
                throw std::runtime_error("Error - cannot access density data");
        }

        NO_DISCARD field_type& density()
        {
            if (isUsable())
                return *rhoPtr();
            else
                throw std::runtime_error("Error - cannot access density data");
        }

        
        void computeDensity()
        {
            rhoPtr()->zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.density();
                std::transform(std::begin(*rhoPtr()), std::end(*rhoPtr()), std::begin(popDensity),
                               std::begin(*rhoPtr()), std::plus<Float>{});
            }
        }


        NO_DISCARD field_type const& massDensity() const
        {
            if (isUsable())
                return sameMasses_ ? *rhoPtr() : *massDensity_; 
            throw std::runtime_error("Error - cannot access density data");
        }

        NO_DISCARD field_type& massDensity()
        {
            if (isUsable())
                return sameMasses_ ? *rhoPtr() : *massDensity_;
            else
                throw std::runtime_error("Error - cannot access density data");
        }


        void computeMassDensity()
        {
            massDensity_->zero();

            for (auto const& pop : populations_)
            {
                // we sum over all nodes contiguously, including ghosts
                // nodes. This is more efficient and easier to code as we don't
                // have to account for the field dimensionality.

                auto& popDensity = pop.density();
                std::transform(
                    std::begin(*massDensity_), std::end(*massDensity_), std::begin(popDensity),
                    std::begin(*massDensity_),
                    [&pop](auto const& n, auto const& pop_n) { return n + pop_n * pop.mass(); });
            }
        }



        void computeBulkVelocity()
        {
            // the bulk velocity is sum(pop_mass * pop_flux) / sum(pop_mass * pop_density)
            // if all populations have the same mass, this is equivalent to sum(pop_flux) /
            // sum(pop_density) sum(pop_density) is rho_ and already known by the time we get here.
            // sum(pop_mass * pop_flux) is massDensity_ and is computed by computeMassDensity() if
            // needed
            if (!sameMasses_)
                computeMassDensity();
            auto const& density = (sameMasses_) ? rhoPtr() : massDensity_;

            bulkVelocity_.zero();
            auto& vx = bulkVelocity_.getComponent(Component::X);
            auto& vy = bulkVelocity_.getComponent(Component::Y);
            auto& vz = bulkVelocity_.getComponent(Component::Z);

            for (auto& pop : populations_)
            {
                // account for mass only if populations have different masses
                std::function<Float(Float, Float)> plus = std::plus<Float>{};
                std::function<Float(Float, Float)> plusMass
                    = [&pop](Float const& v, Float const& f) { return v + f * pop.mass(); };

                auto const& flux    = pop.flux();
                auto&& [fx, fy, fz] = flux();


                std::transform(std::begin(vx), std::end(vx), std::begin(fx), std::begin(vx),
                               sameMasses_ ? plus : plusMass);
                std::transform(std::begin(vy), std::end(vy), std::begin(fy), std::begin(vy),
                               sameMasses_ ? plus : plusMass);
                std::transform(std::begin(vz), std::end(vz), std::begin(fz), std::begin(vz),
                               sameMasses_ ? plus : plusMass);
            }


            std::transform(std::begin(vx), std::end(vx), std::begin(*density), std::begin(vx),
                           std::divides<Float>{});
            std::transform(std::begin(vy), std::end(vy), std::begin(*density), std::begin(vy),
                           std::divides<Float>{});
            std::transform(std::begin(vz), std::end(vz), std::begin(*density), std::begin(vz),
                           std::divides<Float>{});
        }


        NO_DISCARD auto begin() { return std::begin(populations_); }
        NO_DISCARD auto end() { return std::end(populations_); }
        NO_DISCARD auto begin() const { return std::begin(populations_); }
        NO_DISCARD auto end() const { return std::end(populations_); }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------


        struct MomentsProperty
        {
            std::string name;
            typename HybridQuantity::Scalar qty;
        };

        using MomentProperties = std::vector<MomentsProperty>;

        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            if (sameMasses_)
                return {{{densityName(), HybridQuantity::Scalar::rho}}};
            else
                return {{{densityName(), HybridQuantity::Scalar::rho},
                         {massDensityName(), HybridQuantity::Scalar::rho}}};
        }

        
        NO_DISCARD std::vector<ParticlePopulation>& getRunTimeResourcesUserList()
        {
            return populations_;
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------



        virtual std::string plasmaComponentName() const { return ""; }


        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << plasmaComponentName() << "\n";
            ss << "------------------------------------\n";
            ss << "number of populations  : " << nbrPopulations() << "\n";
            for (auto& pop : populations_)
                ss << core::to_str(pop);
            return ss.str();
        }

    protected:

        bool allSameMass_() const
        {
            return all(populations_, [this](auto const& pop) { // arbitrary small diff
                return float_equals(pop.mass(), populations_.front().mass(), /*abs_tol=*/1e-10);
            });
        }


        field_type* rho_{nullptr};
        std::vector<ParticlePopulation> populations_;
        vecfield_type bulkVelocity_;
        field_type* massDensity_{nullptr};

        // note this is set at construction when all populations are added
        // in the future if some populations are dynamically created during the simulation
        // this should be updated accordingly
        bool sameMasses_{false};
    };
} // namespace PHARE::core

#endif
