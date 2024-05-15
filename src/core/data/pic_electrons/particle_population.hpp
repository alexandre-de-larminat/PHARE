#ifndef PHARE_PARTICLE_POPULATION_HPP
#define PHARE_PARTICLE_POPULATION_HPP

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <array>


#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"

namespace PHARE::core
{
    template<typename ParticleArray, typename VecField, typename TensorField, typename GridLayout, 
    template<typename, typename, typename, typename> class PopulationType>
    class ParticlePopulation
    {
    public:
        using field_type                       = typename VecField::field_type;
        using pop_type                         = PopulationType<ParticleArray, VecField, TensorField, GridLayout>;
        static constexpr std::size_t dimension = VecField::dimension;

        ParticlePopulation(initializer::PHAREDict const& initializer)
            : name_{initializer["name"].template to<std::string>()}
            , mass_{initializer["mass"].template to<double>()}
            , particleInitializerInfo_{initializer["particle_initializer"]}
        {
        }

        NO_DISCARD double mass() const { return mass_; }

        NO_DISCARD std::string const& name() const { return name_; }

        NO_DISCARD initializer::PHAREDict const& particleInitializerInfo() const { return particleInitializerInfo_; }

        virtual ParticlesPack<ParticleArray> const* getParticlesPtr() const { return particles_; }


        virtual field_type const& rho() const { return *rho_; }

        virtual bool isUsable() const{ return false; }

        virtual bool isSettable() const { return false; }

        NO_DISCARD field_type const& density() const
        {
            if (isUsable())
            {
                return rho();
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to density field");
            }
        }
        
        NO_DISCARD field_type& density()
        {
            return const_cast<field_type&>(static_cast<const pop_type*>(this)->density());
        }

        
        NO_DISCARD auto nbrParticles() const
        {
            if (isUsable())
            {
                return getParticlesPtr()->domainParticles->size();
            }
            else
            {
                throw std::runtime_error("Error - cannot access to particles");
            }
        }


        NO_DISCARD auto& domainParticles() const
        {
            if (isUsable())
            {
                return *getParticlesPtr()->domainParticles;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }


        NO_DISCARD auto& domainParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<pop_type const*>(this)->domainParticles());
        }


        NO_DISCARD auto& patchGhostParticles() const
        {
            if (isUsable())
            {
                return *getParticlesPtr()->patchGhostParticles;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }

        NO_DISCARD auto& patchGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<pop_type const*>(this)->patchGhostParticles());
        }


        NO_DISCARD auto& levelGhostParticles() const
        {
            if (isUsable())
            {
                return *getParticlesPtr()->levelGhostParticles;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }

        NO_DISCARD auto& levelGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<pop_type const*>(this)->levelGhostParticles());
        }


        NO_DISCARD ParticleArray& levelGhostParticlesOld()
        {
            if (isUsable())
            {
                return *getParticlesPtr()->levelGhostParticlesOld;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }

        NO_DISCARD ParticleArray& levelGhostParticlesNew()
        {
            if (isUsable())
            {
                return *getParticlesPtr()->levelGhostParticlesNew;
            }
            else
            {
                throw std::runtime_error("Error - cannot provide access to particle buffers");
            }
        }


        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        struct MomentsProperty
        {
            std::string name;
            typename HybridQuantity::Scalar qty;
        };

       using MomentProperties = std::array<MomentsProperty, 1>;


        struct ParticleProperty
        {
            std::string name;
        };

        using ParticleProperties = std::array<ParticleProperty, 1>;



        NO_DISCARD ParticleProperties getParticleArrayNames() const { return {{{ name() }}}; }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        virtual std::string getParticleTypeName() const { return ""; }

        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << getParticleTypeName() << "Population\n";
            ss << "------------------------------------\n";
            ss << "name                : " << name() << "\n";
            return ss.str();
        }

    protected:
        std::string name_;
        double mass_;
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
        initializer::PHAREDict const& particleInitializerInfo_;
    };
} // namespace PHARE::core

#endif
