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
    template<typename ParticleArray, typename VecField, typename GridLayout, 
    template<typename, typename, typename> class PopulationType>
    class ParticlePopulation
    {
    public:
        using field_type                       = typename VecField::field_type;
        using pop_type                         = PopulationType<ParticleArray, VecField, GridLayout>;
        static constexpr std::size_t dimension = VecField::dimension;


        virtual std::string const& name() const { return name_; }


        virtual ParticlesPack<ParticleArray> const* getParticlesPtr() const { return particles_; }

        virtual VecField& flux() { return flux_; }
        virtual VecField const& flux() const { return flux_ ; }
        virtual field_type const* rhoPtr() const { return rho_; }
        virtual field_type const& rho() const { return *rho_; }

        virtual bool isUsable() const
        {
        }


        virtual bool isSettable() const
        {
        }

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

        void setBuffer(std::string const& bufferName, ParticlesPack<ParticleArray>* pack)
        {
            if (bufferName == name())
                getParticlesPtr() = pack;
            else
                throw std::runtime_error("Error - invalid particle resource name");
        }


        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        virtual std::string getParticleName() const { return ""; }

        NO_DISCARD std::string to_str()
        {
            std::stringstream ss;
            ss << getParticleName() << "Population\n";
            ss << "------------------------------------\n";
            ss << "name                : " << name() << "\n";
            return ss.str();
        }

    private:
        std::string name_;
        VecField flux_{"Placeholder", HybridQuantity::Vector::V};
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
    };
} // namespace PHARE::core

#endif
