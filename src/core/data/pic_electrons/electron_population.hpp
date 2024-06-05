#ifndef PHARE_ELECTRON_POPULATION_HPP
#define PHARE_ELECTRON_POPULATION_HPP

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"
#include "particle_population.hpp"

namespace PHARE::core
{
    template<typename ParticleArray, typename VecField, typename TensorField, typename GridLayout>
    class ElectronPopulation : public ParticlePopulation<ParticleArray, VecField, TensorField, 
                                                         GridLayout, ElectronPopulation>
    {
    public:
        using field_type                       = typename VecField::field_type;
        using particle_array_type              = ParticleArray;
        using particle_resource_type           = ParticlesPack<ParticleArray>;
        using vecfield_type                    = VecField;

        using Super = ParticlePopulation<ParticleArray, VecField, TensorField, GridLayout, 
                                         ElectronPopulation>;
        using Super::name_;
        using Super::flux_;

        ElectronPopulation(initializer::PHAREDict const& initializer)
            : Super{initializer}
        {
        } 


        NO_DISCARD field_type const& rho() const override { return *rho_; }

        NO_DISCARD ParticlesPack<ParticleArray> const* getParticlesPtr() const override { 
            return particles_; }


        NO_DISCARD bool isUsable() const
        {
            return particles_ != nullptr && rho_ != nullptr && flux_.isUsable();
        }


        NO_DISCARD bool isSettable() const
        {
            return particles_ == nullptr && rho_ == nullptr && flux_.isSettable();
        }



        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        void setBuffer(std::string const& bufferName, ParticlesPack<ParticleArray>* pack)
        {
            if (bufferName == name_)
                particles_ = pack;
            else
                throw std::runtime_error("Error - invalid particle resource name");
        }

        void setBuffer(std::string const& bufferName, field_type* field)
        {
            if (bufferName == name_ + "_rho")
                rho_ = field;
            else
                throw std::runtime_error("Error - invalid density buffer name");
        }



        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(flux_);
        }
        

        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        std::string getParticleName() const override { return "Electron"; }


    private:
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
    };
} // namespace PHARE::core

#endif
