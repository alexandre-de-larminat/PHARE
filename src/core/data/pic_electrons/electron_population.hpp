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
    class ElectronPopulation : public ParticlePopulation<ParticleArray, VecField, TensorField, GridLayout, ElectronPopulation>
    {
    public:
        using field_type                       = typename VecField::field_type;
        using particle_array_type              = ParticleArray;
        using particle_resource_type           = ParticlesPack<ParticleArray>;
        using vecfield_type                    = VecField;

        ElectronPopulation(initializer::PHAREDict const& initializer)
            : ParticlePopulation{initializer["name"].template to<std::string>()}
            , flux_{name_ + "_flux", HybridQuantity::Vector::Ve}
        {
        } 

        NO_DISCARD auto const& particleInitializerInfo() const { return particleInitializerInfo_; }


        NO_DISCARD VecField const& flux() const { return flux_; }
        NO_DISCARD VecField& flux() { return flux_; }

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


        using typename ParticlePopulation<ParticleArray, VecField, TensorField, GridLayout, 
        ElectronPopulation>::MomentProperties;


        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{name_ + "_rho", HybridQuantity::Scalar::rhoE}}};
        }



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
            {
                rho_ = field;
            }
            else
            {
                throw std::runtime_error("Error - invalid density buffer name");
            }
        }



        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(flux_);
        }
        

        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        std::string getParticleTypeName() const override { return "Electron"; }


    private:
        VecField flux_;
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
    };
} // namespace PHARE::core

#endif
