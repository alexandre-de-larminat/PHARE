#ifndef PHARE_ION_POPULATION_HPP
#define PHARE_ION_POPULATION_HPP

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <array>


#include "core/def.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "initializer/data_provider.hpp"
#include "particle_pack.hpp"
#include "core/data/pic_electrons/particle_population.hpp"

namespace PHARE::core
{
    template<typename ParticleArray, typename VecField, typename TensorField, typename GridLayout>
    class IonPopulation: public ParticlePopulation<ParticleArray, VecField, TensorField, GridLayout, IonPopulation>
    {
    public:
        using field_type                       = typename VecField::field_type;
        using particle_array_type              = ParticleArray;
        using particle_resource_type           = ParticlesPack<ParticleArray>;
        using vecfield_type                    = VecField;
        using tensorfield_type                 = TensorField;


        IonPopulation(initializer::PHAREDict const& initializer)
            : name_{initializer["name"].template to<std::string>()}
            , mass_{initializer["mass"].template to<double>()}
            , flux_{name_ + "_flux", HybridQuantity::Vector::V}
            , momentumTensor_{name_ + "_momentumTensor", HybridQuantity::Tensor::M}
            , particleInitializerInfo_{initializer["particle_initializer"]}
        {
        }


        NO_DISCARD double mass() const { return mass_; }

        NO_DISCARD std::string const& name() const override { return name_; }


        NO_DISCARD auto const& particleInitializerInfo() const { return particleInitializerInfo_; }



        NO_DISCARD bool isUsable() const override
        {
            return particles_ != nullptr && rho_ != nullptr && flux_.isUsable()
                   && momentumTensor_.isUsable();
        }


        NO_DISCARD bool isSettable() const override
        {
            return particles_ == nullptr && rho_ == nullptr && flux_.isSettable()
                   && momentumTensor_.isSettable();
        }



        NO_DISCARD VecField const& flux() const override { return flux_; }
        NO_DISCARD VecField& flux() { return flux_; }

        NO_DISCARD TensorField const& momentumTensor() const { return momentumTensor_; }
        NO_DISCARD TensorField& momentumTensor() { return momentumTensor_; }

        NO_DISCARD field_type const& rho() const override { return *rho_; }

        NO_DISCARD ParticlesPack<ParticleArray> const* getParticlesPtr() const override { 
            return particles_; }

        //-------------------------------------------------------------------------
        //                  start the ResourcesUser interface
        //-------------------------------------------------------------------------



        using typename ParticlePopulation<ParticleArray, VecField, TensorField, GridLayout, 
        IonPopulation>::MomentProperties;


        NO_DISCARD MomentProperties getFieldNamesAndQuantities() const
        {
            return {{{name_ + "_rho", HybridQuantity::Scalar::rho}}};
        }




        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(flux_, momentumTensor_);
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

        //-------------------------------------------------------------------------
        //                  ends the ResourcesUser interface
        //-------------------------------------------------------------------------

        std::string getParticleName() const override { return "Ion"; }


    private:
        std::string name_;
        double mass_;
        VecField flux_;
        TensorField momentumTensor_;
        field_type* rho_{nullptr};
        ParticlesPack<ParticleArray>* particles_{nullptr};
        initializer::PHAREDict const& particleInitializerInfo_;
    };
} // namespace PHARE::core

#endif
