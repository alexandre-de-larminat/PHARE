#ifndef PHARE_FERMION_UPDATER_HPP
#define PHARE_FERMION_UPDATER_HPP


#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/projector/projector.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/fermions/fermions.hpp"

#include "initializer/data_provider.hpp"

#include "core/logger.hpp"

#include <cstddef>
#include <memory>


namespace PHARE::core
{

template<typename Ions, typename PICElectrons, typename Electromag, typename VecField, typename GridLayout>
class FermionUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
    using Projector         = PHARE::core::Projector<GridLayout>;
    using ParticleArray     = typename Ions::particle_array_type;
    using Particle_t        = typename ParticleArray::Particle_t;
    using PartIterator      = typename ParticleArray::iterator;
    using ParticleRange     = IndexRange<ParticleArray>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::Pusher<dimension, ParticleRange, Electromag, Interpolator,
                                       BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, ParticleRange, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;
    Projector projector_;

public:
    FermionUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
    }

    void updatePopulations(Ions& ions, PICElectrons& electrons, Electromag const& em, VecField& J, GridLayout const& layout, double dt);

    void updateFermions(Ions& ions, PICElectrons& electrons, GridLayout const& layout);

    void reset()
    {
        // clear memory
        // tmp_particles_ = std::move(ParticleArray{Box{}});
    }


private:
    void updateAndDepositAll_(Ions& ions, PICElectrons& electrons, Electromag const& em, VecField& J, GridLayout const& layout, double dt);


    // dealloced on regridding/load balancing coarsest
    // ParticleArray tmp_particles_{Box{}}; //{std::make_unique<ParticleArray>(Box{})};
};




template<typename Ions, typename PICElectrons, typename Electromag, typename VecField, typename GridLayout>
void FermionUpdater<Ions, PICElectrons, Electromag, VecField, GridLayout>::updatePopulations(Ions& ions, PICElectrons& electrons, 
                                                                 Electromag const& em, VecField& J,
                                                                 GridLayout const& layout, 
                                                                 double dt)
{
    PHARE_LOG_SCOPE("FermionUpdater::updatePopulations");

    resetMoments(ions); //TODO add electrons
    pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

    updateFermions(ions, electrons, layout);
    updateAndDepositAll_(ions, electrons, em, J, layout, dt);
}


template<typename Ions, typename PICElectrons, typename Electromag, typename VecField, typename GridLayout>
void FermionUpdater<Ions, PICElectrons, Electromag, VecField, GridLayout>::updateFermions(Ions& ions, PICElectrons& electrons, 
                                                                      GridLayout const& layout)
{
    ions.computeDensity();
    ions.computeBulkVelocity();
    electrons.computeDensity();
    electrons.computeBulkVelocity();
}



template<typename Ions, typename PICElectrons, typename Electromag, typename VecField, typename GridLayout>
void FermionUpdater<Ions, PICElectrons, Electromag, VecField, GridLayout>::updateAndDepositAll_(Ions& ions, PICElectrons& electrons, 
                                                                    Electromag const& em, VecField& J, 
                                                                    GridLayout const& layout,
                                                                    double dt)
{
    PHARE_LOG_SCOPE("FermionUpdater::updateAndDepositAll_");


    for (auto& pop : ions)
        {        
            auto& domain = pop.domainParticles();
            auto rangeIn  = makeIndexRange(domain);

            auto push = [&](auto&& rangeIn) {
            auto rangeOut = pusher_->move(rangeIn, rangeIn, em, pop.mass(), interpolator_, 
                                          layout, [](auto& rangeIn) {return rangeIn;}, 
                                          [](auto& rangeIn) {return rangeIn;});

            projector_(J, rangeOut, rangeIn, dt);

            };

            push(makeIndexRange(pop.patchGhostParticles()));
            push(makeIndexRange(pop.levelGhostParticles()));
            
            interpolator_(makeIndexRange(domain), pop.density(), pop.flux(), layout);
        }

    for (auto& pop : electrons)
        {        
            auto& domain = pop.domainParticles();
            auto rangeIn  = makeIndexRange(domain);

            auto push = [&](auto&& rangeIn) {
            auto rangeOut = pusher_->move(rangeIn, rangeIn, em, pop.mass(), interpolator_, 
                                          layout, [](auto& rangeIn) {return rangeIn;}, 
                                          [](auto& rangeIn) {return rangeIn;});

            projector_(J, rangeOut, rangeIn, dt);

            };

            push(makeIndexRange(pop.patchGhostParticles()));
            push(makeIndexRange(pop.levelGhostParticles()));
            
            interpolator_(makeIndexRange(domain), pop.density(), pop.flux(), layout);
        }
}


} // namespace PHARE::core

#endif // FERMION_UPDATER_HPP
