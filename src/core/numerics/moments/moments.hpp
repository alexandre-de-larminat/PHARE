#ifndef MOMENTS_HPP
#define MOMENTS_HPP

#include <iterator>

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/data/fermions/fermions.hpp"


namespace PHARE
{
namespace core
{
    template<typename Fermions>
    void resetMoments(Fermions& fermions, int picmodel) //TODO add some if-model-is-x?
    {
        auto& ions = fermions.ions;
        auto& electrons = fermions.electrons;

        for (auto& pop : ions)
        {
            pop.density().zero();
            pop.flux().zero();
        }
        for (auto& pop : electrons)
        {
            pop.density().zero();
            pop.flux().zero();
        }
    }
    template<typename Ions>
    void resetMoments(Ions& ions) //TODO add some if-model-is-x?
    {
        for (auto& pop : ions)
        {
            pop.density().zero();
            pop.flux().zero();
        }
    }

    struct DomainDeposit
    {
    };

    struct PatchGhostDeposit
    {
    };
    struct LevelGhostDeposit
    {
    };


    template<typename Fermions, typename GridLayout, typename DepositTag>
    void depositParticles(Fermions& fermions, GridLayout& layout,
                          Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate,
                          DepositTag, int picmodel)
    {
        auto& ions = fermions.ions;
        auto& electrons = fermions.electrons;

        for (auto& pop : ions)
        {
            auto& density = pop.density();
            auto& flux    = pop.flux();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            {
                auto& partArray = pop.domainParticles();
                interpolate(partArray, density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, PatchGhostDeposit>)
            {
                auto& partArray = pop.patchGhostParticles();
                interpolate(partArray, density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, LevelGhostDeposit>)
            {
                auto& partArray = pop.levelGhostParticlesOld();
                interpolate(partArray, density, flux, layout);
            }
        }
        for (auto& pop : electrons)
        {
            auto& density = pop.density();
            auto& flux    = pop.flux();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            {
                auto& partArray = pop.domainParticles();
                interpolate(partArray, density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, PatchGhostDeposit>)
            {
                auto& partArray = pop.patchGhostParticles();
                interpolate(partArray, density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, LevelGhostDeposit>)
            {
                auto& partArray = pop.levelGhostParticlesOld();
                interpolate(partArray, density, flux, layout);
            }
        }
    }
        template<typename Ions, typename GridLayout, typename DepositTag>
    void depositParticles(Ions& ions, GridLayout& layout,
                          Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate,
                          DepositTag)
    {
        for (auto& pop : ions)
        {
            auto& density = pop.density();
            auto& flux    = pop.flux();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            {
                auto& partArray = pop.domainParticles();
                interpolate(partArray, density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, PatchGhostDeposit>)
            {
                auto& partArray = pop.patchGhostParticles();
                interpolate(partArray, density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, LevelGhostDeposit>)
            {
                auto& partArray = pop.levelGhostParticlesOld();
                interpolate(partArray, density, flux, layout);
            }
        }
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_HPP
