#ifndef MOMENTS_HPP
#define MOMENTS_HPP

#include <iterator>

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/pic_electrons/pic_electrons.hpp"


namespace PHARE
{
namespace core
{
    template<typename Ions, typename PICElectrons>
    void resetMoments(Ions& ions, PICElectrons& electrons)
    {
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
    void resetMoments(Ions& ions)
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


    template<typename Ions, typename PICElectrons, typename GridLayout, typename DepositTag>
    void depositParticles(Ions& ions, PICElectrons& electrons, GridLayout& layout,
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
