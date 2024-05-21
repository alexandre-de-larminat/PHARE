#ifndef PHARE_CORE_NUMERICS_MAXWELL_AMPERE_HPP
#define PHARE_CORE_NUMERICS_MAXWELL_AMPERE_HPP

#include <cstddef>
#include <cstdio>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include "initializer/data_provider.hpp"


namespace PHARE::core
{
template<typename GridLayout>
class MaxwellAmpere : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:

    explicit MaxwellAmpere(PHARE::initializer::PHAREDict const& dict)
        : c_norm{dict["normalized_c"].template to<double>()}
    {
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField const& J, VecField& Enew, double dt)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - MaxwellAmpere - GridLayout not set, cannot proceed to calculate MaxwellAmpere()");

        if (!(B.isUsable() && E.isUsable() && J.isUsable() && Enew.isUsable()))
            throw std::runtime_error("Error - MaxwellAmpere - not all VecField parameters are usable");

        this->dt_ = dt;

        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Ex = E(Component::X);
        auto const& Ey = E(Component::Y);
        auto const& Ez = E(Component::Z);

        auto const& Jx = J(Component::X);
        auto const& Jy = J(Component::Y);
        auto const& Jz = J(Component::Z);

        auto& Exnew = Enew(Component::X);
        auto& Eynew = Enew(Component::Y);
        auto& Eznew = Enew(Component::Z);

        layout_->evalOnBox(Exnew,
                                [&](auto&... args) mutable { ExEq_(Ex, B, Jx, Exnew, args...); });
        layout_->evalOnBox(Eynew,
                                [&](auto&... args) mutable { EyEq_(Ey, B, Jy, Eynew, args...); });
        layout_->evalOnBox(Eznew,
                                [&](auto&... args) mutable { EzEq_(Ez, B, Jz, Eznew, args...); });
    }


private:
    double dt_;
    double c_norm;
    double c2 = c_norm * c_norm;
    double mu0 = 1.2566370614e-6;

    template<typename VecField, typename Field, typename... Indexes>
    void ExEq_(Field const& Ex, VecField const& B, Field const& Jx, Field& Exnew, Indexes const&... ijk) const
    {
        auto const& [_, By, Bz] = B();

        if constexpr (dimension == 1)
        {
            Exnew(ijk...) = Ex(ijk...) - dt_ * c2 * mu0 * Jx(ijk...) ; 
        }
        if constexpr (dimension == 2)
        {
            Exnew(ijk...) = Ex(ijk...) + dt_ * c2 * (layout_->template deriv<Direction::Y>(Bz, {ijk...}) 
                            - Jx(ijk...) );
        }
        if constexpr (dimension == 3)
        {
            Exnew(ijk...) = Ex(ijk...) + dt_ * c2 * (layout_->template deriv<Direction::Y>(Bz, {ijk...})
                            - layout_->template deriv<Direction::Z>(By, {ijk...}) - Jx(ijk...));
        }   
    }

    template<typename VecField, typename Field, typename... Indexes>
    void EyEq_(Field const& Ey, VecField const& B, Field const& Jy, Field& Eynew, Indexes const&... ijk) const
    {
        auto const& [Bx, _, Bz] = B();

        if constexpr (dimension == 1 || dimension == 2)
        {
            Eynew(ijk...) = Ey(ijk...) - dt_ * c2 *( layout_->template deriv<Direction::X>(Bz, {ijk...}) 
                            + mu0 * Jy(ijk...));
        }
        if constexpr (dimension == 3)
        {
            Eynew(ijk...) = Ey(ijk...) + dt_ * c2 * (layout_->template deriv<Direction::Z>(Bx, {ijk...})
                            - layout_->template deriv<Direction::X>(Bz, {ijk...}) - Jy(ijk...));
        }
    }

    template<typename VecField, typename Field, typename... Indexes>
    void EzEq_(Field const& Ez, VecField const& B, Field const& Jz, Field& Eznew, Indexes const&... ijk) const
    {
        auto const& [Bx, By, _] = B();

        if constexpr (dimension == 1)
        {
            Eznew(ijk...) = Ez(ijk...) + dt_ * c2 * (layout_->template deriv<Direction::X>(By, {ijk...}) 
                            - mu0 * Jz(ijk...));
        }

        if constexpr (dimension == 2 || dimension == 3)
        {
            Eznew(ijk...) = Ez(ijk...) + dt_ * c2 * (layout_->template deriv<Direction::X>(By, {ijk...})
                            - layout_->template deriv<Direction::Y>(Bx, {ijk...}) - Jz(ijk...));
        }
    }
};

} // namespace PHARE::core

#endif
