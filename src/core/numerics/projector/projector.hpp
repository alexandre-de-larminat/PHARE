#ifndef PHARE_CORE_NUMERICS_PROJECTOR_HPP
#define PHARE_CORE_NUMERICS_PROJECTOR_HPP

#include <cmath>
#include <iostream>
#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/particles/particle.hpp"

#include "core/numerics/interpolator/interpolator.hpp"


namespace PHARE::core
{

// Calculate weights for a single particle, for use in calculating Esirkepov coefficients

template<typename Particle>
auto E_weights_(Particle const& part)
{   
    constexpr auto interpOrder = 1 ; //TODO. just for testing
    constexpr auto supportpts  = interpOrder + 1;

    constexpr auto dimension = 1; //part.iCell.size();

    using Interpolator_t = PHARE::core::Interpolator<dimension, interpOrder>;
    using Weighter_t     = PHARE::core::Weighter<interpOrder>;
    Weighter_t weighter_;
    std::array<double, supportpts> weights;
    std::vector<std::vector<double>> vec_weights;
    
    for (uint i = 0; i < dimension; ++i)
    {
        auto position = part.iCell[i] + part.delta[i] - 0.5;
        auto startIndex = part.iCell[i] - Interpolator_t::template computeStartLeftShift<QtyCentering, QtyCentering::dual>(part.delta[i]);
        weighter_.computeWeight(position, startIndex, weights);

        // Add 0.0 to the beginning and end of the weights for calculation purposes
        std::vector<double> vec_weights_i(std::begin(weights), std::end(weights));
        vec_weights_i.insert(vec_weights_i.begin(), 0.0);
        vec_weights_i.insert(vec_weights_i.end(), 0.0);

        vec_weights.push_back(vec_weights_i);

    }

    return vec_weights; 
   
} // END E_weights_



template<std::size_t dim>
class ProjectJ
{
};


/** \brief specialization of ParticleToMesh for 1D interpolation of J */
template<>
class ProjectJ<1>
{
public:
    template<typename VecField, typename Particle>
    inline void operator()(VecField& J, Particle const& partIn,
                           Particle const& partOut, double dt)
    {
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const& xStartIndex = partIn.iCell[0] ; // central dual node
        auto const& S0          = E_weights_(partIn)[0];
        auto const& S1          = E_weights_(partOut)[0];
        auto const& order_size  = S0.size();

        // requisite for appropriate centering (greatest weight at node considered)
        int iCorr = order_size/2;
        if (S0[1] > S0[S0.size()-2]) { 
            iCorr -= 1;
        }

        double charge_weight = partIn.charge * partIn.weight;
            
        double cr_p = charge_weight/dt;   // current density in the evaluated dimension, assuming dx=1
        double cry_p_1D = charge_weight*partIn.v[1];   // current density in the y-direction in 1D
        double crz_p_1D2D = charge_weight*partIn.v[2]; // current density in the z-direction in 1D or 2D

        std::vector<double> Jx_p(order_size, 0.);

        std::vector<double> Wl(order_size);
        std::vector<double> Wt(order_size);


        for (uint i = 0; i < order_size; ++i)
        {
            auto x = xStartIndex + i - iCorr; // eg, i from -2 to 2 for 3rd order B-splines.

            Wl[i] = S0[i] - S1[i];
            Wt[i] = 0.5 * ( S0[i] + S1[i] );

            Jx_p[i] = Jx_p[i-1] + cr_p * Wl[i-1];
            Jx(x) += Jx_p[i];
            Jy(x)  += cry_p_1D * Wt[i];
            Jz(x)  += crz_p_1D2D * Wt[i];
        }
    }
}; // END ProjectJ<1> specialization


template<>
class ProjectJ<2>
{
public:
    template<typename VecField, typename Particle>
    inline void operator()(VecField& J, Particle const& partIn,
                           Particle const& partOut, double dt)
    {
        double one_third_ = 1./3.;
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const& xStartIndex = partIn.iCell[0];
        auto const& yStartIndex = partIn.iCell[1];
        auto const& oldWeights  = E_weights_(partIn);
        auto const& newWeights  = E_weights_(partOut);

        auto const& Sx0      = oldWeights[0];
        auto const& Sy0      = oldWeights[1];
        auto const& Sx1      = newWeights[0];
        auto const& Sy1      = newWeights[1];

        auto const& order_size = Sx0.size();

        // requisite for appropriate centering
        int iCorr = order_size/2;

        if (Sx0[1] > Sx0[Sx0.size()-2]) { 
            iCorr -= 1;
        }
        int jCorr = order_size/2;
        if (Sy0[1] > Sy0[Sy0.size()-2]) {
            jCorr -= 1;
        }

        double charge_weight = partIn.charge * partIn.weight; // CHECK weight factors in the cell volume
     
        double cr_p = charge_weight/dt;  // current density in the evaluated dimension
        double crz_p_1D2D = charge_weight*partIn.v[2]; // current density in the z-direction in 1D or 2D

        std::vector<std::vector<double>> Jx_p(order_size, std::vector<double>(order_size, 0.));
        std::vector<std::vector<double>> Jy_p(order_size, std::vector<double>(order_size, 0.));

        std::vector<double> DSx(order_size);
        std::vector<double> DSy(order_size);

        std::vector<std::vector<double>> Wx(order_size, std::vector<double>(order_size));
        std::vector<std::vector<double>> Wy(order_size, std::vector<double>(order_size));
        std::vector<std::vector<double>> Wz(order_size, std::vector<double>(order_size));

        for(auto i = 0u; i < order_size; ++i)
        {
            DSx[i] = Sx1[i] - Sx0[i];
            DSy[i] = Sy1[i] - Sy0[i];
        }


        for(auto i = 0u; i < order_size; ++i)
        {
            for(auto j = 0u; j < order_size; ++j)
            {
                Wx[i][j] = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
                Wy[i][j] = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
                Wz[i][j] = Sx0[i] * Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*Sx0[i]*DSy[j] + one_third_*DSx[i]*DSy[j];
            }
        }


        for (auto i = 0u; i < order_size; ++i)
        {
            for(auto j = 0u; j < order_size; ++j)
            {
                auto x = xStartIndex + i - iCorr; // eg, i from -2 to 2 for 3rd order B-splines.
                auto y = yStartIndex + j - jCorr;

                Jx_p[i][j] = Jx_p[i-1][j] + cr_p * Wx[i-1][j];
                Jx(x, y) += Jx_p[i][j] ;

                Jy_p[i][j] = Jy_p[i][j-1] + cr_p * Wy[i][j-1];
                Jy(x, y) += Jy_p[i][j] ;

                Jz(x, y) += crz_p_1D2D * Wz[i][j];
            }
        }

    }
}; // END ProjectJ<2> specialization




template<>
class ProjectJ<3>
{
public:
    double one_third_ = 1./3.;
    
    template<typename VecField, typename Particle>
    inline void operator()(VecField& J, Particle const& partIn,
                           Particle const& partOut, double dt)
    {
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const& xStartIndex = partIn.iCell[0];
        auto const& yStartIndex = partIn.iCell[1];
        auto const& zStartIndex = partIn.iCell[2];

        auto const& oldWeights  = E_weights_(partIn);
        auto const& newWeights  = E_weights_(partOut);

        auto const& Sx0 = oldWeights[0];
        auto const& Sy0 = oldWeights[1];
        auto const& Sz0 = oldWeights[2];
        auto const& Sx1 = newWeights[0];
        auto const& Sy1 = newWeights[1];
        auto const& Sz1 = newWeights[2];

        auto const& order_size  = Sx0.size();

        // requisite for appropriate centering
        int iCorr = order_size/2;
        if (Sx0[1] > Sx0[Sx0.size()-2]) { 
            iCorr -= 1;
        }
        int jCorr = order_size/2;
        if (Sy0[1] > Sy0[Sy0.size()-2]) {
            jCorr -= 1;
        }
        int kCorr = order_size/2;
        if (Sz0[1] > Sz0[Sz0.size()-2]) {
            kCorr -= 1;
        }

        double charge_weight = partIn.charge * partIn.weight; // CHECK: assumes weight factors in the cell volume
            
        double cr_p = charge_weight/dt; // current density in the evaluated dimension

        std::vector<std::vector<std::vector<double>>> Jx_p(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size, 0.)));

        std::vector<std::vector<std::vector<double>>> Jy_p(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size, 0.)));

        std::vector<std::vector<std::vector<double>>> Jz_p(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size, 0.)));
        
        std::vector<double> DSx(order_size);
        std::vector<double> DSy(order_size);
        std::vector<double> DSz(order_size);

        std::vector<std::vector<std::vector<double>>> Wx(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size)));

        std::vector<std::vector<std::vector<double>>> Wy(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size)));

        std::vector<std::vector<std::vector<double>>> Wz(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size)));

        for(auto i = 0u; i < order_size; ++i)
        {
            DSx[i] = Sx1[i] - Sx0[i];
            DSy[i] = Sy1[i] - Sy0[i];
            DSz[i] = Sz1[i] - Sz0[i];
        }
        
        for(auto i = 0u; i < order_size; ++i)
        {
            for(auto j = 0u; j < order_size; ++j)
            {
                for(auto k = 0u; k < order_size; ++k)
                {
                Wx[i][j] = DSx[i] * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*Sy0[j]*DSz[k] + one_third_*DSy[j]*DSz[k]);
                Wy[i][j] = DSy[j] * (Sx0[i]*Sz0[k] + 0.5*DSx[i]*Sz0[k] + 0.5*Sx0[i]*DSz[k] + one_third_*DSx[i]*DSz[k]);
                Wz[i][j] = DSz[k] * (Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*Sx0[i]*DSy[j] + one_third_*DSx[i]*DSy[j]);
                }
            }
        }

        for (auto i = 0u; i < order_size; ++i)
        {
            for (auto j = 0u; j < order_size; ++j)
            {
                for (auto k = 0u; k < order_size; ++k)
                {
                auto x = xStartIndex + i - iCorr; // eg, i from -2 to 2 for 3rd order B-splines.
                auto y = yStartIndex + j - jCorr;
                auto z = zStartIndex + k - kCorr;

                Jx_p[i][j][k] = Jx_p[i-1][j][k] + cr_p * Wx[i-1][j][k];
                Jx(x, y, z) += Jx_p[i][j][k] ;

                Jy_p[i][j][k] = Jy_p[i][j-1][k] + cr_p * Wy[i][j-1][k];
                Jy(x, y, z) += Jy_p[i][j][k] ;

                Jz_p[i][j][k] = Jz_p[i][j][k-1] + cr_p * Wz[i][j][k-1];
                Jz(x, y, z) += Jz_p[i][j][k] ;

                }
            }
        }
    }
}; // END ProjectJ<3> specialization


template<typename GridLayout>
class Projector : public LayoutHolder<GridLayout>
{

public:
    constexpr static auto interp_order = GridLayout::interp_order;
    constexpr static auto dimension    = GridLayout::dimension;
    ProjectJ<dimension> projectJ_;

    template<typename VecField, typename ParticleRange>
    inline void operator()(VecField& J, ParticleRange& rangeOut, ParticleRange& rangeIn, double dt)
    {
        auto& inParticles = rangeIn.array();
        auto& outParticles = rangeOut.array();

        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
                    ++inIdx, ++outIdx)
                {
                    projectJ_(J, outParticles[outIdx], inParticles[inIdx], dt);
                }   
    }

}; // END Projector

} //namespace PHARE::core 
#endif
