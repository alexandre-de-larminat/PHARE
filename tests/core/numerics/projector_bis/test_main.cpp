#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/projector/projector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"

#include "core/data/particles/particle_array.hpp"

#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"

using namespace PHARE::core;


template<int dim, int interp>
class NDlayout
{
    NDlayout() {}

    using nDL = GridLayout<GridLayoutImplYee<dim, interp>>;

public:
    static nDL create()
    {
        if constexpr (dim == 1)
        {
            return {{{1}}, {{10}}, Point{0.}};
        }
        else if constexpr (dim == 2)
        {
            return {{{1, 0.5}}, {{10, 15}}, Point{0., 0.}};
        }
        else if constexpr (dim == 3)
        {
            return {{{1, 0.5, 2}}, {{10, 15, 20}}, Point{0., 0., 0.}};
        }
    }
};

class DummySelector
{
public:
    template<typename Range>
    Range operator()(Range& particles) const
    {
        return particles;
    }
};

class Electromag
{
};

class Interpolator_
{
    using E_B_tuple = std::tuple<std::array<double, 3>, std::array<double, 3>>;

public:
    template<typename Particle_t, typename Electromag, typename GridLayout>
    auto operator()(Particle_t& particle, Electromag const&, GridLayout&)
    {
        E_B_tuple eb_interop;
        auto& [pE, pB]        = eb_interop;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;

        pEx = 0.01;
        pEy = -0.05;
        pEz = 0.05;
        pBx = 1.;
        pBy = 1.;
        pBz = 1.;

        return eb_interop;
    }
};

template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct EsirkepovTest : public ::testing::Test
{
public:
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee  = GridLayout<GridLayoutImplYee<dim, interp>>;
    GridYee layout = NDlayout<dim, interp>::create();
    
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jz;
    VecField<NdArrayVector<dim>, HybridQuantity> J;
    ParticleArray<dim> particlesIn;
    ParticleArray<dim> particlesOut;
    Projector<GridYee> projector;


    using Pusher_ = BorisPusher<dim, IndexRange<ParticleArray<dim>>, Electromag, Interpolator_,
                                BoundaryCondition<dim, 1>, GridLayout<GridLayoutImplYee<dim, interp>>>;


    std::unique_ptr<Pusher_> pusher;
    double mass;
    double dt;
    Electromag em;
    Interpolator_ interpolator;
    DummySelector selector;

    std::array<double, dim> dxyz;


    EsirkepovTest()
        : Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , J{"J", HybridQuantity::Vector::J}
        , particlesIn{layout.AMRBox()}
        , particlesOut(layout.AMRBox())
        , pusher{std::make_unique<Pusher_>()}
        , mass{1}
        , dt{1}

    {
        for (int i = 5; i < 6; ++i)
        {
            particlesIn.emplace_back(
                Particle{1., 1., ConstArray<int, dim>(i), ConstArray<double, dim>(0.6), {10., 10., 10.}});
            particlesOut.emplace_back(
                Particle{1., 1., ConstArray<int, dim>(i), ConstArray<double, dim>(0.), {10., 10., 10.}});
        }

        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);

        if constexpr (dim == 1)
        {
            auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
            auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

            for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
            {
                Jx(ix) = 0.0;
            }
        }
    }

    ~EsirkepovTest() {}
};

using TupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>>;

TYPED_TEST_SUITE(EsirkepovTest, TupleInfos);



TYPED_TEST(EsirkepovTest, EsirkepovCalculatedOk)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
    auto layout   = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

    auto rangeIn  = makeIndexRange(this->particlesIn);
    auto rangeOut = makeIndexRange(this->particlesOut);
    //std::copy(rangeIn.begin(), rangeIn.end(), rangeOut.begin());

    this->projector.setLayout(layout.get());
    //this->projector(this->J, rangeOut, rangeIn, this->layout, 1.);

    //std::vector<double> dJx;
    
    if constexpr (dim == 1)
    {
        auto psi_X = this->layout.physicalStartIndex(this->Jx, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(this->Jx, Direction::X);
        
        std::vector<double> weight_chargeOld(pei_X, 0.);
        std::vector<double> weight_chargeNew(pei_X, 0.);
        std::vector<double> dRho(pei_X, 0.);
        std::vector<double> dJx(pei_X, 0.);

        std::cout << "Projector" << std::endl;

        //rangeOut = this->pusher->move(rangeIn, rangeOut, this->em, this->mass, this->interpolator, this->layout, this->selector);
        this->projector(this->J, rangeOut, rangeIn, this->layout, 1.);

        std::cout << "Post projector" << std::endl;

        auto& particlesIn = rangeIn.array();
        auto& particlesOut = rangeOut.array();

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            int cell_idx = ix - 2 ;
            std::cout << "cell_idx = " << cell_idx << std::endl;

            dJx[ix] = ( (this->Jx(ix+1) - this->Jx(ix)) / this->layout.meshSize()[0]);

            for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
                        ++inIdx, ++outIdx)
                {
                    auto& currPart = particlesIn[inIdx] ;
                    auto& currPartOut = particlesOut[outIdx] ;


                    if (currPart.iCell[0] == ix -2)
                    {
                        if (currPart.delta[0]<0.5)
                        {
                            cell_idx = ix - 3 ;
                        }
                        else
                        {cell_idx = ix - 2 ;
                        }

                        
                        std::cout << "currPart.iCell[0] = " << currPart.iCell[0] << std::endl;
                        auto distance_to_dual1 =  0.5 + currPart.delta[0] ;
                        auto distance_to_dual2 = 1. - distance_to_dual1 ;
                        //std::cout << "deltaIn = " << currPart.delta[0] << std::endl;
                        //std::cout << "deltaOut = " << currPartOut.delta[0] << std::endl;

                        auto distance_to_dual1_out =  0.5 + currPartOut.delta[0] ;
                        auto distance_to_dual2_out = 1. - distance_to_dual1_out ;
                        //std::cout << "distance_to_dual1 = " << distance_to_dual1 << std::endl;
                        //std::cout << "distance_to_dual1_out = " << distance_to_dual1_out << std::endl;

                        if (cell_idx >= 0)
                        {
                            weight_chargeOld[cell_idx] += currPart.charge * currPart.weight * distance_to_dual1 ;
                            weight_chargeNew[cell_idx] += currPart.charge * currPart.weight * distance_to_dual1_out ;
                            std::cout << "weight_chargeOld[" << cell_idx << "] = " << weight_chargeOld[cell_idx] << std::endl;
                            std::cout << "weight_chargeNew[" << cell_idx << "] = " << weight_chargeNew[cell_idx] << std::endl;
                        }

                        weight_chargeOld[cell_idx + 1] += currPart.charge * currPart.weight * distance_to_dual2 ;
                        weight_chargeNew[cell_idx + 1] += currPart.charge * currPart.weight * distance_to_dual2_out ;
                        std::cout << "weight_chargeOld[" << cell_idx+1 << "] = " << weight_chargeOld[cell_idx+1] << std::endl;
                        std::cout << "weight_chargeNew[" << cell_idx+1 << "] = " << weight_chargeNew[cell_idx+1] << std::endl;

                    }
                }
            
            dRho[cell_idx] = (weight_chargeOld[cell_idx] - weight_chargeNew[cell_idx])/this->dt ;
            if (cell_idx < pei_X)
                dRho[cell_idx+1] = (weight_chargeOld[cell_idx+1] - weight_chargeNew[cell_idx+1])/this->dt ;
        }
        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            std::cout << "dRho[" << ix-2 << "] = " << dRho[ix-2] << std::endl;
            std::cout << "dJx[" << ix << "] = " << dJx[ix] << std::endl;
            std::cout << "Jx[" << ix << "] = " << this->Jx(ix) << std::endl;
        }
        std::cout << "psi_X = " << psi_X << std::endl;
        std::cout << "pei_X = " << pei_X << std::endl;

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            auto cell_idx = ix - 2 ;
            std::cout << "TEST\n" << std::endl;
            std::cout << "Comparing dJx[" << ix << "] = " << -dJx[ix] << " with dRho[" << cell_idx << "] = " << dRho[cell_idx] << std::endl;
            EXPECT_NEAR(-dJx[ix], dRho[cell_idx], 1e-12);
        }
    }

    


}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
