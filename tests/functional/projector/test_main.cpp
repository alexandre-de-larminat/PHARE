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

#include "tests/core/data/field/test_field.hpp"
#include "tests/core/data/vecfield/test_vecfield.hpp"
#include "tests/core/data/gridlayout/gridlayout_test.hpp"

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
            return {{{0.1}}, {{50}}, Point{0.}};
        }
        else if constexpr (dim == 2)
        {
            return {{{0.1, 0.2}}, {{50, 30}}, Point{0., 0.}};
        }
        else if constexpr (dim == 3)
        {
            return {{{0.1, 0.2, 0.3}}, {{50, 30, 40}}, Point{0., 0., 0.}};
        }
    }
};





template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct EsirkepovTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee  = GridLayout<GridLayoutImplYee<dim, interp>>;
    GridYee layout = NDlayout<dim, interp>::create();
    
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jz;
    VecField<NdArrayVector<dim>, HybridQuantity> J;
    Projector<GridYee> projector;

public:
    EsirkepovTest()
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , J{"J", HybridQuantity::Vector::J}
        , particlesIn{layout.AMRBox()}
        , particlesOut(layout.AMRBox())

    {
        particlesIn.emplace_back(
            Particle{1., 1., ConstArray<int, dim>(5), ConstArray<double, dim>(0.), {0., 10., 0.}});
        particlesOut.emplace_back(
            Particle{1., 1., ConstArray<int, dim>(5), ConstArray<double, dim>(0.), {0., 10., 0.}});
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
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<3>, InterpConst<1>>>;

TYPED_TEST_SUITE(EsirkepovTest, TupleInfos);




TYPED_TEST(EsirkepovTest, ShouldBeGivenAGridLayoutPointerToBeOperational)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    auto rangeIn  = makeIndexRange(this->particlesIn);
    auto rangeOut = makeIndexRange(this->particlesOut);
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut.begin());


    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    auto layout = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

    // this->maxwellAmpere.setLayout(layout.get());
    EXPECT_ANY_THROW(this->projector(this->J, rangeOut, rangeIn, layout, 1.)); // because the grid layout is not set (TODO issue #3392)
}





TYPED_TEST(EsirkepovTest, EsirkepovCalculatedOk)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
    auto layout   = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

    auto rangeIn  = makeIndexRange(this->particlesIn);
    auto rangeOut = makeIndexRange(this->particlesOut);
    std::copy(rangeIn.begin(), rangeIn.end(), rangeOut.begin());

    this->projector.setLayout(layout.get());
    this->projector(this->J, rangeOut, rangeIn, layout, 1.);
    
    if constexpr (dim == 1)
    {
        auto psi_X = this->layout.physicalStartIndex(this->Jx, Direction::X);
        auto pei_X = this->layout.physicalEndIndex(this->Jx, Direction::X);

        for (auto ix = psi_X; ix <= pei_X; ++ix)
        {
            dJx = (Jx(ix + 1) - Jx(ix)) / this->layout.meshSize()[0];
            EXPECT_THAT(this->Exnew(ix), ::testing::DoubleNear((expected_MAx[ix]), 1e-12));
        }

    }

    


}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
