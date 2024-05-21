#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <fstream>
#include <memory>


#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/ampere/maxwell_ampere.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"

#include "tests/core/data/field/test_field.hpp"
#include "tests/core/data/vecfield/test_vecfield.hpp"
#include "tests/core/data/gridlayout/gridlayout_test.hpp"

using namespace PHARE::core;




struct GridLayoutMock1D
{
    static const auto dimension = 1u;

    template<auto direction>
    double deriv([[maybe_unused]] FieldMock<1> const& f, [[maybe_unused]] MeshIndex<1u> mi)
    {
        return 0;
    }

    std::size_t physicalStartIndex([[maybe_unused]] FieldMock<1>&, [[maybe_unused]] Direction dir)
    {
        return 0;
    }
    std::size_t physicalEndIndex([[maybe_unused]] FieldMock<1>&, [[maybe_unused]] Direction dir)
    {
        return 0;
    }
};

struct GridLayoutMock2D
{
    static const auto dimension = 2u;

    template<auto direction>
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<2u> mi)
    {
        return 0;
    }

    std::size_t physicalStartIndex([[maybe_unused]] FieldMock<dimension>&,
                                   [[maybe_unused]] Direction dir)
    {
        return 0;
    }
    std::size_t physicalEndIndex([[maybe_unused]] FieldMock<dimension>&,
                                 [[maybe_unused]] Direction dir)
    {
        return 0;
    }
};

struct GridLayoutMock3D
{
    static const auto dimension = 3u;


    template<auto direction>
    double deriv([[maybe_unused]] FieldMock<dimension> const& f, [[maybe_unused]] MeshIndex<3u> mi)
    {
        return 0;
    }

    std::size_t physicalStartIndex([[maybe_unused]] FieldMock<dimension>&,
                                   [[maybe_unused]] Direction dir)
    {
        return 0;
    }
    std::size_t physicalEndIndex([[maybe_unused]] FieldMock<dimension>&,
                                 [[maybe_unused]] Direction dir)
    {
        return 0;
    }
};


PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["normalized_c"] = 1.0;

    return dict;
}


TEST(MaxwellAmpere, canBe1D)
{
    MaxwellAmpere<GridLayoutMock1D> maxwellAmpere{createDict()};
}


TEST(MaxwellAmpere, canBe2D)
{
    MaxwellAmpere<GridLayoutMock2D> maxwellAmpere{createDict()};
}


TEST(MaxwellAmpere, canBe3D)
{
    MaxwellAmpere<GridLayoutMock3D> maxwellAmpere{createDict()};
}


TEST(MaxwellAmpere, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    {
        using GridLayout = GridLayout<GridLayoutImplYee<1, 1>>;
        VecFieldMock<FieldMock<1>> B_1, E_1, J_1, Enew_1;
        MaxwellAmpere<GridLayout> MaxwellAmpere1d{createDict()};
        auto layout1d = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(MaxwellAmpere1d(B_1, E_1, J_1, Enew_1, 1.));
        MaxwellAmpere1d.setLayout(layout1d.get());
    }
/*
    {
        using GridLayout = GridLayout<GridLayoutImplYee<2, 1>>;
        VecFieldMock<FieldMock<2>> B_2, E_2, Bnew_2;
        MaxwellAmpere<GridLayout> MaxwellAmpere2d;
        auto layout2d = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(MaxwellAmpere2d(B_2, E_2, Bnew_2, 1.));
        MaxwellAmpere2d.setLayout(layout2d.get());
    }

    {
        using GridLayout = GridLayout<GridLayoutImplYee<3, 1>>;
        VecFieldMock<FieldMock<3>> B_3, E_3, Bnew_3;
        MaxwellAmpere<GridLayout> MaxwellAmpere3d;
        auto layout3d = std::make_unique<TestGridLayout<GridLayout>>();
        EXPECT_ANY_THROW(MaxwellAmpere3d(B_3, E_3, Bnew_3, 1.));
        MaxwellAmpere3d.setLayout(layout3d.get());
    }*/
}




std::vector<double> read(std::string filename)
{
    std::ifstream readFile(filename);
    assert(readFile.is_open());
    std::vector<double> x;

    std::copy(std::istream_iterator<double>(readFile), std::istream_iterator<double>(),
              std::back_inserter(x));
    return x;
}



class MaxwellAmpere1DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<1, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jz;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Exnew;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Eynew;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Eznew;
    VecField<NdArrayVector<1>, HybridQuantity> B;
    VecField<NdArrayVector<1>, HybridQuantity> E;
    VecField<NdArrayVector<1>, HybridQuantity> J;
    VecField<NdArrayVector<1>, HybridQuantity> Enew;
    MaxwellAmpere<GridLayout<GridLayoutImpl>> maxwellAmpere;

public:
    MaxwellAmpere1DTest()
        : layout{{{0.1}}, {{50}}, Point{0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Ey{"Ey", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Ez{"Ez", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , Exnew{"Exnew", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Eynew{"Eynew", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Eznew{"Eznew", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , J{"J", HybridQuantity::Vector::J}
        , B{"B", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
        , Enew{"Enew", HybridQuantity::Vector::E}
        , maxwellAmpere{createDict()}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
        J.setBuffer("J_x", &Jx);
        J.setBuffer("J_y", &Jy);
        J.setBuffer("J_z", &Jz);
        Enew.setBuffer("Enew_x", &Exnew);
        Enew.setBuffer("Enew_y", &Eynew);
        Enew.setBuffer("Enew_z", &Eznew);
    }
};


/*

class MaxwellAmpere2DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<2, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Exnew;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Eynew;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Eznew;
    VecField<NdArrayVector<2>, HybridQuantity> B;
    VecField<NdArrayVector<2>, HybridQuantity> E;
    VecField<NdArrayVector<2>, HybridQuantity> Enew;
    MaxwellAmpere<GridLayout<GridLayoutImpl>> MaxwellAmpere;

public:
    MaxwellAmpere2DTest()
        : layout{{{0.1, 0.2}}, {{50, 30}}, Point{0., 0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Ey{"Ey", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Ez{"Ez", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Exnew{"Exnew", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , Eynew{"Bynew", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::By)}
        , Eznew{"Bznew", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , B{"B", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
        , Enew{"Enew", HybridQuantity::Vector::E}
        , J{"J", HybridQuantity::Vector::J}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
        Enew.setBuffer("Bnew_x", &Exnew);
        Enew.setBuffer("Bnew_y", &Eynew);
        Enew.setBuffer("Bnew_z", &Eznew);
    }
};




class MaxwellAmpere3DTest : public ::testing::Test
{
protected:
    using GridLayoutImpl = GridLayoutImplYee<3, 1>;
    GridLayout<GridLayoutImpl> layout;
    static constexpr auto interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bxnew;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bynew;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Bznew;
    VecField<NdArrayVector<3>, HybridQuantity> B;
    VecField<NdArrayVector<3>, HybridQuantity> E;
    VecField<NdArrayVector<3>, HybridQuantity> Bnew;
    MaxwellAmpere<GridLayout<GridLayoutImpl>> MaxwellAmpere;

public:
    MaxwellAmpere3DTest()
        : layout{{{0.1, 0.2, 0.3}}, {{50, 30, 40}}, Point{0., 0., 0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Ey{"Ey", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Ez{"Ez", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
        , Bxnew{"Bxnew", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , Bynew{"Bynew", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bznew{"Bznew", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , B{"B", HybridQuantity::Vector::B}
        , E{"E", HybridQuantity::Vector::E}
        , Bnew{"Bnew", HybridQuantity::Vector::B}
    {
        B.setBuffer("B_x", &Bx);
        B.setBuffer("B_y", &By);
        B.setBuffer("B_z", &Bz);
        E.setBuffer("E_x", &Ex);
        E.setBuffer("E_y", &Ey);
        E.setBuffer("E_z", &Ez);
        Bnew.setBuffer("Bnew_x", &Bxnew);
        Bnew.setBuffer("Bnew_y", &Bynew);
        Bnew.setBuffer("Bnew_z", &Bznew);
    }
};


*/

TEST_F(MaxwellAmpere1DTest, MaxwellAmpere1DCalculatedOk)
{
    auto filename_dbxdt = std::string{"dbxdt_yee_1D_order1.txt"};
    auto filename_dbydt = std::string{"dbydt_yee_1D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_1D_order1.txt"};
    auto expected_dbxdt = read(filename_dbxdt);
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(Ey, Point{0.}, ix);

        Ey(ix) = std::cos(2 * M_PI / 5. * point[0]);
        Ez(ix) = std::sin(2 * M_PI / 5. * point[0]);
    }

    auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        auto point = this->layout.fieldNodeCoordinates(By, Point{0.}, ix);

        By(ix) = std::tanh(point[0] - 5. / 2.);
        Bz(ix) = std::tanh(point[0] - 5. / 2.);
    }

    maxwellAmpere.setLayout(&layout);
    maxwellAmpere(B, E, J, Enew, 1.);

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        EXPECT_THAT(Exnew(ix), ::testing::DoubleNear((expected_dbxdt[ix]), 1e-12));
        EXPECT_THAT(Eynew(ix), ::testing::DoubleNear((expected_dbydt[ix]), 1e-12));
        EXPECT_THAT(Eznew(ix), ::testing::DoubleNear((expected_dbzdt[ix]), 1e-12));
    }
}


/*

TEST_F(MaxwellAmpere2DTest, MaxwellAmpere2DCalculatedOk)
{
    auto filename_dbxdt = std::string{"dbxdt_yee_2D_order1.txt"};
    auto filename_dbydt = std::string{"dbydt_yee_2D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_2D_order1.txt"};
    auto expected_dbxdt = read(filename_dbxdt);
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
    auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Ex, Point{0., 0.}, ix, iy);

            Ex(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Ey, Point{0., 0.}, ix, iy);

            Ey(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Ez, Point{0., 0.}, ix, iy);

            Ez(ix, iy) = std::sin(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Bx, Point{0., 0.}, ix, iy);

            Bx(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(By, Point{0., 0.}, ix, iy);

            By(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(Bz, Point{0., 0.}, ix, iy);

            Bz(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
        }
    }

    MaxwellAmpere.setLayout(&layout);
    MaxwellAmpere(B, E, Bnew, 1.);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    auto nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bx);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Bxnew(ix, iy), ::testing::DoubleNear((expected_dbxdt[index_]), 1e-12));
        }
    }

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::By);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Bynew(ix, iy), ::testing::DoubleNear((expected_dbydt[index_]), 1e-12));
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bz);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(Bznew(ix, iy), ::testing::DoubleNear((expected_dbzdt[index_]), 1e-12));
        }
    }
}



TEST_F(MaxwellAmpere3DTest, MaxwellAmpere3DCalculatedOk)
{
    auto filename_dbxdt = std::string{"dbxdt_yee_3D_order1.txt"};
    auto filename_dbydt = std::string{"dbydt_yee_3D_order1.txt"};
    auto filename_dbzdt = std::string{"dbzdt_yee_3D_order1.txt"};
    auto expected_dbxdt = read(filename_dbxdt);
    auto expected_dbydt = read(filename_dbydt);
    auto expected_dbzdt = read(filename_dbzdt);

    auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto gsi_p_Y = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Y);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gsi_d_Y = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Y);
    auto gei_d_Y = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Y);

    auto gsi_p_Z = this->layout.ghostStartIndex(QtyCentering::primal, Direction::Z);
    auto gei_p_Z = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Z);
    auto gsi_d_Z = this->layout.ghostStartIndex(QtyCentering::dual, Direction::Z);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Ex, Point{0., 0., 0.}, ix, iy, iz);

                Ex(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                 * std::cos(2 * M_PI / 6. * point[1])
                                 * std::tanh(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Ey, Point{0., 0., 0.}, ix, iy, iz);

                Ey(ix, iy, iz) = std::tanh(2 * M_PI / 5. * point[0])
                                 * std::sin(2 * M_PI / 6. * point[1])
                                 * std::cos(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Ez, Point{0., 0., 0.}, ix, iy, iz);

                Ez(ix, iy, iz) = std::cos(2 * M_PI / 5. * point[0])
                                 * std::tanh(2 * M_PI / 6. * point[1])
                                 * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Bx, Point{0., 0., 0.}, ix, iy, iz);

                Bx(ix, iy, iz) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.)
                                 * std::tanh(point[2] - 12. / 2.);
            }
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = gsi_d_Z; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(By, Point{0., 0., 0.}, ix, iy, iz);

                By(ix, iy, iz) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.)
                                 * std::tanh(point[2] - 12. / 2.);
            }
        }
    }

    for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
    {
        for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
        {
            for (auto iz = gsi_p_Z; iz <= gei_p_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(Bz, Point{0., 0., 0.}, ix, iy, iz);

                Bz(ix, iy, iz) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.)
                                 * std::tanh(point[2] - 12. / 2.);
            }
        }
    }

    MaxwellAmpere.setLayout(&layout);
    MaxwellAmpere(B, E, Bnew, 1.);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    auto psi_p_Z = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
    auto pei_p_Z = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Z);
    auto psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    auto pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    auto nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bx);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Bxnew(ix, iy, iz),
                            ::testing::DoubleNear((expected_dbxdt[index_]), 1e-12));
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::By);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Bynew(ix, iy, iz),
                            ::testing::DoubleNear((expected_dbydt[index_]), 1e-12));
            }
        }
    }

    nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Bz);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (auto iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(Bznew(ix, iy, iz),
                            ::testing::DoubleNear((expected_dbzdt[index_]), 1e-12));
            }
        }
    }
}*/




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
