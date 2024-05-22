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
            return {{{0.1}}, {{50}}, {0.}};
        }
        else if constexpr (dim == 2)
        {
            return {{{0.1, 0.2}}, {{50, 40}}, {0., 0.}};
        }
        else if constexpr (dim == 3)
        {
            return {{{0.1, 0.2, 0.3}}, {{50, 40, 30}}, {0., 0., 0.}};
        }
    }
};


PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;

    dict["normalized_c"] = 1.0;

    return dict;
}


template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct MaxwellAmpereTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee  = GridLayout<GridLayoutImplYee<dim, interp>>;
    GridYee layout = NDlayout<dim, interp>::create();
    
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Bx;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Bz;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Ex;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Ey;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Ez;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Jz;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Exnew;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Eynew;
    Field<NdArrayVector<dim>, HybridQuantity::Scalar> Eznew;
    VecField<NdArrayVector<dim>, HybridQuantity> B;
    VecField<NdArrayVector<dim>, HybridQuantity> E;
    VecField<NdArrayVector<dim>, HybridQuantity> J;
    VecField<NdArrayVector<dim>, HybridQuantity> Enew;
    MaxwellAmpere<GridYee> maxwellAmpere;

public:
    MaxwellAmpereTest()
        : Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
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

    ~MaxwellAmpereTest() {}
};

using TupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<3>, InterpConst<1>>>;

TYPED_TEST_SUITE(MaxwellAmpereTest, TupleInfos);


TYPED_TEST(MaxwellAmpereTest, ThatOhmHasCtorWithDict)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    MaxwellAmpere<GridYee> maxwellAmpere(createDict());
}


TYPED_TEST(MaxwellAmpereTest, ShouldBeGivenAGridLayoutPointerToBeOperational)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    auto layout = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

    // this->maxwellAmpere.setLayout(layout.get());
    EXPECT_ANY_THROW(this->maxwellAmpere(this->B, this->E, this->J, 
                     this->Enew, 1.)); // because the grid layout is not set (TODO issue #3392)
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




TYPED_TEST(MaxwellAmpereTest, MaxwellAmpere1DCalculatedOk)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    if constexpr (dim == 1)
    {
        auto filename_MAx = std::string{"MAx_yee_1D_order1.txt"};
        auto filename_MAy = std::string{"MAy_yee_1D_order1.txt"};
        auto filename_MAz = std::string{"MAz_yee_1D_order1.txt"};
        auto expected_MAx = read(filename_MAx);
        auto expected_MAy = read(filename_MAy);
        auto expected_MAz = read(filename_MAz);

        using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
        auto layout   = std::make_unique<GridYee>(NDlayout<dim, interp>::create());


        auto gsi_p_X = this->layout.ghostStartIndex(QtyCentering::primal, Direction::X);
        auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            auto point = this->layout.fieldNodeCoordinates(this->Ey, Point{0.}, ix);

            this->Ey(ix) = std::cos(2 * M_PI / 5. * point[0]);
            this->Ez(ix) = std::sin(2 * M_PI / 5. * point[0]);
        }

        auto gsi_d_X = this->layout.ghostStartIndex(QtyCentering::dual, Direction::X);
        auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            auto point = this->layout.fieldNodeCoordinates(this->By, Point{0.}, ix);

            this->By(ix) = std::tanh(point[0] - 5. / 2.);
            this->Bz(ix) = std::tanh(point[0] - 5. / 2.);
        }

        this->maxwellAmpere.setLayout(layout.get());
        this->maxwellAmpere(this->B, this->E, this->J, this->Enew, 1.);

        auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
        auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

        for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
        {
            EXPECT_THAT(this->Exnew(ix), ::testing::DoubleNear((expected_MAx[ix]), 1e-12));
            EXPECT_THAT(this->Eynew(ix), ::testing::DoubleNear((expected_MAy[ix]), 1e-12));
            EXPECT_THAT(this->Eznew(ix), ::testing::DoubleNear((expected_MAz[ix]), 1e-12));
        }
    }

    if constexpr (dim == 2)
    {
        auto filename_MAx = std::string{"MAx_yee_2D_order1.txt"};
        auto filename_MAy = std::string{"MAy_yee_2D_order1.txt"};
        auto filename_MAz = std::string{"MAz_yee_2D_order1.txt"};
        auto expected_MAx = read(filename_MAx);
        auto expected_MAy = read(filename_MAy);
        auto expected_MAz = read(filename_MAz);

        using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
        auto layout   = std::make_unique<GridYee>(NDlayout<dim, interp>::create());

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
                auto point = this->layout.fieldNodeCoordinates(this->Ey, Point{0., 0.}, ix, iy);

                this->Ex(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
            }
        }

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Ey, Point{0., 0.}, ix, iy);

                this->Ey(ix, iy) = std::cos(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
            }
        }

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Ez, Point{0., 0.}, ix, iy);

                this->Ez(ix, iy) = std::sin(2 * M_PI / 5. * point[0]) * std::tanh(2 * M_PI / 6. * point[1]);
            }
        }

        for (auto ix = gsi_p_X; ix <= gei_p_X; ++ix)
        {
            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Bx, Point{0., 0.}, ix, iy);

                this->Bx(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
            }
        }

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            for (auto iy = gsi_p_Y; iy <= gei_p_Y; ++iy)
            {
                auto point = this->layout.fieldNodeCoordinates(this->By, Point{0., 0.}, ix, iy);

                this->By(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
            }
        }

        for (auto ix = gsi_d_X; ix <= gei_d_X; ++ix)
        {
            for (auto iy = gsi_d_Y; iy <= gei_d_Y; ++iy)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Bz, Point{0., 0.}, ix, iy);

                this->Bz(ix, iy) = std::tanh(point[0] - 5. / 2.) * std::tanh(point[1] - 6. / 2.);
            }
        }

        this->maxwellAmpere.setLayout(layout.get());
        this->maxwellAmpere(this->B, this->E, this->J, this->Enew, 1.);

        auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
        auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
        auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

        for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
        {
            for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
            {
                auto index_ = ix * this->layout.allocSize(HybridQuantity::Scalar::Ex)[1] + iy;
                EXPECT_THAT(this->Exnew(ix, iy), ::testing::DoubleNear((expected_MAx[index_]), 1e-12));
                EXPECT_THAT(this->Eynew(ix, iy), ::testing::DoubleNear((expected_MAy[index_]), 1e-12));
                EXPECT_THAT(this->Eznew(ix, iy), ::testing::DoubleNear((expected_MAz[index_]), 1e-12));
            }
        }
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
