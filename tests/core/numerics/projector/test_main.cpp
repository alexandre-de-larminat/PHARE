#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <list>
#include <random>

#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "phare_core.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/numerics/projector/projector.hpp"


using namespace PHARE::core;

template<std::size_t dimension_, std::size_t interp_order_>
struct DummyLayout
{
    static constexpr std::size_t dimension = dimension_;
    static constexpr auto interp_order = interp_order_;
    std::array<unsigned int, dimension> nbrCells_;
    auto nbrCells() const { return nbrCells_; }
    auto AMRBox() const { return PHARE::core::emptyBox<int, dimension>(); }
};

struct DummyField
{
};

template<std::size_t dimension_>
struct DummyVecField
{
    static constexpr std::size_t dimension = dimension_;
    using field_type                       = DummyField;
    DummyVecField(std::string name, HybridQuantity::Vector /*v*/) { (void)name; }
    bool isUsable() const { return false; }
    bool isSettable() const { return true; }
};



template<std::size_t dimension_, std::size_t interp_order_>
auto E_weights_computer()
{
    static const int dim = dimension_;
    static const int interp_order = interp_order_;

    using Particle = typename ParticleArray<dim>::Particle_t;

    static const int nbr_tests = 1000;

    std::array<double, nbr_tests> normalizedPositions;
    std::vector<double> weightsSums(nbr_tests*dim);
    DummyLayout<dim, interp_order> layout;

    std::vector<std::vector<std::vector<double>>> weights_(nbr_tests, std::vector<std::vector<double>>(dim, std::vector<double>(interp_order + 3)));

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(5, 10);


    // generate the nbr_tests random (normalized) particle position
    std::generate(std::begin(normalizedPositions), std::end(normalizedPositions),
                      [&dis, &gen]() { return dis(gen); });
 
    // now for each random position, calculate
    // the start index and the N(interp_order) weights
    for (auto i = 0u; i < nbr_tests*dim ; i+=(1*dim))
    {
        auto smol_i = i/dim;
        auto icell = static_cast<int>(normalizedPositions[smol_i]);
        auto delta = normalizedPositions[smol_i] - icell;
        Particle particle{1., 1., ConstArray<int, dim>(icell), ConstArray<double, dim>(delta), {0., 10., 0.}};
        weights_[smol_i] = E_weights_(particle, layout);

        weightsSums[i] = std::accumulate(std::begin(weights_[smol_i][0]), std::end(weights_[smol_i][0]), 0.);
        if (dim > 1)
        {
            weightsSums[i+1] = std::accumulate(std::begin(weights_[smol_i][1]), std::end(weights_[smol_i][1]), 0.);
        }
        if (dim > 2)
        {
            weightsSums[i+2] = std::accumulate(std::begin(weights_[smol_i][2]), std::end(weights_[smol_i][2]), 0.);
        }   
    } 
  
    return weightsSums;
}


TEST(Weights_for_Esirkepov1D, ComputesWeightThatSumIsOne)
{
    auto weightsSums1 = E_weights_computer<1, 1>();
    auto equalsOne = [](double sum) { return std::abs(sum - 1.) < 1e-10; };
    EXPECT_TRUE(std::all_of(std::begin(weightsSums1), std::end(weightsSums1), equalsOne));

    auto weightsSums2 = E_weights_computer<1, 2>();
    EXPECT_TRUE(std::all_of(std::begin(weightsSums2), std::end(weightsSums2), equalsOne));

    auto weightsSums3 = E_weights_computer<1, 3>();
    EXPECT_TRUE(std::all_of(std::begin(weightsSums3), std::end(weightsSums3), equalsOne));
}

TEST(Weights_for_Esirkepov2D, ComputesWeightThatSumIsOne)
{
    auto weightsSums1 = E_weights_computer<2, 1>();
    auto equalsOne = [](double sum) { return std::abs(sum - 1.) < 1e-10; };
    EXPECT_TRUE(std::all_of(std::begin(weightsSums1), std::end(weightsSums1), equalsOne));

    auto weightsSums2 = E_weights_computer<2, 2>();
    EXPECT_TRUE(std::all_of(std::begin(weightsSums2), std::end(weightsSums2), equalsOne));

    auto weightsSums3 = E_weights_computer<2, 3>();
    EXPECT_TRUE(std::all_of(std::begin(weightsSums3), std::end(weightsSums3), equalsOne));
}

TEST(Weights_for_Esirkepov3D, ComputesWeightThatSumIsOne)
{
    auto weightsSums1 = E_weights_computer<3, 1>();
    auto equalsOne = [](double sum) { return std::abs(sum - 1.) < 1e-10; };
    EXPECT_TRUE(std::all_of(std::begin(weightsSums1), std::end(weightsSums1), equalsOne));

    auto weightsSums2 = E_weights_computer<3, 2>();
    EXPECT_TRUE(std::all_of(std::begin(weightsSums2), std::end(weightsSums2), equalsOne));

    auto weightsSums3 = E_weights_computer<3, 3>();
    EXPECT_TRUE(std::all_of(std::begin(weightsSums3), std::end(weightsSums3), equalsOne));
}



template<typename DummyLayout>
class AProjector : public ::testing::Test
{
public:
    static constexpr auto dim    = DummyLayout::dimension;
    static constexpr auto interp_order = DummyLayout::interp_order;
    AProjector()
    : particlesIn{layout.AMRBox()}
    , particlesOut(layout.AMRBox())
    {
        particlesIn.emplace_back(
            Particle{1., 1., ConstArray<int, dim>(5), ConstArray<double, dim>(0.), {0., 10., 0.}});
        particlesOut.emplace_back(
            Particle{1., 1., ConstArray<int, dim>(5), ConstArray<double, dim>(0.1), {0., 10., 0.}});

        this->projector.Projector(vecfield, particlesOut, particlesIn, layout, dt);
    }

protected:

    using Particle = typename ParticleArray<dim>::Particle_t;
    DummyLayout layout;
    ParticleArray<dim> particlesIn;
    ParticleArray<dim> particlesOut;
    DummyVecField<dim> vecfield;
    Projector<DummyLayout> projector;
    double dt{0.0001};
};

using AProjector1D = AProjector<DummyLayout<1,1>>;
using AProjector2D = AProjector<DummyLayout<2,1>>;
using AProjector3D = AProjector<DummyLayout<3,1>>;



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
