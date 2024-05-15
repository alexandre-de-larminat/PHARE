#include <type_traits>




#include "core/data/pic_electrons/electron_population.hpp"
#include "core/data/particles/particle_array.hpp"
#include "initializer/data_provider.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;
using namespace PHARE::initializer;



struct DummyField
{
};


struct DummyVecField
{
    static constexpr std::size_t dimension = 1;
    using field_type                       = DummyField;
    DummyVecField(std::string name, HybridQuantity::Vector /*v*/) { (void)name; }
    bool isUsable() const { return false; }
    bool isSettable() const { return true; }
};

struct DummyTensorField
{
    static constexpr std::size_t dimension = 1;
    DummyTensorField(std::string name, HybridQuantity::Tensor /*v*/) { (void)name; }
    bool isUsable() const { return false; }
    bool isSettable() const { return true; }
};

struct DummyParticleInitializer
{
};


struct DummyLayout
{
};

PHAREDict getDict()
{
    PHAREDict dict;
    dict["name"]                         = std::string{"electrons"};
    dict["mass"]                         = 1.;
    dict["particle_initializer"]["name"] = std::string{"DummyParticleInitializer"};
    return dict;
}

struct AnElectronPopulation : public ::testing::Test
{
    ElectronPopulation<ParticleArray<1>, DummyVecField, DummyTensorField, DummyLayout> electrons{
        getDict()};
    virtual ~AnElectronPopulation();
};

AnElectronPopulation::~AnElectronPopulation() {}



TEST_F(AnElectronPopulation, hasAMass)
{
    EXPECT_DOUBLE_EQ(1., electrons.mass());
}




TEST_F(AnElectronPopulation, hasAName)
{
    EXPECT_EQ("electrons", electrons.name());
}




TEST_F(AnElectronPopulation, isNonUsableUponConstruction)
{
    EXPECT_FALSE(electrons.isUsable());
}




TEST_F(AnElectronPopulation, isSettableIfNonUsable)
{
    if (!electrons.isUsable())
    {
        EXPECT_TRUE(electrons.isSettable());
    }
}




TEST_F(AnElectronPopulation, throwsIfOneWantsToAccessParticleBuffersWhileNotUsable)
{
    EXPECT_ANY_THROW(auto& p = electrons.domainParticles());
    EXPECT_ANY_THROW(auto& p = electrons.patchGhostParticles());
    EXPECT_ANY_THROW(auto& p = electrons.levelGhostParticles());
}




TEST_F(AnElectronPopulation, isResourceUserAndHasGetParticleArrayNamesOK)
{
    auto bufferNames = electrons.getParticleArrayNames();
    EXPECT_EQ(1, bufferNames.size());
    EXPECT_EQ(electrons.name(), bufferNames[0].name);
}



TEST_F(AnElectronPopulation, isResourceUserAndHasFieldNamesAndQuantitiesOK)
{
    auto fieldProperties = electrons.getFieldNamesAndQuantities();
    EXPECT_EQ(electrons.name() + std::string{"_rho"}, fieldProperties[0].name);
    EXPECT_EQ(HybridQuantity::Scalar::rhoE, fieldProperties[0].qty);
}



TEST_F(AnElectronPopulation, hasAVecFieldSubResource)
{
    [[maybe_unused]] DummyVecField const& vf
        = std::get<0>(electrons.getCompileTimeResourcesUserList());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
