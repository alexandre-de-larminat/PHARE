#include <type_traits>



#include "core/data/fermions/electron_population.hpp"
#include "core/data/fermions/pic_electrons.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "initializer/data_provider.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

using namespace PHARE::core;

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using MaxwellianParticleInitializer1D = MaxwellianParticleInitializer<ParticleArray<1>, GridYee1D>;



class thePICElectrons : public ::testing::Test
{
protected:
    using VecField1D       = VecField<NdArrayVector<1>, HybridQuantity>;
    using SymTensorField1D = SymTensorField<NdArrayVector<1>, HybridQuantity>;
    using InitFunctionT    = PHARE::initializer::InitFunction<1>;

    using IonPopulation1D
        = IonPopulation<ParticleArray<1>, VecField1D, SymTensorField1D, GridYee1D>;
    PICElectrons<IonPopulation1D, GridYee1D> pic_electrons;

    PHARE::initializer::PHAREDict createPICElectronsDict()
    {
        PHARE::initializer::PHAREDict dict;
        dict["pic_electrons"]["nbrPopulations"] = std::size_t{2};
        dict["pic_electrons"]["pop0"]["name"]   = std::string{"protons"};
        dict["pic_electrons"]["pop0"]["mass"]   = 1.;
        dict["pic_electrons"]["pop0"]["particle_initializer"]["name"]
            = std::string{"MaxwellianParticleInitializer"};
        dict["pic_electrons"]["pop0"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT>(density);

        dict["pic_electrons"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT>(vx);

        dict["pic_electrons"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT>(vy);

        dict["pic_electrons"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT>(vz);


        dict["pic_electrons"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT>(vthx);

        dict["pic_electrons"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT>(vthy);

        dict["pic_electrons"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT>(vthz);


        dict["pic_electrons"]["pop0"]["particle_initializer"]["nbrPartPerCell"] = int{100};
        dict["pic_electrons"]["pop0"]["particle_initializer"]["charge"]         = -1.;
        dict["pic_electrons"]["pop0"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

        dict["pic_electrons"]["pop1"]["name"] = std::string{"alpha"};
        dict["pic_electrons"]["pop1"]["mass"] = 1.;
        dict["pic_electrons"]["pop1"]["particle_initializer"]["name"]
            = std::string{"MaxwellianParticleInitializer"};
        dict["pic_electrons"]["pop1"]["particle_initializer"]["density"]
            = static_cast<InitFunctionT>(density);

        dict["pic_electrons"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
            = static_cast<InitFunctionT>(vx);

        dict["pic_electrons"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
            = static_cast<InitFunctionT>(vy);

        dict["pic_electrons"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
            = static_cast<InitFunctionT>(vz);


        dict["pic_electrons"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
            = static_cast<InitFunctionT>(vthx);

        dict["pic_electrons"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
            = static_cast<InitFunctionT>(vthy);

        dict["pic_electrons"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
            = static_cast<InitFunctionT>(vthz);


        dict["pic_electrons"]["pop1"]["particle_initializer"]["nbrPartPerCell"] = int{100};
        dict["pic_electrons"]["pop1"]["particle_initializer"]["charge"]         = -1.;
        dict["pic_electrons"]["pop1"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

        return dict;
    }


    thePICElectrons()
        : pic_electrons{createPICElectronsDict()["pic_electrons"]}
    {
    }

public:
    ~thePICElectrons();
};

thePICElectrons::~thePICElectrons() {}


TEST_F(thePICElectrons, areAContainerOfIonPopulations)
{
    //
    for (auto& pop : pic_electrons)
    {
        (void)pop;
    }
}




TEST_F(thePICElectrons, areNotUsableUponConstruction)
{
    EXPECT_FALSE(pic_electrons.isUsable());
}




TEST_F(thePICElectrons, areSettableUponConstruction)
{
    EXPECT_TRUE(pic_electrons.isSettable());
}




TEST_F(thePICElectrons, throwIfAccessingDensityWhileNotUsable)
{
    EXPECT_ANY_THROW(auto& n = pic_electrons.density());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
