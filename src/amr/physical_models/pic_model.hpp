#ifndef PHARE_PIC_MODEL_HPP
#define PHARE_PIC_MODEL_HPP

#include <string>
#include <cmath>

#include "initializer/data_provider.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/models/pic_state.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"
#include "amr/resources_manager/resources_manager.hpp"
#include "amr/messengers/pic_messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/def.hpp"

#include "core/data/fermions/fermions.hpp"

namespace PHARE::solver
{
/**
 * @brief The PICModel class is a concrete implementation of a IPhysicalModel. The class
 * holds a PICState and a ResourcesManager.
 */
template<typename GridLayoutT, typename Electromag, typename Fermions, typename AMR_Types>
class PICModel : public IPhysicalModel<AMR_Types>
{
public:

    using type_list = PHARE::core::type_list<GridLayoutT, Electromag, Fermions, AMR_Types>;
    using Interface = IPhysicalModel<AMR_Types>;
    using amr_types = AMR_Types;
    using patch_t   = typename AMR_Types::patch_t;
    using level_t   = typename AMR_Types::level_t;

    static const std::string model_name;
    using gridlayout_type              = GridLayoutT;
    using electromag_type              = Electromag;
    using vecfield_type                = typename Electromag::vecfield_type;
    using field_type                   = typename vecfield_type::field_type;
    using ions_type                    = typename Fermions::ions_type;
    using fluid_electrons_type         = typename Fermions::fluid_electrons_type;
    using fermions_type                = Fermions;
    using particle_array_type          = typename Fermions::particle_array_type;
    using hybrid_model_type            = HybridModel<GridLayoutT, Electromag, ions_type, fluid_electrons_type, AMR_Types>;
    using resources_manager_type       = amr::ResourcesManager<gridlayout_type>;
    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::interp_order;

    using ParticleInitializerFactory
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;

    core::PICState<Electromag, Fermions> state; 
    std::shared_ptr<resources_manager_type> resourcesManager;



    virtual void initialize(level_t& level) override;

    /**
     * @brief allocate uses the ResourcesManager to allocate PICState physical quantities on
     * the given Patch at the given allocateTime
     */
    virtual void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
    }

    // TODO check purpose
    auto patch_data_ids() const { return resourcesManager->restart_patch_data_ids(*this); }

    /**
     * @brief fillMessengerInfo describes which variables of the model are to be initialized or
     * filled at ghost nodes.
     */
    virtual void
    fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;

    NO_DISCARD auto setOnPatch(patch_t& patch)
    {
        return resourcesManager->setOnPatch(patch, *this);
    }

   
    explicit PICModel(PHARE::initializer::PHAREDict const& dict,
                      std::shared_ptr<resources_manager_type> const& _resourcesManager)
                    : IPhysicalModel<AMR_Types>{model_name}
                    , state{dict}
                    , resourcesManager{std::move(_resourcesManager)}
    {
    }

    virtual ~PICModel() override = default;

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesUserList() const { return std::forward_as_tuple(state); }

    NO_DISCARD auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(state); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;

};

template<typename GridLayoutT, typename Electromag, typename Fermions, typename AMR_Types>
const std::string PICModel<GridLayoutT, Electromag, Fermions, AMR_Types>::model_name = "PICModel";

//-------------------------------------------------------------------------
//                             definitions
//-------------------------------------------------------------------------


template<typename GridLayoutT, typename Electromag, typename Fermions,
         typename AMR_Types>
void PICModel<GridLayoutT, Electromag, Fermions, AMR_Types>::initialize(level_t& level)
{
    printf("PICModel::initialize\n");
    using InitFunctionT = PHARE::initializer::InitFunction<dimension>;
    double T = 0.0;

    for (auto& patch : level)
    {
        // first initialize the ions
        printf("Declaring layout\n");
        auto layout = amr::layoutFromPatch<gridlayout_type>(*patch);
        auto& ions = state.ions;
        auto& electrons = state.pic_electrons;
        auto _      = this->resourcesManager->setOnPatch(*patch, state.electromag, ions, electrons);

        for (auto& pop : ions)
        {
            printf("Initializing ions\n");
            auto const& info         = pop.particleInitializerInfo();
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);         
        }
        // then the electrons
        for (auto& pop : electrons)
        {
            printf("Initializing electrons\n");
            auto const& info         = pop.particleInitializerInfo(); 
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);
        }
        printf("Initializing electromag\n");
        state.electromag.initialize(layout);
    }
    printf("PICModel::initialize done, register for restarts\n");
    resourcesManager->registerForRestarts(*this);
}

// TODO adapt to PIC messenger
template<typename GridLayoutT, typename Electromag, typename Fermions, typename AMR_Types>
void PICModel<GridLayoutT, Electromag, Fermions, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    printf("model::fillMessengerInfo\n");
    auto& picInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);
    printf("THIS PRINTS??????\n");
    auto& ions = state.ions;

    picInfo.modelMagnetic        = core::VecFieldNames{state.electromag.B};
    picInfo.modelElectric        = core::VecFieldNames{state.electromag.E};
    picInfo.modelIonDensity      = ions.densityName();
    picInfo.modelIonBulkVelocity = core::VecFieldNames{ions.velocity()};
    picInfo.modelCurrent         = core::VecFieldNames{state.J};

    picInfo.initElectric.emplace_back(core::VecFieldNames{state.electromag.E});
    picInfo.initMagnetic.emplace_back(core::VecFieldNames{state.electromag.B});

    picInfo.ghostElectric.push_back(picInfo.modelElectric);
    picInfo.ghostMagnetic.push_back(picInfo.modelMagnetic);
    picInfo.ghostCurrent.push_back(core::VecFieldNames{state.J});
    picInfo.ghostBulkVelocity.push_back(picInfo.modelIonBulkVelocity);


    auto transform_ = [](auto& ions, auto& inserter) {
        std::transform(std::begin(ions), std::end(ions), std::back_inserter(inserter),
                       [](auto const& pop) { return pop.name(); });
    };
    transform_(ions, picInfo.interiorParticles);
    transform_(ions, picInfo.levelGhostParticlesOld);
    transform_(ions, picInfo.levelGhostParticlesNew);
    transform_(ions, picInfo.patchGhostParticles);
}


/* 
Used to define certain behaviors at compile. Used in diagnotics and restart. Useless unless those
files are modified, which will happen at some point.
*/
template<typename... Args>
PICModel<Args...> pic_model_from_type_list(core::type_list<Args...>);

template<typename TypeList>
struct type_list_to_pic_model
{
    using type = decltype(pic_model_from_type_list(std::declval<TypeList>()));
};

template<typename TypeList>
using type_list_to_pic_model_t = typename type_list_to_pic_model<TypeList>::type;




} // namespace PHARE::solver

#endif
