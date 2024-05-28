#ifndef PHARE_LEVEL_INITIALIZER_FACTORY_HPP
#define PHARE_LEVEL_INITIALIZER_FACTORY_HPP

#include "hybrid_level_initializer.hpp"
#include "pic_level_initializer.hpp"
#include "level_initializer.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"

#include <memory>
#include <string>

namespace PHARE
{
namespace solver
{
    template<typename HybridModel, typename PICModel>
    class LevelInitializerFactory
    {
        using AMRTypes = typename PICModel::amr_types;//EDITED

    public:
        NO_DISCARD static std::unique_ptr<LevelInitializer<AMRTypes>>
        create(std::string modelName, PHARE::initializer::PHAREDict const& dict)
        {
            if (modelName == "HybridModel")
            {
                return std::make_unique<HybridLevelInitializer<HybridModel>>(dict);
            }
            if (modelName == "PICModel")
            {
                return std::make_unique<PICLevelInitializer<PICModel>>(dict);
            }
            return nullptr;
        }
    };

} // namespace solver
} // namespace PHARE



#endif
