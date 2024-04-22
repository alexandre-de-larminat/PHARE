#ifndef PHARE_TAGGER_FACTORY_HPP
#define PHARE_TAGGER_FACTORY_HPP

#include <string>
#include <memory>

#include "tagger.hpp"
#include "hybrid_tagger.hpp"
#include "hybrid_tagger_strategy.hpp"
#include "default_hybrid_tagger_strategy.hpp"
#include "pic_tagger.hpp"
#include "pic_tagger_strategy.hpp"
#include "default_pic_tagger_strategy.hpp"
#include "core/def.hpp"

namespace PHARE::amr
{
template<typename PHARE_T>
class TaggerFactory
{
public:
    TaggerFactory() = delete;
    NO_DISCARD static std::unique_ptr<Tagger> make(std::string modelName, std::string methodName);
};

template<typename PHARE_T>
std::unique_ptr<Tagger> TaggerFactory<PHARE_T>::make(std::string modelName, std::string methodName)
{
    if (modelName == "HybridModel")
    {
        using HybridModel = typename PHARE_T::HybridModel_t;
        using HT          = HybridTagger<HybridModel>;

        if (methodName == "default")
        {
            using HTS = DefaultHybridTaggerStrategy<HybridModel>;
            return std::make_unique<HT>(std::make_unique<HTS>());
        }
    }
    
    if (modelName == "PICModel")
    {
        using PICModel = typename PHARE_T::PICModel_t;
        using PT          = PICTagger<PICModel>;

        if (methodName == "default")
        {
            using PTS = DefaultPICTaggerStrategy<PICModel>;
            return std::make_unique<PT>(std::make_unique<PTS>());
        }
    }
    return nullptr;
}


} // namespace PHARE::amr




#endif
