#ifndef PIC_TAGGER_STRATEGY_HPP
#define PIC_TAGGER_STRATEGY_HPP

namespace PHARE::amr
{

template<typename PICModel>
class PICTaggerStrategy
{
    using gridlayout_type = typename PICModel::gridlayout_type;

public:
    virtual void tag(PICModel& model, gridlayout_type const& layout, int* tags) const = 0;
    virtual ~PICTaggerStrategy()                                                      = 0;
};

template<typename PICModel>
PICTaggerStrategy<PICModel>::~PICTaggerStrategy()
{
}
}

#endif // HYBRID_TAGGER_STRATEGY_HPP
