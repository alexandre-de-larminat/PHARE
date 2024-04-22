#ifndef PHARE_PIC_MESSENGER_INFO_HPP
#define PHARE_PIC_MESSENGER_INFO_HPP

#include "messenger_info.hpp"



namespace PHARE::amr
{
    class PICMessengerInfo : public IMessengerInfo
    {
    public:
        virtual ~PICMessengerInfo() = default;
    };

} // namespace PHARE::amr
#endif
