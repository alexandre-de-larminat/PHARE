#ifndef PHARE_PIC_MESSENGER_INFO_HPP
#define PHARE_PIC_MESSENGER_INFO_HPP

#include "messenger_info.hpp"



namespace PHARE::amr
{
    class PICMessengerInfo : public HybridMessengerInfo
    {

    public:

        core::VecFieldNames modelElectronBulkVelocity;
        std::string modelElectronDensity;

        virtual ~PICMessengerInfo() = default;
    };

} // namespace PHARE::amr
#endif
