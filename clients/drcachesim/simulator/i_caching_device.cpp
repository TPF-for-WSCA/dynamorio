#include "i_caching_device.h"

I_caching_device_t::access_update(int block_idx, int way)
{
    cache_->access_update(block_idx, way);
}