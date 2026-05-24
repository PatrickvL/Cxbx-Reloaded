#ifndef CXBXR_MCPX_DSP_QEMU_BSWAP_H
#define CXBXR_MCPX_DSP_QEMU_BSWAP_H

#include <stdint.h>
#include <string.h>

static inline uint32_t ldl_le_p(const void* ptr)
{
    uint32_t value = 0;
    memcpy(&value, ptr, sizeof(value));
    return value;
}

static inline void stl_le_p(void* ptr, uint32_t value)
{
    memcpy(ptr, &value, sizeof(value));
}

#endif
