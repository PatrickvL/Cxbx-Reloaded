#ifndef CXBXR_MCPX_DSP_DEBUG_H
#define CXBXR_MCPX_DSP_DEBUG_H

#ifndef DEBUG_DSP
#define DEBUG_DSP 0
#endif

#define TRACE_DSP_DISASM 0
#define TRACE_DSP_DISASM_REG 0
#define TRACE_DSP_DISASM_MEM 0

#define DPRINTF(fmt, ...) \
    do { \
        if (DEBUG_DSP) fprintf(stderr, fmt, ##__VA_ARGS__); \
    } while (0)

#endif
