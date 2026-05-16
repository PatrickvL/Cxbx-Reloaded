#ifndef NV2A_PGRAPH_BACKEND_H
#define NV2A_PGRAPH_BACKEND_H

struct NV2AState;

// Backend interface for PGRAPH rendering.
// The PGRAPH method handler dispatches draw/state operations through these
// function pointers, allowing different rendering backends (D3D11, Vulkan, etc.)
// to be plugged in without changing PGRAPH logic.
struct PgraphBackend {
    void (*draw)(NV2AState *d);
    void (*draw_state_update)(NV2AState *d);
    void (*draw_clear)(NV2AState *d);
    void (*draw_patch)(NV2AState *d);
    void (*flip_stall)(NV2AState *d);
    void (*zpass_begin)(NV2AState *d);
    void (*zpass_end)(NV2AState *d);
    void (*zpass_collect)(NV2AState *d);
    void (*launch_transform_program)(NV2AState *d, unsigned int program_start);
};

// Global backend instance — set once at init by the active renderer.
extern PgraphBackend g_pgraph_backend;

#endif // NV2A_PGRAPH_BACKEND_H
