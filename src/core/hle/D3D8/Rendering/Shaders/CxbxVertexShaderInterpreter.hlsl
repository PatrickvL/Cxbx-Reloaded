// CxbxVertexShaderInterpreter.hlsl
//
// DX11 / SM5 Xbox NV2A vertex shader interpreter ubershader.
//
// Instead of recompiling each Xbox vertex shader program into host HLSL,
// this single precompiled shader interprets the raw NV2A microcode at
// runtime.  All state comes from two StructuredBuffers:
//   - g_PGRegs (t12): shared PGRAPH register array (same buffer as PS)
//     provides CHEOPS_PROGRAM_START to locate the active program slot.
//   - g_XFPR (t5): NV2A XFPR (Transform Program RAM) — 136 × uint4.
//     On real hardware this is on-chip XF SRAM behind the RDI interface.
// Vertex constants (c0–c191) remain in the existing cbuffer b0.
//
// Architecture mirrors the register combiner interpreter:
//   - C++ uploads raw data to StructuredBuffers
//   - This shader loops over instruction slots, decodes fields, executes ops
//   - No CPU-side D3DCompile, no shader cache, no async compilation
//
// Reference: nv2a_vsh_cpu (import/nv2a_vsh_cpu/src/nv2a_vsh_emulator.c)

#include "CxbxVertexShaderCommon.hlsli"
#include "CxbxVertexFetch.hlsli"

// Xbox constant registers
#define X_D3DSCM_CORRECTION 96
#define X_D3DVS_CONSTREG_COUNT 192
uniform float4 C[X_D3DVS_CONSTREG_COUNT] : register(c0);

#include "CxbxScreenspaceTransform.hlsli"
#include "CxbxPGRAPHRegs.hlsli"
#include "CxbxVertexShaderInterpreterState.hlsli"
#include "CxbxNV2AMathHelpers.hlsli"
#include "CxbxNV2AVshOps.hlsli"

// ============================================================
// Per-invocation interpreter state.
//
// HLSL `static` globals are per-thread (each vertex gets its own copy,
// reset at shader entry). Keeping state here instead of in locals avoids
// passing 12 + 16 + 16 float4 arrays through every helper call, which
// is the single biggest driver of compile time and helper-inlining cost
// in this shader.
//
// Guard register optimisation
// ---------------------------
// s_r[] is padded with one zero-valued "guard" element at the start and end.
// Every runtime index is biased by +1 and clamped to [0, arraySize-1], so
// out-of-bounds accesses (e.g. r13-r15) silently land on a guard that reads
// as zero — matching NV2A hardware behaviour without a conditional branch.
//
// Constants (c0-c191) are read directly from the cbuffer C[] to avoid a
// massive 192×float4 per-thread indexable temp array that would destroy GPU
// occupancy.  Out-of-range reads (via a0.x) return zero via clamp+select.
// A small 4-entry context-write cache handles the rare case where a vertex
// program writes back to constant registers (NV2A context writes, out_orb==0).
// ============================================================
#define X_D3DVS_TEMPREG_COUNT  12
#define S_R_GUARD_SIZE (X_D3DVS_TEMPREG_COUNT + 2) // 14: [0]=lo guard, [1..12]=r0-r11, [13]=hi guard
#define GUARD_BIAS 1

static float4 s_r[S_R_GUARD_SIZE];  // temporary registers r0-r11 + 2 zero guards
static float4 s_v[16];      // input (vertex attribute) registers v0-v15
static float4 s_oRegs[16];  // output registers (indexed by OUTPUT_REG_*)
static int    s_a0;         // address register

// Context register write cache — tiny associative cache for the rare case
// where a vertex program writes back to constant registers (out_orb==0).
// NV2A programs can write to c[] so subsequent reads see updated values.
// Only 4 slots needed: no known Xbox title writes more than 2-3 constants
// per vertex program.  If all slots are full, the oldest is evicted (FIFO).
// ENABLE_CTX_WRITES can be set to 0 at compile time to produce a variant
// that skips context write support entirely (faster for the common case).
#ifndef ENABLE_CTX_WRITES
#define ENABLE_CTX_WRITES 1
#endif
#if ENABLE_CTX_WRITES
#define CTX_CACHE_SIZE 4
static float4 s_ctx_val[CTX_CACHE_SIZE];
static uint   s_ctx_idx[CTX_CACHE_SIZE]; // 0xFFFFFFFF = empty slot
static uint   s_ctx_count;               // number of valid entries (0..CTX_CACHE_SIZE)
#endif

// Named indices into s_oRegs[].  These match the NV2A output address
// encoding: the 4-bit out_address field maps directly to these slots.
// The layout mirrors the NV2A vertex attribute (input) slot numbering
// used for v[] registers (see CxbxFixedFunctionVertexShader.hlsl and
// VertexShader.cpp OReg_Name[]):
//   0=oPos/position  3=oD0/diffuse    4=oD1/specular   5=oFog/fogCoord
//   6=oPts/pointSize 7=oB0/backDiff   8=oB1/backSpec   9-12=oT0-oT3/texcoord0-3
// Slots 1-2 (weight/normal on the input side) are unused for output.
// Slot 15 is a0.x in the OReg encoding (handled via s_a0/ARL, not s_oRegs).
#define OUTPUT_REG_OPOS  0u  // position
#define OUTPUT_REG_OD0   3u  // front diffuse colour
#define OUTPUT_REG_OD1   4u  // front specular colour
#define OUTPUT_REG_OFOG  5u  // fog coordinate (only .x used by rasterizer)
#define OUTPUT_REG_OPTS  6u  // point sprite size
#define OUTPUT_REG_OB0   7u  // back diffuse colour
#define OUTPUT_REG_OB1   8u  // back specular colour
#define OUTPUT_REG_OT0   9u  // texture coordinate 0
#define OUTPUT_REG_OT1  10u  // texture coordinate 1
#define OUTPUT_REG_OT2  11u  // texture coordinate 2
#define OUTPUT_REG_OT3  12u  // texture coordinate 3
#define OUTPUT_REG_A0X  15u  // a0.x (writes go through s_a0/ARL, not write_output)

// ============================================================
// Swizzle helper: rearrange float4 components by packed index.
// Packed format: bits [7:6]=X [5:4]=Y [3:2]=Z [1:0]=W.
// Direct float4 indexing (no temp scalar array).
// ============================================================
float4 apply_swizzle(float4 v, uint swz)
{
    return float4(v[(swz >> 6) & 3], v[(swz >> 4) & 3],
                  v[(swz >> 2) & 3], v[ swz       & 3]);
}

// ============================================================
// Fetch an input source register value (A, B, or C input).
// ============================================================
float4 fetch_input(uint mux, uint r_idx, uint v_idx, uint const_idx,
                   uint swz, bool is_neg, bool use_a0x)
{
    float4 raw;

    if (mux == VSI_MUX_R) {
        // r0-r11; r12 aliases oPos; >12 is NV2A-undefined → guard reads zero.
        raw = (r_idx == 12) ? s_oRegs[OUTPUT_REG_OPOS]
                            : s_r[clamp((int)r_idx + GUARD_BIAS, 0, S_R_GUARD_SIZE - 1)];
    }
    else if (mux == VSI_MUX_V) {
        raw = s_v[v_idx & 0xF];
    }
    else {
        // Constant register c0..c191, optionally offset by a0.
        // Read from cbuffer C[] (fast path) with context-write cache overlay.
        // Out-of-range indices return zero.
        int c_index = (int)(const_idx & 0xFF) + s_a0 * (int)use_a0x;
        if (c_index >= 0 && c_index < X_D3DVS_CONSTREG_COUNT) {
            // Check context write cache first (most recent writes win)
            raw = C[c_index];
#if ENABLE_CTX_WRITES
            uint ci = (uint)c_index;
            [unroll] for (uint k = 0; k < CTX_CACHE_SIZE; k++) {
                if (k < s_ctx_count && s_ctx_idx[k] == ci)
                    raw = s_ctx_val[k];
            }
#endif
        } else {
            raw = float4(0, 0, 0, 0);
        }
    }

    float4 sw = apply_swizzle(raw, swz);
    // Branchless negate: multiply by ±1.  Compiles to a single XOR on the
    // sign bits when the compiler sees the pattern (no actual mul emitted).
    float neg = is_neg ? -1.0f : 1.0f;
    return sw * neg;
}

// ============================================================
// Write `src` into `dest` respecting the 4-bit component mask.
// Four independent scalar bit-tests compile to four `movc`s without
// FXC having to build a bool4 from a uint constant first.
// ============================================================
void write_masked(inout float4 dest, float4 src, uint mask)
{
    if (mask & VSI_MASK_X) dest.x = src.x;
    if (mask & VSI_MASK_Y) dest.y = src.y;
    if (mask & VSI_MASK_Z) dest.z = src.z;
    if (mask & VSI_MASK_W) dest.w = src.w;
}

// ============================================================
// Write to a temporary register (r0-r11) or oPos (r12).
// Indices > 12 are undefined on NV2A and silently ignored.
// ============================================================
void write_r(uint dest, float4 result, uint mask)
{
    if (dest == 12)     write_masked(s_oRegs[OUTPUT_REG_OPOS], result, mask);
    else if (dest < 12) write_masked(s_r[dest + GUARD_BIAS],  result, mask);
}

// ============================================================
// Write to context register cache (constant write-back).
// If the index already exists in the cache, update in-place.
// Otherwise allocate the next free slot (FIFO eviction when full).
// ============================================================
#if ENABLE_CTX_WRITES
void write_ctx(uint addr, float4 result, uint mask)
{
    // Find existing entry or allocate a new slot
    uint slot = s_ctx_count < CTX_CACHE_SIZE ? s_ctx_count : (CTX_CACHE_SIZE - 1);
    bool found = false;
    [unroll] for (uint k = 0; k < CTX_CACHE_SIZE; k++) {
        if (k < s_ctx_count && s_ctx_idx[k] == addr) {
            slot = k;
            found = true;
        }
    }
    if (!found) {
        // New slot — seed with cbuffer value so partial-mask writes
        // preserve the unwritten components (e.g. MOV c58.xy keeps .zw).
        s_ctx_idx[slot] = addr;
        s_ctx_val[slot] = C[addr];
        s_ctx_count = min(s_ctx_count + 1, CTX_CACHE_SIZE);
    }
    write_masked(s_ctx_val[slot], result, mask);
}
#endif // ENABLE_CTX_WRITES

// NV2A-accurate multiply and dot product helpers are in CxbxNV2AMathHelpers.hlsli

// ============================================================
// MAC unit operations
//
// [branch]: All vertices in a draw call execute the same instruction,
// so opcode is wavefront-uniform.  [branch] prevents FXC from
// flattening the switch (which would evaluate ALL 12 cases per
// iteration) and generates a single indexed jump instead.
// ============================================================
float4 exec_mac(uint opcode, float4 a, float4 b, float4 c_in)
{
    [branch] switch (opcode) {
        case VSI_MAC_MOV: return mac_mov(a);
        case VSI_MAC_MUL: return mac_mul(a, b);
        case VSI_MAC_ADD: return mac_add(a, c_in);
        case VSI_MAC_MAD: return mac_mad(a, b, c_in);
        case VSI_MAC_DP3: return mac_dp3(a, b);
        case VSI_MAC_DPH: return mac_dph(a, b);
        case VSI_MAC_DP4: return mac_dp4(a, b);
        case VSI_MAC_DST: return mac_dst(a, b);
        case VSI_MAC_MIN: return mac_min(a, b);
        case VSI_MAC_MAX: return mac_max(a, b);
        case VSI_MAC_SLT: return mac_slt(a, b);
        case VSI_MAC_SGE: return mac_sge(a, b);
        default: return float4(0, 0, 0, 0);
    }
}

// ============================================================
// ILU unit operations
// ============================================================

// [branch]: Same rationale as exec_mac — opcode is wavefront-uniform.
float4 exec_ilu(uint opcode, float4 c_in)
{
    [branch] switch (opcode) {
        case VSI_ILU_MOV: return ilu_mov(c_in);
        case VSI_ILU_RCP: return ilu_rcp(c_in);
        case VSI_ILU_RCC: return ilu_rcc(c_in);
        case VSI_ILU_RSQ: return ilu_rsq(c_in);
        case VSI_ILU_EXP: return ilu_exp(c_in);
        case VSI_ILU_LOG: return ilu_log(c_in);
        case VSI_ILU_LIT: return ilu_lit(c_in);
        default: return float4(0, 0, 0, 0);
    }
}

// ============================================================
// Fog output register write helper
//
// NV2A remaps writes to oFog so that the most significant masked
// component ends up in .x (the only component the rasterizer reads).
// For example, a write with mask ,y puts the .y result into oFog.x.
// This matches xemu's fog_mask_str remapping table.
// ============================================================
void write_fog_output(float4 result, uint mask)
{
    // Highest set bit in mask (x=8,y=4,z=2,w=1) selects component → .x.
    // Ternary chain compiles to a cascade of movc (branchless).
    s_oRegs[OUTPUT_REG_OFOG].x = (mask & VSI_MASK_X) ? result.x
                                : (mask & VSI_MASK_Y) ? result.y
                                : (mask & VSI_MASK_Z) ? result.z
                                :                        result.w;
}

// ============================================================
// Write to output register, with fog mask remapping.
// ============================================================
void write_output(uint out_address, float4 result, uint mask)
{
    uint addr = out_address & 0xF;
    if (addr == OUTPUT_REG_OFOG)
        write_fog_output(result, mask);
    else
        write_masked(s_oRegs[addr], result, mask);
}

// ============================================================
// Main vertex shader entry point
// ============================================================
VS_OUTPUT main(const VS_INPUT xIn)
{
    // Output register defaults.
    //   0=oPos, 1-2=unused, 3=oD0, 4=oD1, 5=oFog, 6=oPts, 7=oB0, 8=oB1,
    //   9-12=oT0-oT3, 13-15=padding so (out_address & 0xF) can never
    //   land on an unmapped slot and write into other state.
    s_oRegs[OUTPUT_REG_OPOS] = float4(0, 0, 0, 1);
    s_oRegs[1]               = float4(0, 0, 0, 0); // unused
    s_oRegs[2]               = float4(0, 0, 0, 0); // unused
    s_oRegs[OUTPUT_REG_OD0]  = float4(1, 1, 1, 1); // TODO: xemu uses (0,0,0,1)
    s_oRegs[OUTPUT_REG_OD1]  = float4(0, 0, 0, 1); // specular must default black — white saturates V1R0_SUM (Water)
    s_oRegs[OUTPUT_REG_OFOG] = float4(1, 1, 1, 1); // TODO: xemu uses (0,0,0,1) — (0,0,0,1) turns labels white via fog (Water)
    s_oRegs[OUTPUT_REG_OPTS] = float4(0, 0, 0, 0);
    s_oRegs[OUTPUT_REG_OB0]  = float4(1, 1, 1, 1); // TODO: xemu uses (0,0,0,1)
    s_oRegs[OUTPUT_REG_OB1]  = float4(0, 0, 0, 1); // back-face specular, same reasoning as oD1
    s_oRegs[OUTPUT_REG_OT0]  = float4(0, 0, 0, 1);
    s_oRegs[OUTPUT_REG_OT1]  = float4(0, 0, 0, 1);
    s_oRegs[OUTPUT_REG_OT2]  = float4(0, 0, 0, 1);
    s_oRegs[OUTPUT_REG_OT3]  = float4(0, 0, 0, 1);
    s_oRegs[13]              = float4(0, 0, 0, 0); // padding
    s_oRegs[14]              = float4(0, 0, 0, 0); // padding
    s_oRegs[OUTPUT_REG_A0X]  = float4(0, 0, 0, 0); // a0.x (writes go through s_a0/ARL, not here)

    s_a0 = 0;

    // Zero r0-r11 plus guard slots (Xbox semantics: all temp regs start at zero).
    [unroll] for (uint ri = 0; ri < S_R_GUARD_SIZE; ri++) s_r[ri] = float4(0, 0, 0, 0);

    // Context write cache starts empty.
#if ENABLE_CTX_WRITES
    s_ctx_count = 0;
    [unroll] for (uint ci = 0; ci < CTX_CACHE_SIZE; ci++) {
        s_ctx_idx[ci] = 0xFFFFFFFFu;
        s_ctx_val[ci] = float4(0, 0, 0, 0);
    }
#endif

    // Constants are read from cbuffer C[] with a small context-write cache
    // overlay for the rare case where a vertex program writes back to c[].

    // Populate v0-v15 directly into the static array.
    FetchAllAttributes(ResolveVertexIndex(xIn.vertexId), s_v);

    // ============================================================
    // Instruction execution loop
    // Read program start address from PGRAPH register CSV0_C.
    // The shader loops from startSlot until FLD_FINAL or max slots.
    // ============================================================
    uint startSlot = (PG_UINT(NV_PGRAPH_CSV0_C) >> NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START_SHIFT)
                   & NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START_MASK;

    [loop]
    for (uint pc = 0; pc < XFPR_LENGTH && (startSlot + pc) < XFPR_LENGTH; pc++) {
        uint4 inst = g_XFPR[startSlot + pc];
        // inst.x = SubToken 0 (unused by fields)
        // inst.y = SubToken 1
        // inst.z = SubToken 2
        // inst.w = SubToken 3

        uint dw1 = inst.y;
        uint dw2 = inst.z;
        uint dw3 = inst.w;

        // Opcode fields
        uint ilu_op   = (dw1 >> VSI_FLD_ILU_SHIFT) & VSI_FLD_ILU_MASK;
        uint mac_op   = (dw1 >> VSI_FLD_MAC_SHIFT) & VSI_FLD_MAC_MASK;
        bool is_final = ((dw3 >> VSI_FLD_FINAL_BIT3) & 1) != 0;

        bool has_mac = (mac_op != VSI_MAC_NOP);
        bool has_ilu = (ilu_op != VSI_ILU_NOP);

        // Skip decode entirely when both units are idle (padding slots)
        if (!has_mac && !has_ilu) {
            if (is_final) break;
            continue;
        }

        // -- Register indices
        uint const_idx = (dw1 >> VSI_FLD_CONST_SHIFT) & VSI_FLD_CONST_MASK;
        uint v_idx     = (dw1 >> VSI_FLD_V_SHIFT)     & VSI_FLD_V_MASK;

        // -- Input A (dw1 + dw2)
        uint a_mux = (dw2 >> VSI_FLD_A_MUX_SHIFT) & VSI_FLD_A_MUX_MASK;
        uint a_reg = (dw2 >> VSI_FLD_A_R_SHIFT)   & VSI_FLD_A_R_MASK;
        bool a_neg = ((dw1 >> VSI_FLD_A_NEG_BIT1) & 1) != 0;
        uint a_swz = dw1 & 0xFF;          // Packed XYZW: [7:6]=X [5:4]=Y [3:2]=Z [1:0]=W

        // -- Input B (dw2)
        uint b_mux = (dw2 >> VSI_FLD_B_MUX_SHIFT) & VSI_FLD_B_MUX_MASK;
        uint b_reg = (dw2 >> VSI_FLD_B_R_SHIFT)   & VSI_FLD_B_R_MASK;
        bool b_neg = ((dw2 >> VSI_FLD_B_NEG_BIT2) & 1) != 0;
        uint b_swz = (dw2 >> 17) & 0xFF;  // Packed XYZW from bits [24:17]

        // -- Input C (dw2 + dw3)
        uint c_mux   = (dw3 >> VSI_FLD_C_MUX_SHIFT3)    & VSI_FLD_C_MUX_MASK;
        uint c_r_hi  = (dw2 >> VSI_FLD_C_R_HIGH_SHIFT2) & VSI_FLD_C_R_HIGH_MASK;
        uint c_r_lo  = (dw3 >> VSI_FLD_C_R_LOW_SHIFT3)  & VSI_FLD_C_R_LOW_MASK;
        uint c_reg   = (c_r_hi << 2) | c_r_lo;
        bool c_neg   = ((dw2 >> VSI_FLD_C_NEG_BIT2) & 1) != 0;
        uint c_swz   = (dw2 >> 2) & 0xFF; // Packed XYZW from bits [9:2]

        // -- Output fields
        uint out_mac_mask = (dw3 >> VSI_FLD_OUT_MAC_MASK_SHIFT) & VSI_FLD_OUT_MAC_MASK_MASK;
        uint out_r_addr   = (dw3 >> VSI_FLD_OUT_R_SHIFT)        & VSI_FLD_OUT_R_MASK;
        uint out_ilu_mask = (dw3 >> VSI_FLD_OUT_ILU_MASK_SHIFT) & VSI_FLD_OUT_ILU_MASK_MASK;
        uint out_o_mask   = (dw3 >> VSI_FLD_OUT_O_MASK_SHIFT)   & VSI_FLD_OUT_O_MASK_MASK;
        bool out_orb      = ((dw3 >> VSI_FLD_OUT_ORB_BIT3) & 1) != 0; // 0=context, 1=output
        uint out_address  = (dw3 >> VSI_FLD_OUT_ADDRESS_SHIFT)  & VSI_FLD_OUT_ADDRESS_MASK;
        uint out_mux      = (dw3 >> VSI_FLD_OUT_MUX_BIT3) & 1;        // 0=MAC, 1=ILU
        bool use_a0x      = ((dw3 >> VSI_FLD_A0X_BIT3) & 1) != 0;

        bool is_paired     = has_mac && has_ilu;
        bool mac_is_output = (out_mux == 0);
        bool do_out        = (out_o_mask != 0) && out_orb;
#if ENABLE_CTX_WRITES
        bool do_ctx        = (out_o_mask != 0) && !out_orb;
#endif

        // ========================================================
        // Snapshot inputs BEFORE either unit writes back.
        // Critical for correctness: when paired, ILU must see the
        // pre-MAC value of its source register, not the post-MAC one.
        //
        // [branch] on has_mac / has_ilu: these bools derive from the
        // instruction word, so every vertex in the wave takes the
        // same path.  [branch] prevents FXC from flattening (which
        // would evaluate both fetch+execute paths every iteration).
        // ========================================================
        float4 in_a = float4(0, 0, 0, 0);
        float4 in_b = float4(0, 0, 0, 0);

        [branch] if (has_mac) {
            in_a = fetch_input(a_mux, a_reg, v_idx, const_idx, a_swz, a_neg, use_a0x);
            in_b = fetch_input(b_mux, b_reg, v_idx, const_idx, b_swz, b_neg, use_a0x);
        }

        // C input is used by both MAC and ILU; fetch unconditionally.
        float4 in_c = fetch_input(c_mux, c_reg, v_idx, const_idx, c_swz, c_neg, use_a0x);

        // ========================================================
        // Execute MAC
        // ========================================================
        [branch] if (has_mac) {
            [branch] if (mac_op == VSI_MAC_ARL) {
                s_a0 = mac_arl(in_a);
            }
            else {
                float4 mac_result = exec_mac(mac_op, in_a, in_b, in_c);

                if (out_mac_mask != 0)
                    write_r(out_r_addr, mac_result, out_mac_mask);

                if (mac_is_output && do_out)
                    write_output(out_address, mac_result, out_o_mask);

#if ENABLE_CTX_WRITES
                if (mac_is_output && do_ctx && out_address < X_D3DVS_CONSTREG_COUNT)
                    write_ctx(out_address, mac_result, out_o_mask);
#endif
            }
        }

        // ========================================================
        // Execute ILU
        // ========================================================
        [branch] if (has_ilu) {
            float4 ilu_result = exec_ilu(ilu_op, in_c);

            // When paired, ILU always writes to r1; otherwise shares out_r_addr with MAC.
            uint ilu_r_dest = is_paired ? 1 : out_r_addr;
            if (out_ilu_mask != 0)
                write_r(ilu_r_dest, ilu_result, out_ilu_mask);

            if (!mac_is_output && do_out)
                write_output(out_address, ilu_result, out_o_mask);

#if ENABLE_CTX_WRITES
            if (!mac_is_output && do_ctx && out_address < X_D3DVS_CONSTREG_COUNT)
                write_ctx(out_address, ilu_result, out_o_mask);
#endif
        }

        if (is_final) break;
    }

    // ============================================================
    // Copy to output struct (same footer as other vertex shaders).
    // Footer expects these named variables in scope.
    // ============================================================
    float4 oPos = s_oRegs[OUTPUT_REG_OPOS];
    float4 oD0  = s_oRegs[OUTPUT_REG_OD0];
    float4 oD1  = s_oRegs[OUTPUT_REG_OD1];
    float4 oFog = s_oRegs[OUTPUT_REG_OFOG];
    float4 oPts = s_oRegs[OUTPUT_REG_OPTS];
    float4 oB0  = s_oRegs[OUTPUT_REG_OB0];
    float4 oB1  = s_oRegs[OUTPUT_REG_OB1];
    float4 oT0  = s_oRegs[OUTPUT_REG_OT0];
    float4 oT1  = s_oRegs[OUTPUT_REG_OT1];
    float4 oT2  = s_oRegs[OUTPUT_REG_OT2];
    float4 oT3  = s_oRegs[OUTPUT_REG_OT3];

    VS_OUTPUT xOut;
#include "CxbxVertexOutputFooter.hlsli"

    return xOut;
}
