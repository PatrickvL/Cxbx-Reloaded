// VertexShaderInterpreter.hlsl
//
// DX11 / SM5 Xbox NV2A vertex shader interpreter ubershader.
//
// Instead of recompiling each Xbox vertex shader program into host HLSL,
// this single precompiled shader interprets the raw NV2A microcode at
// runtime. The 128-bit instruction slots are uploaded to cbuffer b3;
// vertex constants (c0–c191) are in the existing cbuffer b0.
//
// Architecture mirrors the register combiner interpreter:
//   - C++ uploads raw Xbox microcode bytes to a constant buffer
//   - This shader loops over instruction slots, decodes fields, executes ops
//   - No CPU-side D3DCompile, no shader cache, no async compilation
//
// Reference: nv2a_vsh_cpu (import/nv2a_vsh_cpu/src/nv2a_vsh_emulator.c)

#include "CxbxVertexShaderCommon.hlsli"
#include "CxbxVertexFetch.hlsli"

// Xbox constant registers (same as in CxbxVertexShaderTemplate.hlsl)
#define X_D3DSCM_CORRECTION 96
#define X_D3DVS_CONSTREG_COUNT 192
uniform float4 C[X_D3DVS_CONSTREG_COUNT] : register(c0);

#include "CxbxScreenspaceTransform.hlsli"
#include "VertexShaderInterpreterState.hlsli"

// ============================================================
// Swizzle helper: rearrange float4 components by packed index
// Packed format: bits [7:6]=X [5:4]=Y [3:2]=Z [1:0]=W
// ============================================================
float4 apply_swizzle(float4 v, uint swz)
{
    float arr[4] = { v.x, v.y, v.z, v.w };
    return float4(arr[(swz >> 6) & 3], arr[(swz >> 4) & 3],
                  arr[(swz >> 2) & 3], arr[swz & 3]);
}

// ============================================================
// Fetch an input source register value
// ============================================================
float4 fetch_input(
    uint mux, uint r_idx, uint v_idx, uint const_idx,
    uint swz, bool is_neg, bool use_a0x, int a0,
    float4 r[12], float4 oPos_reg, float4 v_regs[16])
{
    float4 raw;

    if (mux == VSI_MUX_R) {
        // Temporary register r0-r11; r12 aliases oPos
        if (r_idx == 12)
            raw = oPos_reg;
        else if (r_idx < 12)
            raw = r[r_idx];
        else
            raw = float4(0, 0, 0, 0);
    }
    else if (mux == VSI_MUX_V) {
        // Vertex input register v0-v15
        raw = v_regs[v_idx & 0xF];
    }
    else {
        // Constant register c0-c191
        // The Xbox encoding collapses to: (const_idx & 0xFF) maps to 0..191
        int c_index = (int)(const_idx & 0xFF);
        if (use_a0x)
            c_index += a0;
        raw = (c_index >= 0 && c_index < (int)X_D3DVS_CONSTREG_COUNT)
            ? C[c_index] : float4(0, 0, 0, 0);
    }

    float4 swizzled = apply_swizzle(raw, swz);
    return is_neg ? -swizzled : swizzled;
}

// ============================================================
// Write result to register with writemask
// ============================================================
void write_masked(inout float4 dest, float4 src, uint mask)
{
    dest = float4(
        (mask & VSI_MASK_X) ? src.x : dest.x,
        (mask & VSI_MASK_Y) ? src.y : dest.y,
        (mask & VSI_MASK_Z) ? src.z : dest.z,
        (mask & VSI_MASK_W) ? src.w : dest.w
    );
}

// ============================================================
// Write result to a temporary register (r0-r11) or oPos (r12)
// Indices > 12 are undefined on NV2A and silently ignored.
// ============================================================
void write_r(uint dest, inout float4 r[12], inout float4 oPos, float4 result, uint mask)
{
    if (dest == 12)     write_masked(oPos,    result, mask);
    else if (dest < 12) write_masked(r[dest], result, mask);
}

// ============================================================
// MAC unit operations
// ============================================================
float4 exec_mac(uint opcode, float4 a, float4 b, float4 c_in)
{
    switch (opcode) {
        case VSI_MAC_MOV: return a;
        case VSI_MAC_MUL: return a * b;
        case VSI_MAC_ADD: return a + c_in;
        case VSI_MAC_MAD: return a * b + c_in;
        case VSI_MAC_DP3: { float d = dot(a.xyz, b.xyz); return float4(d, d, d, d); }
        case VSI_MAC_DPH: { float d = dot(a.xyz, b.xyz) + b.w; return float4(d, d, d, d); }
        case VSI_MAC_DP4: { float d = dot(a, b); return float4(d, d, d, d); }
        case VSI_MAC_DST: return float4(1.0, a.y * b.y, a.z, b.w);
        case VSI_MAC_MIN: return min(a, b);
        case VSI_MAC_MAX: return max(a, b);
        case VSI_MAC_SLT: return 1.0 - step(b, a);  // 1 where a < b
        case VSI_MAC_SGE: return step(b, a);           // 1 where a >= b
        case VSI_MAC_ARL: return a; // ARL result stored to a0 by caller
        default: return float4(0, 0, 0, 0);
    }
}

// ============================================================
// ILU unit operations
// ============================================================

// Floor with bias (matching CxbxVertexShaderTemplate.hlsl)
#define BIAS 0.001

float vsi_floor(float src)
{
    return floor(src + BIAS);
}

float4 exec_ilu(uint opcode, float4 c_in)
{
    float s = c_in.x; // Scalar input

    switch (opcode) {
        case VSI_ILU_MOV: return c_in;
        case VSI_ILU_RCP:
        case VSI_ILU_RCC: {
            float rv = 1.0 / s;
            rv = (rv >= 0)
                ? clamp(rv, 5.42101e-020f, 1.84467e+019f)
                : clamp(rv, -1.84467e+019f, -5.42101e-020f);
            return float4(rv, rv, rv, rv);
        }
        case VSI_ILU_RSQ: {
            float a = abs(s);
            float r = rsqrt(a);
            return float4(r, r, r, r);
        }
        case VSI_ILU_EXP: {
            float fl = vsi_floor(s);
            return float4(exp2(fl), s - fl, exp2(s), 1.0);
        }
        case VSI_ILU_LOG: {
            float ex;
            float mantissa = frexp(s, ex);
            float z = log2(s);
            return float4(ex, mantissa, z, 1.0);
        }
        case VSI_ILU_LIT: {
            float diffuse = c_in.x;
            float blinn = c_in.y;
            float specPower = clamp(c_in.w, -(128.0 - 1.0/256.0), 128.0 - 1.0/256.0);
            float litZ = (diffuse > 0 && blinn > 0) ? pow(abs(blinn), specPower) : 0;
            return float4(1.0, max(0.0, diffuse), litZ, 1.0);
        }
        default: return float4(0, 0, 0, 0);
    }
}

// ============================================================
// Main vertex shader entry point
// ============================================================
VS_OUTPUT main(const VS_INPUT xIn)
{
    // Output registers: sparse array indexed by NV2A output address
    // 0=oPos, 1-2=unused, 3=oD0, 4=oD1, 5=oFog, 6=oPts, 7=oB0, 8=oB1, 9-12=oT0-oT3
    float4 oRegs[13];
    oRegs[0]  = float4(0, 0, 0, 1); // oPos
    oRegs[1]  = float4(0, 0, 0, 0); // unused
    oRegs[2]  = float4(0, 0, 0, 0); // unused
    oRegs[3]  = float4(0, 0, 0, 1); // oD0
    oRegs[4]  = float4(0, 0, 0, 1); // oD1
    oRegs[5]  = float4(1, 1, 1, 1); // oFog
    oRegs[6]  = float4(0, 0, 0, 0); // oPts
    oRegs[7]  = float4(0, 0, 0, 1); // oB0
    oRegs[8]  = float4(0, 0, 0, 1); // oB1
    oRegs[9]  = float4(0, 0, 0, 1); // oT0
    oRegs[10] = float4(0, 0, 0, 1); // oT1
    oRegs[11] = float4(0, 0, 0, 1); // oT2
    oRegs[12] = float4(0, 0, 0, 1); // oT3

    // Address register
    int a0 = 0;

    // Temporary registers r0-r11
    float4 r[12];
    [unroll] for (uint ri = 0; ri < 12; ri++) r[ri] = float4(0, 0, 0, 0);

    // Input registers v0-v15
    float4 v_regs[16];
    {
        float4 v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;
#include "CxbxVertexInputLoad.hlsli"
        v_regs[0]=v0; v_regs[1]=v1; v_regs[2]=v2; v_regs[3]=v3;
        v_regs[4]=v4; v_regs[5]=v5; v_regs[6]=v6; v_regs[7]=v7;
        v_regs[8]=v8; v_regs[9]=v9; v_regs[10]=v10; v_regs[11]=v11;
        v_regs[12]=v12; v_regs[13]=v13; v_regs[14]=v14; v_regs[15]=v15;
    }

    // ============================================================
    // Instruction execution loop
    // ============================================================
    uint instCount = min(InstructionCount, VSI_MAX_SLOTS);

    [loop]
    for (uint pc = 0; pc < instCount; pc++) {
        uint4 inst = Instructions[pc];
        // inst.x = SubToken 0 (unused by fields)
        // inst.y = SubToken 1
        // inst.z = SubToken 2
        // inst.w = SubToken 3

        uint dw1 = inst.y;
        uint dw2 = inst.z;
        uint dw3 = inst.w;

        // Decode opcodes
        uint ilu_op = (dw1 >> VSI_FLD_ILU_SHIFT) & VSI_FLD_ILU_MASK;
        uint mac_op = (dw1 >> VSI_FLD_MAC_SHIFT) & VSI_FLD_MAC_MASK;

        // Skip decode entirely when both units are idle (padding slots)
        if (mac_op == VSI_MAC_NOP && ilu_op == VSI_ILU_NOP) {
            if (((dw3 >> VSI_FLD_FINAL_BIT3) & 1) != 0) break;
            continue;
        }

        // Decode register indices
        uint const_idx = (dw1 >> VSI_FLD_CONST_SHIFT) & VSI_FLD_CONST_MASK;
        uint v_idx     = (dw1 >> VSI_FLD_V_SHIFT)     & VSI_FLD_V_MASK;

        // Input A (SubToken 1 + SubToken 2)
        uint a_mux   = (dw2 >> VSI_FLD_A_MUX_SHIFT) & VSI_FLD_A_MUX_MASK;
        uint a_reg   = (dw2 >> VSI_FLD_A_R_SHIFT) & VSI_FLD_A_R_MASK;
        bool a_neg   = ((dw1 >> VSI_FLD_A_NEG_BIT1) & 1) != 0;
        uint a_swz   = dw1 & 0xFF; // Packed XYZW swizzle: bits [7:6]=X [5:4]=Y [3:2]=Z [1:0]=W

        // Input B (SubToken 2)
        uint b_mux   = (dw2 >> VSI_FLD_B_MUX_SHIFT) & VSI_FLD_B_MUX_MASK;
        uint b_reg   = (dw2 >> VSI_FLD_B_R_SHIFT) & VSI_FLD_B_R_MASK;
        bool b_neg   = ((dw2 >> VSI_FLD_B_NEG_BIT2) & 1) != 0;
        uint b_swz   = (dw2 >> 17) & 0xFF; // Packed XYZW swizzle from bits [24:17]

        // Input C (SubToken 2 + SubToken 3)
        uint c_mux     = (dw3 >> VSI_FLD_C_MUX_SHIFT3) & VSI_FLD_C_MUX_MASK;
        uint c_r_high  = (dw2 >> VSI_FLD_C_R_HIGH_SHIFT2) & VSI_FLD_C_R_HIGH_MASK;
        uint c_r_low   = (dw3 >> VSI_FLD_C_R_LOW_SHIFT3) & VSI_FLD_C_R_LOW_MASK;
        uint c_reg     = (c_r_high << 2) | c_r_low;
        bool c_neg     = ((dw2 >> VSI_FLD_C_NEG_BIT2) & 1) != 0;
        uint c_swz     = (dw2 >> 2) & 0xFF; // Packed XYZW swizzle from bits [9:2]

        // Output fields
        uint out_mac_mask = (dw3 >> VSI_FLD_OUT_MAC_MASK_SHIFT) & VSI_FLD_OUT_MAC_MASK_MASK;
        uint out_r_addr   = (dw3 >> VSI_FLD_OUT_R_SHIFT) & VSI_FLD_OUT_R_MASK;
        uint out_ilu_mask = (dw3 >> VSI_FLD_OUT_ILU_MASK_SHIFT) & VSI_FLD_OUT_ILU_MASK_MASK;
        uint out_o_mask   = (dw3 >> VSI_FLD_OUT_O_MASK_SHIFT) & VSI_FLD_OUT_O_MASK_MASK;
        bool out_orb      = ((dw3 >> VSI_FLD_OUT_ORB_BIT3) & 1) != 0; // 0=context, 1=output
        uint out_address  = (dw3 >> VSI_FLD_OUT_ADDRESS_SHIFT) & VSI_FLD_OUT_ADDRESS_MASK;
        uint out_mux      = (dw3 >> VSI_FLD_OUT_MUX_BIT3) & 1; // 0=MAC, 1=ILU
        bool use_a0x      = ((dw3 >> VSI_FLD_A0X_BIT3) & 1) != 0;
        bool is_final     = ((dw3 >> VSI_FLD_FINAL_BIT3) & 1) != 0;

        bool is_paired = (mac_op != VSI_MAC_NOP) && (ilu_op != VSI_ILU_NOP);

        // ============================================================
        // Snapshot inputs before executing (prevents order-dependent behavior)
        // MAC uses inputs A, B, C; ILU uses input C (same parameters)
        // ============================================================
        float4 in_a, in_b, in_c;

        if (mac_op != VSI_MAC_NOP) {
            in_a = fetch_input(a_mux, a_reg, v_idx, const_idx, a_swz, a_neg, use_a0x, a0, r, oRegs[0], v_regs);
            in_b = fetch_input(b_mux, b_reg, v_idx, const_idx, b_swz, b_neg, use_a0x, a0, r, oRegs[0], v_regs);
            in_c = fetch_input(c_mux, c_reg, v_idx, const_idx, c_swz, c_neg, use_a0x, a0, r, oRegs[0], v_regs);
        }
        else if (ilu_op != VSI_ILU_NOP) {
            // ILU-only: C-input not yet fetched
            in_c = fetch_input(c_mux, c_reg, v_idx, const_idx, c_swz, c_neg, use_a0x, a0, r, oRegs[0], v_regs);
        }

        // ============================================================
        // Execute MAC operation
        // ============================================================
        if (mac_op != VSI_MAC_NOP) {
            float4 mac_result = exec_mac(mac_op, in_a, in_b, in_c);

            // ARL writes to address register
            if (mac_op == VSI_MAC_ARL) {
                a0 = (int)vsi_floor(mac_result.x);
            }
            else {
                // Write to R register (unless paired and R=1, which is reserved for ILU)
                uint mac_r_dest = out_r_addr;
                if (!(is_paired && mac_r_dest == 1) && out_mac_mask != 0)
                    write_r(mac_r_dest, r, oRegs[0], mac_result, out_mac_mask);

                // Write to output register (if MAC is the output source)
                if (out_mux == 0 && out_o_mask != 0 && out_orb)
                    write_masked(oRegs[out_address & 0xF], mac_result, out_o_mask);
            }
        }

        // ============================================================
        // Execute ILU operation
        // ============================================================
        if (ilu_op != VSI_ILU_NOP) {
            float4 ilu_result = exec_ilu(ilu_op, in_c);

            // ILU writes to R register
            // When paired, ILU always writes to R1
            uint ilu_r_dest = is_paired ? 1 : out_r_addr;
            if (out_ilu_mask != 0)
                write_r(ilu_r_dest, r, oRegs[0], ilu_result, out_ilu_mask);

            // Write to output register (if ILU is the output source)
            if (out_mux == 1 && out_o_mask != 0 && out_orb)
                write_masked(oRegs[out_address & 0xF], ilu_result, out_o_mask);
        }

        // Stop at the final instruction
        if (is_final)
            break;
    }

    // ============================================================
    // Copy to output struct (same footer as CxbxVertexShaderTemplate.hlsl)
    // ============================================================
    // Unpack named outputs for the footer (expects named variables in scope)
    float4 oPos = oRegs[0];
    float4 oD0  = oRegs[3];
    float4 oD1  = oRegs[4];
    float4 oFog = oRegs[5];
    float4 oPts = oRegs[6];
    float4 oB0  = oRegs[7];
    float4 oB1  = oRegs[8];
    float4 oT0  = oRegs[9];
    float4 oT1  = oRegs[10];
    float4 oT2  = oRegs[11];
    float4 oT3  = oRegs[12];

    VS_OUTPUT xOut;
#include "CxbxVertexOutputFooter.hlsli"

    return xOut;
}
