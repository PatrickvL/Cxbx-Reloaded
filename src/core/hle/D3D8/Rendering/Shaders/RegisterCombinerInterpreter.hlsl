// RegisterCombinerInterpreter.hlsl — SM5.0 Xbox Pixel Shader (Register Combiner) Interpreter
//
// This ubershader interprets Xbox NV2A register combiner programs at runtime,
// reading the pixel shader definition from a constant buffer containing the
// relevant PGRAPH registers (D3DRS_PS* render states).
//
// Ported from the SM3.0 RegisterCombinerInterpreted.fx (PatrickvL/HLSL_PS2 branch).
// Key SM5 changes:
//   - Native uint32 bitwise operations replace fmod()-based bit extraction
//   - Texture2D/3D/Cube.Sample() replace tex2D/3D/CUBE
//   - cbuffer replaces SetPixelShaderConstantF
//   - Integer switch/case for all decode paths

// ---------------------------------------------------------------
// Constant buffer: Xbox pixel shader register combiner state
// Uploaded from D3D__RenderState[] (PGRAPH combiner registers).
// Each uint4 holds a packed 32-bit register value in .x, with the
// remaining .yzw available for additional data.
// ---------------------------------------------------------------
cbuffer CxbxPSState : register(b1)
{
    // Combiner inputs/outputs (8 stages)
    uint4 PSAlphaInputs[8];                  // c0-c7:   D3DRS_PSALPHAINPUTS[0..7]
    uint4 PSFinalCombinerInputsABCD;         // c8:      D3DRS_PSFINALCOMBINERINPUTSABCD
    uint4 PSFinalCombinerInputsEFG;          // c9:      D3DRS_PSFINALCOMBINERINPUTSEFG
    float4 PSConstant0[8];                   // c10-c17: D3DRS_PSCONSTANT0[0..7] (RGBA float colors)
    float4 PSConstant1[8];                   // c18-c25: D3DRS_PSCONSTANT1[0..7] (RGBA float colors)
    uint4 PSAlphaOutputs[8];                 // c26-c33: D3DRS_PSALPHAOUTPUTS[0..7]
    uint4 PSRGBInputs[8];                    // c34-c41: D3DRS_PSRGBINPUTS[0..7]
    uint4 PSCompareMode;                     // c42:     D3DRS_PSCOMPAREMODE
    float4 PSFinalCombinerConstant[2];       // c43-c44: D3DRS_PSFINALCOMBINERCONSTANT0/1
    uint4 PSRGBOutputs[8];                   // c45-c52: D3DRS_PSRGBOUTPUTS[0..7]
    uint4 PSCombinerCount;                   // c53:     D3DRS_PSCOMBINERCOUNT
    uint4 PSTextureModes;                    // c54:     D3DRS_PSTEXTUREMODES
    uint4 PSDotMapping;                      // c55:     D3DRS_PSDOTMAPPING
    uint4 PSInputTexture;                    // c56:     D3DRS_PSINPUTTEXTURE
    float4 PSColorSign[4];                   // c57-c60: per-stage color sign conversion
};

// ---------------------------------------------------------------
// Texture samplers
// ---------------------------------------------------------------
Texture2D    g_Tex2D[4] : register(t0);
Texture3D    g_Tex3D[4] : register(t4);
TextureCube  g_TexCube[4] : register(t8);
SamplerState g_Sampler[4] : register(s0);

// ---------------------------------------------------------------
// PS_REGISTER constants
// ---------------------------------------------------------------
#define PS_REGISTER_ZERO        0
#define PS_REGISTER_DISCARD     0
#define PS_REGISTER_C0          1
#define PS_REGISTER_C1          2
#define PS_REGISTER_FOG         3
#define PS_REGISTER_V0          4
#define PS_REGISTER_V1          5
#define PS_REGISTER_T0          8
#define PS_REGISTER_T1          9
#define PS_REGISTER_T2          10
#define PS_REGISTER_T3          11
#define PS_REGISTER_R0          12
#define PS_REGISTER_R1          13
#define PS_REGISTER_V1R0_SUM    14
#define PS_REGISTER_EF_PROD     15
#define PS_REGISTER_COUNT       16

// Input mapping constants
#define PS_INPUTMAPPING_UNSIGNED_IDENTITY  0u
#define PS_INPUTMAPPING_UNSIGNED_INVERT    1u
#define PS_INPUTMAPPING_EXPAND_NORMAL      2u
#define PS_INPUTMAPPING_EXPAND_NEGATE      3u
#define PS_INPUTMAPPING_HALFBIAS_NORMAL    4u
#define PS_INPUTMAPPING_HALFBIAS_NEGATE    5u
#define PS_INPUTMAPPING_SIGNED_IDENTITY    6u
#define PS_INPUTMAPPING_SIGNED_NEGATE      7u

// Channel select
#define PS_CHANNEL_RGB   0u
#define PS_CHANNEL_ALPHA 1u

// Output mapping flags
#define PS_COMBINEROUTPUT_BIAS_BIT       0x08u
#define PS_COMBINEROUTPUT_SCALE_MASK     0x30u
#define PS_COMBINEROUTPUT_SHIFTLEFT_1    0x10u
#define PS_COMBINEROUTPUT_SHIFTLEFT_2    0x20u
#define PS_COMBINEROUTPUT_SHIFTRIGHT_1   0x30u

// Texture modes
#define PS_TEXTUREMODES_NONE                  0u
#define PS_TEXTUREMODES_PROJECT2D             1u
#define PS_TEXTUREMODES_PROJECT3D             2u
#define PS_TEXTUREMODES_CUBEMAP               3u
#define PS_TEXTUREMODES_PASSTHRU              4u
#define PS_TEXTUREMODES_CLIPPLANE             5u
#define PS_TEXTUREMODES_BUMPENVMAP            6u
#define PS_TEXTUREMODES_BUMPENVMAP_LUM        7u
#define PS_TEXTUREMODES_BRDF                  8u
#define PS_TEXTUREMODES_DOT_ST                9u
#define PS_TEXTUREMODES_DOT_ZW                10u
#define PS_TEXTUREMODES_DOT_RFLCT_DIFF        11u
#define PS_TEXTUREMODES_DOT_RFLCT_SPEC        12u
#define PS_TEXTUREMODES_DOT_STR_3D            13u
#define PS_TEXTUREMODES_DOT_STR_CUBE          14u
#define PS_TEXTUREMODES_DPNDNT_AR             15u
#define PS_TEXTUREMODES_DPNDNT_GB             16u
#define PS_TEXTUREMODES_DOTPRODUCT            17u
#define PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST  18u

// Stage identifiers
#define STAGE_FINAL_COMBINER       8u
#define STAGE_FINAL_COMBINER_ABCD  9u
#define MAX_COMBINER_STAGE_COUNT   8u

// ---------------------------------------------------------------
// Interpreter state
// ---------------------------------------------------------------
struct ps_state
{
    float4 Regs[16]; // All registers indexed by PS_REGISTER_*
    uint   stage;
    bool   FlagMuxMsb;
    bool   FlagUniqueC0;
    bool   FlagUniqueC1;
};

// ---------------------------------------------------------------
// Bitfield extraction helpers (native integer ops in SM5)
// ---------------------------------------------------------------
uint ExtractBits(uint value, uint offset, uint count)
{
    return (value >> offset) & ((1u << count) - 1u);
}

// ---------------------------------------------------------------
// Register access
// ---------------------------------------------------------------
float4 GetReg(ps_state state, uint index)
{
    return state.Regs[index & 0xF];
}

void SetReg(inout ps_state state, uint index, float4 value)
{
    if (index != PS_REGISTER_DISCARD)
        state.Regs[index & 0xF] = value;
}

// ---------------------------------------------------------------
// Input mapping
// ---------------------------------------------------------------
float4 ApplyInputMapping(uint mapping, float4 value)
{
    switch (mapping)
    {
    case PS_INPUTMAPPING_UNSIGNED_IDENTITY:
        return max(0.0, value);
    case PS_INPUTMAPPING_UNSIGNED_INVERT:
        return 1.0 - clamp(value, 0.0, 1.0);
    case PS_INPUTMAPPING_EXPAND_NORMAL:
        return 2.0 * max(0.0, value) - 1.0;
    case PS_INPUTMAPPING_EXPAND_NEGATE:
        return -2.0 * max(0.0, value) + 1.0;
    case PS_INPUTMAPPING_HALFBIAS_NORMAL:
        return max(0.0, value) - 0.5;
    case PS_INPUTMAPPING_HALFBIAS_NEGATE:
        return 0.5 - max(0.0, value);
    case PS_INPUTMAPPING_SIGNED_IDENTITY:
        return value;
    case PS_INPUTMAPPING_SIGNED_NEGATE:
    default:
        return -value;
    }
}

// ---------------------------------------------------------------
// Output mapping
// ---------------------------------------------------------------
float4 ApplyOutputMapping(uint outputFlags, float4 value)
{
    float bias = (outputFlags & PS_COMBINEROUTPUT_BIAS_BIT) ? -0.5 : 0.0;
    float scale;
    uint scaleBits = outputFlags & PS_COMBINEROUTPUT_SCALE_MASK;
    if (scaleBits == PS_COMBINEROUTPUT_SHIFTLEFT_1)
        scale = 2.0;
    else if (scaleBits == PS_COMBINEROUTPUT_SHIFTLEFT_2)
        scale = 4.0;
    else if (scaleBits == PS_COMBINEROUTPUT_SHIFTRIGHT_1)
        scale = 0.5;
    else
        scale = 1.0;

    return clamp((value + bias) * scale, -1.0, 1.0);
}

// ---------------------------------------------------------------
// Decode a combiner input byte and fetch the mapped register value
// ---------------------------------------------------------------
float4 GetCombinerInput(ps_state state, uint inputByte, bool isAlpha)
{
    uint regIndex = inputByte & 0xF;
    uint channelBit = (inputByte >> 4) & 1u;
    uint mapping = (inputByte >> 5) & 7u;

    float4 value = (float4)0;

    // Resolve C0/C1 per-stage constants
    if (regIndex == PS_REGISTER_C0)
    {
        if (state.stage < STAGE_FINAL_COMBINER)
            value = state.FlagUniqueC0 ? PSConstant0[state.stage] : PSConstant0[0];
        else
            value = PSFinalCombinerConstant[0];
    }
    else if (regIndex == PS_REGISTER_C1)
    {
        if (state.stage < STAGE_FINAL_COMBINER)
            value = state.FlagUniqueC1 ? PSConstant1[state.stage] : PSConstant1[0];
        else
            value = PSFinalCombinerConstant[1];
    }
    else if (regIndex == PS_REGISTER_FOG)
    {
        float4 fog = GetReg(state, PS_REGISTER_FOG);
        if (state.stage < STAGE_FINAL_COMBINER)
            value = float4(fog.rgb, 1.0); // Only RGB available in combiner stages
        else
            value = float4(0, 0, 0, fog.a); // Only alpha available in final combiner
    }
    else
    {
        value = GetReg(state, regIndex);
    }

    // Channel select
    if (channelBit) // PS_CHANNEL_ALPHA
        value = value.aaaa;
    else
        value = value.rgbb; // PS_CHANNEL_RGB (blue replicated to alpha for .b channel)

    return ApplyInputMapping(mapping, value);
}

// ---------------------------------------------------------------
// Color combiner stage (processes either RGB or Alpha)
// ---------------------------------------------------------------
void DoCombinerStage(inout ps_state state, bool isAlpha)
{
    // Fetch the 32-bit packed input/output registers for this stage
    uint inputs = isAlpha ? PSAlphaInputs[state.stage].x : PSRGBInputs[state.stage].x;
    uint outputs = isAlpha ? PSAlphaOutputs[state.stage].x : PSRGBOutputs[state.stage].x;

    // Decode four input bytes from the packed 32-bit value
    uint A_byte = (inputs >> 24) & 0xFF;
    uint B_byte = (inputs >> 16) & 0xFF;
    uint C_byte = (inputs >> 8) & 0xFF;
    uint D_byte = inputs & 0xFF;

    // Fetch input values
    float4 A = GetCombinerInput(state, A_byte, isAlpha);
    float4 B = GetCombinerInput(state, B_byte, isAlpha);
    float4 C = GetCombinerInput(state, C_byte, isAlpha);
    float4 D = GetCombinerInput(state, D_byte, isAlpha);

    // Decode output register indices and flags
    uint outRegCD = outputs & 0xF;
    uint outRegAB = (outputs >> 4) & 0xF;
    uint outRegABCD = (outputs >> 8) & 0xF;
    bool cdDot = ((outputs >> 12) & 1u) != 0;
    bool abDot = ((outputs >> 13) & 1u) != 0;
    bool mux = ((outputs >> 14) & 1u) != 0;
    uint outputMapping = (outputs >> 15) & 0x38; // bits [17:15] → scale+bias
    // Reconstruct output mapping byte: bit3=bias from bit15, bits[5:4]=scale from bits[17:16]
    uint outMapByte = ((outputs >> 12) & 0x38); // Extract bits 15-17 shifted to positions 3-5

    // Recalculate output mapping from raw bits
    bool biasBit = ((outputs >> 15) & 1u) != 0;
    uint scaleBits = (outputs >> 16) & 3u;
    outMapByte = (biasBit ? 0x08u : 0u) | (scaleBits << 4);

    bool cdBlueToAlpha = ((outputs >> 18) & 1u) != 0;
    bool abBlueToAlpha = ((outputs >> 19) & 1u) != 0;

    // Calculate AB and CD products (or dot products for RGB)
    float4 AB_value, CD_value;
    if (!isAlpha && abDot)
        AB_value = (float4)dot(A.rgb, B.rgb);
    else
        AB_value = A * B;

    if (!isAlpha && cdDot)
        CD_value = (float4)dot(C.rgb, D.rgb);
    else
        CD_value = C * D;

    // Calculate AB_CD (sum or mux)
    float4 ABCD_value;
    if (mux)
    {
        float r0a = GetReg(state, PS_REGISTER_R0).a;
        if (state.FlagMuxMsb)
            ABCD_value = (r0a >= 0.5) ? CD_value : AB_value;
        else
            ABCD_value = (frac(r0a * 255.0) >= 0.5) ? CD_value : AB_value; // LSB test
    }
    else
    {
        ABCD_value = AB_value + CD_value;
    }

    // Write outputs
    if (isAlpha)
    {
        // Alpha: write single component
        float abOut = ApplyOutputMapping(outMapByte, AB_value).a;
        float cdOut = ApplyOutputMapping(outMapByte, CD_value).a;
        float abcdOut = ApplyOutputMapping(outMapByte, ABCD_value).a;

        if (outRegAB != PS_REGISTER_DISCARD) {
            float4 r = GetReg(state, outRegAB);
            r.a = abOut;
            SetReg(state, outRegAB, r);
        }
        if (outRegCD != PS_REGISTER_DISCARD) {
            float4 r = GetReg(state, outRegCD);
            r.a = cdOut;
            SetReg(state, outRegCD, r);
        }
        if (outRegABCD != PS_REGISTER_DISCARD) {
            float4 r = GetReg(state, outRegABCD);
            r.a = abcdOut;
            SetReg(state, outRegABCD, r);
        }
    }
    else
    {
        // RGB: write .rgb, optionally blue-to-alpha
        float4 abMapped = ApplyOutputMapping(outMapByte, AB_value);
        float4 cdMapped = ApplyOutputMapping(outMapByte, CD_value);
        float4 abcdMapped = ApplyOutputMapping(outMapByte, ABCD_value);

        if (outRegAB != PS_REGISTER_DISCARD) {
            float4 r = GetReg(state, outRegAB);
            r.rgb = abMapped.rgb;
            if (abBlueToAlpha) r.a = abMapped.b;
            SetReg(state, outRegAB, r);
        }
        if (outRegCD != PS_REGISTER_DISCARD) {
            float4 r = GetReg(state, outRegCD);
            r.rgb = cdMapped.rgb;
            if (cdBlueToAlpha) r.a = cdMapped.b;
            SetReg(state, outRegCD, r);
        }
        if (outRegABCD != PS_REGISTER_DISCARD) {
            float4 r = GetReg(state, outRegABCD);
            r.rgb = abcdMapped.rgb;
            SetReg(state, outRegABCD, r);
        }
    }
}

// ---------------------------------------------------------------
// Final combiner
// ---------------------------------------------------------------
float4 DoFinalCombiner(inout ps_state state)
{
    uint efg = PSFinalCombinerInputsEFG.x;
    uint abcd = PSFinalCombinerInputsABCD.x;

    // If both are zero, just return R0
    if (efg == 0 && abcd == 0)
        return GetReg(state, PS_REGISTER_R0);

    // Decode E, F, G inputs and settings
    uint E_byte = (efg >> 24) & 0xFF;
    uint F_byte = (efg >> 16) & 0xFF;
    uint G_byte = (efg >> 8) & 0xFF;
    uint settings = efg & 0xFF;

    state.stage = STAGE_FINAL_COMBINER;

    float3 E_value = GetCombinerInput(state, E_byte, false).rgb;
    float3 F_value = GetCombinerInput(state, F_byte, false).rgb;
    float G_alpha = GetCombinerInput(state, G_byte, true).a;

    // EF product
    float3 EF_PROD = E_value * F_value;
    SetReg(state, PS_REGISTER_EF_PROD, float4(EF_PROD, 1.0));

    // V1+R0 sum with optional complement and clamping
    float4 R0 = GetReg(state, PS_REGISTER_R0);
    float3 R0_rgb = R0.rgb;
    if (settings & 0x20) // COMPLEMENT_R0
        R0_rgb = 1.0 - R0_rgb;

    float3 V1_rgb = GetReg(state, PS_REGISTER_V1).rgb;
    if (settings & 0x40) // COMPLEMENT_V1
        V1_rgb = 1.0 - V1_rgb;

    float3 V1R0_SUM = V1_rgb + R0_rgb;
    if (settings & 0x80) // CLAMP_SUM
        V1R0_SUM = saturate(V1R0_SUM);

    SetReg(state, PS_REGISTER_V1R0_SUM, float4(V1R0_SUM, 1.0));

    // Decode A, B, C, D inputs
    uint A_byte = (abcd >> 24) & 0xFF;
    uint B_byte = (abcd >> 16) & 0xFF;
    uint C_byte = (abcd >> 8) & 0xFF;
    uint D_byte = abcd & 0xFF;

    state.stage = STAGE_FINAL_COMBINER_ABCD;

    float4 A = GetCombinerInput(state, A_byte, false);
    float4 B = GetCombinerInput(state, B_byte, false);
    float4 C = GetCombinerInput(state, C_byte, false);
    float4 D = GetCombinerInput(state, D_byte, false);

    // Final output: RGB = A*B + (1-A)*C + D, Alpha = G
    float4 result;
    result.rgb = min(lerp(C.rgb, B.rgb, A.rgb) + D.rgb, 1.0);
    result.a = G_alpha;

    return result;
}

// ---------------------------------------------------------------
// Color sign conversion
// ---------------------------------------------------------------
float4 PerformColorSign(float4 colorSign, float4 t)
{
    // >0 means convert unsigned [0,1] to signed [-1,+1]
    // <0 means convert signed [-1,1] to unsigned [0,1]
    if (colorSign.r > 0) t.r = t.r * 2.0 - 1.0;
    if (colorSign.g > 0) t.g = t.g * 2.0 - 1.0;
    if (colorSign.b > 0) t.b = t.b * 2.0 - 1.0;
    if (colorSign.a > 0) t.a = t.a * 2.0 - 1.0;
    if (colorSign.r < 0) t.r = (t.r + 1.0) * 0.5;
    if (colorSign.g < 0) t.g = (t.g + 1.0) * 0.5;
    if (colorSign.b < 0) t.b = (t.b + 1.0) * 0.5;
    if (colorSign.a < 0) t.a = (t.a + 1.0) * 0.5;
    return t;
}

// ---------------------------------------------------------------
// Texture fetch (all texture addressing modes)
// ---------------------------------------------------------------
void FetchTexture(inout ps_state state, uint ts, float4 texCoords, uint texMode)
{
    float4 texValue = float4(0, 0, 0, 1);

    switch (texMode)
    {
    case PS_TEXTUREMODES_NONE:
        return;
    case PS_TEXTUREMODES_PROJECT2D:
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], texCoords.xy / texCoords.w);
        break;
    case PS_TEXTUREMODES_PROJECT3D:
        texValue = g_Tex3D[ts].Sample(g_Sampler[ts], texCoords.xyz / texCoords.w);
        break;
    case PS_TEXTUREMODES_CUBEMAP:
        texValue = g_TexCube[ts].Sample(g_Sampler[ts], texCoords.xyz);
        break;
    case PS_TEXTUREMODES_PASSTHRU:
        texValue = texCoords;
        break;
    case PS_TEXTUREMODES_CLIPPLANE:
        if (texCoords.x < 0 || texCoords.y < 0 || texCoords.z < 0 || texCoords.w < 0)
            discard;
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], texCoords.xy);
        break;
    case PS_TEXTUREMODES_BUMPENVMAP:
    case PS_TEXTUREMODES_BUMPENVMAP_LUM:
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], texCoords.xy);
        break;
    case PS_TEXTUREMODES_DPNDNT_AR:
    {
        float4 prev = GetReg(state, PS_REGISTER_T0 + ts - 1);
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], float2(prev.a, prev.r));
        break;
    }
    case PS_TEXTUREMODES_DPNDNT_GB:
    {
        float4 prev = GetReg(state, PS_REGISTER_T0 + ts - 1);
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], float2(prev.g, prev.b));
        break;
    }
    case PS_TEXTUREMODES_DOTPRODUCT:
    {
        float4 src = GetReg(state, PS_REGISTER_T0 + ts - 1);
        texValue = float4(dot(texCoords.xyz, src.xyz), 0, 0, 0);
        break;
    }
    case PS_TEXTUREMODES_DOT_ST:
    {
        float4 prevDot = GetReg(state, PS_REGISTER_T0 + ts - 1);
        float4 src = GetReg(state, PS_REGISTER_T0 + ts - 2);
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], float2(src.x, prevDot.x));
        break;
    }
    case PS_TEXTUREMODES_DOT_ZW:
    {
        float4 src = GetReg(state, PS_REGISTER_T0 + ts - 1);
        texValue = float4(0, 0, dot(texCoords.xyz, src.xyz), 1);
        break;
    }
    case PS_TEXTUREMODES_DOT_RFLCT_DIFF:
    {
        float4 src = GetReg(state, PS_REGISTER_T0 + ts - 1);
        float3 n = normalize(texCoords.xyz);
        float3 refl = n * dot(n, src.xyz) * 2.0 - src.xyz;
        texValue = g_TexCube[ts].Sample(g_Sampler[ts], refl);
        break;
    }
    case PS_TEXTUREMODES_DOT_RFLCT_SPEC:
    case PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST:
    {
        float3 n = normalize(texCoords.xyz);
        float3 eye = float3(0, 0, 1); // TODO: proper eye vector
        float3 refl = reflect(-eye, n);
        texValue = g_TexCube[ts].Sample(g_Sampler[ts], refl);
        break;
    }
    case PS_TEXTUREMODES_DOT_STR_3D:
    {
        float4 dot2 = GetReg(state, PS_REGISTER_T0 + ts - 2);
        float4 dot1 = GetReg(state, PS_REGISTER_T0 + ts - 1);
        float4 src = GetReg(state, PS_REGISTER_T0 + ts);
        float dv = dot(texCoords.xyz, src.xyz);
        texValue = g_Tex3D[ts].Sample(g_Sampler[ts], float3(dot2.x, dot1.x, dv));
        break;
    }
    case PS_TEXTUREMODES_DOT_STR_CUBE:
    {
        float4 dot2 = GetReg(state, PS_REGISTER_T0 + ts - 2);
        float4 dot1 = GetReg(state, PS_REGISTER_T0 + ts - 1);
        float4 src = GetReg(state, PS_REGISTER_T0 + ts);
        float dv = dot(texCoords.xyz, src.xyz);
        texValue = g_TexCube[ts].Sample(g_Sampler[ts], float3(dot2.x, dot1.x, dv));
        break;
    }
    case PS_TEXTUREMODES_BRDF:
    default:
        texValue = g_Tex2D[ts].Sample(g_Sampler[ts], texCoords.xy);
        break;
    }

    // Apply color sign conversion
    texValue = PerformColorSign(PSColorSign[ts], texValue);

    SetReg(state, PS_REGISTER_T0 + ts, texValue);
}

// ---------------------------------------------------------------
// Pixel shader input structure (matches VS_OUTPUT from vertex shaders)
// ---------------------------------------------------------------
struct PS_INPUT
{
    float4 Position : SV_Position;
    float4 oD0 : COLOR0;     // Diffuse
    float4 oD1 : COLOR1;     // Specular
    float4 oB0 : COLOR2;     // Back diffuse
    float4 oB1 : COLOR3;     // Back specular
    float  oFog : FOG;
    float  oPts : PSIZE;
    float4 oT0 : TEXCOORD0;
    float4 oT1 : TEXCOORD1;
    float4 oT2 : TEXCOORD2;
    float4 oT3 : TEXCOORD3;
};

// ---------------------------------------------------------------
// Main pixel shader entry point
// ---------------------------------------------------------------
float4 main(PS_INPUT input) : SV_Target
{
    ps_state state = (ps_state)0;

    // Decode PSTextureModes (packed 5 bits per stage)
    uint tmRaw = PSTextureModes.x;
    uint texModes[4];
    texModes[0] = ExtractBits(tmRaw, 0, 5);
    texModes[1] = ExtractBits(tmRaw, 5, 5);
    texModes[2] = ExtractBits(tmRaw, 10, 5);
    texModes[3] = ExtractBits(tmRaw, 15, 5);

    // Fetch textures
    FetchTexture(state, 0, input.oT0, texModes[0]);
    FetchTexture(state, 1, input.oT1, texModes[1]);
    FetchTexture(state, 2, input.oT2, texModes[2]);
    FetchTexture(state, 3, input.oT3, texModes[3]);

    // Initialize register state
    SetReg(state, PS_REGISTER_FOG, float4(input.oB0.rgb, input.oFog));
    SetReg(state, PS_REGISTER_V0, input.oD0);
    SetReg(state, PS_REGISTER_V1, input.oD1);
    SetReg(state, PS_REGISTER_R0, float4(1, 1, 1, GetReg(state, PS_REGISTER_T0).a));

    // Decode combiner count and flags
    uint ccRaw = PSCombinerCount.x;
    uint numStages = ccRaw & 0xF;
    state.FlagMuxMsb = ((ccRaw >> 8) & 1u) != 0;
    state.FlagUniqueC0 = ((ccRaw >> 12) & 1u) != 0;
    state.FlagUniqueC1 = ((ccRaw >> 16) & 1u) != 0;

    numStages = clamp(numStages, 1u, 8u);

    // Execute combiner stages
    [loop] for (uint i = 0; i < numStages; i++)
    {
        state.stage = i;
        DoCombinerStage(state, false); // RGB
        DoCombinerStage(state, true);  // Alpha
    }

    // Final combiner
    return DoFinalCombiner(state);
}
