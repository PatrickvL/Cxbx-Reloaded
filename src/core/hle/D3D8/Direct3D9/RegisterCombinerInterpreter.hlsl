//
// Cxbx-Reloaded Xbox Pixel Shader Interpreter, written in Direct3D 9.0 Shader Model 3.0 HLSL
//
// This shader interprets Xbox Pixel Shaders, as declared in the Xbox Direct3D__RenderState registers.
// These registers must be fed into this shader though host SetPixelShaderConstantF(),
// so that they can be accessed from the corresponding constant registers.
//
// Since Shader Model doesn't support bitwise operations nor native integers,
// the below code expects 32 bit input data to be offered split up in four parts.
// The .r part must contain bits  0.. 7 of any input, expressed as a float in the range [0-255].
// The .g part must contain bits  8..15 of any input, expressed as a float in the range [0-255].
// The .b part must contain bits 16..23 of any input, expressed as a float in the range [0-255].
// The .a part must contain bits 24..31 of any input, expressed as a float in the range [0-255].
// Using this representation, it's possible to extract bits from floats via fmod(),
// as described in http://theinstructionlimit.com/encoding-boolean-flags-into-a-float-in-hlsl
// Also, see https://en.wikipedia.org/wiki/Single-precision_floating-point_format
// Test fmod() online via https://www.functions-online.com/fmod.html
// Note, that this only applies to bitfield registers; color constant registers are expressed
// as four floats, one for each color channel, in the range [0.0-1.0].

//
// Utility functions
//

#define SUPPORTS_INDIRECT_INDEX // When set, register value accesses indexes an array instead of accessing separate fields // [125 slots] when unset?
//#define ASSERT_VALID_ACCESSES // When set, checks for valid accesses are compiled (doesn't result in additional code)
//#define AVOID_INVALID_ACCESSES // When set, invalid accesses are avoided (although assertions may have failed already)

typedef uint byte_t; // Indicates a float in the range [0.0-255.0] (always positive with zero fraction)

#ifdef ASSERT_VALID_ACCESSES
void assert(bool condition)
{
	// TODO
}
#endif

// Input float is read from one of the D3DRS_* 4 byte floats, which MUST already be in range [0.0-255.0]!
// Given this type of input, do as little as possible to get this compiling without errors.
byte_t float_to_byte(uniform float value) // TODO : Test this thoroughly, aiming for the fastest correct implementation that compiles without problems
{
#ifdef ASSERT_VALID_ACCESSES
    assert(value >= 0.0f);
    assert(value <= 255.0f);
#endif
    //byte_t byte_result = value; // error X3548: in ps_3_0 uints can only be used with known-positive values, use int if possible
    //byte_t byte_result = fmod(value, 256u); // error X3548: in ps_3_0 uints can only be used with known-positive values, use int if possible
    //byte_t byte_result = fmod(abs(value), 256u); // error X3548: in ps_3_0 uints can only be used with known-positive values, use int if possible

    byte_t byte_result = max(value, 255u); // [623 slots] Compiles in Visual Studio, not in AMD GPU ShaderAnalyzer
    //byte_t byte_result = max(abs(value), 255); // [624 slots] Compiles in Visual Studio, not in AMD GPU ShaderAnalyzer
    //byte_t byte_result = abs(value); // [626 slots] Compiles in Visual Studio, not in AMD GPU ShaderAnalyzer
    //byte_t byte_result = abs(fmod(value, 256)); // [652 slots] Probably most correct. Compiles in Visual Studio, not in AMD GPU ShaderAnalyzer

    //byte_t byte_result = max(fmod(value, 256u), 255u); // [13 slots] Visual Studio incorrectly optimizes nearly everything away
#ifdef ASSERT_VALID_ACCESSES
    assert(byte_result >= 0u);
    assert(byte_result <= 255u);
#endif
    return byte_result;
}

// byte functions
byte_t get_byte_0_from_float4(uniform float4 value)
{
    return float_to_byte(value.r);
}

byte_t get_byte_1_from_float4(uniform float4 value)
{
    return float_to_byte(value.g);
}

byte_t get_byte_2_from_float4(uniform float4 value)
{
    return float_to_byte(value.b);
}

byte_t get_byte_3_from_float4(uniform float4 value)
{
    return float_to_byte(value.a);
}

// bit functions
bool get_bit_0_from_byte(uniform byte_t value)
{
    return fmod(value, 2.0f) >= 1.0f;
}

#if 0 // unused for now
bool get_bit_1_from_byte(uniform byte_t value)
{
    return fmod(value, 4.0f) >= 2.0f;
}
#endif

bool get_bit_2_from_byte(uniform byte_t value)
{
    return fmod(value, 8.0f) >= 4.0f;
}

bool get_bit_3_from_byte(uniform byte_t value)
{
    return fmod(value, 16.0f) >= 8.0f;
}

bool get_bit_4_from_byte(uniform byte_t value)
{
    return fmod(value, 32.0f) >= 16.0f;
}

bool get_bit_5_from_byte(uniform byte_t value)
{
    return fmod(value, 64.0f) >= 32.0f;
}

bool get_bit_6_from_byte(uniform byte_t value)
{
    return fmod(value, 128.0f) >= 64.0f;
}

bool get_bit_7_from_byte(uniform byte_t value)
{
    return fmod(value, 256.0f) >= 128.0f;
}

// nibble functions
byte_t get_nibble_0_from_byte(uniform byte_t value)
{
    byte_t nibble_0_result = value;
    //nibble_0_result = abs(nibble_0_result);
    nibble_0_result = fmod(nibble_0_result, 16u);
    //nibble_0_result = max(nibble_0_result, 15u);
#ifdef ASSERT_VALID_ACCESSES
    assert(nibble_0_result >= 0u);
    assert(nibble_0_result <= 15u);
#endif
    return nibble_0_result;
}

byte_t get_nibble_1_from_byte(uniform byte_t value)
{
    byte_t nibble_1_result = value;
    //nibble_1_result = abs(nibble_1_result);
    nibble_1_result = nibble_1_result / 16u;
    //nibble_1_result = max(nibble_1_result, 15u);
#ifdef ASSERT_VALID_ACCESSES
    assert(nibble_1_result >= 0u);
    assert(nibble_1_result <= 15u);
#endif
    return nibble_1_result;
}

//
// Xbox Pixel Shader definition
//

/*---------------------------------------------------------------------------------*/
/*  Color combiners - The following members of the D3DPixelShaderDef structure     */
/*  define the state for the eight stages of color combiners:                      */
/*      PSCombinerCount - Number of stages                                         */
/*      PSAlphaInputs[8] - Inputs for alpha portion of each stage                  */
/*      PSRGBInputs[8] - Inputs for RGB portion of each stage                      */
/*      PSConstant0[8] - Constant 0 for each stage                                 */
/*      PSConstant1[8] - Constant 1 for each stage                                 */
/*      PSFinalCombinerConstant0 - Constant 0 for final combiner                   */
/*      PSFinalCombinerConstant1 - Constant 1 for final combiner                   */
/*      PSAlphaOutputs[8] - Outputs for alpha portion of each stage                */
/*      PSRGBOutputs[8] - Outputs for RGB portion of each stage                    */
/*---------------------------------------------------------------------------------*/

// =========================================================================================================
// PSCombinerCount
// --------.--------.--------.----xxxx // number of combiners (1-8)
// --------.--------.-------x.-------- // mux bit (0= LSB, 1= MSB)
// --------.--------.---x----.-------- // separate C0
// --------.-------x.--------.-------- // separate C1

// #define PS_COMBINERCOUNT(count, flags) (((flags)<<8)|(count))
// count is 1-8, flags contains one or more values from PS_COMBINERCOUNTFLAGS

// enum PS_COMBINERCOUNTFLAGS
//    PS_COMBINERCOUNT_MUX_LSB = 0x0000L, // mux on r0.a lsb
//    PS_COMBINERCOUNT_MUX_MSB = 0x0001L, // mux on r0.a msb
//    PS_COMBINERCOUNT_SAME_C0 = 0x0000L, // c0 same in each stage
//    PS_COMBINERCOUNT_UNIQUE_C0 = 0x0010L, // c0 unique in each stage
//    PS_COMBINERCOUNT_SAME_C1 = 0x0000L, // c1 same in each stage
//    PS_COMBINERCOUNT_UNIQUE_C1 = 0x0100L // c1 unique in each stage

// =========================================================================================================
// PSRGBInputs[0-7]
// PSAlphaInputs[0-7]
// PSFinalCombinerInputsABCD
// PSFinalCombinerInputsEFG
// --------.--------.--------.----xxxx // D register
// --------.--------.--------.---x---- // D channel (0= RGB/BLUE, 1= ALPHA)
// --------.--------.--------.xxx----- // D input mapping
// --------.--------.----xxxx.-------- // C register
// --------.--------.---x----.-------- // C channel (0= RGB/BLUE, 1= ALPHA)
// --------.--------.xxx-----.-------- // C input mapping
// --------.----xxxx.--------.-------- // B register
// --------.---x----.--------.-------- // B channel (0= RGB/BLUE, 1= ALPHA)
// --------.xxx-----.--------.-------- // B input mapping
// ----xxxx.--------.--------.-------- // A register
// ---x----.--------.--------.-------- // A channel (0= RGB/BLUE, 1= ALPHA)
// xxx-----.--------.--------.-------- // A input mapping

// examples:
//
// shader.PSRGBInputs[3]= PS_COMBINERINPUTS(
//     PS_REGISTER_T0 | PS_INPUTMAPPING_EXPAND_NORMAL     | PS_CHANNEL_RGB,
//     PS_REGISTER_C0 | PS_INPUTMAPPING_UNSIGNED_IDENTITY | PS_CHANNEL_ALPHA,
//     PS_REGISTER_ZERO,
//     PS_REGISTER_ZERO);
//
// shader.PSFinalCombinerInputsABCD= PS_COMBINERINPUTS(
//     PS_REGISTER_T0     | PS_INPUTMAPPING_UNSIGNED_IDENTITY | PS_CHANNEL_ALPHA,
//     PS_REGISTER_ZERO   | PS_INPUTMAPPING_EXPAND_NORMAL     | PS_CHANNEL_RGB,
//     PS_REGISTER_EFPROD | PS_INPUTMAPPING_UNSIGNED_INVERT   | PS_CHANNEL_RGB,
//     PS_REGISTER_ZERO);
//
// PS_FINALCOMBINERSETTING is set in 4th field of PSFinalCombinerInputsEFG with PS_COMBINERINPUTS
// example:
//
// shader.PSFinalCombinerInputsEFG= PS_COMBINERINPUTS(
//     PS_REGISTER_R0 | PS_INPUTMAPPING_UNSIGNED_IDENTITY | PS_CHANNEL_RGB,
//     PS_REGISTER_R1 | PS_INPUTMAPPING_UNSIGNED_IDENTITY | PS_CHANNEL_RGB,
//     PS_REGISTER_R1 | PS_INPUTMAPPING_UNSIGNED_IDENTITY | PS_CHANNEL_BLUE,
//    PS_FINALCOMBINERSETTING_CLAMP_SUM | PS_FINALCOMBINERSETTING_COMPLEMENT_R0);

//#define PS_COMBINERINPUTS(a,b,c,d) (((a)<<24)|((b)<<16)|((c)<<8)|(d))
// For PSFinalCombinerInputsEFG,
//     a,b,c contain a value from PS_REGISTER, PS_CHANNEL, and PS_INPUTMAPPING for input E,F, and G
//     d contains values from PS_FINALCOMBINERSETTING
// For all other inputs,
//     a,b,c,d each contain a value from PS_REGISTER, PS_CHANNEL, and PS_INPUTMAPPING

// The input can have the following mappings applied :
//
// PS_INPUTMAPPING_UNSIGNED_IDENTITY : y =  1*max(0,x) - 0.0  ?  OR  abs(x)                        ?
// PS_INPUTMAPPING_UNSIGNED_INVERT   : y = -1*max(0,x) + 1.0  ?  OR  1 - x   MAYBE EVEN  1-abs(x)  ?
// PS_INPUTMAPPING_EXPAND_NORMAL     : y =  2*max(0,x) - 1.0
// PS_INPUTMAPPING_EXPAND_NEGATE     : y = -2*max(0,x) + 1.0
// PS_INPUTMAPPING_HALFBIAS_NORMAL   : y =  1*max(0,x) - 0.5
// PS_INPUTMAPPING_HALFBIAS_NEGATE   : y = -1*max(0,x) + 0.5
// PS_INPUTMAPPING_SIGNED_IDENTITY   : y =  1*      x  - 0.0
// PS_INPUTMAPPING_SIGNED_NEGATE     : y = -1*      x  + 0.0
//
// (Note : I don't know for sure if the max() operation mentioned above is indeed what happens,
// as there's no further documentation available on this. Native Direct3D can clamp with the
// '_sat' instruction modifier, but that's not really the same as these Xbox1 input mappings.)
//
// When the input register is PS_ZERO, the above mappings result in the following constants:
//
// PS_REGISTER_NEGATIVE_ONE      (PS_INPUTMAPPING_EXPAND_NORMAL on zero)   : y = -1.0
// PS_REGISTER_NEGATIVE_ONE_HALF (PS_INPUTMAPPING_HALFBIAS_NORMAL on zero) : y = -0.5
// PS_REGISTER_ZERO itself                                                 : y =  0.0
// PS_REGISTER_ONE_HALF          (PS_INPUTMAPPING_HALFBIAS_NEGATE on zero) : y =  0.5
// PS_REGISTER_ONE               (PS_INPUTMAPPING_UNSIGNED_INVERT on zero) : y =  1.0
// (Note : It has no define, but PS_INPUTMAPPING_EXPAND_NEGATE on zero results in ONE too!)
static const byte_t PS_INPUTMAPPING_UNSIGNED_IDENTITY = 0x00L; //     max(0,x)     OK for final combiner: y = abs(x)
static const byte_t PS_INPUTMAPPING_UNSIGNED_INVERT   = 0x20L; // 1 - max(0,x)     OK for final combiner: y = 1 - x
static const byte_t PS_INPUTMAPPING_EXPAND_NORMAL     = 0x40L; // 2*max(0,x) - 1   invalid for final combiner
static const byte_t PS_INPUTMAPPING_EXPAND_NEGATE     = 0x60L; // 1 - 2*max(0,x)   invalid for final combiner
static const byte_t PS_INPUTMAPPING_HALFBIAS_NORMAL   = 0x80L; // max(0,x) - 1/2   invalid for final combiner
static const byte_t PS_INPUTMAPPING_HALFBIAS_NEGATE   = 0xa0L; // 1/2 - max(0,x)   invalid for final combiner
static const byte_t PS_INPUTMAPPING_SIGNED_IDENTITY   = 0xc0L; //  x               invalid for final combiner
static const byte_t PS_INPUTMAPPING_SIGNED_NEGATE     = 0xe0L; // -x               invalid for final combiner

// Pixel Shader Registers
static const byte_t PS_REGISTER_ZERO = 0x00L; // r
static const byte_t PS_REGISTER_DISCARD = 0x00L; // w
static const byte_t PS_REGISTER_C0 = 0x01L; // r
static const byte_t PS_REGISTER_C1 = 0x02L; // r
static const byte_t PS_REGISTER_FOG = 0x03L; // r
static const byte_t PS_REGISTER_V0 = 0x04L; // r/w
static const byte_t PS_REGISTER_V1 = 0x05L; // r/w
// Filling the PS_REGISTER_* gap here with V2 and V3 is just guesswork,
// they could be unsupported, read-only, special-purpose, or whatever else.
// Needs investigation.
//static const byte_t PS_REGISTER_V2 =                0x06L; // unknown r/w state
//static const byte_t PS_REGISTER_V3 =                0x07L; // unknown r/w state
static const byte_t PS_REGISTER_T0 = 0x08L; // r/w
static const byte_t PS_REGISTER_T1 = 0x09L; // r/w
static const byte_t PS_REGISTER_T2 = 0x0aL; // r/w
static const byte_t PS_REGISTER_T3 = 0x0bL; // r/w
static const byte_t PS_REGISTER_R0 = 0x0cL; // r/w
static const byte_t PS_REGISTER_R1 = 0x0dL; // r/w
static const byte_t PS_REGISTER_V1R0_SUM = 0x0eL; // r
static const byte_t PS_REGISTER_EF_PROD = 0x0fL; // r
//static const byte_t PS_REGISTER_ONE = PS_REGISTER_ZERO | PS_INPUTMAPPING_UNSIGNED_INVERT; // 0x20 OK for final combiner
//static const byte_t PS_REGISTER_NEGATIVE_ONE = PS_REGISTER_ZERO | PS_INPUTMAPPING_EXPAND_NORMAL; // 0x40 invalid for final combiner
//static const byte_t PS_REGISTER_ONE_HALF = PS_REGISTER_ZERO | PS_INPUTMAPPING_HALFBIAS_NEGATE; // 0xa0 invalid for final combiner
//static const byte_t PS_REGISTER_NEGATIVE_ONE_HALF = PS_REGISTER_ZERO | PS_INPUTMAPPING_HALFBIAS_NORMAL; // 0x80 invalid for final combiner

static const byte_t PS_CHANNEL_RGB = 0x00; // used as RGB source
static const byte_t PS_CHANNEL_BLUE = 0x00; // used as ALPHA source
static const byte_t PS_CHANNEL_ALPHA = 0x10; // used as RGB or ALPHA source

static const byte_t PS_FINALCOMBINERSETTING_CLAMP_SUM = 0x80; // V1+R0 sum clamped to [0;1]
static const byte_t PS_FINALCOMBINERSETTING_COMPLEMENT_V1 = 0x40; // unsigned invert mapping  (1 - v1) is used as an input to the sum rather than v1
static const byte_t PS_FINALCOMBINERSETTING_COMPLEMENT_R0 = 0x20; // unsigned invert mapping  (1 - r0) is used as an input to the sum rather than r0

// =========================================================================================================
// PSRGBOutputs[0-7]
// PSAlphaOutputs[0-7]
// --------.--------.--------.----xxxx // CD register
// --------.--------.--------.xxxx---- // AB register
// --------.--------.----xxxx.-------- // SUM register
// --------.--------.---x----.-------- // CD output (0= multiply, 1= dot product)
// --------.--------.--x-----.-------- // AB output (0= multiply, 1= dot product)
// --------.--------.-x------.-------- // AB_CD mux/sum select (0= sum, 1= mux)
// --------.------xx.x-------.-------- // Output mapping
// --------.-----x--.--------.-------- // CD blue to alpha
// --------.----x---.--------.-------- // AB blue to alpha
//
//#define PS_COMBINEROUTPUTS(ab,cd,mux_sum,flags) (((flags)<<12)|((mux_sum)<<8)|((ab)<<4)|(cd))
// ab,cd,mux_sum contain a value from PS_REGISTER
// flags contains values from PS_COMBINEROUTPUT
static const byte_t PS_COMBINEROUTPUT_CD_MULTIPLY = 0x00L;
static const byte_t PS_COMBINEROUTPUT_CD_DOT_PRODUCT = 0x01L; // RGB only
static const byte_t PS_COMBINEROUTPUT_AB_MULTIPLY = 0x00L;
static const byte_t PS_COMBINEROUTPUT_AB_DOT_PRODUCT = 0x02L; // RGB only
static const byte_t PS_COMBINEROUTPUT_AB_CD_SUM = 0x00L; // 3rd output is AB+CD
static const byte_t PS_COMBINEROUTPUT_AB_CD_MUX = 0x04L; // 3rd output is MUX(AB;CD) based on R0.a
static const byte_t PS_COMBINEROUTPUT_IDENTITY = 0x00L; // y = x
static const byte_t PS_COMBINEROUTPUT_BIAS = 0x08L; // y = x - 0.5
static const byte_t PS_COMBINEROUTPUT_SHIFTLEFT_1 = 0x10L; // y = x*2
static const byte_t PS_COMBINEROUTPUT_SHIFTLEFT_1_BIAS = 0x18L; // y = (x - 0.5)*2
static const byte_t PS_COMBINEROUTPUT_SHIFTLEFT_2 = 0x20L; // y = x*4
static const byte_t PS_COMBINEROUTPUT_SHIFTLEFT_2_BIAS  = 0x28L; // y = (x - 0.5)*4 // TODO : Unknown if this actually exists, needs verification
static const byte_t PS_COMBINEROUTPUT_SHIFTRIGHT_1 = 0x30L; // y = x/2
static const byte_t PS_COMBINEROUTPUT_SHIFTRIGHT_1_BIAS = 0x38L; // y = (x - 0.5)/2 // TODO : Unknown if this actually exists, needs verification
static const byte_t PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA = 0x40L; // RGB only
static const byte_t PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA = 0x80L; // RGB only
// AB_CD register output must be DISCARD if either AB_DOT_PRODUCT or CD_DOT_PRODUCT are set

// =========================================================================================================
// PSC0Mapping
// PSC1Mapping
// --------.--------.--------.----xxxx // offset of D3D constant for stage 0
// --------.--------.--------.xxxx---- // offset of D3D constant for stage 1
// --------.--------.----xxxx.-------- // offset of D3D constant for stage 2
// --------.--------.xxxx----.-------- // offset of D3D constant for stage 3
// --------.----xxxx.--------.-------- // offset of D3D constant for stage 4
// --------.xxxx----.--------.-------- // offset of D3D constant for stage 5
// ----xxxx.--------.--------.-------- // offset of D3D constant for stage 6
// xxxx----.--------.--------.-------- // offset of D3D constant for stage 7

//#define PS_CONSTANTMAPPING(s0,s1,s2,s3,s4,s5,s6,s7) \
//     (((DWORD)(s0)&0xf)<< 0) | (((DWORD)(s1)&0xf)<< 4) | \
//     (((DWORD)(s2)&0xf)<< 8) | (((DWORD)(s3)&0xf)<<12) | \
//     (((DWORD)(s4)&0xf)<<16) | (((DWORD)(s5)&0xf)<<20) | \
//     (((DWORD)(s6)&0xf)<<24) | (((DWORD)(s7)&0xf)<<28)
// s0-s7 contain the offset of the D3D constant that corresponds to the
// c0 or c1 constant in stages 0 through 7.  These mappings are only used in
// SetPixelShaderConstant().

// =========================================================================================================
// PSFinalCombinerConstants
// --------.--------.--------.----xxxx // offset of D3D constant for C0
// --------.--------.--------.xxxx---- // offset of D3D constant for C1
// --------.--------.-------x.-------- // Adjust texture flag

//#define PS_FINALCOMBINERCONSTANTS(c0,c1,flags) (((DWORD)(flags) << 8) | ((DWORD)(c0)&0xf)<< 0) | (((DWORD)(c1)&0xf)<< 4)
// c0 and c1 contain the offset of the D3D constant that corresponds to the
// constants in the final combiner.  These mappings are only used in
// SetPixelShaderConstant().  Flags contains values from PS_GLOBALFLAGS

// enum PS_GLOBALFLAGS
// if this flag is set, the texture mode for each texture stage is adjusted as follows:
//     if set texture is a cubemap,
//         change PS_TEXTUREMODES_PROJECT2D to PS_TEXTUREMODES_CUBEMAP
//         change PS_TEXTUREMODES_PROJECT3D to PS_TEXTUREMODES_CUBEMAP
//         change PS_TEXTUREMODES_DOT_STR_3D to PS_TEXTUREMODES_DOT_STR_CUBE
//     if set texture is a volume texture,
//         change PS_TEXTUREMODES_PROJECT2D to PS_TEXTUREMODES_PROJECT3D
//         change PS_TEXTUREMODES_CUBEMAP to PS_TEXTUREMODES_PROJECT3D
//         change PS_TEXTUREMODES_DOT_STR_CUBE to PS_TEXTUREMODES_DOT_STR_3D
//     if set texture is neither cubemap or volume texture,
//         change PS_TEXTUREMODES_PROJECT3D to PS_TEXTUREMODES_PROJECT2D
//         change PS_TEXTUREMODES_CUBEMAP to PS_TEXTUREMODES_PROJECT2D
//
//    PS_GLOBALFLAGS_NO_TEXMODE_ADJUST = 0x0000L, // don"t adjust texture modes
//    PS_GLOBALFLAGS_TEXMODE_ADJUST = 0x0001L, // adjust texture modes according to set texture

// Cxbx constants
static const byte_t MAX_COMBINER_STAGE_COUNT = 8u;
static const byte_t STAGE_FINAL_COMBINER = MAX_COMBINER_STAGE_COUNT; // Used for E,F,G inputs
static const byte_t STAGE_FINAL_COMBINER_ABCD = STAGE_FINAL_COMBINER + 1u;
static const float DEFAULT_RGB = 0.0f;
static const float DEFAULT_ALPHA = 1.0f;

#ifdef INPUT_BUFFER // TODO : Expand this for Shader Model 4 and up :
// Xbox Pixel Shader definition (a byte_t array) D3D__RenderState member indices :
#else // Shader Model 3 :

//
// Xbox Pixel Shader definition (32 bit dword values, taken from "D3D__RenderState", split in float4 per byte where needed)
//

// Since shader model 3.0 seems to have no input-buffer feature, and using texture-sampling might interpolate values, use host pixel shader constants to access the Xbox renderstate input :
nointerpolation float4 D3DRS_PSALPHAINPUTS[8]           : register(c0); // 4 byte floats
nointerpolation float4 D3DRS_PSFINALCOMBINERINPUTSABCD  : register(c8); // 4 byte floats
nointerpolation float4 D3DRS_PSFINALCOMBINERINPUTSEFG   : register(c9); // 4 byte floats
nointerpolation float4 D3DRS_PSCONSTANT0[8]             : register(c10);
nointerpolation float4 D3DRS_PSCONSTANT1[8]             : register(c18);
nointerpolation float4 D3DRS_PSALPHAOUTPUTS[8]          : register(c26); // 4 byte floats
nointerpolation float4 D3DRS_PSRGBINPUTS[8]             : register(c34); // 4 byte floats
//nointerpolation float4 D3DRS_PSCOMPAREMODE              : register(c42); // 4 byte floats - TODO unimplemented
nointerpolation float4 D3DRS_PSFINALCOMBINERCONSTANT[2] : register(c43);
nointerpolation float4 D3DRS_PSRGBOUTPUTS[8]            : register(c45); // 4 byte floats
nointerpolation float4 D3DRS_PSCOMBINERCOUNT            : register(c53); // 4 byte floats
//nointerpolation float4 D3DRS_PSTEXTUREMODES             : register(c54); // 4 byte floats - CxbxUpdateActivePixelShader_HLSL() reads D3DRS_PSTEXTUREMODES (at render state offset 136), and puts it here in c54, which is named D3DRS_PS_RESERVED in Xbox D3D
//nointerpolation float4 D3DRS_PSDOTMAPPING               : register(c55); // 4 byte floats - TODO unimplemented
//nointerpolation float4 D3DRS_PSINPUTTEXTURE             : register(c56); // 4 byte floats - TODO unimplemented

#endif

//
// Our Xbox Pixel Shader Interpreter state
//
typedef
struct _ps_state
{
	// Xbox NV2A pixel shader register values, in order of PS_REGISTER_* :
#ifndef SUPPORTS_INDIRECT_INDEX
	// nointerpolation float4 ZERO;
    // nointerpolation float4 C[2];    // Constant registers, read-only
    nointerpolation float4 FOG;      // Note, FOG ALPHA is only available in final combiner
    nointerpolation float4 V[2];     // Vertex color registers, read/write (V[0] = diffuse, V[1] = specular
	// nointerpolation float4 Unknown[2];
    nointerpolation float4 T[4];     // Texture registers, initially sampled or (0,0,0,1), after texture-addressing they're read/write just like temporary registers
    nointerpolation float4 R[2];     // Temporary registers, read/write (R[0].a (R0_ALPHA) is initialized to T[0].s (T0_ALPHA) in stage 0); Final result is in R[0]
    nointerpolation float3 V1R0_SUM; // Note : V1R0_SUM and EF_PROD are only available in final combiner (A,B,C,D inputs only)
    nointerpolation float3 EF_PROD;  // Note : V1R0_SUM_ALPHA and EF_PROD_ALPHA are not available, hence the float3 type here
#else
    nointerpolation float4 RegisterValues[16]; // Room for 16 registers : PS_REGISTER_ZERO (0) upto PS_REGISTER_EF_PROD (15)
#endif
	// Flags, used in get_input_register_as_float4() :
    nointerpolation byte_t stage; // Currently active stage (0 upto 7, 8 and 9 are for final combiner)
	nointerpolation bool FlagMuxMsb;   // Mux on r0.a lsb or msb, 0 = PS_COMBINERCOUNT_MUX_LSB, 1 = PS_COMBINERCOUNT_MUX_MSB
	nointerpolation bool FlagUniqueC0; // C[0] unique in each stage, 0 = PS_COMBINERCOUNT_SAME_C0, 1 = PS_COMBINERCOUNT_UNIQUE_C0
    nointerpolation bool FlagUniqueC1; // C[1] unique in each stage, 0 = PS_COMBINERCOUNT_SAME_C1, 1 = PS_COMBINERCOUNT_UNIQUE_C1

} ps_state;

byte_t mask_register(uniform byte_t value)
{
    byte_t mask_register_result = value;
    mask_register_result = fmod(mask_register_result, PS_REGISTER_EF_PROD + 1u);
    mask_register_result = max(mask_register_result, PS_REGISTER_EF_PROD); // max() avoids error X3500: array reference cannot be used as an l-value; not natively addressable
#ifdef ASSERT_VALID_ACCESSES
    assert(mask_register_result >= PS_REGISTER_ZERO);
    assert(mask_register_result <= PS_REGISTER_EF_PROD);
#endif
    return mask_register_result;
}

byte_t mask_inputmapping(uniform byte_t value)
{
    byte_t mask_inputmapping_result = value;
    mask_inputmapping_result = mask_inputmapping_result / 2u; // remove least significant bit
    mask_inputmapping_result = fmod(mask_inputmapping_result, 8u); // == (PS_INPUTMAPPING_SIGNED_NEGATE/2) + 1;
    mask_inputmapping_result = mask_inputmapping_result * 2u; // shift bits back to their original place
#ifdef ASSERT_VALID_ACCESSES
    assert(mask_inputmapping_result >= PS_INPUTMAPPING_UNSIGNED_IDENTITY);
    assert(mask_inputmapping_result <= PS_INPUTMAPPING_SIGNED_NEGATE);
#endif
    return mask_inputmapping_result;
}

byte_t mask_outputmapping(uniform byte_t byte_1, uniform byte_t byte_2)
{
    byte_t mask_outputmapping_result = byte_2;
    mask_outputmapping_result = fmod(mask_outputmapping_result, 4u); // only use the lowest 2 bits of byte 2
    mask_outputmapping_result = mask_outputmapping_result * PS_COMBINEROUTPUT_SHIFTLEFT_1; // shift them left 4 bits
    if (get_bit_7_from_byte(byte_1)) // when bit 7 of byte 1 is set
        mask_outputmapping_result = mask_outputmapping_result + PS_COMBINEROUTPUT_BIAS; // include the BIAS bit in the result

#ifdef ASSERT_VALID_ACCESSES
	assert(mask_outputmapping_result >= PS_COMBINEROUTPUT_IDENTITY);
    assert(mask_outputmapping_result <= PS_COMBINEROUTPUT_SHIFTRIGHT_1_BIAS);
#endif
    return mask_outputmapping_result;
}

// Read without any post-processing nor validation
float4 get_plain_register_as_float4(uniform ps_state state, uniform byte_t register_index)
{
#ifndef SUPPORTS_INDIRECT_INDEX
    switch (register_index)
    {
		case PS_REGISTER_FOG:
			return state.FOG;
        case PS_REGISTER_V0:
            return state.V[0];
        case PS_REGISTER_V1:
            return state.V[1];
        case PS_REGISTER_T0:
            return state.T[0];
        case PS_REGISTER_T1:
            return state.T[1];
        case PS_REGISTER_T2:
            return state.T[2];
        case PS_REGISTER_T3:
            return state.T[3];
		case PS_REGISTER_R0:
            return state.R[0];
        case PS_REGISTER_R1:
            return state.R[1];
        case PS_REGISTER_V1R0_SUM:
            return float4(state.V1R0_SUM, DEFAULT_ALPHA);
        case PS_REGISTER_EF_PROD:
            return float4(state.EF_PROD, DEFAULT_ALPHA);
		default: // This also handles PS_REGISTER_ZERO
            return 0; // Avoids error X3507: 'get_plain_register_as_float4': Not all control paths return a value
    }
#else
    return state.RegisterValues[register_index];
#endif
}

// Write without any pre-processing and minimal validation
void set_plain_register_as_float4(inout ps_state state, uniform byte_t register_index, uniform float4 value)
{
#ifndef SUPPORTS_INDIRECT_INDEX
    switch (register_index)
    {
		case PS_REGISTER_FOG:
		    state.FOG = value;
            break;
        case PS_REGISTER_V0:
            state.V[0] = value;
            break;
        case PS_REGISTER_V1:
            state.V[1] = value;
            break;
        case PS_REGISTER_T0:
            state.T[0] = value;
            break;
        case PS_REGISTER_T1:
            state.T[1] = value;
            break;
        case PS_REGISTER_T2:
            state.T[2] = value;
            break;
        case PS_REGISTER_T3:
            state.T[3] = value;
            break;
		case PS_REGISTER_R0:
            state.R[0] = value;
            break;
        case PS_REGISTER_R1:
            state.R[1] = value;
            break;
        case PS_REGISTER_V1R0_SUM:
			state.V1R0_SUM = value.rgb;
			break;
        case PS_REGISTER_EF_PROD:
			state.EF_PROD = value.rgb;
			break;
    }
#else
	// Avoid writing to PS_REGISTER_ZERO
    if (register_index == PS_REGISTER_DISCARD)
        return;

    state.RegisterValues[register_index] = value;
#endif
}

float4 apply_input_mapping(uniform byte_t input_mapping, uniform float4 value)
{
    switch (input_mapping)
    {
        case PS_INPUTMAPPING_UNSIGNED_IDENTITY:  // = 0x00L : y = max(0,x)       =  1*max(0,x) - 0.0
            return         abs(value); // TODO : Verify, abs() might need to be replaced by max() operation?
        case PS_INPUTMAPPING_UNSIGNED_INVERT:    // = 0x20L : y = 1 - max(0,x)   = -1*max(0,x) + 1.0
            return  1.0f - abs(value); // TODO : Verify, abs() might need to do be removed, or replaced by max() operation?
        case PS_INPUTMAPPING_EXPAND_NORMAL:      // = 0x40L : y = 2*max(0,x) - 1 =  2*max(0,x) - 1.0
            return  2.0f * max(0.0f, value) - 1.0f;
        case PS_INPUTMAPPING_EXPAND_NEGATE:      // = 0x60L : y = 1 - 2*max(0,x) = -2*max(0,x) + 1.0
            return -2.0f * max(0.0f, value) + 1.0f;
        case PS_INPUTMAPPING_HALFBIAS_NORMAL:    // = 0x80L : y = max(0,x) - 1/2 =  1*max(0,x) - 0.5
            return  1.0f * max(0.0f, value) - 0.5f;
        case PS_INPUTMAPPING_HALFBIAS_NEGATE:    // = 0xa0L : y = 1/2 - max(0,x) = -1*max(0,x) + 0.5
            return -1.0f * max(0.0f, value) + 0.5f;
        default: // This also handles
		// case PS_INPUTMAPPING_SIGNED_IDENTITY: // = 0xc0L : y = x              =  1*      x  - 0.0
            return value;            
        case PS_INPUTMAPPING_SIGNED_NEGATE:      // = 0xe0L : y = -x             = -1*      x  + 0.0
            return -value;
    }
/* TODO : Instead of the above, mimick apply_output_mapping by determining bias and scale upfront
	// Determine bias
    float input_bias;
    if (get_bit_?_from_byte(input_mapping))
        input_bias = 0.5f;
    else
        input_bias = 0.0f;

	// Determine scale
    float input_scale;

	// TODO : Apply input biasing and scaling (and possibly clamping?) in that order (or clamp first?) :
    return input_scale * (value - input_bias);
*/
}

float4 apply_output_mapping(uniform byte_t output_mapping, uniform float4 value)
{
	// Determine bias
    float output_bias;
    if (get_bit_3_from_byte(output_mapping))
        output_bias = 0.5f; // Initial PS_COMBINEROUTPUT_BIAS (or any other PS_COMBINEROUTPUT_*_BIAS)
    else
        output_bias = 0.0f; // Initial PS_COMBINEROUTPUT_IDENTITY (or any other PS_COMBINEROUTPUT_* without _BIAS)

	// Determine scale
    float output_scale;
    if (get_bit_4_from_byte(output_mapping))
    {
        if (get_bit_5_from_byte(output_mapping))
            output_scale = 0.5f; // Apply PS_COMBINEROUTPUT_SHIFTRIGHT_1 (or PS_COMBINEROUTPUT_SHIFTRIGHT_1_BIAS)
        else
            output_scale = 4.0f; // Apply PS_COMBINEROUTPUT_SHIFTLEFT_2 (or PS_COMBINEROUTPUT_SHIFTLEFT_2_BIAS)
    }
    else
    {
        if (get_bit_5_from_byte(output_mapping))
            output_scale = 2.0f; // Apply PS_COMBINEROUTPUT_SHIFTLEFT_1 (or PS_COMBINEROUTPUT_SHIFTLEFT_1_BIAS)
        else
            output_scale = 1.0f; // Keep PS_COMBINEROUTPUT_IDENTITY or PS_COMBINEROUTPUT_BIAS
    }

	// Apply output biasing, scaling and clamping (in that order) :
    return clamp(output_scale * (value - output_bias), -1.0f, 1.0f);
}

float4 get_input_register_as_float4(uniform ps_state state, uniform byte_t reg_byte, uniform bool is_alpha = false)
{
    float4 input_result = 0;

    byte_t register_index = mask_register(reg_byte);
    byte_t input_mapping = mask_inputmapping(reg_byte);

    if (state.stage < STAGE_FINAL_COMBINER)
	{
		// Below are invalid in color combiners (only final combiner) :
        switch (register_index)
        {
			case PS_REGISTER_FOG:
				// Only FOG.RGB may be read in color combiner stages :
                if (is_alpha)
                {
#ifdef ASSERT_VALID_ACCESSES
					assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
		            return input_result;
#endif
                }
                break;
            case PS_REGISTER_V1R0_SUM:
            case PS_REGISTER_EF_PROD:
#ifdef ASSERT_VALID_ACCESSES
				assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
                return input_result;
#else
                break;
#endif
		}
	}
    else // (state.stage >= STAGE_FINAL_COMBINER)
    {
		// Below are only valid for final combiner :
        switch (register_index)
        {
            case PS_REGISTER_FOG:
				// Only FOG.Alpha may be read in final combiner stages :
                if (!is_alpha)
                {
#ifdef ASSERT_VALID_ACCESSES
					assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
		            return input_result;
#endif
                }
                break;
			case PS_REGISTER_V1R0_SUM:
			case PS_REGISTER_EF_PROD:
				if (state.stage == STAGE_FINAL_COMBINER_ABCD)
                    break; // Allowed for A,B,C,D inputs, forbidden for E,F,G inputs
				else
                {
#ifdef ASSERT_VALID_ACCESSES
					assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
		            return input_result;
#endif
                }
				break;
        }

		// Below are invalid for final combiner :
        switch (input_mapping)
        {
            case PS_INPUTMAPPING_EXPAND_NORMAL:
            case PS_INPUTMAPPING_HALFBIAS_NORMAL:
            case PS_INPUTMAPPING_SIGNED_IDENTITY:
#ifdef ASSERT_VALID_ACCESSES
                assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
                input_mapping = PS_INPUTMAPPING_UNSIGNED_IDENTITY;
#endif				
                break;
            case PS_INPUTMAPPING_EXPAND_NEGATE:
            case PS_INPUTMAPPING_HALFBIAS_NEGATE:
            case PS_INPUTMAPPING_SIGNED_NEGATE:
#ifdef ASSERT_VALID_ACCESSES
                assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
                input_mapping = PS_INPUTMAPPING_UNSIGNED_INVERT;
#endif				
                break;
        }
    }

    switch (register_index)
    {
        case PS_REGISTER_C0:
            if (state.stage < STAGE_FINAL_COMBINER)
            {
                if (state.FlagUniqueC0)
                    input_result = D3DRS_PSCONSTANT0[state.stage]; // PS_COMBINERCOUNT_UNIQUE_C0
                else
                    input_result = D3DRS_PSCONSTANT0[0u]; // PS_COMBINERCOUNT_SAME_C0
            }
            else
				input_result = D3DRS_PSFINALCOMBINERCONSTANT[0u];
            break;
        case PS_REGISTER_C1:
            if (state.stage < STAGE_FINAL_COMBINER)
            {
                if (state.FlagUniqueC1)
                    input_result = D3DRS_PSCONSTANT1[state.stage]; // PS_COMBINERCOUNT_UNIQUE_C1
                else
                    input_result = D3DRS_PSCONSTANT1[0u]; // PS_COMBINERCOUNT_SAME_C1
            }
            else
                input_result = D3DRS_PSFINALCOMBINERCONSTANT[1u];
            break;
        case PS_REGISTER_FOG:{
                float4 FOG_value = get_plain_register_as_float4(state, PS_REGISTER_FOG);
                if (state.stage < STAGE_FINAL_COMBINER)
                    input_result = float4(FOG_value.rgb, DEFAULT_ALPHA);
                else
                    input_result = float4(DEFAULT_RGB, DEFAULT_RGB, DEFAULT_RGB, FOG_value.a);
            }
            break;
		default:
            input_result = get_plain_register_as_float4(state, register_index);
            break;
    }

    bool use_channel_alpha = get_bit_4_from_byte(reg_byte);
    if (use_channel_alpha)
        input_result = input_result.aaaa; // PS_CHANNEL_ALPHA
    else
        input_result = input_result.rgbb; // PS_CHANNEL_BLUE and PS_CHANNEL_RGB

    return apply_input_mapping(input_mapping, input_result);
}

void set_output_register_rgb(
	inout ps_state state,
	uniform byte_t reg_byte,
	uniform byte_t output_mapping,
	uniform bool flag_BlueToAlpha,
	uniform float3 rgb_value)
{
    byte_t register_index = mask_register(reg_byte);

	// Below are not writeable :
    switch (register_index)
    {
        case PS_REGISTER_C0:
        case PS_REGISTER_C1:
        case PS_REGISTER_FOG:
        case 6: // Unknown r/w state
        case 7: // Unknown r/w state
        case PS_REGISTER_V1R0_SUM:
        case PS_REGISTER_EF_PROD:
#ifdef ASSERT_VALID_ACCESSES
		    assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
            return;
#else
            break;
#endif
    }

    float4 new_value = apply_output_mapping(output_mapping, float4(rgb_value, 0.0f)); // Rely on the optimizer to remove useless operations (otherwise, implement apply_output_mapping_rgb)

	// Mix in the new RGB channels into the destination register :
    if (flag_BlueToAlpha)
        new_value = float4(new_value.rgb, new_value.b); // Put blue channel in alpha (this is the only way RGB can set alpha)
    else
        new_value = float4(new_value.rgb, get_plain_register_as_float4(state, register_index).a); // Keep destination alpha channel intact

    set_plain_register_as_float4(state, register_index, new_value);
}

void set_output_register_alpha(
	inout ps_state state,
	uniform byte_t reg_byte,
	uniform byte_t output_mapping,
	uniform float alpha_value)
{
    byte_t register_index = mask_register(reg_byte);

	// Below are not writeable :
    switch (register_index)
    {
        case PS_REGISTER_C0:
        case PS_REGISTER_C1:
        case PS_REGISTER_FOG:
        case 6: // Unknown r/w state
        case 7: // Unknown r/w state
        case PS_REGISTER_V1R0_SUM:
        case PS_REGISTER_EF_PROD:
#ifdef ASSERT_VALID_ACCESSES
		    assert(false);
#endif
#ifdef AVOID_INVALID_ACCESSES
            return;
#else
            break;
#endif
    }

    float new_alpha = apply_output_mapping(output_mapping, float4(DEFAULT_RGB, DEFAULT_RGB, DEFAULT_RGB, alpha_value)).a; // Rely on the optimizer to remove useless operations (otherwise, implement apply_output_mapping_alpha)

	// Only mix in the new alpha value into the destination register :
    float4 new_value = float4(get_plain_register_as_float4(state, register_index).rgb, new_alpha);

    set_plain_register_as_float4(state, register_index, new_value);
}

void do_color_combiner_stage(inout ps_state state, uniform bool is_alpha)
{
	// Fetch input register bytes :
    byte_t A_input_reg_byte = is_alpha ? get_byte_0_from_float4(D3DRS_PSALPHAINPUTS[state.stage]) : get_byte_0_from_float4(D3DRS_PSRGBINPUTS[state.stage]);
    byte_t B_input_reg_byte = is_alpha ? get_byte_1_from_float4(D3DRS_PSALPHAINPUTS[state.stage]) : get_byte_1_from_float4(D3DRS_PSRGBINPUTS[state.stage]);
    byte_t C_input_reg_byte = is_alpha ? get_byte_2_from_float4(D3DRS_PSALPHAINPUTS[state.stage]) : get_byte_2_from_float4(D3DRS_PSRGBINPUTS[state.stage]);
    byte_t D_input_reg_byte = is_alpha ? get_byte_3_from_float4(D3DRS_PSALPHAINPUTS[state.stage]) : get_byte_3_from_float4(D3DRS_PSRGBINPUTS[state.stage]);

	// Fetch input (source) register values :
    float4 A_value = get_input_register_as_float4(state, A_input_reg_byte, is_alpha); // A.k.a. 's0'
    float4 B_value = get_input_register_as_float4(state, B_input_reg_byte, is_alpha); // A.k.a. 's1'
    float4 C_value = get_input_register_as_float4(state, C_input_reg_byte, is_alpha); // A.k.a. 's2'
    float4 D_value = get_input_register_as_float4(state, D_input_reg_byte, is_alpha); // A.k.a. 's3'

	// Fetch output bytes :
    byte_t output_byte_0 = is_alpha ? get_byte_0_from_float4(D3DRS_PSALPHAOUTPUTS[state.stage]) : get_byte_0_from_float4(D3DRS_PSRGBOUTPUTS[state.stage]);
    byte_t output_byte_1 = is_alpha ? get_byte_1_from_float4(D3DRS_PSALPHAOUTPUTS[state.stage]) : get_byte_1_from_float4(D3DRS_PSRGBOUTPUTS[state.stage]);
    byte_t output_byte_2 = is_alpha ? get_byte_2_from_float4(D3DRS_PSALPHAOUTPUTS[state.stage]) : get_byte_2_from_float4(D3DRS_PSRGBOUTPUTS[state.stage]);

	// Decode already required flags from output byte 1 :
    bool combiner_flag_CD_DotProduct = get_bit_4_from_byte(output_byte_1); // PS_COMBINEROUTPUT_CD_DOT_PRODUCT = 0x01L; // RGB only
    bool combiner_flag_AB_DotProduct = get_bit_5_from_byte(output_byte_1); // PS_COMBINEROUTPUT_AB_DOT_PRODUCT = 0x02L; // RGB only
    bool combiner_flag_ABCD_Mux = get_bit_6_from_byte(output_byte_1); // PS_COMBINEROUTPUT_AB_CD_MUX = 0x04L; // 3rd output is MUX(AB;CD) based on R0.a

	// xmma (mul/mul/sum)     > calculating : d0=s0*s1,     d1=s2*s3,     d2=             s2*s3 + s0*s1
	// xmmc (mul/mul/mux)     > calculating : d0=s0*s1,     d1=s2*s3,     d2=(r0.a>0.5) ? s2*s3 : s0*s1
	// xdm  (dot/mul/discard) > calculating : d0=s0 dot s1, d1=s2*s3
	// xdd  (dot/dot/discard) > calculating : d0=s0 dot s1, d1=s2 dot s3

	// Calculate AB side, doing either DOT or MUL :
    float4 AB_value = combiner_flag_AB_DotProduct ? dot(A_value, B_value) : A_value * B_value;
	// Calculate CD side, doing either DOT or MUL :
    float4 CD_value = combiner_flag_CD_DotProduct ? dot(C_value, D_value) : C_value * D_value;

	// Perform the AB,CD operation (either MUX or SUM) :
    float4 AB_CD_value;
    if (combiner_flag_ABCD_Mux) // PS_COMBINEROUTPUT_AB_CD_MUX
    {
        float R0_a_value = get_plain_register_as_float4(state, PS_REGISTER_R0).a; // TODO : Normalize? How?
        if (state.FlagMuxMsb)
            AB_CD_value = (R0_a_value > 0.5f) ? CD_value : AB_value; // PS_COMBINERCOUNT_MUX_MSB
        else
            AB_CD_value = get_bit_0_from_byte(R0_a_value) ? CD_value : AB_value; // PS_COMBINERCOUNT_MUX_LSB
    }
    else // PS_COMBINEROUTPUT_AB_CD_SUM
        AB_CD_value = AB_value + CD_value;

	// Decode further output register numbers :
    byte_t combiner_output_reg_AB = get_nibble_0_from_byte(output_byte_0); // A.k.a. 'd0'
    byte_t combiner_output_reg_CD = get_nibble_1_from_byte(output_byte_0); // A.k.a. 'd1'
    byte_t combiner_output_reg_AB_CD = get_nibble_0_from_byte(output_byte_1); // A.k.a. 'd2'
    byte_t combiner_output_mapping = mask_outputmapping(output_byte_1, output_byte_2);

	// AB_CD register output must be DISCARD if either AB_DOT_PRODUCT or CD_DOT_PRODUCT are set
    if (combiner_flag_AB_DotProduct || combiner_flag_CD_DotProduct)
    {
#ifdef ASSERT_VALID_ACCESSES
		assert(combiner_output_reg_AB_CD == PS_REGISTER_DISCARD);
#endif
#ifdef AVOID_INVALID_ACCESSES
		// Fix-up, just to avoid issues (TODO : Investigate if real hardware does this too?)
        combiner_output_reg_AB_CD = PS_REGISTER_DISCARD;
#endif
    }

	// Store resulting values in output registers :
    if (is_alpha)
    {
        set_output_register_alpha(state, combiner_output_reg_AB, combiner_output_mapping, AB_value.a);
        set_output_register_alpha(state, combiner_output_reg_CD, combiner_output_mapping, CD_value.a);
        set_output_register_alpha(state, combiner_output_reg_AB_CD, combiner_output_mapping, AB_CD_value.a);
    }
    else
    {
		// Decode the last two output register flags (for RGB only) :
		bool flag_CD_BlueToAlpha = get_bit_2_from_byte(output_byte_2); // PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA = 0x40L; // RGB only
        bool flag_AB_BlueToAlpha = get_bit_3_from_byte(output_byte_2); // PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA = 0x80L; // RGB only

        set_output_register_rgb(state, combiner_output_reg_AB, combiner_output_mapping, flag_AB_BlueToAlpha, AB_value.rgb);
        set_output_register_rgb(state, combiner_output_reg_CD, combiner_output_mapping, flag_CD_BlueToAlpha, CD_value.rgb);
        set_output_register_rgb(state, combiner_output_reg_AB_CD, combiner_output_mapping, /*flag_BlueToAlpha=*/false, AB_CD_value.rgb);
    }
}

float4 do_final_combiner(inout ps_state state)
{
	// Fetch R0 register value (in final combiner used are input, otherwise returned) :
    float4 R0_value = get_plain_register_as_float4(state, PS_REGISTER_R0);

	// TODO : Is the following actually required? D3D seems to always set the final combiner register...
	// Check if the final combiner doesn't need to run :
    if (!any(D3DRS_PSFINALCOMBINERINPUTSABCD) && !any(D3DRS_PSFINALCOMBINERINPUTSEFG))
        return R0_value;

	// Decompose the final combiner inputs into separate bytes :
    byte_t E_reg_byte = get_byte_0_from_float4(D3DRS_PSFINALCOMBINERINPUTSEFG);
    byte_t F_reg_byte = get_byte_1_from_float4(D3DRS_PSFINALCOMBINERINPUTSEFG);
    byte_t G_reg_byte = get_byte_2_from_float4(D3DRS_PSFINALCOMBINERINPUTSEFG);
    byte_t setting_byte = get_byte_3_from_float4(D3DRS_PSFINALCOMBINERINPUTSEFG);

	// Tell get_input_register_as_float4() we're now in the final combiner stage for E,F,G :
    state.stage = STAGE_FINAL_COMBINER;

	// Fetch the E and F register values :
    float3 E_value = get_input_register_as_float4(state, E_reg_byte).rgb; // A.k.a. 's4'
    float3 F_value = get_input_register_as_float4(state, F_reg_byte).rgb; // A.k.a. 's5'
    float G_alpha = get_input_register_as_float4(state, G_reg_byte, /*is_alpha=*/true).a; // A.k.a. 's6' // Note : Set is_alpha argument merely avoids FOG errors

	// Calculate the product of E*F and set the result in it's dedicated register :
    float3 EF_PROD_value = E_value * F_value;
    set_plain_register_as_float4(state, PS_REGISTER_EF_PROD, float4(EF_PROD_value, DEFAULT_ALPHA));

	// Do optional complement operation on R0 register value :
    if (get_bit_5_from_byte(setting_byte)) // PS_FINALCOMBINERSETTING_COMPLEMENT_R0
    {
		// unsigned invert mapping  (1 - r0) is used as an input to the sum rather than r0
        R0_value = 1.0f - R0_value; // TODO : Is abs() needed here?
    }

	// Fetch V1 register value (including an optional complement operation) :
    float3 V1_value = get_plain_register_as_float4(state, PS_REGISTER_V1).rgb;
    if (get_bit_6_from_byte(setting_byte)) // PS_FINALCOMBINERSETTING_COMPLEMENT_V1
    {
		// unsigned invert mapping  (1 - v1) is used as an input to the sum rather than v1
        V1_value = 1.0f - V1_value; // TODO : Is abs() needed here?
    }

	// Calculate V1 + R0 value (including an optional clamping operation) :
    float3 V1R0_SUM_value = V1_value + R0_value.rgb;
    if (get_bit_7_from_byte(setting_byte)) // PS_FINALCOMBINERSETTING_CLAMP_SUM
    {
		// V1+R0 sum clamped to [0,1]
        V1R0_SUM_value = clamp(V1R0_SUM_value, 0.0f, 1.0f);
    }

	// Set the V1+R0 result in it's dedicated register :
    set_plain_register_as_float4(state, PS_REGISTER_V1R0_SUM, float4(V1R0_SUM_value, DEFAULT_ALPHA));

	// Decompose the final combiner inputs into separate bytes :
    byte_t A_reg_byte = get_byte_0_from_float4(D3DRS_PSFINALCOMBINERINPUTSABCD);
    byte_t B_reg_byte = get_byte_1_from_float4(D3DRS_PSFINALCOMBINERINPUTSABCD);
    byte_t C_reg_byte = get_byte_2_from_float4(D3DRS_PSFINALCOMBINERINPUTSABCD);
    byte_t D_reg_byte = get_byte_3_from_float4(D3DRS_PSFINALCOMBINERINPUTSABCD);

	// Tell get_input_register_as_float4() we're now in the final combiner stage for A,B,C,D :
    state.stage = STAGE_FINAL_COMBINER_ABCD;

	// Fetch the A, B, C and D register values :
    float4 A_value = get_input_register_as_float4(state, A_reg_byte); // A.k.a. 's0'
    float4 B_value = get_input_register_as_float4(state, B_reg_byte); // A.k.a. 's1'
    float4 C_value = get_input_register_as_float4(state, C_reg_byte); // A.k.a. 's2'
    float4 D_value = get_input_register_as_float4(state, D_reg_byte); // A.k.a. 's3'

	// Calculate output RGB = a*b + (1-a)*c + d :
    float4 FinalOutput;
    FinalOutput.rgb = lerp(C_value.rgb, B_value.rgb, A_value.rgb) + D_value.rgb;

	// Fetch the output alpha channel from the G register value :
    FinalOutput.a = G_alpha;

    return FinalOutput;
}

texture tex0;
texture tex1;
texture tex2;
texture tex3;
sampler samp2d; // TODO : Implement all possible samplers

struct VS_OUTPUT // TODO : Complete and pass this in through code (and/or vertex shader?)
{
//	float4 Position : SV_POSITION; // vertex position 
	float4 Diffuse : COLOR0; // vertex diffuse color (note that COLOR0 is clamped from 0..1)
	float4 Specular : COLOR0; // vertex specular color (note that COLOR0 is clamped from 0..1)
	float4 Fog : COLOR0;
    float4 TextureCoords[4] : TEXCOORD0; // vertex texture coords 
} In;

float4 main() : COLOR
{
    nointerpolation ps_state state = (ps_state) 0; // Clearing like this avoids error X3508: 'do_color_combiner_stage': output parameter 'state' not completely initialized
	// Note, the above also sets state.RegisterValues[PS_REGISTER_ZERO] to 0!

	// Process all four texture stages :
    for (byte_t texture_stage = 0u; texture_stage < 4u; texture_stage++)
    {
		// TODO : Fully implement texture fetch (preferrably in this shader, including conversion of X_D3DFMT_P8 and other unsupported textures formats)
        float4 texture_value = tex2D(samp2d, (float2) In.TextureCoords[texture_stage]); // or default float4(0,0,0,1)

        byte_t texture_register = PS_REGISTER_T0 + texture_stage;
        texture_register = max(texture_register, PS_REGISTER_T3); // Avoids warning X3550: array reference cannot be used as an l-value; not natively addressable, forcing loop to unroll 
        set_plain_register_as_float4(state, texture_register, texture_value);
    }

	// Calculate initial register values :
    float4 FOG_value = In.Fog; // TODO : Is this correct, or use float4(1.0f, 1.0f, 1.0f, DEFAULT_ALPHA)?
    float4 V0_value = In.Diffuse; // TODO : Is this correct?
    float4 V1_value = In.Specular; // TODO : Is this correct?
    float4 R0_value = float4(1.0f, 1.0f, 1.0f, DEFAULT_ALPHA); // TODO : Is this correct?

	// Set initial register values :
	set_plain_register_as_float4(state, PS_REGISTER_FOG, FOG_value);
    set_plain_register_as_float4(state, PS_REGISTER_V0, V0_value);
    set_plain_register_as_float4(state, PS_REGISTER_V1, V1_value);
    set_plain_register_as_float4(state, PS_REGISTER_R0, R0_value);

	// Decode the global combiner count and flags :
    byte_t combiner_byte_0 = get_byte_0_from_float4(D3DRS_PSCOMBINERCOUNT);
    byte_t combiner_byte_1 = get_byte_1_from_float4(D3DRS_PSCOMBINERCOUNT);
    byte_t combiner_byte_2 = get_byte_2_from_float4(D3DRS_PSCOMBINERCOUNT);
    byte_t num_stages = combiner_byte_0;
    state.FlagMuxMsb = get_bit_0_from_byte(combiner_byte_1); // PS_COMBINERCOUNT_MUX_MSB = 0x0001L, // mux on r0.a msb
    state.FlagUniqueC0 = get_bit_4_from_byte(combiner_byte_1); // PS_COMBINERCOUNT_UNIQUE_C0 = 0x0010L, // c0 unique in each stage
    state.FlagUniqueC1 = get_bit_0_from_byte(combiner_byte_2); // PS_COMBINERCOUNT_UNIQUE_C1 = 0x0100L // c1 unique in each stage

#ifdef ASSERT_VALID_ACCESSES
    assert(num_stages >= 1);
    assert(num_stages <= 8);
#endif
#ifdef AVOID_INVALID_ACCESSES
    num_stages = min(num_stages, 1u);
    num_stages = max(num_stages, 8u);
#endif

	// Loop over at most 8 combiner stages (start counting from zero)
    for (byte_t combiner_stage = 0u; combiner_stage < MAX_COMBINER_STAGE_COUNT; combiner_stage++)
    {
        if (combiner_stage >= num_stages)
            break;

        state.stage = combiner_stage; // tell do_color_combiner_stage() and get_input_register_as_float4() the currently active stage
        do_color_combiner_stage(state, false); // for RGB
        do_color_combiner_stage(state, true); // for Alpha
    }

    return do_final_combiner(state);
}
