//
// Cxbx-Reloaded Xbox Pixel Shader Interpreter, written in Direct3D 9.0 Shader Model 3.0
//
// This shader interprets Xbox Pixel Shaders, as declared in the Xbox Direct3D__RenderState registers.
// These registers must be fed into this shader though host SetPixelShaderConstantF(),
// so that they can be accessed from the corresponding constant registers.
//
// Since Shader Model doesn't support bitwise operations nor native integers,
// the below code expects 32 bit input data to be offered split up in two parts.
// The Lower16 part must contain the lower 16 bits of any input, expressed as a float.
// The Upper16 part must contain the upper 16 bits of any input, expressed as a float.
// Using this representation, it's possible to extract bits from floats via fmod(),
// as described in http://theinstructionlimit.com/encoding-boolean-flags-into-a-float-in-hlsl
// Also, see https://en.wikipedia.org/wiki/Single-precision_floating-point_format

//
// Utility functions
//

#define AVOID_INVALID_ACCESSES
//#define ASSERT_VALID_ACCESSES

typedef float floor_t; // Indicates an already floor()'ed float
typedef float floor1_t; // Indicates an already floor()'ed float, limited to lowest 8 bits
typedef float floor2_t; // Indicates an already floor()'ed float, limited to lowest 16 bits

#ifdef ASSERT_VALID_ACCESSES
void assert(bool condition)
{
	// TODO
}
#endif

bool get_bit_0_from_value(const floor_t value)
{
    return fmod(value, 2) == 1;
}

bool get_bit_1_from_value(const floor_t value)
{
    return fmod(value, 4) >= 2;
}

bool get_bit_2_from_value(const floor_t value)
{
    return fmod(value, 8) >= 4;
}

bool get_bit_3_from_value(const floor_t value)
{
    return fmod(value, 16) >= 8;
}

bool get_bit_4_from_value(const floor_t value)
{
    return fmod(value, 32) >= 16;
}

bool get_bit_5_from_value(const floor_t value)
{
    return fmod(value, 64) >= 32;
}

bool get_bit_6_from_value(const floor_t value)
{
    return fmod(value, 128) >= 64;
}

bool get_bit_7_from_value(const floor_t value)
{
    return fmod(value, 256) >= 128;
}

bool get_bit_8_from_value(const floor_t value)
{
    return fmod(value, 512) >= 256;
}

bool get_bit_9_from_value(const floor_t value)
{
    return fmod(value, 1024) >= 512;
}

bool get_bit_10_from_value(const floor_t value)
{
    return fmod(value, 2048) >= 1024;
}

bool get_bit_11_from_value(const floor_t value)
{
    return fmod(value, 4096) >= 2048;
}

bool get_bit_12_from_value(const floor_t value)
{
    return fmod(value, 8192) >= 4096;
}

bool get_bit_13_from_value(const floor_t value)
{
    return fmod(value, 16384) >= 8192;
}

bool get_bit_14_from_value(const floor_t value)
{
    return fmod(value, 32768) >= 16384;
}

bool get_bit_15_from_value(const floor_t value)
{
    return fmod(value, 65536) >= 32768;
}

bool get_bit_16_from_value(const floor_t value)
{
    return fmod(value, 131072) >= 65536;
}

bool get_bit_17_from_value(const floor_t value)
{
    return fmod(value, 262144) >= 131072;
}

bool get_bit_18_from_value(const floor_t value)
{
    return fmod(value, 524288) >= 262144;
}

bool get_bit_19_from_value(const floor_t value)
{
    return fmod(value, 1048576) >= 524288;
}

bool get_bit_20_from_value(const floor_t value)
{
    return fmod(value, 2097152) >= 1048576;
}

bool get_bit_21_from_value(const floor_t value)
{
    return fmod(value, 4194304) >= 2097152;
}

bool get_bit_22_from_value(const floor_t value)
{
    return fmod(value, 8388608) >= 4194304;
}

uint1 get_byte_0_from_word(const floor2_t value)
{
#ifdef ASSERT_VALID_ACCESSES
    assert(value <= 65535);
#endif
    return (uint1) abs(fmod(value, 256));
}

uint1 get_byte_1_from_word(const floor2_t value)
{
#ifdef ASSERT_VALID_ACCESSES
    assert(value <= 65535);
#endif
    return (uint) abs(value / 256); // Cast avoids warning X3556: integer divides may be much slower, try using uints if possible.
}

uint1 get_nibble_0_from_byte(const floor1_t value)
{
#ifdef ASSERT_VALID_ACCESSES
    assert(value <= 255);
#endif
    return (uint1) abs(fmod(value, 16));
}

uint1 get_nibble_1_from_byte(const floor1_t value)
{
#ifdef ASSERT_VALID_ACCESSES
    assert(value <= 255);
#endif
    return (uint) abs(value / 16); // Cast avoids warning X3556: integer divides may be much slower, try using uints if possible.
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
// PS_INPUTMAPPING_UNSIGNED_IDENTITY : y =  1*max(0,x) + 0.0  ?  OR  abs(x)                        ?
// PS_INPUTMAPPING_UNSIGNED_INVERT   : y = -1*max(0,x) + 1.0  ?  OR  1 - x   MAYBE EVEN  1-abs(x)  ?
// PS_INPUTMAPPING_EXPAND_NORMAL     : y =  2*max(0,x) - 1.0
// PS_INPUTMAPPING_EXPAND_NEGATE     : y = -2*max(0,x) + 1.0
// PS_INPUTMAPPING_HALFBIAS_NORMAL   : y =  1*max(0,x) - 0.5
// PS_INPUTMAPPING_HALFBIAS_NEGATE   : y = -1*max(0,x) + 0.5
// PS_INPUTMAPPING_SIGNED_IDENTITY   : y =  1*      x  + 0.0
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
static const uint PS_INPUTMAPPING_UNSIGNED_IDENTITY = 0x00L; // max(0;x)         OK for final combiner: y = abs(x)
static const uint PS_INPUTMAPPING_UNSIGNED_INVERT =   0x20L; // 1 - max(0;x)     OK for final combiner: y = 1 - x
static const uint PS_INPUTMAPPING_EXPAND_NORMAL =     0x40L; // 2*max(0;x) - 1   invalid for final combiner
static const uint PS_INPUTMAPPING_EXPAND_NEGATE =     0x60L; // 1 - 2*max(0;x)   invalid for final combiner
static const uint PS_INPUTMAPPING_HALFBIAS_NORMAL =   0x80L; // max(0;x) - 1/2   invalid for final combiner
static const uint PS_INPUTMAPPING_HALFBIAS_NEGATE =   0xa0L; // 1/2 - max(0;x)   invalid for final combiner
static const uint PS_INPUTMAPPING_SIGNED_IDENTITY =   0xc0L; // x                invalid for final combiner
static const uint PS_INPUTMAPPING_SIGNED_NEGATE =     0xe0L; // -x               invalid for final combiner

// Pixel Shader Registers
static const uint PS_REGISTER_ZERO =              0x00L; // r
static const uint PS_REGISTER_DISCARD =           0x00L; // w
static const uint PS_REGISTER_C0 =                0x01L; // r
static const uint PS_REGISTER_C1 =                0x02L; // r
static const uint PS_REGISTER_FOG =               0x03L; // r
static const uint PS_REGISTER_V0 =                0x04L; // r/w
static const uint PS_REGISTER_V1 =                0x05L; // r/w
// Filling the PS_REGISTER_* gap here with V2 and V3 is just guesswork,
// they could be unsupported, read-only, special-purpose, or whatever else.
// Needs investigation.
//static const uint PS_REGISTER_V2 =                0x06L; // unknown r/w state
//static const uint PS_REGISTER_V3 =                0x07L; // unknown r/w state
static const uint PS_REGISTER_T0 =                0x08L; // r/w
static const uint PS_REGISTER_T1 =                0x09L; // r/w
static const uint PS_REGISTER_T2 =                0x0aL; // r/w
static const uint PS_REGISTER_T3 =                0x0bL; // r/w
static const uint PS_REGISTER_R0 =                0x0cL; // r/w
static const uint PS_REGISTER_R1 =                0x0dL; // r/w
static const uint PS_REGISTER_V1R0_SUM =          0x0eL; // r
static const uint PS_REGISTER_EF_PROD =           0x0fL; // r
static const uint PS_REGISTER_ONE =               PS_REGISTER_ZERO | PS_INPUTMAPPING_UNSIGNED_INVERT; // 0x20 OK for final combiner
static const uint PS_REGISTER_NEGATIVE_ONE =      PS_REGISTER_ZERO | PS_INPUTMAPPING_EXPAND_NORMAL;   // 0x40 invalid for final combiner
static const uint PS_REGISTER_ONE_HALF =          PS_REGISTER_ZERO | PS_INPUTMAPPING_HALFBIAS_NEGATE; // 0xa0 invalid for final combiner
static const uint PS_REGISTER_NEGATIVE_ONE_HALF = PS_REGISTER_ZERO | PS_INPUTMAPPING_HALFBIAS_NORMAL; // 0x80 invalid for final combiner
static const uint PS_REGISTER_CXBX_PROD =         PS_REGISTER_ZERO | PS_INPUTMAPPING_SIGNED_IDENTITY; // Cxbx internal use

static const uint PS_CHANNEL_RGB =   0x00; // used as RGB source
static const uint PS_CHANNEL_BLUE =  0x00; // used as ALPHA source
static const uint PS_CHANNEL_ALPHA = 0x10; // used as RGB or ALPHA source

static const uint PS_FINALCOMBINERSETTING_CLAMP_SUM =     0x80; // V1+R0 sum clamped to [0;1]
static const uint PS_FINALCOMBINERSETTING_COMPLEMENT_V1 = 0x40; // unsigned invert mapping  (1 - v1) is used as an input to the sum rather than v1
static const uint PS_FINALCOMBINERSETTING_COMPLEMENT_R0 = 0x20; // unsigned invert mapping  (1 - r0) is used as an input to the sum rather than r0

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
static const uint PS_COMBINEROUTPUT_IDENTITY             = 0x00L; // y = x
static const uint PS_COMBINEROUTPUT_BIAS                 = 0x08L; // y = x - 0.5
static const uint PS_COMBINEROUTPUT_SHIFTLEFT_1          = 0x10L; // y = x*2
static const uint PS_COMBINEROUTPUT_SHIFTLEFT_1_BIAS     = 0x18L; // y = (x - 0.5)*2
static const uint PS_COMBINEROUTPUT_SHIFTLEFT_2          = 0x20L; // y = x*4
// static const uint PS_COMBINEROUTPUT_SHIFTLEFT_2_BIAS  = 0x28L; // y = (x - 0.5)*4
static const uint PS_COMBINEROUTPUT_SHIFTRIGHT_1         = 0x30L; // y = x/2
// static const uint PS_COMBINEROUTPUT_SHIFTRIGHT_1_BIAS = 0x38L; // y = (x - 0.5)/2
static const uint PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA     = 0x80L; // RGB only
static const uint PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA     = 0x40L; // RGB only
static const uint PS_COMBINEROUTPUT_AB_MULTIPLY          = 0x00L;
static const uint PS_COMBINEROUTPUT_AB_DOT_PRODUCT       = 0x02L; // RGB only
static const uint PS_COMBINEROUTPUT_CD_MULTIPLY          = 0x00L;
static const uint PS_COMBINEROUTPUT_CD_DOT_PRODUCT       = 0x01L; // RGB only
static const uint PS_COMBINEROUTPUT_AB_CD_SUM            = 0x00L; // 3rd output is AB+CD
static const uint PS_COMBINEROUTPUT_AB_CD_MUX            = 0x04L; // 3rd output is MUX(AB;CD) based on R0.a
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
static const uint MAX_COMBINER_STAGE_COUNT = 8;
static const uint FINAL_COMBINER_STAGE = MAX_COMBINER_STAGE_COUNT; // Used for E,F,G inputs
static const uint FINAL_COMBINER_STAGE_ABCD = FINAL_COMBINER_STAGE + 1;
static const float DEFAULT_ALPHA = 1.0f;

#ifdef INPUT_BUFFER // TODO : Expand this for Shader Model 4 and up :
// Xbox Pixel Shader definition (a uint array) D3D__RenderState member indices :
#else // Shader Model 3 :

//
// Xbox Pixel Shader definition (32 bit dword/float values, taken from "D3D__RenderState")
//

// Since shader model 3.0 seems to have no input-buffer feature, and using texture-sampling might interpolate values, use host pixel shader constants to access the Xbox renderstate input :
float D3DRS_PSALPHAINPUTS_Lower16            [8] : register(c0);
float D3DRS_PSFINALCOMBINERINPUTSABCD_Lower16    : register(c8);
float D3DRS_PSFINALCOMBINERINPUTSEFG_Lower16     : register(c9);
float D3DRS_PSCONSTANT0                      [8] : register(c10); // Read as float, so no split upper16/lower16
float D3DRS_PSCONSTANT1                      [8] : register(c18); // Read as float, so no split upper16/lower16
float D3DRS_PSALPHAOUTPUTS_Lower16           [8] : register(c26);
float D3DRS_PSRGBINPUTS_Lower16              [8] : register(c34);
float D3DRS_PSCOMPAREMODE_Lower16                : register(c42); // Read as float, so no split upper16/lower16
float D3DRS_PSFINALCOMBINERCONSTANT          [2] : register(c43);
float D3DRS_PSRGBOUTPUTS_Lower16             [8] : register(c45);
float D3DRS_PSCOMBINERCOUNT_Lower16              : register(c53);
float D3DRS_PS_RESERVED_Lower16                  : register(c54); // Dxbx note : This takes the slot of X_D3DPIXELSHADERDEF.PSTextureModes); set by D3DDevice_SetRenderState_LogicOp?
float D3DRS_PSDOTMAPPING_Lower16                 : register(c55);
float D3DRS_PSINPUTTEXTURE_Lower16               : register(c56);
// Keep a few indices free between lower16 and upper16 renderstate input constants, in case we need to pass in more data.
float D3DRS_PSALPHAINPUTS_Upper16            [8] : register(c64);
float D3DRS_PSFINALCOMBINERINPUTSABCD_Upper16    : register(c72);
float D3DRS_PSFINALCOMBINERINPUTSEFG_Upper16     : register(c73);
//float D3DRS_PSCONSTANT0_Upper16              [8] : register(c74);
//float D3DRS_PSCONSTANT1_Upper16              [8] : register(c82);
float D3DRS_PSALPHAOUTPUTS_Upper16           [8] : register(c90);
float D3DRS_PSRGBINPUTS_Upper16              [8] : register(c98);
float D3DRS_PSCOMPAREMODE_Upper16                : register(c106);
//float D3DRS_PSFINALCOMBINERCONSTANT_Upper16  [2] : register(c107);
float D3DRS_PSRGBOUTPUTS_Upper16             [8] : register(c109);
//float D3DRS_PSCOMBINERCOUNT_Upper16              : register(c117); // unused
float D3DRS_PS_RESERVED_Upper16                  : register(c118); // Dxbx note : This takes the slot of X_D3DPIXELSHADERDEF.PSTextureModes); set by D3DDevice_SetRenderState_LogicOp?
float D3DRS_PSDOTMAPPING_Upper16                 : register(c119);
float D3DRS_PSINPUTTEXTURE_Upper16               : register(c120);

#endif

//
// Our Xbox Pixel Shader Interpreter state
//
typedef
struct
{
	// Flags, used in get_input_register_as_float4() :
    uint stage;        // Currently active stage (0 upto 7, 8 is for final combiner)
    bool FlagMuxMsb;   // Mux on r0.a lsb or msb, 0 = PS_COMBINERCOUNT_MUX_LSB, 1 = PS_COMBINERCOUNT_MUX_MSB
    bool FlagUniqueC0; // C[0] unique in each stage, 0 = PS_COMBINERCOUNT_SAME_C0, 1 = PS_COMBINERCOUNT_UNIQUE_C0
    bool FlagUniqueC1; // C[1] unique in each stage, 0 = PS_COMBINERCOUNT_SAME_C1, 1 = PS_COMBINERCOUNT_UNIQUE_C1

	// Xbox NV2A pixel shader register values, in order of PS_REGISTER_* :
	// float4 ZERO;
    // float4 C[2]; // Constant registers, read-only
    // float4 FOG;  // Note, FOG ALPHA is only available in final combiner
    // float4 V[2]; // Vertex color registers, read/write (V[0] = diffuse, V[1] = specular
	// float4 Unknown[2];
    // float4 T[4]; // Texture registers, initially sampled or (0,0,0,1), after texture-addressing they're read/write just like temporary registers
    // float4 R[2]; // Temporary registers, read/write (R[0].a (R0_ALPHA) is initialized to T[0].s (T0_ALPHA) in stage 0); Final result is in R[0]
    // float3 V1R0_SUM;	// Note : V1R0_SUM and EF_PROD are only available in final combiner (A,B,C,D inputs only)
    // float3 EF_PROD;  // Note : V1R0_SUM_ALPHA and EF_PROD_ALPHA are not available, hence the float3 type here
    float4 RegisterValues[16]; // Room for 16 registers : PS_REGISTER_ZERO (0) upto PS_REGISTER_EF_PROD (15)
} ps_state;

uint1 mask_register(const uint1 value)
{
    return fmod(abs(value), 16); // == PS_REGISTER_EF_PROD + 1
}

uint1 mask_inputmapping(const uint1 value)
{
    floor1_t Result = abs(value);
    Result = Result / 2;
    Result = fmod(Result, 8); // == (PS_INPUTMAPPING_SIGNED_NEGATE/2) + 1;
    Result = Result + Result;
    return (uint1) Result;
}

float4 get_input_register_as_float4(const ps_state state, const uint1 reg_byte, const bool is_alpha = false)
{
    float4 Result = 0;

    uint register_index = mask_register(reg_byte);
    uint input_mapping = mask_inputmapping(reg_byte);
#ifdef AVOID_INVALID_ACCESSES
    if (state.stage >= FINAL_COMBINER_STAGE)
	{
		// Below are only valid for final combiner :
        switch (register_index)
		{
            case PS_REGISTER_V1R0_SUM:
			case PS_REGISTER_EF_PROD:
				if (state.stage == FINAL_COMBINER_STAGE_ABCD)
					; // Allowed for A,B,C,D inputs, forbidden for E,F,G inputs
                else
                {
#ifdef ASSERT_VALID_ACCESSES
					assert(false);
#endif
                    return Result;
                }
        }

		// Below are invalid for final combiner :
        switch (input_mapping)
        {
            case PS_INPUTMAPPING_EXPAND_NORMAL:
            case PS_INPUTMAPPING_EXPAND_NEGATE:
            case PS_INPUTMAPPING_HALFBIAS_NORMAL:
            case PS_INPUTMAPPING_HALFBIAS_NEGATE:
            case PS_INPUTMAPPING_SIGNED_IDENTITY:
            case PS_INPUTMAPPING_SIGNED_NEGATE:
	#ifdef ASSERT_VALID_ACCESSES
                assert(false);
	#endif
                return Result;
        }
    }

#endif
    switch (register_index)
    {
        case PS_REGISTER_C0:
            if (state.stage >= FINAL_COMBINER_STAGE)
                Result = D3DRS_PSFINALCOMBINERCONSTANT[0];
            else if (state.FlagUniqueC0)
                Result = D3DRS_PSCONSTANT0[state.stage];
            else
                Result = D3DRS_PSCONSTANT0[0];
            break;
        case PS_REGISTER_C1:
            if (state.stage >= FINAL_COMBINER_STAGE)
                Result = D3DRS_PSFINALCOMBINERCONSTANT[1];
            else if (state.FlagUniqueC1)
                Result = D3DRS_PSCONSTANT1[state.stage];
            else
                Result = D3DRS_PSCONSTANT1[1];
            break;
        case PS_REGISTER_FOG: {
                float4 FOG_value;
                FOG_value = state.RegisterValues[PS_REGISTER_FOG];
                if (state.stage >= FINAL_COMBINER_STAGE)
                {
                    if (is_alpha)
                        Result = FOG_value.a;
                    else
                        Result = FOG_value; // TODO : Verify this!
                }
                else
                    Result = float4(FOG_value.rgb, DEFAULT_ALPHA);
            }
            break;
        default:
            Result = state.RegisterValues[register_index];
            break;
    }

    switch (input_mapping)
    {
        case PS_INPUTMAPPING_UNSIGNED_IDENTITY: // = 0x00L : y = max(0,x)       =  1*max(0,x) + 0.0
            Result = abs(Result);
            break;
        case PS_INPUTMAPPING_UNSIGNED_INVERT:   // = 0x20L : y = 1 - max(0,x)   = -1*max(0,x) + 1.0
            Result = 1 - Result;
            break;
        case PS_INPUTMAPPING_EXPAND_NORMAL:     // = 0x40L : y = 2*max(0,x) - 1 =  2*max(0,x) - 1.0
            Result = 2 * max(0, Result) - 1.0;
            break;
        case PS_INPUTMAPPING_EXPAND_NEGATE:     // = 0x60L : y = 1 - 2*max(0,x) = -2*max(0,x) + 1.0
            Result = -2 * max(0, Result) + 1.0;
            break;
        case PS_INPUTMAPPING_HALFBIAS_NORMAL:   // = 0x80L : y = max(0,x) - 1/2 =  1*max(0,x) - 0.5
            Result = 1 * max(0, Result) - 0.5;
            break;
        case PS_INPUTMAPPING_HALFBIAS_NEGATE:   // = 0xa0L : y = 1/2 - max(0,x) = -1*max(0,x) + 0.5
            Result = -1 * max(0, Result) + 0.5;
            break;
        case PS_INPUTMAPPING_SIGNED_IDENTITY:   // = 0xc0L : y = x              =  1*      x  + 0.0
            Result = Result;
            break;
        case PS_INPUTMAPPING_SIGNED_NEGATE:     // = 0xe0L : y = -x             = -1*      x  + 0.0
            Result = -Result;
            break;
    }

    bool use_channel_alpha = get_bit_4_from_value(reg_byte);
    if (use_channel_alpha)
        Result = Result.a; // PS_CHANNEL_ALPHA = 0x10; // used as RGB or ALPHA source
    else if (is_alpha)
        Result = Result.b; // PS_CHANNEL_BLUE = 0x00; // used as ALPHA source
    else
        Result = float4(Result.rgb, DEFAULT_ALPHA); // PS_CHANNEL_RGB: // = 0x00; // used as RGB source

    return Result;
}

void set_output_register_rgb(inout ps_state state, const uint1 reg_byte, const float4 value)
{
    uint register_index = mask_register(reg_byte);

	// Below are not writeable :
    switch (register_index)
    {
        case PS_REGISTER_DISCARD:
            return;
#ifdef AVOID_INVALID_ACCESSES
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
            return;
#endif
    }

    state.RegisterValues[register_index].rgb = value.rgb; // TODO : Fix error X3500: array reference cannot be used as an l-value; not natively addressable
}

void set_output_register_alpha(inout ps_state state, const uint1 reg_byte, const float4 value)
{
    uint register_index = mask_register(reg_byte);

	// Below are not writeable :
    switch (register_index)
    {
        case PS_REGISTER_DISCARD:
            return;
#ifdef AVOID_INVALID_ACCESSES
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
            return;
#endif
    }

    state.RegisterValues[register_index].a = value.a; // TODO : Fix error X3500: array reference cannot be used as an l-value; not natively addressable
}

void do_color_combiner_stage(inout ps_state state, const bool is_alpha)
{
	// Decode input register words :
    floor2_t AB_input_regs_word = is_alpha ? floor(D3DRS_PSALPHAINPUTS_Lower16[state.stage]) : floor(D3DRS_PSRGBINPUTS_Lower16[state.stage]);
    floor2_t CD_input_regs_word = is_alpha ? floor(D3DRS_PSALPHAINPUTS_Upper16[state.stage]) : floor(D3DRS_PSRGBINPUTS_Upper16[state.stage]);
    floor1_t A_input_reg_byte = get_byte_0_from_word(AB_input_regs_word);
    floor1_t B_input_reg_byte = get_byte_1_from_word(AB_input_regs_word);
    floor1_t C_input_reg_byte = get_byte_0_from_word(CD_input_regs_word);
    floor1_t D_input_reg_byte = get_byte_1_from_word(CD_input_regs_word);

	// Fetch input (source) register values :
    float4 A_value = get_input_register_as_float4(state, A_input_reg_byte, is_alpha); // A.k.a. 's0'
    float4 B_value = get_input_register_as_float4(state, B_input_reg_byte, is_alpha); // A.k.a. 's1'
    float4 C_value = get_input_register_as_float4(state, C_input_reg_byte, is_alpha); // A.k.a. 's2'
    float4 D_value = get_input_register_as_float4(state, D_input_reg_byte, is_alpha); // A.k.a. 's3'

	// Decode output words :
    floor2_t output_word0 = is_alpha ? floor(D3DRS_PSALPHAOUTPUTS_Lower16[state.stage]) : floor(D3DRS_PSRGBOUTPUTS_Lower16[state.stage]);
    floor2_t output_word1 = is_alpha ? floor(D3DRS_PSALPHAOUTPUTS_Upper16[state.stage]) : floor(D3DRS_PSRGBOUTPUTS_Upper16[state.stage]);
    floor1_t output_byte0 = get_byte_0_from_word(output_word0);
    floor1_t output_byte1 = get_byte_1_from_word(output_word0);
    floor1_t output_byte2 = get_byte_0_from_word(output_word1);

	// Decode output register numbers and flags :
    uint1 combiner_output_reg_AB = get_nibble_0_from_byte(output_byte0); // A.k.a. 'd0'
    uint1 combiner_output_reg_CD = get_nibble_1_from_byte(output_byte0); // A.k.a. 'd1'
    uint1 combiner_output_reg_AB_CD = get_nibble_0_from_byte(output_byte1); // A.k.a. 'd2'
    bool combiner_flag_CD_DotProduct = get_bit_4_from_value(output_byte1); // PS_COMBINEROUTPUT_CD_DOT_PRODUCT = 0x01L; // RGB only
    bool combiner_flag_AB_DotProduct = get_bit_5_from_value(output_byte1); // PS_COMBINEROUTPUT_AB_DOT_PRODUCT = 0x02L; // RGB only
    bool combiner_flag_ABCD_Mux = get_bit_6_from_value(output_byte1); // PS_COMBINEROUTPUT_AB_CD_MUX = 0x04L; // 3rd output is MUX(AB;CD) based on R0.a
	/* TODO : Fetch and interpret 3 output mapping bits
	static const uint PS_COMBINEROUTPUT_IDENTITY             = 0x00L; // y = x
	static const uint PS_COMBINEROUTPUT_BIAS                 = 0x08L; // y = x - 0.5
	static const uint PS_COMBINEROUTPUT_SHIFTLEFT_1          = 0x10L; // y = x*2
	static const uint PS_COMBINEROUTPUT_SHIFTLEFT_1_BIAS     = 0x18L; // y = (x - 0.5)*2
	static const uint PS_COMBINEROUTPUT_SHIFTLEFT_2          = 0x20L; // y = x*4
	// static const uint PS_COMBINEROUTPUT_SHIFTLEFT_2_BIAS  = 0x28L; // y = (x - 0.5)*4
	static const uint PS_COMBINEROUTPUT_SHIFTRIGHT_1         = 0x30L; // y = x/2
	// static const uint PS_COMBINEROUTPUT_SHIFTRIGHT_1_BIAS = 0x38L; // y = (x - 0.5)/2
	*/
    bool combiner_flag_CD_BlueToAlpha = get_bit_2_from_value(output_byte2); // PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA = 0x40L; // RGB only
    bool combiner_flag_AB_BlueToAlpha = get_bit_3_from_value(output_byte2); // PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA = 0x80L; // RGB only

	// AB_CD register output must be DISCARD if either AB_DOT_PRODUCT or CD_DOT_PRODUCT are set
    if (combiner_flag_AB_DotProduct || combiner_flag_CD_DotProduct)
    {
#ifdef AVOID_INVALID_ACCESSES
		// Fix-up, just to avoid issues (TODO : Investigate if real hardware does this too?)
        combiner_output_reg_AB_CD = PS_REGISTER_DISCARD;
#else
	#ifdef ASSERT_VALID_ACCESSES
		assert(combiner_output_reg_AB_CD == PS_REGISTER_DISCARD);
	#endif
#endif
    }

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
        float R0_a_value = state.RegisterValues[PS_REGISTER_R0].a;
        if (state.FlagMuxMsb)
            AB_CD_value = (R0_a_value > 0.5) ? CD_value : AB_value; // PS_COMBINERCOUNT_MUX_MSB
        else
            AB_CD_value = get_bit_0_from_value(R0_a_value) ? CD_value : AB_value; // PS_COMBINERCOUNT_MUX_LSB
    }
    else // PS_COMBINEROUTPUT_AB_CD_SUM
        AB_CD_value = AB_value + CD_value;

	// Store resulting values in output registers :
	// TODO : Handle output_mapping in set_output_register_*()
    // uint output_mapping = mask_outputmapping(output_byte1, output_byte2);
    // TODO : Apply combiner_flag_AB_BlueToAlpha and combiner_flag_CD_BlueToAlpha
    if (is_alpha)
    {
        set_output_register_alpha(state, combiner_output_reg_AB, AB_value);
        set_output_register_alpha(state, combiner_output_reg_CD, CD_value);
        set_output_register_alpha(state, combiner_output_reg_AB_CD, AB_CD_value);
    }
    else
    {
        set_output_register_rgb(state, combiner_output_reg_AB, AB_value);
        set_output_register_rgb(state, combiner_output_reg_CD, CD_value);
        set_output_register_rgb(state, combiner_output_reg_AB_CD, AB_CD_value);
    }
}

float4 do_final_combiner(inout ps_state state)
{
	// Fetch R0 register value (in final combiner used are input, otherwise returned) :
    float4 R0_value = state.RegisterValues[PS_REGISTER_R0];

	// Get the final combiner words :
    floor2_t AB_regs_word = floor(D3DRS_PSFINALCOMBINERINPUTSABCD_Lower16);
    floor2_t CD_regs_word = floor(D3DRS_PSFINALCOMBINERINPUTSABCD_Upper16);
    floor2_t EF_regs_word = floor(D3DRS_PSFINALCOMBINERINPUTSEFG_Lower16);
    floor2_t GS_regs_word = floor(D3DRS_PSFINALCOMBINERINPUTSEFG_Upper16);

	// Check if the final combiner needs to run at all :
    if (AB_regs_word > 0 || CD_regs_word > 0 || EF_regs_word > 0 || GS_regs_word > 0)
    {
		// Decompose the final combiner words into separate bytes :
        uint1 A_reg_byte = get_byte_0_from_word(AB_regs_word);
        uint1 B_reg_byte = get_byte_1_from_word(AB_regs_word);
        uint1 C_reg_byte = get_byte_0_from_word(CD_regs_word);
        uint1 D_reg_byte = get_byte_1_from_word(CD_regs_word);
        uint1 E_reg_byte = get_byte_0_from_word(EF_regs_word);
        uint1 F_reg_byte = get_byte_1_from_word(EF_regs_word);
        uint1 G_reg_byte = get_byte_0_from_word(GS_regs_word);
        uint1 setting_byte = get_byte_1_from_word(GS_regs_word);

		// Tell get_input_register_as_float4() we're now in the final combiner stage for E,F,G :
        state.stage = FINAL_COMBINER_STAGE;

		// Fetch the E and F register values :
        float3 E_value = (float3) get_input_register_as_float4(state, E_reg_byte);
        float3 F_value = (float3) get_input_register_as_float4(state, F_reg_byte);

		// Calculate the product of E*F and set the result in it's dedicated register :
        float3 EF_PROD_value = E_value * F_value;
        state.RegisterValues[PS_REGISTER_EF_PROD] = float4(EF_PROD_value, DEFAULT_ALPHA);

		// Do optional complement operation on R0 register value :
        if (get_bit_5_from_value(setting_byte)) // PS_FINALCOMBINERSETTING_COMPLEMENT_R0
        {
			// unsigned invert mapping  (1 - r0) is used as an input to the sum rather than r0
            R0_value = 1 - R0_value;
        }

		// Fetch V1 register value (including an optional complement operation) :
        float3 V1_value;
        V1_value = (float3) state.RegisterValues[PS_REGISTER_V1];
        if (get_bit_6_from_value(setting_byte)) // PS_FINALCOMBINERSETTING_COMPLEMENT_V1
        {
			// unsigned invert mapping  (1 - v1) is used as an input to the sum rather than v1
            V1_value = 1 - V1_value;
        }

		// Calculate V1 + R0 value (including an optional clamping operation) :
        float3 V1R0_SUM_value = V1_value + R0_value.rgb;
        if (get_bit_7_from_value(setting_byte)) // PS_FINALCOMBINERSETTING_CLAMP_SUM
        {
			// V1+R0 sum clamped to [0,1]
            V1R0_SUM_value = clamp(V1R0_SUM_value, 0, 1);
        }

		// Set the result in it's dedicated register
        state.RegisterValues[PS_REGISTER_V1R0_SUM] = float4(V1R0_SUM_value, DEFAULT_ALPHA);

		// Tell get_input_register_as_float4() we're now in the final combiner stage for A,B,C,D :
        state.stage = FINAL_COMBINER_STAGE_ABCD;
		// Fetch the A, B, C and D register values :
        float4 A_value = get_input_register_as_float4(state, A_reg_byte);
        float4 B_value = get_input_register_as_float4(state, B_reg_byte);
        float4 C_value = get_input_register_as_float4(state, C_reg_byte);
        float4 D_value = get_input_register_as_float4(state, D_reg_byte);

		// Calculate output RGB = a*b + (1-a)*c + d :
        float4 FinalOutput;
        FinalOutput.rgb = lerp(B_value.rgb, C_value.rgb, A_value.rgb) + D_value.rgb;

		// Fetch the output alpha channel from the G register value :
        float4 G_value = get_input_register_as_float4(state, G_reg_byte, /*is_alpha=*/true); // TODO : Is is_alpha correct?
        FinalOutput.a = G_value.a;

        return FinalOutput;
    }
    else
        return R0_value;
}

float4 main() : SV_TARGET
{
    ps_state state = (ps_state) 0; // Clearing like this avoids error X3508: 'do_color_combiner_stage': output parameter 'state' not completely initialized
	// Note, the above also sets state.RegisterValues[PS_REGISTER_ZERO] to 0.0

	// TODO : Is the following initialization of R0 and FOG correct?
    state.RegisterValues[PS_REGISTER_R0] = float4(1.0f, 1.0f, 1.0f, DEFAULT_ALPHA);
    state.RegisterValues[PS_REGISTER_FOG] = float4(1.0f, 1.0f, 1.0f, DEFAULT_ALPHA);    

	// Decode the global combiner count and flags :
    floor2_t CombinerCountWord = floor(D3DRS_PSCOMBINERCOUNT_Lower16);
    state.FlagUniqueC1 = get_bit_16_from_value(CombinerCountWord); // PS_COMBINERCOUNT_UNIQUE_C1 = 0x0100L // c1 unique in each stage
    state.FlagUniqueC0 = get_bit_12_from_value(CombinerCountWord); // PS_COMBINERCOUNT_UNIQUE_C0 = 0x0010L, // c0 unique in each stage
    state.FlagMuxMsb = get_bit_8_from_value(CombinerCountWord); // PS_COMBINERCOUNT_MUX_MSB = 0x0001L, // mux on r0.a msb

    uint num_stages = get_byte_0_from_word(CombinerCountWord);
#ifdef ASSERT_VALID_ACCESSES
    assert(num_stages >= 1);
    assert(num_stages <= 8);
#endif

	// TODO : Implement texture fetch (preferrably in this shader, including conversion of X_D3DFMT_P8 and other unsupported textures formats)

	// Loop over at most 8 stages (start counting from zero)
    for (uint stage = 0; stage < MAX_COMBINER_STAGE_COUNT; stage++)
    {
		// Help loop-unrolling (which doesn't work with an uncertain limit) by breaking ourselves
        if (stage >= num_stages)
        {
            break;
        }

        state.stage = stage; // tell get_input_register_as_float4() the currently active stage
        do_color_combiner_stage(state, /*rgb*/false);
        do_color_combiner_stage(state, /*is_alpha=*/true);
    }

    return do_final_combiner(state);
}
