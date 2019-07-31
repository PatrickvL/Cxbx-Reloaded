//
// Utility functions
//
void assert(bool condition)
{
	// TODO
}

bool get_bit_0_from_value(const int value)
{
    return fmod(value, 2) == 1;
}

bool get_bit_1_from_value(const int value)
{
    return fmod(value, 4) >= 2;
}

bool get_bit_2_from_value(const int value)
{
    return fmod(value, 8) >= 4;
}

bool get_bit_3_from_value(const int value)
{
    return fmod(value, 16) >= 8;
}

bool get_bit_4_from_value(const int value)
{
    return fmod(value, 32) >= 16;
}

bool get_bit_5_from_value(const int value)
{
    return fmod(value, 64) >= 32;
}

bool get_bit_6_from_value(const int value)
{
    return fmod(value, 128) >= 64;
}

bool get_bit_7_from_value(const int value)
{
    return fmod(value, 256) >= 128;
}

bool get_bit_8_from_value(const int value)
{
    return fmod(value, 512) >= 256;
}

bool get_bit_9_from_value(const int value)
{
    return fmod(value, 1024) >= 512;
}

bool get_bit_10_from_value(const int value)
{
    return fmod(value, 2048) >= 1024;
}

bool get_bit_11_from_value(const int value)
{
    return fmod(value, 4096) >= 2048;
}

bool get_bit_12_from_value(const int value)
{
    return fmod(value, 8192) >= 4096;
}

bool get_bit_13_from_value(const int value)
{
    return fmod(value, 16384) >= 8192;
}

bool get_bit_14_from_value(const int value)
{
    return fmod(value, 32768) >= 16384;
}

bool get_bit_15_from_value(const int value)
{
    return fmod(value, 65536) >= 32768;
}

bool get_bit_16_from_value(const int value)
{
    return fmod(value, 131072) >= 65536;
}

bool get_bit_17_from_value(const int value)
{
    return fmod(value, 262144) >= 131072;
}

bool get_bit_18_from_value(const int value)
{
    return fmod(value, 524288) >= 262144;
}

bool get_bit_19_from_value(const int value)
{
    return fmod(value, 1048576) >= 524288;
}

bool get_bit_20_from_value(const int value)
{
    return fmod(value, 2097152) >= 1048576;
}

bool get_bit_21_from_value(const int value)
{
    return fmod(value, 4194304) >= 2097152;
}

bool get_bit_22_from_value(const int value)
{
    return fmod(value, 8388608) >= 4194304;
}

bool get_bit_from_value(const int value, const int bit_nr) // valid bit_nr : 0..22
{
    assert(bit_nr >= 0);
    assert(bit_nr <= 22);

	// Based on http://theinstructionlimit.com/encoding-boolean-flags-into-a-float-in-hlsl
	// and https://en.wikipedia.org/wiki/Single-precision_floating-point_format
    switch (bit_nr)
    {
        case 0:
            return get_bit_0_from_value(value);
        case 1:
            return get_bit_1_from_value(value);
        case 2:
            return get_bit_2_from_value(value);
        case 3:
            return get_bit_3_from_value(value);
        case 4:
            return get_bit_4_from_value(value);
        case 5:
            return get_bit_5_from_value(value);
        case 6:
            return get_bit_6_from_value(value);
        case 7:
            return get_bit_7_from_value(value);
        case 8:
            return get_bit_8_from_value(value);
        case 9:
            return get_bit_9_from_value(value);
        case 10:
            return get_bit_10_from_value(value);
        case 11:
            return get_bit_11_from_value(value);
        case 12:
            return get_bit_12_from_value(value);
        case 13:
            return get_bit_13_from_value(value);
        case 14:
            return get_bit_14_from_value(value);
        case 15:
            return get_bit_15_from_value(value);
        case 16:
            return get_bit_16_from_value(value);
        case 17:
            return get_bit_17_from_value(value);
        case 18:
            return get_bit_18_from_value(value);
        case 19:
            return get_bit_19_from_value(value);
        case 20:
            return get_bit_20_from_value(value);
        case 21:
            return get_bit_21_from_value(value);
        case 22:
            return get_bit_22_from_value(value);
        default:
            assert(false);
            return 0;
    }
}

//
// Xbox Pixel Shader Constants
//

// Xbox Pixel Shader definition (a dword array) D3D__RenderState member indices :
static const dword X_D3DRS_PSALPHAINPUTS0 = 0u;
static const dword X_D3DRS_PSALPHAINPUTS1 = 1u;
static const dword X_D3DRS_PSALPHAINPUTS2 = 2u;
static const dword X_D3DRS_PSALPHAINPUTS3 = 3u;
static const dword X_D3DRS_PSALPHAINPUTS4 = 4u;
static const dword X_D3DRS_PSALPHAINPUTS5 = 5u;
static const dword X_D3DRS_PSALPHAINPUTS6 = 6u;
static const dword X_D3DRS_PSALPHAINPUTS7 = 7u;
static const dword X_D3DRS_PSFINALCOMBINERINPUTSABCD = 8u;
static const dword X_D3DRS_PSFINALCOMBINERINPUTSEFG = 9u;
static const dword X_D3DRS_PSCONSTANT0_0 = 10u;
static const dword X_D3DRS_PSCONSTANT0_1 = 11u;
static const dword X_D3DRS_PSCONSTANT0_2 = 12u;
static const dword X_D3DRS_PSCONSTANT0_3 = 13u;
static const dword X_D3DRS_PSCONSTANT0_4 = 14u;
static const dword X_D3DRS_PSCONSTANT0_5 = 15u;
static const dword X_D3DRS_PSCONSTANT0_6 = 16u;
static const dword X_D3DRS_PSCONSTANT0_7 = 17u;
static const dword X_D3DRS_PSCONSTANT1_0 = 18u;
static const dword X_D3DRS_PSCONSTANT1_1 = 19u;
static const dword X_D3DRS_PSCONSTANT1_2 = 20u;
static const dword X_D3DRS_PSCONSTANT1_3 = 21u;
static const dword X_D3DRS_PSCONSTANT1_4 = 22u;
static const dword X_D3DRS_PSCONSTANT1_5 = 23u;
static const dword X_D3DRS_PSCONSTANT1_6 = 24u;
static const dword X_D3DRS_PSCONSTANT1_7 = 25u;
static const dword X_D3DRS_PSALPHAOUTPUTS0 = 26u;
static const dword X_D3DRS_PSALPHAOUTPUTS1 = 27u;
static const dword X_D3DRS_PSALPHAOUTPUTS2 = 28u;
static const dword X_D3DRS_PSALPHAOUTPUTS3 = 29u;
static const dword X_D3DRS_PSALPHAOUTPUTS4 = 30u;
static const dword X_D3DRS_PSALPHAOUTPUTS5 = 31u;
static const dword X_D3DRS_PSALPHAOUTPUTS6 = 32u;
static const dword X_D3DRS_PSALPHAOUTPUTS7 = 33u;
static const dword X_D3DRS_PSRGBINPUTS0 = 34u;
static const dword X_D3DRS_PSRGBINPUTS1 = 35u;
static const dword X_D3DRS_PSRGBINPUTS2 = 36u;
static const dword X_D3DRS_PSRGBINPUTS3 = 37u;
static const dword X_D3DRS_PSRGBINPUTS4 = 38u;
static const dword X_D3DRS_PSRGBINPUTS5 = 39u;
static const dword X_D3DRS_PSRGBINPUTS6 = 40u;
static const dword X_D3DRS_PSRGBINPUTS7 = 41u;
static const dword X_D3DRS_PSCOMPAREMODE = 42u;
static const dword X_D3DRS_PSFINALCOMBINERCONSTANT0 = 43u;
static const dword X_D3DRS_PSFINALCOMBINERCONSTANT1 = 44u;
static const dword X_D3DRS_PSRGBOUTPUTS0 = 45u;
static const dword X_D3DRS_PSRGBOUTPUTS1 = 46u;
static const dword X_D3DRS_PSRGBOUTPUTS2 = 47u;
static const dword X_D3DRS_PSRGBOUTPUTS3 = 48u;
static const dword X_D3DRS_PSRGBOUTPUTS4 = 49u;
static const dword X_D3DRS_PSRGBOUTPUTS5 = 50u;
static const dword X_D3DRS_PSRGBOUTPUTS6 = 51u;
static const dword X_D3DRS_PSRGBOUTPUTS7 = 52u;
static const dword X_D3DRS_PSCOMBINERCOUNT = 53u;
static const dword X_D3DRS_PS_RESERVED = 54u; // Dxbx note : This takes the slot of X_D3DPIXELSHADERDEF.PSTextureModesu; set by D3DDevice_SetRenderState_LogicOp?
static const dword X_D3DRS_PSDOTMAPPING = 55u;
static const dword X_D3DRS_PSINPUTTEXTURE = 56u;

// The input can have the following mappings applied :
//
// PS_INPUTMAPPING_UNSIGNED_IDENTITY : y = max(0,x)       =  1*max(0,x) + 0.0
// PS_INPUTMAPPING_UNSIGNED_INVERT   : y = 1 - max(0,x)   = -1*max(0,x) + 1.0
// PS_INPUTMAPPING_EXPAND_NORMAL     : y = 2*max(0,x) - 1 =  2*max(0,x) - 1.0
// PS_INPUTMAPPING_EXPAND_NEGATE     : y = 1 - 2*max(0,x) = -2*max(0,x) + 1.0
// PS_INPUTMAPPING_HALFBIAS_NORMAL   : y = max(0,x) - 1/2 =  1*max(0,x) - 0.5
// PS_INPUTMAPPING_HALFBIAS_NEGATE   : y = 1/2 - max(0,x) = -1*max(0,x) + 0.5
// PS_INPUTMAPPING_SIGNED_IDENTITY   : y = x              =  1*      x  + 0.0
// PS_INPUTMAPPING_SIGNED_NEGATE     : y = -x             = -1*      x  + 0.0
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
static const dword PS_INPUTMAPPING_UNSIGNED_IDENTITY = 0x00L; // max(0;x)         OK for final combiner: y = abs(x)
static const dword PS_INPUTMAPPING_UNSIGNED_INVERT =   0x20L; // 1 - max(0;x)     OK for final combiner: y = 1 - x
static const dword PS_INPUTMAPPING_EXPAND_NORMAL =     0x40L; // 2*max(0;x) - 1   invalid for final combiner
static const dword PS_INPUTMAPPING_EXPAND_NEGATE =     0x60L; // 1 - 2*max(0;x)   invalid for final combiner
static const dword PS_INPUTMAPPING_HALFBIAS_NORMAL =   0x80L; // max(0;x) - 1/2   invalid for final combiner
static const dword PS_INPUTMAPPING_HALFBIAS_NEGATE =   0xa0L; // 1/2 - max(0;x)   invalid for final combiner
static const dword PS_INPUTMAPPING_SIGNED_IDENTITY =   0xc0L; // x                invalid for final combiner
static const dword PS_INPUTMAPPING_SIGNED_NEGATE =     0xe0L; // -x               invalid for final combiner

// Pixel Shader Registers
static const dword PS_REGISTER_ZERO =              0x00L; // r
static const dword PS_REGISTER_DISCARD =           0x00L; // w
static const dword PS_REGISTER_C0 =                0x01L; // r
static const dword PS_REGISTER_C1 =                0x02L; // r
static const dword PS_REGISTER_FOG =               0x03L; // r
static const dword PS_REGISTER_V0 =                0x04L; // r/w
static const dword PS_REGISTER_V1 =                0x05L; // r/w
static const dword PS_REGISTER_T0 =                0x08L; // r/w
static const dword PS_REGISTER_T1 =                0x09L; // r/w
static const dword PS_REGISTER_T2 =                0x0aL; // r/w
static const dword PS_REGISTER_T3 =                0x0bL; // r/w
static const dword PS_REGISTER_R0 =                0x0cL; // r/w
static const dword PS_REGISTER_R1 =                0x0dL; // r/w
static const dword PS_REGISTER_V1R0_SUM =          0x0eL; // r
static const dword PS_REGISTER_EF_PROD =           0x0fL; // r
static const dword PS_REGISTER_ONE =               PS_REGISTER_ZERO | PS_INPUTMAPPING_UNSIGNED_INVERT; // 0x20 OK for final combiner
static const dword PS_REGISTER_NEGATIVE_ONE =      PS_REGISTER_ZERO | PS_INPUTMAPPING_EXPAND_NORMAL;   // 0x40 invalid for final combiner
static const dword PS_REGISTER_ONE_HALF =          PS_REGISTER_ZERO | PS_INPUTMAPPING_HALFBIAS_NEGATE; // 0xa0 invalid for final combiner
static const dword PS_REGISTER_NEGATIVE_ONE_HALF = PS_REGISTER_ZERO | PS_INPUTMAPPING_HALFBIAS_NORMAL; // 0x80 invalid for final combiner
static const dword PS_REGISTER_CXBX_PROD =         PS_REGISTER_ZERO | PS_INPUTMAPPING_SIGNED_IDENTITY; // Cxbx internal use

static const dword PS_CHANNEL_RGB =   0x00; // used as RGB source
static const dword PS_CHANNEL_BLUE =  0x00; // used as ALPHA source
static const dword PS_CHANNEL_ALPHA = 0x10; // used as RGB or ALPHA source

static const dword PS_FINALCOMBINERSETTING_CLAMP_SUM =     0x80; // V1+R0 sum clamped to [0;1]
static const dword PS_FINALCOMBINERSETTING_COMPLEMENT_V1 = 0x40; // unsigned invert mapping  (1 - v1) is used as an input to the sum rather than v1
static const dword PS_FINALCOMBINERSETTING_COMPLEMENT_R0 = 0x20; // unsigned invert mapping  (1 - r0) is used as an input to the sum rather than r0

//
// Xbox Pixel Shader definition (a 32 bit dword/float array) buffer, taken from "D3D__RenderState"
//
extern const Buffer<float> g_XboxPSDef_Full32;
extern const Buffer<float> g_XboxPSDef_Upper16; // TODO : Either Upper16 or Lower16 could be read from g_XboxPSDef_Full32
extern const Buffer<float> g_XboxPSDef_Lower16;

float get_render_state(const int render_state)
{
    return g_XboxPSDef_Full32.Load(render_state);
}

int get_render_state_low_word(const int render_state)
{
    return g_XboxPSDef_Lower16.Load(render_state);
}

int get_render_state_high_word(const int render_state)
{
    return g_XboxPSDef_Upper16.Load(render_state);
}

// Retrieves a field sized 1 (up to 23) bits from the specified render state
// Note : Shader Model 3 lacks support for integer inputs, so we have to interpret the input as floats here.
// Also, since 32-bit floats have only 23 fraction bits, we split the input in a low and high 16-bit part.
//
// start_bit starts at 0 for least significant bit, and ends on 31 for most significant bit
// bit_count must ly between 1 and 32-start_bit (both inclusive)
int get_render_state_bitfield(const dword render_state, const int start_bit, int bit_count)
{
    assert(start_bit >= 0);
    assert(start_bit <= 31);
    assert(bit_count > 0);
    assert(start_bit + bit_count <= 32);

    int Result = 0;

    int current_bit = start_bit;
    int value = get_render_state_low_word(render_state);

	while (current_bit < 16) {
        Result = Result + Result + get_bit_from_value(value, current_bit++);
		if (bit_count-- == 0)
            return Result;
    }

    value = get_render_state_high_word(render_state);
	while (bit_count-- > 0) {
		Result = Result + Result + get_bit_from_value(value, current_bit++);
	}

	return Result;
}

//
// Our Xbox Pixel Shader Interpreter state
//
typedef struct
{
	// Flags, used in get_input_register_as_float4() :
    int stage;         // Currently active stage
    bool is_fc;        // Set when final combiner stage is active
    bool FlagMuxMsb;   // Mux on r0.a lsb or msb, 0 = PS_COMBINERCOUNT_MUX_LSB, 1 = PS_COMBINERCOUNT_MUX_MSB
    bool FlagUniqueC0; // C[0] unique in each stage, 0 = PS_COMBINERCOUNT_SAME_C0, 1 = PS_COMBINERCOUNT_UNIQUE_C0
    bool FlagUniqueC1; // C[1] unique in each stage, 0 = PS_COMBINERCOUNT_SAME_C1, 1 = PS_COMBINERCOUNT_UNIQUE_C1

	// Xbox NV2A pixel shader register values :
#ifdef EXPLICIT_REGISTERS
    float4 C[8]; // Constant registers, read-only
    float4 R[2]; // Temporary registers, read/write (R[0].a (R0_ALPHA) is initialized to T[0].s (T0_ALPHA) in stage 0); Final result is in R[0]
    float4 T[4]; // Texture registers, initially sampled or (0,0,0,1), after texture-addressing they're read/write just like temporary registers
    float4 V[2]; // Vertex color registers, read/write (V[0] = diffuse, V[1] = specular
    float4 FOG;  // Note, FOG ALPHA is only available in final combiner
	// V1R0_SUM and EF_PROD are only available in final combiner (A,B,C,D inputs only)
	// V1R0_SUM_ALPHA and EF_PROD_ALPHA are not available, hence the float3 type here :
    float3 V1R0_SUM;	
    float3 EF_PROD;
#else
    float4 RegisterValues[16];
#endif

	// The final output value :
    float4 RGBout;
} ps_state;

float4 get_input_register_as_float4(inout ps_state state, int reg_bits, bool is_alpha = false)
{
    float4 Result = 0;

    int register_index = reg_bits & 0xF;

    if (state.is_fc) {
		// Below are only valid for final combiner :
		// TODO : Limit to A,B,C,D inputs
        switch (register_index) {
            case PS_REGISTER_V1R0_SUM:
			case PS_REGISTER_EF_PROD:
		        assert(false);
        }
	}

    switch (register_index) {
        case PS_REGISTER_ZERO:
            Result = 0;
            break;
		case PS_REGISTER_C0:
            if (state.is_fc)
                Result = get_render_state(X_D3DRS_PSFINALCOMBINERCONSTANT0);
			else
				if (state.FlagUniqueC0)
					Result = get_render_state(X_D3DRS_PSCONSTANT0_0 + state.stage);
				else
#ifdef EXPLICIT_REGISTERS
					Result = state.C[0];
#else
                Result = state.RegisterValues[PS_REGISTER_C0];
#endif
            break;
		case PS_REGISTER_C1:
            if (state.is_fc)
                Result = get_render_state(X_D3DRS_PSFINALCOMBINERCONSTANT1);
            else
				if (state.FlagUniqueC1)
					Result = get_render_state(X_D3DRS_PSCONSTANT1_0 + state.stage);
				else
#ifdef EXPLICIT_REGISTERS
					Result = state.C[1];
#else
                Result = state.RegisterValues[PS_REGISTER_C1];
#endif
            break;
		case PS_REGISTER_FOG: {
#ifdef EXPLICIT_REGISTERS
            float4 Fog = stage.FOG;
#else
            float4 Fog = state.RegisterValues[PS_REGISTER_FOG];
#endif
            if (state.is_fc)
                if (is_alpha)
                    Result = Fog.a;
                else
                    Result = Fog; // TODO : Verify this!
            else
                Result = Fog.rgb;
			}
            break;
#ifdef EXPLICIT_REGISTERS
        case PS_REGISTER_V0:
            Result = state.V[0];
            break;
        case PS_REGISTER_V1:
            Result = state.V[1];
            break;
        case PS_REGISTER_T0:
            Result = state.T[0];
            break;
        case PS_REGISTER_T1:
            Result = state.T[1];
            break;
        case PS_REGISTER_T2:
            Result = state.T[2];
            break;
        case PS_REGISTER_T3:
            Result = state.T[3];
            break;
        case PS_REGISTER_R0:
            Result = state.R[0];
            break;
        case PS_REGISTER_R1:
            Result = state.R[1];
            break;
		// Below are only valid for final combiner :
		// TODO : Limit to A,B,C,D inputs
        case PS_REGISTER_V1R0_SUM:
            Result = state.V1R0_SUM;
            break;
        case PS_REGISTER_EF_PROD:
			Result = state.EF_PROD;
            break;
#else
		default:
            Result = state.RegisterValues[register_index];
            break;
#endif
    }

    switch (reg_bits & 0xE0)
    {
        case PS_INPUTMAPPING_UNSIGNED_IDENTITY: // = 0x00L : y = max(0,x)       =  1*max(0,x) + 0.0
            Result = abs(Result);
		case PS_INPUTMAPPING_UNSIGNED_INVERT:   // = 0x20L : y = 1 - max(0,x)   = -1*max(0,x) + 1.0
            Result = 1 - Result;
            break;
		// Below are invalid for final combiner :
		case PS_INPUTMAPPING_EXPAND_NORMAL:     // = 0x40L : y = 2*max(0,x) - 1 =  2*max(0,x) - 1.0
            if (state.is_fc)
                assert(false);
			else
				Result =  2 * max(0, Result) - 1.0;
            break;
		case PS_INPUTMAPPING_EXPAND_NEGATE:     // = 0x60L : y = 1 - 2*max(0,x) = -2*max(0,x) + 1.0
            if (state.is_fc)
                assert(false);
            else
				Result = -2 * max(0, Result) + 1.0;
            break;
		case PS_INPUTMAPPING_HALFBIAS_NORMAL:   // = 0x80L : y = max(0,x) - 1/2 =  1*max(0,x) - 0.5
            if (state.is_fc)
                assert(false);
            else
				Result =  1 * max(0, Result) - 0.5;
            break;
		case PS_INPUTMAPPING_HALFBIAS_NEGATE:   // = 0xa0L : y = 1/2 - max(0,x) = -1*max(0,x) + 0.5
            if (state.is_fc)
                assert(false);
            else
				Result = -1 * max(0, Result) + 0.5;
            break;
		case PS_INPUTMAPPING_SIGNED_IDENTITY:   // = 0xc0L : y = x              =  1*      x  + 0.0
            if (state.is_fc)
                assert(false);
            else
	            Result =             Result;
            break;
        case PS_INPUTMAPPING_SIGNED_NEGATE:     // = 0xe0L : y = -x             = -1*      x  + 0.0
            if (state.is_fc)
                assert(false);
            else
		        Result = -           Result;
            break;
	}

    switch (reg_bits & 0x10)
    {
		case PS_CHANNEL_RGB:// = 0x00; // used as RGB source
            if (is_alpha) // PS_CHANNEL_BLUE = 0x00; // used as ALPHA source
                Result = Result.b;
			else
                Result = Result.rgb;
			break;
		case PS_CHANNEL_ALPHA: // = 0x10; // used as RGB or ALPHA source
            Result = Result.a;
            break;
    }

    return Result;
}

void do_color_combiner_stage(inout ps_state state, const bool is_alpha)
{
    int input_low = get_render_state_low_word(is_alpha ? X_D3DRS_PSALPHAINPUTS0 : X_D3DRS_PSRGBINPUTS0 + state.stage);
    int input_high = get_render_state_high_word(is_alpha ? X_D3DRS_PSALPHAINPUTS0 : X_D3DRS_PSRGBINPUTS0 + state.stage);
    int output_low = get_render_state_low_word(is_alpha ? X_D3DRS_PSALPHAOUTPUTS0 : X_D3DRS_PSRGBOUTPUTS0 + state.stage);
    int output_high = get_render_state_high_word(is_alpha ? X_D3DRS_PSALPHAOUTPUTS0 : X_D3DRS_PSRGBOUTPUTS0 + state.stage);

    float4 input1_value = get_input_register_as_float4(state, input_low & 0xFF, is_alpha);
    float4 input2_value = get_input_register_as_float4(state, input_low >> 8, is_alpha);
	// TODO
}

void do_final_combiner(inout ps_state state)
{
    int AB = get_render_state_low_word(X_D3DRS_PSFINALCOMBINERINPUTSABCD);
    int CD = get_render_state_high_word(X_D3DRS_PSFINALCOMBINERINPUTSABCD);
    int EF = get_render_state_low_word(X_D3DRS_PSFINALCOMBINERINPUTSEFG);
    int GS = get_render_state_high_word(X_D3DRS_PSFINALCOMBINERINPUTSEFG); // S = settings
    if (AB || CD || EF || GS) {
        state.is_fc = true; // tell get_input_register_as_float4() we're now in the final combiner stage

        float4 val_e = get_input_register_as_float4(state, EF & 0xFF);
        float4 val_f = get_input_register_as_float4(state, EF >> 8);
#ifdef EXPLICIT_REGISTERS
        state.EF_PROD = val_e * val_f;
#else
        state.RegisterValues[PS_REGISTER_EF_PROD] = val_e * val_f;
#endif

        int Setting = GS / 256; // TODO : Fix warning X3556: integer divides may be much slower, try using uints if possible.

#ifdef EXPLICIT_REGISTERS
        float4 V1 = state.V[1];
#else
        float4 V1 = state.RegisterValues[PS_REGISTER_V1];
#endif
        if (get_bit_6_from_value(Setting)) { // PS_FINALCOMBINERSETTING_COMPLEMENT_V1= 0x40
			// unsigned invert mapping  (1 - v1) is used as an input to the sum rather than v1
            V1 = 1 - V1;
        }

#ifdef EXPLICIT_REGISTERS
        float4 R0 = state.R[0];
#else
        float4 R0 = state.RegisterValues[PS_REGISTER_R0];
#endif
        if (get_bit_5_from_value(Setting)) { // PS_FINALCOMBINERSETTING_COMPLEMENT_R0= 0x20
			// unsigned invert mapping  (1 - r0) is used as an input to the sum rather than r0
            R0 = 1 - R0;
        }

        float4 V1R0_SUM = V1 + R0;
        if (get_bit_7_from_value(Setting)) { // PS_FINALCOMBINERSETTING_CLAMP_SUM=     0x80
			// V1+R0 sum clamped to [0,1]
            V1R0_SUM = clamp(V1R0_SUM, 0, 1);
        }
#ifdef EXPLICIT_REGISTERS
        state.V1R0_SUM = V1R0_SUM;
#else
        state.RegisterValues[PS_REGISTER_V1R0_SUM] = V1R0_SUM;
#endif

        float4 val_a = get_input_register_as_float4(state, AB & 0xFF);
        float4 val_b = get_input_register_as_float4(state, AB >> 8);
        float4 val_c = get_input_register_as_float4(state, CD & 0xFF);
        float4 val_d = get_input_register_as_float4(state, CD >> 8);
        state.RGBout.rgb = lerp(val_b.rgb, val_c.rgb, val_a.rgb) + val_d.rgb;

        float4 val_g = get_input_register_as_float4(state, GS & 0xFF, /*is_alpha=*/true);
        state.RGBout.a = val_g.a;
    }
}

float4 main() : SV_TARGET
{
    ps_state state = (ps_state) 0; // Clearing like this avoids error X3508: 'do_color_combiner_stage': output parameter 'state' not completely initialized
    state.RGBout = float4(1.0f, 1.0f, 1.0f, 1.0f);

	// Decode the global combiner count and flags :
    state.FlagUniqueC1 = get_render_state_bitfield(X_D3DRS_PSCOMBINERCOUNT, 16, 1); // PS_COMBINERCOUNT_UNIQUE_C1
    state.FlagUniqueC0 = get_render_state_bitfield(X_D3DRS_PSCOMBINERCOUNT, 12, 1); // PS_COMBINERCOUNT_UNIQUE_C0
    state.FlagMuxMsb = get_render_state_bitfield(X_D3DRS_PSCOMBINERCOUNT, 8, 1); // PS_COMBINERCOUNT_MUX_MSB
    int num_stages = get_render_state_bitfield(X_D3DRS_PSCOMBINERCOUNT, 0, 8);    
    assert(num_stages >= 1);
    assert(num_stages <= 8);

	// TODO : Implement texture fetch (use a separate shader?)

	// TODO : Should we set C[0] and C[1] here, and with what value?

    state.is_fc = false;
    for (int i = 0; i < 8; i++) { // loop over at most 8 stages (start counting from zero)
		// Help loop-unrolling (which doesn't work with an uncertain limit) by breaking ourselves
        if (i >= num_stages) {
            break;
        }

        state.stage = i; // tell get_input_register_as_float4() the currently active stage
        do_color_combiner_stage(state, /*rgb*/false);
        do_color_combiner_stage(state, /*is_alpha=*/true);
    }

    do_final_combiner(state);

    return state.RGBout;
}
