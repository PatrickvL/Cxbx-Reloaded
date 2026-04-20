#include "CxbxVertexShaderCommon.hlsli"
#ifdef CXBX_IA_BYPASS
#include "CxbxVertexFetch.hlsli"
#endif

#define X_D3DSCM_CORRECTION                 96 // Add 96 to arrive at the range 0..191 (instead of -96..95)
#define X_D3DVS_CONSTREG_COUNT              192

// Xbox constant registers
uniform float4 C[X_D3DVS_CONSTREG_COUNT] : register(c0);

#include "CxbxScreenspaceTransform.hlsli"

// Overloaded casts, assuring all inputs are treated as float4
float4 _tof4(float  src) { return float4(src, src, src, src); }
float4 _tof4(float2 src) { return src.xyyy; }
float4 _tof4(float3 src) { return src.xyzz; }
float4 _tof4(float4 src) { return src; }
float4 _ssss(float s)    { return float4(s, s, s, s); } // a scalar output replicated across a 4-component vector
#define _scalar(src) _tof4(src).x /* a scalar input */

float4 c(int register_number)
{
	// Map Xbox [-96, 95] to Host [0, 191]
	// Account for Xbox's negative constant indexes
    register_number += X_D3DSCM_CORRECTION;

    // Like Xbox, out-of-range indices are guaranteed to be zero in HLSL
    // so no need to bounds check negative numbers
    // if (register_number < 0)
    //    return 0;

    // If the index is too large, set it to -1 so 0 is returned
    // Note: returning 0 directly requires many more instructions
    if (register_number >= X_D3DVS_CONSTREG_COUNT) // X_D3DVS_CONSTREG_COUNT
        register_number = -1;

    return C[register_number];
}

// Due to rounding differences with the Xbox (and increased precision on PC?)
// some titles produce values just below the threshold of the next integer.
// We can add a small bias to make sure it's bumped over the threshold
// Test Case: Azurik (divides indexes 755, then scales them back in the vertex shader)
#define BIAS 0.001
// NOTE : Was 0.0001, unlike xqemu

// 2.14.1.11  Vertex Program Floating Point Requirements
// The floor operations used by the ARL and EXP instructions must
// operate identically.  Specifically, the EXP instruction's floor(t.x)
// intermediate result must exactly match the integer stored in the
// address register by the ARL instruction.
float x_floor(float src)
{
	return floor(src + BIAS);
}

// https://xboxdevwiki.net/NV2A/Vertex_Shader
// https://www.khronos.org/registry/OpenGL/extensions/NV/NV_vertex_program.txt
// https://www.khronos.org/registry/OpenGL/extensions/NV/NV_vertex_program1_1.txt

// Functions for MAC ('Multiply And Accumulate') opcodes

// 2.14.1.10.1  ARL: Address Register Load
// The address register should be floored
#define x_arl(dest, mask, src0) dest.mask = x_floor(_tof4(src0).x).mask

// 2.14.1.10.2  MOV: Move
#define x_mov(dest, mask, src0) dest.mask = (_tof4(src0)).mask

// 2.14.1.10.3  MUL: Multiply
#define x_mul(dest, mask, src0, src1) dest.mask = (_tof4(src0) * _tof4(src1)).mask

// 2.14.1.10.4  ADD: Add
#define x_add(dest, mask, src0, src1) dest.mask = (_tof4(src0) + _tof4(src1)).mask

// 2.14.1.10.5  MAD: Multiply and Add
#define x_mad(dest, mask, src0, src1, src2) dest.mask = (_tof4(src0) * _tof4(src1) + _tof4(src2)).mask

// 2.14.1.10.8  DP3: Three-Component Dot Product
#define x_dp3(dest, mask, src0, src1) dest.mask = _ssss(dot(_tof4(src0).xyz, _tof4(src1).xyz)).mask

//  2.14.1.10.9  DP4: Four-Component Dot Product
#define x_dp4(dest, mask, src0, src1) dest.mask = _ssss(dot(_tof4(src0), _tof4(src1))).mask

// 2.14.1.10.10  DST: Distance Vector
#define x_dst(dest, mask, src0, src1) dest.mask = dst(_tof4(src0), _tof4(src1)).mask /* equals { dest.x = 1; dest.y = src0.y * src1.y; dest.z = src0.z; dest.w = src1.w; } */

// 2.14.1.10.11  MIN: Minimum
#define x_min(dest, mask, src0, src1) dest.mask = min(_tof4(src0), _tof4(src1)).mask

// 2.14.1.10.12  MAX: Maximum
#define x_max(dest, mask, src0, src1) dest.mask = max(_tof4(src0), _tof4(src1)).mask

// 2.14.1.10.13  SLT: Set On Less Than
#define x_slt(dest, mask, src0, src1) dest.mask = _slt(_tof4(src0), _tof4(src1)).mask
float4 _slt(float4 src0, float4 src1)
{
	float4 dest;
	dest.x = (src0.x < src1.x) ? 1 : 0;
	dest.y = (src0.y < src1.y) ? 1 : 0;
	dest.z = (src0.z < src1.z) ? 1 : 0;
	dest.w = (src0.w < src1.w) ? 1 : 0;
	return dest;
}

// 2.14.1.10.14  SGE: Set On Greater or Equal Than
#define x_sge(dest, mask, src0, src1) dest.mask = _sge(_tof4(src0), _tof4(src1)).mask
float4 _sge(float4 src0, float4 src1)
{
	float4 dest;
	dest.x = (src0.x >= src1.x) ? 1 : 0;
	dest.y = (src0.y >= src1.y) ? 1 : 0;
	dest.z = (src0.z >= src1.z) ? 1 : 0;
	dest.w = (src0.w >= src1.w) ? 1 : 0;
	return dest;
}

// 2.14.1.10.18  DPH: Homogeneous Dot Product
#define x_dph(dest, mask, src0, src1) dest.mask = _ssss(_dph(_tof4(src0), _tof4(src1))).mask
float _dph(float4 src0, float4 src1)
{
	return dot(src0.xyz, src1.xyz) + src1.w;
}

// Xbox ILU Functions

// 2.14.1.10.7  RSQ: Reciprocal Square Root
#define x_rsq(dest, mask, src0) dest.mask = _ssss(_rsq(_scalar(src0))).mask
float _rsq(float src)
{
	float a = abs(src);
#if 0 // TODO : Enable
	if (a == 1) return 1;
	if (a == 0) return 1.#INF;
#endif
	return rsqrt(a);
}

// 2.14.1.10.15  EXP: Exponential Base 2
#define x_expp(dest, mask, src0) dest.mask = _expp(_scalar(src0)).mask
float4 _expp(float src)
{
    float floor_src = x_floor(src);

    float4 dest;
    dest.x = exp2(floor_src);
    dest.y = src - floor_src;
    dest.z = exp2(src);
    dest.w = 1;

	return dest;
}

// 2.14.1.10.16  LOG: Logarithm Base 2
#define x_logp(dest, mask, src0) dest.mask = _logp(_scalar(src0)).mask
float4 _logp(float src)
{
    float4 dest;
#if 0 // TODO : Enable
	float t = abs(src);
	if (t != 0) {
		if (t == 1.#INF) {
			dest.x = 1.#INF;
			dest.y = 1;
			dest.z = 1.#INF;
		} else {
#endif
			float exponent;
			float mantissa = frexp(src/* + BIAS*/, /*out*/exponent);
			float z = log2(src);
			dest.x = exponent;
			dest.y = mantissa;
			dest.z = z;
#if 0
		}
	} else {
		dest.x = -1.#INF;
		dest.y = 1;
		dest.z = -1.#INF;
	}
#endif
    dest.w = 1;    
	return dest;
}

// 2.14.1.10.17  LIT: Light Coefficients
#define x_lit(dest, mask, src) dest.mask = _lit(_tof4(src)).mask
float4 _lit(float4 src0)
{
	const float epsilon = 1.0f / 256.0f;

	float diffuse = src0.x;
	float blinn = src0.y;
	float specPower = clamp(src0.w, -(128 - epsilon), (128 - epsilon));

	float4 dest;
	dest.x = 1;
	dest.y = max(0, diffuse);
	dest.z = (diffuse > 0) && (blinn > 0) ? pow(blinn, specPower) : 0;
	dest.w = 1;

	return dest;
}

// 2.14.1.10.19  RCC: Reciprocal Clamped
#define x_rcc(dest, mask, src0) dest.mask = _ssss(_rcc(_scalar(src0))).mask
float _rcc(float src)
{
	// Calculate the reciprocal
	float r = 1 / src;

	// Clamp
	return (r >= 0)
		? clamp(r,  5.42101e-020f,  1.84467e+019f)  // the IEEE 32-bit binary values 0x1F800000 and 0x5F800000
		: clamp(r, -1.84467e+019f, -5.42101e-020f); // the IEEE 32-bit binary values 0xDF800000 and 0x9F800000
}

// 2.14.1.10.6  RCP: Reciprocal
#define x_rcp(dest, mask, src0) dest.mask = _ssss(_rcp(_scalar(src0))).mask
float _rcp(float src)
{
	// OpenGL/NVidia extension definition
#if 0 // TODO : Enable?
	if (src == 1) return 1;
	if (src == 0) return 1.#INF;
	return 1 / src;
#endif
	// Forward to Xbox clamped reciprocal
	// So we have defined behaviour with rcp(0)
	// This prevents issues with XYZRHW modes
	// where the w component may be 0
	return _rcc(src);
}

VS_OUTPUT main(const VS_INPUT xIn)
{
	// Output variables
	float4 oPos, oD0, oD1, oB0, oB1, oT0, oT1, oT2, oT3;
	oPos = oD0 = oD1 = oB0 = oB1 = oT0 = oT1 = oT2 = oT3 = float4(0, 0, 0, 1); // Pre-initialize w component of outputs to 1

	// Single component outputs
	float4 oFog, oPts; // x is write-only on Xbox. Use float4 as some games use incorrect masks
	oFog = 1; // Default to no fog. Test case: Lego Star Wars II
	oPts = 0;

	// Address (index) register
	int1 a0 = 0;

	// Temporary registers
	float4 r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = r11 = float4(0, 0, 0, 0);
	#define r12 oPos // oPos and r12 are two ways of accessing the same register on Xbox

	// Input registers
	float4 v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

#ifdef CXBX_IA_BYPASS
	// IA bypass: fetch all attributes from ByteAddressBuffer using SV_VertexID
	{
		uint xboxVtxIdx = ResolveVertexIndex(xIn.vertexId);
		float4 vArr[16];
		FetchAllAttributes(xboxVtxIdx, vArr);
		v0=vArr[0]; v1=vArr[1]; v2=vArr[2]; v3=vArr[3];
		v4=vArr[4]; v5=vArr[5]; v6=vArr[6]; v7=vArr[7];
		v8=vArr[8]; v9=vArr[9]; v10=vArr[10]; v11=vArr[11];
		v12=vArr[12]; v13=vArr[13]; v14=vArr[14]; v15=vArr[15];
	}
#elif defined(CXBX_USE_D3D11)
	// D3D11: The input assembler delivers correct values for all 16 attributes.
	// Streamed attributes come from their real vertex buffer slots.
	// Non-streamed attributes come from the zero-stride defaults buffer (slot 16)
	// which holds the NV2A's sticky inline_value[] registers.
	#define init_v(i) v##i = xIn.v[i];
	// Note : unroll manually instead of for-loop, because of the ## concatenation
	init_v( 0); init_v( 1); init_v( 2); init_v( 3);
	init_v( 4); init_v( 5); init_v( 6); init_v( 7);
	init_v( 8); init_v( 9); init_v(10); init_v(11);
	init_v(12); init_v(13); init_v(14); init_v(15);
#else
	// D3D9: Lerp between vertex data and constant buffer defaults
	float vRegisterDefaultFlags[16] = (float[16])vRegisterDefaultFlagsPacked;
	#define init_v(i) v##i = lerp(xIn.v[i], vRegisterDefaultValues[i], vRegisterDefaultFlags[i]);
	// Note : unroll manually instead of for-loop, because of the ## concatenation
	init_v( 0); init_v( 1); init_v( 2); init_v( 3);
	init_v( 4); init_v( 5); init_v( 6); init_v( 7);
	init_v( 8); init_v( 9); init_v(10); init_v(11);
	init_v(12); init_v(13); init_v(14); init_v(15);
#endif

	// Temp variable for paired VS instruction
	float4 temp;

	// Xbox shader program will be inserted here
	// <XBOX SHADER PROGRAM GOES HERE>
	// End Xbox shader program

	// Copy variables to output struct
    VS_OUTPUT xOut;
    
	xOut.oPos = reverseScreenspaceTransform(oPos);
	xOut.oD0 = saturate(oD0);
	xOut.oD1 = saturate(oD1);
	xOut.oFog = oFog.x; // Don't clamp here: table fog modes pass a distance (can be > 1); PS clamps the final fog factor
	xOut.oPts = oPts.x;
	xOut.oB0 = saturate(oB0);
	xOut.oB1 = saturate(oB1);
	// Apply TEXCOORDINDEX remapping: NV2A texture units route interpolated
	// texcoords to stages based on D3DTSS_TEXCOORDINDEX. In D3D11 we do this
	// in the VS since there's no hardware texcoord routing post-interpolation.
	{
		float4 texcoordSets[4] = { oT0, oT1, oT2, oT3 };
		oT0 = texcoordSets[(int)xboxTexCoordIndex.x];
		oT1 = texcoordSets[(int)xboxTexCoordIndex.y];
		oT2 = texcoordSets[(int)xboxTexCoordIndex.z];
		oT3 = texcoordSets[(int)xboxTexCoordIndex.w];
	}
	// Scale textures (TODO : or should we apply this to the input register values?)
	xOut.oT0 = oT0 / xboxTextureScale[0];
	xOut.oT1 = oT1 / xboxTextureScale[1];
	xOut.oT2 = oT2 / xboxTextureScale[2];
	xOut.oT3 = oT3 / xboxTextureScale[3];

	return xOut;
}
