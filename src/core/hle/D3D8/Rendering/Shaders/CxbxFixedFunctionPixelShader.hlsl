#include "CxbxFixedFunctionPixelShader.hlsli"

#include "CxbxPixelShaderInput.hlsli"

// Texture objects and samplers — DX11 SM5 style.
// Register layout matches all other pixel shaders:
//   t0-t3  = Texture2D,  t4-t7  = Texture3D,  t8-t11 = TextureCube
//   s0-s3  = shared SamplerState
Texture2D   Tex2D_0   : register(t0);
Texture2D   Tex2D_1   : register(t1);
Texture2D   Tex2D_2   : register(t2);
Texture2D   Tex2D_3   : register(t3);
Texture3D   Tex3D_0   : register(t4);
Texture3D   Tex3D_1   : register(t5);
Texture3D   Tex3D_2   : register(t6);
Texture3D   Tex3D_3   : register(t7);
TextureCube TexCube_0 : register(t8);
TextureCube TexCube_1 : register(t9);
TextureCube TexCube_2 : register(t10);
TextureCube TexCube_3 : register(t11);

SamplerState Samp0   : register(s0);
SamplerState Samp1   : register(s1);
SamplerState Samp2   : register(s2);
SamplerState Samp3   : register(s3);

static const float4 WarningColor = float4(0, 1, 1, 1); // Returned when unhandled scenario is encountered

#include "CxbxPixelShaderFunctions.hlsli"

uniform FixedFunctionPixelShaderState state : register(c0);

// These 'D3DTA' texture argument values
// may be used during each texture stage
struct TextureArgs {
	float4 CURRENT;
	float4 TEXTURE;
	float4 DIFFUSE;
	float4 SPECULAR;
	float4 TEMP;
	float4 TFACTOR;
};

static float4 TexCoords[4];

// When creating an instance of the fixed function shader
// we string-replace the assignment below with a value
// The define keeps the shader compilable without the replacement
#define TEXTURE_SAMPLE_TYPE {SAMPLE_2D, SAMPLE_2D, SAMPLE_2D, SAMPLE_2D};
static int TextureSampleType[4] = TEXTURE_SAMPLE_TYPE;

// D3D11: Use native integer bitwise operations for D3DTA flag extraction
bool HasFlag(int value, int flag) {
	return (value & flag) != 0;
}

float4 GetArg(int arg, TextureArgs ctx) {
	// https://docs.microsoft.com/en-us/windows/win32/direct3d9/d3dta
	bool alphaReplicate = (arg & X_D3DTA_ALPHAREPLICATE) != 0;
	bool complement = (arg & X_D3DTA_COMPLEMENT) != 0;
	arg = arg & 0xF; // Extract base selector (bits 0-3)

	float4 o;

	if (arg == X_D3DTA_DIFFUSE)
		o = ctx.DIFFUSE;
	if (arg == X_D3DTA_CURRENT)
		o = ctx.CURRENT;
	if (arg == X_D3DTA_TEXTURE)
		o = ctx.TEXTURE;
	if (arg == X_D3DTA_TFACTOR)
		o = ctx.TFACTOR;
	if (arg == X_D3DTA_SPECULAR)
		o = ctx.SPECULAR;
	if (arg == X_D3DTA_TEMP)
		o = ctx.TEMP;

	if (alphaReplicate)
		return o.aaaa;
	else if (complement)
		return 1 - o;
	else
		return o;
}

float4 ExecuteTextureOp(int op, float4 arg1, float4 arg2, float4 arg0, TextureArgs ctx, PsTextureStageState stage) {
	// https://docs.microsoft.com/en-us/windows/win32/direct3d9/d3dtextureop

	// Note : When we use separate "if"'s here instead of below "else if"'s,
	// D3DCompile may stackoverflow at runtime
	if (op == X_D3DTOP_SELECTARG1)
		return arg1;
	else if (op == X_D3DTOP_SELECTARG2)
		return arg2;
	else if (op == X_D3DTOP_MODULATE)
		return arg1 * arg2;
	else if (op == X_D3DTOP_MODULATE2X)
		return 2 * (arg1 * arg2);
	else if (op == X_D3DTOP_MODULATE4X)
		return 4 * (arg1 * arg2);
	else if (op == X_D3DTOP_ADD)
		return arg1 + arg2;
	else if (op == X_D3DTOP_ADDSIGNED)
		return arg1 + arg2 - 0.5;
	else if (op == X_D3DTOP_ADDSIGNED2X)
		return 2 * (arg1 + arg2 - 0.5);
	else if (op == X_D3DTOP_SUBTRACT)
		return arg1 - arg2;
	else if (op == X_D3DTOP_ADDSMOOTH)
		return arg1 + arg2 * (1 - arg1);
	else if (op == X_D3DTOP_BLENDDIFFUSEALPHA)
		return arg1 * ctx.DIFFUSE.a + arg2 * (1 - ctx.DIFFUSE.a);
	else if (op == X_D3DTOP_BLENDCURRENTALPHA)
		return arg1 * ctx.CURRENT.a + arg2 * (1 - ctx.CURRENT.a);
	else if (op == X_D3DTOP_BLENDTEXTUREALPHA)
		return arg1 * ctx.TEXTURE.a + arg2 * (1 - ctx.TEXTURE.a);
	else if (op == X_D3DTOP_BLENDFACTORALPHA)
		return arg1 * ctx.TFACTOR.a + arg2 * (1 - ctx.TFACTOR.a);
	else if (op == X_D3DTOP_BLENDTEXTUREALPHAPM)
		return arg1 + arg2 * (1 - ctx.TEXTURE.a);
	else if (op == X_D3DTOP_PREMODULATE)
		return arg1; // Note this also multiplies the next stage's CURRENT by its texture
	else if (op == X_D3DTOP_MODULATEALPHA_ADDCOLOR)
		return float4(arg1.rgb + arg1.a * arg2.rgb, 1);
	else if (op == X_D3DTOP_MODULATECOLOR_ADDALPHA)
		return float4(arg1.rgb * arg2.rgb + arg1.a, 1);
	else if (op == X_D3DTOP_MODULATEINVALPHA_ADDCOLOR)
		return float4((1 - arg1.a) * arg2.rgb + arg1.rgb, 1);
	else if (op == X_D3DTOP_MODULATEINVCOLOR_ADDALPHA)
		return float4((1 - arg1.rgb) * arg2.rgb + arg1.a, 1);
	else if (op == X_D3DTOP_DOTPRODUCT3) {
		// Test case: PerPixelLighting
		// DOT_PRODUCT3 expands inputs from [0,1] to [-1,1] range before the dot product
		arg1.rgb = (arg1.rgb - 0.5) * 2;
		arg2.rgb = (arg2.rgb - 0.5) * 2;
		return saturate(dot(arg1.rgb, arg2.rgb));
	}
	// Note arg0 below is arg1 in D3D docs
	// since it becomes the first argument for operations supporting 3 arguments...
	else if (op == X_D3DTOP_MULTIPLYADD)
		return arg0 + arg1 * arg2;
	else if (op == X_D3DTOP_LERP)
		return arg0 * arg1 + (1 - arg0) * arg2;
	else if (op >= X_D3DTOP_BUMPENVMAP) { // Also handles X_D3DTOP_BUMPENVMAPLUMINANCE
		arg1 = ctx.CURRENT; // Bump mapping uses the previous stage's result (CURRENT) as the perturbation source
		arg1.xy = (arg1.xy - 0.5) * 2; // Expand perturbation from [0,1] to [-1,1]
		// Note : default component order .xyzw is identical to .rgba, and z (red) should not be scaled here, as it's used for luminance which is an unsigned input.
		return float4(
			arg1.x * stage.BUMPENVMAT00 + arg1.y * stage.BUMPENVMAT10,
			arg1.x * stage.BUMPENVMAT01 + arg1.y * stage.BUMPENVMAT11,
			(op == X_D3DTOP_BUMPENVMAPLUMINANCE) ? arg1.z * stage.BUMPENVLSCALE + stage.BUMPENVLOFFSET : 1,
			1);
	}
	// Something is amiss... we should have returned by now!
	return WarningColor;
}

TextureArgs ExecuteTextureStage(
	int i,
	TextureArgs ctx,
	PsTextureHardcodedState s,
	int previousOp
)
{
	// Early exit if this stage is disabled (and therefore all further stages are too)
	if (s.COLOROP == X_D3DTOP_DISABLE)
		return ctx;

	PsTextureStageState stage = state.stages[i];

	// Fetch the texture coordinates for this stage
	float3 TexCoord = TexCoords[i].xyz;

	// Bump environment mapping special case
	if (previousOp >= X_D3DTOP_BUMPENVMAP) { // Also handles X_D3DTOP_BUMPENVMAPLUMINANCE
		// Assume U, V, L is in CURRENT
		// Add U', V', to the texture coordinates
		// https://docs.microsoft.com/en-us/windows/win32/direct3d9/bump-mapping-formulas
		TexCoord.xy += ctx.CURRENT.xy;
	}

	// Sample the texture
	float4 t;
	int type = TextureSampleType[i];
	if (type == SAMPLE_NONE)
		t = 1; // Test case JSRF
	else if (type == SAMPLE_2D) {
		// Use switch to dispatch to individual texture/sampler pairs
		switch (i) {
			case 0: t = Tex2D_0.Sample(Samp0, TexCoord.xy); break;
			case 1: t = Tex2D_1.Sample(Samp1, TexCoord.xy); break;
			case 2: t = Tex2D_2.Sample(Samp2, TexCoord.xy); break;
			case 3: t = Tex2D_3.Sample(Samp3, TexCoord.xy); break;
		}
	}
	else if (type == SAMPLE_3D) {
		switch (i) {
			case 0: t = Tex3D_0.Sample(Samp0, TexCoord.xyz); break;
			case 1: t = Tex3D_1.Sample(Samp1, TexCoord.xyz); break;
			case 2: t = Tex3D_2.Sample(Samp2, TexCoord.xyz); break;
			case 3: t = Tex3D_3.Sample(Samp3, TexCoord.xyz); break;
		}
	}
	else if (type == SAMPLE_CUBE) {
		switch (i) {
			case 0: t = TexCube_0.Sample(Samp0, TexCoord.xyz); break;
			case 1: t = TexCube_1.Sample(Samp1, TexCoord.xyz); break;
			case 2: t = TexCube_2.Sample(Samp2, TexCoord.xyz); break;
			case 3: t = TexCube_3.Sample(Samp3, TexCoord.xyz); break;
		}
	}

	// Apply texture format channel fixup (D3D11: luminance replication, channel swizzle)
	t = ApplyTexFmtFixup(t, stage.TEXFMTFIXUP);

	// Bump environment mapping with luminance special case
	if (previousOp == X_D3DTOP_BUMPENVMAPLUMINANCE) {
		// Multiply sampled texture rgb values by L'
		t.rgb *= ctx.CURRENT.z;
	}

	// TODO : Figure out in which order the following operations should be performed :
	t = PerformColorSign(stage.COLORSIGN, t);
	t = PerformColorKeyOp(stage.COLORKEYOP, stage.COLORKEYCOLOR, t);
	PerformAlphaKill(stage.ALPHAKILL, t);

	// Assign the final value for TEXTURE
	ctx.TEXTURE = t;

	// Premodulate special case
	if (previousOp == X_D3DTOP_PREMODULATE) {
		ctx.CURRENT *= ctx.TEXTURE;
	}

	// Get arguments for the texture operation
	// Almost all operate on 2 arguments, Arg1 and Arg2
	// Arg0 is a third argument that seems to have been tacked on
	// for MULTIPLYADD and LERP

	// Colour operation arguments
	float4 cArg1 = GetArg(s.COLORARG1, ctx);
	float4 cArg2 = GetArg(s.COLORARG2, ctx);
	float4 cArg0 = GetArg(s.COLORARG0, ctx);

	// Alpha operation arguments
	float4 aArg1 = GetArg(s.ALPHAARG1, ctx);
	float4 aArg2 = GetArg(s.ALPHAARG2, ctx);
	float4 aArg0 = GetArg(s.ALPHAARG0, ctx);

	// Execute texture operation
	// ALPHAOP == X_D3DTOP_DISABLE is undefined behaviour
	// Using an intermediate value matches known cases...
	// Test case: DoA:Xtreme (menu water), GTA III (logos), Crash Wrath of Cortex (relics UI)
	static float4 value = 1;
	value.rgb = ExecuteTextureOp(s.COLOROP, cArg1, cArg2, cArg0, ctx, stage).rgb;
	if (s.ALPHAOP != X_D3DTOP_DISABLE) {
		value.a = ExecuteTextureOp(s.ALPHAOP, aArg1, aArg2, aArg0, ctx, stage).a;
	}

	// Save the result
	// Note RESULTARG should either be CURRENT or TEMP
	// But some titles seem to set it to DIFFUSE
	// Use CURRENT for anything other than TEMP
	// Test case: DoA 3
	if (s.RESULTARG == X_D3DTA_TEMP)
		ctx.TEMP = value;
	else
		ctx.CURRENT = value;

	return ctx;
}

float4 main(const PS_INPUT input) : SV_Target {
    // iFog is already a computed fog factor from the VS (EXP/EXP2/LINEAR/passthrough),
    // linearly interpolated by the rasterizer. NV2A clamps to [0,1] before use.
    float fogFactor = state.FogEnable
        ? saturate(input.iFog)
        : 1.0f;

    // Map input texture coordinates to an array, for indexing purposes
    TexCoords[0] = input.iT0;
    TexCoords[1] = input.iT1;
    TexCoords[2] = input.iT2;
    TexCoords[3] = input.iT3;

	// Each stage is passed and returns
	// a set of texture arguments
	// And will usually update the CURRENT value
	TextureArgs ctx;

	// The CURRENT register
	// Default to the diffuse value
	// TODO determine whether to use the front or back colours
	// and set them here
	ctx.CURRENT = input.iD0;
	ctx.DIFFUSE = input.iD0;
	ctx.SPECULAR = input.iD1;
	// The TEMP register
	// Default to 0
	ctx.TEMP = float4(0, 0, 0, 0);
	ctx.TFACTOR = state.TextureFactor;

	PsTextureHardcodedState stages[4];
	stages[0].COLOROP = X_D3DTOP_DISABLE;
	stages[1].COLOROP = X_D3DTOP_DISABLE;
	stages[2].COLOROP = X_D3DTOP_DISABLE;
	stages[3].COLOROP = X_D3DTOP_DISABLE;

	// Define stages
	// https://docs.microsoft.com/en-us/windows/win32/direct3d9/d3dtexturestagestatetype
	// We'll find comment below and insert the definitions after it
	// STAGE DEFINITIONS
	// END STAGE DEFINITIONS

	// Run each stage
	int previousOp = -1;
	for (int i = 0; i < 4; i++) {

		ctx = ExecuteTextureStage(
			i,
			ctx,
			stages[i],
			previousOp
		);

		previousOp = stages[i].COLOROP;
	}

	// Add fog if enabled
	if (state.FogEnable) {
		ctx.CURRENT.rgb = lerp(state.FogColor.rgb, ctx.CURRENT.rgb, fogFactor);
	}

	// Add specular if enabled
	if (state.SpecularEnable) {
		ctx.CURRENT.rgb += ctx.SPECULAR.rgb;
	}

	PerformAlphaTest(state.AlphaTest.xyz, ctx.CURRENT.a);

	// Output whatever is in current at the end
	return ctx.CURRENT;
}
