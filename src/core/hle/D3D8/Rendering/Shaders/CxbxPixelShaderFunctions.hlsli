// CxbxPixelShaderFunctions.hlsli — pure math helper functions for pixel shaders
//
// Shared by CxbxPixelShaderTemplate.hlsl (compiled PS), CxbxFixedFunctionPixelShader.hlsl,
// and CxbxRegisterCombinerInterpreter.hlsl (RC interpreter ubershader).
//
// This file must NOT declare any I/O structs, samplers, or texture objects
// so that it can be included from shaders with different I/O declarations.

#ifndef CXBX_PIXEL_SHADER_FUNCTIONS_HLSLI
#define CXBX_PIXEL_SHADER_FUNCTIONS_HLSLI

// Integer steering type: D3D11 uses native int, D3D9 used float for SM3 compat
#ifndef CXBX_STEERING_INT
#define CXBX_STEERING_INT int
#endif

// Color sign conversion (Xbox X_D3DTSS_COLORSIGN extension)
// Per channel: >0 = expand [0,1]→[-1,1]; <0 = contract [-1,1]→[0,1]; 0 = identity
float4 PerformColorSign(const float4 ColorSign, float4 t)
{
	// Vectorized: compiler emits movc per channel, no per-component branches
	float4 expand   = t * 2.0f - 1.0f;   // unsigned_to_signed
	float4 contract = t * 0.5f + 0.5f;   // signed_to_unsigned
	bool4  pos = ColorSign > 0.0f;
	bool4  neg = ColorSign < 0.0f;
	return pos ? expand : (neg ? contract : t);
}

// Color key operation (D3DTCOLORKEYOP)
// Compare at 8-bit precision: bilinear-filtered samples rarely hit exact
// float equality, but Xbox hardware compares at native texel precision (8-bit).
float4 PerformColorKeyOp(const CXBX_STEERING_INT ColorKeyOp, const float4 ColorKeyColor, float4 t)
{
	if (ColorKeyOp == 0) // _DISABLE
		return t;

	uint4 tI = (uint4)(saturate(t)              * 255.0f + 0.5f);
	uint4 kI = (uint4)(saturate(ColorKeyColor)  * 255.0f + 0.5f);
	if (any(tI != kI))
		return t; // No match

	if (ColorKeyOp == 1) // _ALPHA
		return float4(t.rgb, 0);

	if (ColorKeyOp == 2) // _RGBA
		return (float4)0;

	if (ColorKeyOp == 3) // _KILL
		clip(-1);

	return t;
}

// Alpha kill (D3DTALPHAKILL_ENABLE)
void PerformAlphaKill(const CXBX_STEERING_INT AlphaKill, float4 t)
{
	if (AlphaKill)
		if (t.a == 0)
			clip(-1);
}

// Alpha test (D3D11 has no fixed-function alpha test)
// NV2A quantizes both alpha output and reference to 8-bit before comparing.
// alphaTest.x = AlphaTestEnable, .y = AlphaRef, .z = AlphaFunc (D3DCMPFUNC)
void PerformAlphaTest(const float3 alphaTest, float alpha)
{
	[branch] if (alphaTest.x != 0.0f) {
		uint alphaVal  = (uint)(saturate(alpha)       * 255.0f + 0.5f);
		uint alphaRefI = (uint)(saturate(alphaTest.y) * 255.0f + 0.5f);
		int  alphaFunc = (int)alphaTest.z;
		// D3DCMPFUNC: 1=NEVER,2=LESS,3=EQUAL,4=LESSEQUAL,5=GREATER,6=NOTEQUAL,7=GREATEREQUAL,8=ALWAYS
		bool alphaPass;
		switch (alphaFunc) {
			case 1:  alphaPass = false;                    break; // NEVER
			case 2:  alphaPass = (alphaVal <  alphaRefI);  break; // LESS
			case 3:  alphaPass = (alphaVal == alphaRefI);  break; // EQUAL
			case 4:  alphaPass = (alphaVal <= alphaRefI);  break; // LESSEQUAL
			case 5:  alphaPass = (alphaVal >  alphaRefI);  break; // GREATER
			case 6:  alphaPass = (alphaVal != alphaRefI);  break; // NOTEQUAL
			case 7:  alphaPass = (alphaVal >= alphaRefI);  break; // GREATEREQUAL
			default: alphaPass = true;                     break; // 8 = ALWAYS
		}
		if (!alphaPass) clip(-1);
	}
}

// Texture format channel fixup (D3D11: luminance replication, channel swizzle)
// fixup: 0=identity, 1=.gbar, 2=.abgr, 3=luminance, 4=alpha-luminance, 5=opaque-alpha
float4 ApplyTexFmtFixup(float4 t, CXBX_STEERING_INT fixup)
{
	[branch] if (fixup != 0) {
		if (fixup == 1) return t.gbar;                       // B8G8R8A8 uploaded as R8G8B8A8
		if (fixup == 2) return t.abgr;                       // A8B8G8R8 uploaded as R8G8B8A8
		if (fixup == 3) return float4(t.r, t.r, t.r, t.a);   // Luminance: R→(R,R,R,A)
		if (fixup == 4) return float4(t.r, t.r, t.r, t.g);   // Alpha-luminance: RG→(R,R,R,G)
		if (fixup == 5) return float4(t.rgb, 1.0f);           // Opaque-alpha: X8R8G8B8/X1R5G5B5
	}
	return t;
}

// Fog factor computation (caller must check fogEnable before calling)
// fogTableMode: 0=NONE (vertex fog passthrough), 1=EXP, 2=EXP2, 3=LINEAR
float CalculateFogFactor(CXBX_STEERING_INT fogTableMode, float fogDensity, float fogStart, float fogEnd, float fogDepth)
{
	if (fogTableMode == 1)      // EXP
		return 1.0f / exp(fogDepth * fogDensity);
	else if (fogTableMode == 2) // EXP2
		return 1.0f / exp(pow(fogDepth * fogDensity, 2));
	else if (fogTableMode == 3) // LINEAR — per spec, clamped to [0,1]
		return saturate((fogEnd - fogDepth) / (fogEnd - fogStart));
	else                        // 0 = NONE (vertex fog passthrough, already [0,1])
		return fogDepth;
}

#endif // CXBX_PIXEL_SHADER_FUNCTIONS_HLSLI
