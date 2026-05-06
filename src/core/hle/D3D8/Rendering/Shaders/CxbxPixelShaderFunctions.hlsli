// CxbxPixelShaderFunctions.hlsli — pure math helper functions for pixel shaders
//
// Shared by CxbxRegisterCombinerInterpreter.hlsl (RC interpreter ubershader).
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
// alphaTest.x = AlphaTestEnable, .y = AlphaRef, .z = AlphaFunc (NV2A PGRAPH encoding)
void PerformAlphaTest(const float3 alphaTest, float alpha)
{
	[branch] if (alphaTest.x != 0.0f) {
		uint alphaVal  = (uint)(saturate(alpha)       * 255.0f + 0.5f);
		uint alphaRefI = (uint)(saturate(alphaTest.y) * 255.0f + 0.5f);
		int  alphaFunc = (int)alphaTest.z;
		// NV2A PGRAPH encoding: 0=NEVER,1=LESS,2=EQUAL,3=LEQUAL,4=GREATER,5=NOTEQUAL,6=GEQUAL,7=ALWAYS
		bool alphaPass;
		switch (alphaFunc) {
			case 0:  alphaPass = false;                    break; // NEVER
			case 1:  alphaPass = (alphaVal <  alphaRefI);  break; // LESS
			case 2:  alphaPass = (alphaVal == alphaRefI);  break; // EQUAL
			case 3:  alphaPass = (alphaVal <= alphaRefI);  break; // LEQUAL
			case 4:  alphaPass = (alphaVal >  alphaRefI);  break; // GREATER
			case 5:  alphaPass = (alphaVal != alphaRefI);  break; // NOTEQUAL
			case 6:  alphaPass = (alphaVal >= alphaRefI);  break; // GEQUAL
			default: alphaPass = true;                     break; // 7 = ALWAYS
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

// Dot mapping conversion (PS_DOTMAPPING modes 0-7).
// Remaps a texture register value for use in a dot product calculation.
// Pure function — no cbuffer dependency; both the compiled PS template and
// the RC interpreter ubershader call this with a compile-time or runtime mode.
// Math matches xemu's HW-verified implementation (psh.c sign1/sign2/sign3).
//
// mode: 0=ZERO_TO_ONE, 1=MINUS1_TO_1_D3D, 2=MINUS1_TO_1_GL, 3=MINUS1_TO_1,
//       4=HILO_1, 5=HILO_HEMISPHERE_D3D, 6=HILO_HEMISPHERE_GL, 7=HILO_HEMISPHERE
float3 ApplyDotMapping(uint mode, float4 src)
{
    if (mode == 0u)
        return src.rgb;

    float3 b = round(saturate(src.rgb) * 255.0f);

    if (mode == 1u) // D3D: (byte - 128) / 127
        return (b - 128.0f) / 127.0f;

    if (mode == 2u) { // GL: two's complement with +0.5 bias
        float3 s = (b >= 128.0f) ? (b - 255.5f) : (b + 0.5f);
        return s / 127.5f;
    }

    if (mode == 3u) { // Generic two's complement
        float3 s = (b >= 128.0f) ? (b - 256.0f) : b;
        return s / 127.0f;
    }

    // HILO modes: reconstruct two 16-bit values from ARGB channels.
    // Channel order follows xemu HW verification: HI = (A<<8|R), LO = (G<<8|B)
    {
        float4 c = round(saturate(src) * 255.0f);
        float H = c.a * 256.0f + c.r;  // 0..65535
        float L = c.g * 256.0f + c.b;  // 0..65535

        if (mode == 4u) // HILO_1: unsigned [0,1], Z=1
            return float3(H / 65535.0f, L / 65535.0f, 1.0f);

        // HILO_HEMISPHERE modes 5-7: signed H,L with Z = sqrt(1 - H² - L²)
        float Hs, Ls;
        if (mode == 5u) { // D3D: (val - 32768) / 32767
            Hs = (H - 32768.0f) / 32767.0f;
            Ls = (L - 32768.0f) / 32767.0f;
        } else if (mode == 6u) { // GL: two's complement with +0.5 bias
            Hs = (H >= 32768.0f) ? (H - 65535.5f) / 32767.5f : (H + 0.5f) / 32767.5f;
            Ls = (L >= 32768.0f) ? (L - 65535.5f) / 32767.5f : (L + 0.5f) / 32767.5f;
        } else { // mode 7: Generic two's complement
            Hs = (H >= 32768.0f) ? (H - 65536.0f) / 32767.0f : H / 32767.0f;
            Ls = (L >= 32768.0f) ? (L - 65536.0f) / 32767.0f : L / 32767.0f;
        }
        return float3(Hs, Ls, sqrt(max(0.0f, 1.0f - Hs*Hs - Ls*Ls)));
    }
}

// Shadow depth comparison — pure function, no register reads.
// shadowFunc: 3-bit NV2A SHADOW_ZFUNC (0=NEVER,1=LESS,...,7=ALWAYS)
// sampled: value sampled from the shadow/depth texture
// ref: reference depth (typically projected Z coordinate)
float4 ApplyShadowCompare(uint shadowFunc, float4 sampled, float ref)
{
    if (shadowFunc == 0u) return sampled;       // NEVER — pass through
    if (shadowFunc == 7u) return 1.0f.xxxx;     // ALWAYS

    float depth = sampled.r;
    // Bitmask: bit0=LESS, bit1=EQUAL, bit2=GREATER (matches HW encoding)
    uint cmp = (depth <  ref ? 1u : 0u)
             | (depth == ref ? 2u : 0u)
             | (depth >  ref ? 4u : 0u);
    return ((cmp & shadowFunc) != 0u ? 1.0f : 0.0f).xxxx;
}

// Clip plane / compare mode — pure function, no register reads.
// clipBits: 4-bit per-stage mask from NV_PGRAPH_SHADERCLIPMODE
//   Bit set = discard if coord >= 0 (keep negative half-space)
//   Bit clear = discard if coord < 0 (keep positive half-space)
void ApplyCompareMode(uint clipBits, float4 coords)
{
    bool k0 = (clipBits & 1u) != 0u ? (coords.x >= 0.0f) : (coords.x < 0.0f);
    bool k1 = (clipBits & 2u) != 0u ? (coords.y >= 0.0f) : (coords.y < 0.0f);
    bool k2 = (clipBits & 4u) != 0u ? (coords.z >= 0.0f) : (coords.z < 0.0f);
    bool k3 = (clipBits & 8u) != 0u ? (coords.w >= 0.0f) : (coords.w < 0.0f);
    if (k0 || k1 || k2 || k3)
        discard;
}

#endif // CXBX_PIXEL_SHADER_FUNCTIONS_HLSLI
