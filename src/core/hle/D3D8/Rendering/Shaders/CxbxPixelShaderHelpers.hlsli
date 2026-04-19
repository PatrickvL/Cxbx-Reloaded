#ifndef CXBX_PIXEL_SHADER_HELPERS_HLSLI
#define CXBX_PIXEL_SHADER_HELPERS_HLSLI

static const float4 WarningColor = float4(0, 1, 1, 1); // Returned when unhandled scenario is encountered

#define unsigned_to_signed(x) (((x) * 2) - 1) // Shifts range from [0..1] to [-1..1] (just like s_bx2)
#define signed_to_unsigned(x) (((x) + 1) / 2) // Shifts range from [-1..1] to [0..1]

float4 PerformColorSign(const float4 ColorSign, float4 t)
{
	// Per color channel, based on the ColorSign setting :
	// either keep the value range as-is (when ColorSign is zero)
	// or convert from [0..1] to [-1..+1] (when ColorSign is more than zero, often used for bumpmaps),
	// or convert from [-1..1] to [0..1] (when ColorSign is less than zero):
	if (ColorSign.r > 0) t.r = unsigned_to_signed(t.r);
	if (ColorSign.g > 0) t.g = unsigned_to_signed(t.g);
	if (ColorSign.b > 0) t.b = unsigned_to_signed(t.b);
	if (ColorSign.a > 0) t.a = unsigned_to_signed(t.a);
	if (ColorSign.r < 0) t.r = signed_to_unsigned(t.r);
	if (ColorSign.g < 0) t.g = signed_to_unsigned(t.g);
	if (ColorSign.b < 0) t.b = signed_to_unsigned(t.b);
	if (ColorSign.a < 0) t.a = signed_to_unsigned(t.a);
	// TODO : Instead of the above, create a mirror texture with a host format that has identical component layout, but with all components signed.
	// Then, in here, when any component has to be read as signed, sample the signed texture (ouch : with what dimension and coordinate?!)
	// and replace the components that we read from the unsigned texture, but which have to be signed, with the signed components read from the signed mirror texture.
	// This way, texture filtering can still be allowed, as that would be performed separately over the unsigned vs unsigned textures (so no mixing between the two).

	return t;
}

float4 PerformColorKeyOp(const float ColorKeyOp, const float4 ColorKeyColor, float4 t)
{
	// Handle all D3DTCOLORKEYOP_ modes :
	if (ColorKeyOp == 0) // = _DISABLE
		return t; // No color-key checking

	if (any(t - ColorKeyColor))
		return t; // Cxbx assumption : On color mismatch, simply return the input. TODO : This might require a more elaborate operation? (Like "when any of the texels were filtered with a non-zero weight", whatever that means)

	if (ColorKeyOp == 1) // = _ALPHA
		return float4(t.rgb, 0);

	if (ColorKeyOp == 2) // = _RGBA
		return 0;

	if (ColorKeyOp == 3) // = _KILL
		clip(-1);

	// Undefined ColorKeyOp mode
	return WarningColor;
}

void PerformAlphaKill(const float AlphaKill, float4 t)
{
	if (AlphaKill)
		if (t.a == 0)
			clip(-1);
}

// D3D11: Alpha test (D3D11 has no fixed-function alpha test)
// NV2A quantizes both alpha output and reference to 8-bit before comparing,
// so we replicate that to avoid float precision issues with == and !=.
// alphaTest.x = AlphaTestEnable, .y = AlphaRef, .z = AlphaFunc (D3DCMPFUNC)
void PerformAlphaTest(const float3 alphaTest, float alpha)
{
	if (alphaTest.x) {
		int alphaVal = (int)round(saturate(alpha) * 255);
		int alphaRefI = (int)round(saturate(alphaTest.y) * 255);
		int alphaFunc = (int)alphaTest.z;
		// D3DCMPFUNC: 1=NEVER,2=LESS,3=EQUAL,4=LESSEQUAL,5=GREATER,6=NOTEQUAL,7=GREATEREQUAL,8=ALWAYS
		bool alphaPass = (alphaFunc == 8); // ALWAYS
		if (alphaFunc == 1) alphaPass = false;                          // NEVER
		if (alphaFunc == 2) alphaPass = (alphaVal < alphaRefI);         // LESS
		if (alphaFunc == 3) alphaPass = (alphaVal == alphaRefI);        // EQUAL
		if (alphaFunc == 4) alphaPass = (alphaVal <= alphaRefI);        // LESSEQUAL
		if (alphaFunc == 5) alphaPass = (alphaVal > alphaRefI);         // GREATER
		if (alphaFunc == 6) alphaPass = (alphaVal != alphaRefI);        // NOTEQUAL
		if (alphaFunc == 7) alphaPass = (alphaVal >= alphaRefI);        // GREATEREQUAL
		if (!alphaPass) clip(-1);
	}
}

#endif // CXBX_PIXEL_SHADER_HELPERS_HLSLI
