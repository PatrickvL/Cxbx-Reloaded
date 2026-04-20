// C++ / HLSL shared state block for fixed function support
#ifdef  __cplusplus
#pragma once

#include <DirectXMath.h> // for XMFLOAT2, XMFLOAT3, XMFLOAT4, XMMATRIX
using namespace DirectX;
#define float2 XMFLOAT2
#define float3 xbox::X_D3DVECTOR // XMFLOAT3
#define float4 XMFLOAT4
#define float4x4 XMMATRIX

#include <array> // for std::array<>
#define arr(name, type, length) std::array<type, length> name

// A float that occupies a full 16-byte constant register (matching Xbox register layout)
#define PADDED_FLOAT(name) alignas(16) float name
// A float3 that occupies a full 16-byte constant register (matching Xbox register layout)
#define PADDED_FLOAT3(name) alignas(16) float3 name
// An integer steering value that occupies a full 16-byte constant register
#define PADDED_INT(name) alignas(16) int32_t name
#define CXBX_STEERING_INT int32_t

#else
// HLSL
#define arr(name, type, length) type name[length]
#define alignas(x)
#define const static

// Ensure each field occupies a full float4 constant register in SM5 cbuffer packing,
// matching the C++ alignas(16) layout for raw Xbox register data upload.
#define PADDED_FLOAT(name) float name; float3 _pad_##name
#define PADDED_FLOAT3(name) float3 name; float _pad_##name
#define PADDED_INT(name) int name; int3 _pad_##name
#define CXBX_STEERING_INT int

#endif //  __cplusplus

namespace FixedFunctionVertexShader {
	// Fog depth is taken from the vertex, rather than generated
	const CXBX_STEERING_INT FOG_DEPTH_NONE = 0;
	// Fog depth is the output Z coordinate
    const CXBX_STEERING_INT FOG_DEPTH_Z = 1;
	// Fog depth is based on the output W coordinate (1 / W)
    const CXBX_STEERING_INT FOG_DEPTH_W = 2;
	// Fog depth is based distance of the vertex from the eye position
    const CXBX_STEERING_INT FOG_DEPTH_RANGE = 3;

    // https://docs.microsoft.com/en-us/windows/win32/direct3d9/fog-formulas
    const CXBX_STEERING_INT FOG_TABLE_NONE = 0;
    const CXBX_STEERING_INT FOG_TABLE_EXP = 1;
    const CXBX_STEERING_INT FOG_TABLE_EXP2 = 2;
    const CXBX_STEERING_INT FOG_TABLE_LINEAR = 3;
}

// Shared HLSL structures
// Taking care with packing rules
// In VS_3_0, packing works in mysterious ways
// * Structs inside arrays are not packed
// * Floats can't be packed at all (?)
// We don't get documented packing until vs_4_0

struct Transforms {
    float4x4 View; // 0
    float4x4 Projection; // 1
    arr(Texture, float4x4, 4); // 2, 3, 4, 5
	// World matrices are 6, 7, 8, 9
	// But we use combined WorldView matrices in the shader
    arr(WorldView, float4x4, 4); 
	arr(WorldViewInverseTranspose, float4x4, 4);
};

// See D3DLIGHT
struct Light {
    // TODO in vs_4_0+ when floats are packable
    // Change colour values to float3
    // And put something useful in the alpha slot
    float4 Diffuse;
    float4 Specular;

    // Viewspace light position
    PADDED_FLOAT3(PositionV);
    PADDED_FLOAT(Range);

    // Viewspace light direction (normalized)
    PADDED_FLOAT3(DirectionVN);
    PADDED_INT(Type); // 1=Point, 2=Spot, 3=Directional

    PADDED_FLOAT3(Attenuation);
    PADDED_FLOAT(Falloff);

    PADDED_FLOAT(CosHalfPhi);
    // cos(theta/2) - cos(phi/2)
    PADDED_FLOAT(SpotIntensityDivisor);
};

struct Material {
    float4 Diffuse;
    float4 Ambient;
    float4 Specular;
    float4 Emissive;

    PADDED_FLOAT(Power);
};

struct Modes {
    PADDED_INT(AmbientMaterialSource);
    PADDED_INT(DiffuseMaterialSource);
    PADDED_INT(SpecularMaterialSource);
    PADDED_INT(EmissiveMaterialSource);

    PADDED_INT(BackAmbientMaterialSource);
    PADDED_INT(BackDiffuseMaterialSource);
    PADDED_INT(BackSpecularMaterialSource);
    PADDED_INT(BackEmissiveMaterialSource);

    PADDED_INT(Lighting);
    PADDED_INT(TwoSidedLighting);
//  PADDED_INT(SpecularEnable);
    PADDED_INT(LocalViewer);

/// PADDED_INT(ColorVertex);
    PADDED_INT(VertexBlend_NrOfMatrices);
    PADDED_INT(VertexBlend_CalcLastWeight); // Could be a bool in higer shader models
    PADDED_INT(NormalizeNormals);
};

struct PointSprite {
    PADDED_FLOAT(PointSize);
    PADDED_FLOAT(PointSize_Min);
    PADDED_FLOAT(PointSize_Max);
//  PADDED_FLOAT(PointScaleEnable);
    PADDED_FLOAT(XboxRenderTargetHeight);
    PADDED_FLOAT3(PointScaleABC);
    PADDED_FLOAT(RenderUpscaleFactor);
};

struct TextureState {
    PADDED_INT(TextureTransformFlagsCount);
    PADDED_INT(TextureTransformFlagsProjected);
    PADDED_INT(TexCoordIndex);
    PADDED_INT(TexCoordIndexGen);
};

struct Fog {
    PADDED_INT(Enable);
    PADDED_INT(DepthMode);
    PADDED_INT(TableMode);
    PADDED_FLOAT(Density); // EXP fog density
    PADDED_FLOAT(Start); // LINEAR fog start
    PADDED_FLOAT(End); // LINEAR fog end
};

// Vertex lighting
// Both frontface and backface lighting can be calculated
struct TwoSidedColor
{
	PADDED_FLOAT3(Front);
	PADDED_FLOAT3(Back);
};

struct FixedFunctionVertexShaderState {
    alignas(16) Transforms Transforms;
    alignas(16) arr(Lights, Light, 8);
    alignas(16) TwoSidedColor TotalLightsAmbient;
    alignas(16) arr(Materials, Material, 2);
    alignas(16) Modes Modes;
    alignas(16) Fog Fog;
    alignas(16) arr(TextureStates, TextureState, 4);
    alignas(16) PointSprite PointSprite;
    alignas(16) float4 TexCoordComponentCount;
};

#ifdef  __cplusplus
#undef float2
#undef float3
#undef float4
#undef float4x4
#else // HLSL
#undef alignas
#undef const
#endif //  __cplusplus

#undef arr
#undef PADDED_FLOAT
#undef PADDED_FLOAT3
#undef PADDED_INT
#undef CXBX_STEERING_INT

#undef arr
#undef PADDED_FLOAT
#undef PADDED_FLOAT3
