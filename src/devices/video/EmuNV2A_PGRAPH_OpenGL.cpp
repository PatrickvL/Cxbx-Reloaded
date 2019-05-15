// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
// ******************************************************************
// *
// *  This file is part of the Cxbx project.
// *
// *  Cxbx and Cxbe are free software; you can redistribute them
// *  and/or modify them under the terms of the GNU General Public
// *  License as published by the Free Software Foundation; either
// *  version 2 of the license, or (at your option) any later version.
// *
// *  This program is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// *  GNU General Public License for more details.
// *
// *  You should have recieved a copy of the GNU General Public License
// *  along with this program; see the file COPYING.
// *  If not, write to the Free Software Foundation, Inc.,
// *  59 Temple Place - Suite 330, Bostom, MA 02111-1307, USA.
// *
// *  This file is heavily based on code from XQEMU
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a_pgraph.c
// *  Copyright (c) 2012 espes
// *  Copyright (c) 2015 Jannik Vogel
// *  Copyright (c) 2018 Matt Borgerson
// *
// *  Contributions for Cxbx-Reloaded
// *  Copyright (c) 2017-2018 Luke Usher <luke.usher@outlook.com>
// *  Copyright (c) 2018 Patrick van Logchem <pvanlogchem@gmail.com>
// *
// *  All rights reserved
// *
// ******************************************************************

// Xbox uses 4 KiB pages
#define TARGET_PAGE_BITS 12
#define TARGET_PAGE_SIZE (1 << TARGET_PAGE_BITS)
#define TARGET_PAGE_MASK ~(TARGET_PAGE_SIZE - 1)
#define TARGET_PAGE_ALIGN(addr) (((addr) + TARGET_PAGE_SIZE - 1) & TARGET_PAGE_MASK)

enum FormatEncoding {
	linear = 0,
	swizzled, // for all NV097_SET_TEXTURE_FORMAT_*_SZ_*
	compressed // for all NV097_SET_TEXTURE_FORMAT_*_DXT*
};

typedef struct ColorFormatInfo {
    unsigned int bytes_per_pixel; // Derived from the total number of channel bits
    FormatEncoding encoding;
    GLint gl_internal_format;
    GLenum gl_format; // == 0 for compressed formats
    GLenum gl_type;
    GLint *gl_swizzle_mask; // == nullptr when gl_internal_format, gl_format and gl_type are sufficient
} ColorFormatInfo;

// Resulting gl_internal_format, gl_format and gl_type values, for formats handled by convert_texture_data()
#define GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV

static GLint gl_swizzle_mask_0RG1[4] = { GL_ZERO, GL_RED, GL_GREEN, GL_ONE };
static GLint gl_swizzle_mask_111R[4] = { GL_ONE, GL_ONE, GL_ONE, GL_RED };
static GLint gl_swizzle_mask_ARGB[4] = { GL_ALPHA, GL_RED, GL_GREEN, GL_BLUE };
static GLint gl_swizzle_mask_BGRA[4] = { GL_BLUE, GL_GREEN, GL_RED, GL_ALPHA };
static GLint gl_swizzle_mask_GGGR[4] = { GL_GREEN, GL_GREEN, GL_GREEN, GL_RED };
static GLint gl_swizzle_mask_R0G1[4] = { GL_RED, GL_ZERO, GL_GREEN, GL_ONE };
static GLint gl_swizzle_mask_RRR1[4] = { GL_RED, GL_RED, GL_RED, GL_ONE };
static GLint gl_swizzle_mask_RRRG[4] = { GL_RED, GL_RED, GL_RED, GL_GREEN };
static GLint gl_swizzle_mask_RRRR[4] = { GL_RED, GL_RED, GL_RED, GL_RED };

// Note : Avoid designated initializers to facilitate C++ builds
static const ColorFormatInfo kelvin_color_format_map[256] = {
    //0x00 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_Y8] =
        {1, swizzled, GL_R8, GL_RED, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_RRR1},
    //0x01 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_AY8] =
        {1, swizzled, GL_R8, GL_RED, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_RRRR},
    //0x02 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_A1R5G5B5] =
        {2, swizzled, GL_RGB5_A1, GL_BGRA, GL_UNSIGNED_SHORT_1_5_5_5_REV},
    //0x03 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_X1R5G5B5] =
        {2, swizzled, GL_RGB5, GL_BGRA, GL_UNSIGNED_SHORT_1_5_5_5_REV},
    //0x04 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_A4R4G4B4] =
        {2, swizzled, GL_RGBA4, GL_BGRA, GL_UNSIGNED_SHORT_4_4_4_4_REV},
    //0x05 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R5G6B5] =
        {2, swizzled, GL_RGB565, GL_RGB, GL_UNSIGNED_SHORT_5_6_5},
    //0x06 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_A8R8G8B8] =
        {4, swizzled, GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV},
    //0x07 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_X8R8G8B8] =
        {4, swizzled, GL_RGB8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV},
	//0x08 [?] =
		{},
	//0x09 [?] =
		{},
	//0x0A [?] =
		{},

    /* paletted texture */
    //0x0B [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_I8_A8R8G8B8] = // See convert_texture_data
        {1, swizzled, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT},

    //0x0C [NV097_SET_TEXTURE_FORMAT_COLOR_L_DXT1_A1R5G5B5] =
        {4, compressed, GL_COMPRESSED_RGBA_S3TC_DXT1_EXT, 0, GL_RGBA},
	//0x0D [?] =
		{},
    //0x0E [NV097_SET_TEXTURE_FORMAT_COLOR_L_DXT23_A8R8G8B8] =
        {4, compressed, GL_COMPRESSED_RGBA_S3TC_DXT3_EXT, 0, GL_RGBA},
    //0x0F [NV097_SET_TEXTURE_FORMAT_COLOR_L_DXT45_A8R8G8B8] =
        {4, compressed, GL_COMPRESSED_RGBA_S3TC_DXT5_EXT, 0, GL_RGBA},
    //0x10 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A1R5G5B5] =
        {2, linear, GL_RGB5_A1, GL_BGRA, GL_UNSIGNED_SHORT_1_5_5_5_REV},
    //0x11 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R5G6B5] =
        {2, linear, GL_RGB565, GL_RGB, GL_UNSIGNED_SHORT_5_6_5},
    //0x12 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8R8G8B8] =
        {4, linear, GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV},
    //0x13 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_Y8] =
        {1, linear, GL_R8, GL_RED, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_RRR1},
	//0x14 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_SY8] =
        {1, linear, GL_R8_SNORM, GL_RED, GL_BYTE,
         gl_swizzle_mask_RRR1}, // TODO : Verify
	//0x15 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_X7SY9] = // See convert_texture_data FIXME
		{2, linear, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify
	//0x16 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R8B8] =
        {2, linear, GL_RG8, GL_RG, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_R0G1}, // TODO : Verify
	//0x17 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_G8B8] =
        {2, linear, GL_RG8, GL_RG, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_0RG1}, // TODO : Verify
	//0x18 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_SG8SB8] =
        {2, linear, GL_RG8_SNORM, GL_RG, GL_BYTE,
         gl_swizzle_mask_0RG1}, // TODO : Verify

    //0x19 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_A8] =
        {1, swizzled, GL_R8, GL_RED, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_111R},
    //0x1A [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_A8Y8] =
        {2, swizzled, GL_RG8, GL_RG, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_GGGR},
    //0x1B [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_AY8] =
        {1, linear, GL_R8, GL_RED, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_RRRR},
    //0x1C [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_X1R5G5B5] =
        {2, linear, GL_RGB5, GL_BGRA, GL_UNSIGNED_SHORT_1_5_5_5_REV},
    //0x1D [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A4R4G4B4] =
        {2, linear, GL_RGBA4, GL_BGRA, GL_UNSIGNED_SHORT_4_4_4_4_REV}, // TODO : Verify this is truely linear
    //0x1E [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_X8R8G8B8] =
        {4, linear, GL_RGB8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV},
    //0x1F [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8] =
        {1, linear, GL_R8, GL_RED, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_111R},
    //0x20 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8Y8] =
        {2, linear, GL_RG8, GL_RG, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_GGGR},
	//0x21 [?] =
		{},
	//0x22 [?] =
		{},
	//0x23 [?] =
		{},
	//0x24 [NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_CR8YB8CB8YA8] = // See convert_texture_data calling ____UYVYToARGBRow_C
		{2, linear, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify
	//0x25 [NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_YB8CR8YA8CB8] = // See convert_texture_data calling ____YUY2ToARGBRow_C
		{2, linear, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify
	//0x26 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8CR8CB8Y8] = // See convert_texture_data FIXME
		{2, linear, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify

    //0x27 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R6G5B5] = // See convert_texture_data calling __R6G5B5ToARGBRow_C
        {2, swizzled, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify
    //0x28 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_G8B8] =
        {2, swizzled, GL_RG8, GL_RG, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_0RG1},
    //0x29 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R8B8] =
        {2, swizzled, GL_RG8, GL_RG, GL_UNSIGNED_BYTE,
         gl_swizzle_mask_R0G1},
	//0x2A [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_X8_Y24_FIXED] =
		{4, swizzled, GL_DEPTH24_STENCIL8, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8}, // TODO : Verify
	//0x2B [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_X8_Y24_FLOAT] =
		{4, swizzled, GL_DEPTH24_STENCIL8, GL_DEPTH_STENCIL, GL_FLOAT_32_UNSIGNED_INT_24_8_REV}, // TODO : Verify
	//0x2C [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_Y16_FIXED] =
		{2, swizzled, GL_DEPTH_COMPONENT16, GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT}, // TODO : Verify
	//0x2D [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_Y16_FLOAT] =
		{2, swizzled, GL_DEPTH_COMPONENT16, GL_DEPTH_COMPONENT, GL_FLOAT}, // TODO : Verify


    /* TODO: format conversion */
    //0x2E [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_X8_Y24_FIXED] =
        {4, linear, GL_DEPTH24_STENCIL8, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8},
	//0x2F [?] =
		{},
    //0x30 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_Y16_FIXED] =
        {2, linear, GL_DEPTH_COMPONENT16, GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT},
	//0x31 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_Y16_FLOAT] =
        {2, linear, GL_DEPTH_COMPONENT16, GL_DEPTH_COMPONENT, GL_FLOAT}, // TODO : Verify
	//0x32 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_Y16] =
        {2, swizzled, GL_R16, GL_RED, GL_UNSIGNED_SHORT, // TODO : Verify
         gl_swizzle_mask_RRR1},
	//0x33 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_YB16YA16] =
        {4, swizzled, GL_RG16, GL_RG, GL_UNSIGNED_SHORT, // TODO : Verify
         gl_swizzle_mask_RRRG},
	//0x34 [NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_A4V6YB6A4U6YA6] = // TODO : handle in convert_texture_data
		{2, linear, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify
    //0x35 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_Y16] =
        {2, linear, GL_R16, GL_RED, GL_UNSIGNED_SHORT,
         gl_swizzle_mask_RRR1},
	//0x36 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_YB16YA16] =
        {4, linear, GL_RG16, GL_RG, GL_UNSIGNED_SHORT, // TODO : Verify
         gl_swizzle_mask_RRRG},
	//0x37 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R6G5B5] = // See convert_texture_data calling __R6G5B5ToARGBRow_C
        {2, linear, GL_CONVERT_TEXTURE_DATA_RESULTING_FORMAT}, // TODO : Verify
	//0x38 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R5G5B5A1] =
		{2, swizzled, GL_RGB5_A1, GL_RGBA, GL_UNSIGNED_SHORT_5_5_5_1}, // TODO : Verify
	//0x39 [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R4G4B4A4] =
		{2, swizzled, GL_RGBA, GL_RGBA, GL_UNSIGNED_SHORT_4_4_4_4}, // TODO : Verify
    //0x3A [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_A8B8G8R8] =
        {4, swizzled, GL_RGBA8, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8_REV}, // TODO : Verify
	//0x3B [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_B8G8R8A8] =
		{4, swizzled, GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8}, // TODO : Verify

    //0x3C [NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R8G8B8A8] =
        {4, swizzled, GL_RGBA8, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8},
	//0x3D [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R5G5B5A1] =
		{2, linear, GL_RGB5_A1, GL_RGBA, GL_UNSIGNED_SHORT_5_5_5_1}, // TODO : Verify
	//0x3E [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R4G4B4A4] =
        {2, linear, GL_RGBA, GL_RGBA, GL_UNSIGNED_SHORT_4_4_4_4}, // TODO : Verify

    //0x3F [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8B8G8R8] =
        {4, linear, GL_RGBA8, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8_REV}, // TODO : Verify
    //0x40 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_B8G8R8A8] =
        {4, linear, GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8}, // TODO : Verify
    //0x41 [NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R8G8B8A8] =
        {4, linear, GL_RGBA8, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8}, // TODO : Verify
};

static void OpenGL_pgraph_update_shader_constants(PGRAPHState *pg, ShaderBinding *binding, bool binding_changed, bool vertex_program, bool fixed_function);
static void OpenGL_pgraph_bind_shaders(PGRAPHState *pg);
static void OpenGL_pgraph_update_surface_part(NV2AState *d, bool upload, bool color);
static void OpenGL_pgraph_update_surface(NV2AState *d, bool upload, bool color_write, bool zeta_write);
static void OpenGL_pgraph_bind_textures(NV2AState *d);
static void OpenGL_pgraph_update_memory_buffer(NV2AState *d, hwaddr addr, hwaddr size, bool f);
static void OpenGL_pgraph_bind_vertex_attributes(NV2AState *d, unsigned int num_elements, bool inline_data, unsigned int inline_stride);
static unsigned int OpenGL_pgraph_bind_inline_array(NV2AState *d);

static int OpenGL_upload_texture(GLenum gl_target, const TextureShape s, const uint8_t *texture_data, const uint8_t *palette_data);
static TextureBinding* OpenGL_generate_texture(const TextureShape s, const uint8_t *texture_data, const uint8_t *palette_data);
static gpointer OpenGL_texture_key_retrieve(gpointer key, gpointer user_data, GError **error);
static void OpenGL_texture_binding_destroy(gpointer data);

void OpenGL_draw_end(NV2AState *d); // forward declaration

void OpenGL_draw_arrays(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);
	assert(pg->shader_binding);

	OpenGL_pgraph_bind_vertex_attributes(d, pg->draw_arrays_max_count,
		false, 0);
	glMultiDrawArrays(pg->shader_binding->gl_primitive_mode,
		pg->gl_draw_arrays_start,
		pg->gl_draw_arrays_count,
		pg->draw_arrays_length);

	OpenGL_draw_end(d);
}

void OpenGL_draw_inline_buffer(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);
	assert(pg->shader_binding);

	for (unsigned int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
		VertexAttribute *vertex_attribute = &pg->vertex_attributes[i];

		if (vertex_attribute->inline_buffer) {

			glBindBuffer(GL_ARRAY_BUFFER,
				vertex_attribute->gl_inline_buffer);
			glBufferData(GL_ARRAY_BUFFER,
				pg->inline_buffer_length
				* sizeof(float) * 4,
				vertex_attribute->inline_buffer,
				GL_DYNAMIC_DRAW);

			/* Clear buffer for next batch */
			g_free(vertex_attribute->inline_buffer);
			vertex_attribute->inline_buffer = NULL;

			glVertexAttribPointer(i, 4, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(i);
		}
		else {
			glDisableVertexAttribArray(i);
			glVertexAttrib4fv(i, vertex_attribute->inline_value);
		}
	}

	glDrawArrays(pg->shader_binding->gl_primitive_mode,
		0, pg->inline_buffer_length);

	OpenGL_draw_end(d);
}

void OpenGL_draw_inline_array(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);
	assert(pg->shader_binding);

	unsigned int index_count = OpenGL_pgraph_bind_inline_array(d);
	glDrawArrays(pg->shader_binding->gl_primitive_mode,
		0, index_count);

	OpenGL_draw_end(d);
}

void OpenGL_draw_inline_elements(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);
	assert(pg->shader_binding);

	uint32_t max_element = 0;
	uint32_t min_element = (uint32_t)-1;
	for (unsigned int i = 0; i < pg->inline_elements_length; i++) {
		max_element = MAX(pg->inline_elements[i], max_element);
		min_element = MIN(pg->inline_elements[i], min_element);
	}
	OpenGL_pgraph_bind_vertex_attributes(d, max_element + 1, false, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pg->gl_element_buffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		pg->inline_elements_length * sizeof(pg->inline_elements[0]),
		pg->inline_elements,
		GL_DYNAMIC_DRAW);

	glDrawRangeElements(pg->shader_binding->gl_primitive_mode,
		min_element, max_element,
		pg->inline_elements_length,
		GL_UNSIGNED_SHORT, // Cxbx-Reloaded TODO : Restore GL_UNSIGNED_INT once HLE_draw_inline_elements can draw using uint32_t
		(void*)0);

	OpenGL_draw_end(d);
}

void OpenGL_draw_state_update(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

	NV2A_GL_DGROUP_BEGIN("NV097_SET_BEGIN_END: 0x%x", pg->primitive_mode);

	uint32_t control_0 = pg->regs[NV_PGRAPH_CONTROL_0];
	uint32_t control_1 = pg->regs[NV_PGRAPH_CONTROL_1];

	bool depth_test = control_0
		& NV_PGRAPH_CONTROL_0_ZENABLE;
	bool stencil_test = control_1
		& NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE;

	OpenGL_pgraph_update_surface(d, true, true, depth_test || stencil_test);

	bool alpha = control_0 & NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE;
	bool red = control_0 & NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE;
	bool green = control_0 & NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE;
	bool blue = control_0 & NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE;
	glColorMask(red, green, blue, alpha);
	glDepthMask(!!(control_0 & NV_PGRAPH_CONTROL_0_ZWRITEENABLE));
	glStencilMask(GET_MASK(control_1,
		NV_PGRAPH_CONTROL_1_STENCIL_MASK_WRITE));

	uint32_t blend = pg->regs[NV_PGRAPH_BLEND];
	if (blend & NV_PGRAPH_BLEND_EN) {
		static const GLenum pgraph_blend_factor_map[] = {
			GL_ZERO,
			GL_ONE,
			GL_SRC_COLOR,
			GL_ONE_MINUS_SRC_COLOR,
			GL_SRC_ALPHA,
			GL_ONE_MINUS_SRC_ALPHA,
			GL_DST_ALPHA,
			GL_ONE_MINUS_DST_ALPHA,
			GL_DST_COLOR,
			GL_ONE_MINUS_DST_COLOR,
			GL_SRC_ALPHA_SATURATE,
			0,
			GL_CONSTANT_COLOR,
			GL_ONE_MINUS_CONSTANT_COLOR,
			GL_CONSTANT_ALPHA,
			GL_ONE_MINUS_CONSTANT_ALPHA,
		};

		glEnable(GL_BLEND);
		uint32_t sfactor = GET_MASK(blend,
			NV_PGRAPH_BLEND_SFACTOR);
		uint32_t dfactor = GET_MASK(blend,
			NV_PGRAPH_BLEND_DFACTOR);
		assert(sfactor < ARRAY_SIZE(pgraph_blend_factor_map));
		assert(dfactor < ARRAY_SIZE(pgraph_blend_factor_map));
		glBlendFunc(pgraph_blend_factor_map[sfactor],
			pgraph_blend_factor_map[dfactor]);

		{
			static const GLenum pgraph_blend_equation_map[] = {
				GL_FUNC_SUBTRACT,
				GL_FUNC_REVERSE_SUBTRACT,
				GL_FUNC_ADD,
				GL_MIN,
				GL_MAX,
				GL_FUNC_REVERSE_SUBTRACT,
				GL_FUNC_ADD,
			};

			uint32_t equation = GET_MASK(blend,
				NV_PGRAPH_BLEND_EQN);
			assert(equation < ARRAY_SIZE(pgraph_blend_equation_map));
			glBlendEquation(pgraph_blend_equation_map[equation]);
		}

/*
			static const GLenum pgraph_blend_logicop_map[] = {
				GL_CLEAR,
				GL_AND,
				GL_AND_REVERSE,
				GL_COPY,
				GL_AND_INVERTED,
				GL_NOOP,
				GL_XOR,
				GL_OR,
				GL_NOR,
				GL_EQUIV,
				GL_INVERT,
				GL_OR_REVERSE,
				GL_COPY_INVERTED,
				GL_OR_INVERTED,
				GL_NAND,
				GL_SET,
			};


*/
		uint32_t blend_color = pg->regs[NV_PGRAPH_BLENDCOLOR];
		glBlendColor(((blend_color >> 16) & 0xFF) / 255.0f, /* red */
			((blend_color >> 8) & 0xFF) / 255.0f,  /* green */
			(blend_color & 0xFF) / 255.0f,         /* blue */
			((blend_color >> 24) & 0xFF) / 255.0f);/* alpha */
	}
	else {
		glDisable(GL_BLEND);
	}

	/* Face culling */
	uint32_t setupraster = pg->regs[NV_PGRAPH_SETUPRASTER];
	if (setupraster
		& NV_PGRAPH_SETUPRASTER_CULLENABLE) {
		static const GLenum pgraph_cull_face_map[] = {
			0,
			GL_FRONT,
			GL_BACK,
			GL_FRONT_AND_BACK
		};

		uint32_t cull_face = GET_MASK(setupraster,
			NV_PGRAPH_SETUPRASTER_CULLCTRL);
		assert(cull_face < ARRAY_SIZE(pgraph_cull_face_map));
		glCullFace(pgraph_cull_face_map[cull_face]);
		glEnable(GL_CULL_FACE);
	}
	else {
		glDisable(GL_CULL_FACE);
	}

	/* Front-face select */
	glFrontFace(setupraster
		& NV_PGRAPH_SETUPRASTER_FRONTFACE
		? GL_CCW : GL_CW);

	/* Polygon offset */
	/* FIXME: GL implementation-specific, maybe do this in VS? */
	if (setupraster &
		NV_PGRAPH_SETUPRASTER_POFFSETFILLENABLE) {
		glEnable(GL_POLYGON_OFFSET_FILL);
	}
	else {
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	if (setupraster &
		NV_PGRAPH_SETUPRASTER_POFFSETLINEENABLE) {
		glEnable(GL_POLYGON_OFFSET_LINE);
	}
	else {
		glDisable(GL_POLYGON_OFFSET_LINE);
	}
	if (setupraster &
		NV_PGRAPH_SETUPRASTER_POFFSETPOINTENABLE) {
		glEnable(GL_POLYGON_OFFSET_POINT);
	}
	else {
		glDisable(GL_POLYGON_OFFSET_POINT);
	}
	if (setupraster &
		(NV_PGRAPH_SETUPRASTER_POFFSETFILLENABLE |
			NV_PGRAPH_SETUPRASTER_POFFSETLINEENABLE |
			NV_PGRAPH_SETUPRASTER_POFFSETPOINTENABLE)) {
		GLfloat zfactor = *(float*)&pg->regs[NV_PGRAPH_ZOFFSETFACTOR];
		GLfloat zbias = *(float*)&pg->regs[NV_PGRAPH_ZOFFSETBIAS];
		glPolygonOffset(zfactor, zbias);
	}

	/* Depth testing */
	if (depth_test) {
		static const GLenum pgraph_depth_func_map[] = {
			GL_NEVER,
			GL_LESS,
			GL_EQUAL,
			GL_LEQUAL,
			GL_GREATER,
			GL_NOTEQUAL,
			GL_GEQUAL,
			GL_ALWAYS,
		};

		glEnable(GL_DEPTH_TEST);

		uint32_t depth_func = GET_MASK(control_0,
			NV_PGRAPH_CONTROL_0_ZFUNC);
		assert(depth_func < ARRAY_SIZE(pgraph_depth_func_map));
		glDepthFunc(pgraph_depth_func_map[depth_func]);
	}
	else {
		glDisable(GL_DEPTH_TEST);
	}

	if (stencil_test) {
		static const GLenum pgraph_stencil_func_map[] = {
			GL_NEVER,
			GL_LESS,
			GL_EQUAL,
			GL_LEQUAL,
			GL_GREATER,
			GL_NOTEQUAL,
			GL_GEQUAL,
			GL_ALWAYS,
		};

		glEnable(GL_STENCIL_TEST);

		uint32_t stencil_func = GET_MASK(control_1,
			NV_PGRAPH_CONTROL_1_STENCIL_FUNC);
		uint32_t stencil_ref = GET_MASK(control_1,
			NV_PGRAPH_CONTROL_1_STENCIL_REF);
		uint32_t func_mask = GET_MASK(control_1,
			NV_PGRAPH_CONTROL_1_STENCIL_MASK_READ);
		uint32_t control2 = pg->regs[NV_PGRAPH_CONTROL_2];
		uint32_t op_fail = GET_MASK(control2,
			NV_PGRAPH_CONTROL_2_STENCIL_OP_FAIL);
		uint32_t op_zfail = GET_MASK(control2,
			NV_PGRAPH_CONTROL_2_STENCIL_OP_ZFAIL);
		uint32_t op_zpass = GET_MASK(control2,
			NV_PGRAPH_CONTROL_2_STENCIL_OP_ZPASS);

		assert(stencil_func < ARRAY_SIZE(pgraph_stencil_func_map));

		static const GLenum pgraph_stencil_op_map[] = {
			0,
			GL_KEEP,
			GL_ZERO,
			GL_REPLACE,
			GL_INCR,
			GL_DECR,
			GL_INVERT,
			GL_INCR_WRAP,
			GL_DECR_WRAP,
		};

		assert(op_fail < ARRAY_SIZE(pgraph_stencil_op_map));
		assert(op_zfail < ARRAY_SIZE(pgraph_stencil_op_map));
		assert(op_zpass < ARRAY_SIZE(pgraph_stencil_op_map));

		glStencilFunc(
			pgraph_stencil_func_map[stencil_func],
			stencil_ref,
			func_mask);

		glStencilOp(
			pgraph_stencil_op_map[op_fail],
			pgraph_stencil_op_map[op_zfail],
			pgraph_stencil_op_map[op_zpass]);

	}
	else {
		glDisable(GL_STENCIL_TEST);
	}

	/* Dither */
	/* FIXME: GL implementation dependent */
	if (control_0 &
		NV_PGRAPH_CONTROL_0_DITHERENABLE) {
		glEnable(GL_DITHER);
	}
	else {
		glDisable(GL_DITHER);
	}

	OpenGL_pgraph_bind_shaders(pg);
	OpenGL_pgraph_bind_textures(d);

	//glDisableVertexAttribArray(NV2A_VERTEX_ATTR_DIFFUSE);
	//glVertexAttrib4f(NV2A_VERTEX_ATTR_DIFFUSE, 1.0f, 1.0f, 1.0f, 1.0f);

	unsigned int width, height;
	pgraph_get_surface_dimensions(pg, &width, &height);
	pgraph_apply_anti_aliasing_factor(pg, &width, &height);
	glViewport(0, 0, width, height);

	/* Visibility testing */
	if (pg->zpass_pixel_count_enable) {
		GLuint gl_query;
		glGenQueries(1, &gl_query);
		pg->gl_zpass_pixel_count_query_count++;
		pg->gl_zpass_pixel_count_queries = (GLuint*)g_realloc(
			pg->gl_zpass_pixel_count_queries,
			sizeof(GLuint) * pg->gl_zpass_pixel_count_query_count);
		pg->gl_zpass_pixel_count_queries[
			pg->gl_zpass_pixel_count_query_count - 1] = gl_query;
		glBeginQuery(GL_SAMPLES_PASSED, gl_query);
	}
}

void OpenGL_draw_end(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

	/* End of visibility testing */
	if (pg->zpass_pixel_count_enable) {
		glEndQuery(GL_SAMPLES_PASSED);
	}

	NV2A_GL_DGROUP_END();
}

void OpenGL_draw_clear(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

	NV2A_DPRINTF("---------PRE CLEAR ------\n");
	GLbitfield gl_mask = 0;

	bool write_color = (pg->clear_surface & NV097_CLEAR_SURFACE_COLOR);
	bool write_zeta =
		(pg->clear_surface & (NV097_CLEAR_SURFACE_Z | NV097_CLEAR_SURFACE_STENCIL));

	if (write_zeta) {
		uint32_t clear_zstencil =
			d->pgraph.regs[NV_PGRAPH_ZSTENCILCLEARVALUE];
		GLint gl_clear_stencil;
		GLfloat gl_clear_depth;

		/* FIXME: Put these in some lookup table */
		const float f16_max = 511.9375f;
		/* FIXME: 7 bits of mantissa unused. maybe use full buffer? */
		const float f24_max = 3.4027977E38f;

		switch (pg->surface_shape.zeta_format) {
		case NV097_SET_SURFACE_FORMAT_ZETA_Z16: {
			if (pg->clear_surface & NV097_CLEAR_SURFACE_Z) {
				gl_mask |= GL_DEPTH_BUFFER_BIT;
				uint16_t z = clear_zstencil & 0xFFFF;
				if (pg->surface_shape.z_format) {
					gl_clear_depth = convert_f16_to_float(z) / f16_max;
					assert(false); /* FIXME: Untested */
				}
				else {
					gl_clear_depth = z / (float)0xFFFF;
				}
			}
			break;
		}
		case NV097_SET_SURFACE_FORMAT_ZETA_Z24S8: {
			if (pg->clear_surface & NV097_CLEAR_SURFACE_STENCIL) {
				gl_mask |= GL_STENCIL_BUFFER_BIT;
				gl_clear_stencil = clear_zstencil & 0xFF;
			}
			if (pg->clear_surface & NV097_CLEAR_SURFACE_Z) {
				gl_mask |= GL_DEPTH_BUFFER_BIT;
				uint32_t z = clear_zstencil >> 8;
				if (pg->surface_shape.z_format) {
					gl_clear_depth = convert_f24_to_float(z) / f24_max;
					assert(false); /* FIXME: Untested */
				}
				else {
					gl_clear_depth = z / (float)0xFFFFFF;
				}
			}
			break;
		}
		default:
			fprintf(stderr, "Unknown zeta surface format: 0x%x\n", pg->surface_shape.zeta_format);
			assert(false);
			break;
		}

		if (gl_mask & GL_DEPTH_BUFFER_BIT) {
			glDepthMask(GL_TRUE);
			glClearDepth(gl_clear_depth);
		}

		if (gl_mask & GL_STENCIL_BUFFER_BIT) {
			glStencilMask(0xff);
			glClearStencil(gl_clear_stencil);
		}
	}

	if (write_color) {
		gl_mask |= GL_COLOR_BUFFER_BIT;
		glColorMask((pg->clear_surface & NV097_CLEAR_SURFACE_R)
			? GL_TRUE : GL_FALSE,
			(pg->clear_surface & NV097_CLEAR_SURFACE_G)
			? GL_TRUE : GL_FALSE,
			(pg->clear_surface & NV097_CLEAR_SURFACE_B)
			? GL_TRUE : GL_FALSE,
			(pg->clear_surface & NV097_CLEAR_SURFACE_A)
			? GL_TRUE : GL_FALSE);
		uint32_t clear_color = d->pgraph.regs[NV_PGRAPH_COLORCLEARVALUE];

		/* Handle RGB */
		GLfloat red, green, blue;
		switch (pg->surface_shape.color_format) {
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1R5G5B5_Z1R5G5B5:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1R5G5B5_O1R5G5B5:
			red = ((clear_color >> 10) & 0x1F) / 31.0f;
			green = ((clear_color >> 5) & 0x1F) / 31.0f;
			blue = (clear_color & 0x1F) / 31.0f;
			assert(false); /* Untested */
			break;
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_R5G6B5:
			red = ((clear_color >> 11) & 0x1F) / 31.0f;
			green = ((clear_color >> 5) & 0x3F) / 63.0f;
			blue = (clear_color & 0x1F) / 31.0f;
			break;
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X8R8G8B8_Z8R8G8B8:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X8R8G8B8_O8R8G8B8:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_Z1A7R8G8B8:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_O1A7R8G8B8:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_A8R8G8B8:
			red = ((clear_color >> 16) & 0xFF) / 255.0f;
			green = ((clear_color >> 8) & 0xFF) / 255.0f;
			blue = (clear_color & 0xFF) / 255.0f;
			break;
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_B8:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_G8B8:
			/* Xbox D3D doesn't support clearing those */
		default:
			red = 1.0f;
			green = 0.0f;
			blue = 1.0f;
			fprintf(stderr, "CLEAR_SURFACE for color_format 0x%x unsupported",
				pg->surface_shape.color_format);
			assert(false);
			break;
		}

		/* Handle alpha */
		GLfloat alpha;
		switch (pg->surface_shape.color_format) {
			/* FIXME: CLEAR_SURFACE seems to work like memset, so maybe we
			*        also have to clear non-alpha bits with alpha value?
			*        As GL doesn't own those pixels we'd have to do this on
			*        our own in xbox memory.
			*/
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_Z1A7R8G8B8:
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_O1A7R8G8B8:
			alpha = ((clear_color >> 24) & 0x7F) / 127.0f;
			assert(false); /* Untested */
			break;
		case NV097_SET_SURFACE_FORMAT_COLOR_LE_A8R8G8B8:
			alpha = ((clear_color >> 24) & 0xFF) / 255.0f;
			break;
		default:
			alpha = 1.0f;
			break;
		}

		glClearColor(red, green, blue, alpha);
	}

	if (gl_mask) {
		OpenGL_pgraph_update_surface(d, true, write_color, write_zeta);

		glEnable(GL_SCISSOR_TEST);

		unsigned int xmin = GET_MASK(pg->regs[NV_PGRAPH_CLEARRECTX],
			NV_PGRAPH_CLEARRECTX_XMIN);
		unsigned int xmax = GET_MASK(pg->regs[NV_PGRAPH_CLEARRECTX],
			NV_PGRAPH_CLEARRECTX_XMAX);
		unsigned int ymin = GET_MASK(pg->regs[NV_PGRAPH_CLEARRECTY],
			NV_PGRAPH_CLEARRECTY_YMIN);
		unsigned int ymax = GET_MASK(pg->regs[NV_PGRAPH_CLEARRECTY],
			NV_PGRAPH_CLEARRECTY_YMAX);

		unsigned int scissor_x = xmin;
		unsigned int scissor_y = pg->surface_shape.clip_height - ymax - 1;

		unsigned int scissor_width = xmax - xmin + 1;
		unsigned int scissor_height = ymax - ymin + 1;

		pgraph_apply_anti_aliasing_factor(pg, &scissor_x, &scissor_y);
		pgraph_apply_anti_aliasing_factor(pg, &scissor_width, &scissor_height);

		/* FIXME: Should this really be inverted instead of ymin? */
		glScissor(scissor_x, scissor_y, scissor_width, scissor_height);

		/* FIXME: Respect window clip?!?! */

		NV2A_DPRINTF("------------------CLEAR 0x%x %d,%d - %d,%d  %x---------------\n",
			parameter, xmin, ymin, xmax, ymax, d->pgraph.regs[NV_PGRAPH_COLORCLEARVALUE]);

		/* Dither */
		/* FIXME: Maybe also disable it here? + GL implementation dependent */
		if (pg->regs[NV_PGRAPH_CONTROL_0] &
			NV_PGRAPH_CONTROL_0_DITHERENABLE) {
			glEnable(GL_DITHER);
		}
		else {
			glDisable(GL_DITHER);
		}

		glClear(gl_mask);

		glDisable(GL_SCISSOR_TEST);
	}

	pgraph_set_surface_dirty(pg, write_color, write_zeta);
}

void OpenGL_pgraph_report_get(NV2AState *d)
{
	unsigned int i;

	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

	/* FIXME: Multisampling affects this (both: OGL and Xbox GPU),
	 *        not sure if CLEARs also count
	 */
	 /* FIXME: What about clipping regions etc? */
	for (i = 0; i < pg->gl_zpass_pixel_count_query_count; i++) {
		GLuint gl_query_result;
		glGetQueryObjectuiv(pg->gl_zpass_pixel_count_queries[i],
			GL_QUERY_RESULT,
			&gl_query_result);
		pg->zpass_pixel_count_result += gl_query_result;
	}
	if (pg->gl_zpass_pixel_count_query_count) {
		glDeleteQueries(pg->gl_zpass_pixel_count_query_count,
			pg->gl_zpass_pixel_count_queries);
	}
	pg->gl_zpass_pixel_count_query_count = 0;
}

void OpenGL_pgraph_report_clear(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

	/* FIXME: Does this have a value in parameter? Also does this (also?) modify
	 *        the report memory block?
	 */
	if (pg->gl_zpass_pixel_count_query_count) {
		glDeleteQueries(pg->gl_zpass_pixel_count_query_count,
			pg->gl_zpass_pixel_count_queries);
	}
	pg->gl_zpass_pixel_count_query_count = 0;
}

void OpenGL_init_pgraph_plugins()
{
	pgraph_draw_arrays = OpenGL_draw_arrays;
	pgraph_draw_inline_buffer = OpenGL_draw_inline_buffer;
	pgraph_draw_inline_array = OpenGL_draw_inline_array;
	pgraph_draw_inline_elements = OpenGL_draw_inline_elements;
	pgraph_draw_state_update = OpenGL_draw_state_update;
	pgraph_draw_clear = OpenGL_draw_clear;
	pgraph_update_surface = OpenGL_pgraph_update_surface;
	pgraph_report_get = OpenGL_pgraph_report_get;
	pgraph_report_clear = OpenGL_pgraph_report_clear;
}

void OpenGL_init(PGRAPHState *pg)
{
	int i;

	/* attach OpenGL render plugins */
	OpenGL_init_pgraph_plugins();

	/* fire up opengl */

    pg->gl_context = glo_context_create();
    assert(pg->gl_context);

#ifdef DEBUG_NV2A_GL
    glEnable(GL_DEBUG_OUTPUT);
#endif

    glextensions_init();

    /* DXT textures */
    assert(glo_check_extension("GL_EXT_texture_compression_s3tc"));
    /*  Internal RGB565 texture format */
    assert(glo_check_extension("GL_ARB_ES2_compatibility"));

    GLint max_vertex_attributes;
    glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &max_vertex_attributes);
    assert(max_vertex_attributes >= NV2A_VERTEXSHADER_ATTRIBUTES);

    glGenFramebuffers(1, &pg->gl_framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, pg->gl_framebuffer);

    /* need a valid framebuffer to start with */
    glGenTextures(1, &pg->gl_color_buffer);
    glBindTexture(GL_TEXTURE_2D, pg->gl_color_buffer);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 640, 480,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, pg->gl_color_buffer, 0);

    assert(glCheckFramebufferStatus(GL_FRAMEBUFFER)
            == GL_FRAMEBUFFER_COMPLETE);

    //glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

#ifdef USE_TEXTURE_CACHE
    pg->texture_cache = g_lru_cache_new_full(
        0,
        NULL,
        texture_key_destroy,
        0,
        NULL,
		OpenGL_texture_binding_destroy,
        texture_key_hash,
        texture_key_equal,
        OpenGL_texture_key_retrieve,
        NULL,
        NULL
        );

    g_lru_cache_set_max_size(pg->texture_cache, 512);
#endif

#ifdef USE_SHADER_CACHE
    pg->shader_cache = g_hash_table_new(shader_hash, shader_equal);
#endif

    for (i=0; i<NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
        glGenBuffers(1, &pg->vertex_attributes[i].gl_converted_buffer);
        glGenBuffers(1, &pg->vertex_attributes[i].gl_inline_buffer);
    }
    glGenBuffers(1, &pg->gl_inline_array_buffer);
    glGenBuffers(1, &pg->gl_element_buffer);

    glGenBuffers(1, &pg->gl_memory_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, pg->gl_memory_buffer);
    glBufferData(GL_ARRAY_BUFFER,
                 g_SystemMaxMemory, // Was : d->vram_size,
                 NULL,
                 GL_DYNAMIC_DRAW);

    glGenVertexArrays(1, &pg->gl_vertex_array);
    glBindVertexArray(pg->gl_vertex_array);

//    assert(glGetError() == GL_NO_ERROR);

    glo_set_current(NULL);
}

void OpenGL_destroy(PGRAPHState *pg)
{
	glo_set_current(pg->gl_context);

	if (pg->gl_color_buffer) {
		glDeleteTextures(1, &pg->gl_color_buffer);
	}
	if (pg->gl_zeta_buffer) {
		glDeleteTextures(1, &pg->gl_zeta_buffer);
	}
	glDeleteFramebuffers(1, &pg->gl_framebuffer);

	// TODO: clear out shader cached
	// TODO: clear out texture cache

	glo_set_current(NULL);

	glo_context_destroy(pg->gl_context);
}

static void OpenGL_pgraph_update_shader_constants(PGRAPHState *pg,
                                           ShaderBinding *binding,
                                           bool binding_changed,
                                           bool vertex_program,
                                           bool fixed_function)
{
	assert(pg->opengl_enabled);

	unsigned int i, j;

    /* update combiner constants */
    for (i = 0; i<= 8; i++) {
        uint32_t constant[2];
        if (i == 8) {
            /* final combiner */
            constant[0] = pg->regs[NV_PGRAPH_SPECFOGFACTOR0];
            constant[1] = pg->regs[NV_PGRAPH_SPECFOGFACTOR1];
        } else {
            constant[0] = pg->regs[NV_PGRAPH_COMBINEFACTOR0 + i * 4];
            constant[1] = pg->regs[NV_PGRAPH_COMBINEFACTOR1 + i * 4];
        }

        for (j = 0; j < 2; j++) {
            GLint loc = binding->psh_constant_loc[i][j];
            if (loc != -1) {
                float value[4];
                value[0] = (float) ((constant[j] >> 16) & 0xFF) / 255.0f;
                value[1] = (float) ((constant[j] >> 8) & 0xFF) / 255.0f;
                value[2] = (float) (constant[j] & 0xFF) / 255.0f;
                value[3] = (float) ((constant[j] >> 24) & 0xFF) / 255.0f;

                glUniform4fv(loc, 1, value);
            }
        }
    }
    if (binding->alpha_ref_loc != -1) {
        float alpha_ref = GET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
                                   NV_PGRAPH_CONTROL_0_ALPHAREF) / 255.0f;
        glUniform1f(binding->alpha_ref_loc, alpha_ref);
    }


    /* For each texture stage */
    for (i = 0; i < NV2A_MAX_TEXTURES; i++) {
        // char name[32];
        GLint loc;

        /* Bump luminance only during stages 1 - 3 */
        if (i > 0) {
            loc = binding->bump_mat_loc[i];
            if (loc != -1) {
                glUniformMatrix2fv(loc, 1, GL_FALSE, pg->bump_env_matrix[i - 1]);
            }
            loc = binding->bump_scale_loc[i];
            if (loc != -1) {
                glUniform1f(loc, *(float*)&pg->regs[
                                NV_PGRAPH_BUMPSCALE1 + (i - 1) * 4]);
            }
            loc = binding->bump_offset_loc[i];
            if (loc != -1) {
                glUniform1f(loc, *(float*)&pg->regs[
                            NV_PGRAPH_BUMPOFFSET1 + (i - 1) * 4]);
            }
        }

    }

    if (binding->fog_color_loc != -1) {
        uint32_t fog_color = pg->regs[NV_PGRAPH_FOGCOLOR];
        glUniform4f(binding->fog_color_loc,
                    GET_MASK(fog_color, NV_PGRAPH_FOGCOLOR_RED) / 255.0f,
                    GET_MASK(fog_color, NV_PGRAPH_FOGCOLOR_GREEN) / 255.0f,
                    GET_MASK(fog_color, NV_PGRAPH_FOGCOLOR_BLUE) / 255.0f,
                    GET_MASK(fog_color, NV_PGRAPH_FOGCOLOR_ALPHA) / 255.0f);
    }
    if (binding->fog_param_loc[0] != -1) {
        glUniform1f(binding->fog_param_loc[0],
                    *(float*)&pg->regs[NV_PGRAPH_FOGPARAM0]);
    }
    if (binding->fog_param_loc[1] != -1) {
        glUniform1f(binding->fog_param_loc[1],
                    *(float*)&pg->regs[NV_PGRAPH_FOGPARAM1]);
    }

    float zclip_max = *(float*)&pg->regs[NV_PGRAPH_ZCLIPMAX];
    float zclip_min = *(float*)&pg->regs[NV_PGRAPH_ZCLIPMIN];

    if (fixed_function) {
        /* update lighting constants */
        struct {
            uint32_t* v;
            bool* dirty;
            GLint* locs;
            size_t len;
        } lighting_arrays[] = {
			// TODO : Change pointers into offset_of(), so this variable can become static
            {&pg->ltctxa[0][0], &pg->ltctxa_dirty[0], binding->ltctxa_loc, NV2A_LTCTXA_COUNT},
            {&pg->ltctxb[0][0], &pg->ltctxb_dirty[0], binding->ltctxb_loc, NV2A_LTCTXB_COUNT},
            {&pg->ltc1[0][0], &pg->ltc1_dirty[0], binding->ltc1_loc, NV2A_LTC1_COUNT},
        };

        for (i=0; i<ARRAY_SIZE(lighting_arrays); i++) {
            uint32_t *lighting_v = lighting_arrays[i].v;
            bool *lighting_dirty = lighting_arrays[i].dirty;
            GLint *lighting_locs = lighting_arrays[i].locs;
            size_t lighting_len = lighting_arrays[i].len;
            for (j=0; j<lighting_len; j++) {
                if (!lighting_dirty[j] && !binding_changed) continue;
                GLint loc = lighting_locs[j];
                if (loc != -1) {
                    glUniform4fv(loc, 1, (const GLfloat*)&lighting_v[j*4]);
                }
                lighting_dirty[j] = false;
            }
        }


        for (i = 0; i < NV2A_MAX_LIGHTS; i++) {
            GLint loc;
            loc = binding->light_infinite_half_vector_loc[i];
            if (loc != -1) {
                glUniform3fv(loc, 1, pg->light_infinite_half_vector[i]);
            }
            loc = binding->light_infinite_direction_loc[i];
            if (loc != -1) {
                glUniform3fv(loc, 1, pg->light_infinite_direction[i]);
            }

            loc = binding->light_local_position_loc[i];
            if (loc != -1) {
                glUniform3fv(loc, 1, pg->light_local_position[i]);
            }
            loc = binding->light_local_attenuation_loc[i];
            if (loc != -1) {
                glUniform3fv(loc, 1, pg->light_local_attenuation[i]);
            }
        }

        /* estimate the viewport by assuming it matches the surface ... */
        //FIXME: Get surface dimensions?
        float m11 = 0.5f * pg->surface_shape.clip_width;
        float m22 = -0.5f * pg->surface_shape.clip_height;
        float m33 = zclip_max - zclip_min;
        //float m41 = m11;
        //float m42 = -m22;
        float m43 = zclip_min;
        //float m44 = 1.0f;

        if (m33 == 0.0f) {
            m33 = 1.0f;
        }
        float invViewport[16] = {
            1.0f/m11, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f/m22, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f/m33, 0.0f,
            -1.0f, 1.0f, -m43/m33, 1.0f
        };

        if (binding->inv_viewport_loc != -1) {
            glUniformMatrix4fv(binding->inv_viewport_loc,
                               1, GL_FALSE, &invViewport[0]);
        }

    }

    /* update vertex program constants */
    for (i=0; i<NV2A_VERTEXSHADER_CONSTANTS; i++) {
        if (!pg->vsh_constants_dirty[i] && !binding_changed) continue;

        GLint loc = binding->vsh_constant_loc[i];
        //assert(loc != -1);
        if (loc != -1) {
            glUniform4fv(loc, 1, (const GLfloat*)pg->vsh_constants[i]);
        }
        pg->vsh_constants_dirty[i] = false;
    }

    if (binding->surface_size_loc != -1) {
        glUniform2f(binding->surface_size_loc, (GLfloat)pg->surface_shape.clip_width,
                    (GLfloat)pg->surface_shape.clip_height);
    }

    if (binding->clip_range_loc != -1) {
        glUniform2f(binding->clip_range_loc, zclip_min, zclip_max);
    }

}

static void OpenGL_pgraph_bind_shaders(PGRAPHState *pg)
{
	assert(pg->opengl_enabled);

	unsigned int i, j;

	uint32_t csv0_d = pg->regs[NV_PGRAPH_CSV0_D];
    bool vertex_program = GET_MASK(csv0_d,
                                   NV_PGRAPH_CSV0_D_MODE) == 2;

    bool fixed_function = GET_MASK(csv0_d,
                                   NV_PGRAPH_CSV0_D_MODE) == 0;

	uint32_t csv0_c = pg->regs[NV_PGRAPH_CSV0_C];
    int program_start = GET_MASK(csv0_c,
                                 NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START);

    NV2A_GL_DGROUP_BEGIN("%s (VP: %s FFP: %s)", __func__,
                         vertex_program ? "yes" : "no",
                         fixed_function ? "yes" : "no");

	ShaderBinding* old_binding = pg->shader_binding;

	ShaderState state;
	/* register combiner stuff */
	state.psh.window_clip_exclusive = pg->regs[NV_PGRAPH_SETUPRASTER]
                                       & NV_PGRAPH_SETUPRASTER_WINDOWCLIPTYPE,
	state.psh.combiner_control = pg->regs[NV_PGRAPH_COMBINECTL];
	state.psh.shader_stage_program = pg->regs[NV_PGRAPH_SHADERPROG];
	state.psh.other_stage_input = pg->regs[NV_PGRAPH_SHADERCTL];
	state.psh.final_inputs_0 = pg->regs[NV_PGRAPH_COMBINESPECFOG0];
	state.psh.final_inputs_1 = pg->regs[NV_PGRAPH_COMBINESPECFOG1];
	uint32_t control0 = pg->regs[NV_PGRAPH_CONTROL_0];

	state.psh.alpha_test = control0
		& NV_PGRAPH_CONTROL_0_ALPHATESTENABLE;
	state.psh.alpha_func = (enum PshAlphaFunc)GET_MASK(control0,
		NV_PGRAPH_CONTROL_0_ALPHAFUNC);

    /* fixed function stuff */
	state.skinning = (enum VshSkinning)GET_MASK(csv0_d,
		NV_PGRAPH_CSV0_D_SKIN);
	state.lighting = csv0_c
		& NV_PGRAPH_CSV0_C_LIGHTING;
	state.normalization = csv0_c
		& NV_PGRAPH_CSV0_C_NORMALIZATION_ENABLE;

	state.fixed_function = fixed_function;

    /* vertex program stuff */
	state.vertex_program = vertex_program;
	state.z_perspective = control0
		& NV_PGRAPH_CONTROL_0_Z_PERSPECTIVE_ENABLE;

    /* geometry shader stuff */
	state.primitive_mode = (enum ShaderPrimitiveMode)pg->primitive_mode;
	state.polygon_front_mode = (enum ShaderPolygonMode)GET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
		NV_PGRAPH_SETUPRASTER_FRONTFACEMODE);
	state.polygon_back_mode = (enum ShaderPolygonMode)GET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
		NV_PGRAPH_SETUPRASTER_BACKFACEMODE);

    state.program_length = 0;
    memset(state.program_data, 0, sizeof(state.program_data));

    if (vertex_program) {
        // copy in vertex program tokens
        for (i = program_start; i < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH; i++) {
            uint32_t *cur_token = (uint32_t*)&pg->program_data[i];
            memcpy(&state.program_data[state.program_length],
                   cur_token,
                   VSH_TOKEN_SIZE * sizeof(uint32_t));
            state.program_length++;

            if (vsh_get_field(cur_token, FLD_FINAL)) {
                break;
            }
        }
    }

    /* Texgen */
    for (i = 0; i < 4; i++) {
        unsigned int reg = (i < 2) ? NV_PGRAPH_CSV1_A : NV_PGRAPH_CSV1_B;
        for (j = 0; j < 4; j++) {
            unsigned int masks[] = {
                // NOTE: For some reason, Visual Studio thinks NV_PGRAPH_xxxx is signed integer. (possible bug?)
                (i % 2U) ? (unsigned int)NV_PGRAPH_CSV1_A_T1_S : (unsigned int)NV_PGRAPH_CSV1_A_T0_S,
                (i % 2U) ? (unsigned int)NV_PGRAPH_CSV1_A_T1_T : (unsigned int)NV_PGRAPH_CSV1_A_T0_T,
                (i % 2U) ? (unsigned int)NV_PGRAPH_CSV1_A_T1_R : (unsigned int)NV_PGRAPH_CSV1_A_T0_R,
                (i % 2U) ? (unsigned int)NV_PGRAPH_CSV1_A_T1_Q : (unsigned int)NV_PGRAPH_CSV1_A_T0_Q
            };
            state.texgen[i][j] = (enum VshTexgen)GET_MASK(pg->regs[reg], masks[j]);
        }
    }

    /* Fog */
    state.fog_enable = pg->regs[NV_PGRAPH_CONTROL_3]
                           & NV_PGRAPH_CONTROL_3_FOGENABLE;
    if (state.fog_enable) {
        /*FIXME: Use CSV0_D? */
        state.fog_mode = (enum VshFogMode)GET_MASK(pg->regs[NV_PGRAPH_CONTROL_3],
                                  NV_PGRAPH_CONTROL_3_FOG_MODE);
        state.foggen = (enum VshFoggen)GET_MASK(csv0_d,
                                NV_PGRAPH_CSV0_D_FOGGENMODE);
    } else {
        /* FIXME: Do we still pass the fogmode? */
        state.fog_mode = (enum VshFogMode)0;
        state.foggen = (enum VshFoggen)0;
    }

    /* Texture matrices */
    for (i = 0; i < 4; i++) {
        state.texture_matrix_enable[i] = pg->texture_matrix_enable[i];
    }

    /* Lighting */
    if (state.lighting) {
        for (i = 0; i < NV2A_MAX_LIGHTS; i++) {
            state.light[i] = (enum VshLight)GET_MASK(csv0_d,
                                      NV_PGRAPH_CSV0_D_LIGHT0 << (i * 2));
        }
    }

    /* Window clip
     *
     * Optimization note: very quickly check to ignore any repeated or zero-size
     * clipping regions. Note that if region number 7 is valid, but the rest are
     * not, we will still add all of them. Clip regions seem to be typically
     * front-loaded (meaning the first one or two regions are populated, and the
     * following are zeroed-out), so let's avoid adding any more complicated
     * masking or copying logic here for now unless we discover a valid case.
     */
    assert(!state.psh.window_clip_exclusive); /* FIXME: Untested */
    state.psh.window_clip_count = 0;
    uint32_t last_x = 0, last_y = 0;

    for (i = 0; i < 8; i++) {
        const uint32_t x = pg->regs[NV_PGRAPH_WINDOWCLIPX0 + i * 4];
        const uint32_t y = pg->regs[NV_PGRAPH_WINDOWCLIPY0 + i * 4];
        const uint32_t x_min = GET_MASK(x, NV_PGRAPH_WINDOWCLIPX0_XMIN);
        const uint32_t x_max = GET_MASK(x, NV_PGRAPH_WINDOWCLIPX0_XMAX);
        const uint32_t y_min = GET_MASK(y, NV_PGRAPH_WINDOWCLIPY0_YMIN);
        const uint32_t y_max = GET_MASK(y, NV_PGRAPH_WINDOWCLIPY0_YMAX);

        /* Check for zero width or height clipping region */
        if ((x_min == x_max) || (y_min == y_max)) {
            continue;
        }

        /* Check for in-order duplicate regions */
        if ((x == last_x) && (y == last_y)) {
            continue;
        }

        NV2A_DPRINTF("Clipping Region %d: min=(%d, %d) max=(%d, %d)\n",
            i, x_min, y_min, x_max, y_max);

        state.psh.window_clip_count = i + 1;
        last_x = x;
        last_y = y;
    }

    for (i = 0; i < 8; i++) {
        state.psh.rgb_inputs[i] = pg->regs[NV_PGRAPH_COMBINECOLORI0 + i * 4];
        state.psh.rgb_outputs[i] = pg->regs[NV_PGRAPH_COMBINECOLORO0 + i * 4];
        state.psh.alpha_inputs[i] = pg->regs[NV_PGRAPH_COMBINEALPHAI0 + i * 4];
        state.psh.alpha_outputs[i] = pg->regs[NV_PGRAPH_COMBINEALPHAO0 + i * 4];
        //constant_0[i] = pg->regs[NV_PGRAPH_COMBINEFACTOR0 + i * 4];
        //constant_1[i] = pg->regs[NV_PGRAPH_COMBINEFACTOR1 + i * 4];
    }

    for (i = 0; i < 4; i++) {
        state.psh.rect_tex[i] = false;
        bool enabled = pg->regs[NV_PGRAPH_TEXCTL0_0 + i*4]
                         & NV_PGRAPH_TEXCTL0_0_ENABLE;
        unsigned int color_format =
            GET_MASK(pg->regs[NV_PGRAPH_TEXFMT0 + i*4],
                     NV_PGRAPH_TEXFMT0_COLOR);

        if (enabled && kelvin_color_format_map[color_format].encoding == linear) {
            state.psh.rect_tex[i] = true;
        }

        for (j = 0; j < 4; j++) {
            state.psh.compare_mode[i][j] =
                (pg->regs[NV_PGRAPH_SHADERCLIPMODE] >> (4 * i + j)) & 1;
        }
        state.psh.alphakill[i] = pg->regs[NV_PGRAPH_TEXCTL0_0 + i*4]
                               & NV_PGRAPH_TEXCTL0_0_ALPHAKILLEN;
    }

#ifdef USE_SHADER_CACHE
	ShaderBinding* cached_shader = (ShaderBinding*)g_hash_table_lookup(pg->shader_cache, &state);

    if (cached_shader) {
        pg->shader_binding = cached_shader;
    } else {
#endif
        pg->shader_binding = generate_shaders(state);

#ifdef USE_SHADER_CACHE
        /* cache it */
        ShaderState *cache_state = (ShaderState *)g_malloc(sizeof(*cache_state));
        memcpy(cache_state, &state, sizeof(*cache_state));
        g_hash_table_insert(pg->shader_cache, cache_state,
                            (gpointer)pg->shader_binding);
    }
#endif

    bool binding_changed = (pg->shader_binding != old_binding);

    glUseProgram(pg->shader_binding->gl_program);

    /* Clipping regions */
    for (i = 0; i < state.psh.window_clip_count; i++) {
        if (pg->shader_binding->clip_region_loc[i] == -1) {
            continue;
        }

        uint32_t x   = pg->regs[NV_PGRAPH_WINDOWCLIPX0 + i * 4];
        GLuint x_min = GET_MASK(x, NV_PGRAPH_WINDOWCLIPX0_XMIN);
        GLuint x_max = GET_MASK(x, NV_PGRAPH_WINDOWCLIPX0_XMAX);

        /* Adjust y-coordinates for the OpenGL viewport: translate coordinates
         * to have the origin at the bottom-left of the surface (as opposed to
         * top-left), and flip y-min and y-max accordingly.
         */
        uint32_t y   = pg->regs[NV_PGRAPH_WINDOWCLIPY0 + i * 4];
        GLuint y_min = (pg->surface_shape.clip_height - 1) -
                       GET_MASK(y, NV_PGRAPH_WINDOWCLIPY0_YMAX);
        GLuint y_max = (pg->surface_shape.clip_height - 1) -
                       GET_MASK(y, NV_PGRAPH_WINDOWCLIPY0_YMIN);

        pgraph_apply_anti_aliasing_factor(pg, &x_min, &y_min);
        pgraph_apply_anti_aliasing_factor(pg, &x_max, &y_max);

        glUniform4i(pg->shader_binding->clip_region_loc[i],
                    x_min, y_min, x_max, y_max);
    }

    OpenGL_pgraph_update_shader_constants(pg, pg->shader_binding, binding_changed,
                                   vertex_program, fixed_function);

    NV2A_GL_DGROUP_END();
}

static void OpenGL_pgraph_update_surface_part(NV2AState *d, bool upload, bool color)
{
    PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

    unsigned int width, height;
    pgraph_get_surface_dimensions(pg, &width, &height);
    pgraph_apply_anti_aliasing_factor(pg, &width, &height);

    Surface *surface;
    hwaddr dma_address;
    GLuint *gl_buffer;
    unsigned int bytes_per_pixel;
    GLint gl_internal_format;
    GLenum gl_format, gl_type, gl_attachment;

    if (color) {
		typedef struct SurfaceColorFormatInfo {
			unsigned int bytes_per_pixel;
			GLint gl_internal_format;
			GLenum gl_format;
			GLenum gl_type;
		} SurfaceColorFormatInfo;

		// Note : Avoid designated initializers to facilitate C++ builds
		static const SurfaceColorFormatInfo kelvin_surface_color_format_map[16] = {
			//0x00 [?] = 
				{},
			//0x01 [NV097_SET_SURFACE_FORMAT_COLOR_LE_X1R5G5B5_Z1R5G5B5] =
				{2, GL_RGB5_A1, GL_BGRA, GL_UNSIGNED_SHORT_1_5_5_5_REV},
			//0x02 [NV097_SET_SURFACE_FORMAT_COLOR_LE_X1R5G5B5_O1R5G5B5] =
				{},
			//0x03 [NV097_SET_SURFACE_FORMAT_COLOR_LE_R5G6B5] =
				{2, GL_RGB565, GL_RGB, GL_UNSIGNED_SHORT_5_6_5},
			//0x04 [NV097_SET_SURFACE_FORMAT_COLOR_LE_X8R8G8B8_Z8R8G8B8] =
				{4, GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV},
			//0x05 [NV097_SET_SURFACE_FORMAT_COLOR_LE_X8R8G8B8_O8R8G8B8] =
				{},
			//0x06 [NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_Z1A7R8G8B8] =
				{},
			//0x07 [NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_O1A7R8G8B8] =
				{},
			//0x08 [NV097_SET_SURFACE_FORMAT_COLOR_LE_A8R8G8B8] =
				{4, GL_RGBA8, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV},
			//0x09 [NV097_SET_SURFACE_FORMAT_COLOR_LE_B8] =
				{}, // TODO : {1, GL_R8, GL_RED, GL_UNSIGNED_BYTE}, // PatrickvL guesstimate
			//0x0A [NV097_SET_SURFACE_FORMAT_COLOR_LE_G8B8] =
				{}, // TODO : {2, GL_RG8, GL_RG, GL_UNSIGNED_BYTE}, // PatrickvL guesstimate
			//0x0B [?] = 
				{},
			//0x0C [?] = 
				{},
			//0x0D [?] = 
				{},
			//0x0E [?] = 
				{},
			//0x0F [?] = 
				{}
		};

        surface = &pg->surface_color;
        dma_address = pg->dma_color;
        gl_buffer = &pg->gl_color_buffer;

        assert(pg->surface_shape.color_format != 0);
        assert(pg->surface_shape.color_format
                < ARRAY_SIZE(kelvin_surface_color_format_map));
        SurfaceColorFormatInfo f =
            kelvin_surface_color_format_map[pg->surface_shape.color_format];
        if (f.bytes_per_pixel == 0) {
            fprintf(stderr, "nv2a: unimplemented color surface format 0x%x\n",
                    pg->surface_shape.color_format);
            abort();
        }

        bytes_per_pixel = f.bytes_per_pixel;
        gl_internal_format = f.gl_internal_format;
        gl_format = f.gl_format;
        gl_type = f.gl_type;
        gl_attachment = GL_COLOR_ATTACHMENT0;

    } else {
        surface = &pg->surface_zeta;
        dma_address = pg->dma_zeta;
        gl_buffer = &pg->gl_zeta_buffer;

        assert(pg->surface_shape.zeta_format != 0);
        switch (pg->surface_shape.zeta_format) {
        case NV097_SET_SURFACE_FORMAT_ZETA_Z16:
            bytes_per_pixel = 2;
            gl_format = GL_DEPTH_COMPONENT;
            gl_attachment = GL_DEPTH_ATTACHMENT;
            if (pg->surface_shape.z_format) {
                gl_type = GL_HALF_FLOAT;
                gl_internal_format = GL_DEPTH_COMPONENT32F;
            } else {
                gl_type = GL_UNSIGNED_SHORT;
                gl_internal_format = GL_DEPTH_COMPONENT16;
            }
            break;
        case NV097_SET_SURFACE_FORMAT_ZETA_Z24S8:
            bytes_per_pixel = 4;
            gl_format = GL_DEPTH_STENCIL;
            gl_attachment = GL_DEPTH_STENCIL_ATTACHMENT;
            if (pg->surface_shape.z_format) {
                assert(false);
                gl_type = GL_FLOAT_32_UNSIGNED_INT_24_8_REV;
                gl_internal_format = GL_DEPTH32F_STENCIL8;
            } else {
                gl_type = GL_UNSIGNED_INT_24_8;
                gl_internal_format = GL_DEPTH24_STENCIL8;
            }
            break;
        default:
            assert(false);
            break;
        }
    }


    DMAObject dma = nv_dma_load(d, dma_address);
    /* There's a bunch of bugs that could cause us to hit this function
     * at the wrong time and get a invalid dma object.
     * Check that it's sane. */
    assert(dma.dma_class == NV_DMA_IN_MEMORY_CLASS);

    assert(dma.address + surface->offset != 0);
    assert(surface->offset <= dma.limit);
    assert(surface->offset + surface->pitch * height <= dma.limit + 1);

    hwaddr data_len;
    uint8_t *data = (uint8_t*)nv_dma_map(d, dma_address, &data_len);

    /* TODO */
    // assert(pg->surface_clip_x == 0 && pg->surface_clip_y == 0);

    bool swizzle = (pg->surface_type == NV097_SET_SURFACE_FORMAT_TYPE_SWIZZLE);

    uint8_t *buf = data + surface->offset;
    if (swizzle) {
        buf = (uint8_t*)g_malloc(height * surface->pitch);
    }

    bool dirty = surface->buffer_dirty;
    if (color) {
#if 1
		// HACK: Always mark as dirty
        dirty |= 1;
#else
        dirty |= memory_region_test_and_clear_dirty(d->vram,
                                               dma.address + surface->offset,
                                               surface->pitch * height,
                                               DIRTY_MEMORY_NV2A);
#endif
    }
    if (upload && dirty) {
        /* surface modified (or moved) by the cpu.
         * copy it into the opengl renderbuffer */
        // TODO: Why does this assert?
		//assert(!surface->draw_dirty);
		assert(surface->pitch % bytes_per_pixel == 0);

        if (swizzle) {
            unswizzle_rect(data + surface->offset,
                           width, height,
                           buf,
                           surface->pitch,
                           bytes_per_pixel);
        }

		if (!color) {
			/* need to clear the depth_stencil and depth attachment for zeta */
			glFramebufferTexture2D(GL_FRAMEBUFFER,
									GL_DEPTH_ATTACHMENT,
									GL_TEXTURE_2D,
									0, 0);
			glFramebufferTexture2D(GL_FRAMEBUFFER,
									GL_DEPTH_STENCIL_ATTACHMENT,
									GL_TEXTURE_2D,
									0, 0);
		}

		glFramebufferTexture2D(GL_FRAMEBUFFER,
								gl_attachment,
								GL_TEXTURE_2D,
								0, 0);

		if (*gl_buffer) {
			glDeleteTextures(1, gl_buffer);
			*gl_buffer = 0;
		}

		glGenTextures(1, gl_buffer);
		glBindTexture(GL_TEXTURE_2D, *gl_buffer);

		/* This is VRAM so we can't do this inplace! */
		uint8_t *flipped_buf = (uint8_t*)g_malloc(width * height * bytes_per_pixel);
		unsigned int irow;
		for (irow = 0; irow < height; irow++) {
			memcpy(&flipped_buf[width * (height - irow - 1)
										* bytes_per_pixel],
					&buf[surface->pitch * irow],
					width * bytes_per_pixel);
		}

		glTexImage2D(GL_TEXTURE_2D, 0, gl_internal_format,
						width, height, 0,
						gl_format, gl_type,
						flipped_buf);

		g_free(flipped_buf);

		glFramebufferTexture2D(GL_FRAMEBUFFER,
								gl_attachment,
								GL_TEXTURE_2D,
								*gl_buffer, 0);

		assert(glCheckFramebufferStatus(GL_FRAMEBUFFER)
			== GL_FRAMEBUFFER_COMPLETE);

        if (color) {
            OpenGL_pgraph_update_memory_buffer(d, dma.address + surface->offset,
                                        surface->pitch * height, true);
        }
        surface->buffer_dirty = false;

#ifdef DEBUG_NV2A
        uint8_t *out = data + surface->offset + 64;
        NV2A_DPRINTF("upload_surface %s 0x%" HWADDR_PRIx " - 0x%" HWADDR_PRIx ", "
                      "(0x%" HWADDR_PRIx " - 0x%" HWADDR_PRIx ", "
                        "%d %d, %d %d, %d) - %x %x %x %x\n",
            color ? "color" : "zeta",
            dma.address, dma.address + dma.limit,
            dma.address + surface->offset,
            dma.address + surface->pitch * height,
            pg->surface_shape.clip_x, pg->surface_shape.clip_y,
            pg->surface_shape.clip_width,
            pg->surface_shape.clip_height,
            surface->pitch,
            out[0], out[1], out[2], out[3]);
#endif
    }

    if (!upload && surface->draw_dirty) {
		/* read the opengl framebuffer into the surface */

		glo_readpixels(gl_format, gl_type,
						bytes_per_pixel, surface->pitch,
						width, height,
						buf);

//		assert(glGetError() == GL_NO_ERROR);

        if (swizzle) {
            swizzle_rect(buf,
                         width, height,
                         data + surface->offset,
                         surface->pitch,
                         bytes_per_pixel);
        }

        // memory_region_set_client_dirty(d->vram,
        //                                dma.address + surface->offset,
        //                                surface->pitch * height,
        //                                DIRTY_MEMORY_VGA);

        if (color) {
            OpenGL_pgraph_update_memory_buffer(d, dma.address + surface->offset,
                                        surface->pitch * height, true);
        }

        surface->draw_dirty = false;
        surface->write_enabled_cache = false;

#ifdef DEBUG_NV2A
        uint8_t *out = data + surface->offset + 64;
        NV2A_DPRINTF("read_surface %s 0x%" HWADDR_PRIx " - 0x%" HWADDR_PRIx ", "
                      "(0x%" HWADDR_PRIx " - 0x%" HWADDR_PRIx ", "
                        "%d %d, %d %d, %d) - %x %x %x %x\n",
            color ? "color" : "zeta",
            dma.address, dma.address + dma.limit,
            dma.address + surface->offset,
            dma.address + surface->pitch * pg->surface_shape.clip_height,
            pg->surface_shape.clip_x, pg->surface_shape.clip_y,
            pg->surface_shape.clip_width, pg->surface_shape.clip_height,
            surface->pitch,
            out[0], out[1], out[2], out[3]);
#endif
    }

    if (swizzle) {
        g_free(buf);
    }
}

static void OpenGL_pgraph_update_surface(NV2AState *d, bool upload,
	bool color_write, bool zeta_write)
{
    PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

    pg->surface_shape.z_format = GET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
                                          NV_PGRAPH_SETUPRASTER_Z_FORMAT);

    /* FIXME: Does this apply to CLEARs too? */
    color_write = color_write && pgraph_get_color_write_enabled(pg);
    zeta_write = zeta_write && pgraph_get_zeta_write_enabled(pg);

    if (upload && pgraph_get_framebuffer_dirty(pg)) {
        assert(!pg->surface_color.draw_dirty);
        assert(!pg->surface_zeta.draw_dirty);

        pg->surface_color.buffer_dirty = true;
        pg->surface_zeta.buffer_dirty = true;

		glFramebufferTexture2D(GL_FRAMEBUFFER,
								GL_COLOR_ATTACHMENT0,
								GL_TEXTURE_2D,
								0, 0);

		if (pg->gl_color_buffer) {
			glDeleteTextures(1, &pg->gl_color_buffer);
			pg->gl_color_buffer = 0;
		}

		glFramebufferTexture2D(GL_FRAMEBUFFER,
								GL_DEPTH_ATTACHMENT,
								GL_TEXTURE_2D,
								0, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER,
								GL_DEPTH_STENCIL_ATTACHMENT,
								GL_TEXTURE_2D,
								0, 0);

		if (pg->gl_zeta_buffer) {
			glDeleteTextures(1, &pg->gl_zeta_buffer);
			pg->gl_zeta_buffer = 0;
		}
	}

    memcpy(&pg->last_surface_shape, &pg->surface_shape,
            sizeof(SurfaceShape));

    if ((color_write || (!upload && pg->surface_color.write_enabled_cache))
        && (upload || pg->surface_color.draw_dirty)) {
		OpenGL_pgraph_update_surface_part(d, upload, true);
    }


    if ((zeta_write || (!upload && pg->surface_zeta.write_enabled_cache))
        && (upload || pg->surface_zeta.draw_dirty)) {
		OpenGL_pgraph_update_surface_part(d, upload, false);
    }

}

static void OpenGL_pgraph_bind_textures(NV2AState *d)
{
    int i;
    PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

    NV2A_GL_DGROUP_BEGIN("%s", __func__);

    for (i=0; i<NV2A_MAX_TEXTURES; i++) {

        uint32_t ctl_0 = pg->regs[NV_PGRAPH_TEXCTL0_0 + i*4];
        uint32_t ctl_1 = pg->regs[NV_PGRAPH_TEXCTL1_0 + i*4];
        uint32_t fmt = pg->regs[NV_PGRAPH_TEXFMT0 + i*4];
        uint32_t filter = pg->regs[NV_PGRAPH_TEXFILTER0 + i*4];
        uint32_t address =  pg->regs[NV_PGRAPH_TEXADDRESS0 + i*4];
        uint32_t palette =  pg->regs[NV_PGRAPH_TEXPALETTE0 + i*4];

        bool enabled = ctl_0 & NV_PGRAPH_TEXCTL0_0_ENABLE;
        unsigned int min_mipmap_level =
            GET_MASK(ctl_0, NV_PGRAPH_TEXCTL0_0_MIN_LOD_CLAMP);
        unsigned int max_mipmap_level =
            GET_MASK(ctl_0, NV_PGRAPH_TEXCTL0_0_MAX_LOD_CLAMP);

        unsigned int pitch =
            GET_MASK(ctl_1, NV_PGRAPH_TEXCTL1_0_IMAGE_PITCH);

        bool dma_select =
            fmt & NV_PGRAPH_TEXFMT0_CONTEXT_DMA;
        bool cubemap =
            fmt & NV_PGRAPH_TEXFMT0_CUBEMAPENABLE;
        unsigned int dimensionality =
            GET_MASK(fmt, NV_PGRAPH_TEXFMT0_DIMENSIONALITY);
        unsigned int color_format = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_COLOR);
        unsigned int levels = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_MIPMAP_LEVELS);
        unsigned int log_width = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_BASE_SIZE_U);
        unsigned int log_height = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_BASE_SIZE_V);
        unsigned int log_depth = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_BASE_SIZE_P);

        unsigned int rect_width =
            GET_MASK(pg->regs[NV_PGRAPH_TEXIMAGERECT0 + i*4],
                     NV_PGRAPH_TEXIMAGERECT0_WIDTH);
        unsigned int rect_height =
            GET_MASK(pg->regs[NV_PGRAPH_TEXIMAGERECT0 + i*4],
                     NV_PGRAPH_TEXIMAGERECT0_HEIGHT);
#ifdef DEBUG_NV2A
        unsigned int lod_bias =
            GET_MASK(filter, NV_PGRAPH_TEXFILTER0_MIPMAP_LOD_BIAS);
#endif
        unsigned int min_filter = GET_MASK(filter, NV_PGRAPH_TEXFILTER0_MIN);
        unsigned int mag_filter = GET_MASK(filter, NV_PGRAPH_TEXFILTER0_MAG);

        unsigned int addru = GET_MASK(address, NV_PGRAPH_TEXADDRESS0_ADDRU);
        unsigned int addrv = GET_MASK(address, NV_PGRAPH_TEXADDRESS0_ADDRV);
        unsigned int addrp = GET_MASK(address, NV_PGRAPH_TEXADDRESS0_ADDRP);

		bool border_source_color = (fmt & NV_PGRAPH_TEXFMT0_BORDER_SOURCE); // != NV_PGRAPH_TEXFMT0_BORDER_SOURCE_TEXTURE;
        uint32_t border_color = pg->regs[NV_PGRAPH_BORDERCOLOR0 + i*4];

        unsigned int offset = pg->regs[NV_PGRAPH_TEXOFFSET0 + i*4];

        bool palette_dma_select =
            palette & NV_PGRAPH_TEXPALETTE0_CONTEXT_DMA;
        unsigned int palette_length_index =
            GET_MASK(palette, NV_PGRAPH_TEXPALETTE0_LENGTH);
        unsigned int palette_offset =
            palette & NV_PGRAPH_TEXPALETTE0_OFFSET;

        unsigned int palette_length = 0;
        switch (palette_length_index) {
        case NV_PGRAPH_TEXPALETTE0_LENGTH_256: palette_length = 256; break;
        case NV_PGRAPH_TEXPALETTE0_LENGTH_128: palette_length = 128; break;
        case NV_PGRAPH_TEXPALETTE0_LENGTH_64: palette_length = 64; break;
        case NV_PGRAPH_TEXPALETTE0_LENGTH_32: palette_length = 32; break;
        default: assert(false); break;
        }

        /* Check for unsupported features */
        assert(!(filter & NV_PGRAPH_TEXFILTER0_ASIGNED));
        assert(!(filter & NV_PGRAPH_TEXFILTER0_RSIGNED));
        assert(!(filter & NV_PGRAPH_TEXFILTER0_GSIGNED));
        assert(!(filter & NV_PGRAPH_TEXFILTER0_BSIGNED));

        glActiveTexture(GL_TEXTURE0 + i);
        if (!enabled) {
            glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
            glBindTexture(GL_TEXTURE_RECTANGLE, 0);
            glBindTexture(GL_TEXTURE_1D, 0);
            glBindTexture(GL_TEXTURE_2D, 0);
            glBindTexture(GL_TEXTURE_3D, 0);
            continue;
        }

        if (!pg->texture_dirty[i] && pg->texture_binding[i]) {
            glBindTexture(pg->texture_binding[i]->gl_target,
                          pg->texture_binding[i]->gl_texture);
            continue;
        }

        NV2A_DPRINTF(" texture %d is format 0x%x, off 0x%x (r %d, %d or %d, %d, %d; %d%s),"
                        " filter %x %x, levels %d-%d %d bias %d\n",
                     i, color_format, offset,
                     rect_width, rect_height,
                     1 << log_width, 1 << log_height, 1 << log_depth,
                     pitch,
                     cubemap ? "; cubemap" : "",
                     min_filter, mag_filter,
                     min_mipmap_level, max_mipmap_level, levels,
                     lod_bias);

        assert(color_format < ARRAY_SIZE(kelvin_color_format_map));
        ColorFormatInfo f = kelvin_color_format_map[color_format];
        if (f.bytes_per_pixel == 0) {
            fprintf(stderr, "nv2a: unimplemented texture color format 0x%x\n",
                    color_format);
            abort();
        }

        unsigned int width, height, depth;
        if (f.encoding == linear) {
            assert(dimensionality == 2);
            width = rect_width;
            height = rect_height;
            depth = 1;
        } else {
            width = 1 << log_width;
            height = 1 << log_height;
            depth = 1 << log_depth;

            /* FIXME: What about 3D mipmaps? */
            levels = MIN(levels, max_mipmap_level + 1);
            if (f.encoding == swizzled) {
                /* Discard mipmap levels that would be smaller than 1x1.
                 * FIXME: Is this actually needed?
                 *
                 * >> Level 0: 32 x 4
                 *    Level 1: 16 x 2
                 *    Level 2: 8 x 1
                 *    Level 3: 4 x 1
                 *    Level 4: 2 x 1
                 *    Level 5: 1 x 1
                 */
                levels = MIN(levels, MAX(log_width, log_height) + 1);
            } else {
                /* OpenGL requires DXT textures to always have a width and
                 * height a multiple of 4. The Xbox and DirectX handles DXT
                 * textures smaller than 4 by padding the reset of the block.
                 *
                 * See:
                 * https://msdn.microsoft.com/en-us/library/windows/desktop/bb204843(v=vs.85).aspx
                 * https://msdn.microsoft.com/en-us/library/windows/desktop/bb694531%28v=vs.85%29.aspx#Virtual_Size
                 *
                 * Work around this for now by discarding mipmap levels that
                 * would result in too-small textures. A correct solution
                 * will be to decompress these levels manually, or add texture
                 * sampling logic.
                 *
                 * >> Level 0: 64 x 8
                 *    Level 1: 32 x 4
                 *    Level 2: 16 x 2 << Ignored
                 * >> Level 0: 16 x 16
                 *    Level 1: 8 x 8
                 *    Level 2: 4 x 4 << OK!
                 */
                if (log_width < 2 || log_height < 2) {
                    /* Base level is smaller than 4x4... */
                    levels = 1;
                } else {
                    levels = MIN(levels, MIN(log_width, log_height) - 1);
                }
            }
            assert(levels > 0);
        }

        hwaddr dma_len;
        uint8_t *texture_data;
        if (dma_select) {
            texture_data = (uint8_t*)nv_dma_map(d, pg->dma_b, &dma_len);
        } else {
            texture_data = (uint8_t*)nv_dma_map(d, pg->dma_a, &dma_len);
        }
        assert(offset < dma_len);
        texture_data += offset;

        hwaddr palette_dma_len;
        uint8_t *palette_data;
        if (palette_dma_select) {
            palette_data = (uint8_t*)nv_dma_map(d, pg->dma_b, &palette_dma_len);
        } else {
            palette_data = (uint8_t*)nv_dma_map(d, pg->dma_a, &palette_dma_len);
        }
        assert(palette_offset < palette_dma_len);
        palette_data += palette_offset;

        NV2A_DPRINTF(" - 0x%tx\n", texture_data - d->vram_ptr);

        size_t length = 0;
        if (f.encoding == linear) {
            assert(cubemap == false);
            assert(dimensionality == 2);
            length = height * pitch;
        } else {
            if (dimensionality >= 2) {
                unsigned int w = width, h = height;
                unsigned int level;
                if (f.encoding == swizzled) {
                    for (level = 0; level < levels; level++) {
                        w = MAX(w, 1); h = MAX(h, 1);
                        length += w * h * f.bytes_per_pixel;
                        w /= 2;
                        h /= 2;
                    }
                } else {
                    /* Compressed textures are a bit different */
                    unsigned int block_size;
                    if (f.gl_internal_format ==
                            GL_COMPRESSED_RGBA_S3TC_DXT1_EXT) {
                        block_size = 8;
                    } else {
                        block_size = 16;
                    }

                    for (level = 0; level < levels; level++) {
                        w = MAX(w, 4); h = MAX(h, 4);
                        length += w/4 * h/4 * block_size;
                        w /= 2; h /= 2;
                    }
                }
                if (cubemap) {
                    assert(dimensionality == 2);
                    length *= 6;
                }
                if (dimensionality >= 3) {
                    length *= depth;
                }
            }
        }

		TextureShape state;
		state.cubemap = cubemap;
		state.dimensionality = dimensionality;
		state.color_format = color_format;
		state.levels = levels;
		state.width = width;
		state.height = height;
		state.depth = depth;
		state.min_mipmap_level = min_mipmap_level;
		state.max_mipmap_level = max_mipmap_level;
        state.pitch = pitch;

#ifdef USE_TEXTURE_CACHE
		TextureKey key;
		key.state = state;
		key.data_hash = fast_hash(texture_data, length, 5003)
			^ fnv_hash(palette_data, palette_length);
		key.texture_data = texture_data;
		key.palette_data = palette_data;

        gpointer cache_key = g_malloc(sizeof(TextureKey));
        memcpy(cache_key, &key, sizeof(TextureKey));

        GError *err;
        TextureBinding *binding = (TextureBinding *)g_lru_cache_get(pg->texture_cache, cache_key, &err);
        assert(binding);
        binding->refcnt++;
#else
        TextureBinding *binding = OpenGL_generate_texture(state,
                                                   texture_data, palette_data);
#endif

        glBindTexture(binding->gl_target, binding->gl_texture);


        if (f.encoding == linear) {
            /* sometimes games try to set mipmap min filters on linear textures.
             * this could indicate a bug... */
            switch (min_filter) {
            case NV_PGRAPH_TEXFILTER0_MIN_BOX_NEARESTLOD:
            case NV_PGRAPH_TEXFILTER0_MIN_BOX_TENT_LOD:
                min_filter = NV_PGRAPH_TEXFILTER0_MIN_BOX_LOD0;
                break;
            case NV_PGRAPH_TEXFILTER0_MIN_TENT_NEARESTLOD:
            case NV_PGRAPH_TEXFILTER0_MIN_TENT_TENT_LOD:
                min_filter = NV_PGRAPH_TEXFILTER0_MIN_TENT_LOD0;
                break;
            }
        }

		{
			static const GLenum pgraph_texture_min_filter_map[] = {
				0,
				GL_NEAREST,
				GL_LINEAR,
				GL_NEAREST_MIPMAP_NEAREST,
				GL_LINEAR_MIPMAP_NEAREST,
				GL_NEAREST_MIPMAP_LINEAR,
				GL_LINEAR_MIPMAP_LINEAR,
				GL_LINEAR, /* TODO: Convolution filter... */
			};

			assert(min_filter < ARRAY_SIZE(pgraph_texture_min_filter_map));
			glTexParameteri(binding->gl_target, GL_TEXTURE_MIN_FILTER,
				pgraph_texture_min_filter_map[min_filter]);
		}

		{
			static const GLenum pgraph_texture_mag_filter_map[] = {
				0,
				GL_NEAREST,
				GL_LINEAR,
				0,
				GL_LINEAR /* TODO: Convolution filter... */
			};

			assert(mag_filter < ARRAY_SIZE(pgraph_texture_mag_filter_map));
			glTexParameteri(binding->gl_target, GL_TEXTURE_MAG_FILTER,
				pgraph_texture_mag_filter_map[mag_filter]);
		}
        /* Texture wrapping */
		{
			static const GLenum pgraph_texture_addr_map[] = {
				0,
				GL_REPEAT,
				GL_MIRRORED_REPEAT,
				GL_CLAMP_TO_EDGE,
				GL_CLAMP_TO_BORDER,
				// GL_CLAMP
			};

			assert(addru < ARRAY_SIZE(pgraph_texture_addr_map));
			glTexParameteri(binding->gl_target, GL_TEXTURE_WRAP_S,
				pgraph_texture_addr_map[addru]);
			if (dimensionality > 1) {
				assert(addrv < ARRAY_SIZE(pgraph_texture_addr_map));
				glTexParameteri(binding->gl_target, GL_TEXTURE_WRAP_T,
					pgraph_texture_addr_map[addrv]);
			}
			if (dimensionality > 2) {
				assert(addrp < ARRAY_SIZE(pgraph_texture_addr_map));
				glTexParameteri(binding->gl_target, GL_TEXTURE_WRAP_R,
					pgraph_texture_addr_map[addrp]);
			}
		}

        /* FIXME: Only upload if necessary? [s, t or r = GL_CLAMP_TO_BORDER] */
        if (border_source_color) {
            GLfloat gl_border_color[] = {
                /* FIXME: Color channels might be wrong order */
                ((border_color >> 16) & 0xFF) / 255.0f, /* red */
                ((border_color >> 8) & 0xFF) / 255.0f,  /* green */
                (border_color & 0xFF) / 255.0f,         /* blue */
                ((border_color >> 24) & 0xFF) / 255.0f  /* alpha */
            };
            glTexParameterfv(binding->gl_target, GL_TEXTURE_BORDER_COLOR,
                gl_border_color);
        }

        if (pg->texture_binding[i]) {
			OpenGL_texture_binding_destroy(pg->texture_binding[i]);
        }
        pg->texture_binding[i] = binding;
        pg->texture_dirty[i] = false;
    }
    NV2A_GL_DGROUP_END();
}

static void OpenGL_pgraph_update_memory_buffer(NV2AState *d, hwaddr addr, hwaddr size,
                                        bool f)
{
	glBindBuffer(GL_ARRAY_BUFFER, d->pgraph.gl_memory_buffer);

	hwaddr end = TARGET_PAGE_ALIGN(addr + size);
	addr &= TARGET_PAGE_MASK;

	assert(end < d->vram_size);

    // if (f || memory_region_test_and_clear_dirty(d->vram,
    //                                             addr,
    //                                             end - addr,
    //                                             DIRTY_MEMORY_NV2A)) {
		glBufferSubData(GL_ARRAY_BUFFER, addr, end - addr, d->vram_ptr + addr);
    // }

//		auto error = glGetError();
//		assert(error == GL_NO_ERROR);
}

static void OpenGL_pgraph_bind_vertex_attributes(NV2AState *d,
                                          unsigned int num_elements,
                                          bool inline_data,
                                          unsigned int inline_stride)
{
	unsigned int i, j;
    PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

    if (inline_data) {
        NV2A_GL_DGROUP_BEGIN("%s (num_elements: %d inline stride: %d)",
                             __func__, num_elements, inline_stride);
    } else {
        NV2A_GL_DGROUP_BEGIN("%s (num_elements: %d)", __func__, num_elements);
    }

    for (i=0; i<NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
        VertexAttribute *vertex_attribute = &pg->vertex_attributes[i];
        if (vertex_attribute->count) {
            uint8_t *data;
            unsigned int in_stride;
            if (inline_data && vertex_attribute->needs_conversion) {
                data = (uint8_t*)pg->inline_array
                        + vertex_attribute->inline_array_offset;
                in_stride = inline_stride;
            } else {
                hwaddr dma_len;
                if (vertex_attribute->dma_select) {
                    data = (uint8_t*)nv_dma_map(d, pg->dma_vertex_b, &dma_len);
                } else {
                    data = (uint8_t*)nv_dma_map(d, pg->dma_vertex_a, &dma_len);
                }

                assert(vertex_attribute->offset < dma_len);
                data += vertex_attribute->offset;

                in_stride = vertex_attribute->stride;
            }

            if (vertex_attribute->needs_conversion) {
                NV2A_DPRINTF("converted %d\n", i);

                unsigned int out_stride = vertex_attribute->converted_size
                                        * vertex_attribute->converted_count;

                if (num_elements > vertex_attribute->converted_elements) {
                    vertex_attribute->converted_buffer = (uint8_t*)g_realloc(
                        vertex_attribute->converted_buffer,
                        num_elements * out_stride);
                }

                for (j=vertex_attribute->converted_elements; j<num_elements; j++) {
                    uint8_t *in = data + j * in_stride;
                    uint8_t *out = vertex_attribute->converted_buffer + j * out_stride;

                    switch (vertex_attribute->format) {
                    case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_CMP: {
                        uint32_t p = ldl_le_p((uint32_t*)in);
                        float *xyz = (float*)out;
                        xyz[0] = ((int32_t)(((p >>  0) & 0x7FF) << 21) >> 21)
                                                                      / 1023.0f;
                        xyz[1] = ((int32_t)(((p >> 11) & 0x7FF) << 21) >> 21)
                                                                      / 1023.0f;
                        xyz[2] = ((int32_t)(((p >> 22) & 0x3FF) << 22) >> 22)
                                                                       / 511.0f;
                        break;
                    }
                    default:
                        assert(false);
                        break;
                    }
                }


                glBindBuffer(GL_ARRAY_BUFFER, vertex_attribute->gl_converted_buffer);
                if (num_elements != vertex_attribute->converted_elements) {
                    glBufferData(GL_ARRAY_BUFFER,
                                 num_elements * out_stride,
                                 vertex_attribute->converted_buffer,
                                 GL_DYNAMIC_DRAW);
                    vertex_attribute->converted_elements = num_elements;
                }


                glVertexAttribPointer(i,
                    vertex_attribute->converted_count,
                    vertex_attribute->gl_type,
                    vertex_attribute->gl_normalize,
                    out_stride,
                    0);
            } else if (inline_data) {
                glBindBuffer(GL_ARRAY_BUFFER, pg->gl_inline_array_buffer);
                glVertexAttribPointer(i,
                                      vertex_attribute->gl_count,
                                      vertex_attribute->gl_type,
                                      vertex_attribute->gl_normalize,
                                      inline_stride,
                                      (void*)(uintptr_t)vertex_attribute->inline_array_offset);
            } else {
                hwaddr addr = data - d->vram_ptr;
                OpenGL_pgraph_update_memory_buffer(d, addr,
                                            num_elements * vertex_attribute->stride,
                                            false);
                glVertexAttribPointer(i,
                    vertex_attribute->gl_count,
                    vertex_attribute->gl_type,
                    vertex_attribute->gl_normalize,
                    vertex_attribute->stride,
                    (void*)(uint64_t)(addr));
            }
            glEnableVertexAttribArray(i);
        } else {
            glDisableVertexAttribArray(i);

            glVertexAttrib4fv(i, vertex_attribute->inline_value);
        }
    }
    NV2A_GL_DGROUP_END();
}

static unsigned int OpenGL_pgraph_bind_inline_array(NV2AState *d)
{
	int i;

    PGRAPHState *pg = &d->pgraph;

	assert(pg->opengl_enabled);

    unsigned int offset = 0;
    for (i=0; i<NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
        VertexAttribute *vertex_attribute = &pg->vertex_attributes[i];
        if (vertex_attribute->count) {
            vertex_attribute->inline_array_offset = offset;

            NV2A_DPRINTF("bind inline vertex_attribute %d size=%d, count=%d\n",
                i, vertex_attribute->size, vertex_attribute->count);
            offset += vertex_attribute->size * vertex_attribute->count;
            assert(offset % 4 == 0);
        }
    }

    unsigned int vertex_size = offset;


    unsigned int index_count = pg->inline_array_length*4 / vertex_size;

    NV2A_DPRINTF("draw inline array %d, %d\n", vertex_size, index_count);

    glBindBuffer(GL_ARRAY_BUFFER, pg->gl_inline_array_buffer);
    glBufferData(GL_ARRAY_BUFFER, pg->inline_array_length*4, pg->inline_array,
                 GL_DYNAMIC_DRAW);

	OpenGL_pgraph_bind_vertex_attributes(d, index_count, true, vertex_size);

    return index_count;
}

/* returns the format of the output, either identical to the input format, or the converted format - see converted_format */
static int OpenGL_upload_texture(GLenum gl_target,
                              const TextureShape s,
                              const uint8_t *texture_data,
                              const uint8_t *palette_data)
{
	// assert(pg->opengl_enabled); // Lacks in-scope pg
    int resulting_format = s.color_format;
    ColorFormatInfo f = kelvin_color_format_map[s.color_format];

    switch(gl_target) {
    case GL_TEXTURE_1D:
        assert(false);
        break;
    case GL_TEXTURE_RECTANGLE: {
        /* Can't handle strides unaligned to pixels */
        assert(s.pitch % f.bytes_per_pixel == 0);
        glPixelStorei(GL_UNPACK_ROW_LENGTH,
                      s.pitch / f.bytes_per_pixel);

		uint8_t *unswizzled = NULL;
        if (f.encoding == swizzled) { // TODO : Verify this works correctly
            unswizzled = (uint8_t*)g_malloc(s.height * s.pitch);
            unswizzle_rect(texture_data, s.width, s.height,
                            unswizzled, s.pitch, f.bytes_per_pixel);
        }
        uint8_t *converted = convert_texture_data(s.color_format, unswizzled ? unswizzled : texture_data,
                                                  palette_data,
                                                  s.width, s.height, 1,
                                                  s.pitch, 0);

        resulting_format = converted ? converted_format : s.color_format;
        ColorFormatInfo cf = kelvin_color_format_map[resulting_format];
        glTexImage2D(gl_target, 0, cf.gl_internal_format,
                     s.width, s.height, 0,
                     cf.gl_format, cf.gl_type,
                     converted ? converted : unswizzled ? unswizzled : texture_data);

        if (converted) {
          g_free(converted);
        }
		if (unswizzled) {
			g_free(unswizzled);
        }

        glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
        break;
    }
    case GL_TEXTURE_2D:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_X:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_X:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_Y:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_Y:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_Z:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_Z: {

        unsigned int width = s.width, height = s.height;

        unsigned int level;
        for (level = 0; level < s.levels; level++) {
            if (f.encoding == compressed) {

                width = MAX(width, 4); height = MAX(height, 4);

                unsigned int block_size;
                if (f.gl_internal_format == GL_COMPRESSED_RGBA_S3TC_DXT1_EXT) {
                    block_size = 8;
                } else {
                    block_size = 16;
                }

                glCompressedTexImage2D(gl_target, level, f.gl_internal_format,
                                       width, height, 0,
                                       width/4 * height/4 * block_size,
                                       texture_data);

                texture_data += width/4 * height/4 * block_size;
            } else {

                width = MAX(width, 1); height = MAX(height, 1);

                unsigned int pitch = width * f.bytes_per_pixel;
                uint8_t *unswizzled = NULL;
                if (f.encoding == swizzled) {
                    unswizzled = (uint8_t*)g_malloc(height * pitch);
                    unswizzle_rect(texture_data, width, height,
                                   unswizzled, pitch, f.bytes_per_pixel);
                }

                uint8_t *converted = convert_texture_data(s.color_format, unswizzled ? unswizzled : texture_data,
                                                          palette_data,
                                                          width, height, 1,
                                                          pitch, 0);

                resulting_format = converted ? converted_format : s.color_format;
                ColorFormatInfo cf = kelvin_color_format_map[resulting_format];
                glTexImage2D(gl_target, level, cf.gl_internal_format,
                             width, height, 0,
                             cf.gl_format, cf.gl_type,
                             converted ? converted : unswizzled ? unswizzled : texture_data);

                if (converted) {
                    g_free(converted);
                }
                if (unswizzled) {
                    g_free(unswizzled);
                }

                texture_data += pitch * height;
            }

            width /= 2;
            height /= 2;
        }

        break;
    }
    case GL_TEXTURE_3D: {

        unsigned int width = s.width, height = s.height, depth = s.depth;

        unsigned int level;
        for (level = 0; level < s.levels; level++) {

            unsigned int row_pitch = width * f.bytes_per_pixel;
            unsigned int slice_pitch = row_pitch * height;
            uint8_t *unswizzled = NULL;
            if (f.encoding == swizzled) {
                unswizzled = (uint8_t*)g_malloc(slice_pitch * depth);
                unswizzle_box(texture_data, width, height, depth, unswizzled,
                               row_pitch, slice_pitch, f.bytes_per_pixel);
            }
            uint8_t *converted = convert_texture_data(s.color_format, unswizzled ? unswizzled : texture_data,
                                                      palette_data,
                                                      width, height, depth,
                                                      row_pitch, slice_pitch);

            resulting_format = converted ? converted_format : s.color_format;
            ColorFormatInfo cf = kelvin_color_format_map[resulting_format];
            glTexImage3D(gl_target, level, cf.gl_internal_format,
                         width, height, depth, 0,
                         cf.gl_format, cf.gl_type,
                         converted ? converted : unswizzled ? unswizzled : texture_data);

            if (converted) {
                g_free(converted);
            }
            if (unswizzled) {
                g_free(unswizzled);
            }

            texture_data += width * height * depth * f.bytes_per_pixel;

            width /= 2;
            height /= 2;
            depth /= 2;
        }
        break;
    }
    default:
        assert(false);
        break;
    }
	return resulting_format;
}

static TextureBinding* OpenGL_generate_texture(const TextureShape s,
                                        const uint8_t *texture_data,
                                        const uint8_t *palette_data)
{
	// assert(pg->opengl_enabled); // Lacks in-scope pg

    ColorFormatInfo f = kelvin_color_format_map[s.color_format];

    /* Create a new opengl texture */
    GLuint gl_texture;
    glGenTextures(1, &gl_texture);

    GLenum gl_target;
    if (s.cubemap) {
        assert(f.encoding != linear);
        assert(s.dimensionality == 2);
        gl_target = GL_TEXTURE_CUBE_MAP;
    } else {
        if (f.encoding == linear) { /* FIXME : Include compressed too? (!= swizzled) */
            /* linear textures use unnormalised texcoords.
             * GL_TEXTURE_RECTANGLE_ARB conveniently also does, but
             * does not allow repeat and mirror wrap modes.
             *  (or mipmapping, but xbox d3d says 'Non swizzled and non
             *   compressed textures cannot be mip mapped.')
             * Not sure if that'll be an issue. */

            /* FIXME: GLSL 330 provides us with textureSize()! Use that? */
            gl_target = GL_TEXTURE_RECTANGLE;
            assert(s.dimensionality == 2);
        } else {
            switch(s.dimensionality) {
            case 1: gl_target = GL_TEXTURE_1D; break;
            case 2: gl_target = GL_TEXTURE_2D; break;
            case 3: gl_target = GL_TEXTURE_3D; break;
            default:
                assert(false);
                break;
            }
        }
    }

    glBindTexture(gl_target, gl_texture);

    NV2A_GL_DLABEL(GL_TEXTURE, gl_texture,
                   "format: 0x%02X%s, %d dimensions%s, width: %d, height: %d, depth: %d",
                   s.color_format, (f.encoding == linear) ? "" : (f.encoding == swizzled) ? " (SZ)" : " (DXT)", // compressed
                   s.dimensionality, s.cubemap ? " (Cubemap)" : "",
                   s.width, s.height, s.depth);

    /* Linear textures don't support mipmapping */
    if (f.encoding != linear) {
        glTexParameteri(gl_target, GL_TEXTURE_BASE_LEVEL,
            s.min_mipmap_level);
        glTexParameteri(gl_target, GL_TEXTURE_MAX_LEVEL,
            s.levels - 1);
    }

	/* Set this before calling OpenGL_upload_texture() to prevent potential conversions */
    if (f.gl_swizzle_mask) {
        glTexParameteriv(gl_target, GL_TEXTURE_SWIZZLE_RGBA,
                         f.gl_swizzle_mask);
    }

    if (gl_target == GL_TEXTURE_CUBE_MAP) {

        size_t length = 0;
        unsigned int w = s.width, h = s.height;
        unsigned int level;
        for (level = 0; level < s.levels; level++) {
            /* FIXME: This is wrong for compressed textures and textures with 1x? non-square mipmaps */
            length += w * h * f.bytes_per_pixel;
            w /= 2;
            h /= 2;
        }

        OpenGL_upload_texture(GL_TEXTURE_CUBE_MAP_POSITIVE_X,
                          s, texture_data + 0 * length, palette_data);
        OpenGL_upload_texture(GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
                          s, texture_data + 1 * length, palette_data);
        OpenGL_upload_texture(GL_TEXTURE_CUBE_MAP_POSITIVE_Y,
                          s, texture_data + 2 * length, palette_data);
        OpenGL_upload_texture(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
                          s, texture_data + 3 * length, palette_data);
        OpenGL_upload_texture(GL_TEXTURE_CUBE_MAP_POSITIVE_Z,
                          s, texture_data + 4 * length, palette_data);
        OpenGL_upload_texture(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,
                          s, texture_data + 5 * length, palette_data);
    } else {
        OpenGL_upload_texture(gl_target, s, texture_data, palette_data);
    }

    TextureBinding* ret = (TextureBinding *)g_malloc(sizeof(TextureBinding));
    ret->gl_target = gl_target;
    ret->gl_texture = gl_texture;
    ret->refcnt = 1;
    return ret;
}

static gpointer OpenGL_texture_key_retrieve(gpointer key, gpointer user_data, GError **error)
{
    const TextureKey *k = (const TextureKey *)key;
    TextureBinding *v = OpenGL_generate_texture(k->state,
                                         k->texture_data,
                                         k->palette_data);
    if (error != NULL) {
        *error = NULL;
    }
    return v;
}

static void OpenGL_texture_binding_destroy(gpointer data)
{
    TextureBinding *binding = (TextureBinding *)data;

	// assert(pg->opengl_enabled); // Lacks in-scope pg

    assert(binding->refcnt > 0);
    binding->refcnt--;
    if (binding->refcnt == 0) {
        glDeleteTextures(1, &binding->gl_texture);
        g_free(binding);
    }
}
