# gen_cso_header.cmake — Convert a compiled shader binary (.cso) to a C byte array header.
# Usage: cmake -DCSO_FILE=path/to/shader.cso -DHEADER_FILE=path/to/output.h
#              -DARRAY_NAME=shader_name -P gen_cso_header.cmake

file(READ "${CSO_FILE}" CSO_DATA HEX)
string(LENGTH "${CSO_DATA}" CSO_HEX_LEN)
math(EXPR CSO_BYTE_COUNT "${CSO_HEX_LEN} / 2")

# Build unsigned char array from hex pairs
set(ARRAY_BODY "")
set(COL 0)
math(EXPR LAST_BYTE "${CSO_BYTE_COUNT} - 1")
foreach(I RANGE 0 ${LAST_BYTE})
    math(EXPR POS "${I} * 2")
    string(SUBSTRING "${CSO_DATA}" ${POS} 2 BYTE_HEX)
    string(APPEND ARRAY_BODY "0x${BYTE_HEX},")
    math(EXPR COL "${COL} + 1")
    if(COL EQUAL 16)
        string(APPEND ARRAY_BODY "\n    ")
        set(COL 0)
    endif()
endforeach()

file(WRITE "${HEADER_FILE}"
"// Auto-generated from ${CSO_FILE} — do not edit.\n"
"#pragma once\n"
"#include <cstddef>\n"
"static const unsigned char ${ARRAY_NAME}[] = {\n"
"    ${ARRAY_BODY}\n"
"};\n"
"static const size_t ${ARRAY_NAME}_size = ${CSO_BYTE_COUNT};\n"
)
