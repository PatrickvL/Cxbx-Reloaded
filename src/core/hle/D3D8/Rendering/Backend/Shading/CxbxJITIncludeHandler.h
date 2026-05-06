// CxbxJITIncludeHandler.h — ID3DInclude that resolves #include directives
// by loading files relative to the emulator executable's "hlsl\" subdirectory.
// Used by both the VS and PS JIT compilers.

#pragma once

#include <d3dcompiler.h>
#include <string>
#include <windows.h>

class CxbxJITIncludeHandler : public ID3DInclude
{
public:
    CxbxJITIncludeHandler()
    {
        char exePath[MAX_PATH] = {};
        GetModuleFileNameA(nullptr, exePath, MAX_PATH);
        char* lastSlash = strrchr(exePath, '\\');
        if (lastSlash) *(lastSlash + 1) = '\0';
        m_shaderDir = std::string(exePath) + "hlsl\\";
    }

    STDMETHOD(Open)(D3D_INCLUDE_TYPE /*IncludeType*/, LPCSTR pFileName,
        LPCVOID /*pParentData*/, LPCVOID* ppData, UINT* pBytes) override
    {
        std::string fullPath = m_shaderDir + pFileName;
        HANDLE hFile = CreateFileA(fullPath.c_str(), GENERIC_READ, FILE_SHARE_READ,
            nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (hFile == INVALID_HANDLE_VALUE)
            return E_FAIL;

        DWORD fileSize = GetFileSize(hFile, nullptr);
        if (fileSize == INVALID_FILE_SIZE) { CloseHandle(hFile); return E_FAIL; }

        char* pData = new (std::nothrow) char[fileSize];
        if (!pData) { CloseHandle(hFile); return E_OUTOFMEMORY; }

        DWORD bytesRead = 0;
        BOOL ok = ReadFile(hFile, pData, fileSize, &bytesRead, nullptr);
        CloseHandle(hFile);
        if (!ok || bytesRead != fileSize) { delete[] pData; return E_FAIL; }

        *ppData = pData;
        *pBytes = bytesRead;
        return S_OK;
    }

    STDMETHOD(Close)(LPCVOID pData) override
    {
        delete[] static_cast<const char*>(pData);
        return S_OK;
    }

private:
    std::string m_shaderDir;
};
