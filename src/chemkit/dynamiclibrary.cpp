/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "dynamiclibrary.h"

#include <boost/algorithm/string.hpp>

#if defined(CHEMKIT_OS_UNIX)
#include <dlfcn.h>
#elif defined(CHEMKIT_OS_WIN32)
#include <windows.h>
#endif

// Setup typedef for LibraryHandle
#if defined(CHEMKIT_OS_UNIX)
typedef void* LibraryHandle;
#elif defined(CHEMKIT_OS_WIN32)
typedef HMODULE LibraryHandle;
#else
typedef void* LibraryHandle;
#endif

namespace chemkit {

// === DynamicLibraryPrivate =============================================== //
class DynamicLibraryPrivate
{
public:
    std::string fileName;
    std::string errorString;
    LibraryHandle handle;
};

// === DynamicLibrary ====================================================== //
/// \class DynamicLibrary dynamiclibrary.h chemkit/dynamiclibrary.h
/// \ingroup chemkit
/// \internal
/// \brief The DynamicLibrary class represents a dynamically loaded
///        library.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new dynamic library object.
DynamicLibrary::DynamicLibrary()
    : d(new DynamicLibraryPrivate)
{
    d->handle = 0;
}

/// Creates a new dynamic library object with \p fileName.
DynamicLibrary::DynamicLibrary(const std::string &fileName)
    : d(new DynamicLibraryPrivate)
{
    d->handle = 0;
    d->fileName = fileName;
}

/// Destroys the dynamic library object. This will call close() if
/// the library is currently open.
DynamicLibrary::~DynamicLibrary()
{
    close();

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the file name for the library to \p fileName.
void DynamicLibrary::setFileName(const std::string &fileName)
{
    d->fileName = fileName;
}

/// Returns the file name of the library.
std::string DynamicLibrary::fileName() const
{
    return d->fileName;
}

// --- Loading and Unloading ----------------------------------------------- //
/// Opens the library. Returns \c false if an error occurs.
bool DynamicLibrary::open()
{
    // close library if currently open
    if(d->handle){
        close();
    }

    // open library
#if defined(CHEMKIT_OS_UNIX)
    LibraryHandle handle = dlopen(d->fileName.c_str(), RTLD_LAZY);
    if(!handle){
        setErrorString(dlerror());
        return false;
    }

    d->handle = handle;
    return true;

#elif defined(CHEMKIT_OS_WIN32)
    LibraryHandle handle = LoadLibrary(d->fileName.c_str());
    if(!handle){
        return false;
    }

    d->handle = handle;
    return true;

#else
    return false;
#endif
}

/// Unloads and closes the library.
void DynamicLibrary::close()
{
    if(!d->handle){
        return;
    }

#if defined(CHEMKIT_OS_UNIX)
    dlclose(d->handle);
#elif defined(CHEMKIT_OS_WIN32)
    FreeLibrary(d->handle);
#endif

    d->handle = 0;
}

// --- Symbol Resolution --------------------------------------------------- //
/// Resolves \p symbol and returns its address. Returns \c 0 if the
/// symbol does not exist.
void* DynamicLibrary::resolve(const std::string &symbol)
{
    // open library if not already opened
    if(!d->handle){
        bool ok = open();

        if(!ok){
            return 0;
        }
    }

    // lookup and return the symbol's address
#if defined(CHEMKIT_OS_UNIX)
    return dlsym(d->handle, symbol.c_str());

#elif defined(CHEMKIT_OS_WIN32)
    return reinterpret_cast<void *>(GetProcAddress(d->handle, symbol.c_str()));

#else
    return 0;
#endif
}

/// Resolves \p symbol and returns its address as a function
/// pointer. Returns \c 0 if the symbol does not exist.
DynamicLibrary::Function DynamicLibrary::resolveFunction(const std::string &symbol)
{
    void *address = resolve(symbol);
    if(!address){
        return 0;
    }

    union FunctionPointer {
        void *address;
        Function function;
    };

    FunctionPointer pointer;
    pointer.address = address;

    return pointer.function;
}

// --- Error Handling ------------------------------------------------------ //
void DynamicLibrary::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occurred.
std::string DynamicLibrary::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns \c true if \p fileName ends with the native operating
/// systems dynamic library suffix (e.g. ".so" on *nix and ".dll"
/// on Windows).
bool DynamicLibrary::isLibrary(const std::string &fileName)
{
#if defined(CHEMKIT_OS_UNIX)
 #if defined(CHEMKIT_OS_MAC)
    return boost::algorithm::ends_with(fileName, ".dylib");
 #else
    return boost::algorithm::ends_with(fileName, ".so");
 #endif
#elif defined(CHEMKIT_OS_WIN32)
    return boost::algorithm::ends_with(fileName, ".dll");
#else
    return false;
#endif
}

} // end chemkit namespace
