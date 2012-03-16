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

#ifndef CHEMKIT_CHEMKIT_H
#define CHEMKIT_CHEMKIT_H

#include <cstddef>

#include <boost/preprocessor/stringize.hpp>

#include "config.h"

#if defined(CHEMKIT_OS_UNIX)
    #define CHEMKIT_DECL_IMPORT
    #define CHEMKIT_DECL_EXPORT  __attribute__((visibility("default")))
#elif defined(CHEMKIT_OS_WIN32)
    #define CHEMKIT_DECL_IMPORT  __declspec(dllimport)
    #define CHEMKIT_DECL_EXPORT  __declspec(dllexport)
#endif

#ifdef CHEMKIT_LIBRARY
    #define CHEMKIT_EXPORT CHEMKIT_DECL_EXPORT
#else
    #define CHEMKIT_EXPORT CHEMKIT_DECL_IMPORT
#endif

/// A string containing the version number of the chemkit
/// library (e.g. "1.2").
#define CHEMKIT_VERSION_STRING BOOST_PP_STRINGIZE(CHEMKIT_VERSION_MAJOR) "." \
                               BOOST_PP_STRINGIZE(CHEMKIT_VERSION_MINOR)

#define CHEMKIT_UNUSED(variable) (void) variable

/// This macro marks a class as not copyable. It should be used in
/// the private section of a class's declaration.
#define CHEMKIT_DISABLE_COPY(Class) \
    Class(const Class &); \
    Class &operator=(const Class &);

namespace chemkit {

/// Typedef for a real number.
typedef double Real;

} // end chemkit namespace

#endif // CHEMKIT_CHEMKIT_H
