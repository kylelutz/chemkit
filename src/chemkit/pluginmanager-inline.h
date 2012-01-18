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

#ifndef CHEMKIT_PLUGINMANAGER_INLINE_H
#define CHEMKIT_PLUGINMANAGER_INLINE_H

#include "pluginmanager.h"

#include <typeinfo>

namespace chemkit {

// --- Plugin Classes ------------------------------------------------------ //
/// Creates and returns a new instance of a plugin class \p T from
/// \p pluginName. Returns \c 0 if \p pluginName is not found.
///
/// The ownership of the returned object is passed to the caller.
template<class T>
inline T* PluginManager::createPluginClass(const std::string &pluginName) const
{
    Function function = pluginClassFunction(typeid(T).name(), pluginName);

    if(function){
        return static_cast<T*>(function());
    }

    return 0;
}

/// Returns a vector of the names of the plugins registered for the
/// class \p T.
template<class T>
inline std::vector<std::string> PluginManager::pluginClassNames() const
{
    return pluginClassNames(typeid(T).name());
}

} // end chemkit namespace

#endif // CHEMKIT_PLUGINMANAGER_INLINE_H
