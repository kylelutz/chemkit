/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
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
    typename T::CreateFunction function = reinterpret_cast<typename T::CreateFunction>(pluginClassFunction(typeid(T).name(), pluginName));

    if(function){
        return function();
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
