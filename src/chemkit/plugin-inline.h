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

#ifndef CHEMKIT_PLUGIN_INLINE_H
#define CHEMKIT_PLUGIN_INLINE_H

#include "plugin.h"

#include <typeinfo>

#include "pluginmanager.h"

namespace chemkit {

/// Registers a new plugin class with \p name and \p function.
template<class T>
inline bool Plugin::registerPluginClass(const std::string &name, typename T::CreateFunction function)
{
    return PluginManager::instance()->registerPluginClass(typeid(T).name(), name, reinterpret_cast<PluginManager::Function>(function));
}

/// Unregisters the plugin class with \p name.
template<class T>
inline bool Plugin::unregisterPluginClass(const std::string &name)
{
    return PluginManager::instance()->unregisterPluginClass(typeid(T).name(), name);
}

} // end chemkit namespace

#endif // CHEMKIT_PLUGIN_INLINE_H
