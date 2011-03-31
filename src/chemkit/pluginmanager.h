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

#ifndef CHEMKIT_PLUGINMANAGER_H
#define CHEMKIT_PLUGINMANAGER_H

#include "chemkit.h"

#include <string>

#include <QtCore>

namespace chemkit {

class Plugin;
class PluginManagerPrivate;

class CHEMKIT_EXPORT PluginManager : public QObject
{
    Q_OBJECT

    public:
        // enumerations
        typedef void (*Function)();

        // properties
        Plugin* plugin(const std::string &name) const;
        QList<Plugin *> plugins() const;
        int pluginCount() const;

        // plugin loading
        bool loadPlugin(const std::string &fileName);
        void loadPlugins(const std::string &directory);
        void loadDefaultPlugins();
        bool unloadPlugin(Plugin *plugin);
        bool unloadPlugin(const std::string &name);

        // plugin classes
        template<class T> T* createPluginClass(const std::string &pluginName) const;
        template<class T> std::vector<std::string> pluginClassNames() const;

        // error handling
        std::string errorString() const;

        // static methods
        static PluginManager* instance();

    Q_SIGNALS:
        void pluginLoaded(const chemkit::Plugin *plugin);
        void pluginUnloaded(const chemkit::Plugin *plugin);

    private:
        PluginManager();
        ~PluginManager();

        void setErrorString(const std::string &errorString);
        bool registerPluginClass(const std::string &className, const std::string &pluginName, Function function);
        bool unregisterPluginClass(const std::string &className, const std::string &pluginName);
        std::vector<std::string> pluginClassNames(const std::string &className) const;
        Function pluginClassFunction(const std::string &className, const std::string &pluginName) const;

        Q_DISABLE_COPY(PluginManager);

        friend class Plugin;

    private:
        PluginManagerPrivate* const d;
};

} // end chemkit namespace

#include "pluginmanager-inline.h"

#endif // CHEMKIT_PLUGINMANAGER_H
