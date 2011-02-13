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

#include <QtCore>

namespace chemkit {

class Plugin;
class PluginManagerPrivate;

class CHEMKIT_EXPORT PluginManager : public QObject
{
    Q_OBJECT

    public:
        // properties
        Plugin* plugin(const QString &name);
        const Plugin* plugin(const QString &name) const;
        QList<Plugin *> plugins();
        QList<const Plugin *> plugins() const;
        int pluginCount() const;

        // plugin loading
        bool loadPlugin(const QString &fileName);
        void loadPlugins(const QString &directory);
        void loadDefaultPlugins();
        bool unloadPlugin(Plugin *plugin);
        bool unloadPlugin(const QString &name);

        // error handling
        QString errorString() const;

        // static methods
        static PluginManager* instance();

    signals:
        void pluginLoaded(const chemkit::Plugin *plugin);
        void pluginUnloaded(const chemkit::Plugin *plugin);

    private:
        PluginManager();
        ~PluginManager();

        void setErrorString(const QString &errorString);

        Q_DISABLE_COPY(PluginManager);

    private:
        PluginManagerPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PLUGINMANAGER_H
