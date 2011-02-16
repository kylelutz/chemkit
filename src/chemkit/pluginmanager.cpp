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

#include "pluginmanager.h"

#include "plugin.h"

namespace chemkit {

// === PluginManagerPrivate ================================================ //
class PluginManagerPrivate
{
    public:
        QList<Plugin *> plugins;
        QString errorString;
        bool defaultPluginsLoaded;
};

// === PluginManager ======================================================= //
/// \class PluginManager pluginmanager.h chemkit/pluginmanager.h
/// \ingroup chemkit
/// \brief The PluginManager class manages the loading and unloading
///        of plugins.
///
/// \see Plugin

// --- Construction and Destruction ---------------------------------------- //
PluginManager::PluginManager()
    : QObject(),
      d(new PluginManagerPrivate)
{
    d->defaultPluginsLoaded = false;
}

PluginManager::~PluginManager()
{
    foreach(Plugin *plugin, d->plugins){
        delete plugin;
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the plugin with \p name. Returns \c 0 if no plugin with
// \p name is loaded.
Plugin* PluginManager::plugin(const QString &name) const
{
    foreach(Plugin *plugin, d->plugins){
        if(plugin->name() == name){
            return plugin;
        }
    }

    return 0;
}

/// Returns a list of all the loaded plugins.
QList<Plugin *> PluginManager::plugins() const
{
    return d->plugins;
}

/// Returns the number of loaded plugins.
int PluginManager::pluginCount() const
{
    return d->plugins.size();
}

// --- Plugin Loading ------------------------------------------------------ //
/// Loads a plugin from \p fileName. Returns \c false if an error
/// occurs.
bool PluginManager::loadPlugin(const QString &fileName)
{
    QPluginLoader plugin(fileName);

    Plugin *instance = qobject_cast<Plugin *>(plugin.instance());
    if(!instance){
        qDebug() << "Failed to load plugin (" << fileName << "): " << plugin.errorString();
        return false;
    }

    instance->setFileName(fileName);

    d->plugins.append(instance);

    emit pluginLoaded(instance);

    return true;
}

/// Loads all plugins from \p directory.
void PluginManager::loadPlugins(const QString &directory)
{
    QDir dir(directory);

    if(!dir.exists()){
        return;
    }

    foreach(const QString &fileName, dir.entryList(QDir::Files)){
        if(QLibrary::isLibrary(fileName)){
            loadPlugin(dir.filePath(fileName));
        }
    }
}

void PluginManager::loadDefaultPlugins()
{
    if(d->defaultPluginsLoaded){
        return;
    }

    // list of directories to load plugins from
    QStringList directories;

    // add default plugin directory
#if defined(Q_OS_LINUX)
    directories.append("/usr/share/chemkit/plugins/");
#elif defined(Q_OS_WIN32)
    QSettings registry("HKEY_LOCAL_MACHINE\\Software\\chemkit", QSettings::NativeFormat);
    directories.append(registry.value("PluginPath").toString());
#endif

    // add directory from the CHEMKIT_PLUGIN_PATH environment variable
    QProcessEnvironment environment = QProcessEnvironment::systemEnvironment();
    QString path = environment.value("CHEMKIT_PLUGIN_PATH");
    if(!path.isEmpty()){
        directories.append(path);
    }

    // load plugins from each directory
    foreach(const QString &directory, directories){
        loadPlugins(directory);
    }

    d->defaultPluginsLoaded = true;
}

/// Unloads the plugin.
bool PluginManager::unloadPlugin(Plugin *plugin)
{
    emit pluginUnloaded(plugin);

    return false;
}

/// Unloads the plugin with \p name.
bool PluginManager::unloadPlugin(const QString &name)
{
    Q_UNUSED(name);

    return false;
}

// --- Error Handling ------------------------------------------------------ //
void PluginManager::setErrorString(const QString &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
QString PluginManager::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the instance of the plugin manager.
PluginManager* PluginManager::instance()
{
    static PluginManager singleton;

    return &singleton;
}

// --- Signals ------------------------------------------------------------- //
/// \fn void PluginManager::pluginLoaded(const chemkit::Plugin *plugin)
///
/// This signal is emitted when a new plugin is loaded.

/// \fn void PluginManager::pluginUnloaded(const chemkit::Plugin *plugin)
///
/// This signal is emitted when a plugin is unloaded.

} // end chemkit namespace
