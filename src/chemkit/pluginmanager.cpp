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

#include <map>
#include <boost/algorithm/string/case_conv.hpp>

#include "plugin.h"

namespace chemkit {

// === PluginManagerPrivate ================================================ //
class PluginManagerPrivate
{
    public:
        QList<Plugin *> plugins;
        std::string errorString;
        bool defaultPluginsLoaded;
        std::map<std::string, std::map<std::string, PluginManager::Function> > pluginClasses;
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
    Q_FOREACH(Plugin *plugin, d->plugins){
        delete plugin;
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the plugin with \p name. Returns \c 0 if no plugin with
// \p name is loaded.
Plugin* PluginManager::plugin(const std::string &name) const
{
    Q_FOREACH(Plugin *plugin, d->plugins){
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

    instance->setFileName(fileName.toStdString());

    d->plugins.append(instance);

    Q_EMIT pluginLoaded(instance);

    return true;
}

/// Loads all plugins from \p directory.
void PluginManager::loadPlugins(const QString &directory)
{
    QDir dir(directory);

    if(!dir.exists()){
        return;
    }

    Q_FOREACH(const QString &fileName, dir.entryList(QDir::Files)){
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
    Q_FOREACH(const QString &directory, directories){
        loadPlugins(directory);
    }

    d->defaultPluginsLoaded = true;
}

/// Unloads the plugin.
bool PluginManager::unloadPlugin(Plugin *plugin)
{
    Q_EMIT pluginUnloaded(plugin);

    return false;
}

/// Unloads the plugin with \p name.
bool PluginManager::unloadPlugin(const QString &name)
{
    Q_UNUSED(name);

    return false;
}

// --- Error Handling ------------------------------------------------------ //
void PluginManager::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
std::string PluginManager::errorString() const
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

// --- Internal Methods ---------------------------------------------------- //
/// Registers a new plugin function for \p className and
/// \p pluginName.
bool PluginManager::registerPluginClass(const std::string &className, const std::string &pluginName, Function function)
{
    std::map<std::string, Function> &classPlugins = d->pluginClasses[className];

    // use lower case plugin name
    std::string lowerCasePluginName = boost::algorithm::to_lower_copy(pluginName);

    // prevent overwriting of previously registered plugins
    if(classPlugins.find(lowerCasePluginName) != classPlugins.end()){
        return false;
    }

    // add plugin class
    classPlugins[lowerCasePluginName] = function;

    return true;
}

/// Unregisters a plugin function for \p className and \p pluginName.
bool PluginManager::unregisterPluginClass(const std::string &className, const std::string &pluginName)
{
    std::map<std::string, Function> &classPlugins = d->pluginClasses[className];

    // use lower case plugin name
    std::string lowerCasePluginName = boost::algorithm::to_lower_copy(pluginName);

    // remove plugin class
    return classPlugins.erase(lowerCasePluginName) > 0;
}

/// Returns a vector of strings containing the names of registered
/// plugins for \p className.
std::vector<std::string> PluginManager::pluginClassNames(const std::string &className) const
{
    // ensure default plugins are loaded
    const_cast<PluginManager *>(this)->loadDefaultPlugins();

    const std::map<std::string, Function> &classPlugins = d->pluginClasses[className];

    std::vector<std::string> names;
    for(std::map<std::string, Function>::const_iterator i = classPlugins.begin(); i != classPlugins.end(); ++i){
        names.push_back(i->first);
    }

    return names;
}

/// Returns the registered function for the given \p className and \p pluginName.
PluginManager::Function PluginManager::pluginClassFunction(const std::string &className, const std::string &pluginName) const
{
    // ensure default plugins are loaded
    const_cast<PluginManager *>(this)->loadDefaultPlugins();

    // use lower case plugin name
    std::string lowerCasePluginName = boost::algorithm::to_lower_copy(pluginName);

    const std::map<std::string, Function> &classPlugins = d->pluginClasses[className];

    std::map<std::string, Function>::const_iterator location = classPlugins.find(lowerCasePluginName);
    if(location == classPlugins.end()){
        return 0;
    }

    return location->second;
}

} // end chemkit namespace
