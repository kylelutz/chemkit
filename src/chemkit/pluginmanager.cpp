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

#include "pluginmanager.h"

#include <map>
#include <boost/algorithm/string/case_conv.hpp>

#include "plugin.h"
#include "foreach.h"

namespace chemkit {

// === PluginManagerPrivate ================================================ //
class PluginManagerPrivate
{
    public:
        std::vector<Plugin *> plugins;
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
    foreach(Plugin *plugin, d->plugins){
        if(plugin->name() == name){
            return plugin;
        }
    }

    return 0;
}

/// Returns a list of all the loaded plugins.
const std::vector<Plugin *>& PluginManager::plugins() const
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
bool PluginManager::loadPlugin(const std::string &fileName)
{
    QPluginLoader plugin(QString::fromStdString(fileName));

    Plugin *instance = qobject_cast<Plugin *>(plugin.instance());
    if(!instance){
        qDebug() << "Failed to load plugin (" << fileName.c_str() << "): " << plugin.errorString();
        return false;
    }

    instance->setFileName(fileName);

    d->plugins.push_back(instance);

    Q_EMIT pluginLoaded(instance);

    return true;
}

/// Loads all plugins from \p directory.
void PluginManager::loadPlugins(const std::string &directory)
{
    QDir dir(QString::fromStdString(directory));

    if(!dir.exists()){
        return;
    }

    Q_FOREACH(const QString &fileName, dir.entryList(QDir::Files)){
        if(QLibrary::isLibrary(fileName)){
            loadPlugin(dir.filePath(fileName).toStdString());
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
    directories.append(CHEMKIT_INSTALL_PREFIX "/share/chemkit/plugins/");
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
        loadPlugins(directory.toStdString());
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
bool PluginManager::unloadPlugin(const std::string &name)
{
    CHEMKIT_UNUSED(name);

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
