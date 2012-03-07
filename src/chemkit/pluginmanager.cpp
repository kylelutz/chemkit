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
#include <cstdlib>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "plugin.h"
#include "foreach.h"
#include "dynamiclibrary.h"

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
    : d(new PluginManagerPrivate)
{
    d->defaultPluginsLoaded = false;
}

PluginManager::~PluginManager()
{
    d->pluginClasses.clear();

    foreach(Plugin *plugin, d->plugins){
        DynamicLibrary *library = plugin->library();
        delete plugin;
        delete library;
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
    DynamicLibrary *library = new DynamicLibrary;
    library->setFileName(fileName);
    bool ok = library->open();
    if(!ok){
        std::cerr << "PluginManager: Error: Failed to load plugin: (" << fileName << "):"
                  << library->errorString() << std::endl;

        return false;
    }

    typedef Plugin* (*InitFunction)();
    InitFunction initFunction = reinterpret_cast<InitFunction>(library->resolveFunction("chemkit_plugin_init"));
    if(!initFunction){
        std::cerr << "PluginManager: Error: Failed to load plugin: (" << fileName << "):"
                  << "Plugin contains no init() function." << std::endl;

        return false;
    }

    Plugin *plugin = initFunction();
    if(!plugin){
        std::cerr << "PluginManager: Error: Failed to load plugin: (" << fileName << "):"
                  << "Calling the plugin's init() function failed." << std::endl;

        return false;
    }

    plugin->setLibrary(library);

    d->plugins.push_back(plugin);

    return true;
}

/// Loads all plugins from \p directory.
void PluginManager::loadPlugins(const std::string &directory)
{
    boost::filesystem::path dir(directory);

    if(!boost::filesystem::exists(dir)){
        return;
    }

    for(boost::filesystem::directory_iterator iter(dir); iter != boost::filesystem::directory_iterator(); ++iter){
        std::string fileName = boost::filesystem::path(iter->path().filename()).string();

        if(DynamicLibrary::isLibrary(fileName)){
            loadPlugin(iter->path().string());
        }
    }
}

void PluginManager::loadDefaultPlugins()
{
    if(d->defaultPluginsLoaded){
        return;
    }

    // list of directories to load plugins from
    std::vector<std::string> directories;

    // add default plugin directory
#if defined(CHEMKIT_OS_UNIX)
    directories.push_back(CHEMKIT_INSTALL_PREFIX "/lib/chemkit/plugins/");
#endif

    // add directory from the CHEMKIT_PLUGIN_PATH environment variable
    const char *path = getenv("CHEMKIT_PLUGIN_PATH");
    if(path){
        directories.push_back(path);
    }

    // load plugins from each directory
    foreach(const std::string &directory, directories){
        loadPlugins(directory);
    }

    d->defaultPluginsLoaded = true;
}

/// Unloads the plugin.
bool PluginManager::unloadPlugin(Plugin *plugin)
{
    if(!plugin){
        return false;
    }

    d->plugins.erase(std::remove(d->plugins.begin(), d->plugins.end(), plugin));
    DynamicLibrary *library = plugin->library();
    delete plugin;
    delete library;

    return true;
}

/// Unloads the plugin with \p name.
bool PluginManager::unloadPlugin(const std::string &name)
{
    return unloadPlugin(plugin((name)));
}

// --- Error Handling ------------------------------------------------------ //
void PluginManager::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occurred.
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

// --- Internal Methods ---------------------------------------------------- //
/// Registers a new plugin function for \p className and
/// \p pluginName.
bool PluginManager::registerPluginClass(const std::string &className, const std::string &pluginName, boost::function<void* ()> function)
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
boost::function<void* ()> PluginManager::pluginClassFunction(const std::string &className, const std::string &pluginName) const
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
