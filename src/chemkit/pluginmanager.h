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

#ifndef CHEMKIT_PLUGINMANAGER_H
#define CHEMKIT_PLUGINMANAGER_H

#include "chemkit.h"

#include <string>
#include <vector>

#include <boost/function.hpp>

namespace chemkit {

class Plugin;
class PluginManagerPrivate;

class CHEMKIT_EXPORT PluginManager
{
public:
    // enumerations
    typedef boost::function<void* ()> Function;

    // properties
    Plugin* plugin(const std::string &name) const;
    const std::vector<Plugin *>& plugins() const;
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

private:
    PluginManager();
    ~PluginManager();

    void setErrorString(const std::string &errorString);
    bool registerPluginClass(const std::string &className, const std::string &pluginName, Function function);
    bool unregisterPluginClass(const std::string &className, const std::string &pluginName);
    std::vector<std::string> pluginClassNames(const std::string &className) const;
    Function pluginClassFunction(const std::string &className, const std::string &pluginName) const;

    CHEMKIT_DISABLE_COPY(PluginManager)

    friend class Plugin;

private:
    PluginManagerPrivate* const d;
};

} // end chemkit namespace

#include "pluginmanager-inline.h"

#endif // CHEMKIT_PLUGINMANAGER_H
