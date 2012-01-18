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

#include "plugin.h"

#include <boost/filesystem.hpp>

#include "foreach.h"
#include "dynamiclibrary.h"

namespace chemkit {

// === PluginPrivate ======================================================= //
class PluginPrivate
{
public:
    std::string name;
    DynamicLibrary *library;
    std::vector<std::pair<std::string, std::string> > pluginClasses;
};

// === Plugin ============================================================== //
/// \class Plugin plugin.h chemkit/plugin.h
/// \ingroup chemkit
/// \brief The Plugin class is the base class for all dynamically
///        loaded chemkit plugins.
///
/// \sa PluginManager

// --- Construction and Destruction ---------------------------------------- //
Plugin::Plugin(const std::string &name)
    : d(new PluginPrivate)
{
    d->name = name;
    d->library = 0;
}

Plugin::~Plugin()
{
    // unregister plugin classes
    std::pair<std::string, std::string> pluginClass;
    foreach(pluginClass, d->pluginClasses){
        PluginManager::instance()->unregisterPluginClass(pluginClass.second,
                                                         pluginClass.first);
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the plugin.
std::string Plugin::name() const
{
    return d->name;
}

std::string Plugin::fileName() const
{
    if(!d->library){
        return std::string();
    }

    return d->library->fileName();
}

std::string Plugin::dataPath() const
{
    if(!d->library){
        return std::string();
    }

    boost::filesystem::path path(d->library->fileName());

    return (path.parent_path() / "data" / d->name / "/").string();
}

// --- Internal Methods ---------------------------------------------------- //
void Plugin::setLibrary(DynamicLibrary *library)
{
    d->library = library;
}

DynamicLibrary* Plugin::library() const
{
    return d->library;
}

void Plugin::addClassRegistration(const std::string &name, const std::string &className)
{
    d->pluginClasses.push_back(std::make_pair(name, className));
}

void Plugin::removeClassRegistration(const std::string &name, const std::string &className)
{
    std::pair<std::string, std::string> value(name, className);

    d->pluginClasses.erase(std::remove(d->pluginClasses.begin(),
                                       d->pluginClasses.end(),
                                       value),
                           d->pluginClasses.end());
}

} // end chemkit namespace
