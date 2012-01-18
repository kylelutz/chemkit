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

#ifndef CHEMKIT_PLUGIN_H
#define CHEMKIT_PLUGIN_H

#include "chemkit.h"

#include <string>

#include <boost/function.hpp>
#include <boost/lambda/construct.hpp>

namespace chemkit {

class PluginPrivate;
class DynamicLibrary;

class CHEMKIT_EXPORT Plugin
{
public:
    // properties
    std::string name() const;
    std::string fileName() const;
    std::string dataPath() const;

protected:
    // construction and destruction
    Plugin(const std::string &name);
    virtual ~Plugin();

    template<class T> bool registerPluginClass(const std::string &name, boost::function<T* ()> function);
    template<class T> bool unregisterPluginClass(const std::string &name);

private:
    void setLibrary(DynamicLibrary *library);
    DynamicLibrary* library() const;
    void addClassRegistration(const std::string &name, const std::string &className);
    void removeClassRegistration(const std::string &name, const std::string &className);

    friend class PluginManager;

private:
    PluginPrivate* const d;
};

} // end chemkit namespace

/// Export a plugin.
#define CHEMKIT_EXPORT_PLUGIN(name, className) \
    extern "C" CHEMKIT_DECL_EXPORT chemkit::Plugin* chemkit_plugin_init() \
    { \
        return new className; \
    }

/// Registers a plugin class with \p name.
///
/// This method must be called within the constructor of a
/// Plugin derived class.
#define CHEMKIT_REGISTER_PLUGIN_CLASS(name, baseClass, pluginClass) \
    registerPluginClass<baseClass>(name, boost::lambda::new_ptr<pluginClass>())

#include "plugin-inline.h"

#endif // CHEMKIT_PLUGIN_H
