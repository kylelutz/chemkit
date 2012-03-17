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

#include "lineformat.h"

#include <map>

#include <boost/format.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/algorithm/string.hpp>

#include "molecule.h"
#include "variantmap.h"
#include "pluginmanager.h"

namespace chemkit {

// === LineFormatPrivate =================================================== //
class LineFormatPrivate
{
public:
    std::string name;
    std::string errorString;
    VariantMap options;
};

// === LineFormat ========================================================== //
/// \class LineFormat lineformat.h chemkit/lineformat.h
/// \ingroup chemkit
/// \brief The LineFormat class provides a generic interface for
///        chemical line formats.
///
/// A list of supported line formats is available at:
/// http://wiki.chemkit.org/Features#Line_Formats

// --- Construction and Destruction ---------------------------------------- //
LineFormat::LineFormat(const std::string &name)
    : d(new LineFormatPrivate)
{
    d->name = boost::algorithm::to_lower_copy(name);
}

/// Destroys the line format object.
LineFormat::~LineFormat()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the line format.
std::string LineFormat::name() const
{
    return d->name;
}

// --- Options ------------------------------------------------------------- //
/// Sets an option for the line format.
void LineFormat::setOption(const std::string &name, const Variant &value)
{
    d->options[name] = value;
}

/// Returns the value of an option for the line format.
Variant LineFormat::option(const std::string &name) const
{
    VariantMap::const_iterator location = d->options.find(name);
    if(location != d->options.end()){
        return location->second;
    }
    else{
        return defaultOption(name);
    }
}

Variant LineFormat::defaultOption(const std::string &name) const
{
    CHEMKIT_UNUSED(name);

    return Variant();
}

// --- Input and Output ---------------------------------------------------- //
/// Reads and returns the molecule represented by the given
/// \p formula. Returns \c 0 if \p formula could not be
/// read.
Molecule* LineFormat::read(const std::string &formula)
{
    CHEMKIT_UNUSED(formula);

    setErrorString((boost::format("'%s' read not supported.") % name()).str());
    return 0;
}

/// Write and return the formula of a molecule.
std::string LineFormat::write(const Molecule *molecule)
{
    CHEMKIT_UNUSED(molecule);

    setErrorString((boost::format("'%s' write not supported.") % name()).str());
    return std::string();
}

// --- Error Handling ------------------------------------------------------ //
void LineFormat::setErrorString(const std::string &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occurred.
std::string LineFormat::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new line format object.
LineFormat* LineFormat::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<LineFormat>(name);
}

/// Returns a list of all the supported line formats.
std::vector<std::string> LineFormat::formats()
{
    return PluginManager::instance()->pluginClassNames<LineFormat>();
}

/// Converts \p formula in \p inputFormatName to \p outputFormatName.
std::string LineFormat::convert(const std::string &formula,
                                const std::string &inputFormatName,
                                const std::string &outputFormatName)
{
    // create input format
    boost::scoped_ptr<LineFormat> inputFormat(create(inputFormatName));
    if(!inputFormat){
        return std::string();
    }

    // read input formula
    boost::scoped_ptr<Molecule> molecule(inputFormat->read(formula));
    if(!molecule){
        return std::string();
    }

    // create output format
    boost::scoped_ptr<LineFormat> outputFormat(create(outputFormatName));
    if(!outputFormat){
        return std::string();
    }

    // return output formula
    return outputFormat->write(molecule.get());
}

} // end chemkit namespace
