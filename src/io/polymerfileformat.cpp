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

#include "polymerfileformat.h"

#include <boost/format.hpp>

#include <chemkit/pluginmanager.h>

namespace chemkit {

// === PolymerFileFormatPrivate ============================================ //
class PolymerFileFormatPrivate
{
public:
    std::string name;
    std::string errorString;
};

// === PolymerFileFormat =================================================== //
/// \class PolymerFileFormat polymerfileformat.h chemkit/polymerfileformat.h
/// \ingroup chemkit-io
/// \brief The PolymerFileFormat class represents a polymer file
///        format.
///
/// A list of supported polymer file formats is available at:
/// http://wiki.chemkit.org/Features#Polymer_File_Formats

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new polymer file format with \p name.
PolymerFileFormat::PolymerFileFormat(const std::string &name)
    : d(new PolymerFileFormatPrivate)
{
    d->name = name;
}

/// Destroys the polymer file format object.
PolymerFileFormat::~PolymerFileFormat()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the file format.
std::string PolymerFileFormat::name() const
{
    return d->name;
}

// --- Input and Output ---------------------------------------------------- //
/// Read the data from \p input into \p file.
bool PolymerFileFormat::read(std::istream &input, PolymerFile *file)
{
    CHEMKIT_UNUSED(input);
    CHEMKIT_UNUSED(file);

    setErrorString((boost::format("'%1' reading not supported.") % name()).str());
    return false;
}

/// Read the data from \p input into \p file.
///
/// \internal
bool PolymerFileFormat::readMappedFile(const boost::iostreams::mapped_file_source &input,
                                       PolymerFile *file)
{
    CHEMKIT_UNUSED(input);
    CHEMKIT_UNUSED(file);

    setErrorString((boost::format("'%s' mapped file reading not supported.") % name()).str());
    return false;
}

/// Write the contents of \p file to \p output.
bool PolymerFileFormat::write(const PolymerFile *file, std::ostream &output)
{
    CHEMKIT_UNUSED(file);
    CHEMKIT_UNUSED(output);

    setErrorString((boost::format("'%1' writing not supported.") % name()).str());
    return false;
}

// --- Error Handling ------------------------------------------------------ //
void PolymerFileFormat::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occurred.
std::string PolymerFileFormat::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new polymer file format with \p name. Returns \c 0 if
/// \p name is invalid.
PolymerFileFormat* PolymerFileFormat::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<PolymerFileFormat>(name);
}

/// Returns a list of available polymer file formats.
std::vector<std::string> PolymerFileFormat::formats()
{
    return PluginManager::instance()->pluginClassNames<PolymerFileFormat>();
}

} // end chemkit namespace
