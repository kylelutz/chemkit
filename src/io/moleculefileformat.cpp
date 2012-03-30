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

#include "moleculefileformat.h"

#include <map>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/variantmap.h>
#include <chemkit/pluginmanager.h>

namespace chemkit {

// === MoleculeFileFormatPrivate =========================================== //
class MoleculeFileFormatPrivate
{
public:
    std::string name;
    std::string errorString;
    VariantMap options;
};

// === MoleculeFileFormat ================================================== //
/// \class MoleculeFileFormat moleculefileformat.h chemkit/moleculefileformat.h
/// \ingroup chemkit-io
/// \brief The MoleculeFileFormat class represents a molecule file
///        format.
///
/// The MoleculeFileFormat class allows read and write access to a
/// molecule file's data. This class only deals with interpreting a
/// file format. To access the molecules contained in a file use the
/// MoleculeFile class.
///
/// A list of supported molecule file formats is available at:
/// http://wiki.chemkit.org/Features#Molecule_File_Formats
///
/// \see MoleculeFile, PolymerFileFormat

// --- Construction and Destruction ---------------------------------------- //
/// Construct a molecule file format.
MoleculeFileFormat::MoleculeFileFormat(const std::string &name)
    : d(new MoleculeFileFormatPrivate)
{
    d->name = boost::algorithm::to_lower_copy(name);
}

/// Destroys a molecule file format.
MoleculeFileFormat::~MoleculeFileFormat()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the format.
std::string MoleculeFileFormat::name() const
{
    return d->name;
}

// --- Options ------------------------------------------------------------- //
/// Sets an option for the format.
void MoleculeFileFormat::setOption(const std::string &name, const Variant &value)
{
    d->options[name] = value;
}

/// Returns the option for the format.
Variant MoleculeFileFormat::option(const std::string &name) const
{
    VariantMap::iterator element = d->options.find(name);
    if(element != d->options.end()){
        return element->second;
    }
    else{
        return defaultOption(name);
    }
}

Variant MoleculeFileFormat::defaultOption(const std::string &name) const
{
    CHEMKIT_UNUSED(name);

    return Variant();
}

// --- Input and Output ---------------------------------------------------- //
/// Read the data from \p input into \p file.
bool MoleculeFileFormat::read(std::istream &input, MoleculeFile *file)
{
    CHEMKIT_UNUSED(input);
    CHEMKIT_UNUSED(file);

    setErrorString((boost::format("'%s' reading not supported.") % name()).str());
    return false;
}

/// Read the data from \p input into \p file.
///
/// \internal
bool MoleculeFileFormat::readMappedFile(const boost::iostreams::mapped_file_source &input,
                                        MoleculeFile *file)
{
    CHEMKIT_UNUSED(input);
    CHEMKIT_UNUSED(file);

    setErrorString((boost::format("'%s' mapped file reading not supported.") % name()).str());
    return false;
}

/// Write the contents of \p file to \p output.
bool MoleculeFileFormat::write(const MoleculeFile *file, std::ostream &output)
{
    CHEMKIT_UNUSED(file);
    CHEMKIT_UNUSED(output);

    setErrorString((boost::format("'%s' writing not supported.") % name()).str());
    return false;
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a string describing the last error that occurred.
void MoleculeFileFormat::setErrorString(const std::string &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occurred.
std::string MoleculeFileFormat::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new molecule file format.
MoleculeFileFormat* MoleculeFileFormat::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<MoleculeFileFormat>(name);
}

/// Returns a list of all supported file formats.
std::vector<std::string> MoleculeFileFormat::formats()
{
    return PluginManager::instance()->pluginClassNames<MoleculeFileFormat>();
}

} // end chemkit namespace
