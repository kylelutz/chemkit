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

#include "moleculefileformat.h"

#include <map>
#include <boost/algorithm/string/case_conv.hpp>

#include "pluginmanager.h"

namespace chemkit {

// === MoleculeFileFormatPrivate =========================================== //
class MoleculeFileFormatPrivate
{
    public:
        std::string name;
        std::string errorString;
        std::map<std::string, QVariant> options;
};

// === MoleculeFileFormat ================================================== //
/// \class MoleculeFileFormat moleculefileformat.h chemkit/moleculefileformat.h
/// \ingroup chemkit
/// \brief The MoleculeFileFormat class represents a molecule file
///        format.
///
/// The MoleculeFileFormat class allows read and write access to a
/// molecule file's data. This class only deals with interpreting a
/// file format. To access the molecules contained in a file use the
/// MoleculeFile class.
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
void MoleculeFileFormat::setOption(const std::string &name, const QVariant &value)
{
    d->options[name] = value;
}

/// Returns the option for the format.
QVariant MoleculeFileFormat::option(const std::string &name) const
{
    std::map<std::string, QVariant>::iterator element = d->options.find(name);
    if(element != d->options.end()){
        return element->second;
    }

    return QVariant();
}

// --- Input and Output ---------------------------------------------------- //
/// Read from iodev into file.
bool MoleculeFileFormat::read(QIODevice *iodev, MoleculeFile *file)
{
    Q_UNUSED(iodev);
    Q_UNUSED(file);

    setErrorString(QString("'%1' reading not supported.").arg(name().c_str()).toStdString());
    return false;
}

/// Write the contents of the file to iodev.
bool MoleculeFileFormat::write(const MoleculeFile *file, QIODevice *iodev)
{
    Q_UNUSED(file);
    Q_UNUSED(iodev);

    setErrorString(QString("'%1' writing not supported.").arg(name().c_str()).toStdString());
    return false;
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a string describing the last error that occured.
void MoleculeFileFormat::setErrorString(const std::string &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
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
