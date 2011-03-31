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

#include "chemicalfileformat.h"

#include <map>
#include <boost/algorithm/string/case_conv.hpp>

#include "pluginmanager.h"

namespace chemkit {

// === ChemicalFileFormatPrivate =========================================== //
class ChemicalFileFormatPrivate
{
    public:
        std::string name;
        QString errorString;
        std::map<std::string, QVariant> options;
};

// === ChemicalFileFormat ================================================== //
/// \class ChemicalFileFormat chemicalfileformat.h chemkit/chemicalfileformat.h
/// \ingroup chemkit
/// \brief The ChemicalFileFormat class represents a chemical file
///        format.
///
/// The ChemicalFileFormat class allows read and write access to a
/// chemical file's data. This class only deals with interpreting a
/// file format. To access the molecules contained in a file use the
/// ChemicalFile class.
///
/// \see ChemicalFile, PolymerFileFormat

// --- Construction and Destruction ---------------------------------------- //
/// Construct a chemical file format.
ChemicalFileFormat::ChemicalFileFormat(const std::string &name)
    : d(new ChemicalFileFormatPrivate)
{
    d->name = boost::algorithm::to_lower_copy(name);
}

/// Destroys a chemical file format.
ChemicalFileFormat::~ChemicalFileFormat()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the format.
std::string ChemicalFileFormat::name() const
{
    return d->name;
}

// --- Options ------------------------------------------------------------- //
/// Sets an option for the format.
void ChemicalFileFormat::setOption(const std::string &name, const QVariant &value)
{
    d->options[name] = value;
}

/// Returns the option for the format.
QVariant ChemicalFileFormat::option(const std::string &name) const
{
    std::map<std::string, QVariant>::iterator element = d->options.find(name);
    if(element != d->options.end()){
        return element->second;
    }

    return QVariant();
}

// --- Input and Output ---------------------------------------------------- //
/// Read from iodev into file.
bool ChemicalFileFormat::read(QIODevice *iodev, ChemicalFile *file)
{
    Q_UNUSED(iodev);
    Q_UNUSED(file);

    setErrorString(QString("'%1' reading not supported.").arg(name().c_str()));
    return false;
}

/// Write the contents of the file to iodev.
bool ChemicalFileFormat::write(const ChemicalFile *file, QIODevice *iodev)
{
    Q_UNUSED(file);
    Q_UNUSED(iodev);

    setErrorString(QString("'%1' writing not supported.").arg(name().c_str()));
    return false;
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a string describing the last error that occured.
void ChemicalFileFormat::setErrorString(const QString &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
QString ChemicalFileFormat::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new chemical file format.
ChemicalFileFormat* ChemicalFileFormat::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<ChemicalFileFormat>(name);
}

/// Returns a list of all supported file formats.
std::vector<std::string> ChemicalFileFormat::formats()
{
    return PluginManager::instance()->pluginClassNames<ChemicalFileFormat>();
}

} // end chemkit namespace
