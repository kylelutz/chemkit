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

#include "pluginmanager.h"

namespace chemkit {

namespace {

QHash<QString, ChemicalFileFormat::CreateFunction> pluginFormats;

} // end anonymous namespace

// === ChemicalFileFormatPrivate =========================================== //
class ChemicalFileFormatPrivate
{
    public:
        QString name;
        QString errorString;
        QHash<QString, QVariant> options;
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
/// \see ChemicalFile, BiochemicalFileFormat

// --- Construction and Destruction ---------------------------------------- //
/// Construct a chemical file format.
ChemicalFileFormat::ChemicalFileFormat(const QString &name)
    : d(new ChemicalFileFormatPrivate)
{
    d->name = name.toLower();
}

/// Destroys a chemical file format.
ChemicalFileFormat::~ChemicalFileFormat()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the format.
QString ChemicalFileFormat::name() const
{
    return d->name;
}

// --- Options ------------------------------------------------------------- //
/// Sets an option for the format.
void ChemicalFileFormat::setOption(const QString &name, const QVariant &value)
{
    d->options[name] = value;
}

/// Returns the option for the format.
QVariant ChemicalFileFormat::option(const QString &name) const
{
    return d->options.value(name);
}

// --- Input and Output ---------------------------------------------------- //
/// Read from iodev into file.
bool ChemicalFileFormat::read(QIODevice *iodev, ChemicalFile *file)
{
    Q_UNUSED(iodev);
    Q_UNUSED(file);

    setErrorString(QString("'%1' reading not supported.").arg(name()));
    return false;
}

/// Write the contents of the file to iodev.
bool ChemicalFileFormat::write(const ChemicalFile *file, QIODevice *iodev)
{
    Q_UNUSED(file);
    Q_UNUSED(iodev);

    setErrorString(QString("'%1' writing not supported.").arg(name()));
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
ChemicalFileFormat* ChemicalFileFormat::create(const QString &name)
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    CreateFunction createFunction = pluginFormats.value(name.toLower());
    if(createFunction)
        return createFunction();

    return 0;
}

/// Returns a list of all supported file formats.
QStringList ChemicalFileFormat::formats()
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    return pluginFormats.keys();
}

void ChemicalFileFormat::registerFormat(const QString &name, CreateFunction function)
{
    pluginFormats[name.toLower()] = function;
}

void ChemicalFileFormat::unregisterFormat(const QString &name, CreateFunction function)
{
    Q_UNUSED(name);
    Q_UNUSED(function);
}

} // end chemkit namespace
