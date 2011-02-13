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

#include "biochemicalfileformat.h"

#include "pluginmanager.h"

namespace chemkit {

namespace {

QHash<QString, BiochemicalFileFormat::CreateFunction> pluginFormats;

} // end anonymous namespace

// === BiochemicalFileFormatPrivate ======================================== //
class BiochemicalFileFormatPrivate
{
    public:
        QString name;
        QString errorString;
};

// === BiochemicalFileFormat =============================================== //
/// \class BiochemicalFileFormat biochemicalfileformat.h chemkit/biochemicalfileformat.h
/// \ingroup chemkit
/// \brief The BiochemicalFileFormat class represents a biochemical
///        file format.
///
/// The BiochemicalFileFormat class reads and writes biochemical
/// files. This class on handles interpreting the data in the
/// biochemical file. To access the proteins and nucleic acids the
/// file contains use the BiochemicalFile class.
///
/// \see BiochemicalFile, ChemicalFileFormat

// --- Construction and Destruction ---------------------------------------- //
/// Construct a new biochemical file format.
BiochemicalFileFormat::BiochemicalFileFormat(const QString &name)
    : d(new BiochemicalFileFormatPrivate)
{
    d->name = name;
}

/// Destroys a biochemical file format.
BiochemicalFileFormat::~BiochemicalFileFormat()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the file format.
QString BiochemicalFileFormat::name() const
{
    return d->name;
}

// --- Input and Output ---------------------------------------------------- //
/// Read from iodev into file.
bool BiochemicalFileFormat::read(QIODevice *iodev, BiochemicalFile *file)
{
    Q_UNUSED(iodev);
    Q_UNUSED(file);

    setErrorString(QString("'%1' reading not supported.").arg(name()));
    return false;
}

/// Write from file into iodev.
bool BiochemicalFileFormat::write(const BiochemicalFile *file, QIODevice *iodev)
{
    Q_UNUSED(file);
    Q_UNUSED(iodev);

    setErrorString(QString("'%1' writing not supported.").arg(name()));
    return false;
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a error string describing the last error that occured.
void BiochemicalFileFormat::setErrorString(const QString &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
QString BiochemicalFileFormat::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new biochemical file format.
BiochemicalFileFormat* BiochemicalFileFormat::create(const QString &name)
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    CreateFunction createFunction = pluginFormats.value(name.toLower());
    if(createFunction)
        return createFunction();

    return 0;
}

/// Returns a list of all supported biochemical file formats.
QStringList BiochemicalFileFormat::formats()
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    return pluginFormats.keys();
}

void BiochemicalFileFormat::registerFormat(const QString &name, CreateFunction function)
{
    pluginFormats[name.toLower()] = function;
}

void BiochemicalFileFormat::unregisterFormat(const QString &name, CreateFunction function)
{
    Q_UNUSED(name);
    Q_UNUSED(function);
}

} // end chemkit namespace
