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

#include "polymerfileformat.h"

#include "pluginmanager.h"

namespace chemkit {

namespace {

QHash<QString, PolymerFileFormat::CreateFunction> pluginFormats;

} // end anonymous namespace

// === PolymerFileFormatPrivate ============================================ //
class PolymerFileFormatPrivate
{
    public:
        QString name;
        QString errorString;
};

// === PolymerFileFormat =================================================== //
/// \class PolymerFileFormat polymerfileformat.h chemkit/polymerfileformat.h
/// \ingroup chemkit
/// \brief The PolymerFileFormat class represents a polymer file
///        format.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new polymer file format with \p name.
PolymerFileFormat::PolymerFileFormat(const QString &name)
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
QString PolymerFileFormat::name() const
{
    return d->name;
}

// --- Input and Output ---------------------------------------------------- //
/// Reads a file from \p iodev.
bool PolymerFileFormat::read(QIODevice *iodev, PolymerFile *file)
{
    Q_UNUSED(iodev);
    Q_UNUSED(file);

    setErrorString(QString("'%1' reading not supported.").arg(name()));
    return false;
}

/// Writes a file to \p iodev.
bool PolymerFileFormat::write(const PolymerFile *file, QIODevice *iodev)
{
    Q_UNUSED(file);
    Q_UNUSED(iodev);

    setErrorString(QString("'%1' writing not supported.").arg(name()));
    return false;
}

// --- Error Handling ------------------------------------------------------ //
void PolymerFileFormat::setErrorString(const QString &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
QString PolymerFileFormat::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new polymer file format with \p name. Returns \c 0 if
/// \p name is invalid.
PolymerFileFormat* PolymerFileFormat::create(const QString &name)
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    CreateFunction createFunction = pluginFormats.value(name.toLower());
    if(createFunction){
        return createFunction();
    }

    return 0;
}

/// Returns a list of available polymer file formats.
QStringList PolymerFileFormat::formats()
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    return pluginFormats.keys();
}

void PolymerFileFormat::registerFormat(const QString &name, CreateFunction function)
{
    pluginFormats.insert(name.toLower(), function);
}

void PolymerFileFormat::unregisterFormat(const QString &name, CreateFunction function)
{
    if(pluginFormats.value(name.toLower()) == function){
        pluginFormats.remove(name.toLower());
    }
}

} // end chemkit namespace
