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

// === PolymerFileFormatPrivate ============================================ //
class PolymerFileFormatPrivate
{
    public:
        std::string name;
        std::string errorString;
};

// === PolymerFileFormat =================================================== //
/// \class PolymerFileFormat polymerfileformat.h chemkit/polymerfileformat.h
/// \ingroup chemkit
/// \brief The PolymerFileFormat class represents a polymer file
///        format.

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
/// Reads a file from \p iodev.
bool PolymerFileFormat::read(QIODevice *iodev, PolymerFile *file)
{
    Q_UNUSED(iodev);
    Q_UNUSED(file);

    setErrorString(QString("'%1' reading not supported.").arg(name().c_str()).toStdString());
    return false;
}

/// Writes a file to \p iodev.
bool PolymerFileFormat::write(const PolymerFile *file, QIODevice *iodev)
{
    Q_UNUSED(file);
    Q_UNUSED(iodev);

    setErrorString(QString("'%1' writing not supported.").arg(name().c_str()).toStdString());
    return false;
}

// --- Error Handling ------------------------------------------------------ //
void PolymerFileFormat::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
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
