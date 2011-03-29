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

#include "lineformat.h"

#include "pluginmanager.h"

namespace chemkit {

namespace {

QHash<QString, LineFormat::CreateFunction> pluginFormats;

} // end anonymous namespace

// === LineFormatPrivate =================================================== //
class LineFormatPrivate
{
    public:
        std::string name;
        std::string errorString;
        QHash<QString, QVariant> options;
};

// === LineFormat ========================================================== //
/// \class LineFormat lineformat.h chemkit/lineformat.h
/// \ingroup chemkit
/// \brief The LineFormat class provides a generic interface for
///        chemical line formats.
///
/// The following line formats are supported in chemkit:
///   - \c inchi
///   - \c inchikey
///   - \c mcdl
///   - \c formula
///   - \c smiles

// --- Construction and Destruction ---------------------------------------- //
LineFormat::LineFormat(const std::string &name)
    : d(new LineFormatPrivate)
{
    d->name = QString(name.c_str()).toLower().toStdString();
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
void LineFormat::setOption(const std::string &name, const QVariant &value)
{
    d->options[name.c_str()] = value;
}

/// Returns the value of an option for the line format.
QVariant LineFormat::option(const std::string &name) const
{
    return d->options.value(name.c_str(), defaultOption(name));
}

QVariant LineFormat::defaultOption(const std::string &name) const
{
    Q_UNUSED(name);

    return QVariant();
}

// --- Input and Output ---------------------------------------------------- //
/// Reads \p formula and adds its contents to \p molecule. Returns
/// \c false if \p formula could not be read.
bool LineFormat::read(const std::string &formula, Molecule *molecule)
{
    Q_UNUSED(formula);
    Q_UNUSED(molecule);

    setErrorString(QString("'%1' read not supported.").arg(name().c_str()).toStdString());
    return false;
}

/// Reads and returns the molecule represented by the given
/// \p formula. Returns \c 0 if \p formula could not be
/// read.
Molecule* LineFormat::read(const std::string &formula)
{
    Molecule *molecule = new Molecule;

    bool ok = read(formula, molecule);
    if(!ok){
        delete molecule;
        return 0;
    }

    return molecule;
}

/// Write and return the formula of a molecule.
std::string LineFormat::write(const Molecule *molecule)
{
    Q_UNUSED(molecule);

    setErrorString(QString("'%1' write not supported.").arg(name().c_str()).toStdString());
    return std::string();
}

// --- Error Handling ------------------------------------------------------ //
void LineFormat::setErrorString(const std::string &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
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

} // end chemkit namespace
