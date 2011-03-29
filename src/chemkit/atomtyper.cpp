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

#include "atomtyper.h"

#include "atom.h"
#include "foreach.h"
#include "molecule.h"
#include "pluginmanager.h"

namespace chemkit {

// === AtomTyperPrivate ==================================================== //
class AtomTyperPrivate
{
    public:
        std::string name;
        const Molecule *molecule;
};

// === AtomTyper =========================================================== //
/// \class AtomTyper atomtyper.h chemkit/atomtyper.h
/// \ingroup chemkit
/// \brief The AtomTyper class assigns symbolic types to atoms.
///
/// To create atom typer objects use the AtomTyper::create() method.

// --- Construction and Destruction ---------------------------------------- //
AtomTyper::AtomTyper(const std::string &name)
    : d(new AtomTyperPrivate)
{
    d->name = name;
    d->molecule = 0;
}

/// Destroys the atom typer object.
AtomTyper::~AtomTyper()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the atom typer.
std::string AtomTyper::name() const
{
    return d->name;
}

/// Sets the molecule for the atom typer to \p molecule.
void AtomTyper::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    assignTypes(molecule);
}

/// Returns the molecule for the atom typer.
const Molecule* AtomTyper::molecule() const
{
    return d->molecule;
}

// --- Types --------------------------------------------------------------- //
/// Returns the symbolic type for the atom at \p index.
QVariant AtomTyper::type(int index) const
{
    Q_UNUSED(index);

    return QVariant();
}

/// Returns the symbolic type for \p atom.
QVariant AtomTyper::type(const Atom *atom) const
{
    return type(atom->index());
}

/// Returns the symbolic type for the atom at \p index as an integer.
int AtomTyper::typeNumber(int index) const
{
    return type(index).toInt();
}

/// Returns the symbolic type for \p atom as an integer.
int AtomTyper::typeNumber(const Atom *atom) const
{
    return type(atom).toInt();
}

/// Returns the symbolic type for the atom at \p index as a string.
std::string AtomTyper::typeString(int index) const
{
    return type(index).toString().toStdString();
}

/// Returns the symbolic type for \p atom as a string.
std::string AtomTyper::typeString(const Atom *atom) const
{
    return type(atom).toString().toStdString();
}

// --- Static Methods ------------------------------------------------------ //
/// Creates and returns a new atom typer with \p name. Returns \c 0
/// if \p name is invalid.
AtomTyper* AtomTyper::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<AtomTyper>(name);
}

/// Returns a list of names of all the available atom typers.
std::vector<std::string> AtomTyper::typers()
{
    return PluginManager::instance()->pluginClassNames<AtomTyper>();
}

// --- Internal Methods ---------------------------------------------------- //
void AtomTyper::assignTypes(const Molecule *molecule)
{
    Q_UNUSED(molecule);
}

} // end chemkit namespace
