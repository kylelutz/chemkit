
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
}

/// Returns the molecule for the atom typer.
const Molecule* AtomTyper::molecule() const
{
    return d->molecule;
}

// --- Types --------------------------------------------------------------- //
/// Returns the symbolic type for \p atom.
Variant AtomTyper::type(const Atom *atom) const
{
    CHEMKIT_UNUSED(atom);

    return Variant();
}

/// Returns the symbolic type for \p atom as an integer.
int AtomTyper::typeNumber(const Atom *atom) const
{
    return type(atom).toInt();
}

/// Returns the symbolic type for \p atom as a string.
std::string AtomTyper::typeString(const Atom *atom) const
{
    return type(atom).toString();
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

} // end chemkit namespace
