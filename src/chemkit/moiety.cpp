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

#include "moiety.h"

#include "atom.h"

namespace chemkit {

// === MoietyPrivate ======================================================= //
class MoietyPrivate
{
    public:
        QList<const Atom *> atoms;
};

// === Moiety ============================================================== //
/// \class Moiety moiety.h chemkit/moiety.h
/// \ingroup chemkit
/// \brief The Moiety class represents a group of atoms in a
///        molecule.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty moiety object.
Moiety::Moiety()
    : d(new MoietyPrivate)
{
}

/// Creates a new moiety object containing \p atoms.
Moiety::Moiety(const QList<const Atom *> &atoms)
    : d(new MoietyPrivate)
{
    d->atoms = atoms;
}

/// Destroys the moiety object.
Moiety::~Moiety()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the moiety.
int Moiety::size() const
{
    return atomCount();
}

/// Returns \c true if the moiety contains no atoms.
bool Moiety::isEmpty() const
{
    return size() == 0;
}

/// Returns the molecule that the moiety is a part of.
const Molecule* Moiety::molecule() const
{
    return d->atoms[0]->molecule();
}

// --- Structure ----------------------------------------------------------- //
/// Returns the atom at \p index.
const Atom* Moiety::atom(int index) const
{
    return d->atoms[index];
}

/// Returns a list of the atoms in the moiety.
QList<const Atom *> Moiety::atoms() const
{
    return d->atoms;
}

/// Returns the number of atoms in the moiety.
int Moiety::atomCount() const
{
    return d->atoms.size();
}

// --- Operators ----------------------------------------------------------- //
Moiety& Moiety::operator=(const Moiety &moiety)
{
    d->atoms = moiety.d->atoms;

    return *this;
}

} // end chemkit namespace
