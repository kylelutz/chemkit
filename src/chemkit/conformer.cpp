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

#include "conformer.h"

#include "molecule.h"

namespace chemkit {

// === ConformerPrivate ==================================================== //
class ConformerPrivate
{
    public:
        const Molecule *molecule;
        QHash<const Atom *, Point3> coordinates;
};

// === Conformer =========================================================== //
/// \class Conformer conformer.h chemkit/conformer.h
/// \ingroup chemkit
/// \brief The Conformer class represents an alternative set of atomic
///        coordinates for a molecule.
///
/// Conformer objects are created using the Molecule::addConformer()
/// method and destroyed with the Molecule::removeConformer() method.

// --- Construction and Destruction ---------------------------------------- //
Conformer::Conformer(const Molecule *molecule)
    : d(new ConformerPrivate)
{
    d->molecule = molecule;

    Q_FOREACH(const Atom *atom, molecule->atoms()){
        setPosition(atom, atom->position());
    }
}

Conformer::~Conformer()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the molecule for the conformer.
const Molecule* Conformer::molecule() const
{
    return d->molecule;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the coordinates for \p atom to \p position.
void Conformer::setPosition(const Atom *atom, const Point3 &position)
{
    d->coordinates[atom] = position;
}

/// Returns the position of the atom in the conformer.
Point3 Conformer::position(const Atom *atom) const
{
    return d->coordinates.value(atom, atom->position());
}

} // end chemkit namespace
