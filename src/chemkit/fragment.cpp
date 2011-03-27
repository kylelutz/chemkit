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

#include "fragment.h"

#include "molecule.h"

namespace chemkit {

// === Fragment ============================================================ //
/// \class Fragment fragment.h chemkit/fragment.h
/// \ingroup chemkit
/// \brief The Fragment class represents a group of connected atoms in
///        a molecule.
///
/// Fragment objects are returned from the various fragment perception
/// methods such as Molecule::fragments() and Atom::fragment().

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new fragment that contains all the atoms attached to
/// \p root.
Fragment::Fragment(Atom *root)
{
    Q_ASSERT(root->m_fragment == 0);

    QList<Atom *> row;
    row.append(root);

    while(!row.isEmpty()){
        QList<Atom *> nextRow;

        Q_FOREACH(Atom *atom, row){
            if(!atom->m_fragment){
                atom->m_fragment = this;
                m_atoms.append(atom);

                Q_FOREACH(Atom *neighbor, atom->neighbors()){
                    nextRow.append(neighbor);
                }
            }
        }

        row = nextRow;
    }
}

/// Destroys the fragment object.
Fragment::~Fragment()
{
}

// --- Structure ----------------------------------------------------------- //
/// Returns a list of all the bonds in the fragment.
QList<Bond *> Fragment::bonds() const
{
    QList<Bond *> bonds;

    Q_FOREACH(Atom *atom, atoms()){
        Q_FOREACH(Bond *bond, atom->bonds()){
            if(!bonds.contains(bond)){
                bonds.append(bond);
            }
        }
    }

    return bonds;
}

/// Returns the number of bonds in the fragment.
int Fragment::bondCount() const
{
    return bonds().size();
}

/// Returns \c true if the fragment contains the bond.
bool Fragment::contains(const Bond *bond) const
{
    return bond->atom1()->fragment() == this;
}

} // end chemkit namespace
