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

#include "foreach.h"
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

    std::vector<Atom *> row;
    row.push_back(root);

    while(!row.empty()){
        std::vector<Atom *> nextRow;

        foreach(Atom *atom, row){
            if(!atom->m_fragment){
                atom->m_fragment = this;
                m_atoms.push_back(atom);

                foreach(Atom *neighbor, atom->neighbors()){
                    nextRow.push_back(neighbor);
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
std::vector<Bond *> Fragment::bonds() const
{
    std::vector<Bond *> bonds;

    foreach(Atom *atom, m_atoms){
        foreach(Bond *bond, atom->bonds()){
            if(std::find(bonds.begin(), bonds.end(), bond) == bonds.end()){
                bonds.push_back(bond);
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
