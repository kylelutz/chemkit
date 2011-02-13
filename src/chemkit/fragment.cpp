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

        foreach(Atom *atom, row){
            if(!atom->m_fragment){
                atom->m_fragment = this;
                m_atoms.append(atom);

                foreach(Atom *neighbor, atom->neighbors()){
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

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the fragment.
int Fragment::size() const
{
    return atomCount();
}

/// Returns the molecule the fragment is a part of.
Molecule* Fragment::molecule()
{
    return m_atoms[0]->molecule();
}

/// \overload
const Molecule* Fragment::molecule() const
{
    return const_cast<const Atom *>(m_atoms[0])->molecule();
}

// --- Structure ----------------------------------------------------------- //
/// Returns a list of all the atoms in the fragment.
QList<Atom *> Fragment::atoms()
{
    return m_atoms;
}

/// \overload
QList<const Atom *> Fragment::atoms() const
{
    QList<const Atom *> atoms;

    foreach(const Atom *atom, m_atoms){
        atoms.append(atom);
    }

    return atoms;
}

/// Returns the number of atoms in the fragment.
int Fragment::atomCount() const
{
    return atoms().size();
}

/// Returns \c true if the fragment contains the atom.
bool Fragment::contains(const Atom *atom) const
{
    return m_atoms.contains(const_cast<Atom *>(atom));
}

/// Returns a list of all the bonds in the fragment.
QList<Bond *> Fragment::bonds()
{
    QList<Bond *> bonds;

    foreach(Atom *atom, atoms()){
        foreach(Bond *bond, atom->bonds()){
            if(!bonds.contains(bond)){
                bonds.append(bond);
            }
        }
    }

    return bonds;
}

/// \overload
QList<const Bond *> Fragment::bonds() const
{
    QList<const Bond *> bonds;

    foreach(const Bond *bond, const_cast<Fragment *>(this)->bonds()){
        bonds.append(bond);
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
