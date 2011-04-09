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

#ifndef CHEMKIT_RING_INLINE_H
#define CHEMKIT_RING_INLINE_H

#include "ring.h"

#include <algorithm>

#include "atom.h"
#include "bond.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the ring.
inline int Ring::size() const
{
    return atomCount();
}

/// Returns the molecule the ring is a part of.
inline Molecule* Ring::molecule() const
{
    return m_atoms[0]->molecule();
}

/// Returns the fragment the ring is a part of.
inline Fragment* Ring::fragment() const
{
    return m_atoms[0]->fragment();
}

// --- Structure ----------------------------------------------------------- //
/// Returns the atom at \p index in the ring.
inline Atom* Ring::atom(int index) const
{
    return m_atoms[index];
}

/// Returns the atoms in the ring.
inline std::vector<Atom *> Ring::atoms() const
{
    return m_atoms;
}

/// Returns the number of atoms in the ring.
inline int Ring::atomCount() const
{
    return atoms().size();
}

/// Returns \c true if the ring contains atom.
inline bool Ring::contains(const Atom *atom) const
{
    return std::find(m_atoms.begin(), m_atoms.end(), atom) != m_atoms.end();
}

/// Returns \c true if the ring contains bond.
inline bool Ring::contains(const Bond *bond) const
{
    return contains(bond->atom1()) && contains(bond->atom2());
}

/// Returns \c true if the ring contains an atom with atomicNumber.
inline bool Ring::contains(const Element &element) const
{
    for(std::vector<Atom *>::const_iterator iter = m_atoms.begin(); iter != m_atoms.end(); ++iter){
        if((*iter)->is(element)){
            return true;
        }
    }

    return false;
}

} // end chemkit namespace

#endif // CHEMKIT_RING_INLINE_H
