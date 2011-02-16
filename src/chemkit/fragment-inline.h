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

#ifndef CHEMKIT_FRAGMENT_INLINE_H
#define CHEMKIT_FRAGMENT_INLINE_H

#include "fragment.h"

#include "atom.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the fragment.
inline int Fragment::size() const
{
    return atomCount();
}

/// Returns the molecule the fragment is a part of.
inline Molecule* Fragment::molecule() const
{
    return m_atoms[0]->molecule();
}

// --- Structure ----------------------------------------------------------- //
/// Returns a list of all the atoms in the fragment.
inline QList<Atom *> Fragment::atoms() const
{
    return m_atoms;
}

/// Returns the number of atoms in the fragment.
inline int Fragment::atomCount() const
{
    return m_atoms.size();
}

/// Returns \c true if the fragment contains the atom.
inline bool Fragment::contains(const Atom *atom) const
{
    return m_atoms.contains(const_cast<Atom *>(atom));
}

} // end chemkit namespace

#endif // CHEMKIT_FRAGMENT_INLINE_H
