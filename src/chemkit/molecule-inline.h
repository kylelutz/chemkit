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

#ifndef CHEMKIT_MOLECULE_INLINE_H
#define CHEMKIT_MOLECULE_INLINE_H

#include "molecule.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the molecule.
inline int Molecule::size() const
{
    return atomCount();
}

/// Returns \c true if the molecule contains no atoms (i.e.
/// size() == 0).
inline bool Molecule::isEmpty() const
{
    return size() == 0;
}

// --- Structure ----------------------------------------------------------- //
/// Returns a list of all the atoms in the molecule.
inline QList<Atom *> Molecule::atoms()
{
    return m_atoms;
}

/// \overload
inline QList<const Atom *> Molecule::atoms() const
{
    QList<const Atom *> atoms;

    foreach(const Atom *atom, m_atoms){
        atoms.append(atom);
    }

    return atoms;
}

/// Returns the number of atoms in the molecule.
inline int Molecule::atomCount() const
{
    return m_atoms.size();
}

/// Returns the atom at \p index.
inline Atom* Molecule::atom(int index)
{
    return m_atoms.value(index, 0);
}

/// \overload
inline const Atom* Molecule::atom(int index) const
{
    return m_atoms.value(index, 0);
}

} // end chemkit namespace

#endif // CHEMKIT_MOLECULE_INLINE_H
