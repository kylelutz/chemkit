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

#ifndef CHEMKIT_ATOM_INLINE_H
#define CHEMKIT_ATOM_INLINE_H

#include "atom.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the atom's element.
inline Element Atom::element() const
{
    return m_element;
}

/// Returns the molecule the atom is a part of.
inline Molecule* Atom::molecule() const
{
    return m_molecule;
}

/// Returns \c true if the atom's element is the same as \p element.
///
/// For example, to check if an atom is carbon or hydrogen:
/// \code
/// if(atom->is(Atom::Carbon) || atom->is(Atom::Hydrogen))
///     // yes
/// else
///     //no
/// \endcode
inline bool Atom::is(const Element &element) const
{
    return m_element == element;
}

} // end chemkit namespace

#endif // CHEMKIT_ATOM_INLINE_H
