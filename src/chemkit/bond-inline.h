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

#ifndef CHEMKIT_BOND_INLINE_H
#define CHEMKIT_BOND_INLINE_H

#include "bond.h"

#include "atom.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the first atom in the bond.
inline Atom* Bond::atom1() const
{
    return m_atom1;
}

/// Returns the second atom in the bond.
inline Atom* Bond::atom2() const
{
    return m_atom2;
}

/// Returns the bond order.
inline int Bond::order() const
{
    return m_order;
}

/// Returns the molecule the bond is a part of.
inline Molecule* Bond::molecule() const
{
    return m_atom1->molecule();
}

} // end chemkit namespace

#endif // CHEMKIT_BOND_INLINE_H
