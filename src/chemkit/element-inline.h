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

#ifndef CHEMKIT_ELEMENT_INLINE_H
#define CHEMKIT_ELEMENT_INLINE_H

#include "element.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the atomic number of the element.
inline int Element::atomicNumber() const
{
    return m_atomicNumber;
}

/// Returns \c true if the element is valid.
inline bool Element::isValid() const
{
    return m_atomicNumber != 0;
}

// --- Operators ----------------------------------------------------------- //
inline bool Element::operator==(const Element &element) const
{
    return m_atomicNumber == element.m_atomicNumber;
}

} // end chemkit namespace

#endif // CHEMKIT_ELEMENT_INLINE_H
