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

#ifndef CHEMKIT_VECTOR_INLINE_H
#define CHEMKIT_VECTOR_INLINE_H

#include "vector.h"

namespace chemkit {

// === Vector ============================================================== //
/// \class Vector vector.h chemkit/vector.h
/// \ingroup chemkit
/// \brief The Vector class represents a direction in
///        three-dimensional space.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new vector containing (\c 0, \c 0, \c 0).
inline Vector::Vector()
    : GenericVector<Float>()
{
}

/// Create a new vector containing (\p x, \p y, \p z).
inline Vector::Vector(Float x, Float y, Float z)
    : GenericVector<Float>(x, y, z)
{
}

inline Vector::Vector(const GenericVector<Float> &vector)
    : GenericVector<Float>(vector)
{
}

inline Vector::Vector(const StaticVector<Float, 3> &vector)
    : GenericVector<Float>(vector)
{
}

} // end chemkit namespace

#endif // CHEMKIT_VECTOR_INLINE_H
