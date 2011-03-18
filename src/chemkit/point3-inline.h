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

#ifndef CHEMKIT_POINT3_INLINE_H
#define CHEMKIT_POINT3_INLINE_H

#include "point3.h"

namespace chemkit {

// === Point3 ============================================================== //
/// \class Point3 point3.h chemkit/point3.h
/// \ingroup chemkit
/// \brief The Point3 class represent a location in three-dimensional
///        space.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new point at (\c 0, \c 0, \c 0).
inline Point3::Point3()
    : GenericPoint<Float>()
{
}

/// Create a new point at (\p x, \p y, \p z).
inline Point3::Point3(Float x, Float y, Float z)
    : GenericPoint<Float>(x, y, z)
{
}

inline Point3::Point3(const GenericPoint<Float> &point)
    : GenericPoint<Float>(point)
{
}

inline Point3::Point3(const StaticVector<Float, 3> &vector)
    : GenericPoint<Float>(vector)
{
}

} // end chemkit namespace

#endif // CHEMKIT_POINT3_INLINE_H
