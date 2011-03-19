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

#ifndef CHEMKIT_QUATERNION_INLINE_H
#define CHEMKIT_QUATERNION_INLINE_H

#include "quaternion.h"

#include "constants.h"

namespace chemkit {

// === Quaternion ========================================================== //
/// \class Quaternion quaternion.h chemkit/quaternion.h
/// \ingroup chemkit
/// \brief The Quaternion class represents a quaternion.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new quaternion with imaginary components (x, y, z) and
/// real component r.
inline Quaternion::Quaternion(Float x, Float y, Float z, Float r)
    : GenericQuaternion<Float>(x, y, z, r)
{
    value(0) = x;
    value(1) = y;
    value(2) = z;
    value(3) = r;
}

/// Creates a new quaternion from point with real component r.
inline Quaternion::Quaternion(const Point3 &point, Float r)
    : GenericQuaternion<Float>(point.x(), point.y(), point.z(), r)
{
}

/// Creates a new quaternion from vector with real component r.
inline Quaternion::Quaternion(const Vector3 &vector, Float r)
    : GenericQuaternion<Float>(vector.x(), vector.y(), vector.z(), r)
{
}

inline Quaternion::Quaternion(const GenericQuaternion<Float> &quaternion)
    : GenericQuaternion<Float>(quaternion)
{
}

inline Quaternion::Quaternion(const StaticVector<Float, 4> &quaternion)
    : GenericQuaternion<Float>(quaternion)
{
}

// === Quaternionf ========================================================= //
/// \class Quaternionf quaternion.h chemkit/quaternion.h
/// \ingroup chemkit
/// \brief The Quaternionf class represents a quaternion.

// --- Construction and Destruction ---------------------------------------- //
inline Quaternionf::Quaternionf(float x, float y, float z, float r)
    : GenericQuaternion<float>(x, y, z, r)
{
}

inline Quaternionf::Quaternionf(const Point3f &point, float r)
    : GenericQuaternion<float>(point.x(), point.y(), point.z(), r)
{
}

inline Quaternionf::Quaternionf(const Vector3f &vector, float r)
    : GenericQuaternion<float>(vector.x(), vector.y(), vector.z(), r)
{
}

inline Quaternionf::Quaternionf(const GenericQuaternion<float> &quaternion)
    : GenericQuaternion<float>(quaternion)
{
}

inline Quaternionf::Quaternionf(const StaticVector<float, 4> &quaternion)
    : GenericQuaternion<float>(quaternion)
{
}

} // end chemkit namespace

#endif // CHEMKIT_QUATERNION_INLINE_H
