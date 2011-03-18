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
inline Quaternion::Quaternion(const Vector &vector, Float r)
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

// --- Properties ---------------------------------------------------------- //
/// Returns the imaginary part of the quaternion as a point.
/// Equivalent to Point3(x(), y(), z()).
inline Point3 Quaternion::toPoint3() const
{
    return Point3(x(), y(), z());
}

/// Returns the imaginary part of the quaternion as a vector.
/// Equivalent to Vector(x(), y(), z()).
inline Vector Quaternion::toVector() const
{
    return Vector(x(), y(), z());
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a quaternion that represents the rotation of \p angle degrees
/// around \p axis.
inline Quaternion  Quaternion::rotation(const Vector &axis, Float angle)
{
    return rotationRadians(axis, angle * chemkit::constants::DegreesToRadians);
}

/// Returns a quaternion that represents the rotation of \p angle radians
/// around \p axis.
inline Quaternion Quaternion::rotationRadians(const Vector &axis, Float angle)
{
    return Quaternion(axis.x() * sin(angle/2.0),
                      axis.y() * sin(angle/2.0),
                      axis.z() * sin(angle/2.0),
                      cos(angle/2.0));
}

/// Returns a new point that is the result of rotating the point
/// angle degrees around axis.
inline Point3 Quaternion::rotate(const Point3 &point, const Vector &axis, Float angle)
{
    return rotateRadians(point, axis, angle * chemkit::constants::DegreesToRadians);
}

/// Returns a new point that is the result of rotating the point
/// angle radians around axis.
inline Point3 Quaternion::rotateRadians(const Point3 &point, const Vector &axis, Float angle)
{
    Quaternion p(point.x(), point.y(), point.z(), 0);
    Quaternion q = rotationRadians(axis, angle);
    Quaternion qc = q.conjugate();

    Quaternion qp = q.multiply(p);
    Quaternion r = qp.multiply(qc);

    return r.toPoint3();
}

/// Returns a new vector that is the result of rotating the vector
/// angle degrees around axis.
inline Vector Quaternion::rotate(const Vector &vector, const Vector &axis, Float angle)
{
    return rotateRadians(vector, axis, angle * chemkit::constants::DegreesToRadians);
}

/// Returns a new vector that is the result of rotating the vector
/// angle radians areound axis.
inline Vector Quaternion::rotateRadians(const Vector &vector, const Vector &axis, Float angle)
{
    Quaternion p(vector.x(), vector.y(), vector.z(), 0);
    Quaternion q  = rotationRadians(axis, angle);
    Quaternion qc = q.conjugate();

    Quaternion qp = q.multiply(p);
    Quaternion r = qp.multiply(qc);

    return r.toVector();
}

} // end chemkit namespace

#endif // CHEMKIT_QUATERNION_INLINE_H
