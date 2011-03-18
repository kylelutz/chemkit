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

#ifndef CHEMKIT_QUATERNIONG_INLINE_H
#define CHEMKIT_QUATERNIONG_INLINE_H

#include "quaterniong.h"

#include <chemkit/constants.h>

namespace chemkit {

// === Quaterniong ========================================================= //
/// \class Quaterniong quaterniong.h chemkit/quateriong.h
/// \ingroup chemkit-graphics
/// \brief The Quaterniong class represents a quaternion.

// --- Construction and Destruction ---------------------------------------- //    
inline Quaterniong::Quaterniong(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z, GraphicsFloat r)
    : GenericQuaternion<GraphicsFloat>(x, y, z, r)
{
}

inline Quaterniong::Quaterniong(const Point3g &point, GraphicsFloat r)
    : GenericQuaternion<GraphicsFloat>(point.x(), point.y(), point.z(), r)
{
}

inline Quaterniong::Quaterniong(const Vector3g &vector, GraphicsFloat r)
    : GenericQuaternion<GraphicsFloat>(vector.x(), vector.y(), vector.z(), r)
{
}

inline Quaterniong::Quaterniong(const GenericQuaternion<GraphicsFloat> &quaternion)
    : GenericQuaternion<GraphicsFloat>(quaternion)
{
}

inline Quaterniong::Quaterniong(const StaticVector<GraphicsFloat, 4> &quaternion)
    : GenericQuaternion<GraphicsFloat>(quaternion)
{
}

// --- Properties ---------------------------------------------------------- //
inline Point3g Quaterniong::toPoint3() const
{
    return Point3g(x(), y(), z());
}

inline Vector3g Quaterniong::toVector3() const
{
    return Vector3g(x(), y(), z());
}

// --- Static Methods ------------------------------------------------------ //
inline Quaterniong Quaterniong::rotation(const Vector3g &axis, GraphicsFloat angle)
{
    return rotationRadians(axis, angle * chemkit::constants::DegreesToRadians);
}

inline Quaterniong Quaterniong::rotationRadians(const Vector3g &axis, GraphicsFloat angle)
{
    return Quaterniong(axis.x() * sin(angle/2.0),
                       axis.y() * sin(angle/2.0),
                       axis.z() * sin(angle/2.0),
                       cos(angle/2.0));
}

inline Point3g Quaterniong::rotate(const Point3g &point, const Vector3g &axis, GraphicsFloat angle)
{
    return rotateRadians(point, axis, angle * constants::DegreesToRadians);
}

inline Point3g Quaterniong::rotateRadians(const Point3g &point, const Vector3g &axis, GraphicsFloat angle)
{
    Quaterniong p(point.x(), point.y(), point.z(), 0);
    Quaterniong q = rotationRadians(axis, angle);
    Quaterniong qc = q.conjugate();

    Quaterniong qp = q.multiply(p);
    Quaterniong r = qp.multiply(qc);

    return r.toPoint3();
}

inline Vector3g Quaterniong::rotate(const Vector3g &vector, const Vector3g &axis, GraphicsFloat angle)
{
    return rotateRadians(vector, axis, angle * constants::DegreesToRadians);
}

inline Vector3g Quaterniong::rotateRadians(const Vector3g &vector, const Vector3g &axis, GraphicsFloat angle)
{
    Quaterniong p(vector.x(), vector.y(), vector.z(), 0);
    Quaterniong q = rotationRadians(axis, angle);
    Quaterniong qc = q.conjugate();

    Quaterniong qp = q.multiply(p);
    Quaterniong r = qp.multiply(qc);

    return r.toVector3();
}

} // end chemkit namespace

#endif // CHEMKIT_QUATERNIONG_INLINE_H
