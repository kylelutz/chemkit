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

#ifndef CHEMKIT_GENERICPOINT_INLINE_H
#define CHEMKIT_GENERICPOINT_INLINE_H

#include "genericpoint.h"

namespace chemkit {

// === GenericPoint ======================================================== //
/// \class GenericPoint genericpoint.h chemkit/genericpoint.h
/// \ingroup chemkit
/// \brief The GenericPoint class provides a template for a point in
///        three-dimensional space.
///
/// The GenericPoint template has one parameter:
///     - \b T: The coordinate type.
///
/// \see Point

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty generic point.
template<typename T>
inline GenericPoint<T>::GenericPoint()
    : StaticVector<T, 3>()
{
}

/// Creates a new generic point with components (\p x, \p y, \p z).
template<typename T>
inline GenericPoint<T>::GenericPoint(T x, T y, T z)
    : StaticVector<T, 3>()
{
    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
}

template<typename T>
inline GenericPoint<T>::GenericPoint(const StaticVector<float, 3> &vector)
    : StaticVector<T, 3>(vector)
{
}

template<typename T>
inline GenericPoint<T>::GenericPoint(const StaticVector<double, 3> &vector)
    : StaticVector<T, 3>(vector)
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the x component of the point.
template<typename T>
inline T GenericPoint<T>::x() const
{
    return (*this)[0];
}

/// Returns the y component of the point.
template<typename T>
inline T GenericPoint<T>::y() const
{
    return (*this)[1];
}

/// Returns the z component of the point.
template<typename T>
inline T GenericPoint<T>::z() const
{
    return (*this)[2];
}

/// Moves the point by (\p dx, \p dy, \p dz).
template<typename T>
inline void GenericPoint<T>::moveBy(T dx, T dy, T dz)
{
    (*this)[0] += dx;
    (*this)[1] += dy;
    (*this)[2] += dz;
}

/// Moves to point by \p vector.
template<typename T>
inline void GenericPoint<T>::moveBy(const StaticVector<T, 3> &vector)
{
    *this += vector;
}

/// Removes the point by \p distance along \p direction.
template<typename T>
inline void GenericPoint<T>::moveBy(T distance, const StaticVector<T, 3> &direction)
{
    *this += direction.normalized().scaled(distance);
}

/// Returns a new point moved by (\p dx, \p dy, \p dz).
template<typename T>
inline GenericPoint<T> GenericPoint<T>::movedBy(T dx, T dy, T dz) const
{
    GenericPoint<T> point = *this;
    point.moveBy(dx, dy, dz);
    return point;
}

/// Returns a new point moved by \p vector.
template<typename T>
inline GenericPoint<T> GenericPoint<T>::movedBy(const StaticVector<T, 3> &vector) const
{
    GenericPoint<T> point = *this;
    point.moveBy(vector);
    return point;
}

/// Returns a new point moved by \p distance along \p direction.
template<typename T>
inline GenericPoint<T> GenericPoint<T>::movedBy(T distance, const StaticVector<T, 3> &direction) const
{
    GenericPoint<T> point = *this;
    point.moveBy(distance, direction);
    return point;
}

// --- Math ---------------------------------------------------------------- //
/// Returns the distance between the point and \p point.
template<typename T>
inline T GenericPoint<T>::distance(const GenericPoint<T> &point) const
{
    return distance(*this, point);
}

/// Returns a point midway between the point and \p point.
template<typename T>
inline GenericPoint<T> GenericPoint<T>::midpoint(const GenericPoint<T> &point) const
{
    return midpoint(*this, point);
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the distance between points \p a and \p b.
template<typename T>
inline T GenericPoint<T>::distance(const GenericPoint<T> &a, const GenericPoint<T> &b)
{
    return (b - a).length();
}

/// Returns the sqaure distance between points \p a and \p b.
template<typename T>
inline T GenericPoint<T>::distanceSquared(const GenericPoint<T> &a, const GenericPoint<T> &b)
{
    return (b - a).lengthSquared();
}

/// Returns the angle between the vectors (\p a, \p b) and (\p b,
/// \p c). Angle is in Degrees.
template<typename T>
inline T GenericPoint<T>::angle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c)
{
    return angleRadians(a, b, c) * chemkit::constants::RadiansToDegrees;
}

/// Returns the angle between the vectors (\p a, \p b) and (\p b,
/// \p c). Angle is Radians.
template<typename T>
inline T GenericPoint<T>::angleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c)
{
    StaticVector<T, 3> ba = b - a;
    StaticVector<T, 3> bc = b - c;

    return acos(ba.dot(bc) / (ba.length() * bc.length()));
}

/// Returns the torsion angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Degrees.
template<typename T>
inline T GenericPoint<T>::torsionAngle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d)
{
    return torsionAngleRadians(a, b, c, d) * chemkit::constants::RadiansToDegrees;
}

/// Returns the torsion angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Degrees.
template<typename T>
inline T GenericPoint<T>::torsionAngleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d)
{
    StaticVector<T, 3> ab = b - a;
    StaticVector<T, 3> bc = c - b;
    StaticVector<T, 3> cd = d - c;

    return atan2(bc.length() * ab.dot(bc.cross(cd)), ab.cross(bc).dot(bc.cross(cd)));
}

/// Returns the wilson angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Degrees.
template<typename T>
inline T GenericPoint<T>::wilsonAngle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d)
{
    return wilsonAngleRadians(a, b, c, d) * chemkit::constants::RadiansToDegrees;
}

/// Returns the wilson angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Radians.
template<typename T>
inline T GenericPoint<T>::wilsonAngleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d)
{
    StaticVector<T, 3> normal = ((b - a).cross(c - b)).normalized();

    StaticVector<T, 3> bd = d - b;
    bd.normalize();

    Float angle = acos(bd.dot(normal));

    return (chemkit::constants::Pi * 0.5) - angle;
}

/// Returns the midpoint between \p a and \p b.
template<typename T>
inline GenericPoint<T> GenericPoint<T>::midpoint(const GenericPoint<T> &a, const GenericPoint<T> &b)
{
    return (a + b).scaled(0.5);
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICPOINT_INLINE_H
