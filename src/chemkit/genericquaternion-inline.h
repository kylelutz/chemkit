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

#ifndef CHEMKIT_GENERICQUATERNION_INLINE_H
#define CHEMKIT_GENERICQUATERNION_INLINE_H

#include "genericquaternion.h"

namespace chemkit {

// === GenericQuaternion =================================================== //
/// \class GenericQuaternion genericquaternion.h chemkit/genericquaternion.h
/// \ingroup chemkit
/// \brief The GenericQuaternion class provides a generic template
///        for quaternions.
///
/// The GenericQuaternion template has one parameter:
///     - \b T: The coordinate type.
///
/// \see Quaternion

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, null quaternion.
template<typename T>
inline GenericQuaternion<T>::GenericQuaternion()
    : StaticVector<T, 4>()
{
}

/// Creates a new quaternion with imaginary components (x, y, z) and
/// real component r.
template<typename T>
inline GenericQuaternion<T>::GenericQuaternion(T x, T y, T z, T r)
    : StaticVector<T, 4>()
{
    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
    (*this)[3] = r;
}

/// Creates a new generic quaternion as a copy of \p quaternion.
template<typename T>
inline GenericQuaternion<T>::GenericQuaternion(const StaticVector<T, 4> &quaternion)
    : StaticVector<T, 4>(quaternion)
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the \c x component.
template<typename T>
inline T GenericQuaternion<T>::x() const
{
    return (*this)[0];
}

/// Returns the \c y component.
template<typename T>
inline T GenericQuaternion<T>::y() const
{
    return (*this)[1];
}

/// Returns the \c z component.
template<typename T>
inline T GenericQuaternion<T>::z() const
{
    return (*this)[2];
}

/// Returns the \c r (real) component.
template<typename T>
inline T GenericQuaternion<T>::r() const
{
    return (*this)[3];
}

// --- Math ---------------------------------------------------------------- //
/// Returns the product of the quaternion mulitplied by \p quaternion.
///
/// The expression \c p.multiply(q) will return the result of
/// \f$ p \cdot q \f$.
template<typename T>
inline GenericQuaternion<T> GenericQuaternion<T>::multiply(const GenericQuaternion<T> &quaternion) const
{
    return multiply(*this, quaternion);
}

/// Returns the conjugate of the quaternion.
///
/// The expression \c q.conjugate() will return \f$ q^{*} \f$.
template<typename T>
inline GenericQuaternion<T> GenericQuaternion<T>::conjugate() const
{
    return GenericQuaternion<T>(-x(), -y(), -z(), r());
}

// --- Operators ----------------------------------------------------------- //
template<typename T>
inline GenericQuaternion<T> GenericQuaternion<T>::operator*(const GenericQuaternion<T> &quaternion) const
{
    return multiply(quaternion);
}

// --- Static Methods ------------------------------------------------------ //
template<typename T>
inline GenericQuaternion<T> GenericQuaternion<T>::multiply(const GenericQuaternion<T> &a, const GenericQuaternion<T> &b)
{
    T x = (a.r()*b.x() + a.x()*b.r() + a.y()*b.z() - a.z()*b.y());
    T y = (a.r()*b.y() + a.y()*b.r() + a.z()*b.x() - a.x()*b.z());
    T z = (a.r()*b.z() + a.z()*b.r() + a.x()*b.y() - a.y()*b.x());
    T r = (a.r()*b.r() - a.x()*b.x() - a.y()*b.y() - a.z()*b.z());

    return GenericQuaternion<T>(x, y, z, r);
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICQUATERNION_INLINE_H
