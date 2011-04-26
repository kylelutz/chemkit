/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
/// Returns the \c r (real) component.
template<typename T>
inline T GenericQuaternion<T>::r() const
{
    return (*this)[3];
}

/// Returns the imaginary part of the quaternion as a point.
template<typename T>
inline GenericPoint<T> GenericQuaternion<T>::toPoint3() const
{
    return GenericPoint<T>(StaticVector<T, 4>::x(),
                           StaticVector<T, 4>::y(),
                           StaticVector<T, 4>::z());
}

/// Returns the imaginary part of the quaternion as a vector.
template<typename T>
inline GenericVector<T> GenericQuaternion<T>::toVector3() const
{
    return GenericVector<T>(StaticVector<T, 4>::x(),
                            StaticVector<T, 4>::y(),
                            StaticVector<T, 4>::z());
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
    return GenericQuaternion<T>(-StaticVector<T, 4>::x(),
                                -StaticVector<T, 4>::y(),
                                -StaticVector<T, 4>::z(),
                                r());
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

// --- Static Methods ------------------------------------------------------ //
template<typename T>
inline GenericQuaternion<T> GenericQuaternion<T>::rotation(const GenericVector<T> &axis, T angle)
{
    return rotationRadians(axis, angle * chemkit::constants::DegreesToRadians);
}

template<typename T>
inline GenericQuaternion<T> GenericQuaternion<T>::rotationRadians(const GenericVector<T> &axis, T angle)
{
    return GenericQuaternion<T>(axis.x() * sin(angle/2.0),
                                axis.y() * sin(angle/2.0),
                                axis.z() * sin(angle/2.0),
                                cos(angle/2.0));
}

template<typename T>
inline GenericPoint<T> GenericQuaternion<T>::rotate(const GenericPoint<T> &point, const GenericVector<T> &axis, T angle)
{
    return rotateRadians(point, axis, angle * constants::DegreesToRadians);
}

template<typename T>
inline GenericPoint<T> GenericQuaternion<T>::rotateRadians(const GenericPoint<T> &point, const GenericVector<T> &axis, T angle)
{
    GenericQuaternion<T> p(point.x(), point.y(), point.z(), 0);
    GenericQuaternion<T> q = rotationRadians(axis, angle);
    GenericQuaternion<T> qc = q.conjugate();

    GenericQuaternion<T> qp = q.multiply(p);
    GenericQuaternion<T> r = qp.multiply(qc);

    return r.toPoint3();
}

template<typename T>
inline GenericVector<T> GenericQuaternion<T>::rotate(const GenericVector<T> &vector, const GenericVector<T> &axis, T angle)
{
    return rotateRadians(vector, axis, angle * constants::DegreesToRadians);
}

template<typename T>
inline GenericVector<T> GenericQuaternion<T>::rotateRadians(const GenericVector<T> &vector, const GenericVector<T> &axis, T angle)
{
    GenericQuaternion<T> p(vector.x(), vector.y(), vector.z(), 0);
    GenericQuaternion<T> q = rotationRadians(axis, angle);
    GenericQuaternion<T> qc = q.conjugate();

    GenericQuaternion<T> qp = q.multiply(p);
    GenericQuaternion<T> r = qp.multiply(qc);

    return r.toVector3();
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICQUATERNION_INLINE_H
