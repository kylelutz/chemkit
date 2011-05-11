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
/// \see Point3

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
    return (*this - point).length();
}

/// Returns a point midway between the point and \p point.
template<typename T>
inline GenericPoint<T> GenericPoint<T>::midpoint(const GenericPoint<T> &point) const
{
    return (*this + point).scaled(0.5);
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICPOINT_INLINE_H
