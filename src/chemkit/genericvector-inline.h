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

#ifndef CHEMKIT_GENERICVECTOR_INLINE_H
#define CHEMKIT_GENERICVECTOR_INLINE_H

#include "genericvector.h"

#include <cstdlib>

namespace chemkit {

// === GenericVector ======================================================= //
/// \class GenericVector genericvector.h chemkit/genericvector.h
/// \ingroup chemkit
/// \brief The GenericVector class provides a generic template for
///        vectors in three-dimensional space.
///
/// The GenericVector template has one parameter:
///     - \b T: The coordinate type.
///
/// \see Vector3

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new generic vector.
template<typename T>
inline GenericVector<T>::GenericVector()
    : StaticVector<T, 3>()
{
}

/// Creates a new generic vector containing components (\p x, \p y,
/// \p z).
template<typename T>
inline GenericVector<T>::GenericVector(T x, T y, T z)
    : StaticVector<T, 3>()
{
    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
}

template<typename T>
inline GenericVector<T>::GenericVector(const StaticVector<float, 3> &vector)
    : StaticVector<T, 3>(vector)
{
}

template<typename T>
inline GenericVector<T>::GenericVector(const StaticVector<double, 3> &vector)
    : StaticVector<T, 3>(vector)
{
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a unit vector along the x-axis. (\c 1, \c 0, \c 0).
template<typename T>
inline GenericVector<T> GenericVector<T>::X()
{
    return GenericVector<T>(1, 0, 0);
}

/// Returns a unit vector along to y-axis. (\c 0, \c 1, \c 0).
template<typename T>
inline GenericVector<T> GenericVector<T>::Y()
{
    return GenericVector<T>(0, 1, 0);
}

/// Returns a unit vector along the z-axis. (\c 0, \c 0, \c 1).
template<typename T>
inline GenericVector<T> GenericVector<T>::Z()
{
    return GenericVector<T>(0, 0, 1);
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICVECTOR_INLINE_H
