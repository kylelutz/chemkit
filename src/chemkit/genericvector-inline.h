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
inline GenericVector<T> GenericVector<T>::UnitX()
{
    return GenericVector<T>(1, 0, 0);
}

/// Returns a unit vector along to y-axis. (\c 0, \c 1, \c 0).
template<typename T>
inline GenericVector<T> GenericVector<T>::UnitY()
{
    return GenericVector<T>(0, 1, 0);
}

/// Returns a unit vector along the z-axis. (\c 0, \c 0, \c 1).
template<typename T>
inline GenericVector<T> GenericVector<T>::UnitZ()
{
    return GenericVector<T>(0, 0, 1);
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICVECTOR_INLINE_H
