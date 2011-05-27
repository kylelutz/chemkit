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

#ifndef CHEMKIT_STATICVECTOR_INLINE_H
#define CHEMKIT_STATICVECTOR_INLINE_H

#include "staticvector.h"

#include <cmath>
#include <limits>
#include <cstdlib>

#ifndef QT_NO_DEBUG_STREAM
#include <QDebug>
#endif

namespace chemkit {

// === StaticVector ======================================================== //
/// \class StaticVector staticvector.h chemkit/staticvector.h
/// \ingroup chemkit
/// \brief The StaticVector template class implements a fixed-size vector.
///
/// The StaticVector template has two parameters:
///     - \b T: The data type.
///     - \b N: The size of the vector.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty vector.
template<typename T, int N>
inline StaticVector<T, N>::StaticVector()
{
    for(int i = 0; i < N; i++){
        m_data[i] = 0;
    }
}

template<typename T, int N>
inline StaticVector<T, N>::StaticVector(const StaticVector<float, N> &vector)
{
    for(int i = 0; i < N; i++){
        m_data[i] = vector[i];
    }
}

template<typename T, int N>
inline StaticVector<T, N>::StaticVector(const StaticVector<double, N> &vector)
{
    for(int i = 0; i < N; i++){
        m_data[i] = vector[i];
    }
}

// --- Properties ---------------------------------------------------------- //
/// Returns the x component of the vector.
template<typename T, int N>
inline T StaticVector<T, N>::x() const
{
    return N > 0 ? (*this)[0] : 0;
}

/// Returns the y component of the vector.
template<typename T, int N>
inline T StaticVector<T, N>::y() const
{
    return N > 1 ? (*this)[1] : 0;
}

/// Returns the z component of the vector.
template<typename T, int N>
inline T StaticVector<T, N>::z() const
{
    return N > 2 ? (*this)[2] : 0;
}

/// Returns the w component of the vector.
template<typename T, int N>
inline T StaticVector<T, N>::w() const
{
    return N > 3 ? (*this)[3] : 0;
}

/// Returns the size of the vector.
template<typename T, int N>
inline int StaticVector<T, N>::size() const
{
    return N;
}

/// Returns \c true if the vector contains all zeros.
template<typename T, int N>
inline bool StaticVector<T, N>::isZero() const
{
    for(int i = 0; i < N; i++){
        if(m_data[i] != 0.0){
            return false;
        }
    }

    return true;
}

/// Returns the data for the vector.
template<typename T, int N>
inline T* StaticVector<T, N>::data()
{
    return m_data;
}

/// \overload
template<typename T, int N>
inline const T* StaticVector<T, N>::data() const
{
    return m_data;
}

// --- Math ---------------------------------------------------------------- //
/// Returns the sum of the vector and \p vector.
///
/// The expression \c a.add(b) will return the result of
/// \f$ a + b \f$.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::add(const StaticVector<T, N> &vector) const
{
    StaticVector<T, N> product;

    for(int i = 0; i < N; i++){
        product[i] = m_data[i] + vector[i];
    }

    return product;
}

/// Returns the difference between the vector and \p vector.
///
/// The expression \c a.subtract(b) will return the result of
/// \f$ a - b \f$.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::subtract(const StaticVector<T, N> &vector) const
{
    StaticVector<T, N> product;

    for(int i = 0; i < N; i++){
        product[i] = m_data[i] - vector[i];
    }

    return product;
}

/// Returns the dot product of the vector with \p vector.
///
/// Given a vector \c a, the expression \c a.dot(b) will return the
/// result of \f$ a \cdot b \f$.
template<typename T, int N>
inline T StaticVector<T, N>::dot(const StaticVector<T, N> &vector) const
{
    T product = 0;

    for(int i = 0; i < N; i++){
        product += m_data[i] * vector[i];
    }

    return product;
}

/// Returns the cross product of the vector with \p vector.
/// Assumes that the vectors have three components (x, y, z).
///
/// The expression a.cross(b) will return the result of
/// \f$ a \times b \f$.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::cross(const StaticVector<T, N> &vector) const
{
    StaticVector<T, N> product;

    const StaticVector<T, N> &a = *this;
    const StaticVector<T, N> &b = vector;
    product[0] = a[1] * b[2] - a[2] * b[1];
    product[1] = a[2] * b[0] - a[0] * b[2];
    product[2] = a[0] * b[1] - a[1] * b[0];

    return product;
}

/// Returns the scalar triple product of the vector with \p a
/// and \p b.
///
/// The expression \c a.scalarTriple(b, c) will return the result of
/// \f$ a \cdot (b \times c) \f$.
template<typename T, int N>
inline T StaticVector<T, N>::scalarTriple(const StaticVector<T, N> &a, const StaticVector<T, N> &b) const
{
    return dot(a.cross(b));
}

/// Returns the norm of the vector.
///
/// The expression \c a.norm() will return the result of
/// \f$ \left|a\right| \f$.
template<typename T, int N>
inline T StaticVector<T, N>::norm() const
{
    T product = 0;

    for(int i = 0; i < N; i++){
        product += m_data[i] * m_data[i];
    }

    return sqrt(product);
}

/// Returns the squared norm of the vector.
///
/// The expression \c a.squaredNorm() will return the result of
/// \f$ \left|a\right|^{2} \f$.
template<typename T, int N>
inline T StaticVector<T, N>::squaredNorm() const
{
    T product = 0;

    for(int i = 0; i < N; i++){
        product += m_data[i] * m_data[i];
    }

    return product;
}

/// Normalizes the vector.
template<typename T, int N>
inline void StaticVector<T, N>::normalize()
{
    *this *= 1 / norm();
}

/// Returns a normalized version of the vector.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::normalized() const
{
    StaticVector<T, N> vector = *this;
    vector.normalize();
    return vector;
}

/// Returns the angle between the vector and \p vector. Angle is in degrees.
///
/// The expression \c a.angle(b) will return the result of
/// \f$ cos^{-1}(\frac{a \cdot b}{\left|a\right| \cdot \left|b\right|}) \cdot \frac{180}{\pi} \f$.
template<typename T, int N>
inline T StaticVector<T, N>::angle(const StaticVector<T, N> &vector) const
{
    return angleRadians(vector) * chemkit::constants::RadiansToDegrees;
}

/// Returns the angle between the vector and \p vector. Angle is in radians.
///
/// The expression \c a.angleRadians(b) will return the result of
/// \f$ cos^{-1}(\frac{a \cdot b}{\left|a\right| \cdot \left|b\right|}) \f$
template<typename T, int N>
inline T StaticVector<T, N>::angleRadians(const StaticVector<T, N> &vector) const
{
    return acos(dot(vector) / (norm() * vector.norm()));
}

// --- Operators ----------------------------------------------------------- //
/// Returns the value of the component at \p index.
///
/// \see value()
template<typename T, int N>
inline T StaticVector<T, N>::operator[](int index) const
{
    return m_data[index];
}

/// Returns the value of the component at \p index.
///
/// \see value()
template<typename T, int N>
inline T& StaticVector<T, N>::operator[](int index)
{
    return m_data[index];
}

/// Returns the sum of the vector and \p vector.
///
/// \see add()
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator+(const StaticVector<T, N> &vector) const
{
    return add(vector);
}

/// Returns the product of \c -1 and the vector.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator-() const
{
    return *this * -1;
}

/// Returns the difference between the vector and \p vector.
///
/// \see subtract()
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator-(const StaticVector<T, N> &vector) const
{
    return subtract(vector);
}

/// Returns the dot product of the vector and \p vector.
///
/// \see dot()
template<typename T, int N>
inline T StaticVector<T, N>::operator*(const StaticVector<T, N> &vector) const
{
    return dot(vector);
}

/// Returns the vector multiplied by \p scalar.
///
/// \see scaled()
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator*(T scalar) const
{
    StaticVector<T, N> product;

    for(int i = 0; i < N; i++){
        product[i] = m_data[i] * scalar;
    }

    return product;
}

/// Returns the vector divided by \p scalar.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator/(T scalar) const
{
    return *this * (1 / scalar);
}

/// Returns the cross product of the vector and \p vector.
///
/// \see cross()
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator^(const StaticVector<T, N> &vector) const
{
    return cross(vector);
}

/// Adds \p vector to the vector.
template<typename T, int N>
inline StaticVector<T, N>& StaticVector<T, N>::operator+=(const StaticVector<T, N> &vector)
{
    *this = *this + vector;
    return *this;
}

/// Subtracts \p vector from the vector.
template<typename T, int N>
inline StaticVector<T, N>& StaticVector<T, N>::operator-=(const StaticVector<T, N> &vector)
{
    *this = *this - vector;
    return *this;
}

/// Multiplies the vector by \p scalar.
template<typename T, int N>
inline StaticVector<T, N>& StaticVector<T, N>::operator*=(T scalar)
{
    *this = *this * scalar;
    return *this;
}

/// Divides the vector by \p scalar.
template<typename T, int N>
inline StaticVector<T, N>& StaticVector<T, N>::operator/=(T scalar)
{
    *this = *this / scalar;
    return *this;
}

/// Returns \c true if the vector is equal to \p vector.
template<typename T, int N>
inline bool StaticVector<T, N>::operator==(const StaticVector<T, N> &vector) const
{
    for(int i = 0; i < N; i++){
        if(std::abs(m_data[i] - vector[i]) > std::numeric_limits<T>::epsilon()){
            return false;
        }
    }

    return true;
}

template<typename T, int N>
inline CommaInitializer<T> StaticVector<T, N>::operator<<(const T value)
{
    m_data[0] = value;

    return CommaInitializer<T>(m_data, 1, N);
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a unit vector with a random direction.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::randomUnitVector()
{
    StaticVector<T, N> vector;

    for(int i = 0; i < N; i++){
        vector[i] = qrand() - (RAND_MAX / 2);
    }

    return vector.normalized();
}

// --- Related Functions --------------------------------------------------- //
/// Returns the product of \p scalar and \p vector.
///
/// \related StaticVector
template<typename T, int N>
inline StaticVector<T, N> operator*(T scalar, const StaticVector<T, N> &vector)
{
    return vector * scalar;
}

#ifndef QT_NO_DEBUG_STREAM
template<typename T, int N>
inline QDebug operator<<(QDebug debug, const StaticVector<T, N> &vector)
{
    debug.nospace() << "(";
    for(int i = 0; i < N; i++){
        if(i > 0){
            debug << ", ";
        }

        debug.nospace() << vector[i];
    }
    debug.nospace() << ")";

    return debug.space();
}
#endif // QT_NO_DEBUG_STREAM

} // end chemkit namespace

#endif // CHEMKIT_STATICVECTOR_INLINE_H
