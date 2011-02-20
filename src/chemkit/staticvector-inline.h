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

#ifndef CHEMKIT_STATICVECTOR_INLINE_H
#define CHEMKIT_STATICVECTOR_INLINE_H

#include "staticvector.h"

#include <cmath>

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

/// Creates a new vector from \p data. The \p data pointer must
/// contain \p size number of items.
template<typename T, int N>
inline StaticVector<T, N>::StaticVector(const T *data, int size)
{
    for(int i = 0; i < size; i++){
        m_data[i] = data[i];
    }
}

template<typename T, int N>
inline StaticVector<T, N>::StaticVector(const StaticVector<float, N> &vector)
{
    for(int i = 0; i < N; i++){
        m_data[i] = vector.value(i);
    }
}

template<typename T, int N>
inline StaticVector<T, N>::StaticVector(const StaticVector<double, N> &vector)
{
    for(int i = 0; i < N; i++){
        m_data[i] = vector.value(i);
    }
}

// --- Properties ---------------------------------------------------------- //
/// Sets the value of the component at \p index to \p value.
template<typename T, int N>
inline void StaticVector<T, N>::setValue(int index, T value)
{
    m_data[index] = value;
}

/// Returns the value of the component at \p index.
template<typename T, int N>
inline T StaticVector<T, N>::value(int index) const
{
    return m_data[index];
}

/// \overload
template<typename T, int N>
inline T& StaticVector<T, N>::value(int index)
{
    return m_data[index];
}

/// Returns the size of the vector.
template<typename T, int N>
inline int StaticVector<T, N>::size() const
{
    return N;
}

/// Returns \c true if the vector contains all zeros.
template<typename T, int N>
inline bool StaticVector<T, N>::isNull() const
{
    for(int i = 0; i < N; i++){
        if(value(i) != 0.0){
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
        product[i] = value(i) + vector[i];
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
        product[i] = value(i) - vector[i];
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
        product += value(i) * vector[i];
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

/// Returns the vector triple product of the vector with \p a
/// and \p b.
///
/// The expression \c a.vectorTriple(b, c) will return the result of
/// \f$ a \times (b \times c) \f$.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::vectorTriple(const StaticVector<T, N> &a, const StaticVector<T, N> &b) const
{
    return cross(a.cross(b));
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
        product += value(i) * value(i);
    }

    return sqrt(product);
}

/// Returns the squared norm of the vector.
///
/// The expression \c a.normSquared() will return the result of
/// \f$ \left|a\right|^{2} \f$.
template<typename T, int N>
inline T StaticVector<T, N>::normSquared() const
{
    T product = 0;

    for(int i = 0; i < N; i++){
        product += value(i) * value(i);
    }

    return product;
}

/// Returns the length of the vector.
template<typename T, int N>
inline T StaticVector<T, N>::length() const
{
    return norm();
}

/// Returns the length squared of the vector.
template<typename T, int N>
inline T StaticVector<T, N>::lengthSquared() const
{
    return normSquared();
}

/// Normalizes the vector.
template<typename T, int N>
inline void StaticVector<T, N>::normalize()
{
    scale(1 / norm());
}

/// Returns a normalized version of the vector.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::normalized() const
{
    StaticVector<T, N> vector = *this;
    vector.normalize();
    return vector;
}

/// Scales the vector by the scalar \p value.
template<typename T, int N>
inline void StaticVector<T, N>::scale(T scalar)
{
    for(int i = 0; i < N; i++){
        value(i) *= scalar;
    }
}

/// Returns the vector scaled by \p value.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::scaled(T scalar) const
{
    StaticVector<T, N> vector = *this;
    vector.scale(scalar);
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
    return value(index);
}

/// Returns the value of the component at \p index.
///
/// \see value()
template<typename T, int N>
inline T& StaticVector<T, N>::operator[](int index)
{
    return value(index);
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
    return scaled(-1);
}

/// Returns the difference between the vector and \p vector.
///
/// \see subtract()
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator-(const StaticVector<T, N> &vector) const
{
    return subtract(vector);
}

/// Returns the vector multiplied by \p scalar.
///
/// \see scaled()
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator*(T scalar) const
{
    return scaled(scalar);
}

/// Returns the vector divided by \p scalar.
template<typename T, int N>
inline StaticVector<T, N> StaticVector<T, N>::operator/(T scalar) const
{
    return scaled(1 / scalar);
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
        if(!qFuzzyCompare(value(i), vector[i])){
            return false;
        }
    }

    return true;
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
