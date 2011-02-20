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

#ifndef CHEMKIT_STATICVECTOR_H
#define CHEMKIT_STATICVECTOR_H

#include "chemkit.h"

#include "blas.h"
#include "lapack.h"
#include "constants.h"
#include "commainitializer.h"

namespace chemkit {

template<typename T, int N>
class StaticVector
{
    public:
        // construction and destruction
        StaticVector();
        StaticVector(const T *data, int size = N);
        StaticVector(const StaticVector<float, N> &vector);
        StaticVector(const StaticVector<double, N> &vector);

        // properties
        void setValue(int index, T value);
        T value(int index) const;
        T& value(int index);
        int size() const;
        bool isNull() const;
        T* data();
        const T* data() const;

        // math
        StaticVector<T, N> add(const StaticVector<T, N> &vector) const;
        StaticVector<T, N> subtract(const StaticVector<T, N> &vector) const;
        T dot(const StaticVector<T, N> &vector) const;
        StaticVector<T, N> cross(const StaticVector<T, N> &vector) const;
        T scalarTriple(const StaticVector<T, N> &a, const StaticVector<T, N> &b) const;
        StaticVector<T, N> vectorTriple(const StaticVector<T, N> &a, const StaticVector<T, N> &b) const;
        T norm() const;
        T normSquared() const;
        T length() const;
        T lengthSquared() const;
        void normalize();
        StaticVector<T, N> normalized() const;
        void scale(T scalar);
        StaticVector<T, N> scaled(T scalar) const;
        T angle(const StaticVector<T, N> &vector) const;
        T angleRadians(const StaticVector<T, N> &vector) const;

        // operators
        T operator[](int index) const;
        T& operator[](int index);
        StaticVector<T, N> operator+(const StaticVector<T, N> &vector) const;
        StaticVector<T, N> operator-() const;
        StaticVector<T, N> operator-(const StaticVector<T, N> &vector) const;
        StaticVector<T, N> operator*(T scalar) const;
        StaticVector<T, N> operator/(T scalar) const;
        StaticVector<T, N>& operator+=(const StaticVector<T, N> &vector);
        StaticVector<T, N>& operator-=(const StaticVector<T, N> &vector);
        StaticVector<T, N>& operator*=(T scalar);
        StaticVector<T, N>& operator/=(T scalar);
        bool operator==(const StaticVector<T, N> &vector) const;
        CommaInitializer<T> operator<<(const T value);

    private:
        T m_data[N];
};

} // end chemkit namespace

#include "staticvector-inline.h"

#endif // CHEMKIT_STATICVECTOR_H
