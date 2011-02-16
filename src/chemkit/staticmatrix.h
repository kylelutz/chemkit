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

#ifndef CHEMKIT_STATICMATRIX_H
#define CHEMKIT_STATICMATRIX_H

#include "staticvector.h"

#include "blas.h"
#include "lapack.h"
#include "commainitializer.h"

namespace chemkit {

template<typename T, int R, int C>
class StaticMatrix
{
    public:
        // construction and destruction
        StaticMatrix();
        StaticMatrix(const T *data);

        // properties
        void setValue(int row, int column, T value);
        T value(int row, int column) const;
        T& value(int row, int column);
        int rowCount() const;
        int columnCount() const;
        T* data();
        const T* data() const;
        void fill(T value);

        // math
        T trace() const;
        T determinant() const;
        void invert();
        StaticMatrix<T, R, C> inverted() const;
        StaticMatrix<T, R, R> multiply(const StaticMatrix<T, C, R> &matrix) const;

        // operators
        T operator()(int row, int column) const;
        T& operator()(int row, int column);
        CommaInitializer<T> operator=(const T value);
        bool operator==(const StaticMatrix<T, R, C> &matrix);

        // static methods
        static StaticMatrix<T, R, C> identity();

    private:
        T m_data[R*C];
};

template<typename T, int N>
class StaticMatrix<T, N, N>
{
    public:
        // construction and destruction
        StaticMatrix();
        StaticMatrix(const T *data);

        // properties
        void setValue(int row, int column, T value);
        T value(int row, int column) const;
        T& value(int row, int column);
        int rowCount() const;
        int columnCount() const;
        T* data();
        const T* data() const;
        void fill(T value);

        // math
        T trace() const;
        T determinant() const;
        void invert();
        StaticMatrix<T, N, N> inverted() const;
        StaticVector<T, N> multiply(const StaticVector<T, N> &vector) const;
        StaticMatrix<T, N, N> multiply(const StaticMatrix<T, N, N> &matrix) const;

        // decompositions
        void svd(StaticMatrix<T, N, N> *u, StaticVector<T, N> *s, StaticMatrix<T, N, N> *v) const;

        // operators
        T operator()(int row, int column) const;
        T& operator()(int row, int column);
        StaticMatrix<T, N, N> operator*(const StaticMatrix<T, N, N> &matrix) const;
        StaticMatrix<T, N, N>& operator*=(const StaticMatrix<T, N, N> &matrix);
        CommaInitializer<T> operator=(const T value);
        bool operator==(const StaticMatrix<T, N, N> &matrix);

        // static methods
        static StaticMatrix<T, N, N> identity();

    private:
        T m_data[N*N];
};

} // end chemkit namespace

#include "staticmatrix-inline.h"

#endif // CHEMKIT_STATICMATRIX_H
