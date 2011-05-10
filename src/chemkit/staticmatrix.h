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

#ifndef CHEMKIT_STATICMATRIX_H
#define CHEMKIT_STATICMATRIX_H

#include "staticvector.h"

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
        StaticMatrix<T, R, C> add(const StaticMatrix<T, R, C> &matrix) const;
        StaticMatrix<T, R, C> subtract(const StaticMatrix<T, R, C> &matrix) const;
        StaticMatrix<T, R, C> multiply(T scalar) const;
        StaticMatrix<T, R, R> multiply(const StaticMatrix<T, C, R> &matrix) const;

        // operators
        T operator()(int row, int column) const;
        T& operator()(int row, int column);
        StaticMatrix<T, R, C> operator+(const StaticMatrix<T, R, C> &matrix) const;
        StaticMatrix<T, R, C> operator-(const StaticMatrix<T, R, C> &matrix) const;
        StaticMatrix<T, R, C> operator*(T scalar) const;
        CommaInitializer<T> operator<<(const T value);
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
        StaticMatrix<T, N, N> add(const StaticMatrix<T, N, N> &matrix) const;
        StaticMatrix<T, N, N> subtract(const StaticMatrix<T, N, N> &matrix) const;
        StaticVector<T, N> multiply(const StaticVector<T, N> &vector) const;
        StaticMatrix<T, N, N> multiply(T scalar) const;
        StaticMatrix<T, N, N> multiply(const StaticMatrix<T, N, N> &matrix) const;

        // decompositions
        void svd(StaticMatrix<T, N, N> *u, StaticVector<T, N> *s, StaticMatrix<T, N, N> *v) const;

        // operators
        T operator()(int row, int column) const;
        T& operator()(int row, int column);
        StaticMatrix<T, N, N> operator+(const StaticMatrix<T, N, N> &matrix) const;
        StaticMatrix<T, N, N> operator-(const StaticMatrix<T, N, N> &matrix) const;
        StaticMatrix<T, N, N> operator*(T scalar) const;
        StaticMatrix<T, N, N> operator*(const StaticMatrix<T, N, N> &matrix) const;
        StaticMatrix<T, N, N>& operator*=(const StaticMatrix<T, N, N> &matrix);
        CommaInitializer<T> operator<<(const T value);
        bool operator==(const StaticMatrix<T, N, N> &matrix);

        // static methods
        static StaticMatrix<T, N, N> identity();

    private:
        T m_data[N*N];
};

} // end chemkit namespace

#include "staticmatrix-inline.h"

#endif // CHEMKIT_STATICMATRIX_H
