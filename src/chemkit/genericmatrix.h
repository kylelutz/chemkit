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

#ifndef CHEMKIT_GENERICMATRIX_H
#define CHEMKIT_GENERICMATRIX_H

#include "chemkit.h"

#include "commainitializer.h"

namespace chemkit {

template<typename T>
class GenericMatrix
{
    public:
        // construction and destruction
        GenericMatrix();
        GenericMatrix(int rowCount, int columnCount);
        GenericMatrix(const GenericMatrix<T> &matrix);
        ~GenericMatrix();

        // properties
        void setValue(int row, int column, T value);
        T value(int row, int column) const;
        T& value(int row, int column);
        void setRowCount(int rowCount);
        int rowCount() const;
        void setColumnCount(int columnCount);
        int columnCount() const;
        void resize(int rowCount, int columnCount);
        int size() const;
        T* data();
        const T* data() const;
        void fill(T value);

        // math
        T trace() const;
        GenericMatrix<T> add(const GenericMatrix<T> &matrix) const;
        GenericMatrix<T> subtract(const GenericMatrix<T> &matrix) const;
        GenericMatrix<T> multiply(const T scalar) const;
        GenericMatrix<T> multiply(const GenericMatrix<T> &matrix) const;

        // operators
        T operator()(int row, int column) const;
        T& operator()(int row, int column);
        GenericMatrix<T> operator+(const GenericMatrix<T> &matrix) const;
        GenericMatrix<T> operator-(const GenericMatrix<T> &matrix) const;
        GenericMatrix<T> operator*(const T scalar) const;
        GenericMatrix<T> operator*(const GenericMatrix<T> &matrix) const;
        GenericMatrix<T>& operator=(const GenericMatrix<T> &matrix);
        CommaInitializer<T> operator=(const T value);
        bool operator==(const GenericMatrix<T> &matrix) const;

        // static methods
        static GenericMatrix<T> ones(int rowCount, int columnCount);
        static GenericMatrix<T> zeros(int rowCount, int columnCount);
        static GenericMatrix<T> random(int rowCount, int columnCount);
        static GenericMatrix<T> identity(int rowCount, int columnCount);

    private:
        T* m_data;
        int m_rowCount;
        int m_columnCount;
};

} // end chemkit namespace

#include "genericmatrix-inline.h"

#endif // CHEMKIT_GENERICMATRIX_H
