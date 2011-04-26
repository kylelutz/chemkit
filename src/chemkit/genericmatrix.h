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
