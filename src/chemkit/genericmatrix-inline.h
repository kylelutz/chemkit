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

#ifndef CHEMKIT_GENERICMATRIX_INLINE_H
#define CHEMKIT_GENERICMATRIX_INLINE_H

#include "genericmatrix.h"

#include <cmath>
#include <limits>
#include <cstdlib>

#include "blas.h"
#include "lapack.h"

namespace chemkit {

// === GenericMatrix ======================================================= //
/// \class GenericMatrix genericmatrix.h chemkit/genericmatrix.h
/// \ingroup chemkit
/// \brief The GenericMatrix template class represents a resizable
///        matrix.
///
/// The GenericMatrix template has one parameter:
///     - \b T: The value type.
///
/// \see Matrix

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty matrix with \c 0 rows and \c 0 columns.
template<typename T>
inline GenericMatrix<T>::GenericMatrix()
{
    m_data = 0;
    m_rowCount = 0;
    m_columnCount = 0;
}

/// Creates a new, empty matrix with \p rowCount rows and
/// \p columnCount columns.
template<typename T>
inline GenericMatrix<T>::GenericMatrix(int rowCount, int columnCount)
{
    m_data = new T[rowCount * columnCount];
    m_rowCount = rowCount;
    m_columnCount = columnCount;

    fill(0);
}

/// Creates a new matrix as a copy of \p matrix.
template<typename T>
inline GenericMatrix<T>::GenericMatrix(const GenericMatrix<T> &matrix)
{
    m_data = new T[matrix.rowCount() * matrix.columnCount()];
    m_rowCount = matrix.rowCount();
    m_columnCount = matrix.columnCount();

    for(int i = 0; i < m_rowCount * m_columnCount; i++){
        m_data[i] = matrix.m_data[i];
    }
}

/// Destroys the matrix.
template<typename T>
inline GenericMatrix<T>::~GenericMatrix()
{
    delete[] m_data;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the value at \p row and \p column to \p value.
template<typename T>
inline void GenericMatrix<T>::setValue(int row, int column, T value)
{
    m_data[column*m_rowCount + row] = value;
}

/// Returns the value at \p row and \p column.
template<typename T>
inline T GenericMatrix<T>::value(int row, int column) const
{
    return m_data[column*m_rowCount + row];
}

/// Returns a reference to the value at \p row and \p column.
template<typename T>
inline T& GenericMatrix<T>::value(int row, int column)
{
    return m_data[column*m_rowCount + row];
}

/// Sets the number of rows to \p rowCount.
///
/// \see resize().
template<typename T>
inline void GenericMatrix<T>::setRowCount(int rowCount)
{
    resize(rowCount, m_columnCount);
}

/// Returns the number of rows in the matrix.
template<typename T>
inline int GenericMatrix<T>::rowCount() const
{
    return m_rowCount;
}

/// Sets the number of columns to \p columnCount.
///
/// \see resize().
template<typename T>
inline void GenericMatrix<T>::setColumnCount(int columnCount)
{
    resize(m_rowCount, columnCount);
}

/// Returns the number of columns in the matrix.
template<typename T>
inline int GenericMatrix<T>::columnCount() const
{
    return m_columnCount;
}

/// Resizes the matrix to contain \p rowCount number of rows and
/// \p columnCount number of columns.
template<typename T>
inline void GenericMatrix<T>::resize(int rowCount, int columnCount)
{
    // allocate new space
    T *data = new T[rowCount * columnCount];

    // copy old data
    for(int i = 0; i < qMin(rowCount, m_rowCount); i++){
        for(int j = 0; j < qMin(columnCount, m_columnCount); j++){
            data[j*rowCount+i] = value(i, j);
        }
    }

    // deallocate old space
    delete[] m_data;

    // set new data
    m_rowCount = rowCount;
    m_columnCount = columnCount;
    m_data = data;
}

/// Returns the size of the matrix. The size is equal to the product
/// of the number of rows and columns.
template<typename T>
inline int GenericMatrix<T>::size() const
{
    return m_rowCount * m_columnCount;
}

/// Returns the data in the matrix. The data is stored in
/// column-major order.
template<typename T>
inline T* GenericMatrix<T>::data()
{
    return m_data;
}

/// \overload
template<typename T>
inline const T* GenericMatrix<T>::data() const
{
    return m_data;
}

/// Fills the matrix with \p value.
template<typename T>
inline void GenericMatrix<T>::fill(const T value)
{
    for(int i = 0; i < m_rowCount * m_columnCount; i++){
        m_data[i] = value;
    }
}

// --- Math ---------------------------------------------------------------- //
/// Returns the trace of the matrix. The trace is the sum of the
/// values along the main diagonal.
///
/// The expression \c A.trace() will return the result of
/// \f$tr(A)\f$.
template<typename T>
inline T GenericMatrix<T>::trace() const
{
    int size = qMin(m_rowCount, m_columnCount);

    T trace = 0;

    for(int i = 0; i < size; i++){
        trace += value(i, i);
    }

    return trace;
}

/// Returns the sum of the matrix and \p matrix.
///
/// The expression \c A.add(B) will return the result of
/// \f$ A + B \f$.
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::add(const GenericMatrix<T> &matrix) const
{
    GenericMatrix<T> product(m_rowCount, m_columnCount);

    for(int i = 0; i < m_rowCount * m_columnCount; i++){
        product.m_data[i] = m_data[i] + matrix.m_data[i];
    }

    return product;
}

/// Returns the difference of the matrix and \p matrix.
///
/// The expression \c A.subtract(B) will return the result of
/// \f$ A - B \f$.
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::subtract(const GenericMatrix<T> &matrix) const
{
    GenericMatrix<T> product(m_rowCount, m_columnCount);

    for(int i = 0; i < m_rowCount * m_columnCount; i++){
        product.m_data[i] = m_data[i] - matrix.m_data[i];
    }

    return product;
}

/// Returns the product of the matrix and \p scalar.
///
/// The expression \c A.multiply(r) will return the result of
/// \f$ r \times A \f$.
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::multiply(T scalar) const
{
    GenericMatrix<T> product(m_rowCount, m_columnCount);

    for(int i = 0; i < size(); i++){
        product.m_data[i] = scalar * m_data[i];
    }

    return product;
}

/// Returns the product of the matrix and \p matrix.
///
/// The expression \c A.multiply(B) will return the result of
/// \f$ A \times B \f$.
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::multiply(const GenericMatrix<T> &matrix) const
{
    GenericMatrix<T> product(rowCount(), matrix.columnCount());

    chemkit::blas::gemm(rowCount(),
                        matrix.columnCount(),
                        columnCount(),
                        m_data,
                        false,
                        matrix.data(),
                        false,
                        product.data());

    return product;
}

// --- Operators ----------------------------------------------------------- //
/// Returns the value at \p row and \p column.
template<typename T>
inline T GenericMatrix<T>::operator()(int row, int column) const
{
    return value(row, column);
}

/// Returns a reference to the value at \p row and \p column.
template<typename T>
inline T& GenericMatrix<T>::operator()(int row, int column)
{
    return value(row, column);
}

/// Returns the sum of the matrix and \p matrix.
///
/// \see add()
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::operator+(const GenericMatrix<T> &matrix) const
{
    return add(matrix);
}

/// Returns the difference of the matrix and \p matrix.
///
/// \see subtract()
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::operator-(const GenericMatrix<T> &matrix) const
{
    return subtract(matrix);
}

/// Returns the product of the matrix and \p scalar.
///
/// \see multiply()
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::operator*(const T scalar) const
{
    return multiply(scalar);
}

/// Returns the product of the matrix and \p matrix.
///
/// \see multiply()
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::operator*(const GenericMatrix<T> &matrix) const
{
    return multiply(matrix);
}

/// Sets the matrix equal to \p matrix.
template<typename T>
inline GenericMatrix<T>& GenericMatrix<T>::operator=(const GenericMatrix<T> &matrix)
{
    delete[] m_data;

    m_rowCount = matrix.rowCount();
    m_columnCount = matrix.columnCount();
    m_data = new T[m_rowCount * m_columnCount];

    for(int i = 0; i < m_rowCount * m_columnCount; i++){
        m_data[i] = matrix.m_data[i];
    }

    return *this;
}

template<typename T>
inline CommaInitializer<T> GenericMatrix<T>::operator=(const T value)
{
    m_data[0] = value;

    return CommaInitializer<T>(m_data, m_rowCount, m_columnCount);
}

/// Returns \c true if the matrix equals \p matrix.
template<typename T>
inline bool GenericMatrix<T>::operator==(const GenericMatrix<T> &matrix) const
{
    if(rowCount() != matrix.rowCount()){
        return false;
    }
    else if(columnCount() != matrix.columnCount()){
        return false;
    }

    for(int i = 0; i < m_rowCount * m_columnCount; i++){
        if(std::abs(m_data[i] - matrix.m_data[i]) > std::numeric_limits<T>::epsilon()){
            return false;
        }
    }

    return true;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a matrix filled with \c 1.
///
/// For example, the following code:
/// \code
/// GenericMatrix<int> matrix = GenericMatrix<int>::ones(2, 3);
/// \endcode
/// Will return the following matrix:
/** \f[
///   \left[
///   {
///     \begin{array}{ccc}
///       1 & 1 & 1 \\
///       1 & 1 & 1 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::ones(int rowCount, int columnCount)
{
    GenericMatrix<T> matrix(rowCount, columnCount);
    fill(1);
    return matrix;
}

/// Returns a matrix filled with \c 0.
///
/// For example, the following code:
/// \code
/// GenericMatrix<int> matrix = GenericMatrix<int>::zeros(2, 3);
/// \endcode
/// Will return the following matrix:
/** \f[
///   \left[
///   {
///     \begin{array}{ccc}
///       0 & 0 & 0 \\
///       0 & 0 & 0 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::zeros(int rowCount, int columnCount)
{
    GenericMatrix<T> matrix(rowCount, columnCount);
    fill(0);
    return matrix;
}

/// Returns a matrix containing random values.
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::random(int rowCount, int columnCount)
{
    GenericMatrix<T> matrix(rowCount, columnCount);

    for(int i = 0; i < size(); i++){
        matrix.m_data[i] = qrand() - (RAND_MAX/2);
    }

    return matrix;
}

/// Returns a matrix with \c 1 on the main diagonal.
///
/// For example, the following code:
/// \code
/// GenericMatrix<int> matrix = GenericMatrix<int>::identity(4, 4);
/// \endcode
/// Will return the following matrix:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       1 & 0 & 0 & 0 \\
///       0 & 1 & 0 & 0 \\
///       0 & 0 & 1 & 0 \\
///       0 & 0 & 0 & 1 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
template<typename T>
inline GenericMatrix<T> GenericMatrix<T>::identity(int rowCount, int columnCount)
{
    GenericMatrix<T> matrix(rowCount, columnCount);

    int size = qMin(rowCount, columnCount);
    for(int i = 0; i < size; i++){
        matrix(i, i) = 1;
    }

    return matrix;
}

// --- Related Functions --------------------------------------------------- //
/// Returns the product of \p scalar and \p matrix.
///
/// \related GenericMatrix
template<typename T>
inline GenericMatrix<T> operator*(const T scalar, const GenericMatrix<T> &matrix)
{
    return matrix * scalar;
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICMATRIX_INLINE_H
