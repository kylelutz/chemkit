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

#ifndef CHEMKIT_STATICMATRIX_INLINE_H
#define CHEMKIT_STATICMATRIX_INLINE_H

#include "staticmatrix.h"

#include <cmath>
#include <limits>
#include <cstdlib>

#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/Eigen>

namespace chemkit {

// === StaticMatrix ======================================================== //
/// \class StaticMatrix staticmatrix.h chemkit/staticmatrix.h
/// \ingroup chemkit
/// \brief The StaticMatrix template class implements a fixed-size
///        matrix.
///
/// The StaticMatrix template has three parameters:
///     - \b T: The data type.
///     - \b R: The number of rows in the matrix.
///     - \b C: The number of columns in the matrix.
///
/// The StaticMatrix class supports initialization via the comma
/// operator:
/// \code
/// StaticMatrix<int, 2, 4> matrix;
///
/// matrix << 1, 2, 3, 4,
///           5, 6, 7, 8;
/// \endcode

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty matrix.
///
/// The following code:
/// \code
/// StaticMatrix<int, 2, 4> matrix;
/// \endcode
/// Will create the following matrix:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       0 & 0 & 0 & 0 \\
///       0 & 0 & 0 & 0 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
template<typename T, int R, int C>
inline StaticMatrix<T, R, C>::StaticMatrix()
{
    fill(0);
}

// --- Properties ---------------------------------------------------------- //
/// Sets the value at \p row and \p column to \p value.
template<typename T, int R, int C>
inline void StaticMatrix<T, R, C>::setValue(int row, int column, T value)
{
    m_data[column*R + row] = value;
}

/// Returns the value at \p row and \p column in the matrix.
template<typename T, int R, int C>
inline T StaticMatrix<T, R, C>::value(int row, int column) const
{
    return m_data[column*R + row];
}

/// Returns a reference to the value at \p row and \p column in the
/// matrix.
template<typename T, int R, int C>
inline T& StaticMatrix<T, R, C>::value(int row, int column)
{
    return m_data[column*R + row];
}

/// Returns the number of rows in the matrix.
///
/// The number of rows is equal to the parameter \b R.
template<typename T, int R, int C>
inline int StaticMatrix<T, R, C>::rowCount() const
{
    return R;
}

/// Returns the number of columns in the matrix.
///
/// The number of columns is equivalent to the parameter \b C.
template<typename T, int R, int C>
inline int StaticMatrix<T, R, C>::columnCount() const
{
    return C;
}

/// Returns the data for the matrix. The data is in column-major
/// order.
template<typename T, int R, int C>
inline T* StaticMatrix<T, R, C>::data()
{
    return m_data;
}

/// \overload
template<typename T, int R, int C>
inline const T* StaticMatrix<T, R, C>::data() const
{
    return m_data;
}

/// Fills the matrix with \p value.
template<typename T, int R, int C>
inline void StaticMatrix<T, R, C>::fill(T value)
{
    for(int i = 0; i < R*C; i++){
        m_data[i] = value;
    }
}

// --- Math --------------------------------------------------------------------- //
/// Returns the trace of the matrix.
///
/// This method is only applicable to sqaure matrices. For
/// non-square matrices \c 0 is returned.
template<typename T, int R, int C>
inline T StaticMatrix<T, R, C>::trace() const
{
    return 0;
}

/// Returns the determinant of the matrix.
///
/// This method is only applicable to sqaure matrices. For
/// non-square matrices \c 0 is returned.
template<typename T, int R, int C>
inline T StaticMatrix<T, R, C>::determinant() const
{
    return 0;
}

/// Inverts the matrix in-place.
///
/// This method is only applicable to sqaure matrices. For
/// non-square matrices this method does nothing.
template<typename T, int R, int C>
inline void StaticMatrix<T, R, C>::invert()
{
}

/// Returns the inverse of the matrix.
///
/// This method is only applicable to sqaure matrices. For
/// non-square matrices this method does nothing and the original
/// matrix is returned.
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::inverted() const
{
    StaticMatrix<T, R, C> matrix = *this;

    return matrix;
}

/// Returns the sum of the matrix and \p matrix.
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::add(const StaticMatrix<T, R, C> &matrix) const
{
    StaticMatrix<T, R, C> sum;

    for(int i = 0; i < R * C; i++){
        sum.m_data[i] = m_data[i] + matrix.m_data[i];
    }

    return sum;
}

/// Returns the difference of the matrix and \p matrix.
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::subtract(const StaticMatrix<T, R, C> &matrix) const
{
    StaticMatrix<T, R, C> difference;

    for(int i = 0; i < R * C; i++){
        difference.m_data[i] = m_data[i] - matrix.m_data[i];
    }

    return difference;
}

/// Multiplies values in the matrix by \p scalar.
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::multiply(T scalar) const
{
    StaticMatrix<T, R, C> product;

    for(int i = 0; i < R * C; i++){
        product.m_data[i] = scalar * m_data[i];
    }

    return product;
}

/// Multiplies the matrix by \p matrix.
template<typename T, int R, int C>
inline StaticMatrix<T, R, R> StaticMatrix<T, R, C>::multiply(const StaticMatrix<T, C, R> &matrix) const
{
    Eigen::Matrix<T, R, C> a;
    for(int i = 0; i < R; i++){
        for(int j = 0; j < C; j++){
            a(i, j) = value(i, j);
        }
    }

    Eigen::Matrix<T, C, R> b;
    for(int i = 0; i < C; i++){
        for(int j = 0; j < R; j++){
            b(i, j) = matrix(i, j);
        }
    }

    Eigen::Matrix<T, R, R> c = a * b;

    StaticMatrix<T, R, R> product;
    for(int i = 0; i < R; i++){
        for(int j = 0; j < R; j++){
            product(i, j) = c(i, j);
        }
    }

    return product;
}

// --- Operators ---------------------------------------------------------------- //
/// Returns the value at \p row and \p column in the matrix.
///
/// \see value()
template<typename T, int R, int C>
inline T StaticMatrix<T, R, C>::operator()(int row, int column) const
{
    return value(row, column);
}

/// Returns a reference to the value at \p row and \p column in the
/// matrix.
///
/// \see value()
template<typename T, int R, int C>
inline T& StaticMatrix<T, R, C>::operator()(int row, int column)
{
    return value(row, column);
}

/// Returns the sum of the matrix and \p matrix.
///
/// \see add()
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::operator+(const StaticMatrix<T, R, C> &matrix) const
{
    return add(matrix);
}

/// Returns the difference of the matrix and \p matrix.
///
/// \see subtract()
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::operator-(const StaticMatrix<T, R, C> &matrix) const
{
    return subtract(matrix);
}

/// Multiplies the values in the matrix by \p scalar.
///
/// \see multiply()
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::operator*(T scalar) const
{
    return multiply(scalar);
}

template<typename T, int R, int C>
inline CommaInitializer<T> StaticMatrix<T, R, C>::operator<<(const T value)
{
    m_data[0] = value;

    return CommaInitializer<T>(m_data, R, C);
}

/// Returns \c true if the matrix is equal to \p matrix.
template<typename T, int R, int C>
inline bool StaticMatrix<T, R, C>::operator==(const StaticMatrix<T, R, C> &matrix)
{
    for(int i = 0; i < R*C; i++){
        if(std::abs(m_data[i] - matrix.m_data[i]) > std::numeric_limits<T>::epsilon()){
            return false;
        }
    }

    return true;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the identity matrix.
///
/// The following code:
/// \code
/// StaticMatrix<int, 4, 4>::identity()
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
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> StaticMatrix<T, R, C>::identity()
{
    StaticMatrix<T, R, C> matrix;

    for(int i = 0; i < C; i++){
        matrix(i, i) = 1;
    }

    return matrix;
}

// --- Related Functions --------------------------------------------------- //
/// Multiplies the matrix by \p scalar.
///
/// \related StaticMatrix
template<typename T, int R, int C>
inline StaticMatrix<T, R, C> operator*(T scalar, const StaticMatrix<T, R, C> &matrix)
{
    return matrix * scalar;
}

// === StaticMatrix<T, N, N> =============================================== //
// --- Construction and Destruction ---------------------------------------- //
template<typename T, int N>
inline StaticMatrix<T, N, N>::StaticMatrix()
{
    fill(0);
}

template<typename T, int N>
inline StaticMatrix<T, N, N>::StaticMatrix(const T *data)
{
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            m_data[j*N+i] = data[i*N+j];
        }
    }
}

// --- Properties ---------------------------------------------------------- //
template<typename T, int N>
inline void StaticMatrix<T, N, N>::setValue(int row, int column, T value)
{
    m_data[column*N + row] = value;
}

template<typename T, int N>
inline T StaticMatrix<T, N, N>::value(int row, int column) const
{
    return m_data[column*N + row];
}

template<typename T, int N>
inline T& StaticMatrix<T, N, N>::value(int row, int column)
{
    return m_data[column*N + row];
}

template<typename T, int N>
inline int StaticMatrix<T, N, N>::rowCount() const
{
    return N;
}

template<typename T, int N>
inline int StaticMatrix<T, N, N>::columnCount() const
{
    return N;
}

template<typename T, int N>
inline T* StaticMatrix<T, N, N>::data()
{
    return m_data;
}

template<typename T, int N>
inline const T* StaticMatrix<T, N, N>::data() const
{
    return m_data;
}

template<typename T, int N>
inline void StaticMatrix<T, N, N>::fill(T value)
{
    for(int i = 0; i < N*N; i++){
        m_data[i] = value;
    }
}

// --- Math ---------------------------------------------------------------- //
template<typename T, int N>
inline T StaticMatrix<T, N, N>::trace() const
{
    T trace = 0;

    for(int i = 0; i < N; i++){
        trace += value(i, i);
    }

    return trace;
}

template<typename T, int N>
inline T StaticMatrix<T, N, N>::determinant() const
{
    Eigen::Matrix<T, N, N> matrix;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            matrix(i, j) = value(i, j);
        }
    }

    return matrix.determinant();
}

template<>
inline double StaticMatrix<double, 1, 1>::determinant() const
{
    return value(0, 0);
}

template<>
inline double StaticMatrix<double, 2, 2>::determinant() const
{
    const StaticMatrix<double, 2, 2> &m = *this;

    return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
}

template<>
inline double StaticMatrix<double, 3, 3>::determinant() const
{
    const StaticMatrix<double, 3, 3> &m = *this;

    return m(0, 0) * m(1, 1) * m(2, 2) +
           m(0, 1) * m(1, 2) * m(2, 0) +
           m(0, 2) * m(1, 0) * m(2, 1) -
           m(0, 2) * m(1, 1) * m(2, 0) -
           m(0, 1) * m(1, 0) * m(2, 2) -
           m(0, 0) * m(1, 2) * m(2, 1);
}

template<>
inline double StaticMatrix<double, 4, 4>::determinant() const
{
    const StaticMatrix<double, 4, 4> &m = *this;

    // minor 0, 0
    double m0 = m(1, 1) * m(2, 2) * m(3, 3) +
                m(1, 2) * m(2, 3) * m(3, 1) +
                m(1, 3) * m(2, 1) * m(3, 2) -
                m(1, 3) * m(2, 2) * m(3, 1) -
                m(1, 2) * m(2, 1) * m(3, 3) -
                m(1, 1) * m(2, 3) * m(3, 2);

    // minor 0, 1
    double m1 = m(1, 0) * m(2, 2) * m(3, 3) +
                m(1, 2) * m(2, 3) * m(3, 0) +
                m(1, 3) * m(2, 0) * m(3, 2) -
                m(1, 3) * m(2, 2) * m(3, 0) -
                m(1, 2) * m(2, 0) * m(3, 3) -
                m(1, 0) * m(2, 3) * m(3, 2);

    // minor 0, 2
    double m2 = m(1, 0) * m(2, 1) * m(3, 3) +
                m(1, 1) * m(2, 3) * m(3, 0) +
                m(1, 3) * m(2, 0) * m(3, 1) -
                m(1, 3) * m(2, 1) * m(3, 0) -
                m(1, 1) * m(2, 0) * m(3, 3) -
                m(1, 0) * m(2, 3) * m(3, 1);

    // minor 0,3
    double m3 = m(1, 0) * m(2, 1) * m(3, 2) +
                m(1, 1) * m(2, 2) * m(3, 0) +
                m(1, 2) * m(2, 0) * m(3, 1) -
                m(1, 2) * m(2, 1) * m(3, 0) -
                m(1, 1) * m(2, 0) * m(3, 2) -
                m(1, 0) * m(2, 2) * m(3, 1);

    return m(0, 0) * m0 - m(0, 1) * m1 + m(0, 2) * m2 - m(0, 3) * m3;
}

template<typename T, int N>
inline void StaticMatrix<T, N, N>::invert()
{
    // copy matrix
    Eigen::Matrix<T, N, N> matrix;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            matrix(i, j) = value(i, j);
        }
    }

    // calculate inverse
    Eigen::Matrix<T, N, N> inverse = matrix.inverse();

    // copy inverse values to the matrix
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            setValue(i, j, inverse(i, j));
        }
    }
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::inverted() const
{
    StaticMatrix<T, N, N> matrix = *this;
    matrix.invert();
    return matrix;
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::add(const StaticMatrix<T, N, N> &matrix) const
{
    StaticMatrix<T, N, N> sum;

    for(int i = 0; i < N * N; i++){
        sum.m_data[i] = m_data[i] + matrix.m_data[i];
    }

    return sum;
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::subtract(const StaticMatrix<T, N, N> &matrix) const
{
    StaticMatrix<T, N, N> difference;

    for(int i = 0; i < N * N; i++){
        difference.m_data[i] = m_data[i] - matrix.m_data[i];
    }

    return difference;
}

template<typename T, int N>
inline StaticVector<T, N> StaticMatrix<T, N, N>::multiply(const StaticVector<T, N> &vector) const
{
    Eigen::Matrix<T, N, N> m;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            m(i, j) = value(i, j);
        }
    }

    Eigen::Matrix<T, N, 1> v;
    for(int i = 0; i < N; i++){
        v[i] = vector[i];
    }

    Eigen::Matrix<T, N, 1> p = m * v;

    StaticVector<T, N> product;
    for(int i = 0; i < N; i++){
        product[i] = p[i];
    }

    return product;
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::multiply(T scalar) const
{
    StaticMatrix<T, N, N> product;

    for(int i = 0; i < N * N; i++){
        product.m_data[i] = scalar * m_data[i];
    }

    return product;
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::multiply(const StaticMatrix<T, N, N> &matrix) const
{
    Eigen::Matrix<T, N, N> a;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            a(i, j) = value(i, j);
        }
    }

    Eigen::Matrix<T, N, N> b;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            b(i, j) = matrix(i, j);
        }
    }

    Eigen::Matrix<T, N, N> c = a * b;

    StaticMatrix<T, N, N> product;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            product(i, j) = c(i, j);
        }
    }

    return product;
}

// --- Decompositions ------------------------------------------------------ //
/// Calculates the singular value decomposition of the matrix.
template<typename T, int N>
inline void StaticMatrix<T, N, N>::svd(StaticMatrix<T, N, N> *u, StaticVector<T, N> *s, StaticMatrix<T, N, N> *v) const
{
    // copy of matrix
    Eigen::Matrix<T, N, N> matrix;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            matrix(i, j) = value(i, j);
        }
    }

    Eigen::JacobiSVD<Eigen::Matrix<T, N, N> > svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            (*u)(i, j) = svd.matrixU()(i, j);
        }
    }

    for(int i = 0; i < N; i++){
        (*s)[i] = svd.singularValues()[i];
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            (*v)(i, j) = svd.matrixV().transpose()(i, j);
        }
    }
}

// --- Operators ----------------------------------------------------------- //
template<typename T, int N>
inline T StaticMatrix<T, N, N>::operator()(int row, int column) const
{
    return value(row, column);
}

template<typename T, int N>
inline T& StaticMatrix<T, N, N>::operator()(int row, int column)
{
    return value(row, column);
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::operator+(const StaticMatrix<T, N, N> &matrix) const
{
    return add(matrix);
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::operator-(const StaticMatrix<T, N, N> &matrix) const
{
    return subtract(matrix);
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::operator*(T scalar) const
{
    return multiply(scalar);
}

template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::operator*(const StaticMatrix<T, N, N> &matrix) const
{
    return multiply(matrix);
}

template<typename T, int N>
inline StaticMatrix<T, N, N>& StaticMatrix<T, N, N>::operator*=(const StaticMatrix<T, N, N> &matrix)
{
    *this = *this * matrix;
    return *this;
}

template<typename T, int N>
inline CommaInitializer<T> StaticMatrix<T, N, N>::operator<<(const T value)
{
    m_data[0] = value;

    return CommaInitializer<T>(m_data, N, N);
}

template<typename T, int N>
inline bool StaticMatrix<T, N, N>::operator==(const StaticMatrix<T, N, N> &matrix)
{
    for(int i = 0; i < N*N; i++){
        if(std::abs(m_data[i] - matrix.m_data[i]) > std::numeric_limits<T>::epsilon()){
            return false;
        }
    }

    return true;
}

// --- Static Methods ------------------------------------------------------ //
template<typename T, int N>
inline StaticMatrix<T, N, N> StaticMatrix<T, N, N>::identity()
{
    StaticMatrix<T, N, N> matrix;

    for(int i = 0; i < N; i++){
        matrix(i, i) = 1;
    }

    return matrix;
}

// --- Related Functions --------------------------------------------------- //
template<typename T, int N>
inline StaticMatrix<T, N, N> operator*(T scalar, const StaticMatrix<T, N, N> &matrix)
{
    return matrix * scalar;
}

} // end chemkit namespace

#endif // CHEMKIT_STATICMATRIX_INLINE_H
