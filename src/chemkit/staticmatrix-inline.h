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

#ifndef CHEMKIT_STATICMATRIX_INLINE_H
#define CHEMKIT_STATICMATRIX_INLINE_H

#include "staticmatrix.h"

#include <cmath>
#include <limits>
#include <cstdlib>

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
/// \see GenericMatrix

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
    StaticMatrix<T, R, R> product;

    chemkit::blas::gemm(R, R, C, m_data, false, matrix.data(), false, product.data());

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
    // temporary copy of matrix
    StaticMatrix<T, N, N> matrix = *this;

    // lu decomposition
    int info = 0;
    int ipiv[N];
    chemkit::lapack::getrf(matrix.data(), N, N, ipiv, &info);

    // if info is > 0, the matrix is singular
    // and the determinant is 0
    if(info > 0){
        return 0;
    }

    // product of diagonal values
    T determinant = 1;
    for(int i = 0; i < N; i++){
        determinant *= matrix(i, i);

        if(ipiv[i] != i+1)
            determinant = -determinant;
    }

    return determinant;
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
    // compute LU factorization
    int info = 0;
    int ipiv[N];
    chemkit::lapack::getrf(m_data, N, N, ipiv, &info);

    // if info is > 0 the matrix is singular and
    // cannot be inverted
    if(info > 0){
        return;
    }

    // compute inverse
    T work[N*N];
    int lwork = N*N;
    chemkit::lapack::getri(m_data, N, ipiv, work, lwork, &info);
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
    StaticVector<T, N> product;
    chemkit::blas::gemv(m_data, N, N, vector.data(), product.data());
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
    StaticMatrix<T, N, N> product;
    chemkit::blas::gemm(N, N, N, m_data, false, matrix.data(), false, product.data());
    return product;
}

// --- Decompositions ------------------------------------------------------ //
/// Calculates the singular value decomposition of the matrix.
template<typename T, int N>
inline void StaticMatrix<T, N, N>::svd(StaticMatrix<T, N, N> *u, StaticVector<T, N> *s, StaticMatrix<T, N, N> *v) const
{
    // copy of matrix
    StaticMatrix<T, N, N> matrix = *this;

    int lwork = 32;
    T work[32];
    int info;

    chemkit::lapack::gesvd(matrix.data(),
                           3,
                           3,
                           u->data(),
                           v->data(),
                           s->data(),
                           work,
                           lwork,
                           &info);
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
