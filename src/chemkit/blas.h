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

#ifndef CHEMKIT_BLAS_H
#define CHEMKIT_BLAS_H

#include "chemkit.h"

extern "C" {

// SGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
void sgemv_(const char *trans, const int *m, const int *n, const float *alpha,
            const float *a, const int *lda, const float *x, const int *incx,
            const float *beta, float *y, const int *incy);

// DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
            const double *a, const int *lda, const double *x, const int *incx,
            const double *beta, double *y, const int *incy);

// SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
void sgemm_(const char *transA, const char *transB, const int *m,
            const int *n, const int *k, const float *alpha, const float *a,
            const int *lda, const float *b, const int *ldb,
            const float *beta, float *c, const int *ldc);

// DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
void dgemm_(const char *transA, const char *transB, const int *m,
            const int *n, const int *k, const double *alpha, const double *a,
            const int *lda, const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);

} // end extern "C"

namespace chemkit {

/// \ingroup chemkit
/// \brief The %chemkit::%blas namespace contains wrappers around the
///        BLAS library functions.
///
/// The Basic Linear Algebra Subprograms (BLAS) library contains
/// functions for performing basic linear algebra operations such
/// vector and matrix multiplication. The functions contained in this
/// namespace provide wrappers around the BLAS functions which make
/// their usage from C++ more natural.
///
/// The functions contained in this namespace provide low-level
/// primitives used for linear algebra. The following classes make
/// use of these primitives and provide an easier to use interface:
///     - Point
///     - Vector
///     - StaticVector
///     - StaticMatrix
///     - GenericMatrix
namespace blas {

/// Multiply the matrix M by vector V and place the result in P.
inline void gemv(const float *M, int rowCount, int columnCount,
                 const float *V, float *P)
{
    char trans = 'N';
    float alpha = 1;
    float beta = 0;
    int inc = 1;

    sgemv_(&trans, &rowCount, &columnCount, &alpha, M, &rowCount, V, &inc,
           &beta, P, &inc);
}

/// \overload
inline void gemv(const double *M, int rowCount, int columnCount,
                 const double *V, double *P)
{
    char trans = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    int inc = 1;

    dgemv_(&trans, &rowCount, &columnCount, &alpha, M, &rowCount, V, &inc,
           &beta, P, &inc);
}

/// Multiply matrix A by B and place the result in C.
inline void gemm(int m, int n, int k, const float *A, bool transA,
                 const float *B, bool transB, float *C)
{
    char tA = transA ? 'T' : 'N';
    char tB = transB ? 'T' : 'N';
    float alpha = 1.0;
    float beta = 0.0;

    sgemm_(&tA, // transA
           &tB, // transB
           &m, // m
           &n, // n
           &k, // k
           &alpha, // alpha
           A, // A
           transA ? &k : &m, // lda
           B, // B
           transB ? &n : &k, // ldb
           &beta, // beta
           C, // C
           &m); // ldc
}

/// \overload
inline void gemm(int m, int n, int k, const double *A, bool transA,
                 const double *B, bool transB, double *C)
{
    char tA = transA ? 'T' : 'N';
    char tB = transB ? 'T' : 'N';
    double alpha = 1.0;
    double beta = 0.0;

    dgemm_(&tA, // transA
           &tB, // transB
           &m, // m
           &n, // n
           &k, // k
           &alpha, // alpha
           A, // A
           transA ? &k : &m, // lda
           B, // B
           transB ? &n : &k, // ldb
           &beta, // beta
           C, // C
           &m); // ldc
}

} // end blas namespace

} // end chemkit namespace

#endif // CHEMKIT_BLAS_H
