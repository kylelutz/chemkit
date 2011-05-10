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
