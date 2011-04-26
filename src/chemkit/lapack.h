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

#ifndef CHEMKIT_LAPACK_H
#define CHEMKIT_LAPACK_H

extern "C" {

// SGETRF(M, N, A, LDA, IPIV, INFO)
void sgetrf_(int *m, int *n, float *A, int *lda, int *ipiv, int *info);

// DGETRF(M, N, A, LDA, IPIV, INFO)
void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);

// SGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
void sgetri_(int *n, float *A, int *lda, int *ipiv, float *work, int *lwork,
             int *info);

// DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
void dgetri_(int *n, double *A, int *lda, int *ipiv, double *work, int *lwork,
             int *info);

// SGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
void sgesvd_(char *jobu, char *jobvt, int *M, int *N,
             float *A, int *lda, float *S, float *U, int *ldu,
             float *VT, int *ldvt, float *work, int *lwork,
             int *info);

// DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
void dgesvd_(char *jobu, char *jobvt, int *M, int *N,
             double *A, int *lda, double *S, double *U, int *ldu,
             double *VT, int *ldvt, double *work, int *lwork,
             int *info);

} // end extern "C"

namespace chemkit {

/// \ingroup chemkit
/// \brief The %chemkit::%lapack namespace contains wrappers around
///        the LAPACK library functions.
///
/// The Linear Algebra PACKage (LAPACK) library contains functions
/// for performing linear algebra and is built on top of BLAS (see
/// \ref chemkit::blas). The functions contained in this namespace
/// provide wrappers around the LAPACK functions which make their
/// usage from C++ more natural.
///
/// The functions contained in this namespace provide low-level
/// primitives used for linear algebra. The following classes make
/// use of these primitives and provide an easier to use interface:
///     - StaticMatrix
///     - GenericMatrix
namespace lapack {

/// Compute the LU decomposition of a general matrix.
inline void getrf(float *A, int rowCount, int columnCount, int *ipiv,
                  int *info)
{
    sgetrf_(&rowCount, &columnCount, A, &rowCount, ipiv, info);
}

/// \overload
inline void getrf(double *A, int rowCount, int columnCount, int *ipiv,
                  int *info)
{
    dgetrf_(&rowCount, &columnCount, A, &rowCount, ipiv, info);
}

/// Compute inverse of general matrix from its LU decomposition.
inline void getri(float *A, int size, int *ipiv, float *work, int lwork,
                  int *info)
{
    sgetri_(&size, A, &size, ipiv, work, &lwork, info);
}

/// \overload
inline void getri(double *A, int size, int *ipiv, double *work, int lwork,
                  int *info)
{
    dgetri_(&size, A, &size, ipiv, work, &lwork, info);
}

/// Compute the singular value decomposition of a general matrix.
inline void gesvd(float *A, int rowCount, int columnCount, float *U,
                  float *Vt, float *S, float *work, int lwork, int *info)
{
    char jobu = 'A';
    char jobvt = 'A';

    sgesvd_(&jobu,
            &jobvt,
            &rowCount,
            &columnCount,
            A,
            &rowCount,
            S,
            U,
            &rowCount,
            Vt,
            &columnCount,
            work,
            &lwork,
            info);
}

/// \overload
inline void gesvd(double *A, int rowCount, int columnCount, double *U,
                  double *Vt, double *S, double *work, int lwork, int *info)
{
    char jobu = 'A';
    char jobvt = 'A';

    dgesvd_(&jobu,
            &jobvt,
            &rowCount,
            &columnCount,
            A,
            &rowCount,
            S,
            U,
            &rowCount,
            Vt,
            &columnCount,
            work,
            &lwork,
            info);
}

} // end lapack namespace

} // end chemkit namespace

#endif // CHEMKIT_LAPACK_H
