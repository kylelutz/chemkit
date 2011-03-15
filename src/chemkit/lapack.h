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
