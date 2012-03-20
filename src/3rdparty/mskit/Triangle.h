/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* Jacques Leroy (matrix inversion)
-* Thomas Malik (matrix multiplication)
-* Whoever wrote EISPACK
Z* -------------------------------------------------------------------
*/

#ifndef _H_MSKIT_TRIANGLE
#define _H_MSKIT_TRIANGLE

#include "Vector.h"
#include "MSKContext.h"

int *TrianglePointsToSurface(MSKContext * G, float *v, float *vn, int n,
                             float cutoff, int *nTriPtr, int **stripPtr, float *extent,
                             int cavity_mode);

int TriangleDegenerate(float *v1, float *n1, float *v2, float *n2, float *v3, float *n3);

#endif
