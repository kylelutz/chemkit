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
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#ifndef _H_MSKIT_SPHERE
#define _H_MSKIT_SPHERE

#include "Vector.h"
#include "MSKContext.h"

typedef struct SphereRec {
  float *area;
  Vector3f *dot;
  int *StripLen;
  int *Sequence;
  int NStrip, NVertTot;
  int nDot;
  int *Tri;
  int NTri;
  int *Mesh;
  int NMesh;
} SphereRec, *SphereRecPtr;

struct _CSphere {
  SphereRecPtr Sphere[5];
  SphereRec *Array;
};

void SphereInit(MSKContext * G);
void SphereFree(MSKContext * G);

#endif
