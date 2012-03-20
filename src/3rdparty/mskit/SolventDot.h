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

#ifndef _H_MSKIT_SOLVENTDOT
#define _H_MSKIT_SOLVENTDOT

#include "MSKContext.h"

typedef struct {
  int nDot;
  float *dot;
  float *dotNormal;
  int *dotCode;
} SolventDot;

SolventDot *SolventDotNew(MSKContext * G,
	                         float *coord,
	                         SurfaceJobAtomInfo * atom_info,
	                         float probe_radius, SphereRec * sp,
	                         int *present,
	                         int circumscribe,
	                         int surface_solvent, int cavity_cull,
	                         float max_vdw,
	                         int cavity_mode, float cavity_radius, 
	                         float cavity_cutoff);

void SolventDotFree(SolventDot * I);

#endif
