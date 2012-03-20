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

#ifndef _H_MSKIT_SURFACEJOB
#define _H_MSKIT_SURFACEJOB

#include "MSKContext.h"

#ifdef NT
#undef NT
#endif

typedef struct {
  float vdw;
  int flags;
} SurfaceJobAtomInfo;

typedef struct {
  /* input */
  float *coord;
  SurfaceJobAtomInfo *atomInfo;

  float maxVdw;

  int nPresent;
  int *presentVla;

  int solventSphereIndex, sphereIndex;
  int circumscribe;

  int surfaceQuality;
  int surfaceType;
  int surfaceSolvent;

  float probeRadius;
  float pointSep;
  float trimCutoff;
  float trimFactor;

  int cavityCull;
  int cavityMode;
  float cavityRadius;
  float cavityCutoff;

  /* results */
  int N, NT;
  float *V, *VN;
  int *T, *S;

  int oneColorFlag, oneAlphaFlag;
  int oneColor, *VC;
  float oneAlpha, *VA;

} SurfaceJob;


SurfaceJob *SurfaceJobNew(MSKContext * G,
	                         float *coord, SurfaceJobAtomInfo * atom_info,
	                         float max_vdw, float probe_radius,
	                         int surface_quality, int surface_type,
	                         int surface_solvent, int cavity_cull,
	                         int cavity_mode, float cavity_radius,
	                         float cavity_cutoff, float trim_cutoff, float trim_factor);

void SurfaceJobFree(MSKContext * G, SurfaceJob * I);

int SurfaceJobRun(MSKContext * G, SurfaceJob * I);

void SurfaceJobColoring(MSKContext *G, SurfaceJob * I, const int *colors, const float *transp);

void SurfaceJobPurgeResult(MSKContext * G, SurfaceJob * I);

#endif
