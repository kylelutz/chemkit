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

#ifndef _H_MSKIT_MSKCONTEXT
#define _H_MSKIT_MSKCONTEXT

typedef struct _CSphere CSphere;
//typedef struct _CColor CColor;

typedef struct _MSKContext {

  int Ready;                    /* is the program fully initialized and ready to receive
                                 * messages? */
  int Interrupt;                /* set when we are attempting to abort time-consuming calculations */
  
  CSphere *Sphere;

  // Progress
  int Stage;
  int Progress;
  
  struct {
    float cSetting_hash_max;
    int cSetting_triangle_max_passes;
    int cSetting_fit_iterations;
    int cSetting_fit_kabsch;
    double cSetting_fit_tolerance;
  } Settings;

} MSKContext;

MSKContext *MSKContextNew();
void MSKContextFree(MSKContext * G);
void MSKContextClean(MSKContext * G);

#define OrthoBusyStage(G, S)     ((G)->Stage = (S), (G)->Progress = 0)
#define OrthoBusyFast(G, A, E)   ((G)->Progress = (A) * 100 / (E))
#define OrthoStage(G)            ((G)->Stage)
#define OrthoFast(G)             ((G)->Progress)

#define SettingSet(G, K, V)      ((G)->Settings.K = (V))
#define SettingGet(G, K)         ((G)->Settings.K)

#endif