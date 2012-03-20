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

#include "os_std.h"

#include "Base.h"
#include "MemoryDebug.h"
#include "MemoryCache.h"
#include "Sphere.h"

#include "MSKContext.h"

MSKContext *MSKContextNew()
{
  MSKContext *G = Calloc(MSKContext, 1);
  
  if (G) {
    //ColorInit(G);
    SphereInit(G);

    // Map Settings
    SettingSet(G, cSetting_hash_max, 100);
    // Triangle Settings
    SettingSet(G, cSetting_triangle_max_passes, 5);
    SettingSet(G, cSetting_fit_tolerance, 0.0000001F);
    SettingSet(G, cSetting_fit_iterations, 1000);
    SettingSet(G, cSetting_fit_kabsch, 0);
    
    G->Ready = true;
  }
  
  return G;
}

void MSKContextFree(MSKContext * G)
{
  if (G) {
    SphereFree(G);
    //ColorFree(G);
    FreeP(G);
  }
}

void MSKContextClean(MSKContext * G)
{
  G->Ready = true;
  G->Interrupt = false;
  G->Stage = 0;
  G->Progress = 0;
}
