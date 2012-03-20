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

#ifndef _H_MSKIT_ERR
#define _H_MSKIT_ERR

#include "MSKContext.h"

void ErrFatal(MSKContext * G, const char *where, const char *what);
void ErrPointer(MSKContext * G, const char *file, int line);
int ErrMessage(MSKContext * G, const char *where, const char *what);

#define ErrChkPtr(G,p) {if(!p) ErrPointer(G,__FILE__,__LINE__);}

#endif
