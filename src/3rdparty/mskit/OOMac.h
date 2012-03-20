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

#ifndef _H_MSKIT_OOMAC
#define _H_MSKIT_OOMAC

#include "MemoryDebug.h"

#define OOAlloc(G,type) \
type *I;					\
I = (type*)mmalloc(sizeof(type));

#define OOCalloc(G,type) \
type *I;					 \
I = (type*)mcalloc(sizeof(type),1);

#define OOFreeP(ptr) \
{if(ptr) {mfree(ptr);ptr=NULL;}}

#endif
