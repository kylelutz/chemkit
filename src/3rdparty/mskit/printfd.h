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

#ifndef _H_MSKIT_PRINTFD
#define _H_MSKIT_PRINTFD

#ifdef _DEBUG
#define PRINTFD(G,sysmod) { fprintf(stderr, "[" #sysmod "]"
#define ENDFD   ); fflush(stderr); }
#else
static __inline__ void __dummy_printfd(const char *s, ...) {}
#define PRINTFD(G,sysmod) { if(0) { __dummy_printfd(
#define ENDFD ); }}
#endif

#endif