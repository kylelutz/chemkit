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

#ifndef _H_MSKIT_BASE
#define _H_MSKIT_BASE

#include "os_limits.h"
#include "os_types.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef NULL
#define NULL ((void*)0)
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef uchar
#define uchar unsigned char
#endif

#ifndef uint
#define uint unsigned int
#endif

#define MAX_VDW 2.5F            /* this has to go */

#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

#ifndef R_SMALL4
#define R_SMALL4 0.0001F
#endif

#ifndef R_SMALL8
#define R_SMALL8 0.00000001F
#endif

#define MAXLINELEN 1024

#endif
