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

#ifndef _H_MSKIT_OS_PREDEF
#define _H_MSKIT_OS_PREDEF

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Macros used by Fortify source in GCC 4.1.x are incompatible with
   MSKIT's Feedback system... */

#ifdef _FORTIFY_SOURCE
#undef _FORTIFY_SOURCE
#endif


/* Alias-able typedefs */

#ifdef __GNUC__
#if __GNUC__ > 3
typedef int aliased_int __attribute__ ((may_alias));
typedef float aliased_float __attribute__ ((may_alias));
#else
typedef int aliased_int;
typedef float aliased_float;
#endif
#else
typedef int aliased_int;
typedef float aliased_float;
#endif

#ifdef WIN32
#define __inline__ __inline
#endif

#ifdef __linux__
#include <malloc.h>
#else
#include <stddef.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define fmax max
#define fmin min
#pragma warning (disable:4996)
#define snprintf sprintf_s
#endif

#include "ov_types.h"

#endif
