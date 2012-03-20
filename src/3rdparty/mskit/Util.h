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

#ifndef _H_MSKIT_UTIL
#define _H_MSKIT_UTIL

#include "os_predef.h"

#include "MSKContext.h"

void UtilZeroMem(void *ptr, ov_size howMuch);
void UtilCopyMem(void *dst, void *src, ov_size howMuch);
void *UtilArrayCalloc(unsigned int *dim, ov_size ndim, ov_size atom_size);
char *UtilConcat(char *where, char *what);
void UtilNConcat(char *dst, char *str, ov_size n);
void UtilConcatVLA(char **vla, ov_size * cc, char *str);
void UtilNPadVLA(char **vla, ov_size * cc, char *str, ov_size len);
void UtilFillVLA(char **vla, ov_size * cc, char what, ov_size len);
void UtilNCopy(char *dst, char *src, ov_size n);        /* up to N-1 chars */
void UtilNCopyToLower(char *dst, char *src, ov_size n); /* up to N-1 chars */
void UtilCleanStr(char *s);
int UtilCountStringVLA(char *vla);

#if 0
double UtilGetSeconds(MSKContext * G);
int UtilInit(MSKContext * G);
void UtilFree(MSKContext * G);
#endif

typedef int UtilOrderFn(void *array, int l, int r);
void UtilSortIndex(int n, void *array, int *x, UtilOrderFn * fOrdered);

void UtilSemiSortFloatIndex(int n, float *array, int *x, int forward);

void UtilApplySortedIndices(int n, int *x, int rec_size, void *src, void *dst);

void UtilSortInPlace(MSKContext * G, void *array, int nItem, unsigned int itemSize,
                     UtilOrderFn * fOrdered);

void UtilExpandArrayElements(void *src, void *dst, int n_entries, int old_rec_size,
                             int new_rec_size);

typedef int UtilOrderFnGlobals(MSKContext * G, void *array, int l, int r);
void UtilSortIndexGlobals(MSKContext * G, int n, void *array, int *x,
                          UtilOrderFnGlobals * fOrdered);

int UtilShouldWePrintQuantity(int quantity);

#endif
