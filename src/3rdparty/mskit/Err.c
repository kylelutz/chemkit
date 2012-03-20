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
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#include "os_predef.h"
#include "os_std.h"

#include "Err.h"

void ErrFatal(MSKContext * G, const char *where, const char *what)
{
  fprintf(stderr, "%s-Error: %s\n", where, what);
  fflush(stderr);
  exit(EXIT_FAILURE);
}

int ErrMessage(MSKContext * G, const char *where, const char *what)
{
  fprintf(stderr, "%s-Error: %s\n", where, what);
  return (0);
}

void ErrPointer(MSKContext * G, const char *file, int line)
{
  fprintf(stderr, "NULL-POINTER-ERROR: in %s line %i\n", file, line);
  fflush(stderr);
  exit(EXIT_FAILURE);
}
