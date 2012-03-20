/* 
 * COPYRIGHT NOTICE: This file contains original source code from the
 * Jenarix (TM) Library, Copyright (C) 2007-8 by Warren L. Delano of
 * DeLano Scientific LLC, Palo Alto, California, United States.
 * Please see the accompanying LICENSE file for further information.
 * All rights not explicitly granted in that LICENSE file are
 * reserved.  It is unlawful to modify or remove this notice.
 * TRADEMARK NOTICE: Jenarix is a Trademark of DeLano Scientific LLC.
*/

/* MACHINE PROCESSED SOURCE CODE -- DO NOT EDIT */

#ifndef _H_MSKIT_OV_PORT
#define _H_MSKIT_OV_PORT

/* headers */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#ifdef __linux__
#include <malloc.h>
#else
#include <stddef.h>
#endif

#include "ov_defines.h"
#include "ov_status.h"


/* memory management */

#define ov_os_malloc malloc
#define ov_os_calloc calloc
#define ov_os_realloc realloc
#define ov_os_free free
#define ov_os_memmove memmove
#define ov_os_memset memset

/* termination */

#define ov_os_abort abort
#define ov_os_exit exit

#endif
