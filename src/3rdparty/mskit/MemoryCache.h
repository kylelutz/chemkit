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

#ifndef _H_MSKIT_MEMORYCACHE
#define _H_MSKIT_MEMORYCACHE

#define _MemoryCache_OFF

/* memory cache off */

#define VLACacheCheck(G,ptr,type,rec,t,i) VLACheck(ptr,type,rec)
#define VLACacheAlloc(G,type,initSize,t,i) VLAlloc(type,initSize)
#define VLACacheFreeP(G,ptr,t,i,f) VLAFreeP(ptr)
#define VLACacheSize(G,ptr,type,size,t,i) VLASize(ptr,type,size)
#define VLACacheSizeForSure(G,ptr,type,size,t,i) VLASizeForSure(ptr,type,size)
#define VLACacheExpand(G,ptr,rec,thread_index,i) VLAExpand(ptr,rec)

#define MemoryCacheInit(x)
#define MemoryCacheDone(x)
#define MemoryCacheReplaceBlock(G,g,o,n)

#define VLACacheMalloc(G,a,b,c,d,t,i) VLAMalloc(a,b,c,d)
#define VLACacheFree(G,p,t,i,f) VLAFree(p)

#define CacheAlloc(G,type,size,thread,id) (type*)mmalloc(sizeof(type)*(size))
#define CacheCalloc(G,type,size,thread,id) (type*)mcalloc(sizeof(type),size)
#define CacheRealloc(G,ptr,type,size,thread,id) (type*)mrealloc(sizeof(type)*(size))
#define CacheFreeP(G,ptr,thread,id,force) {if(ptr) {mfree(ptr);ptr=NULL;}}

#endif
