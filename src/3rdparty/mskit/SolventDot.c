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

#include "os_predef.h"
#include "os_std.h"

#include "Base.h"
#include "MemoryDebug.h"
#include "OOMac.h"
#include "Map.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Vector.h"
#include "Err.h"

// MSKIT
#include "SurfaceJob.h"
#include "SolventDot.h"

SolventDot *SolventDotNew(MSKContext * G,
	                         float *coord,
	                         SurfaceJobAtomInfo * atom_info,
	                         float probe_radius, SphereRec * sp,
	                         int *present,
	                         int circumscribe,
	                         int surface_solvent, int cavity_cull,
	                         float max_vdw,
	                         int cavity_mode, float cavity_radius, 
	                         float cavity_cutoff)
{
  int ok = true;
  int b;
  float vdw;
  float probe_radius_plus;
  int maxDot = 0;
  int stopDot;
  int n_coord = VLAGetSize(atom_info);
  Vector3f *sp_dot = sp->dot;
  OOCalloc(G, SolventDot);

  /*  printf("%p %p %p %f %p %p %p %d %d %d %f\n",
     G,
     coord,
     atom_info,
     probe_radius,sp,
     extent,present,
     circumscribe, 
     surface_solvent,  cavity_cull,
     max_vdw);
   */

  stopDot = n_coord * sp->nDot + 2 * circumscribe;
  I->dot = VLAlloc(float, (stopDot + 1) * 3);
  I->dotNormal = VLAlloc(float, (stopDot + 1) * 3);
  I->dotCode = VLACalloc(int, stopDot + 1);

  probe_radius_plus = probe_radius * 1.5F;

  I->nDot = 0;
  {
    int dotCnt = 0;
    MapType *map =
      MapNewFlagged(G, max_vdw + probe_radius, coord, n_coord, NULL, present);
    if(G->Interrupt)
      ok = false;
    if(map && ok) {
      float *v = I->dot;
      float *n = I->dotNormal;
      int *dc = I->dotCode;
      int maxCnt = 0;

      MapSetupExpress(map);
      {
        int a;
        int skip_flag;

        SurfaceJobAtomInfo *a_atom_info = atom_info;
        for(a = 0; a < n_coord; a++) {
          OrthoBusyFast(G, a, n_coord);
          if((!present) || (present[a])) {
            register int i;
            float *v0 = coord + 3 * a;
            vdw = a_atom_info->vdw + probe_radius;

            skip_flag = false;

            i = *(MapLocusEStart(map, v0));
            if(i) {
              register int j = map->EList[i++];
              while(j >= 0) {
                SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                if(j > a)       /* only check if this is atom trails */
                  if((!present) || present[j]) {
                    if(j_atom_info->vdw == a_atom_info->vdw) {        /* handle singularities */
                      float *v1 = coord + 3 * j;
                      if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                        skip_flag = true;
                    }
                  }
                j = map->EList[i++];
              }
            }
            if(!skip_flag) {
              for(b = 0; b < sp->nDot; b++) {
                float *sp_dot_b = (float*)(sp_dot + b);
                register int i;
                int flag = true;
                v[0] = v0[0] + vdw * (n[0] = sp_dot_b[0]);
                v[1] = v0[1] + vdw * (n[1] = sp_dot_b[1]);
                v[2] = v0[2] + vdw * (n[2] = sp_dot_b[2]);
                i = *(MapLocusEStart(map, v));
                if(i) {
                  register int j = map->EList[i++];
                  while(j >= 0) {
                    SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                    if((!present) || present[j]) {
                      if(j != a) {
                        skip_flag = false;
                        if(j_atom_info->vdw == a_atom_info->vdw) {      /* handle singularities */
                          float *v1 = coord + 3 * j;
                          if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                            skip_flag = true;
                        }
                        if(!skip_flag)
                          if(within3f(coord + 3 * j, v, j_atom_info->vdw + probe_radius)) {
                            flag = false;
                            break;
                          }
                      }
                    }
                    j = map->EList[i++];
                  }
                }
                if(flag && (dotCnt < stopDot)) {
                  dotCnt++;
                  v += 3;
                  n += 3;
                  dc++;
                  I->nDot++;
                }
              }
            }
            if(dotCnt > maxCnt) {
              maxCnt = dotCnt;
              maxDot = I->nDot - 1;
            }
          }
          a_atom_info++;
        }
      }

      /* for each pair of proximal atoms, circumscribe a circle for their intersection */

      /*    CGOReset(G->DebugCGO); */

      {
        MapType *map2 = NULL;
        if(circumscribe && (!surface_solvent))
          map2 =
            MapNewFlagged(G, 2 * (max_vdw + probe_radius), coord, n_coord, NULL, present);

        if(G->Interrupt)
          ok = false;
        if(map2 && ok) {
          /*        CGOBegin(G->DebugCGO,GL_LINES); */
          int a;
          int skip_flag;

          SurfaceJobAtomInfo *a_atom_info = atom_info;
          MapSetupExpress(map2);
          for(a = 0; a < n_coord; a++) {
            if((!present) || present[a]) {
              register int i;
              float vdw2;
              float *v0 = coord + 3 * a;
              vdw = a_atom_info->vdw + probe_radius;
              vdw2 = vdw * vdw;

              skip_flag = false;

              i = *(MapLocusEStart(map2, v0));
              if(i) {
                register int j = map2->EList[i++];
                while(j >= 0) {
                  SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                  if(j > a)     /* only check if this is atom trails */
                    if((!present) || present[j]) {
                      if(j_atom_info->vdw == a_atom_info->vdw) {      /* handle singularities */
                        float *v2 = coord + 3 * j;
                        if((v0[0] == v2[0]) && (v0[1] == v2[1]) && (v0[2] == v2[2]))
                          skip_flag = true;
                      }
                    }
                  j = map2->EList[i++];
                }
              }

              if(!skip_flag) {
                register int ii = *(MapLocusEStart(map2, v0));
                if(ii) {
                  register int jj = map2->EList[ii++];
                  while(jj >= 0) {
                    SurfaceJobAtomInfo *jj_atom_info = atom_info + jj;
                    float dist;
                    if(jj > a)  /* only check if this is atom trails */
                      if((!present) || present[jj]) {
                        float vdw3 = jj_atom_info->vdw + probe_radius;

                        float *v2 = coord + 3 * jj;
                        dist = (float) diff3f(v0, v2);
                        if((dist > R_SMALL4) && (dist < (vdw + vdw3))) {
                          float vz[3], vx[3], vy[3], vp[3];
                          float tri_a = vdw, tri_b = vdw3, tri_c = dist;
                          float tri_s = (tri_a + tri_b + tri_c) * 0.5F;
                          float area = (float) sqrt1f(tri_s * (tri_s - tri_a) *
                                                      (tri_s - tri_b) * (tri_s - tri_c));
                          float radius = (2 * area) / dist;
                          float adj = (float) sqrt1f(vdw2 - radius * radius);

                          subtract3f(v2, v0, vz);
                          get_system1f3f(vz, vx, vy);

                          copy3f(vz, vp);
                          scale3f(vp, adj, vp);
                          add3f(v0, vp, vp);

                          for(b = 0; b <= circumscribe; b++) {
                            float xcos = (float) cos((b * 2 * cPI) / circumscribe);
                            float ysin = (float) sin((b * 2 * cPI) / circumscribe);
                            float xcosr = xcos * radius;
                            float ysinr = ysin * radius;
                            int flag = true;
                            v[0] = vp[0] + vx[0] * xcosr + vy[0] * ysinr;
                            v[1] = vp[1] + vx[1] * xcosr + vy[1] * ysinr;
                            v[2] = vp[2] + vx[2] * xcosr + vy[2] * ysinr;

                            i = *(MapLocusEStart(map, v));
                            if(i) {
                              register int j = map->EList[i++];
                              while(j >= 0) {
                                SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                                if((!present) || present[j])
                                  if((j != a) && (j != jj)) {
                                    skip_flag = false;
                                    if(a_atom_info->vdw == j_atom_info->vdw) {  /* handle singularities */
                                      float *v1 = coord + 3 * j;
                                      if((v0[0] == v1[0]) &&
                                         (v0[1] == v1[1]) && (v0[2] == v1[2]))
                                        skip_flag = true;
                                    }
                                    if(jj_atom_info->vdw == j_atom_info->vdw) { /* handle singularities */
                                      float *v1 = coord + 3 * j;
                                      if((v2[0] == v1[0]) &&
                                         (v2[1] == v1[1]) && (v2[2] == v1[2]))
                                        skip_flag = true;
                                    }
                                    if(!skip_flag)
                                      if(within3f(coord + 3 * j, v,
                                          j_atom_info->vdw + probe_radius)) {
                                        flag = false;
                                        break;
                                      }
                                  }
                                j = map->EList[i++];
                              }
                            }
                            if(flag && (dotCnt < stopDot)) {
                              float vt0[3], vt2[3];
                              subtract3f(v0, v, vt0);
                              subtract3f(v2, v, vt2);
                              normalize3f(vt0);
                              normalize3f(vt2);
                              add3f(vt0, vt2, n);
                              invert3f(n);
                              normalize3f(n);
                              /*
                                 n[0] = vx[0] * xcos + vy[0] * ysin;
                                 n[1] = vx[1] * xcos + vy[1] * ysin;
                                 n[2] = vx[2] * xcos + vy[2] * ysin;
                               */

                              *dc = 1;  /* mark as exempt */

                              dotCnt++;
                              v += 3;
                              n += 3;
                              dc++;
                              I->nDot++;
                            }
                          }
                        }
                      }
                    jj = map2->EList[ii++];
                  }
                }
              }
            }
            a_atom_info++;
            if(G->Interrupt) {
              ok = false;
              break;
            }
          }
        }
        MapFree(map2);
      }
    }
    MapFree(map);
  }

  if(cavity_mode) {
    float *cavityDot = VLAlloc(float, (stopDot + 1) * 3);
    int nCavityDot = 0;
    int dotCnt = 0;
    if(cavity_radius<0.0F) {
      cavity_radius = - probe_radius * cavity_radius;
    }
    if(cavity_cutoff<0.0F) {
      cavity_cutoff = cavity_radius - cavity_cutoff * probe_radius;
    }
    {
      MapType *map =
        MapNewFlagged(G, max_vdw + cavity_radius, coord, n_coord, NULL, present);
      if(G->Interrupt)
        ok = false;
      if(map && ok) {
        float *v = cavityDot;
        MapSetupExpress(map);
        {
          int a;
          int skip_flag;
          SurfaceJobAtomInfo *a_atom_info = atom_info;
          for(a = 0; a < n_coord; a++) {
            if((!present) || (present[a])) {
              register int i;
              float *v0 = coord + 3 * a;
              vdw = a_atom_info->vdw + cavity_radius;
              skip_flag = false;
              i = *(MapLocusEStart(map, v0));
              if(i) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                  if(j > a)       /* only check if this is atom trails */
                    if((!present) || present[j]) {
                      if(j_atom_info->vdw == a_atom_info->vdw) {        /* handle singularities */
                        float *v1 = coord + 3 * j;
                        if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                          skip_flag = true;
                      }
                    }
                  j = map->EList[i++];
                }
              }
              if(!skip_flag) {
                for(b = 0; b < sp->nDot; b++) {
                  float *sp_dot_b = (float*)(sp_dot + b);
                  register int i;
                  int flag = true;
                  v[0] = v0[0] + vdw * (sp_dot_b[0]);
                  v[1] = v0[1] + vdw * (sp_dot_b[1]);
                  v[2] = v0[2] + vdw * (sp_dot_b[2]);
                  i = *(MapLocusEStart(map, v));
                  if(i) {
                    register int j = map->EList[i++];
                    while(j >= 0) {
                      SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                      if((!present) || present[j]) {
                        if(j != a) {
                          skip_flag = false;
                          if(j_atom_info->vdw == a_atom_info->vdw) {
                            /* handle singularities */
                            float *v1 = coord + 3 * j;
                            if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                              skip_flag = true;
                          }
                          if(!skip_flag) {
                            if(within3f(coord + 3 * j, v, j_atom_info->vdw + cavity_radius)) {
                              flag = false;
                              break;
                            }
                          }
                        }
                      }
                      j = map->EList[i++];
                    }
                  }
                  if(flag && (dotCnt < stopDot)) {
                    v += 3;
                    nCavityDot++;
                    dotCnt++;
                  }
                }
              }
            }
            a_atom_info++;
          }
        }
      }
      MapFree(map);
    }
    {
      int *dot_flag = Calloc(int, I->nDot);
      ErrChkPtr(G, dot_flag);
      {
        MapType *map = MapNew(G, cavity_cutoff, cavityDot, nCavityDot, NULL);
        if(map) {
          MapSetupExpress(map);
          {
            int *p = dot_flag;
            float *v = I->dot;
            int a;
            for(a = 0; a < I->nDot; a++) {
              register int i = *(MapLocusEStart(map, v));
              if(i) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  if(within3f(cavityDot + (3 * j), v, cavity_cutoff)) {
                    *p = true;
                    break;
                  }
                  j = map->EList[i++];
                }
              }
              v += 3;
              p++;
              if(G->Interrupt) {
                ok = false;
                break;
              }
            }
          }
        }
        MapFree(map);
      }
      
      {
        float *v0 = I->dot;
        float *n0 = I->dotNormal;
        int *dc0 = I->dotCode;
        int *p = dot_flag;
        int c = I->nDot;
        float *n = n0;
        float *v = v0;
        int *dc = dc0;
        int a;
        I->nDot = 0;
        for(a = 0; a < c; a++) {
          if(!*(p++)) {
            *(v0++) = *(v++);
            *(n0++) = *(n++);
            *(v0++) = *(v++);
            *(n0++) = *(n++);
            *(v0++) = *(v++);
            *(n0++) = *(n++);
            *(dc0++) = *(dc++);
            I->nDot++;
          } else {
            v += 3;
            n += 3;
          }
        }
      }
      FreeP(dot_flag);
    }
    VLAFreeP(cavityDot);

  } 
  if(ok && (cavity_mode != 1) && (cavity_cull > 0) && 
     (probe_radius > 0.75F) && (!surface_solvent)) {
    int *dot_flag = Calloc(int, I->nDot);
    ErrChkPtr(G, dot_flag);

#if 0
    dot_flag[maxDot] = 1;       /* this guarantees that we have a valid dot */
#endif

    {
      MapType *map = MapNew(G, probe_radius_plus, I->dot, I->nDot, NULL);
      if(map) {
        int flag = true;
        MapSetupExpress(map);
        while(flag) {
          int *p = dot_flag;
          float *v = I->dot;
          int a;
          flag = false;
          for(a = 0; a < I->nDot; a++) {
            if(!dot_flag[a]) {
              register int i = *(MapLocusEStart(map, v));
              int cnt = 0;

              if(i) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  if(j != a) {
                    if(within3f(I->dot + (3 * j), v, probe_radius_plus)) {
                      if(dot_flag[j]) {
                        *p = true;
                        flag = true;
                        break;
                      }
                      cnt++;
                      if(cnt > cavity_cull) {
                        *p = true;
                        flag = true;
                        break;
                      }
                    }
                  }
                  j = map->EList[i++];
                }
              }
            }
            v += 3;
            p++;
          }
          if(G->Interrupt) {
            ok = false;
            break;
          }
        }
      }
      MapFree(map);
    }

    {
      float *v0 = I->dot;
      float *n0 = I->dotNormal;
      int *dc0 = I->dotCode;
      int *p = dot_flag;
      int c = I->nDot;
      float *n = n0;
      float *v = v0;
      int *dc = dc0;
      int a;
      I->nDot = 0;
      for(a = 0; a < c; a++) {
        if(*(p++)) {
          *(v0++) = *(v++);
          *(n0++) = *(n++);
          *(v0++) = *(v++);
          *(n0++) = *(n++);
          *(v0++) = *(v++);
          *(n0++) = *(n++);
          *(dc0++) = *(dc++);
          I->nDot++;
        } else {
          v += 3;
          n += 3;
        }
      }
    }

    FreeP(dot_flag);
  }

  if(!ok) {
    SolventDotFree(I);
    I = NULL;
  }

  return I;
}

void SolventDotFree(SolventDot * I)
{
  if(I) {
    VLAFreeP(I->dot);
    VLAFreeP(I->dotNormal);
    VLAFreeP(I->dotCode);
  }
  OOFreeP(I);
}
