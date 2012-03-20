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

// MSKIT
#include "SurfaceJob.h"
#include "SolventDot.h"

void SurfaceJobPurgeResult(MSKContext * G, SurfaceJob * I)
{
  I->N  = 0;
  I->NT = 0;
  VLAFreeP(I->V);
  VLAFreeP(I->VN);
  VLAFreeP(I->T);
  VLAFreeP(I->S);
  FreeP(I->VC);
  FreeP(I->VA);
  I->oneColorFlag = true;
  I->oneAlphaFlag = true;
  I->oneColor = -1;
  I->oneAlpha = -1.0F;
}

SurfaceJob *SurfaceJobNew(MSKContext * G,
	                         float *coord, SurfaceJobAtomInfo * atom_info,
	                         float max_vdw, float probe_radius,
	                         int surface_quality, int surface_type,
	                         int surface_solvent, int cavity_cull,
	                         int cavity_mode, float cavity_radius,
	                         float cavity_cutoff, float trim_cutoff, float trim_factor)
{
  OOCalloc(G, SurfaceJob);

  if (I) {
  	I->coord = coord;
  	I->atomInfo = atom_info;
  	I->maxVdw = max_vdw;
    I->surfaceQuality = surface_quality;
  	I->surfaceType = surface_type;
  	I->surfaceSolvent = surface_solvent;
  	I->probeRadius = probe_radius;
  	I->trimCutoff = trim_cutoff;
  	I->trimFactor = trim_factor;
  	I->cavityCull = cavity_cull;
  	I->cavityMode = cavity_mode;
  	I->cavityRadius = cavity_radius;
  	I->cavityCutoff = cavity_cutoff;
    I->oneColorFlag = true;
    I->oneAlphaFlag = true;
    I->oneColor = -1;
    I->oneAlpha = -1.0F;

#define SURFACE_QUALITY_BEST_SEP       0.25F
#define SURFACE_QUALITY_NORMAL_SEP     0.5F
#define SURFACE_QUALITY_POOR_SEP       0.85F
#define SURFACE_QUALITY_MISERABLE_SEP  2.0F

    if(surface_quality >= 4) {
      /* totally impractical */
      I->pointSep = SURFACE_QUALITY_BEST_SEP / 4;
      I->sphereIndex = 4;
      I->solventSphereIndex = 4;
      I->circumscribe = 91;
    } else {
      switch (surface_quality) {
        case 3:
          /* nearly impractical */
          I->pointSep = SURFACE_QUALITY_BEST_SEP / 3;
          I->sphereIndex = 4;
          I->solventSphereIndex = 3;
          I->circumscribe = 71;
          break;
        case 2:
          /* nearly perfect */
          I->pointSep = SURFACE_QUALITY_BEST_SEP / 2;
          I->sphereIndex = 3;
          I->solventSphereIndex = 3;
          I->circumscribe = 41;
          break;
        case 1:
          /* good */
          I->pointSep = SURFACE_QUALITY_BEST_SEP;
          I->sphereIndex = 2;
          I->solventSphereIndex = 3;
          I->circumscribe = 40;
          break;
        case 0:
          /* 0 - normal */
          I->pointSep = SURFACE_QUALITY_NORMAL_SEP;
          I->sphereIndex = 1;
          I->solventSphereIndex = 2;
          if (surface_type == 6)
            I->circumscribe = 30;
          break;
        case -1:
          /* -1 poor */
          I->pointSep = SURFACE_QUALITY_POOR_SEP;
          I->sphereIndex = 1;
          I->solventSphereIndex = 2;
          if (surface_type == 6)
            I->circumscribe = 10;
          break;
        case -2:
          /* -2 god awful */
          I->pointSep = SURFACE_QUALITY_POOR_SEP * 1.5F;
          I->sphereIndex = 1;
          I->solventSphereIndex = 1;
          break;
        case -3:
          /* -3 miserable */
          I->pointSep = SURFACE_QUALITY_MISERABLE_SEP;
          I->sphereIndex = 1;
          I->solventSphereIndex = 1;
          break;
        default:
          I->pointSep = SURFACE_QUALITY_MISERABLE_SEP * 1.18F;
          I->sphereIndex = 0;
          I->solventSphereIndex = 1;
          break;
      }
    }

    if(!surface_solvent)
      I->circumscribe = 0;
  }
  
  return I;
}

void SurfaceJobFree(MSKContext * G, SurfaceJob * I)
{
  SurfaceJobPurgeResult(G, I);
  VLAFreeP(I->presentVla);
  VLAFreeP(I->atomInfo);
  VLAFreeP(I->coord);
  OOFreeP(I);
}

int SurfaceJobRun(MSKContext * G, SurfaceJob * I)
{
  int ok = true;
  int MaxN;
  int n_index = VLAGetSize(I->atomInfo);
  int n_present = I->nPresent;
  SphereRec *sp = G->Sphere->Sphere[I->sphereIndex];
  SphereRec *ssp = G->Sphere->Sphere[I->solventSphereIndex];

  MSKContextClean(G);
  G->Ready = false;

  OrthoBusyStage(G, 0);

  SurfaceJobPurgeResult(G, I);

  {
    /* compute limiting storage requirements */
    int tmp = n_present;
    if(tmp < 1)
      tmp = 1;
    if(sp->nDot < ssp->nDot)
      MaxN = tmp * ssp->nDot;
    else
      MaxN = tmp * sp->nDot;
  }

  I->V = VLAlloc(float, (MaxN + 1) * 3);
  I->VN = VLAlloc(float, (MaxN + 1) * 3);

  if(G->Interrupt)
    ok = false;

  if(!(I->V && I->VN) || (!ok)) {       /* bail out point -- try to reduce crashes */
    VLAFreeP(I->V);
    I->V = NULL;
    VLAFreeP(I->VN);
    I->VN = NULL;
    ok = false;
  } else {
    SolventDot *sol_dot = NULL;
    float *v = I->V;
    float *vn = I->VN;
    float probe_radius = I->probeRadius;
    int circumscribe = I->circumscribe;
    int surface_type = I->surfaceType;
    float point_sep = I->pointSep;
    float *I_coord = I->coord;
    int *present_vla = I->presentVla;
    SurfaceJobAtomInfo *I_atom_info = I->atomInfo;

    I->N = 0;

    sol_dot = SolventDotNew(G, I->coord, I->atomInfo, probe_radius,
                            ssp, present_vla,
                            circumscribe, I->surfaceSolvent,
                            I->cavityCull, I->maxVdw,
                            I->cavityMode, I->cavityRadius, I->cavityCutoff);

    if((!sol_dot) || (G->Interrupt))
      ok = false;

    OrthoBusyStage(G, 1);

    if(ok && sol_dot) {
      if(!I->surfaceSolvent) {

        float solv_tole = point_sep * 0.04F;
        float probe_rad_more;
        float probe_rad_less;
        float probe_rad_less2;

        if(probe_radius < (2.5F * point_sep)) { /* minimum probe radius allowed */
          probe_radius = 2.5F * point_sep;
        }

        probe_rad_more = probe_radius * (1.0F + solv_tole);

        switch (surface_type) {
        case 0:                /* solid */
        case 3:
        case 4:
        case 5:
        case 6:
          probe_rad_less = probe_radius;
          break;
        default:
          probe_rad_less = probe_radius * (1.0F - solv_tole);
          break;
        }
        probe_rad_less2 = probe_rad_less * probe_rad_less;

        if(surface_type >= 5) { /* effectively double-weights atom points */
          if(sol_dot->nDot) {
            int a;
            float *v0 = sol_dot->dot;
            float *n0 = sol_dot->dotNormal;
            for(a = 0; a < sol_dot->nDot; a++) {
              scale3f(n0, -probe_radius, v);
              add3f(v0, v, v);
              copy3f(n0, vn);
              v += 3;
              vn += 3;
              n0 += 3;
              v0 += 3;
              I->N++;
            }
          }
        }
        if(G->Interrupt)
          ok = false;
        if(ok) {
          MapType *map, *solv_map;
          map = MapNewFlagged(G, I->maxVdw + probe_rad_more,
                              I->coord, VLAGetSize(I->coord) / 3, NULL, NULL);

          solv_map = MapNew(G, probe_rad_less, sol_dot->dot, sol_dot->nDot, NULL);
          if(map && solv_map) {

            MapSetupExpress(solv_map);
            MapSetupExpress(map);

            if(sol_dot->nDot) {
              int *dc = sol_dot->dotCode;
              Vector3f *dot = Alloc(Vector3f, sp->nDot);
              float *v0, *n0;

              {
                int b;
                for(b = 0; b < sp->nDot; b++) {
                  scale3f(sp->dot[b], probe_radius, dot[b]);
                }
              }
              v0 = sol_dot->dot;
              n0 = sol_dot->dotNormal;
              {
                int a, b;
                float dist2 = probe_rad_less2;
                int sp_nDot = sp->nDot;
                for(a = 0; a < sol_dot->nDot; a++) {
                  if(dc[a] || (surface_type < 6)) {     /* surface type 6 is completely scribed */
                    OrthoBusyFast(G, a, sol_dot->nDot);
                    for(b = 0; b < sp_nDot; b++) {
                      float *dot_b = dot[b];
                      v[0] = v0[0] + dot_b[0];
                      v[1] = v0[1] + dot_b[1];
                      v[2] = v0[2] + dot_b[2];
                      {
                        int flag = true;
                        int ii;
                        if((ii = *(MapLocusEStart(solv_map, v)))) {
                          register float *i_dot = sol_dot->dot;
                          register float dist = probe_rad_less;
                          register int *elist_ii = solv_map->EList + ii;
                          register float v_0 = v[0];
                          register int jj_next, jj = *(elist_ii++);
                          register float v_1 = v[1];
                          register float *v1 = i_dot + 3 * jj;
                          register float v_2 = v[2];
                          while(jj >= 0) {
                            /* huge bottleneck -- optimized for superscaler processors */
                            register float dx = v1[0], dy, dz;
                            jj_next = *(elist_ii++);
                            dx -= v_0;
                            if(jj != a) {
                              dx = (dx < 0.0F) ? -dx : dx;
                              dy = v1[1] - v_1;
                              if(!(dx > dist)) {
                                dy = (dy < 0.0F) ? -dy : dy;
                                dz = v1[2] - v_2;
                                if(!(dy > dist)) {
                                  dx = dx * dx;
                                  dz = (dz < 0.0F) ? -dz : dz;
                                  dy = dy * dy;
                                  if(!(dz > dist)) {
                                    dx = dx + dy;
                                    dz = dz * dz;
                                    if(!(dx > dist2))
                                      if((dx + dz) <= dist2) {
                                        flag = false;
                                        break;
                                      }
                                  }
                                }
                              }
                            }
                            v1 = i_dot + 3 * jj_next;
                            jj = jj_next;
                          }
                        }

                        /* at this point, we have points on the interior of the solvent surface,
                           so now we need to further trim that surface to cover atoms that are present */

                        if(flag) {
                          register int i = *(MapLocusEStart(map, v));
                          if(i) {
                            register int j = map->EList[i++];
                            while(j >= 0) {
                              SurfaceJobAtomInfo *atom_info = I_atom_info + j;
                              if((!present_vla) || present_vla[j]) {
                                if(within3f
                                   (I_coord + 3 * j, v,
                                    atom_info->vdw + probe_rad_more)) {
                                  flag = false;
                                  break;
                                }
                              }
                              j = map->EList[i++];
                            }
                          }
                          if(!flag) {   /* compute the normals */
                            vn[0] = -sp->dot[b][0];
                            vn[1] = -sp->dot[b][1];
                            vn[2] = -sp->dot[b][2];
                            if(I->N < MaxN) {
                              I->N++;
                              v += 3;
                              vn += 3;
                            } else {
                              int v_offset = v - I->V;
                              int vn_offset = vn - I->VN;
                              MaxN = MaxN * 2;
                              VLASize(I->V, float, (MaxN + 1) * 3);
                              VLASize(I->VN, float, (MaxN + 1) * 3);
                              v = I->V + v_offset;
                              vn = I->VN + vn_offset;
                            }
                          }
                        }
                      }
                    }
                  }
                  v0 += 3;
                  n0 += 3;
                  if(G->Interrupt) {
                    ok = false;
                    break;
                  }
                }
              }
              FreeP(dot);
            }
          }
          MapFree(solv_map);
          MapFree(map);
        }
      } else {
        float *v0 = sol_dot->dot;
        float *n0 = sol_dot->dotNormal;
        int a;
        circumscribe = 0;
        if(sol_dot->nDot) {
          if (sol_dot->nDot > MaxN) {
            int v_offset = v - I->V;
            int vn_offset = vn - I->VN;
            MaxN = sol_dot->nDot;
            VLASize(I->V, float, (MaxN + 1) * 3);
            VLASize(I->VN, float, (MaxN + 1) * 3);
            v = I->V + v_offset;
            vn = I->VN + vn_offset;
          }
          for(a = 0; a < sol_dot->nDot; a++) {
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            I->N++;
          }
        }
      }
    }
    SolventDotFree(sol_dot);
    sol_dot = NULL;
    if(G->Interrupt)
      ok = false;
    if(ok) {
      int refine, ref_count = 1;

      if((surface_type == 0) && (circumscribe)) {
        ref_count = 2;          /* these constants need more tuning... */
      }

      for(refine = 0; refine < ref_count; refine++) {

        /* add new vertices in regions where curvature is very high
           or where there are gaps with no points */

        if(ok && I->N && (surface_type == 0) && (circumscribe)) {
          int n_new = 0;

          float neighborhood = 2.6 * point_sep; /* these constants need more tuning... */
          float dot_cutoff = 0.666;
          float insert_cutoff = 1.1 * point_sep;

          float map_cutoff = neighborhood;
          float *new_dot = VLAlloc(float, 1000);

          if(map_cutoff < (2.9 * point_sep)) {  /* these constants need more tuning... */
            map_cutoff = 2.9 * point_sep;
          }

          {
            MapType *map = NULL;
            int a;
            map = MapNew(G, map_cutoff, I->V, I->N, NULL);
            MapSetupExpress(map);
            v = I->V;
            vn = I->VN;
            for(a = 0; a < I->N; a++) {
              register int i = *(MapLocusEStart(map, v));
              if(i) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  if(j > a) {
                    float *v0 = I->V + 3 * j;
                    if(within3f(v0, v, map_cutoff)) {
                      int add_new = false;
                      float *n0 = I->VN + 3 * j;
                      VLACheck(new_dot, float, n_new * 6 + 5);
                      {
                        float *v1 = new_dot + n_new * 6;
                        average3f(v, v0, v1);
                        if((dot_product3f(n0, vn) < dot_cutoff)
                           && (within3f(v0, v, neighborhood)))
                          add_new = true;
                        else {
                          /* if points are too far apart, insert a new one */
                          register int ii = *(MapLocusEStart(map, v1));
                          if(ii) {
                            int found = false;
                            register int jj = map->EList[ii++];
                            while(jj >= 0) {
                              if(jj != j) {
                                float *vv0 = I->V + 3 * jj;
                                if(within3f(vv0, v1, insert_cutoff)) {
                                  found = true;
                                  break;
                                }
                              }
                              jj = map->EList[ii++];
                            }
                            if(!found)
                              add_new = true;
                          }
                        }
                        if(add_new) {
                          /* highly divergent */
                          float *n1 = v1 + 3;
                          n_new++;
                          average3f(vn, n0, n1);
                          normalize3f(n1);
                        }
                      }
                    }
                  }
                  j = map->EList[i++];
                }
              }
              v += 3;
              vn += 3;
            }
            MapFree(map);
          }
          if(n_new) {
            float *n1 = new_dot + 3;
            float *v1 = new_dot;
            VLASize(I->V, float, 3 * (I->N + n_new));
            VLASize(I->VN, float, 3 * (I->N + n_new));
            v = I->V + 3 * I->N;
            vn = I->VN + 3 * I->N;
            I->N += n_new;
            while(n_new--) {
              copy3f(v1, v);
              copy3f(n1, vn);
              v += 3;
              vn += 3;
              v1 += 6;
              n1 += 6;
            }
          }
          VLAFreeP(new_dot);
        }

        if(ok && I->N && (surface_type == 0) && (circumscribe)) {

          float cutoff = 0.5 * probe_radius;

          /* combine scribing with an atom proximity cleanup pass */

          int *dot_flag = Calloc(int, I->N);
          MapType *map =
            MapNewFlagged(G, I->maxVdw + probe_radius, I_coord, n_index, NULL,
                          present_vla);
          int a;
          MapSetupExpress(map);
          v = I->V;
          for(a = 0; a < I->N; a++) {
            register int i = *(MapLocusEStart(map, v));
            if(i) {
              register int j = map->EList[i++];
              while(j >= 0) {
                SurfaceJobAtomInfo *atom_info = I_atom_info + j;
                if((!present_vla) || present_vla[j]) {
                  if(within3f(I_coord + 3 * j, v, atom_info->vdw + cutoff)) {
                    dot_flag[a] = true;
                  }
                }
                j = map->EList[i++];
              }
            }
            v += 3;
            if(G->Interrupt) {
              ok = false;
              break;
            }
          }

          MapFree(map);
          map = NULL;

          if(ok) {
            /* purge unused dots */

            float *v0 = I->V;
            float *vn0 = I->VN;
            int *p = dot_flag;
            int c = I->N;
            int a;
            v = I->V;
            vn = I->VN;
            I->N = 0;
            for(a = 0; a < c; a++) {
              if(*(p++)) {
                *(v0++) = *(v++);
                *(v0++) = *(v++);
                *(v0++) = *(v++);
                *(vn0++) = *(vn++);
                *(vn0++) = *(vn++);
                *(vn0++) = *(vn++);
                I->N++;
              } else {
                v += 3;
                vn += 3;
              }
            }
          }
          FreeP(dot_flag);
        }

        if(ok && I->N) {
          int repeat_flag = true;
          float min_dot = 0.1F;
          int *dot_flag = Alloc(int, I->N);
          while(repeat_flag) {
            repeat_flag = false;

            if(surface_type >= 3) {
              register int jj;
              float dist;
              register float nearest;
              register float min_sep2 = point_sep * point_sep;
              float diff[3];
              {
                int a;
                for(a = 0; a < I->N; a++)
                  dot_flag[a] = 1;
              }
              {
                MapType *map = MapNew(G, point_sep + 0.05F, I->V, I->N, NULL);
                int a;
                MapSetupExpress(map);
                v = I->V;
                vn = I->VN;
                for(a = 0; a < I->N; a++) {
                  if(dot_flag[a]) {
                    register int i = *(MapLocusEStart(map, v));
                    if(i) {
                      register int j = map->EList[i++];
                      jj = I->N;
                      nearest = point_sep + 1.0F;
                      while(j >= 0) {
                        if(j > a) {
                          if(dot_flag[j]) {
                            if(dot_product3f(I->VN + (3 * j), vn) > min_dot) {
                              if(within3fret
                                 (I->V + (3 * j), v, point_sep, min_sep2, diff, &dist)) {
                                repeat_flag = true;
                                if(dist < nearest) {
                                  /* try to be as determinstic as possible
                                     in terms of how we collapse points */
                                  jj = j;
                                  nearest = dist;
                                } else if((j < jj) && (fabs(dist - nearest) < R_SMALL4)) {
                                  jj = j;
                                  nearest = dist;
                                }
                              }
                            }
                          }
                        }
                        j = map->EList[i++];
                      }

                      if(jj < I->N) {
                        dot_flag[jj] = 0;
                        add3f(vn, I->VN + (3 * jj), vn);
                        average3f(I->V + (3 * jj), v, v);
                        repeat_flag = true;
                      }
                    }
                  }
                  v += 3;
                  vn += 3;
                  if(G->Interrupt) {
                    ok = false;
                    break;
                  }
                }
                MapFree(map);
              }
            } else {            /* surface types < 3 */
              int a;
              MapType *map = MapNew(G, -point_sep, I->V, I->N, NULL);
              for(a = 0; a < I->N; a++)
                dot_flag[a] = 1;
              MapSetupExpress(map);
              v = I->V;
              vn = I->VN;
              for(a = 0; a < I->N; a++) {
                if(dot_flag[a]) {
                  register int i = *(MapLocusEStart(map, v));
                  if(i) {
                    register int j = map->EList[i++];
                    while(j >= 0) {
                      if(j != a) {
                        if(dot_flag[j]) {
                          if(within3f(I->V + (3 * j), v, point_sep)) {
                            dot_flag[j] = 0;
                            add3f(vn, I->VN + (3 * j), vn);
                            average3f(I->V + (3 * j), v, v);
                            repeat_flag = true;
                          }
                        }
                      }
                      j = map->EList[i++];
                    }
                  }
                }
                v += 3;
                vn += 3;
                if(G->Interrupt) {
                  ok = false;
                  break;
                }
              }
              MapFree(map);
            }

            if(ok) {
              float *v0 = I->V;
              float *vn0 = I->VN;
              int *p = dot_flag;
              int c = I->N;
              int a;
              v = I->V;
              vn = I->VN;
              I->N = 0;
              for(a = 0; a < c; a++) {
                if(*(p++)) {
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
                  normalize3f(vn);
                  *(vn0++) = *(vn++);
                  *(vn0++) = *(vn++);
                  *(vn0++) = *(vn++);
                  I->N++;
                } else {
                  v += 3;
                  vn += 3;
                }
              }
            }
            if(G->Interrupt) {
              ok = false;
              break;
            }
          }
          FreeP(dot_flag);
        }

        /* now eliminate troublesome vertices in regions of extremely high curvature */

        if(ok && (surface_type != 3) &&
           I->N && (I->trimCutoff > 0.0F) && (I->trimFactor > 0.0F)) {
          float trim_cutoff = I->trimCutoff;
          float trim_factor = I->trimFactor;
          int repeat_flag = true;
          float neighborhood = trim_factor * point_sep;
          float dot_sum;
          int n_nbr;
          int *dot_flag = Alloc(int, I->N);
          if(surface_type == 6) {       /* emprical tweaks */
            trim_factor *= 2.5;
            trim_cutoff *= 1.5;
          }
          while(repeat_flag) {
            int a;
            MapType *map = MapNew(G, neighborhood, I->V, I->N, NULL);
            repeat_flag = false;

            for(a = 0; a < I->N; a++)
              dot_flag[a] = 1;
            MapSetupExpress(map);
            v = I->V;
            vn = I->VN;
            for(a = 0; a < I->N; a++) {
              if(dot_flag[a]) {
                register int i = *(MapLocusEStart(map, v));
                if(i) {
                  register int j = map->EList[i++];
                  n_nbr = 0;
                  dot_sum = 0.0F;
                  while(j >= 0) {
                    if(j != a) {
                      if(dot_flag[j]) {
                        float *v0 = I->V + 3 * j;
                        if(within3f(v0, v, neighborhood)) {
                          float *n0 = I->VN + 3 * j;
                          dot_sum += dot_product3f(n0, vn);
                          n_nbr++;
                        }
                      }
                    }
                    j = map->EList[i++];
                  }

                  if(n_nbr) {
                    dot_sum /= n_nbr;
                    if(dot_sum < trim_cutoff) {
                      dot_flag[a] = false;
                      repeat_flag = true;
                    }
                  }
                }
              }
              v += 3;
              vn += 3;
              if(G->Interrupt) {
                ok = false;
                break;
              }
            }

            if(ok) {
              float *v0 = I->V;
              float *vn0 = I->VN;
              int *p = dot_flag;
              int c = I->N;
              v = I->V;
              vn = I->VN;
              I->N = 0;
              for(a = 0; a < c; a++) {
                if(*(p++)) {
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
                  normalize3f(vn);
                  *(vn0++) = *(vn++);
                  *(vn0++) = *(vn++);
                  *(vn0++) = *(vn++);
                  I->N++;
                } else {
                  v += 3;
                  vn += 3;
                }
              }
            }
            MapFree(map);
            if(G->Interrupt) {
              ok = false;
              break;
            }
          }
          FreeP(dot_flag);
        }
        if(G->Interrupt) {
          ok = false;
          break;
        }
      }
    }

    if(I->N && I->V && I->VN) {
      VLASizeForSure(I->V, float, 3 * I->N);
      VLASizeForSure(I->VN, float, 3 * I->N);
    }

    if(G->Interrupt)
      ok = false;

    OrthoBusyStage(G, 2);

    if(ok && I->N) {
      if(surface_type != 1) {   /* not a dot surface... */
        float cutoff = point_sep * 5.0F;
        if((cutoff > probe_radius) && (!I->surfaceSolvent))
          cutoff = probe_radius;
        I->T = TrianglePointsToSurface(G, I->V, I->VN, I->N, cutoff, &I->NT, &I->S, NULL, 
                                       I->cavityMode);
      }
    } else {
      VLASizeForSure(I->V, float, 1);
      VLASizeForSure(I->VN, float, 1);
    }
  }
  
  G->Ready = true;

  return ok;
}

void SurfaceJobColoring(MSKContext *G, SurfaceJob * I, const int *colors, const float *transp)
{
  SurfaceJobAtomInfo *ai2 = NULL;
  MapType *map;
  int a, i0, i, j;
  int c0, c1, *vc;
  float a0, a1, *va;
  float *v0, *n0;
  float probe_radius;
  float dist;

  /* assert(colors != NULL); */

  probe_radius = I->probeRadius;
  ai2 = I->atomInfo;

  if(I->N) {
    I->oneColorFlag = true;
    I->oneAlphaFlag = true;
    I->oneColor = -1;
    I->oneAlpha = -1.0F;

    if(!I->VC)
      I->VC = Alloc(int, I->N);
    vc = I->VC;

    if (transp) {
      if(!I->VA)
        I->VA = Alloc(float, I->N);
      va = I->VA;
    }

    c0 = -1;
    a0 = -1.0F;

    /* now, assign colors to each point */
    map = MapNewFlagged(G, 2 * I->maxVdw + probe_radius,
                        I->coord, VLAGetSize(I->coord) / 3, NULL, NULL);

    if(map) {
      MapSetupExpress(map);

      for(a = 0; a < I->N; a++) {
        SurfaceJobAtomInfo *ai0 = NULL;
        float minDist = MAXFLOAT;

        v0 = I->V + 3 * a;
        n0 = I->VN + 3 * a;
        i0 = -1;

        /* colors */
        i = *(MapLocusEStart(map, v0));
        if(i) {
          j = map->EList[i++];
          while(j >= 0) {
            ai2 = I->atomInfo + j;
            dist = (float) diff3f(v0, I->coord + j * 3) - ai2->vdw;
            if(dist < minDist) {
              i0 = j;
              ai0 = ai2;
              minDist = dist;
            }
            j = map->EList[i++];
          }
        }

        c1 = (i0 >= 0)
                ? *(colors + i0)
                : -1;
        *(vc++) = c1;

        if(I->oneColorFlag) {
          if(c0 >= 0) {
            if(c0 != c1)
              I->oneColorFlag = false;
          } else
            c0 = c1;
        }

        if (transp) {
          a1 = (i0 >= 0)
                  ? (1.0F - *(transp + i0))
                  : -1.0F;
          *(va++) = a1;

          if(I->oneAlphaFlag) {
            if(a0 >= 0) {
              if(a0 != a1)
                I->oneAlphaFlag = false;
            } else
              a0 = a1;
          }
        }
      }

      MapFree(map);
    }

    if(I->oneAlphaFlag) {
      I->oneAlpha = a0;
      if(I->VA) {
        FreeP(I->VA);
        I->VA = NULL;
      }
    }
    else
      I->oneColorFlag = false;

    if(I->oneColorFlag) {
      I->oneColor = c0;
      if(I->VC) {
        FreeP(I->VC);
        I->VC = NULL;
      }
    }
  }

  /*
     if(surface_color>=0) {
     I->oneColorFlag=true;
     I->oneColor=surface_color;
     }
   */

}
