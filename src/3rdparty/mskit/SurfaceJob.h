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

#ifndef _H_MSKIT_SURFACEJOB
#define _H_MSKIT_SURFACEJOB

#include "MSKContext.h"

#ifdef NT
#undef NT
#endif

typedef struct {
  float vdw;
  int flags;
} SurfaceJobAtomInfo;

typedef struct {
  /* input */
  float *coord;
  SurfaceJobAtomInfo *atomInfo;

  float maxVdw;
  int allVisibleFlag;

  int nPresent;
  int *presentVla;

  int solventSphereIndex, sphereIndex;

  int surfaceType;
  int circumscribe;
  float probeRadius;
  float carveCutoff;
  float *carveVla;

  int surfaceMode;
  int surfaceSolvent;
  int cavityCull;
  float pointSep;
  float trimCutoff;
  float trimFactor;

  int cavityMode;
  float cavityRadius;
  float cavityCutoff;

  /* results */
  float *V, *VN;
  int N, *T, *S, NT;

} SurfaceJob;


SurfaceJob *SurfaceJobNew(MSKContext * G);

void SurfaceJobFree(MSKContext * G, SurfaceJob * I);

int SurfaceJobRun(MSKContext * G, SurfaceJob * I);

void SurfaceJobPurgeResult(MSKContext * G, SurfaceJob * I);

#endif
