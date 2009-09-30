#include <nongui/nonguilib.h>
#include <cmath>
#include <cstdio>

#include "maplib.h"
#include <math/mathlib.h>

using namespace chemlib;

#define X 0
#define Y 1
#define Z 2

static float InclusionRadSq, InclusionCenter[3];


static inline bool any_over(int i1, int i2, int i3, int level) {
  return i1 >= level || i2 >= level || i3 >= level;
}

static inline bool any_under(int i1, int i2, int i3, int level) {
  return i1 < level || i2 < level || i3 < level;
}


static int inside_sphere(float p[3]) {
  register float x, y, z;
  x = p[X] - InclusionCenter[X];
  x *= x;
  y = p[Y] - InclusionCenter[Y];
  y *= y;
  z = p[Z] - InclusionCenter[Z];
  z *= z;
  return (x + y + z <= InclusionRadSq);
}

static bool out_of_memory = false;

void EMapBase::Contour(float center[3], MIAtomList * current) {
  if (map_points.size() < 8) {
    return;
  }
  CurrentAtoms = current;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  float fcenter[3] = {center[0], center[1], center[2]};
  int planedirection, ilevel;
  mapheader->CtoF(&fcenter[0], &fcenter[1], &fcenter[2]);
  fcenter[0] *= (float)mapheader->nx;
  fcenter[1] *= (float)mapheader->ny;
  fcenter[2] *= (float)mapheader->nz;
  float xr = (float)mapheader->nx*settings->Radius/mapheader->a;
  float yr = (float)mapheader->ny*settings->Radius/mapheader->b;
  float zr = (float)mapheader->nz*settings->Radius/mapheader->c;
  xmin = ROUND(fcenter[0]-xr);
  xmax = ROUND(fcenter[0]+xr);
  ymin = ROUND(fcenter[1]-yr);
  ymax = ROUND(fcenter[1]+yr);
  zmin = ROUND(fcenter[2]-zr);
  zmax = ROUND(fcenter[2]+zr);
  edges.clear();
  out_of_memory = false;
  for (planedirection = 0; planedirection <= 2; planedirection++) {
    for (ilevel = 0; ilevel < 5; ilevel++) {
      if (settings->MapLevelOn[ilevel]) {
        if (!out_of_memory) {
          if (contur_sec(planedirection, xmin, xmax, ymin, ymax, zmin, zmax,
                (int) settings->MapLevel[ilevel], settings->MapLevelColor[ilevel], center) == -1) {
            out_of_memory = true;
          }
        }
      }
    }
  }
  SetChanged(false);
  Logger::log("Contoured %d edges", edges.size());
  map_center[0] = center[0];
  map_center[1] = center[1];
  map_center[2] = center[2];
}

int
EMapBase::contur_sec(int planedirection, int xmin, int xmax, int ymin, int ymax,
                     int zmin, int zmax, int level, int color, float center[3]) {
  float* sec;
  /* sec contains the rho values of the section to be
   * contoured.
   */

  double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0, z = 0;
  double c1 = 0, c2 = 0, r1 = 0, r2 = 0;
  float c[2][3], p[3], rint = 0, cint = 0, pint = 0;
  float fx, fy, fz;
  double* row1p = NULL, * row2p = NULL, * col1p = NULL, * col2p = NULL, * plane1p = NULL, * plane2p = NULL;
  float aax[3];
  int* ir1 = NULL, * ir2 = NULL, * ic1 = NULL, * ic2 = NULL, * ip1 = NULL, * ip2 = NULL;
  int* irowp = NULL, * icolp = NULL, * iplanep = NULL;
  int rowdirection, coldirection;
  int rowmax, rowmin, colmax, colmin, planemin, planemax;
  int i, j, line;
  //int ix1, iy1, ix2, iy2;
  int irow, icol, iplane, rh = 0;
  int iy = 0, ix = 0, iz = 0, nxw, nyw, nzw;
  int nrow, ncol, nsec;
  int ut1, ut2, ut3, lt1, lt2, lt3, upper;
  //int ilevel;
  int n = 0, npass = 0, rhomin = 99999, rhomax = -9999;
  PLINE e;

  float *mp=0;
  if (map_points.size())
    mp=&map_points[0];

  int plane_limit, row_limit, col_limit;

  nxw = mapheader->nx;
  nyw = mapheader->ny;
  nzw = mapheader->nz;
  int nx = mapheader->nx;
  int ny = mapheader->ny;
  int nz = mapheader->nz;

  int nxyz=nx*ny*nz;

#ifdef USE_ORTHGRID
  // this code is not currently used by MIFit, but it presumably could be.
  // I'm ifdefing it out for performance reasons.
  if (orthgrid) {
    nxw = (int) (mapheader->a / gridspacing);
    factor = ((float)nxw) / ((float)nx);
    xmin = (int)((float) xmin * factor);
    xmax = (int)((float) xmax * factor);
    nyw = (int) (mapheader->b / gridspacing);
    factor = ((float)nyw) / ((float)ny);
    ymin = (int)((float) ymin * factor);
    ymax = (int)((float) ymax * factor);
    nzw = (int) (mapheader->c / gridspacing);
    factor = ((float)nzw) / ((float)nz);
    zmin = (int)((float) zmin * factor);
    zmax = (int)((float) zmax * factor);
  }
#endif

  if (planedirection == X) {
    planemin = xmin;
    planemax = xmax;
    rowdirection = Y;
    rowmin = ymin;
    rowmax = ymax;
    coldirection = Z;
    colmin = zmin;
    colmax = zmax;

  } else if (planedirection == Y) {
    planemin = ymin;
    planemax = ymax;
    rowdirection = Z;
    rowmin = zmin;
    rowmax = zmax;
    coldirection = X;
    colmin = xmin;
    colmax = xmax;

  } else {
    planemin = zmin;
    planemax = zmax;
    rowdirection = X;
    rowmin = xmin;
    rowmax = xmax;
    coldirection = Y;
    colmin = ymin;
    colmax = ymax;
  }
  ncol = colmax - colmin + 1;
  nrow = rowmax - rowmin + 1;
  nsec = ncol * nrow;
  if (nsec <= 0) {
    if (flog) {
      fprintf(flog, "BUG - nsec<=0 nsec=%d ncol=%d nrow=%d\n",
        nsec, ncol, nrow);
    }
    return (-1);
  }
  if (nsec > (int)section.size()) {
    section.resize(nsec);
    if (nsec > (int)section.size()) {
      Logger::message("Not enough memory! Request smaller area");
      return (-1);
    }
  }
  sec = &section[0];
  switch (rowdirection) {
    case X:
      rint = 1.0f/(float)nxw;
      irowp = &ix;
      row1p = &x1;
      row2p = &x2;
      row_limit=nxw;
      break;
    case Y:
      rint = 1.0f/(float)nyw;
      irowp = &iy;
      row1p = &y1;
      row2p = &y2;
      row_limit=nyw;
      break;
    case Z:
      rint = 1.0f/(float)nzw;
      irowp = &iz;
      row1p = &z1;
      row2p = &z2;
      row_limit=nzw;
      break;
  }
  switch (coldirection) {
    case X:
      cint = 1.0f/(float)nxw;
      icolp = &ix;
      col1p = &x1;
      col2p = &x2;
      col_limit=nxw;
      break;
    case Y:
      cint = 1.0f/(float)nyw;
      icolp = &iy;
      col1p = &y1;
      col2p = &y2;
      col_limit=nyw;
      break;
    case Z:
      cint = 1.0f/(float)nzw;
      icolp = &iz;
      col1p = &z1;
      col2p = &z2;
      col_limit=nzw;
      break;
  }
  switch (planedirection) {
    case X:
      pint = 1.0f/(float)nxw;
      iplanep = &ix;
      plane1p = &x1;
      plane2p = &x2;
      plane_limit=nxw;
      break;
    case Y:
      pint = 1.0f/(float)nyw;
      iplanep = &iy;
      plane1p = &y1;
      plane2p = &y2;
      plane_limit=nyw;
      break;
    case Z:
      pint = 1.0f/(float)nzw;
      iplanep = &iz;
      plane1p = &z1;
      plane2p = &z2;
      plane_limit=nzw;
      break;
  }

#define do_edge() { \
    if ((n = segment(&r1, &c1, &r2, &c2, irow, icol, sec, nrow, upper, level)) == 2) { \
               r1 = (r1+rowmin)*rint; r2 = (r2+rowmin)*rint; c1 = (c1+colmin)*cint; c2 = (c2+colmin)*cint; \
               *plane1p = *plane2p = z; *row1p = r1; *row2p = r2; *col1p = c1; *col2p = c2; \
               e.p1.x = (float)(mapheader->ftoc[0][0]*x1 + mapheader->ftoc[0][1]*y1 + mapheader->ftoc[0][2]*z1); \
               e.p2.x = (float)(mapheader->ftoc[0][0]*x2 + mapheader->ftoc[0][1]*y2 + mapheader->ftoc[0][2]*z2); \
               e.p1.y = (float)(mapheader->ftoc[1][0]*x1 + mapheader->ftoc[1][1]*y1 + mapheader->ftoc[1][2]*z1); \
               e.p2.y = (float)(mapheader->ftoc[1][0]*x2 + mapheader->ftoc[1][1]*y2 + mapheader->ftoc[1][2]*z2); \
               e.p1.z = (float)(mapheader->ftoc[2][0]*x1 + mapheader->ftoc[2][1]*y1 + mapheader->ftoc[2][2]*z1); \
               e.p2.z = (float)(mapheader->ftoc[2][0]*x2 + mapheader->ftoc[2][1]*y2 + mapheader->ftoc[2][2]*z2); \
               e.p1.color = e.p2.color = color; edges.push_back(e); npass++; \
               } else { if (flog){fprintf(flog, "error in segment at %d %d %d  returned %d\n", irow, icol, upper, n);}}}

#define set_ixyz() { *iplanep = iplane; *icolp = icol + colmin; *irowp = irow + rowmin; }

#define get_cart_point() { \
  fx = (float)ix/(float)nxw; fy = (float)iy/(float)nyw; fz = (float)iz/(float)nzw; \
  p[0] = fx*mapheader->ftoc[X][X]+fy*mapheader->ftoc[X][Y]+fz*mapheader->ftoc[X][Z]; \
  p[1] = fx*mapheader->ftoc[Y][X]+fy*mapheader->ftoc[Y][Y]+fz*mapheader->ftoc[Y][Z]; \
  p[2] = fx*mapheader->ftoc[Z][X]+fy*mapheader->ftoc[Z][Y]+fz*mapheader->ftoc[Z][Z]; \
}

  /* draw graph */
  unsigned int nCurrentAtoms = 0;
  if (CurrentAtoms) {
    nCurrentAtoms = (*CurrentAtoms).size();
  }
  if (settings->ContourMethod == MAP_SPHERE) {
    InclusionCenter[0] = center[0];
    InclusionCenter[1] = center[1];
    InclusionCenter[2] = center[2];
    InclusionRadSq = settings->Radius * settings->Radius;
  } else if (settings->ContourMethod == MAP_BBLOB && nCurrentAtoms > 0) {
    InclusionRadSq = settings->BlobRadius * settings->BlobRadius;
    ic1 = (int*)malloc(6 * nCurrentAtoms * sizeof(int));
    ic2 = ic1 + nCurrentAtoms;
    ir1 = ic2 + nCurrentAtoms;
    ir2 = ir1 + nCurrentAtoms;
    ip1 = ir2 + nCurrentAtoms;
    ip2 = ip1 + nCurrentAtoms;
    for (i = 0; (unsigned int)i < nCurrentAtoms; i++) {
      aax[X] = (*CurrentAtoms)[i]->x();
      aax[Y] = (*CurrentAtoms)[i]->y();
      aax[Z] = (*CurrentAtoms)[i]->z();
      for (j = 0; j < 3; j++) {
        c[0][j] = aax[j] - settings->BlobRadius;
        c[1][j] = aax[j] + settings->BlobRadius;
      }
      for (j = 0; j < 2; j++) {
        mapheader->CtoF(&c[j][0], &c[j][1], &c[j][2]);
        c[j][0] *= (float)nxw;
        c[j][1] *= (float)nyw;
        c[j][2] *= (float)nzw;
      }
      ic1[i] = ROUND(std::max(c[0][coldirection], (float)colmin) - colmin);
      ic2[i] = ROUND(std::min(c[1][coldirection], (float)colmax) - colmin);
      ir1[i] = ROUND(std::max(c[0][rowdirection], (float)rowmin) - rowmin);
      ir2[i] = ROUND(std::min(c[1][rowdirection], (float)rowmax) - rowmin);
      ip1[i] = ROUND(std::max(c[0][planedirection], (float)planemin));
      ip2[i] = ROUND(std::min(c[1][planedirection], (float)planemax));
    }
  }

  for (iplane = planemin; iplane <= planemax; iplane++) {
    z = iplane * pint;

    *iplanep = (iplane + 64*plane_limit) % plane_limit;
    
    if (settings->ContourMethod == MAP_SPHERE || (settings->ContourMethod == MAP_BBLOB && nCurrentAtoms > 0)) {
      for (icol = 0; icol < ncol; icol++) {
        for (irow = 0; irow < nrow; irow++) {
          sec[nrow*icol+irow] = 0;
        }
      }
      if ((settings->ContourMethod == MAP_BBLOB && nCurrentAtoms > 0)) {
        for (i = 0; (unsigned int)i < nCurrentAtoms; i++) {
          if ((iplane < ip1[i]) || (iplane > ip2[i])) {
            continue;
          }
          InclusionCenter[X] = (*CurrentAtoms)[i]->x();
          InclusionCenter[Y] = (*CurrentAtoms)[i]->y();
          InclusionCenter[Z] = (*CurrentAtoms)[i]->z();
          for (icol = ic1[i]; icol <= ic2[i]; icol++) {
            for (irow = ir1[i]; irow <= ir2[i]; irow++) {
              set_ixyz();
              get_cart_point();
              if (inside_sphere(p)) {
                sec[nrow*icol+irow] = 1;
              }
            }
          }
        }
      }
    }

    for (icol = 0; icol < ncol; icol++) {
      *icolp = ((icol + colmin)+64*col_limit) % col_limit;

      for (irow = 0; irow < nrow; irow++) {
        *irowp = ((irow + rowmin)+64*row_limit) % row_limit;

      /* interchange irow,icol,iplane
         * depending upon plane to be drawn
         */
        //set_ixyz();
        /*  fill section */
        if (settings->ContourMethod == MAP_SPHERE) {
          // Unfortunately, this method can't take advantage of the
          // performance optimizations we put in for the other methods, so
          // we'll have to recalculate ix, iy, and iz here, check the
          // cartesian coordinates for the point, and *then* adjust ix/y/z
          // to be inside the cell.
          set_ixyz();
          get_cart_point();
          if (!inside_sphere(p)) {
            continue;
          }
          ix = (ix+64*nxw)%nxw;
          iy = (iy+64*nyw)%nyw;
          iz = (iz+64*nzw)%nzw;
        } else if (settings->ContourMethod == MAP_BBLOB && nCurrentAtoms > 0) {
          if (sec[nrow*icol+irow] == 0) {
            continue;
          }
          ix = (ix+64*nxw)%nxw;
          iy = (iy+64*nyw)%nyw;
          iz = (iz+64*nzw)%nzw;
        }
        if (UseNCR && (mapheader->nNCRSymmops > 0)) {
          set_ixyz();
          rh = NCRAverage(ix, iy, iz);
#ifdef USE_ORTHGRID
  // this code is not currently used by MIFit, but it presumably could be.
  // I'm ifdefing it out for performance reasons.
        } else if (orthgrid) {
          fx = (float)ix/(float)nxw;
          fy = (float)iy/(float)nyw;
          fz = (float)iz/(float)nzw;
          rh = ROUND(avgrho(fx, fy, fz));
#endif
        } else {

#ifdef OLDWAY
          // note: this is suprisingly expensive (integer division), and
          // somewhat redundant as only the value for the inner loop
          // needed to change here, so I replaced this with *irowp, *icolp, and
          // *iplanep calculations inside the relevant loops.
          ix = (ix+64*nxw)%nxw;
          iy = (iy+64*nyw)%nyw;
          iz = (iz+64*nzw)%nzw;
#endif
          i = nxw*(nyw*iz + iy) + ix;
          if (i >= 0 && i < nxyz) {
            rh = ROUND(mp[i]);
          } else {
            if (flog) {
              fprintf(flog, "%d out of range %d %d %d\n", i, ix, iy, iz);
            }
          }
        }
        sec[nrow*icol+irow] = (float)rh;
        if (rh > rhomax) {
          rhomax = rh;
        }
        if (rh < rhomin) {
          rhomin = rh;
        }
      }
    }

    if (level < rhomin || level > rhomax) {
      continue;
    }
    if (level < 0) {
      line = 3;                 /* dash negative contours */
    }
    for (icol = 0; icol < ncol-1; icol++) {
      for (irow = 0; irow < nrow-1; irow++) {
        upper = 0;
        ut1 =       (int) sec[irow+  nrow* icol   ];
        ut2 = lt1 = (int) sec[irow+1+nrow* icol   ];
        lt2 = (int) sec[irow+1+nrow*(icol+1)];
        ut3 = lt3 = (int) sec[irow+  nrow*(icol+1)];
        if (any_over(ut1, ut2, ut3, level) &&  any_under(ut1, ut2, ut3, level)) {
          upper = 2;
          do_edge();
        }
        if (any_over(lt1, lt2, lt3, level) && any_under(lt1, lt2, lt3, level)) {
          upper = 1;
          do_edge();
        }
      }
    }
  }   /* iplane loop ends here */
  if (settings->ContourMethod == MAP_BBLOB && nCurrentAtoms > 0) {
    free(ic1);
  }
  //UnLockMap();
  return (npass);
}

/*
 *  returns the x and y coordinates relative to 0,0
 *  for the line segment for a contour triangle
 *  using linear interpolation.
 *  do not call this if you are not sure that the
 *  triangle contains a contour line!
 *  function returns 2 if successful.
 * 0 and 1 are errors and a line should not be drawn!
 */
int EMapBase::segment(double* x1, double* y1, double* x2, double* y2, int i, int j,
                      float* sec, int nrow, int upper, int level) {
  int found = 0;
  static float n1, n2, n3, n4, t;
  /* switch on upper or lower triangle
   * then search through three sides
   * first side found gets x1,y1
   * second side gets x2,y2
   * then return
   *
   * if upper = 2 do upper
   * if upper = 1 do lower
   * if upper = 0  or > 2 then do nothing
   *
   *    A box starting at i,j split into two tringles thus:
   *
   *          1.(i,j)-----2.(i+1,j)
   *           |        _/|
   *           | upper_/  |
   *           |    _/    |
   *           |  _/ lower|
   *           | /        |
   *          3.(i,j+1)---4.(i+1,j+1)
   */
  if (upper == 2) {
    /* upper triangle */
    n1 = sec[i+nrow*j]-level;
    n2 = sec[i+1+nrow*j]-level;
    n3 = sec[i+nrow*(j+1)]-level;
    if (n1 == 0) {
      n1 += 1;
    }
    if (n2 == 0) {
      n2 += 1;
    }
    if (n3 == 0) {
      n3 += 1;
    }
    /*  side 1,2 */
    if (n1 * n2 <= 0) {    /* contour crosses side 1,2 */
      found++;
      /* interpolate x,y */
      *x1 = (float)i+(0-n1)/(n2-n1);
      *y1 = j;
    }
    /*  side 2,3 */
    if (n2 * n3 <= 0) {
      found++;
      if (found == 1) {
        t = (0-n2)/(n3-n2);
        *x1 = (float)i+ 1.0 -t;
        *y1 = (float)j + t;
      } else {
        t = (0-n2)/(n3-n2);
        *x2 = (float)i+ 1.0 -t;
        *y2 = (float)j + t;
        return (found);
      }
    }
    /*  side 1,3  if we got this far it must number 2 */
    if (n1 * n3 <= 0) {
      found++;
      *x2 = i;
      *y2 = (float)j + (0-n1)/(n3-n1);
    }
    return (found);
  }
  if (upper == 1) {
    /* lower triangle */
    n4 = sec[i+1+nrow*(j+1)]-level;
    n2 = sec[i+1+nrow*j]-level;
    n3 = sec[i+nrow*(j+1)]-level;
    if (n4 == 0) {
      n4 += 1;
    }
    if (n2 == 0) {
      n2 += 1;
    }
    if (n3 == 0) {
      n3 += 1;
    }
    /*  side 2,4 */
    if (n2 * n4 <= 0) {
      found++;
      /* interpolate x,y */
      *x1 = i+1.0;
      *y1 = (float)j + (0-n2)/(n4-n2);
    }
    /*  side 3,4 */
    if (n3 * n4 <= 0) {
      found++;
      if (found == 1) {
        *x1 = (float)i + (0-n3)/(n4-n3);
        *y1 = j + 1.0;
      } else {
        *x2 = (float)i + (0-n3)/(n4-n3);
        *y2 = j + 1.0;
        return (found);
      }
    }
    /*  side 2,3 */
    if (n2 * n3 <= 0) {
      found++;
      t = (0-n2)/(n3-n2);
      *x2 = (float)i+ 1.0 -t;
      *y2 = (float)j + t;
    }
    return (found);
  }
  return (0);
}

void CMapHeaderBase::NCRTransform(float* fx, float* fy, float* fz, int is) {
  /* transform to cartesian */
  float cx, cy, cz, tx, ty, tz;
  cx = *fx*ftoc[X][X]+*fy*ftoc[X][Y]+*fz*ftoc[X][Z];
  cy = *fx*ftoc[Y][X]+*fy*ftoc[Y][Y]+*fz*ftoc[Y][Z];
  cz = *fx*ftoc[Z][X]+*fy*ftoc[Z][Y]+*fz*ftoc[Z][Z];
  /* transform by symmop */
  tx = cx*NCRSymmops[is][0]+cy*NCRSymmops[is][1]+cz*NCRSymmops[is][2]+ NCRSymmops[is][9];
  ty = cx*NCRSymmops[is][3]+cy*NCRSymmops[is][4]+cz*NCRSymmops[is][5]+ NCRSymmops[is][10];
  tz = cx*NCRSymmops[is][6]+cy*NCRSymmops[is][7]+cz*NCRSymmops[is][8]+ NCRSymmops[is][11];
  /* transform back to fractional */
  *fx = tx*ctof[X][X]+ty*ctof[X][Y]+tz*ctof[X][Z];
  *fy = tx*ctof[Y][X]+ty*ctof[Y][Y]+tz*ctof[Y][Z];
  *fz = tx*ctof[Z][X]+ty*ctof[Z][Y]+tz*ctof[Z][Z];
}

int EMapBase::NCRAverage(int ix, int iy, int iz) {
  float rh = 0;
  int n = 0;
  int is;
  float fx, fy, fz;
  /* tranform ix, iy, iz to fractional */
  fx = (float)ix/(float)mapheader->nx;
  fy = (float)iy/(float)mapheader->ny;
  fz = (float)iz/(float)mapheader->nz;
  for (is = 0; is < mapheader->nNCRSymmops; is++) {
    mapheader->NCRTransform(&fx, &fy, &fz, is);
    /* get rho at this point by interpolation */
    rh += avgrho(fx, fy, fz);
    n++;
  }
  if (n > 0) {
    rh /= (float)n;
  }
  return ROUND(rh);
}

