#include "nonguilib.h"
#include "chemlib.h"

#include "SURFDOT.h"
#include "SurfaceSphere.h"
#include "ViewPoint.h"

//#include "Xguicryst.h"

using namespace chemlib;


int SurfResult = 1;

SurfaceSphere surfaceSpheres[NTYPES];
float spacing = -1.0F;

int initspheres(float* radii) {
  SurfResult = 1;
  for (unsigned int i = 0; i < NTYPES; i++) {
    SurfaceSphere& sphere = surfaceSpheres[i];
    if (sphere.getRadius() != radii[i] || sphere.getSpacing() != spacing) {
      sphere.build(radii[i], spacing);
    }
  }
  return 1;
}

void clearspheres() {
  for (unsigned int i = 0; i < NTYPES; i++) {
    SurfaceSphere& sphere = surfaceSpheres[i];
    sphere.clearPoints();
  }
  spacing = -1.0F;
}

//FIXME: use better memory memory management; something other than void*

/* surface atom a in the context of atoms b[] */
/* b should be the atoms who are potential neighbours */

long
atomsurf(MIAtom* a, float ra, MIAtom** b, long nb, SURFDOT** dots,
         long ndots, long* maxdots, void*& hglb, float dotsper, float radius_offset) {
  unsigned int i;
  long j;
  MIAtom* a2;
  int itype = 0, jtype;
  float x1, y1, z1, dx, dy, dz;
  float r1, r2, dist;
  int obscured = 0;
  size_t nbytes;
  SurfResult = 1;
  spacing = dotsper;

  float radius[NTYPES];
  for (unsigned int i = 0; i < NTYPES; ++i) {
    radius[i] = MIAtom::MIAtomRadiusForType(i)+radius_offset;
  }
  initspheres(radius);

  if (*maxdots == 0) {
    *maxdots = 300;
    nbytes = *maxdots*sizeof(SURFDOT);
    hglb = malloc(nbytes);
    if (hglb == NULL) {
      SurfResult = 0;
      Logger::message("Out of memory in surface routine - Increase the Angstroms/dot and try again.");
      return (ndots);
    }
    *dots = (SURFDOT*) hglb;
  }
  itype = MIAtom::MIGetAtomTypeFromName(a->name());
  r1 = radius[itype]*ra;

  std::vector<APOINT>& points = surfaceSpheres[itype].getPoints();
  for (i = 0; i < points.size(); i++) {
    x1 = a->x();
    x1 = x1 + points[i].x*ra;
    y1 = a->y();
    y1 = y1 + points[i].y*ra;
    z1 = a->z();
    z1 = z1 + points[i].z*ra;
    obscured = 0;
    for (j = 0; j < nb; j++) {
      a2 = b[j];
      if (a2 != a) {
        jtype = MIAtom::MIGetAtomTypeFromName(a2->name());
        r2 = radius[jtype]*ra;
        r2 = r2 *r2;
        dx = a2->x();
        dx = dx - x1;
        dist = dx*dx;
        if (dist < r2) {
          dy = a2->y();
          dy = dy - y1;
          dist += dy*dy;
          if (dist < r2) {
            dz = a2->z();
            dz = dz - z1;
            dist += dz*dz;
            if (dist < r2) {
              obscured = 1;
              break;
            }
          }
        }
      }
    }
    if (!obscured) {
      /* output this point */
      if (ndots >= *maxdots) {
        *maxdots += 500;
        nbytes = *maxdots*sizeof(SURFDOT);
        void* hnew = malloc(nbytes);
        if (!hnew) {
          SurfResult = 0;
          return ndots;
        }
        void* p = hnew;
        if (p == NULL) {
          SurfResult = 0;
          return ndots;
        }
        memcpy(p, *dots, ndots*sizeof(SURFDOT));
        if (hglb != NULL) {
          free(hglb);
          hglb = NULL;
        }
        hglb = hnew;
        *dots = (SURFDOT*)p;
      }
      (*dots)[ndots].x = x1;
      (*dots)[ndots].y = y1;
      (*dots)[ndots].z = z1;
      (*dots)[ndots].w = 1;
      (*dots)[ndots].color = abs(a->color());
      ndots++;
    }
  }
  return (ndots);
}

long
atomsurfradius(MIAtom* a, float r, SURFDOT** dots,
               long ndots, long* maxdots, void*& hglb, float dotsper) {
  unsigned int i;
  float x1, y1, z1;
  int obscured = 0;
  size_t nbytes;
  SurfResult = 1;
  WaitCursor wait("Surface Atoms");
  SurfaceSphere sphere;
  sphere.build(r, dotsper);
  std::vector<APOINT>& points = sphere.getPoints();
  for (i = 0; i < points.size(); i++) {
    x1 = (float)a->x();
    x1 = x1 + (float) points[i].x;
    y1 = (float)a->y();
    y1 = y1 + (float)points[i].y;
    z1 = (float)a->z();
    z1 = z1 + (float)points[i].z;
    obscured = 0;
    if (!obscured) {
      /* output this point */
      if (ndots >= *maxdots) {
        *maxdots += 500;
        nbytes = *maxdots*sizeof(SURFDOT);
        void* hnew = malloc(nbytes);
        if (!hnew) {
          SurfResult = 0;
          return ndots;
        }
        void* p = hnew;
        if (p == NULL) {
          SurfResult = 0;
          return ndots;
        }
        memcpy(p, *dots, ndots*sizeof(SURFDOT));
        if (hglb != NULL) {
          free(hglb);
          hglb = NULL;
        }
        hglb = hnew;
        *dots = (SURFDOT*)p;
      }
      (*dots)[ndots].x = x1;
      (*dots)[ndots].y = y1;
      (*dots)[ndots].z = z1;
      (*dots)[ndots].w = 1;
      (*dots)[ndots].color = abs(a->color());
      ndots++;
      if (ndots%100 == 0) {
        if (wait.CheckForAbort() == true) {
          return (ndots);
        }
      }
    }
  }
  return (ndots);
}

