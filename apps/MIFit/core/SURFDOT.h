#ifndef mifit_model_SURFDOT_h
#define mifit_model_SURFDOT_h

#include <chemlib/chemlib.h>

#ifdef TESTING
extern float spacing;
#endif

/**
 * A dot surface point.
 */
class SURFDOT {
public:
  float x;
  float y;
  float z;
  short color;
  unsigned short w;
};


//@{
// surface atom a in the context of atoms b[].
// b should be the atoms who are potential neighbours.
//@}
long atomsurf(chemlib::MIAtom* a, float ra, chemlib::MIAtom** b, long nb,
              SURFDOT** dots, long ndots, long* maxdots, void*&, float dotsper, float radius_offset = 0.0f);
//@{
// Surface an atom at radius r.
//@}
long atomsurfradius(chemlib::MIAtom* a, float r, SURFDOT** dots,
                    long ndots, long* maxdots, void*& hglb, float dotsper);
//@{
// Called before a surface calculation.
//@}
int initspheres(float* radii);
//@{
// For clearing spheres if they need to be recalculated.
//@}
void clearspheres();

#endif
