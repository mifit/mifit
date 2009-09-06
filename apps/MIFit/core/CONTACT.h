#ifndef MIFIT_MODEL_CONTACT_H_
#define MIFIT_MODEL_CONTACT_H_

#include <chemlib/chemlib.h>

/**
 * A contact between two atoms. A contact is drawn as a line with a distance printed at the middle..
 */
class CONTACT {
public:
  chemlib::Bond line;
  float d;
  short color;
};

#endif /*MIFIT_MODEL_CONTACT_H_*/
