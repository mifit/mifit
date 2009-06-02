#define NOMINMAX
#include "nonguilib.h"
#include "chemlib.h"
#include "RESIDUE_.h"
#include "utillib.h"
#include "maplib.h"
#include "mathlib.h"

#include "RESIDUE.h"
#include "Cfiles.h"
#include "rotlsq.h"


using namespace std;
using namespace chemlib;




#define X 0
#define Y 1
#define Z 2

int read_colors(const RESIDUE* res, char* buf, int nbuf) {
  char name[100];
  int i = 0;
  if (sscanf(buf, "res %s", name) != 1) {
    return 0;
  }
  if (strcmp(res->name().c_str(), name)) {
    Logger::message("Error: out of sync in reading colors");
  }
  if (strstr(buf, "_") == NULL) {
    // old way for compatibility
    int sign;
    int nc = 4 + strlen(name);
    while (nc < nbuf && buf[nc] != '\0') {
      if (!isspace(buf[nc])) {
        if (buf[nc] == '-') {
          sign = -1; nc++;
        } else {sign = 1;}
        if (i < res->atomCount()) {
          res->atom(i)->setColor((buf[nc] - '0')*sign);
          if (res->atom(i)->color() == 0) {
            res->atom(i)->setColor(10);
          }
          i++;
        } else {break;}
      }
      nc++;
    }
  } else {
    char* pos = strstr(buf, "_");
    int c;
    while (pos != NULL && pos < buf+strlen(buf)-1 && i < res->atomCount()) {
      if (sscanf(pos+1, "%d", &c) == 1) {
        res->atom(i)->setColor(c);
      }
      i++;
      pos = strstr(pos+1, "_");
    }
  }
  return i;
}

int read_radii(const RESIDUE* res, char* buf, int nbuf) {
  char name[20];
  int i = 0;
  if (sscanf(buf, "res %s", name) != 1) {
    return 0;
  }
  if (strcmp(res->name().c_str(), name)) {
    Logger::message("Error: out of sync in reading radii");
  }
  int sign;
  int nc = 4 + strlen(name);
  while (nc < nbuf && buf[nc] != '\0') {
    if (!isspace(buf[nc])) {
      if (buf[nc] == '-') {
        sign = -1; nc++;
      } else {sign = 1;}
      if (i < res->atomCount()) {
        res->atom(i)->set_radius_type((buf[nc] - '0')*sign);
        i++;
      } else {break;}
    }
    nc++;
  }
  return i;
}

string resid(const RESIDUE* res) {
  string id;
  if (res != NULL) {
    char c = (char)(res->chain_id()&255);
    if (c == ' ') {
      id = format("%s %s", res->type().c_str(), res->name().c_str());
    } else {
      id = format("%s %s %c", res->type().c_str(), res->name().c_str(), c);
    }
  }
  return id;
}

static char one[] = {'A', 'C', 'D', 'E', 'F',
                     'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R',
                     'S', 'T', 'V', 'W', 'Y',
                     'O', 'O', 'X', 'U', 'Z',
                     'M', 'Y', 'S', 'T', 'C', 'C'};
static char charge[] = {' ', ' ', '-', '-', ' ',
                        ' ', '+', ' ', '+', ' ',
                        ' ', ' ', ' ', ' ', '+',
                        ' ', ' ', ' ', ' ', ' ',
                        ' ', ' ', '+', '+', '+',
                        ' ', ' ', ' ', ' ', ' ',
                        ' ', ' ', ' ', ' ', ' ', ' '};

char chargetype(char t) {
  // convert one letter residue names to single letter charge type
  for (unsigned int i = 0; i < sizeof(one); i++) {
    if (t == one[i]) {
      return (charge[i]);
    }
  }
  return (' ');
}

float
phi(const RESIDUE* prev, const RESIDUE* res) {
  MIAtom* a1 = NULL;
  MIAtom* a2 = NULL;
  MIAtom* a3 = NULL;
  MIAtom* a4 = NULL;
  MIAtom* a = NULL;
  int i;
  for (i = 0; i < res->atomCount(); i++) {
    a = res->atom(i);
    if (!strcmp(a->name(), "N")) {
      a2 = a;
    }
    if (!strcmp(a->name(), "CA")) {
      a3 = a;
    }
    if (!strcmp(a->name(), "C")) {
      a4 = a;
    }
  }
  for (i = 0; i < prev->atomCount(); i++) {
    if (!strcmp(prev->atom(i)->name(), "C")) {
      a1 = prev->atom(i);
    }
  }
  if (a1 && a2 && a3 && a4) {
    return (CalcAtomTorsion(a1, a2, a3, a4));
  }
  return (-360.0);
}

float
psi(const RESIDUE* res, const RESIDUE* next) {
  MIAtom* a1 = NULL;
  MIAtom* a2 = NULL;
  MIAtom* a3 = NULL;
  MIAtom* a4 = NULL;
  MIAtom* a = NULL;
  int i;
  for (i = 0; i < res->atomCount(); i++) {
    a = res->atom(i);

    if (!strcmp(a->name(), "N")) {
      a1 = a;
    }
    if (!strcmp(a->name(), "CA")) {
      a2 = a;
    }
    if (!strcmp(a->name(), "C")) {
      a3 = a;
    }
  }
  for (i = 0; i < next->atomCount(); i++) {
    if (!strcmp(next->atom(i)->name(), "N")) {
      a4 = next->atom(i);
    }
  }
  if (a1 && a2 && a3 && a4) {
    return (CalcAtomTorsion(a1, a2, a3, a4));
  }
  return (-360.0);
}

bool MoveOnto(const RESIDUE* res, RESIDUE* fitres, int nres) {
  // move res onto fitres
  int ma = 0, mb = 0, m = 0;
  int i;
  double a[4][3], b[4][3];
  double r[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  double v[3];
  double w[4];
  double x, y, z;
  double tx, ty, tz;
  MIAtom* at;
  for (i = 0; i < 4; i++) {
    w[i] = 1.0;
  }
  if (IsPeptide(*res) && IsPeptide(*fitres)) {
    for (i = 0; i < res->atomCount(); i++) {
      if (!strcmp("N", res->atom(i)->name())) {
        ma++; res->atom(i)->getPosition(a[0]);
      }
      if (!strcmp("C", res->atom(i)->name())) {
        ma++; res->atom(i)->getPosition(a[1]);
      }
      if (!strcmp("CA", res->atom(i)->name())) {
        ma++; res->atom(i)->getPosition(a[2]);
      }
      if (!strcmp("CB", res->atom(i)->name())) {
        ma++; res->atom(i)->getPosition(a[3]);
      }
    }
    for (i = 0; i < fitres->atomCount(); i++) {
      if (!strcmp("N", fitres->atom(i)->name())) {
        mb++; fitres->atom(i)->getPosition(b[0]);
      }
      if (!strcmp("C", fitres->atom(i)->name())) {
        mb++; fitres->atom(i)->getPosition(b[1]);
      }
      if (!strcmp("CA", fitres->atom(i)->name())) {
        mb++; fitres->atom(i)->getPosition(b[2]);
      }
      if (!strcmp("CB", fitres->atom(i)->name())) {
        mb++; fitres->atom(i)->getPosition(b[3]);
      }
    }
    if (mb == 3 || ma == 3) { // glycine add O to matrix to make stable
      if ((at = atom_from_name("O", *fitres))) {
        mb++;
        at->getPosition(b[3]);
      }
      if ((at = atom_from_name("O", *res))) {
        ma++;
        at->getPosition(a[3]);
      }
      w[3] = 0.5; //downweight O
    }
    m = std::min(ma, mb);
  }

  if (m==0){ // nothing to superimpose!
	  return false;  
  }

  // find the rotation matrix and translation vector
  if (!rotlsqfit(a, b, w, m, r, v)) {
    // translate so that first atoms match
    v[0] = b[0][0] - a[0][0];
    v[1] = b[0][1] - a[0][1];
    v[2] = b[0][2] - a[0][2];
  }
  orthomatrix(r, r);
  int n = 0;
  while ((res != NULL) && n < nres) {
    for (i = 0; i < res->atomCount(); i++) {
      tx = res->atom(i)->x();
      ty = res->atom(i)->y();
      tz = res->atom(i)->z();
      x = r[0][0]*tx+r[0][1]*ty+r[0][2]*tz
          +v[0];
      y = r[1][0]*tx+r[1][1]*ty+r[1][2]*tz
          +v[1];
      z = r[2][0]*tx+r[2][1]*ty+r[2][2]*tz
          +v[2];
      res->atom(i)->setPosition((float)x, (float)y, (float)z);
    }
    n++;
    res = res->next();
  }
  return true;
}

void symm(float x, float y, float z, float* xp, float* yp, float* zp, float mat[3][4]) {
  *xp = x * mat[X][X] + y * mat[X][Y] + z * mat[X][Z] + mat[X][3];
  *yp = x * mat[Y][X] + y * mat[Y][Y] + z * mat[Y][Z] + mat[Y][3];
  *zp = x * mat[Z][X] + y * mat[Z][Y] + z * mat[Z][Z] + mat[Z][3];
}

bool outside_sphere(float p[3], float radius, float InclusionCenter[3]) {
  float x, y, z;
  float InclusionRadSq = radius * radius;
  x = p[X] - InclusionCenter[X];
  x *= x;
  if (x <= InclusionRadSq) {
    y = p[Y] - InclusionCenter[Y];
    y *= y;
    if (y <= InclusionRadSq) {
      z = p[Z] - InclusionCenter[Z];
      z *= z;
      if (z <= InclusionRadSq) {
        return (x + y + z > InclusionRadSq);
      }
    }
  }
  return true;
}

RESIDUE* SymmResidue(const RESIDUE* Model, CMapHeaderBase* mh, float center[3], float r, int color) {
  /* set unity to true if operator = x,y,z to avoid looking at 0,0,0 position */
  const RESIDUE* model;
  RESIDUE* res = NULL, * newres, * Res = NULL;
  char label[MAXNAME];
  int nlabel;
  int nadd = 0;
  int i, j, k, n = 0;
  float minxyz[3], maxxyz[3];
  float sminxyz[3], smaxxyz[3];
  float fx, fy, fz;
  int ix, iy, iz;
  float x1, y1, z1;
  float cx1, cy1, cz1, cx2, cy2, cz2;
  float pos[3];
  float cx, cy, cz, dx, dy, dz;
  float symmat[3][4];
  float ctof[3][3];
  float ftoc[3][3];
  for (j = 0; j < 3; j++) {
    for (k = 0; k < 3; k++) {
      ctof[j][k] = mh->ctof[j][k];
      ftoc[j][k] = mh->ftoc[j][k];
    }
  }

  minxyz[0] = 9999.0;
  minxyz[1] = 9999.0;
  minxyz[2] = 9999.0;
  maxxyz[0] = -9999.0;
  maxxyz[1] = -9999.0;
  maxxyz[2] = -9999.0;
  cx = center[X];
  cy = center[Y];
  cz = center[Z];
  cx1 = center[X]-r;
  cy1 = center[Y]-r;
  cz1 = center[Z]-r;
  cx2 = center[X]+r;
  cy2 = center[Y]+r;
  cz2 = center[Z]+r;
  transform(ctof, &cx, &cy, &cz);
  transform(ctof, &cx1, &cy1, &cz1);
  transform(ctof, &cx2, &cy2, &cz2);

  for (int symmop = 0; symmop < mh->nsym; symmop++) {
    /* put res at end of list */
    sprintf(label, "#%d", symmop);
    nlabel = strlen(label);

    for (j = 0; j < 3; j++) {
      for (k = 0; k < 4; k++) {
        symmat[j][k] = mh->symops[j][k][symmop];
      }
    }


    /* find min max of model in fractional */
    model = Model;
    while (model != NULL) {
      for (i = 0; i < model->atomCount(); i++) {
        if (model->atom(i)->x() < minxyz[X]) {
          minxyz[X] = model->atom(i)->x();
        }
        if (model->atom(i)->y() < minxyz[Y]) {
          minxyz[Y] = model->atom(i)->y();
        }
        if (model->atom(i)->z() < minxyz[Z]) {
          minxyz[Z] = model->atom(i)->z();
        }
        if (model->atom(i)->x() > maxxyz[X]) {
          maxxyz[X] = model->atom(i)->x();
        }
        if (model->atom(i)->y() > maxxyz[Y]) {
          maxxyz[Y] = model->atom(i)->y();
        }
        if (model->atom(i)->z() > maxxyz[Z]) {
          maxxyz[Z] = model->atom(i)->z();
        }
      }
      model = model->next();
    }
    /* convert to fractional */
    transform(ctof, &minxyz[X], &minxyz[Y], &minxyz[Z]);
    transform(ctof, &maxxyz[X], &maxxyz[Y], &maxxyz[Z]);

    /* find symmetry min max */
    symm(minxyz[X], minxyz[Y], minxyz[Z],
      &sminxyz[X], &sminxyz[Y], &sminxyz[Z], symmat);
    symm(maxxyz[X], maxxyz[Y], maxxyz[Z],
      &smaxxyz[X], &smaxxyz[Y], &smaxxyz[Z], symmat);
    for (i = 0; i < 3; i++) {
      if (sminxyz[i] > smaxxyz[i]) {
        float t = sminxyz[i];
        sminxyz[i] = smaxxyz[i];
        smaxxyz[i] = t;
      }
    }

    /* look at unit cells -1 to +1 around minmax */
    dx = cx - (sminxyz[X]+smaxxyz[X])/2.0 ;
    dy = cy - (sminxyz[Y]+smaxxyz[Y])/2.0 ;
    dz = cz - (sminxyz[Z]+smaxxyz[Z])/2.0 ;
    for (ix = -1; ix <= 1; ix += 1) {
      for (iy = -1; iy <= 1; iy += 1) {
        for (iz = -1; iz <= 1; iz += 1) {
          /*    if (symmop==0 && ix==0 && iy==0 && iz==0) continue; */
          fx = ROUND(dx + ix);
          fy = ROUND(dy + iy);
          fz = ROUND(dz + iz);
          n++;
          model = Model;
          while (model != NULL) {
            /* look for any atom to be within radius-
               if so copy entire residue and move on */
            for (i = 0; i < model->atomCount(); i++) {
              x1 = model->atom(i)->x();
              y1 = model->atom(i)->y();
              z1 = model->atom(i)->z();
              transform(ctof, &x1, &y1, &z1);
              symm(x1, y1, z1, &x1, &y1, &z1, symmat);
              x1 += fx;
              y1 += fy;
              z1 += fz;
              transform(ftoc, &x1, &y1, &z1);
              pos[X] = x1;
              pos[Y] = y1;
              pos[Z] = z1;
              if (outside_sphere(pos, r, center)) {
                continue;
              }
              if (symmop == 0
                  && x1 < model->atom(i)->x()+.1
                  && x1 > model->atom(i)->x()-.1
                  && y1 < model->atom(i)->y()+.1
                  && y1 > model->atom(i)->y()-.1
                  && z1 < model->atom(i)->z()+.1
                  && z1 > model->atom(i)->z()-.1) {
                continue;
              }
              newres = new RESIDUE(*model);
              if (Res == NULL) {
                Res = res = newres;
              } else {
                res = res->insertResidue(newres);
              }
              /* now transform the atoms */
              if ((int)strlen(res->name().c_str()) < MAXNAME-nlabel) {
                res->setName(res->name() + std::string(label));
              }
              for (j = 0; j < res->atomCount(); j++) {
                x1 = res->atom(j)->x();
                y1 = res->atom(j)->y();
                z1 = res->atom(j)->z();
                transform(ctof, &x1, &y1, &z1);
                symm(x1, y1, z1, &x1, &y1, &z1, symmat);
                x1 += fx;
                y1 += fy;
                z1 += fz;
                transform(ftoc, &x1, &y1, &z1);
                res->atom(j)->setPosition(x1, y1, z1);
                res->atom(j)->setColor(color);
                res->atom(j)->setSymmop(symmop);
                res->atom(j)->addType(AtomType::SYMMATOM);
                /* save fx,fy,fz for later */
                res->atom(j)->resetDelta();
                res->atom(j)->addDelta(fx, fy, fz);
                nadd++;
              }
              break;
            }
            model = model->next();
            //      if(wait.CheckForAbort()) model=NULL;
          }
        }
      }
    }
  }
  Logger::log("Built %d symmatoms", nadd);
  return (Res);
}

#undef X
#undef Y
#undef Z




void getchain(unsigned short chain_id, RESIDUE* reslist, RESIDUE*& nter, RESIDUE*& cter) {
  nter = cter = NULL;
  RESIDUE* res = reslist;
  while ((res != NULL) && res->chain_id() != chain_id) {
    res = res->next();
  }
  nter = res;
  while ((res != NULL) && res->chain_id() == chain_id) {
    cter = res;
    res = res->next();
  }
}

int order_ends(RESIDUE*& res1, RESIDUE*& res2, RESIDUE* res) {
  int nres = 0;
  RESIDUE* start = NULL, * end = NULL;
  if (res1 == res2) {
    nres = 1;
    start = res1;
    end = res1;
  } else {
    res = res1;
    nres = 0;
    while (res != NULL) {
      nres++;
      if (res == res2) {
        start = res1;
        end = res2;
        break;
      }
      res = res->next();
    }
    if (!start) {
      res = res2;
      nres = 0;
      while (res != NULL) {
        nres++;
        if (res == res1) {
          start = res2;
          end = res1;
          break;
        }
        res = res->next();
      }
    }
  }
  res1 = start;
  res2 = end;
  return nres;
}

const string chainstring(const RESIDUE* res) {
  char chainid;
  int chainno;
  static string buff;
  chainid = (char)(res->chain_id() & 255);
  const RESIDUE* res2 = res;
  const RESIDUE* r = res;
  while (r != NULL) {
    if (r->chain_id() == res->chain_id()) {
      res2 = r;
    } else {break;}
    r = r->next();
  }
  if (chainid == ' ') {
    chainid = '_';
  }
  chainno = res->chain_id()/256;
  if (res->linkage_type()&PEPTIDE || res->linkage_type()&NUCLEIC) {
    buff = "Chain:";
  } else {
    buff = "Segmt:";
  }
  buff += chainid;
  buff += ".";
  buff += ftoa(chainno);
  buff += ":";
  buff += resid(res);
  buff += "->";
  buff += resid(res2);
  return (buff);
}

RESIDUE* make_res(const MIAtomList& atoms) {
  RESIDUE* res = new RESIDUE();
  res->setAtoms(atoms);
  return res;
}

