#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctype.h>

#include "nonguilib.h"
#include "mathlib.h"
#include "chemlib.h"
#include "RESIDUE_.h"

//all these headers are verified as required
#include "sfcalc.h"
#include "sfcalc_data.h"
#include "maplib.h"
#include "fft.h"


using namespace chemlib;

#define X 0
#define Y 1
#define Z 2

void
sfinit() {
  int i, j;
  static int init = 0;
  double fj, sol, sol2;

  if (init != 0) {
    return;
  }
  /* reciprocal space form factor table */
  ftable = (double**)malloc(MAXFTABLE*sizeof(double*));
  for (i = 0; i < MAXFTABLE; i++) {
    ftable[i] = (double*)malloc(100*sizeof(double));
    for (j = 0; j < 100; j++) {
      fj = j;
      sol = fj*0.005;
      sol2 = sol*sol;
      ftable[i][j] = a1[i]*exp(-b1[i]*sol2) +
                     a2[i]*exp(-b2[i]*sol2) +
                     co[i] +
                     a3[i]*exp(-b3[i]*sol2) +
                     a4[i]*exp(-b4[i]*sol2);
    }
  }
  init = 1;
  srand(1);
}

int sfatom(const MIAtom& atom, CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  static int ir, in, it;
  static double scftmp, sf,  ah, ak, al, phase;
  static double cp, sp;
  int type = 0;
  int nr = 0;
  static double* xpart = NULL, * ypart = NULL, * zpart = NULL, * trans = NULL;
  static double x, y, z, sthol;
  static int init = 0, index;
  static double twopi;
  static double degtor;
  static int nsym;
  static CREFL* oldrefl = NULL;
  static int oldnrefl = 0;

  if (atom.occ() < 0.00001) {
    return 0;
  }
  if (atom.type() & AtomType::DUMMYATOM) {
    return 0;
  }

  nsym = mh->nsym;
  if (init == 0) {
    sfinit();
    twopi = 2.0 * acos(-1.0);
    degtor = acos(-1.0)/180.0;
    init = 1;
  }
  if (oldrefl != &refl[0] || oldnrefl != nrefl) {
    Logger::log("Allocating and initializing reflection data");
    if (xpart != NULL) {
      free(xpart);
    }
    if (ypart != NULL) {
      free(ypart);
    }
    if (zpart != NULL) {
      free(zpart);
    }
    if (trans != NULL) {
      free(trans);
    }
    xpart = (double*)malloc(nrefl*nsym*sizeof(double));
    ypart = (double*)malloc(nrefl*nsym*sizeof(double));
    zpart = (double*)malloc(nrefl*nsym*sizeof(double));
    trans = (double*)malloc(nrefl*nsym*sizeof(double));
    if (trans == NULL) {
      Logger::log("SFCalc: cannot allocate memory!");
      return (0);
    }
    for (ir = 0; ir < nrefl; ir++) {
      ah = refl[ir].ind[0];
      ak = refl[ir].ind[1];
      al = refl[ir].ind[2];
      for (in = 0; in < nsym; in++) {
        trans[nsym*ir+in] = ah * mh->symops[X][3][in] +
                            ak * mh->symops[Y][3][in] +
                            al * mh->symops[Z][3][in];
        xpart[nsym*ir+in] = ah * mh->symops[X][X][in] +
                            ak * mh->symops[Y][X][in] +
                            al * mh->symops[Z][X][in];
        ypart[nsym*ir+in] = ah * mh->symops[X][Y][in] +
                            ak * mh->symops[Y][Y][in] +
                            al * mh->symops[Z][Y][in];
        zpart[nsym*ir+in] = ah * mh->symops[X][Z][in] +
                            ak * mh->symops[Y][Z][in] +
                            al * mh->symops[Z][Z][in];
      }
    }
    oldrefl = &refl[0];
    oldnrefl = nrefl;
  }

  type = ScattIndex(&atom.name()[0], "*");

  /*  convert from cartesian to fractional coordinates */
  x = atom.x()*mh->ctof[X][X]+atom.y()*mh->ctof[X][Y]+atom.z()*mh->ctof[X][Z];
  y = atom.x()*mh->ctof[Y][X]+atom.y()*mh->ctof[Y][Y]+atom.z()*mh->ctof[Y][Z];
  z = atom.x()*mh->ctof[Z][X]+atom.y()*mh->ctof[Z][Y]+atom.z()*mh->ctof[Z][Z];

  for (ir = 0; ir < nrefl; ir++) {
    sthol = refl[ir].sthol;
    if (sthol > 0.5/mh->resmin
        || sthol < 0.5/mh->resmax) {
      continue;
    }
    nr++;
    it = (int)(sthol*200. + 0.5);
    if (it > 99) {
      printf("Error in ftable: out of bounds\n");
      it = 99;
    }
    sf = ftable[type][it];
    scftmp = sf * exp(-sthol*sthol*atom.BValue());
    scftmp *= atom.occ();
    ah = refl[ir].ind[0];
    ak = refl[ir].ind[1];
    al = refl[ir].ind[2];
    for (in = 0; in < nsym; in++) {
      /*
         trans = ah * mh->symops[X][3][in] +
         ak * mh->symops[Y][3][in] +
         al * mh->symops[Z][3][in];

         xpart = ah * mh->symops[X][X][in] +
         ak * mh->symops[Y][X][in] +
         al * mh->symops[Z][X][in];

         ypart = ah * mh->symops[X][Y][in] +
         ak * mh->symops[Y][Y][in] +
         al * mh->symops[Z][Y][in];

         zpart = ah * mh->symops[X][Z][in] +
         ak * mh->symops[Y][Z][in] +
         al * mh->symops[Z][Z][in];
       */
      index = nsym*ir +in;
      phase = twopi*(xpart[index]*x + ypart[index]*y + zpart[index]*z + trans[index]);
      cp = cos(phase)*scftmp;
      sp = sin(phase)*scftmp;
      refl[ir].acalc += (float)cp;
      refl[ir].bcalc += (float)sp;
    }
  }
  return (nr);
}

/* calculate structure factors for a list of residues
 * if init = 0 then sums are zeroed before beginning
 * else they are summed
 */
int sfcalc(RESIDUE*  res, CREFL refl[], int nrefl, CMapHeaderBase* mh, int init) {
  int i, n = 0, nr = 0;
  float natoms = 0.0;
  char buf[200];
  RESIDUE* start;

  start = res;
  if (init == 0) {
    Logger::log("Clearing any old structure factors...");
    for (i = 0; i < nrefl; i++) {
      refl[i].acalc = 0.0;
      refl[i].bcalc = 0.0;
    }
  }

  while (res != NULL) {
    if (!((strcmp(res->type().c_str(), "BND") == 0 && res->name().size() > 0 && res->name()[0] == '#')) ) {
      natoms += res->atomCount();
    }
    res = res->next();
  }
  res = start;
  while (res != NULL) {
    if (!(strcmp(res->type().c_str(), "BND") == 0 && res->name().size() > 0 && res->name()[0] == '#') ) {
      res = res->next();
      continue;
    }
    for (i = 0; i < res->atomCount(); i++) {
      nr = sfatom(*res->atom(i), refl, nrefl, mh);
      if (nr) {
        n++;
      }
    }
    //sprintf(buf,"SFCalc: %0.1f percent done",100.0*(float)n/natoms);
    //leftfooter(buf);
    res = res->next();
  }
  sprintf(buf, "Calculated str factors for %d atoms", n);
  Logger::log(buf);
  return (nr);
}

int sfcalcatom(MIAtom* atoms[], int natoms, CREFL refl[], int nrefl, CMapHeaderBase* mh, int init) {
  int i, n = 0, nr = 0;
  char buf[200];

  if (init == 0) {
    Logger::log("Clearing any old structure factors...");
    for (i = 0; i < nrefl; i++) {
      refl[i].acalc = 0.0;
      refl[i].bcalc = 0.0;
    }
  }
  for (i = 0; i < natoms; i++) {
    nr = sfatom(*atoms[i], refl, nrefl, mh);
    if (nr) {
      n++;
    }
    //sprintf(buf,"SFCalc: %0.1f percent done",(float)n*100.0/(float)natoms);
    //leftfooter(buf);
  }
  sprintf(buf, "Calculated str factors for %d atoms", n);
  Logger::log(buf);
  return (nr);
}


float ComputeScale(std::vector<CREFL>& refl, CMapHeaderBase* mh) {
  float scale = 0.0;
  double fcsum = 0.0, fosum = 0.0;
  for (unsigned int i = 0; i < refl.size(); i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    fosum += refl[i].fo;
    refl[i].fc = (float)sqrt(refl[i].acalc*refl[i].acalc +
                   refl[i].bcalc*refl[i].bcalc);
    fcsum += refl[i].fc;
  }
  scale = (float)(fosum / fcsum);
  return scale;
}

float ComputeScale2(CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  int i;
  float scale = 0.0;
  double fcsum = 0.0, fosum = 0.0;
  float s2, B = 15.0, K = 1.0, sc;
  if (mh->use_bulksolvent) {
    mh->use_bulksolvent = CalcBulkSolvent(refl, nrefl, mh);
    B = mh->Bsolvent;
    K = mh->Ksolvent;
  }
  if (mh->use_bulksolvent == 0) {
    return 1.0;
  }
  for (i = 0; i < nrefl; i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    fosum += refl[i].fo;
    if (mh->use_bulksolvent) {
      s2 = refl[i].sthol * refl[i].sthol;
      sc = (float)(1.0-(K*exp(-B*s2)));
    } else {sc = 1.0;}
    fcsum += refl[i].fc *sc;
  }
  scale = (float)(fosum / fcsum);
  printf("Scale = %0.5f , (%0.2f / %0.2f)\n", scale, fosum, fcsum);
  return (scale);
}

int ApplyScale(std::vector<CREFL>& refl, float scale, CMapHeaderBase* mh) {
  int n = 0;
  for (unsigned int i = 0; i < refl.size(); i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    refl[i].fo /= scale;
    refl[i].sigma /= scale;
    n++;
  }
  return (n);
}

int
ApplyScale2(CREFL refl[], int nrefl, float scale, CMapHeaderBase* mh) {
  int i, n = 0;
  float s2, B, K, sc;
  for (i = 0; i < nrefl; i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    refl[i].fo /= scale;
    refl[i].sigma /= scale;
    n++;
  }
  if (mh->use_bulksolvent) {
    B = mh->Bsolvent;
    K = mh->Ksolvent;
    for (i = 0; i < nrefl; i++) {
      if (refl[i].sthol > 0.5/mh->resmin
          || refl[i].sthol < 0.5/mh->resmax) {
        continue;
      }
      s2 = refl[i].sthol * refl[i].sthol;
      sc = (float)(1.0-(K*exp(-B*s2)));
      refl[i].fc *= sc;
      refl[i].acalc *= sc;
      refl[i].bcalc *= sc;
      refl[i].awhole *= sc;
      refl[i].bwhole *= sc;
      n++;
    }
  }
  return (n);
}

int CalcBulkSolvent(CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  float bestR, r;
  float* B, * K;
  float sumfo = 0.0, sumdif = 0.0;
  int i = 0;
  B = &(mh->Bsolvent);
  K = &(mh->Ksolvent);
  if (*B <= 0.0) {
    *B = 150.0F;
  }
  if (*K <= 0.0) {
    *K = 0.90F;
  }
  for (i = 0; i < nrefl; i++) {
    if (refl[i].sthol > 0.5F/mh->resmin
        || refl[i].sthol < 0.5F/mh->resmax) {
      continue;
    }
    sumfo += refl[i].fo;
    sumdif += (float)fabs(refl[i].fo-refl[i].fc);
  }
  if (sumfo > 0.0) {
    printf("---Calculate Bulk Solvent Correction---\nStart:   R=%0.3f\n", sumdif/sumfo);
  }
  bestR = sumdif/sumfo;
  for (i = 0; i < 10; i++) {
    r = EstimateBulkSolvent(refl, nrefl, B, K, i, mh);
    printf("Cycle %d: R=%0.3f Bsolv=%0.1f Ksolv=%0.3f\n", i+1, r, *B, *K);
    if (r < bestR-0.003) {
      bestR = r;
    } else {
      goto end;
    }
  }
end:
  if (i == 0) {
    printf("R got worse or no better - no correction will be applied\n");
    mh->Ksolvent = 0.0;
    return 0;
  } else {
    return 1;
  }

}

float EstimateBulkSolvent(CREFL refl[], int nrefl, float* B, float* K, int ntimes, CMapHeaderBase* mh) {
  int i, iB, iK;
  float trialB, trialK, bestB, bestK, bestR;
  float scale;
  float fcsum = 0.0, fosum = 0.0;
  float s2;
  bestR = 100.0;
  bestB = *B;
  bestK = *K;
  for (iB = -4; iB <= 5; iB++) {
    trialB = *B + (float)iB*(*B/(5.0F+(float)ntimes));
    for (iK = -3; iK <= 3; iK++) {
      trialK = *K + (float)iK*(*K/(9.0F+(float)ntimes));
      fosum = 0.0;
      fcsum = 0.0;
      for (i = 0; i < nrefl; i++) {
        if (refl[i].sthol > 0.5F/mh->resmin
            || refl[i].sthol < 0.5F/mh->resmax) {
          continue;
        }
        s2 = refl[i].sthol * refl[i].sthol;
        fcsum += refl[i].fc * (float)(1.0-(trialK*exp(-trialB*s2)));
        fosum += refl[i].fo;
      }
      scale = fosum / fcsum;
      fcsum = 0.0;
      fosum = 0.0;
      for (i = 0; i < nrefl; i++) {
        if (refl[i].sthol > 0.5F/mh->resmin
            || refl[i].sthol < 0.5F/mh->resmax) {
          continue;
        }
        s2 = refl[i].sthol * refl[i].sthol;
        fcsum += (float)fabs(refl[i].fo/scale - (refl[i].fc * (float)(1.0-(trialK*exp(-trialB*s2)))));
        fosum += refl[i].fo/scale;
      }
      if (fcsum/fosum < bestR) {
        bestK = trialK;
        bestB = trialB;
        bestR = fcsum/fosum;
      }
    }
  }
  *B = bestB;
  *K = bestK;
  return bestR;
}

int SubtractPartial(CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  int i;
  for (i = 0; i < nrefl; i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    refl[i].apart = refl[i].acalc;
    refl[i].bpart = refl[i].bcalc;
    refl[i].acalc = refl[i].awhole - refl[i].apart;
    refl[i].bcalc = refl[i].bwhole - refl[i].bpart;
  }
  return 1;
}

int GetStatic(CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  int i;
  for (i = 0; i < nrefl; i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    refl[i].astatic = refl[i].awhole - refl[i].apart;
    refl[i].bstatic = refl[i].bwhole - refl[i].bpart;
  }
  return 1;
}

int AddStaticPartial(CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  int i;
  for (i = 0; i < nrefl; i++) {
    if (refl[i].sthol > 0.5/mh->resmin
        || refl[i].sthol < 0.5/mh->resmax) {
      continue;
    }
    refl[i].acalc = refl[i].astatic + refl[i].apart;
    refl[i].bcalc = refl[i].bstatic + refl[i].bpart;
  }
  return 1;
}

float RePhase(std::vector<CREFL>& refl, CMapHeaderBase* mh) {
  int n = 0;
  float phi, sumphi = 0.0, dphi;
  float fosum = 0.0, difsum = 0.0;
  float rfactor;
  double degtor;
  char buf[200];

  degtor = 180.0/acos(-1.0);
  for (unsigned int i = 0; i < refl.size(); i++) {
    if (refl[i].sthol > 0.5F/mh->resmin
        || refl[i].sthol < 0.5F/mh->resmax) {
      continue;
    }
    n++;
    phi = (float)(atan2(refl[i].bcalc, refl[i].acalc) *degtor);
    refl[i].fc = (float)sqrt(refl[i].acalc*refl[i].acalc +
                   refl[i].bcalc*refl[i].bcalc);
    dphi = (float)fabs(phi - refl[i].phi);
    if (dphi > 180.0F) {
      dphi = 360.0F-dphi;
    }
    sumphi += dphi;
    refl[i].phi = phi;
    difsum += (float)fabs(refl[i].fc - refl[i].fo);
    fosum += refl[i].fo;
    /*
       if(i<10)
       printf("%d %d %d %0.3f %0.3f %0.3f a=%0.3f b=%0.3f st=%0.3f\n",refl[i].ind[0],refl[i].ind[1],refl[i].ind[2],refl[i].fo,refl[i].fc,refl[i].phi,refl[i].acalc, refl[i].bcalc, refl[i].sthol);
     */
  }
  /*
     sprintf(buf,"sumphi = %0.3f  nrefl = %d",sumphi,n);
     Logger::log(buf);
   */
  rfactor =  difsum/fosum;
  sprintf(buf, "R-factor = %0.3f", rfactor);
  Logger::log(buf);
  return ( sumphi / (float)n);
}

int randomize(RESIDUE* res, float maxdev) {
  int i;
  int n = 0;
  float rsum = 0.0;
  float dsum = 0.0;
  float dx, dy, dz;
  float r;
  float d = (float)RAND_MAX;
  char buf[100];

  while (res != NULL) {
    for (i = 0; i < res->atomCount(); i++) {
      r = (float)rand();
      r = r/d * maxdev;
      rsum += r*r;
      res->atom(i)->translate(r, 0.0f, 0.0f);
      dx = r;
      n++;
      r = (float)rand();
      r = r/d * maxdev;
      rsum += r*r;
      res->atom(i)->translate(0.0f, r, 0.0f);
      dy = r;
      n++;
      r = (float)rand();
      r = r/d * maxdev;
      rsum += r*r;
      res->atom(i)->translate(0.0f, 0.0f, r);
      dz = r;
      n++;
      dsum += (float)sqrt(dx*dx + dy*dy + dz*dz);
    }
    res = res->next();
  }
  r = (float)sqrt(rsum/(double)n);
  sprintf(buf, "Randomize: RMS fluctuation added to coords = %f", r);
  Logger::log(buf);
  r = (float)sqrt(dsum/3.0/(double)n);
  sprintf(buf, "    resulting in an RMS distance of = %f", r);
  Logger::log(buf);
  return (n);
}

int arrow(LINE* vu, float dx, float dy, float dz, float ox, float oy, float oz, int color) {
  /* draw an arrow along vector dx,dy,dz
   * starting at ox,oy,oz in LINEs */
  float d, theta;
  float a[3], b[3], c[3];
  float mat[4][3];
  float rtodeg = 0.f;
  static int init = 0;
  int i;
  float x, y, z;
  //extern float vectorangle();
  if (init == 0) {
    rtodeg = (float)(180.0/acos(-1.0));
    init = 1;
  }
  d = (float)sqrt(dx*dx+dy*dy+dz*dz);
  /* draw arrow along x */
  /* arrow shaft */
  vu[0].x1 = 0.0;
  vu[0].y1 = 0.0;
  vu[0].z1 = 0.0;
  vu[0].x2 = d;
  vu[0].y2 = 0.0;
  vu[0].z2 = 0.0;
  /* arrow head 0.1 of shaft length */
  vu[1].x1 = d;
  vu[1].y1 = 0.0;
  vu[1].z1 = 0.0;
  vu[1].x2 = d-d*0.1F;
  vu[1].y2 = -d*0.08F;
  vu[1].z2 = 0.0;
  vu[2].x1 = d;
  vu[2].y1 = 0.0;
  vu[2].z1 = 0.0;
  vu[2].x2 = d-d*0.1F;
  vu[2].y2 = d*0.08F;
  vu[2].z2 = 0.0;
  vu[3].x1 = d;
  vu[3].y1 = 0.0;
  vu[3].z1 = 0.0;
  vu[3].x2 = d-d*0.1F;
  vu[3].y2 = 0.0;
  vu[3].z2 = -d*0.08F;
  vu[4].x1 = d;
  vu[4].y1 = 0.0;
  vu[4].z1 = 0.0;
  vu[4].x2 = d-d*0.1F;
  vu[4].y2 = 0.0;
  vu[4].z2 = d*0.08F;
  /* build matrix for rotating from x to (dx,dy,dz) */
  /* find cross vector = axis about which to rotate */
  a[X] = dx;
  a[Y] = dy;
  a[Z] = dz;
  b[X] = d;
  b[Y] = 0.0;
  b[Z] = 0.0;
  cross(a, b, c);
  theta = vectorangle(b, a)*rtodeg;
  initrotate(0.0, 0.0, 0.0, c[0], c[1], c[2], -theta, mat);
  for (i = 0; i < 5; i++) {
    x = vu[i].x1;
    y = vu[i].y1;
    z = vu[i].z1;
    xl_rotate(x, y, z, &(vu[i].x1), &(vu[i].y1), &(vu[i].z1), mat);
    x = vu[i].x2;
    y = vu[i].y2;
    z = vu[i].z2;
    xl_rotate(x, y, z, &(vu[i].x2), &(vu[i].y2), &(vu[i].z2), mat);
    vu[i].x1 += ox;
    vu[i].x2 += ox;
    vu[i].y1 += oy;
    vu[i].y2 += oy;
    vu[i].z1 += oz;
    vu[i].z2 += oz;
    vu[i].color = (short)color;
    vu[i].color2 = (short)color;
  }
  return (5);
}

int read_scattering_factors(char*** name, int** lname, char*** rname, double** a1, double** a2, double** a3, double** a4,
                            double** b1, double** b2, double** b3, double** b4, double** c,
                            double** zeff, double** frano, double** fiano,
                            const std::string& crystal_data_dir,
                            const std::string& molimage_home_dir) {
  /* read in ScatteringFactors file */
  FILE* fp = fopen("ScatteringFactors", "r");
  int i, n = 0, nread;
  char file[1024];
  char buf[1024];
  char string[200], rstring[200];
  float aa1, aa2, aa3, aa4, bb1, bb2, bb3, bb4, cc, zzeff, ffrano, ffiano;
  if (!fp) {
    sprintf(file, "%s/ScatteringFactors", crystal_data_dir.c_str());
    fp = fopen(file, "r");
  } else {strcpy(file, "./ScatteringFactors");}
  if (!fp) {
    sprintf(file, "%s/data/ScatteringFactors", molimage_home_dir.c_str());
    fp = fopen(file, "r");
  }
  if (!fp) {
    Logger::log("Can't open ScatteringFactors file - check installation data directory.");

    return 0;
  }
  sprintf(buf, "Using scattering factors from %s", file);
  Logger::log(buf);
  /* first get number of scattering factors */
  while (fgets(file, sizeof file, fp) != NULL) {
    if (file[0] == '#' || file[0] == '!') {
      continue;
    }
    nread = sscanf(file, "%s%s%f%f%f%f%f%f%f%f%f%f%f%f", string, rstring,
              &aa1, &bb1, &aa2, &bb2, &aa3, &bb3, &aa4, &bb4,
              &cc, &zzeff, &ffrano, &ffiano);
    if (nread != 14) {
      Logger::log("Warning: can't interpret ScatteringFactors record:");
      Logger::log(file);
      continue;
    }
    n++;
  }
  *a1 = (double*)malloc(n*sizeof(double));
  *a2 = (double*)malloc(n*sizeof(double));
  *a3 = (double*)malloc(n*sizeof(double));
  *a4 = (double*)malloc(n*sizeof(double));
  *b1 = (double*)malloc(n*sizeof(double));
  *b2 = (double*)malloc(n*sizeof(double));
  *b3 = (double*)malloc(n*sizeof(double));
  *b4 = (double*)malloc(n*sizeof(double));
  *c = (double*)malloc(n*sizeof(double));
  *zeff = (double*)malloc(n*sizeof(double));
  *frano = (double*)malloc(n*sizeof(double));
  *fiano = (double*)malloc(n*sizeof(double));
  *name = (char**)malloc(sizeof(char*)*n);
  for (i = 0; i < n; i++) {
    (*name)[i] = (char*)malloc(2*sizeof(char)*MAXHNAME);
  }
  *rname = (char**)malloc(sizeof(char*)*n);
  for (i = 0; i < n; i++) {
    (*rname)[i] = (char*)malloc(2*sizeof(char)*MAXHNAME);
  }
  *lname = (int*)malloc(n*sizeof(int));
  if (*lname == NULL) {
    Logger::log("Out of memory in read_scattering_factors");
    return 0;
  }
  rewind(fp);
  n = 0;
  while (fgets(file, sizeof file, fp) != NULL) {
    if (file[0] == '#' || file[0] == '!') {
      continue;
    }
    nread = sscanf(file, "%s%s%f%f%f%f%f%f%f%f%f%f%f%f", string, rstring,
              &aa1, &bb1, &aa2, &bb2, &aa3, &bb3, &aa4, &bb4,
              &cc, &zzeff, &ffrano, &ffiano);
    if (nread != 14) {
      Logger::log("Warning: can't interpret ScatteringFactors record:");
      Logger::log(file);
      continue;
    }
    strncpy((*name)[n], string, 2*MAXHNAME-1);
    strncpy((*rname)[n], rstring, 2*MAXHNAME-1);
    (*lname)[n] = strlen((*name)[n]);
    (*a1)[n] = aa1;
    (*a2)[n] = aa2;
    (*a3)[n] = aa3;
    (*a4)[n] = aa4;
    (*b1)[n] = bb1;
    (*b2)[n] = bb2;
    (*b3)[n] = bb3;
    (*b4)[n] = bb4;
    (*c)[n] = cc;
    (*zeff)[n] = zzeff;
    (*frano)[n] = ffrano;
    (*fiano)[n] = ffiano;
    n++;
  }

  return n;
}

int default_scattering_factors(char*** name, int** lname, char*** rname, double** a1, double** a2, double** a3, double** a4,
                               double** b1, double** b2, double** b3, double** b4, double** c,
                               double** zeff, double** frano, double** fiano) {
  /* load default Scattering Factors */
  int i, n = 0;
  char buf[1024];
  n = sizeof(da1)/sizeof(double);
  sprintf(buf, "Initializing %d default scattering factors", n);
  Logger::log(buf);
  *a1 = (double*)malloc(n*sizeof(double));
  *a2 = (double*)malloc(n*sizeof(double));
  *a3 = (double*)malloc(n*sizeof(double));
  *a4 = (double*)malloc(n*sizeof(double));
  *b1 = (double*)malloc(n*sizeof(double));
  *b2 = (double*)malloc(n*sizeof(double));
  *b3 = (double*)malloc(n*sizeof(double));
  *b4 = (double*)malloc(n*sizeof(double));
  *c = (double*)malloc(n*sizeof(double));
  *zeff = (double*)malloc(n*sizeof(double));
  *frano = (double*)malloc(n*sizeof(double));
  *fiano = (double*)malloc(n*sizeof(double));
  *name = (char**)malloc(sizeof(char*)*n);
  for (i = 0; i < n; i++) {
    (*name)[i] = (char*)malloc(2*sizeof(char)*MAXHNAME);
  }
  *rname = (char**)malloc(sizeof(char*)*n);
  for (i = 0; i < n; i++) {
    (*rname)[i] = (char*)malloc(2*sizeof(char)*MAXHNAME);
  }
  *lname = (int*)malloc(n*sizeof(int));
  if (*lname == NULL) {
    Logger::log("Out of memory in default_scattering_factors");
    return 0;
  }
  for (i = 0; i < n; i++) {
    strncpy((*name)[i], datypes[i], 2*MAXHNAME-1);
    strncpy((*rname)[i], "*", 2*MAXHNAME-1);
    (*lname)[i] = strlen((*name)[i]);
    (*a1)[i] = da1[i];
    (*a2)[i] = da2[i];
    (*a3)[i] = da3[i];
    (*a4)[i] = da4[i];
    (*b1)[i] = db1[i];
    (*b2)[i] = db2[i];
    (*b3)[i] = db3[i];
    (*b4)[i] = db4[i];
    (*c)[i] = dco[i];
    (*zeff)[i] = dzeff[i];
    (*frano)[i] = dfrano[i];
    (*fiano)[i] = dfiano[i];
  }
  /*
     {
     FILE *fp=fopen("ScatteringFactors","w");
     Logger::log("Writing ScatteringFactors file");
     for(i=0;i<n;i++)
     fprintf(fp,"%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t\n",(*name)[i],(*rname)[i],
        da1[i],db1[i],
        da2[i],db2[i],
        da3[i],db3[i],
        da4[i],db4[i],
        dco[i],dzeff[i],
        dfrano[i],dfiano[i]);
     }
   */
  return n;
}

void
free_scattering_factors(char** name, int* lname, char** rname, double* a1, double* a2, double* a3,
                        double* a4, double* b1, double* b2, double* b3, double* b4, double* c,
                        double* zeff, double* frano, double* fiano) {
  int i;
  free(a1);
  free(a2);
  free(a3);
  free(a4);
  free(b1);
  free(b2);
  free(b3);
  free(b4);
  free(c);
  free(zeff);
  free(frano);
  free(fiano);

  for (i = 0; i < MAXFTABLE; i++) {
    free(name[i]);
  }
  free(name);

  for (i = 0; i < MAXFTABLE; i++) {
    free(rname[i]);
  }
  free(rname);

  free(lname);
}

int
ScattIndex(const char* aname, const char* restype) {
  int i, n = -1;

  for (i = 0; i < MAXFTABLE; i++) {
    if (atypes[i][0] == '*' || !strncmp(atypes[i], aname, atypelen[i])) {
      if (rtypes[i][0] == '*' || restype[0] == '*') {
        n = i;
      } else if (!strcmp(restype, rtypes[i])) {
        n = i;
      }
    }
  }
  return n;
}

#define MAXBINS 50
#ifdef MOVEDTOFFTSUBSC

/* compute sigmaA by the method of Read(1986) Acta Cryst A42, 140-149
 * as done in the CCP4 program sigmaa.f
 */
/*  computes m and D from a reflection list
 * a sigmaA map is of the form
 * 2mFo - DFc, alphacalc
 * or a difference sigmaA map
 * mFo - Fc, alphacalc
 * mFo-DFc is stored in acalc field
 * and mFo-Fc is stored in bcalc field
 */
int
SigmaA(refl, nrefl, mhin)
CREFL *refl;
int nrefl;
CMapHeaderBase* mhin;
{
  extern float sim( /* X */);
  float* epsilon;
  char buf[200];
  unsigned char* centric;
  int i, j, k, ibin, nbin[MAXBINS];
  int ih, ik, il;
  /* increased to 3 x 3 x 96 - dem 7/31/97 */
  int iss[864];
  /* increased to 3 x 96 - dem 7/31/97 */
  int its[288];
  float eps, wt;
  int mult, mk, iflg, iflg2;
  float sigman[MAXBINS], sigmap[MAXBINS];
  float sum22[MAXBINS], sum20[MAXBINS], sum02[MAXBINS], sum40[MAXBINS], sum04[MAXBINS];
  float sum22T = 0, sum20T = 0, sum02T = 0, sum40T = 0, sum04T = 0;
  float sumwt[MAXBINS];
  float sigmaA[MAXBINS], signum, sigden;
  float eo, ec, m, Xarg;
  float corr, mtot = 0.0;
  int ntot = 0, numbins = 10, nfom = 0;
#ifdef XVIEW
  float ratio = (float)atof((char*)xv_get(fftpop->sigmaratio, PANEL_VALUE));
#else
  float ratio = 1.0F;
#endif

  /* try to catch nonsense values */
  if (ratio <= 0.0) {
    ratio = 1.0;
  }
  if (ratio > 3.0) {
    ratio = 1.0;
  }
  printf("Ratio of reference to working set is %f\n", ratio);


  if (nrefl < 100) {
    Logger::log("SigmaA: too few reflections");
    return 0;
  }

  /* determine number of bins such that each bin has at least
   * 500 reflections */
  numbins = std::min(MAXBINS, ROUND((float)nrefl/500.0));

  /* sort the reflections by sthol */
  qsort(refl, nrefl, sizeof(CREFL), scompare);

  epsilon = (float*)malloc(sizeof(float)*nrefl);
  centric = (unsigned char*)malloc(sizeof(unsigned char)*nrefl);
  if (!epsilon || !centric) {
    Logger::log("Out of memory in SigmaA: sorry");
    return 0;
  }
  /* first get the centric flag and the epsilon */
  /* set up the symmetry information */
  for (k = 0; k < mhin->nsym; k++) {
    for (j = 0; j < 3; j++) {
      its[j+k*3] =  ROUND(mhin->symops[j][3][k]*24.0);
      for (i = 0; i < 3; i++) {
        iss[i+3*(j+3*k)] = (int)mhin->symops[i][j][k];
      }
    }
  }
  for (i = 0; i < nrefl; i++) {
    ih = refl[i].ind[0];
    ik = refl[i].ind[1];
    il = refl[i].ind[2];
    stdrefl_(&ih, &ik, &il, &mult, &eps, &mk, &iflg,
      &iflg2, iss, its, &mhin->nsym);
    epsilon[i] = eps;
    centric[i] = mk-1;
    /* check to make sure that Fc column
     * does not contain figure-of-merits
     */
    if (refl[i].fc <= 1.0) {
      nfom++;
    }
  }
  if ((float)nfom/(float)nrefl > 0.80) {
    Logger::log("SigmaA: Fc's are figure-of-merits.");
    printf("SigmaA: Fc's are figure-of-merits.\n");
    for (i = 0; i < nrefl; i++) {
      refl[i].fom = refl[i].fc;
    }
    return 0;
  }

  /* zero the sums */
  for (i = 0; i < numbins; i++) {
    sum22[i] = 0;
    sum20[i] = 0;
    sum02[i] = 0;
    sum40[i] = 0;
    sum04[i] = 0;
    sumwt[i] = 0;
    sigman[i] = 0;
    sigmap[i] = 0;
    nbin[i] = 0;
  }

  /* now accumalate sums */
  for (i = 0; i < nrefl; i++) {
    ibin = (i*numbins)/(nrefl);
    wt = 2.0;
    if (centric[i]) {
      wt = 1.0;
    }
    sumwt[ibin] += wt;
    nbin[ibin]++;
    sigman[ibin] += refl[i].fo*refl[i].fo*wt/epsilon[i];
    sigmap[ibin] += refl[i].fc*refl[i].fc*wt/epsilon[i];
  }

  for (ibin = 0; ibin < numbins; ibin++) {
    if (sumwt[ibin]) {
      sigman[ibin] /= sumwt[ibin];
      sigmap[ibin] /= sumwt[ibin];
    }
  }

  for (i = 0; i < nrefl; i++) {
    ibin = (i*numbins)/(nrefl);
    eo = refl[i].fo/sqrt(sigman[ibin]*epsilon[i]);
    ec = refl[i].fc/sqrt(sigmap[ibin]*epsilon[i]);
    eo *= eo;
    ec *= ec;
    sum22[ibin] += eo*ec;
    sum20[ibin] += eo;
    sum02[ibin] += ec;
    sum40[ibin] += eo*eo;
    sum04[ibin] += ec*ec;
  }

  printf("SigmaA by resolution\n  Bin     Resolution   Number     SigmaA   Correlation\n");
  for (ibin = 0; ibin < numbins; ibin++) {
    if (nbin[ibin] > 0) {
      ntot += nbin[ibin];
      sum22T += sum22[ibin];
      sum20T += sum20[ibin];
      sum02T += sum02[ibin];
      sum40T += sum40[ibin];
      sum04T += sum04[ibin];
      signum = nbin[ibin]*sum22[ibin] - sum20[ibin]*sum02[ibin];
      if (signum > 0.0) {
        sigden = sqrt((nbin[ibin]*sum40[ibin]-sum20[ibin]*sum20[ibin])*
                   (nbin[ibin]*sum04[ibin]-sum02[ibin]-sum02[ibin]));
        sigmaA[ibin] = sqrt(signum/sigden)*ratio;
      } else {
        printf("Warning: correlation between E**2's is non-positive for bin %d - sigmaA set to 0.05\n", ibin);
        sigmaA[ibin] = 0.05;
      }
      corr = (nbin[ibin]*sum22[ibin] - sum20[ibin]*sum02[ibin])/
             sqrt((nbin[ibin]*sum40[ibin]-sum20[ibin]*sum20[ibin])*(nbin[ibin]*sum04[ibin]-sum02[ibin]*sum02[ibin]));
      printf("%5d  %6.2f-%6.2f   %6d    %6.5f   %6.5f\n", ibin+1,
        (float)(0.5/refl[ibin*nrefl/numbins].sthol),
        (float)(0.5/refl[std::min((ibin+1)*nrefl/numbins-1, nrefl-1)].sthol),
        nbin[ibin], sigmaA[ibin], corr);
    }
  }

  corr = (ntot*sum22T - sum20T*sum02T)/
         sqrt((ntot*sum40T-sum20T*sum20T)*(ntot*sum04T-sum02T*sum02T));
  sprintf(buf, "SigmaA: Overall correlation is %f for %d refl\n", corr, ntot);
  printf("%s", buf);
  Logger::log(buf);

  /* modify Fo and Fc according to sigmaA */
  for (i = 0; i < nrefl; i++) {
    ibin = i*numbins/nrefl;
    eo = refl[i].fo/sqrt(sigman[ibin]*epsilon[i]);
    ec = refl[i].fc/sqrt(sigman[ibin]*epsilon[i]);
    Xarg = sigmaA[ibin]*2.0*eo*ec/(1.0-sigmaA[ibin]*sigmaA[ibin]);
    if (!centric[i]) {
      m = sim(Xarg);
      refl[i].acalc = (2.0*m*eo-sigmaA[ibin]*ec)*sqrt(sigman[ibin]*epsilon[i]);
      refl[i].bcalc = m*refl[i].fo-refl[i].fc;
    } else {
      m = tanh(Xarg);
      refl[i].acalc = m*refl[i].fo;
      refl[i].bcalc = refl[i].fo-refl[i].fc;
    }
    refl[i].fom = m;
    mtot += m;
    /*
       if(i%100==0)printf("%3d %3d %3d %5.2f %5.2f %d %f\n",refl[i].ind[0],refl[i].ind[1],refl[i].ind[2],refl[i].fo, refl[i].fc, (int)centric[i], epsilon[i]);
     */
  }
  sprintf(buf, "SigmaA: Average figure of merit is %0.3f\n", mtot/(float)ntot);
  Logger::log(buf);
  printf("%s", buf);

  free(epsilon);
  free(centric);
  return nrefl;
}

int scompare(void* in, void* jn) {
  /* use with qsort to sort reflections according to sthol */
  CREFL* i = (CREFL*)in;
  CREFL* j = (CREFL*)jn;
  int k = 1;
  if (i->sthol == j->sthol) {
    return (0);
  }
  if (i->sthol < j->sthol) {
    return (-1);
  }
  return (k);
}

float
sim(x)
float x;
{
  /*
     Calculate Sim & Srinivasan non-centric figure of merit as
     I1(X)/I0(X), where I1 and I0 are the modified 1st and zero
     order bessel functions.
     References: Sim, G. A. (1960) Acta Cryst. 13, 511-512;
               Srinivasan, R. (1966) Acta Cryst. 20, 143-144;
               Abramowitz & Stegun, Handbook of Mathematical Functions, 378.
   */
  float t;
  float sim = 0.0;
  if (fabs(x) > 0.0001) {
    t = fabs(x)/3.75;
    if (t > 1.0) {
      t = 1.0/t;
      sim = ((((((((0.01787654-t*0.00420059)*t+ (-0.02895312))*t+
                  0.02282967)*t+ (-0.01031555))*t+0.00163801)*t+
               (-0.00362018))*t+ (-0.03988024))*t+0.39894228)*
            SIGN(x)/ ((((((((-0.01647633+t*0.00392377)*t+
                            0.02635537)*t+ (-0.02057706))*t+0.00916281)*t+
                         (-0.00157565))*t+0.00225319)*t+0.01328592)*t+0.39894228);
    } else {
      t = t*t;
      sim = ((((((t*0.00032411+0.00301532)*t+0.02658733)*t+
                0.15084934)*t+0.51498869)*t+0.87890594)*t+0.5)*x/
            ((((((t*0.0045813+0.0360768)*t+0.2659732)*t+
                1.2067492)*t+3.0899424)*t+3.5156229)*t+1.0);
    }
  }
  return sim;
}
#endif /*MOVEDTOFFTSUBSC*/

int scaleaniso(CREFL refl[], int nrefl, CMapHeaderBase* mh) {
  float resbin[MAXBINS], rfactbin[MAXBINS];
  float deltabin[MAXBINS];
  int number = 0;
  float diff, scale;
  float sthol, width;
  float fp, fm, delta, wtdel;
  int i, j;
  float ssqmax, ssqmin;
  char buff[200];
  float ah, ak, al;
  float as, bs, cs, cosgs, cosbs, cosas, V;
  float cosa, cosb, cosg, sina, sinb, sing;
  float alpha, beta, gamma;
  float rr[7], sc[7];
  float am[22], v[7], dv[7], diag[7];
  float scbin[7][MAXBINS];
  float dstar, zh, zk, zl;
  long int nv = 6, nm, ising;
  float sumwr;
  float sumn1, sumn2, sumd1;
  int no;
  int loop;
  float or11, or21, or22, or31, or32, or33;
  int fail, ii, iid;
  int numberbins = 1;
  float degtor;
  //int sigma = 0;
  int ih, ik, il;
  degtor = (float)(acos(-1.0)/180.0);

  if (mh->maptype == (int)MIMapType::Fofom || mh->fc_is_fom != 0) {
    Logger::log("Cannot apply anisotropic scaling to Fo*f.o.m phase type");
    return 0;
  }

  for (i = 0; i < numberbins; i++) {
    rfactbin[i] = 0.0;
    deltabin[i] = 0.0;
  }
  nm = nv * ( nv + 1)/2;
  alpha = mh->alpha * degtor;
  beta = mh->beta * degtor;
  gamma = mh->gamma * degtor;
  cosa = (float)cos(alpha);
  cosb = (float)cos(beta);
  cosg = (float)cos(gamma);
  sina = (float)sin(alpha);
  sinb = (float)sin(beta);
  sing = (float)sin(gamma);
  cosas = (cosb*cosg - cosa)/sinb/sing ;
  cosbs = (cosa*cosg - cosb)/sina/sing ;
  cosgs = (cosa*cosb - cosg)/sina/sinb ;
  V = mh->a*mh->b*mh->c * (float)sqrt(1-cosa*cosa-cosb*cosb-cosg*cosg+2.0*cosa*cosg*cosb);
  as = mh->b*mh->c*sina/V;
  bs = mh->a*mh->c*sinb/V;
  cs = mh->a*mh->b*sing/V;
  or11 = sing* (float)sin(acos(cosbs));
  or21 = -cosg*(float)sin(acos(cosbs));
  or22 = (float)sin(acos(cosas));
  or31 = cosbs;
  or32 = cosas;
  or33 = 1.0;

  ssqmax = 0.5F/mh->resmin;
  ssqmax = ssqmax*ssqmax;
  ssqmin = 0.5F/mh->resmax;
  ssqmin = ssqmin*ssqmin;
  width = (ssqmax - ssqmin) / (float)numberbins;
  if (width <= 0.0) {
    sprintf(buff, "Bin width too small %f\n", width);
    Logger::log(buff);
    return (-1);
  }
  printf("---Anisotropic Scaling---------------\n");

  for (i = 0; i < numberbins; i++) {
    /* initialize matrices */
    sc[1] = 1.0;
    sc[2] = 1.0;
    sc[3] = 1.0;
    sc[4] = 0.0;
    sc[5] = 0.0;
    sc[6] = 0.0;
    sumn1 = 0.0;
    sumn2 = 0.0;
    sumd1 = 0.0;
    sumwr = 0.0;
    for (loop = 1; loop <= nm; loop++) {
      am[loop] = 0.0;
    }
    for (loop = 1; loop <= nv; loop++) {
      v[loop] = 0.0;
    }
    /* first pass to find scale factor */
    number = 0;
    no = 0;
    for (j = 0; j < nrefl; j++) {
      ah = (float)refl[j].ind[0];
      ak = (float)refl[j].ind[1];
      al = (float)refl[j].ind[2];
      sthol = refl[j].sthol;
      dstar = sthol/ 0.5F;
      if (sthol*sthol > (ssqmin+i*width)
          && sthol*sthol < (ssqmin+(i+1)*width) ) {
        if (refl[i].fo > 0.0) {
          fp = refl[j].fo;
          fm = refl[j].fc;
          delta = (float)fabs(fp -fm);
          if (delta < (fp+fm)/2.0) {
            zh = ah * as;
            zk = ak * bs;
            zl = al * cs;
            zl = (or31*zh+or32*zk+or33*zl)/dstar;
            zk = (or21*zh+or22*zk)/dstar;
            zh = or11*zh/dstar;
            /*            zh,zk,zl now directional cosines of reciprocal lattice */
            rr[1] = zh*zh;
            rr[2] = zk*zk;
            rr[3] = zl*zl;
            rr[4] = zh*zk*2;
            rr[5] = zh*zl*2;
            rr[6] = zk*zl*2;
            scale = 0.0;
            for (loop = 1; loop <= 6; loop++) {
              scale = scale + sc[loop]*rr[loop];
            }
            delta = fp - scale*fm;
            wtdel = delta;
            for (loop = 1; loop <= 6; loop++) {
              dv[loop] = fm*rr[loop];
            }
            matset(&am[1], &dv[1], &nm, &nv);
            for (loop = 1; loop <= nv; loop++) {
              v[loop] = v[loop]+wtdel*dv[loop];
            }
          }
        }
      }
    }

    ising = 0;
    fail = 0;
    ii = nm;
    iid = nv;
    for (loop = 1; loop <= nv; loop++) {
      diag[loop] = am[loop];
      if (diag[loop] == 0.0) {
        sprintf(buff, "ANISOSCALE: Diagonal element zero for variable %d\n", loop);
        Logger::log(buff);
        ising = 1;
      }
      ii = ii - iid;
      iid = iid -1;
    }
    if (ising == 1) {
      sprintf(buff, "ANISOSCALE: Matrix is singular!! Try decreasing number of bins\n");
      Logger::log(buff);
      return -1;
    }
    smi(am+1, &nm, &nv, &ising);
    if (ising != 0) {
      sprintf(buff, "ANISOSCALE: Matrix is singular?! Try decreasing number of bins\n");
      Logger::log(buff);
      return -1;
    }
    solv(&am[1], &v[1], &dv[1], &diag[1], &nm, &nv, &sumwr);
    for (loop = 1; loop <= 6; loop++) {
      sc[loop] = sc[loop]+ dv[loop];
    }
    printf("Shell %d: Scale Factors:", i+1);
    for (loop = 1; loop <= 6; loop++) {
      scbin[loop][i] = sc[loop];
      printf(" %0.3f", sc[loop]);
    }

    printf("\n    Resltn %0.3f-%0.3f", (float)(0.5/(sqrt(ssqmin+(i+1)*width))), (float)(0.5/(sqrt(ssqmin+i*width))));

    /*  second pass to apply scale factor */
    for (j = 0; j < nrefl; j++) {
      ah = (float)refl[j].ind[0];
      ak = (float)refl[j].ind[1];
      al = (float)refl[j].ind[2];
      sthol = refl[j].sthol;
      dstar = sthol/ 0.5F;
      if (sthol*sthol > (ssqmin+i*width)
          && sthol*sthol < (ssqmin+(i+1)*width) ) {
        if (refl[j].fo > 0.0) {
          ih = (int)ah;
          ik = (int)ak;
          il = (int)al;
          fp = refl[j].fo;
          fm = refl[j].fc;
          delta = (float)fabs(fp -fm);
          sumd1 += fp;
          sumn1 += delta;
          zh = ah * as;
          zk = ak * bs;
          zl = al * cs;
          zl = (or31*zh+or32*zk+or33*zl)/dstar;
          zk = (or21*zh+or22*zk)/dstar;
          zh = or11*zh/dstar;
          /*           zh,zk,zl now directional cosines of reciprocal lattice */
          rr[1] = zh*zh;
          rr[2] = zk*zk;
          rr[3] = zl*zl;
          rr[4] = zh*zk*2;
          rr[5] = zh*zl*2;
          rr[6] = zk*zl*2;
          scale = 0.0;
          for (loop = 1; loop <= 6; loop++) {
            scale = scale + sc[loop]*rr[loop];
          }
          refl[j].fc *= scale;
          refl[j].acalc *= scale;
          refl[j].bcalc *= scale;
          refl[j].awhole *= scale;
          refl[j].bwhole *= scale;
          diff = (float)fabs(fp - refl[j].fc);
          no++;
          sumn2 = sumn2 + diff;
        }
      }
    }
    if (sumd1 > 0.0) {
      rfactbin[i] = sumn2/sumd1;
      deltabin[i] = sumn2/(float)no;
      printf("  Rstart=%0.3f Rscaled=%0.3f  n=%d rfls\n", sumn1/sumd1, rfactbin[i], no);
    } else {
      rfactbin[i] = 0.0;
      deltabin[i] = 0.0;
    }
    resbin[i] = 0.5F/(float)(sqrt(ssqmin + i*width));
    resbin[i+1] = 0.5F/(float)(sqrt(ssqmin + (i+1)*width));
  }
  return (1);
}

unsigned int MIMapInitializeScatteringFactorTables(
  const std::string& crystal_data_dir,
  const std::string& molimage_home_dir) {
  MAXFTABLE = read_scattering_factors(&atypes, &atypelen, &rtypes, &a1, &a2, &a3, &a4, &b1, &b2, &b3, &b4, &co, &zeff, &frano, &fiano,
                crystal_data_dir,
                molimage_home_dir);
  if (MAXFTABLE == 0) {
    MAXFTABLE = default_scattering_factors(&atypes, &atypelen, &rtypes, &a1, &a2, &a3, &a4, &b1, &b2, &b3, &b4, &co, &zeff, &frano, &fiano);
    if (MAXFTABLE == 0) {
      Logger::log("Can't find scattering factors - Structure factor calculations not available");
    }
  }
  return MAXFTABLE;
}

unsigned int MIMapFreeScatteringFactorTables() {
  free_scattering_factors(atypes, atypelen, rtypes, a1, a2, a3, a4, b1, b2, b3, b4, co, zeff, frano, fiano);
  return 0;
}

