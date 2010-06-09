#ifndef SFCALCH
#define SFCALCH

#include <vector>

#include "CREFL.h"
#include "CMapHeaderBase.h"

namespace chemlib
{
    class Residue;
    class MIAtom;
}

#define MAXTYPES 20
/* Find the maximum value for rand(). Use 2^31-1 if not in stdlib.h */
#include <cstdlib>
#ifndef RAND_MAX
#define RAND_MAX 2147483647.0
#endif

//private
float ComputeScale(std::vector<CREFL> &refl, CMapHeaderBase *mh);
int ApplyScale(std::vector<CREFL> &refl, float scale, CMapHeaderBase *mh);
float ComputeScale2(CREFL refl[], int nrefl, CMapHeaderBase *mh);
float RePhase(std::vector<CREFL> &refl, CMapHeaderBase *mh);
int sfcalc(chemlib::Residue *  res, CREFL refl[], int nrefl, CMapHeaderBase *mh, int init);
int AtomDeriv();
int ScattIndex(const char *aname, const char *restype);
int CalcBulkSolvent(CREFL refl[], int nrefl, CMapHeaderBase *mh);
float EstimateBulkSolvent(CREFL refl[], int nrefl, float *B, float *K, int ntimes, CMapHeaderBase *mh);
int sfcalcatom(chemlib::MIAtom *atoms[], int natoms, CREFL refl[], int nrefl, CMapHeaderBase *mh, int init);
void sfinit();


#endif /* SFCALCH */
