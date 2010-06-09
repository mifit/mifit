#include <algorithm>

#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <chemlib/MIMoleculeBase.h>
#include <math/mathlib.h>

#include "EMapBase.h"
#include "MINATOM.h"
#include "PEAK.h"

#include "fssubs.h"  // private to library
#include "fft.h"     // private to library
#include "sfcalc.h"  // private to library
#include "rescalc.h" // private to library
#include "fssubs.h"  // private to library

typedef float real;
typedef struct
{
    real r, i;
} micomplex;


#ifdef _WIN32
#define _MVS
#define i386
#include <umtz/mmtzlib.h>
#undef _MVS
#define strncasecmp strnicmp
#else
#include <umtz/mmtzlib.h>
#endif

using namespace chemlib;
using namespace std;


// keep these all in sync!
const unsigned int MAP_TYPE_COUNT = 11;
static const char *maptypes[MAP_TYPE_COUNT] = {"Fo", "Fc", "2Fo-Fc", "Fo-Fc", "Fo*Fo", "Fo*fom",
                                               "3Fo-2Fc", "5Fo-3Fc", "2mFo-DFc(SigmaA)", "Fo-DFc(SigmaA)",
                                               "Direct FFT"};
static bool regular_maptype[MAP_TYPE_COUNT] = {true, true, true, false, true, true,
                                               true, true, true, false, true };

// end

bool IsRegularMapType(unsigned int maptype)
{
    if (maptype >= MAP_TYPE_COUNT)
        return true;
    return regular_maptype[maptype];
}

const char *StringForMapType(unsigned int maptype)
{
    if (maptype >= MAP_TYPE_COUNT)
        return "";
    return maptypes[maptype];
}

unsigned int MapTypeForString(const std::string &str)
{
    for (unsigned int i = 0; i<MAP_TYPE_COUNT; ++i)
    {
        if (strcasecmp(maptypes[i], str.c_str())==0)
            return i;
    }
    return UINT_MAX;
}



//defined in sfcalc.cpp
extern double Bscale;
extern char **rtypes;
extern char **atypes;
extern int *atypelen;
extern double *a1, *a2, *a3, *a4, *b1, *b2, *b3, *b4, *co;
extern double *fiano;
extern double *frano;
extern double *zeff;
extern int MAXFTABLE;

#define MAXLOOP 25
unsigned int scanrow(char *buf, char strings[MAXLOOP][81])
{
    return sscanf(buf, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                  strings[0], strings[1], strings[2], strings[3], strings[4], strings[5],
                  strings[6], strings[7], strings[8], strings[9], strings[10], strings[11],
                  strings[12], strings[13], strings[14], strings[15], strings[16], strings[17],
                  strings[18], strings[19]);
}

char *symmtrim(char *line)
{
    int l;
    while (isspace(line[0]) || line[0] == '\'')
    {
        line++;
    }
    l = strlen(line)-1;
    while (l > 0 && (line[l] == '\'' || isspace(line[l])))
    {
        line[l] = '\0';
        l--;
    }
    return line;
}

/**
 * An fget that also converts to lowercase and removes leading and trailing spaces.
 */
static char *fgetsl(char *buf, size_t size, FILE *fp)
{
    int i, ret;
    if (fgets(buf, size, fp) == NULL)
    {
        return NULL;
    }
    //strlower(buf);
    ret = strlen(buf);
    for (i = 0; i < ret-1; i++)
    {
        if (isupper(buf[i]))
        {
            buf[i] = tolower(buf[i]);
        }
    }
    i = ret-1;
    /* trim trailing edge white space */
    while (isspace(buf[i]))
    {
        buf[i] = '\0';
        i--;
        if (i < 0)
        {
            break;
        }
    }
    /* trim leading edge white space */
    i = 0;
    while (isspace(buf[0]))
    {
        memcpy(buf, buf+1, size-i-1);
        i++;
        if (i >= (int)(size-1))
        {
            break;
        }
    }
    return buf;
}

/* readin a binary integer with conditional unsigned char-swabbing */
static int map_in_int(FILE *fp, unsigned int swab)
{
    unsigned char in[4];
    fread(&(in[0]), 1, 4, fp);
    if (!swab)
    {
        return (int)((unsigned int)in[0]+((unsigned int)(in[1])<<8)+((unsigned int)(in[2])<<16)+((unsigned int)(in[3])<<24));
    }
    else
    {
        return (int)(((unsigned int)(in[0])>>24)+((unsigned int)(in[1])>>16)+((unsigned int)(in[2])>>8)+((unsigned int)(in[3])));
    }
}

/* readin a binary float with conditional unsigned char-swabbing */
static float map_in_float(FILE *fp, unsigned int swab)
{
    unsigned char in[4];
    union
    {
        unsigned char buf[4];
        float f;
    } t;
    fread(&(in[0]), 1, 4, fp);
    if (!swab)
    {
        t.buf[0] = in[0];
        t.buf[1] = in[1];
        t.buf[2] = in[2];
        t.buf[3] = in[3];
    }
    else
    {
        t.buf[0] = in[3];
        t.buf[1] = in[2];
        t.buf[2] = in[1];
        t.buf[3] = in[0];
    }
    return t.f;
}

#ifdef TRILINEAR
/* trilinear interpolation given 8 grid points and the fractional
 * distance in between.
 */
static float
trilinear(m000, m001, m010, m011, m100, m101, m110, m111, dx1, dy1, dz1)
float m000, m001, m010, m011, m100, m101, m110, m111, dx1, dy1, dz1;
{
    float dx0, dy0, dz0;
    dx0 = 1.0 - dx1;
    dy0 = 1.0 - dy1;
    dz0 = 1.0 - dz1;
    return (
               dx0*dy0*dz0*m000
               +dx0*dy0*dz1*m001
               +dx0*dy1*dz0*m010
               +dx0*dy1*dz1*m011
               +dx1*dy0*dz0*m100
               +dx1*dy0*dz1*m101
               +dx1*dy1*dz0*m110
               +dx1*dy1*dz1*m111 );
}
#endif // ifdef TRILINEAR

static float
pseudospline(float m000, float m001, float m010, float m011, float m100, float m101, float m110, float m111, float dx1, float dy1, float dz1)
{
    float dx0, dy0, dz0;

    /* does a "pseudo-spline interpolation by weighting corners
       and then doing a linear interpolation */
    if (dx1 > 0.5F)
    {
        dx1 = 1.0F - dx1;
        dx1 = 2.0F * dx1 * dx1;
        dx1 = 1.0F - dx1;
    }
    else
    {
        dx1 = 2.0F * dx1 *dx1;
    }

    if (dy1 > 0.5F)
    {
        dy1 = 1.0F - dy1;
        dy1 = 2.0F * dy1 * dy1;
        dy1 = 1.0F - dy1;
    }
    else
    {
        dy1 = 2.0F * dy1 *dy1;
    }

    if (dz1 > 0.5F)
    {
        dz1 = 1.0F - dz1;
        dz1 = 2.0F * dz1 * dz1;
        dz1 = 1.0F - dz1;
    }
    else
    {
        dz1 = 2.0F * dz1 *dz1;
    }

    dx0 = 1.0F - dx1;
    dy0 = 1.0F - dy1;
    dz0 = 1.0F - dz1;

    return (
               m000 * dx0 * dy0 * dz0
               +m001 * dx0 * dy0 * dz1
               +m010 * dx0 * dy1 * dz0
               +m100 * dx1 * dy0 * dz0
               +m011 * dx0 * dy1 * dz1
               +m110 * dx1 * dy1 * dz0
               +m101 * dx1 * dy0 * dz1
               +m111 * dx1 * dy1 * dz1
               );
}

static bool next_line(FILE *fp,  char *buf, size_t n, std::string &ubuf)
{
    while (fgets(buf, n, fp) != NULL)
    {
        ubuf = buf;
        ubuf = MIToUpper(ubuf);
        MIStringTrim(ubuf, false); // trims space from left
        if (strncmp(ubuf.c_str(), "REMARK", 6) == 0 || ubuf.length() == 0)
        {
            continue;
        }
        return true;
    }
    return false;
}

static int
atom_has_density(MIAtom *a /*, RESIDUE *res*/)
{
    /* ignore Hydrogens with 0.0 B-values (i.e. XPLOR) */
    if (!strncmp("H", a->name(), 1) && a->BValue() == 0.0)
    {
        //fprintf(stderr,"Ignoring %s %s %s B=%0.1f\n",res->type(), res->name(),a->name,a->BValue);
        return false;
    }
    if (a->occ() < 0.00001)
    {
        return false;
    }
    if (a->type() & AtomType::DUMMYATOM)
    {
        return false;
    }
    return true;
}

#undef mdex
#define mdex(ix, iy, iz, nx, ny, nz) (nx*(ny*((iz+64*nz)%nz)+(iy+64*ny)%ny)+(ix+64*nx)%nx)


static void rhocoefs(MIAtom *atom, Residue *res, int nx, int ny, int nz, CMapHeaderBase *mh, float *ae, float *be, float *boxrad, float *dmax)
{
    float d, dx, dy, dz;
    float r, r0;
    int i;
    int type;
    float Bv, occ;

    occ = atom->occ();
    Bv = atom->BValue();
    type = 1; /* default = carbon */
    type = ScattIndex(atom->name(), res->type().c_str());

    if (type==-1)
    {
        type = 1;
        Logger::log("Warning, atom type %s unknown, treated as carbon", atom->name());
    }

    /* derive real space form factor from reciprocal space version */
    ae[1] = (float)(occ*a1[type]*pow(sqrt(4.0*PI/(b1[type]+Bv+Bscale)), 3.0));
    ae[2] = (float)(occ*a2[type]*pow(sqrt(4.0*PI/(b2[type]+Bv+Bscale)), 3.0));
    ae[3] = (float)(occ*a3[type]*pow(sqrt(4.0*PI/(b3[type]+Bv+Bscale)), 3.0));
    ae[4] = (float)(occ*a4[type]*pow(sqrt(4.0*PI/(b4[type]+Bv+Bscale)), 3.0));
    ae[5] = (float)(occ*co[type]*pow(sqrt(4.0*PI/(         Bv+Bscale)), 3.0));
    be[1] = (float)(4.0*PI*PI/(b1[type]+Bv+Bscale));
    be[2] = (float)(4.0*PI*PI/(b2[type]+Bv+Bscale));
    be[3] = (float)(4.0*PI*PI/(b3[type]+Bv+Bscale));
    be[4] = (float)(4.0*PI*PI/(b4[type]+Bv+Bscale));
    be[5] = (float)(4.0*PI*PI/(         Bv+Bscale));
    /* calculate box bounds */
    /* the density at d = 0.0 */
    r0 = ae[1] + ae[2] + ae[3] + ae[4] + ae[5];
    /* search for point where density is 1/100 of r0 */
    for (i = 1; i < 100; i++)
    {
        d = (float)i * 0.3f;
        d *= d;
        r = ae[1]*exp(-be[1]*d) + ae[2]*exp(-be[2]*d) + ae[3]*exp(-be[3]*d)
            +ae[4]*exp(-be[4]*d) + ae[5]*exp(-be[5]*d);
        if (r0/r > 500.0)
        {
            break;
        }
    }
    *dmax = (float)i*0.3f;
    dx = dy = dz = *dmax;
    transform(mh->ctof, &dx, &dy, &dz);
    boxrad[0] = fabs(dx)*(float)nx;
    boxrad[1] = fabs(dy)*(float)ny;
    boxrad[2] = fabs(dz)*(float)nz;
    /* check against dmax^2 to save a square root */
    *dmax *= *dmax;
}

#ifdef _WIN32
static double drand48(void)
{
    return ((double)rand()/(double)RAND_MAX);
}

#endif

static void addrho(MIAtom *atom, micomplex *cmap, int nx, int ny, int nz, CMapHeaderBase *mh, Residue *res, float shake_coords = 0.0)
{
    int ix, iy, iz;
    int xl, xu, yl, yu, zl, zu;
    float sx, sy, sz;
    float fx, fy, fz;
    float d, dx, dy, dz;
    float fix, fiy, fiz;
    float ffx, ffy, ffz;
    int isx, isy, isz;
    float r;
    int is;
    float dmax, ax, ay, az;
    float ae1[6], be1[6], boxrad1[3];
    float *ae, *be, *boxrad;

    if (!atom_has_density(atom /*,res*/))
    {
        return;
    }
    ae = ae1;
    be = be1;
    boxrad = boxrad1;
    rhocoefs(atom, res, nx, ny, nz, mh, ae, be, boxrad, &dmax);
    ax = atom->x();
    ay = atom->y();
    az = atom->z();
    if (shake_coords > 0.0)
    {
        //extern double drand48();
        ax += (float)drand48()*shake_coords*2.0f - shake_coords;
        ay += (float)drand48()*shake_coords*2.0f - shake_coords;
        az += (float)drand48()*shake_coords*2.0f - shake_coords;
    }
    /* transform to fractional coords */
    fx = ax;
    fy = ay;
    fz = az;
    transform(mh->ctof, &fx, &fy, &fz);
    fx *= (float)nx;
    fy *= (float)ny;
    fz *= (float)nz;
    xl = ROUND(fx-boxrad[0]);
    xu = ROUND(fx+boxrad[0]);
    yl = ROUND(fy-boxrad[1]);
    yu = ROUND(fy+boxrad[1]);
    zl = ROUND(fz-boxrad[2]);
    zu = ROUND(fz+boxrad[2]);


    /* loop through box */
    for (iz = zl; iz <= zu; iz++)
    {
        ffz = (float)iz/(float)nz;
        for (iy = yl; iy <= yu; iy++)
        {
            ffy = (float)iy/(float)ny;
            for (ix = xl; ix <= xu; ix++)
            {
                ffx = (float)ix/(float)nx;
                fix = ffx;
                fiy = ffy;
                fiz = ffz;
                transform(mh->ftoc, &fix, &fiy, &fiz);

                /* distance squared in A from center of atom */
                dx = fix - ax;
                dy = fiy - ay;
                dz = fiz - az;
                d = dx*dx + dy*dy + dz*dz;

                /* check if within sphere */
                if (d > dmax)
                {
                    continue;
                }

                /* calculate rho at this point.
                 */
                r = ae[1]*exp(-be[1]*d)
                    +ae[2]*exp(-be[2]*d)
                    +ae[3]*exp(-be[3]*d)
                    +ae[4]*exp(-be[4]*d)
                    +ae[5]*exp(-be[5]*d);
                /* add rho to symmetry positions */
                for (is = 0; is < mh->nsym; is++)
                {
                    sx = ffx * mh->symops[0][0][is]
                         +ffy * mh->symops[0][1][is]
                         +ffz * mh->symops[0][2][is]
                         +mh->symops[0][3][is];
                    sy = ffx * mh->symops[1][0][is]
                         +ffy * mh->symops[1][1][is]
                         +ffz * mh->symops[1][2][is]
                         +mh->symops[1][3][is];
                    sz = ffx * mh->symops[2][0][is]
                         +ffy * mh->symops[2][1][is]
                         +ffz * mh->symops[2][2][is]
                         +mh->symops[2][3][is];
                    isx = ROUND(sx * (float)nx);
                    isy = ROUND(sy * (float)ny);
                    isz = ROUND(sz * (float)nz);
                    cmap[mdex(isx, isy, isz, nx, ny, nz)].r += r;
                }
            }
        }
    }
}

static micomplex *buildrho(Residue *reslist, CMapHeaderBase *mh)
{
    micomplex *cmap;
    int nx, ny, nz;
    int i, doingsym;
    Residue *res;
    MIAtom *a;
    int n = 0;
    sfinit();

    Bscale = std::max(exp(mh->resmin+.1)/2.0, 0.0) + 2.0;

    /* first figure out the size of the map given the cell
     * and the resolution desired
     */
    nx = ROUND(3.0*mh->a/mh->resmin);
    ny = ROUND(3.0*mh->b/mh->resmin);
    nz = ROUND(3.0*mh->c/mh->resmin);
    if (nx > 120)
    {
        /* if grid large make coarser and add to Bscale to compensate*/
        nx = ROUND(2.3*mh->a/mh->resmin);
        Bscale += 10.0;
    }
    if (nx > 150)
    {
        nx = ROUND(2.0*mh->a/mh->resmin);
        Bscale += 10.0;
    }
    if (ny > 120)
    {
        ny = ROUND(2.3*mh->b/mh->resmin);
        Bscale += 10.0;
    }
    if (ny > 150)
    {
        ny = ROUND(2.0*mh->b/mh->resmin);
        Bscale += 10.0;
    }
    if (nz > 120)
    {
        nz = ROUND(2.3*mh->c/mh->resmin);
        Bscale += 10.0;
    }
    if (nz > 150)
    {
        nz = ROUND(2.0*mh->c/mh->resmin);
        Bscale += 10.0;
    }

    nx = MIMapFactor(nx, FFT_PRIME, EVEN, 2);
    ny = MIMapFactor(ny, FFT_PRIME, EVEN, 2);
    nz = MIMapFactor(nz, FFT_PRIME, EVEN, 2);
    Logger::log("sfFFT: A B-value of %0.1f will be added", Bscale);
    Logger::log("sfFFT: FFT grid: nx, ny, nz = %d %d %d", nx, ny, nz);

    cmap = (micomplex*)malloc(nx*ny*nz*sizeof(micomplex));
    if (cmap == NULL)
    {
        Logger::log("sfFFT: unable to allocate map space");
        return (NULL);
    }
    for (i = 0; i < nx*ny*nz; i++)
    {
        cmap[i].r = 0.0;
        cmap[i].i = 0.0;
    }

    doingsym = false;
    res = reslist;
    while (res != NULL)
    {
        if (strcmp(res->type().c_str(), "BND") == 0 && res->name().size() > 0 && res->name()[0] == '#')
        {
            res = res->next();
            continue;
        }
        for (i = 0; i < res->atomCount(); i++)
        {
            a = res->atom(i);
            addrho(a, cmap, nx, ny, nz, mh, res);
            n++;
        }
        //printf(".");
        //fflush(stdout);
        res = res->next();
    }

    Logger::log("sfFFT: built rho for %d atoms", n);

    /* copy these to mh for later steps */
    mh->nx = nx;
    mh->ny = ny;
    mh->nz = nz;

    return (cmap);
}

/* invert a P1 map and copy structure factors to refl */
static int InvertMap(micomplex x[], CMapHeaderBase *mh, std::vector<CREFL> &refl)
{
    int j;
    int ix, iy, iz;
    long int d[5];
    long int nx, ny, nz;
    int n = 0;
    float fact;
    float Vfact;
    nx = mh->nx;
    ny = mh->ny;
    nz = mh->nz;

    /*  transform fast dimension */
    /*  transforms on x. */
    d[0] = 2*nz*nx*ny;
    d[1] = 2;
    d[2] = d[0];
    d[3] = d[0];
    d[4] = 2*nx;

    cmplft_((float*)x, (float*)&(x[0].i), &nx, d);

    /*  transform medium dimension */
    /*  calculates fourier transforms on y  */
    d[1] = 2*nx;
    d[2] = 2*nx*ny;
    d[3] = d[1]  ;
    d[4] = 2;

    cmplft_((float*)x, (float*)&(x[0].i), &ny, d);

    /*  transform slow dimension */
    /*  transforms on z. */
    d[1] = 2*nx*ny ;
    d[2] = d[0];
    d[3] = d[1];
    d[4] = 2;

    cmplft_((float*)x, (float*)&(x[0].i), &nz, d);

    Vfact = Volume(mh->a, mh->b, mh->c, mh->alpha, mh->beta, mh->gamma)
            /((float)nx*(float)ny*(float)nz);

    /* copy out the fc values into the refls */
    for (unsigned int i = 0; i < refl.size(); i++)
    {
        if (refl[i].sthol > 0.5/mh->resmin
            || refl[i].sthol < 0.5/mh->resmax)
        {
            continue;
        }
        ix = refl[i].ind[0];
        iy = refl[i].ind[1];
        iz = refl[i].ind[2];
        j = mdex(ix, iy, iz, nx, ny, nz);
        if (j < 0 || j > 2*nx*ny*nz)
        {
            continue;
        }
        /* later steps are used to figure out F and phase */
        refl[i].acalc = x[j].r;
        /* for some reason b = -b in output of cmplft */
        refl[i].bcalc = -x[j].i;
        /* correction for Bscale */
        fact = (float)exp(Bscale*refl[i].sthol*refl[i].sthol);
        /* correction for grid  Volume/ nx*ny*nz */
        fact *= Vfact;
        refl[i].acalc *= fact;
        refl[i].bcalc *= fact;
        n++;
    }
#ifdef DEBUGFFTINV
    for (ix = 0; ix <= nx; ix++)
    {
        printf("\nSection ix = %d\n", ix);
        for (iy = 0; iy <= ny; iy++)
        {
            printf("iy =%5d: ", iy);
            for (iz = 0; iz <= nz; iz++)
            {
                float r;
                j = mdex(ix, iy, iz, nx, ny, nz);
                r = sqrt(x[j].re*x[j].re+x[j].im*x[j].im);
                printf("%5d", (int)r);
            }
            printf("\n");
        }
    }
    Logger::log("sfFFT: %d reflections calculated", n);
#endif // ifdef DEBUGFFTINV
    return (n);
}

static int symm_mh(float x, float y, float z, float *xp, float *yp, float *zp, CMapHeaderBase *mh, int jth)
{
    *xp = x * mh->symops[0][0][jth] + y * mh->symops[0][1][jth] + z * mh->symops[0][2][jth] + mh->symops[0][3][jth];
    *yp = x * mh->symops[1][0][jth] + y * mh->symops[1][1][jth] + z * mh->symops[1][2][jth] + mh->symops[1][3][jth];
    *zp = x * mh->symops[2][0][jth] + y * mh->symops[2][1][jth] + z * mh->symops[2][2][jth] + mh->symops[2][3][jth];
    return 1;
}

EMapBase::EMapBase()
    : visible(true),
      FreeRSet(false),
      FOMsValid(false),
      FosValid(false),
      FcsValid(false),
      PhicsValid(false)
{

    flog = NULL;
    orthgrid = false;
    gridspacing = 1.0;
    UseNCR = false;
    settings = new MapSettingsBase();
    mapheader = new CMapHeaderBase();
    settings->load();
    mapmin = -250;
    mapmax = 300;
    changed = false;
    mapnumber = 1;
    CurrentAtoms = NULL;
    mapGrid = 1;
    version_string  = std::string(__FILE__);
    version_string += std::string(__DATE__);
    _fostr = _fcstr = _fomstr = _phistr = _sigfstr = _freeRstr = "";
}

EMapBase::~EMapBase()
{
    if (flog)
    {
        fclose(flog);
    }
    flog = NULL;
    delete settings;
    delete mapheader;
}

bool EMapBase::IsCif(const char *pathname)
{
    /* open the file and look for a datablock */
    FILE *fp;
    char buf[512];
    int ndata = 0;
    if ((fp = fopen(pathname, "r")) == NULL)
    {
        return false;

    }
    std::string sbuf;
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        sbuf = buf;
        MIStringTrim(sbuf, false); // trims space from left
        sbuf = MIToLower(sbuf);
        /* should we convert to lowercase ? */
        if (strncmp(sbuf.c_str(), "data_", 5) == 0 && sbuf.length() > 5)
        {
            ndata++;
        }
    }
    fclose(fp);
    return ndata != 0;
}

bool EMapBase::InterpretColumnLabelString(const char *s)
{
    std::string str(s);
    char buf[128]; // column names have a hard-coded 32 character limit,
    std::string::size_type pos;

    pos = str.find("FO=");
    if (pos!=std::string::npos && sscanf(&str[pos], "FO=%s", buf)==1)
    {
        _fostr = buf;
        fColumnName = _fostr;
    }

    pos = str.find("FC=");
    if (pos!=std::string::npos && sscanf(&str[pos], "FC=%s", buf)==1)
        _fcstr = buf;

    pos = str.find("FOM=");
    if (pos!=std::string::npos && sscanf(&str[pos], "FOM=%s", buf)==1)
        _fomstr = buf;

    pos = str.find("PHI=");
    if (pos!=std::string::npos && sscanf(&str[pos], "PHI=%s", buf)==1)
        _phistr = buf;

    pos = str.find("SIGF=");
    if (pos!=std::string::npos && sscanf(&str[pos], "SIGF=%s", buf)==1)
        _sigfstr = buf;

    pos = str.find("FREE=");
    if (pos!=std::string::npos && sscanf(&str[pos], "FREE=%s", buf)==1)
        _freeRstr = buf;
    return true;
}

bool EMapBase::LoadMapPhaseFile(const char *s, bool require_fo)
{
    if (IsCif(s))
    {
        return LoadCIFMap(s) != 0;
    }
    else if (IsCCP4MTZFile(s))
    {
        return LoadMTZMap(s, require_fo) != 0;
    }
    else if (IsScaFile(s))
    {
        // disbling support for this, the reader is broken: doesn't properly handle centric reflections
        // LoadMap(s, EMapBase::HKL_scaI) > 0;
        return false;
    }
    else if (IsFinFile(s))
    {
        return LoadMap(s, EMapBase::XtalView_fin) > 0;
    }
    else if (IsRefFile(s))
    {
        // disbling support for this, the reader is broken: doesn't properly handle centric reflections
        // return LoadMap(s, EMapBase::DTrek_ref) > 0;
        return false;
    }

    return LoadMap(s, EMapBase::XtalView_phase) > 0;
}


void GetTypes(std::vector<unsigned int> &types,
              bool has_fo, bool has_fc, bool has_fom, bool has_phi)
{
    types.clear();
    if (has_fo && has_phi)
    {
        types.push_back(MIMapType::DirectFFT);
        types.push_back(MIMapType::Fo);
    }

    if (has_fc && has_phi)
    {
        types.push_back(MIMapType::Fc);
    }

    if (has_fo && has_fc && has_phi)
    {
        types.push_back(MIMapType::TwoFoFc);
        types.push_back(MIMapType::FoFc);
    }

    if (has_fo && has_phi)
        types.push_back(MIMapType::FoFo);

    if (has_fo && has_fom && has_phi)
        types.push_back(MIMapType::Fofom);

    if (has_fo && has_fc && has_phi)
    {
        types.push_back(MIMapType::ThreeFoTwoFc);
        types.push_back(MIMapType::FiveFoThreeFc);
        types.push_back(MIMapType::TwoFoFcSigmaA);
        types.push_back(MIMapType::FoFcSigmaA);
    }

}


bool FindCIFIndices(FILE *fp, int &nrefl_start,
                    int &hindex, int &kindex, int &lindex,
                    int &fofoindex, int &foindex, int &fcindex,
                    int &phsindex, int &sigmaindex, int &sigmasquareindex)
{
    char buf[512];

    int nline = -1;
    int nfound = 0;
    // State: 0 = ignore, 1 = loop, 2 = reflnLoop, 3 = complete
    int state = 0;
    bool reflnLoopFound = false;
    while (state < 3 && fgetsl(buf, sizeof buf, fp) != NULL)
    {
        ++nline;
        if (buf[0] == '\0' || buf[0] == '#')
        {
            continue;
        }
        if (state == 0)
        {
            if (strncmp(buf, "loop_", 5) == 0)
            {
                state = 1;
                continue;
            }
        }
        if (state == 1)
        {
            if (strncmp(buf, "_refln", 5) == 0)
            {
                state = 2;
                reflnLoopFound = true;
            }
            else
            {
                state = 0;
            }
        }
        if (state == 2)
        {
            if (strncmp(buf, "_refln", 5) == 0)
            {
                if (strstr(buf, "index_h"))
                {
                    hindex = nfound;
                }
                else if (strstr(buf, "index_k"))
                {
                    kindex = nfound;
                }
                else if (strstr(buf, "index_l"))
                {
                    lindex = nfound;
                }
                else if (strstr(buf, "f_squared_meas") && !strstr(buf, "sigma"))
                {
                    fofoindex = nfound;
                }
                else if (strstr(buf, "intensity_meas") && !strstr(buf, "sigma"))
                {
                    fofoindex = nfound;
                }
                else if (strstr(buf, "f_meas") && !strstr(buf, "f_meas_sigma"))
                {
                    foindex = nfound;
                }
                else if (strstr(buf, "f_squared_sigma"))
                {
                    sigmasquareindex = nfound;
                }
                else if (strstr(buf, "intensity") && strstr(buf, "sigma"))
                {
                    sigmasquareindex = nfound;
                }
                else if (strstr(buf, "f_meas_sigma"))
                {
                    sigmaindex = nfound;
                }
                else if (strstr(buf, "f_sigma"))
                {
                    sigmaindex = nfound;
                }
                else if (strstr(buf, "f_meas_au") && !strstr(buf, "f_meas_au_sigma"))
                {
                    foindex = nfound;
                }
                else if (strstr(buf, "f_calc"))
                {
                    fcindex = nfound;
                }
                else if (strstr(buf, "phase_calc"))
                {
                    phsindex = nfound;
                }
                nfound++;
                continue;
            }
            else
            {
                nrefl_start = nline;
                state = 3;
            }
        }
    }

    if (!reflnLoopFound)
    {
        Logger::message("CIFFile: Couldn't find a reflection loop - must abort");
        return 0;
    }
    if (nfound < 4)
    {
        Logger::message("CIFFile: Reflection loop has too few items - must abort");
        return 0;
    }
    if (hindex < 0 || kindex < 0 || lindex < 0)
    {
        Logger::message("CIFFile: Reflection loop has no indices - must abort");
        return 0;
    }

    if (fofoindex < 0 && foindex < 0)
    {
        if (fcindex >= 0)
        {
            foindex = fcindex;
            Logger::message("CIFFile: Reflection loop has no Fo item - Fo set to Fc");
        }
        else
        {
            Logger::message("CIFFile: Reflection loop has no Fo item - must abort");
            return 0;
        }
    }
    return 1;
}

bool EMapBase::GetMapTypes(std::vector<unsigned int> &maptypes,
                           std::vector<unsigned int> &with_fc_maptypes)
{

    if (refls.size()==0)
        return false;
    GetTypes(maptypes, true, FcsValid, HasFOMs(), PhicsValid);
    GetTypes(with_fc_maptypes, true, true, HasFOMs(), true);
    return true;
}

bool EMapBase::CanDoMapType(unsigned int type)
{
    std::vector<unsigned int> maptypes;
    std::vector<unsigned int> with_fc_maptypes;
    if (!GetMapTypes(maptypes, with_fc_maptypes))
        return false;

    for (size_t i = 0; i<maptypes.size(); ++i)
        if (maptypes[i]==type)
            return true;
    return false;
}


bool EMapBase::GetPossibleMapTypes(const std::string &pathname,
                                   std::vector<unsigned int> &maptypes,
                                   std::vector<unsigned int> &with_fc_maptypes)
{
    if (!pathname.size())
        return 0;

    //CIF file
    if (IsCif(pathname.c_str()))
    {
        int hindex = (-1), kindex = (-1), lindex = (-1), fofoindex = (-1),
            foindex = (-1), fcindex = (-1), phsindex = (-1),
            sigmaindex = (-1), sigmasquareindex = (-1);
        int nrefl_start = 0;
        FILE *fp = fopen(pathname.c_str(), "r");
        if (!fp)
            return 0;
        if (!FindCIFIndices(fp, nrefl_start,
                            hindex, kindex, lindex,
                            fofoindex, foindex, fcindex, phsindex,
                            sigmaindex, sigmasquareindex))
        {

            fclose(fp);
            return 0;
        }
        fclose(fp);
        GetTypes(maptypes, foindex!=-1, fcindex!=-1, false, phsindex != -1);
        GetTypes(with_fc_maptypes, foindex!=-1, true, false, true);
        return true;
    }

    //Sca file
    if (IsScaFile(pathname.c_str()))
    {
        GetTypes(maptypes, true, false, false, false);
        GetTypes(with_fc_maptypes, true, true, false, true);
        return true;
    }

    //Fin file
    if (IsFinFile(pathname.c_str()))
    {
        GetTypes(maptypes, true, false, false, false);
        GetTypes(with_fc_maptypes, true, true, false, true);
        return true;
    }

    //RefFile
    if (IsRefFile(pathname.c_str()))
    {
        GetTypes(maptypes, true, false, false, false);
        GetTypes(with_fc_maptypes, true, true, false, true);
        return true;
    }

    //MTZ file: note this returns which maptypes are *possible* assuming correct assigment of FO, FC, etc
    if (IsCCP4MTZFile(pathname.c_str()))
    {
        std::vector<std::string> labels;
        std::vector<char> types;
        if (!GetMTZColumnInfo(pathname, labels, types))
            return false;

        bool possibly_has_fo = false;
        bool possibly_has_fc = false;
        bool possibly_has_fom = false;
        bool possibly_has_phase = false;

        for (size_t i = 0; i<types.size(); ++i)
        {
            switch (types[i])
            {
            case 'F':
            case 'G':
            case 'J':
            case 'K':
                possibly_has_fo = true;
                possibly_has_fc = true;
                break;
            case 'W':
                possibly_has_fom = true;
                break;
            case 'P':
                possibly_has_phase = true;
                break;
            default:
                break;
            }
        }

        GetTypes(maptypes, possibly_has_fo, possibly_has_fc, possibly_has_fom, possibly_has_phase);
        GetTypes(with_fc_maptypes, possibly_has_fo, true, possibly_has_fom, true);
        return true;
    }

    GetTypes(maptypes, true, true, false, true);
    GetTypes(with_fc_maptypes, true, true, false, true);
    return true;
}



float EMapBase::CalcRMS()
{
    float sum = 0.0F, rms = 0.0F;
    predictedAsDifferenceMap = false;
    if (map_points.size() == 0)
    {
        return 0.0;
    }
    //unsigned int positivePoints = 0;
    unsigned int negativePoints = 0;
    //unsigned int zeroPoints = 0;
    mapmin = mapmax = map_points[0];
    for (unsigned int i = 0; i < map_points.size(); i++)
    {
        float point = map_points[i];
        sum += point*point;
        if (point > mapmax)
        {
            mapmax = point;
        }
        if (point < mapmin)
        {
            mapmin = point;
        }
        if (point < 0.0)
        {
            ++negativePoints;
        }
    }
    double mapPointRatio = (map_points.size() - negativePoints) / (double)negativePoints;
    if (mapPointRatio > 0.85 )
    {
        predictedAsDifferenceMap = true;
    }
    Logger::debug("map points: %u neg, %u pos, %u total", negativePoints, (map_points.size() - negativePoints), map_points.size());
    Logger::log("map points: ratio %0.3f, diff map %s", mapPointRatio, predictedAsDifferenceMap ? "true" : "false");
    rms = (float)sqrt(sum/(double)map_points.size());
    return rms;
}

void EMapBase::ScaleMap(float rms)
{
    if (rms > 0.0)
    {
        scale = 50.0F/rms;
        for (unsigned int i = 0; i < map_points.size(); i++)
        {
            map_points[i] *= scale;
        }
        mapmax *= scale;
        mapmin *= scale;
    }
    else
    {
        scale = 1.0F;
    }
    Logger::log("The rms value of %f was scaled by %f to put the rms value of the map at 50.0\nMaximum= %0.1f  Minimum=%0.1f", rms, scale, mapmax, mapmin);
}

long EMapBase::LoadMap(const char *pathname, int type)
{
    FILE *fp;
    CREFL refl;
    char buf[512];
    int ih, ik, il, ret;
    float fo, fc, phi, f1, f2, sigf1, sigf2, s1, s2;
    // for D*trek ref file
    int h_col = -1;
    int k_col = -1;
    int l_col = -1;
    int Io_col = -1;
    int sIo_col = -1;
    int Iplus_col = -1;
    int Iminus_col = -1;
    int sIplus_col = -1;
    int sIminus_col = -1;
    int maxcol = -1;
    float c[12];
    if ((fp = fopen(pathname, "r")) == NULL)
    {
        Logger::message("Can't open file! Already open? Doesn't exist?");
        return 0;
    }
    if (!fp)
    {
        return 0;
    }
    else
    {
        pathName = pathname;
    }
    refls.clear();
    // clear all the flag values values
    FreeRSet = false;
    FOMsValid = false;
    FosValid = false;
    FcsValid = false;
    PhicsValid = false;

    memset(&refl, 0, sizeof(refl));
    if (type == EMapBase::HKL_scaI || type == EMapBase::HKL_scaF)
    {
        char SymInfoString[200];
        float a, b, c, al, be, ga;
        // the first 2 lines are header info
        fgets(buf, sizeof buf, fp);
        fgets(buf, sizeof buf, fp);
        fgets(buf, sizeof buf, fp);
        ret = sscanf(buf, "%f%f%f%f%f%f%s", &a, &b, &c, &al, &be, &ga, SymInfoString);
        if (ret == 7)
        {
            mapheader->a = a;
            mapheader->b = b;
            mapheader->c = c;
            mapheader->alpha = al;
            mapheader->beta = be;
            mapheader->gamma = ga;
            mapheader->FindSpacegroup(SymInfoString);
        }
    }
    else if (type == EMapBase::DTrek_ref)
    {
        float a, b, c, al, be, ga;
        char key[sizeof buf];
        int nI, nF, nS, nC;
        int ncol = 0;
        // first line with 4 numbers
        fgets(buf, sizeof buf, fp);
        ret = sscanf(buf, "%d%d%d%d", &nI, &nF, &nS, &nC);
        if (ret < 3)
        {
            Logger::log("EMapBase::LoadMap:Error1: Not a properly formed D*Trek ref file?");
            return 0;
        }
        if (ret == 3)
        {
            nC = 0;
        }
        if (nI+nF > 12)
        {
            Logger::log("EMapBase::LoadMap:Error4: Unsupported format for D*Trek ref file");
            return 0;
        }
        for (int i = 0; i < (nI+nC+nF); i++)
        {
            if (fgets(buf, sizeof buf, fp) == NULL)
            {
                Logger::log("EMapBase::LoadMap:Error2: Not a properly formed D*Trek ref file?");
                return 0;
            }
            sscanf(buf, "%s", key);
            if (strncmp(buf, "CRYSTAL_UNIT_CELL", 17) == 0)
            {
                char *equal = strchr(buf, '=');
                if (equal)
                {
                    ret = sscanf(equal+1, "%f%f%f%f%f%f", &a, &b, &c, &al, &be, &ga);
                    if (ret == 6)
                    {
                        mapheader->a = a;
                        mapheader->b = b;
                        mapheader->c = c;
                        mapheader->alpha = al;
                        mapheader->beta = be;
                        mapheader->gamma = ga;
                    }
                }
            }
            else if (strncmp(buf, "CRYSTAL_SPACEGROUP", 18) == 0)
            {
                char *equal = strchr(buf, '=');
                int spgpno;
                if (equal)
                {
                    ret = sscanf(equal+1, "%d", &spgpno);
                    if (ret == 1)
                    {
                        mapheader->spgpno = spgpno;
                        mapheader->SetSymmOps();
                    }
                }
            }
            else if (strcmp(key, "nH") == 0)
            {
                h_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "nK") == 0)
            {
                k_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "nL") == 0)
            {
                l_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "fIntensity") == 0)
            {
                Io_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "fSigmaI") == 0)
            {
                sIo_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "fIntensity+") == 0)
            {
                Iplus_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "fIntensity-") == 0)
            {
                Iminus_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "fSigmaI+") == 0)
            {
                sIplus_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
            else if (strcmp(key, "fSigmaI-") == 0)
            {
                sIminus_col = ncol;
                if (ncol > maxcol)
                {
                    maxcol = ncol;
                }
                ncol++;
            }
        }
        // check to see if we have enought to read this file
        if (maxcol < 2 || h_col == (-1) || k_col == (-1) || l_col == (-1))
        {
            Logger::log("EMapBase::LoadMap:Error3: Can't find 3 indices in D*Trek ref file");
            return 0;
        }
        if (Io_col == (-1) && Iplus_col == (-1) && Iminus_col == (-1))
        {
            Logger::log("EMapBase::LoadMap:Error3: Can't find intenisty data in D*Trek ref file");
            return 0;
        }
    }
    else if (type == EMapBase::XtalView_phase)
    {
        FcsValid = true;
        PhicsValid = true;
    }
    else if (type == EMapBase::XtalView_fin)
    {
    }
    sthol(1, 1, 1, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 1);
    refls_stholmax = -1.0F;
    refls_stholmin = 999999.0F;
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        if (buf[0] == '#' || buf[0] == ';' || buf[0] == '!')
        {
            continue;
        }
        if (type == EMapBase::XtalView_phase)
        {
            if (sscanf(buf, "%d%d%d%f%f%f", &ih, &ik, &il, &fo, &fc, &phi) == 6)
            {
                refl.ind[0] = ih;
                refl.ind[1] = ik;
                refl.ind[2] = il;
                refl.fo = fo;
                refl.fc = fc;
                refl.phi = phi;
                refl.sthol = sthol(ih, ik, il, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
                if (refl.sthol > refls_stholmax)
                {
                    refls_stholmax = refl.sthol;
                }
                if (refl.sthol < refls_stholmin)
                {
                    refls_stholmin = refl.sthol;
                }
                refls.push_back(refl);
            }
        }
        else if (type == EMapBase::XtalView_fin)
        {
            if (sscanf(buf, "%d%d%d%f%f%f%f", &ih, &ik, &il, &f1, &sigf1, &f2, &sigf2) == 7)
            {
                refl.ind[0] = ih;
                refl.ind[1] = ik;
                refl.ind[2] = il;
                if (f1 == 0.0)
                {
                    f1 = f2;
                    sigf1 = FLT_MAX;
                }
                if (f2 == 0.0)
                {
                    f2 = f1;
                    sigf2 = FLT_MAX;
                }
                if (sigf1 == 0.0)
                {
                    sigf1 = FLT_MAX;
                }
                if (sigf2 == 0.0)
                {
                    sigf2 = FLT_MAX;
                }
                if (f1 == 0.0 && f2 == 0.0)
                {
                    continue;
                }
                refl.fo = (f1+f2)/2.0F;
                refl.sigma = 1.0F/((1.0F/sigf1)+(1.0F/sigf2));
                refl.fc = refl.fo;
                refl.phi = 0.0;
                refl.sthol = sthol(ih, ik, il, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
                if (refl.sthol > refls_stholmax)
                {
                    refls_stholmax = refl.sthol;
                }
                if (refl.sthol < refls_stholmin)
                {
                    refls_stholmin = refl.sthol;
                }
                refls.push_back(refl);
            }
        }
        else if (type == EMapBase::DTrek_ref)
        {
            ret = sscanf(buf, "%f%f%f%f%f%f%f%f%f%f%f%f", &c[0], &c[1], &c[2], &c[3], &c[4], &c[5],
                         &c[6], &c[7], &c[8], &c[9], &c[10], &c[11]);
            if (ret < maxcol+1)
            {
                continue;
            }
            ih = (int)c[h_col];
            ik = (int)c[k_col];
            il = (int)c[l_col];
            if (Io_col == -1)
            {
                if (Iplus_col == -1)
                {
                    f1 = 0.0;
                }
                else
                {
                    f1 = c[Iplus_col];
                }
                if (Iminus_col == -1)
                {
                    f2 = 0.0;
                }
                else
                {
                    f2 = c[Iminus_col];
                }
                if (sIplus_col == -1)
                {
                    s1 = 0.0;
                }
                else
                {
                    s1 = c[sIplus_col];
                }
                if (sIminus_col == -1)
                {
                    s2 = 0.0;
                }
                else
                {
                    s2 = c[sIminus_col];
                }
            }
            else
            {
                f1 = c[Io_col];
                if (sIo_col == -1)
                {
                    s1 = 10*f1;
                }
                f2 = 0.0;
            }
            /* get rid of negative I's */
            if (f1 == 0.0)
            {
                f1 = 0.0;
                s1 = FLT_MAX;
            }
            if (f2 == 0.0)
            {
                f2 = 0.0;
                s2 = FLT_MAX;
            }
            /* skip empty reflections (i.e. negative I's ) */
            if (f1 == 0.0 && f2 == 0.0)
            {
                continue;
            }
            if (s1 <= 0.0 || f1 <= 0.0)
            {
                f1 = 0.0;
                s1 = FLT_MAX;
            }
            else
            {
                f1 = sqrt(fabs(f1));
                s1 = s1/2.0f/f1;
            }
            if (s2 <= 0.0 || f2 <= 0.0)
            {
                f2 = 0.0;
                s2 = FLT_MAX;
            }
            else
            {
                f2 = sqrt(fabs(f2));
                s2 = s2/2.0f/f2;
            }
            refl.ind[0] = ih;
            refl.ind[1] = ik;
            refl.ind[2] = il;
            refl.fo = (f1+f2)/2.0F;
            refl.sigma = 1.0F/((1.0F/s1)+(1.0F/s2));
            refl.fc = refl.fo;
            refl.phi = 0.0;
            refl.sthol = sthol(ih, ik, il, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
            if (refl.sthol > refls_stholmax)
            {
                refls_stholmax = refl.sthol;
            }
            if (refl.sthol < refls_stholmin)
            {
                refls_stholmin = refl.sthol;
            }
            refls.push_back(refl);
        }
        else if (type == EMapBase::HKL_scaI || type == EMapBase::HKL_scaF)
        {
            ret = sscanf(buf, "%d%d%d%f%f%f%f", &ih, &ik, &il, &f1, &s1, &f2, &s2);
            /* only 5 columns scanned - rest were blank
             * sscanf will leave old values so zero to prevent
             * error of using last value */
            if (ret == 5)
            {
                f2 = 0.0;
            }
            /* get rid of negative I's */
            if (f1 == 0.0)
            {
                f1 = 0.0;
                s1 = FLT_MAX;
            }
            if (f2 == 0.0)
            {
                f2 = 0.0;
                s2 = FLT_MAX;
            }
            /* skip empty reflections (i.e. negative I's ) */
            if (f1 == 0.0 && f2 == 0.0)
            {
                continue;
            }
            if (s1 <= 0.0 || f1 < 0.0)
            {
                f1 = 0.0;
                s1 = FLT_MAX;
            }
            else if (type == EMapBase::HKL_scaI)
            {
                f1 = sqrt(fabs(f1));
                s1 = s1/2.0f/f1;
            }
            if (s2 <= 0.0 || f2 < 0.0)
            {
                f2 = 0.0;
                s2 = FLT_MAX;
            }
            else if (type == EMapBase::HKL_scaI)
            {
                f2 = sqrt(fabs(f2));
                s2 = s2/2.0f/f2;
            }
            refl.ind[0] = ih;
            refl.ind[1] = ik;
            refl.ind[2] = il;
            //FIXME: why is this unconditionally averaging? what if f2 was 0?
            refl.fo = (f1+f2)/2.0F;
            refl.sigma = 1.0F/((1.0F/s1)+(1.0F/s2));
            refl.fc = refl.fo;
            refl.phi = 0.0;
            refl.sthol = sthol(ih, ik, il, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
            if (refl.sthol > refls_stholmax)
            {
                refls_stholmax = refl.sthol;
            }
            if (refl.sthol < refls_stholmin)
            {
                refls_stholmin = refl.sthol;
            }
            refls.push_back(refl);
        }
    }
    fclose(fp);
    mapheader->resmax = 0.5F/refls_stholmin;
    mapheader->resmin = 0.5F/refls_stholmax;
    mapheader->nx = MIMapFactor((int)(3.0F*mapheader->a/mapheader->resmin), FFT_PRIME, EVEN, 2);
    mapheader->ny = MIMapFactor((int)(3.0F*mapheader->b/mapheader->resmin), FFT_PRIME, ODDOREVEN, 1);
    mapheader->nz = MIMapFactor((int)(3.0F*mapheader->c/mapheader->resmin), FFT_PRIME, ODDOREVEN, 1);
    Logger::log("Loaded %ld reflections", (long int) refls.size());
    return refls.size();
}

long EMapBase::SaveCCP4Phase(const char *pathname)
{
    mmtzfile fileout;
    mmtz_crystal xtl;
    mmtz_dataset set;
    mmtz_column col;
    long nout = 0;
    float fdata[20];
    int flags[20];
    char hdr[200];
    if (refls.size() == 0)
    {
        Logger::message("SaveCCP4Phase:: Error:: No reflections to write out - aborting");
        return 0;
    }

    if ((fileout  = mmtz_open(pathname, "w")) == NULL)
    {
        Logger::message("SaveCCP4Phase:: Error:: Can't open file. Must Abort\n"
                        "Possible sources: Check permissions, disk full, pathname invalid or to\n"
                        "directory you have no permissions in");
        return 0;
    }

    xtl.cell[0] = mapheader->a;
    xtl.cell[1] = mapheader->b;
    xtl.cell[2] = mapheader->c;
    xtl.cell[3] = mapheader->alpha;
    xtl.cell[4] = mapheader->beta;
    xtl.cell[5] = mapheader->gamma;
    set.wavel = 0;

    /* write the initial headers to a new mtz */
    mmtz_init_headers(fileout, "MIFit created reflections file", xtl.cell);

    sprintf(hdr, "SORT  %3i %3i %3i %3i %3i", 0, 0, 0, 0, 0);
    umtz_add_head(fileout, hdr);

    if (mapheader->nsym > 0)
    {
        if (!mapheader->WriteSymmOpsMTZ(fileout))
        {
            Logger::log("Write MTZ header failed - aborting");
            mmtz_close(fileout);
            return 0;
        }
    }
    int iset;
    iset = 0;
    strcpy(col.label, "H");
    strcpy(col.type, "H");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    strcpy(col.label, "K");
    strcpy(col.type, "H");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    strcpy(col.label, "L");
    strcpy(col.type, "H");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    iset = 1;
    strcpy(col.label, "FP");
    strcpy(col.type, "F");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    strcpy(col.label, "SIGFP");
    strcpy(col.type, "Q");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    strcpy(col.label, "FC");
    strcpy(col.type, "F");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    strcpy(col.label, "FMAP");
    strcpy(col.type, "F");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    strcpy(col.label, "PHIC");
    strcpy(col.type, "P");
    sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
    umtz_add_head(fileout, hdr);

    if (HasFOMs())
    {
        strcpy(col.label, "FOM");
        strcpy(col.type, "W");
        sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
        umtz_add_head(fileout, hdr);
    }

    if (HasFreeRFlag())
    {
        strcpy(col.label, "FreeR_flag");
        strcpy(col.type, "I");
        iset = 0;
        sprintf(hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i  ", col.label, col.type, 0.0, 0.0, iset);
        umtz_add_head(fileout, hdr);
    }

    iset = 0;
    strcpy(xtl.pname, "HKL_base");
    strcpy(xtl.xname, "HKL_base");
    strcpy(set.dname, "HKL_base");
    sprintf(hdr, "PROJECT  %6i %-64s", iset, xtl.pname);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "CRYSTAL  %6i %-64s", iset, xtl.xname);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "DATASET  %6i %-64s", iset, set.dname);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "DCELL    %6i %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", iset,
            xtl.cell[0], xtl.cell[1], xtl.cell[2], xtl.cell[3], xtl.cell[4], xtl.cell[5]);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "DWAVEL   %6i %10.5f", iset, set.wavel);
    umtz_add_head(fileout, hdr);

    iset = 1;
    strcpy(xtl.pname, "Exported_Data");
    strcpy(xtl.xname, "Exported_Data");
    strcpy(set.dname, "Exported_Data");
    sprintf(hdr, "PROJECT  %6i %-64s", iset, xtl.pname);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "CRYSTAL  %6i %-64s", iset, xtl.xname);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "DATASET  %6i %-64s", iset, set.dname);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "DCELL    %6i %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", iset,
            xtl.cell[0], xtl.cell[1], xtl.cell[2], xtl.cell[3], xtl.cell[4], xtl.cell[5]);
    umtz_add_head(fileout, hdr);
    sprintf(hdr, "DWAVEL   %6i %10.5f", iset, set.wavel);
    umtz_add_head(fileout, hdr);

    time_t ltime;
    time(&ltime);
    char buf[1024];
    sprintf(buf, "Created by %s on %s", version_string.c_str(), ctime(&ltime));
    umtz_add_hist(fileout, buf);

    // set flags
    for (int j = 0; j < 10; j++)
    {
        flags[j] = 1;
    }
    int index;
    for (size_t i = 0; i < refls.size(); i++)
    {
        index = 0;
        fdata[index++] = (float)refls[i].ind[0];
        fdata[index++] = (float)refls[i].ind[1];
        fdata[index++] = (float)refls[i].ind[2];
        fdata[index++] = refls[i].fo;
        if (refls[i].sigma == 0.0)
        {
            refls[i].sigma = (float)sqrt(refls[i].fo)/3.0F;
        }
        fdata[index++] = refls[i].sigma;
        fdata[index++] = refls[i].fc;
        fdata[index++] = refls[i].coef;
        fdata[index++] = refls[i].phi;
        if (HasFOMs())
        {
            fdata[index++] = refls[i].fom;
        }
        if (HasFreeRFlag())
        {
            fdata[index++] = refls[i].freeRflag;
        }
        mmtz_add_row(fileout, fdata, flags);
        nout++;
    }
    // done in close
    //umtz_rewrite_headers_ranges(fileout);
    mmtz_close(fileout);
    Logger::log("Wrote out %d reflections to %s", (int)nout, pathname);
    return nout;
}

#ifdef _WIN32
#define NAN 0xfffa5a5a
#endif //Define NAN

bool EMapBase::GetMTZColumnInfo(const std::string &fname,
                                std::vector<std::string> &labels,
                                std::vector<char> &types)
{
    labels.clear();
    types.clear();

    if (!IsCCP4MTZFile(fname.c_str()))
        return false;

    mmtzfile filein;
    mmtz_column col;
    mmtz_crystal xtl;
    mmtz_dataset set;
    unsigned int i;

    if ((filein  = mmtz_open(fname.c_str(), "r")) == NULL)
    {
        return false;
    }
    mmtz_get_setxtl(filein,  1,  &set,  &xtl);
    unsigned int num_cols = mmtz_num_cols(filein);
    for (i = 0; i < num_cols; i++)
    {
        mmtz_get_column(filein, i, &col, &set, &xtl);
        labels.push_back(col.label);
        types.push_back(col.type[0]);
    }
    mmtz_close(filein);
    return true;
}

bool EMapBase::HasColumnLabels(const std::string &fname,
                               const std::vector<std::string> &labels)
{

    std::vector<std::string> col_labels;
    std::vector<char> types;
    if (!GetMTZColumnInfo(fname, col_labels, types))
        return false;

    std::set<std::string> label_set;
    for (size_t i = 0; i < col_labels.size(); ++i)
    {
        label_set.insert(col_labels[i]);
    }

    // now search file_labels for each entry in labels;
    for (size_t i = 0; i < labels.size(); ++i)
    {
        if (label_set.find(labels[i])==label_set.end())
        {
            return false;
        }
    }
    return true;
}

bool EMapBase::GetIndicesFromDefaults(unsigned int num_cols,
                                      mmtz_column_ *col,
                                      int &foindex, int &fcindex, int &fomindex,
                                      int &phsindex, int &sigfindex, int &freeRindex)
{
    if (!_fostr.size() && !_fcstr.size() && !_fomstr.size() && !_phistr.size() && !_sigfstr.size() && !_freeRstr.size())
        return false;

    for (unsigned int i = 0; i<num_cols; ++i)
    {
        std::string label(col[i].label);
        if (label==_fostr)
            foindex = i;
        if (label==_fcstr)
            fcindex = i;
        if (label==_fomstr)
            fomindex = i;
        if (label==_phistr)
            phsindex = i;
        if (label==_sigfstr)
            sigfindex = i;
        if (label==_freeRstr)
            freeRindex = i;
    }
    return true;
}


// use the given values for defaults for column labels
void EMapBase::UseColumnLabels(const std::string &fostr, const std::string &fcstr,
                               const std::string &fomstr, const std::string &phistr,
                               const std::string &sigfstr, const std::string &freeRstr)
{
    _fostr = fostr;
    fColumnName = _fostr;
    _fcstr = fcstr;
    _fomstr = fomstr;
    _phistr = phistr;
    _sigfstr = sigfstr;
    _freeRstr = freeRstr;
}

long EMapBase::LoadMTZMap(const char *pathname, bool require_fo)
{
    mmtzfile filein;
    mmtz_column *col;
    float *fdata;
    int *flags;
    mmtz_crystal xtl;
    mmtz_dataset set;
    char buf[4096];
    unsigned int i;
    int Valm = (int)NAN; //missing value indicator
    std::string s;
    int hindex = (-1), kindex = (-1), lindex = (-1), freeRindex = (-1),
        foindex = (-1), fcindex = (-1), phsindex = (-1), fomindex = (-1), sigfindex = (-1);
    if ((filein  = mmtz_open(pathname, "r")) == NULL)
    {
        return 0;
    }
    //int num_datasets = mmtz_num_datasets (filein);
    int num_cols = mmtz_num_cols(filein);
    col = (mmtz_column*)malloc(sizeof(mmtz_column)*num_cols);
    fdata = (float*)malloc(sizeof(float)*num_cols);
    flags = (int*)malloc(sizeof(int)*num_cols);

    // unit cell
    mmtz_get_setxtl(filein,  1,  &set,  &xtl);
    mapheader->a = xtl.cell[0];
    mapheader->b = xtl.cell[1];
    mapheader->c = xtl.cell[2];
    mapheader->alpha = xtl.cell[3];
    mapheader->beta = xtl.cell[4];
    mapheader->gamma = xtl.cell[5];
    mapheader->crystal_name = xtl.xname;
    Logger::log("Unit cell: %f %f %f %f %f %f",
                mapheader->a,
                mapheader->b,
                mapheader->c,
                mapheader->alpha,
                mapheader->beta,
                mapheader->gamma);
    orthog(mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, mapheader->ctof);
    uinv(mapheader->ctof, mapheader->ftoc);

    umtz_hdr *hpt;
    int ret;
    mapheader->nsym = 0;
    for (hpt = umtz_first_head(filein); hpt < umtz_last_head(filein); hpt++)
    {
        if (umtz_keymatch(hpt->data, "VALM") )
        {
            ret = sscanf(hpt->data, "%*s%s", buf);
            if (ret == 1 && strncmp(buf, "NAN", 3) != 0)
            {
                sscanf(hpt->data, "%*s %i", &Valm);
            }
        }
        else
        if (umtz_keymatch(hpt->data, "SYMINF") )
        {
            ret = sscanf(hpt->data, "%*s%*s%*s%*s%i", &mapheader->spgpno);
            if (ret == 1)
            {
                mapheader->SetSymmOps();
            }
        }
    }

    // spacegroup and symmops
    // only need this if syminf was not found
    if (mapheader->nsym == 0)
    {
        int num_symmops = mmtz_get_num_symops(filein);
        std::string symmstring;
        mapheader->SymopsString.clear();
        for (i = 0; i < (unsigned int)num_symmops; i++)
        {
            mmtz_get_symop(filein, i, buf);
            symmstring = symmtrim(buf);
            mapheader->SymopsString.push_back(symmstring);
        }
        if (num_symmops > 0)
        {
            mapheader->ScanSymops();
        }
        else
        {
            mapheader->nsym = 0;
        }
    }

    // resolve columns
    int nfound = 0;
    int numfc = 0, numphi = 0;
    for (i = 0; i < (unsigned int)num_cols; i++)
    {
        mmtz_get_column(filein, i, &col[i], &set, &xtl);
        if (strcmp("H", col[i].label) == 0)
        {
            hindex = i;
            nfound++;
        }
        if (strcmp("K", col[i].label) == 0)
        {
            kindex = i;
            nfound++;
        }
        if (strcmp("L", col[i].label) == 0)
        {
            lindex = i;
            nfound++;
        }
        if (strcmp("PHIC", col[i].label) == 0)
        {
            phsindex = i;
            nfound++;
        }
        if (strcmp("FC", col[i].label) == 0)
        {
            fcindex = i;
            nfound++;
        }
        if (strcmp("F", col[i].label) == 0)
        {
            foindex = i;
            nfound++;
        }
        if (strcmp("SIGF", col[i].label) == 0)
        {
            sigfindex = i;
            nfound++;
        }
        if (strcmp("FP", col[i].label) == 0)
        {
            foindex = i;
            nfound++;
        }
        if (strcmp("SIGFP", col[i].label) == 0)
        {
            sigfindex = i;
            nfound++;
        }
        if (strcmp("FNAT", col[i].label) == 0)
        {
            foindex = i;
            nfound++;
        }
        if (strcmp("SIGFNAT", col[i].label) == 0)
        {
            sigfindex = i;
            nfound++;
        }
        if (strcmp("FOM", col[i].label) == 0)
        {
            fomindex = i;
            nfound++;
        }
        if (strcmp("FreeR_flag", col[i].label) == 0)
        {
            freeRindex = i;
            nfound++;
        }
        if (strcmp("FREE", col[i].label) == 0)
        {
            freeRindex = i;
            nfound++;
        }
    }
    if (hindex == -1 || kindex == -1 || lindex == -1)
    {
        Logger::message("MtzMapFile: reflection index columns not found - must abort");
        mmtz_close(filein);
        free((void*)col);
        free(fdata);
        free(flags);
        return 0;
    }

    if (!GetIndicesFromDefaults(num_cols, col, foindex, fcindex, fomindex, phsindex, sigfindex, freeRindex)
        && !PromptForColumnLabels(num_cols, col, foindex, fcindex, fomindex, phsindex, sigfindex, freeRindex))
    {
        Logger::message("MtzMapFile: Unable to get label indices for reflection loop - must abort");
        mmtz_close(filein);
        free((void*)col);
        free(fdata);
        free(flags);
        return 0;
    }

    //paranoia: map from num_cols or greater back to -1 (invalid)
    if (foindex>=num_cols) foindex = -1;
    if (fcindex>=num_cols) fcindex = -1;
    if (fomindex>=num_cols) fomindex = -1;
    if (phsindex>=num_cols) phsindex = -1;
    if (sigfindex>=num_cols) sigfindex = -1;
    if (freeRindex>=num_cols) freeRindex = -1;


    _fostr = (foindex<0       ? "" : col[foindex].label);
    _fcstr = (fcindex<0       ? "" : col[fcindex].label);
    _fomstr = (fomindex<0     ? "" : col[fomindex].label);
    _phistr = (phsindex<0     ? "" : col[phsindex].label);
    _sigfstr = (sigfindex<0   ? "" : col[sigfindex].label);
    _freeRstr = (freeRindex<0 ? "" : col[freeRindex].label);

    if (hindex < 0 || kindex < 0 || lindex < 0)
    {
        Logger::message("MtzMapFile: Reflection loop has no indices - must abort");
        return 0;
    }

    if (require_fo && foindex < 0)
    {
        if (fcindex >= 0)
        {
            foindex = fcindex;
            Logger::message("MtzMapFile: Reflection loop has no Fo item - Fo set to Fc");
        }
        else
        {
            Logger::message("MtzMapFile: Reflection loop has no Fo item - must abort");
            return 0;
        }
    }

    FreeRSet =   (freeRindex >=0);
    FOMsValid =  (fomindex >= 0);
    FosValid =   (foindex >= 0);
    FcsValid =   (fcindex >= 0);
    PhicsValid = (phsindex >= 0);

    // finally load data
    CREFL r;
    sthol(1, 1, 1, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 1);
    refls_stholmax = -1.0F;
    refls_stholmin = 999999.0F;
    refls.clear();

    for (i = 0; i < (unsigned int)mmtz_num_rows(filein); i++)
    {
        mmtz_get_row(filein, fdata, flags);
        if (!flags[foindex])
        {
            continue;
        }
        memset(&r, 0, sizeof(r));
        r.ind[0] = (int)fdata[hindex];
        r.ind[1] = (int)fdata[kindex];
        r.ind[2] = (int)fdata[lindex];
        // if index == -1 then user slected "NONE" so skip
        // if flags is false, this reflection is not found
        if (foindex != -1 && flags[foindex])
        {
            r.fo = fdata[foindex];
        }
        if (fcindex != -1 && flags[fcindex])
        {
            r.fc = fdata[fcindex];
            numfc++;
        }
        if (fomindex != -1 && flags[fomindex])
        {
            r.fom = fdata[fomindex];
        }
        if (phsindex != -1 && flags[phsindex])
        {
            r.phi = fdata[phsindex];
            numphi++;
        }
        if (sigfindex != -1 && flags[sigfindex])
        {
            r.sigma = fdata[sigfindex];
        }
        if (freeRindex != -1 && flags[freeRindex])
        {
            r.freeRflag = (short int)fdata[freeRindex];
        }
        r.sthol = sthol(r.ind[0], r.ind[1], r.ind[2], mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
        if (r.sthol > refls_stholmax)
        {
            refls_stholmax = r.sthol;
        }
        if (r.sthol < refls_stholmin)
        {
            refls_stholmin = r.sthol;
        }
        refls.push_back(r);
    }

    FcsValid = (numfc > 0);
    PhicsValid = (numphi > 0);

    free((void*)col);
    free(fdata);
    free(flags);
    mmtz_close(filein);
    mapheader->resmax = 0.5F/refls_stholmin;
    mapheader->resmin = 0.5F/refls_stholmax;
    mapName = pathname;
    pathName = pathname;
    return refls.size();
}




long EMapBase::LoadCIFMap(const char *pathname, int datablock)
{
    FILE *fp;
    char buf[512];
    char strings[MAXLOOP][81];
    bool inloop;
    int ih, ik, il;
    float fo;
    int nstart = 0;
    int nend = 0;
    int nsymm = 0;
    int ncell = 0;
    int nblock = 0;
    int nline = 0;
    int nfound, nrefl_start = -1;
    char *pos;
    if ((fp = fopen(pathname, "r")) == NULL)
    {
        Logger::message("Can't open file! Already open? Doesn't exist?");
        return 0;
    }
    if (!fp)
    {
        return 0;
    }
    else
    {
        pathName = pathname;
    }

    /* find the correct datablock's start and end line numbers*/
    while (fgetsl(buf, sizeof buf, fp) != NULL)
    {
        if (strstr(buf, "data_"))
        {
            if (nblock == datablock)
            {
                nstart = nline;
            }
            nblock++;
        }
        nline++;
        if (nblock == datablock+1)
        {
            nend = nline;
        }
        else if (nblock > datablock+1)
        {
            break;
        }
    }
    if (nstart == 0 && nend == 0)
    {
        Logger::message("ERROR:CIFFile: Data block not found");
        return 0;
    }
    Logger::debug("CIF data block lines %d - %d", nstart, nend);
    nline = 0;
    rewind(fp);
    nfound = 0;
    while (fgetsl(buf, sizeof buf, fp) != NULL)
    {
        if (buf[0] == '#')
        {
            nline++;
            continue;
        }
        if (nline < nstart)
        {
            nline++;
            continue;
        }
        if (nline >= nend)
        {
            break;
        }
        nline++;
        if (strstr(buf, "_cell_") || strstr(buf, "_cell."))
        {
            nfound++;
            pos = strstr(buf, "length_a");
            if (pos)
            {
                mapheader->a = (float)atof(pos+8);
            }
            pos = strstr(buf, "length_b");
            if (pos)
            {
                mapheader->b = (float)atof(pos+8);
            }
            pos = strstr(buf, "length_c");
            if (pos)
            {
                mapheader->c = (float)atof(pos+8);
            }
            pos = strstr(buf, "angle_alpha");
            if (pos)
            {
                mapheader->alpha = (float)atof(pos+11);
            }
            pos = strstr(buf, "angle_beta");
            if (pos)
            {
                mapheader->beta = (float)atof(pos+10);
            }
            pos = strstr(buf, "angle_gamma");
            if (pos)
            {
                mapheader->gamma = (float)atof(pos+11);
            }
        }
        if (nfound == 6)
        {
            break;
        }
    }
    ncell = nfound;
    /* find start of symmetry loop */
    rewind(fp);
    nline = 0;
    nfound = 0;
    inloop = false;
    while (fgetsl(buf, sizeof buf, fp) != NULL)
    {
        if (buf[0] == '#')
        {
            nline++;
            continue;
        }
        if (nline < nstart)
        {
            nline++;
            continue;
        }
        if (nline >= nend)
        {
            break;
        }
        if ((pos = strstr(buf, "int_tables_number")) != NULL)
        {
            /* process space group string and
             * look it up in symmop database */
            char *start;
            std::string SymInfoString;
            if (nsymm == 0)
            {
                start = pos + strlen("int_tables_number")+1;
                int spgpno = atoi(start);
                if (spgpno > 0 && spgpno <= 230)
                {
                    mapheader->spgpno = spgpno;
                    mapheader->SetSymmOps();
                    nsymm = mapheader->nsym;
                    goto endsymm;
                }
            }
        }
        if (strstr(buf, "loop_"))
        {
            inloop = true;
        }
        if (buf[0] == '\0')
        {
            inloop = false;
        }
        if ( (strstr(buf, "_symmetry_equiv_pos_as_xyz")
              || strstr(buf, "_symmetry_equiv.pos_as_xyz")) && inloop)
        {
            char opstring[MISymmop::MAXSYMMOPS][MISymmop::MAXSTRING];
            char symmline[4096];
            symmline[0] = '\0';
            nline++;
            while (fgetsl(buf, sizeof buf, fp) != NULL)
            {
                if (buf[0] == '\0' || buf[0] == '_')
                {
                    break;
                }
                if (nline >= nend)
                {
                    break;
                }
                strcat(symmline, symmtrim(buf));
                strcat(symmline, ";");
                nsymm++;
                nline++;
            }
            /* replace last symmops ; with a . */
            symmline[strlen(symmline)-1] = '.';

            // Seems that SymopsString needs to be set here. Fix for problem Paula
            // was having. PWC 6-7-05
            mapheader->nsym = mapheader->scan_symmops(symmline, mapheader->symops, opstring);
            mapheader->SymopsString.clear();
            for (int i = 0; i < mapheader->nsym; i++)
            {
                mapheader->SymopsString.push_back(std::string(opstring[i]));
            }

            // never works!
            //mapheader->SpgpFromSymmops();
            if (nsymm != mapheader->nsym)
            {
                Logger::message("Warning: CIFFile: error in symm ops translation");
            }
            goto endsymm;
        }
        nline++;
    }
endsymm:
    if (nsymm == 0 || ncell == 0)
    {
        char *start /*,  *rest */;
        char SymInfoString[20];
        int ret, i, n;
        /* look for a commented out CRYST PDB record and process it */
        Logger::log("Searching for a CRYST1 card in comment section...");
        rewind(fp);
        while (fgets(buf, sizeof buf, fp) != NULL)
        {
            if (buf[0] == '#')
            {
                if ((start = strstr(buf, "CRYST1")) != NULL)
                {
                    if (ncell < 6)
                    {
                        ret = sscanf(start,
                                     "%*s%f%f%f%f%f%f",
                                     &mapheader->a, &mapheader->b, &mapheader->c,
                                     &mapheader->alpha, &mapheader->beta,
                                     &mapheader->gamma);
                        if (ret != 6)
                        {
                            goto endcryst;
                        }
                        ncell = 6;
                    }
                    /* process space group string and
                     * look it up in symmop database */
                    if (nsymm == 0)
                    {
                        start += 55; /* i.e column 56 */
                        n = 0;
                        for (i = 0; i < 10; i++)
                        {
                            if (isalnum(start[i]))
                            {
                                SymInfoString[n] = tolower(start[i]);
                                n++;
                            }
                        }
                        SymInfoString[n] = '\0';
                        if (n > 1)
                        {
                            if (mapheader->FindSpacegroup(SymInfoString))
                            {
                                nsymm = mapheader->nsym;
                            }
                        }
                        else
                        {
                            goto endcryst;
                        }
                    }
                }
            }
        }
    }
endcryst:

    /* find start of reflection loop */
    rewind(fp);
    int hindex = (-1), kindex = (-1), lindex = (-1), fofoindex = (-1),
        foindex = (-1), fcindex = (-1), phsindex = (-1),
        sigmaindex = (-1), sigmasquareindex = (-1);

    if (!FindCIFIndices(fp, nrefl_start,
                        hindex, kindex, lindex,
                        fofoindex, foindex, fcindex, phsindex,
                        sigmaindex, sigmasquareindex))
    {
        return 0;
    }

    Logger::debug("index hkl %d %d %d fo %d fc %d sigma %d phs %d", hindex, kindex, lindex,
                  foindex, fcindex, sigmaindex, phsindex);
    Logger::debug("index fofo %d sigmasquare %d", fofoindex, sigmasquareindex);
    Logger::debug("Reflections start on line %d", nrefl_start);

    FreeRSet = false;
    FOMsValid = false;
    FcsValid = (fcindex >= 0);
    PhicsValid = (phsindex >= 0);
    FosValid = (foindex >= 0 || fofoindex >= 0);

    /* find max index */
    int maxindex = std::max(hindex, kindex);
    maxindex = std::max(maxindex, lindex);
    maxindex = std::max(maxindex, fofoindex);
    maxindex = std::max(maxindex, foindex);
    maxindex = std::max(maxindex, fcindex);
    maxindex = std::max(maxindex, phsindex);
    if (maxindex < 3)
    {
        Logger::message("CIFFile: Reflection loop missing needed items");
        return 0;
    }
    if (maxindex >= MAXLOOP-1)
    {
        Logger::message("CIFFile: Reflection loop has too many items in row");
        return 0;
    }
    refls.clear();

    nsymm = mapheader->nsym;
    if (mapheader->a != 0 && mapheader->b != 0 && mapheader->c != 0 && mapheader->alpha != 0 && mapheader->beta != 0
        && mapheader->gamma != 0)
    {
        ncell = 6;
    }
    // if still no crystal info found - get from user
    if (ncell < 6  || (phsindex < 0 && nsymm == 0))
    {
        if (!PromptForCrystal())
        {
            Logger::message("CIFFile: Unable to get crystal information");
            return 0;
        }
    }
    // read in the reflections
    rewind(fp);
    orthog(mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, mapheader->ctof);
    uinv(mapheader->ctof, mapheader->ftoc);
    sthol(1, 1, 1, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 1);
    refls_stholmax = -1.0F;
    refls_stholmin = 999999.0F;
    nline = 0;
    CREFL refl;
    inloop = true;
    while (fgetsl(buf, sizeof buf, fp) != NULL)
    {
        if (nline < nrefl_start)
        {
            nline++;
            continue;
        }
        if (buf[0] == '#')
        {
            nline++;
            continue;
        }
        if (nline > nend)
        {
            goto end;
        }
        if (!inloop && (buf[0] == '\0' || buf[0] == '_'))
        {
            goto end;
        }
        if (buf[0] == '\0')
        {
            nline++;
            continue;
        }
        memset(&refl, 0, sizeof(refl));
        int stringsLength = scanrow(buf, strings);
        if (stringsLength > maxindex)
        {
            if (fofoindex >= 0)
            {
                fo = (float)sqrt(atof(strings[fofoindex]));
            }
            else
            {
                fo = (float)atof(strings[foindex]);
            }
            if (fo > 0.00001)
            {
                refl.ind[0] = ih = atoi(strings[hindex]);
                refl.ind[1] = ik = atoi(strings[kindex]);
                refl.ind[2] = il = atoi(strings[lindex]);
                refl.fo = fo;
                if (fcindex >= 0)
                {
                    refl.fc = (float)atof(strings[fcindex]);
                }
                if (sigmaindex >= 0)
                {
                    refl.sigma = (float)atof(strings[sigmaindex]);
                }
                if (sigmasquareindex >= 0)
                {
                    refl.sigma = 0.5F*(float)atof(strings[sigmasquareindex])/fo;
                }
                if (phsindex >= 0)
                {
                    refl.phi = (float)atof(strings[phsindex]);
                }
                refl.sthol = sthol(ih, ik, il, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
                if (refl.sthol > refls_stholmax)
                {
                    refls_stholmax = refl.sthol;
                }
                if (refl.sthol < refls_stholmin)
                {
                    refls_stholmin = refl.sthol;
                }
                //        Logger::debug("ref hkl %d %d %d fo %f fc %f sigma %f phi %f sthol %f", refl.ind[0], refl.ind[1], refl.ind[2],
                //            refl.fo, refl.fc, refl.sigma, refl.phi, refl.sthol);
                refls.push_back(refl);
            }
        }
        else
        {
            goto end;
        }
        nline++;
    }
end:
    fclose(fp);
    mapheader->resmax = 0.5F/refls_stholmin;
    mapheader->resmin = 0.5F/refls_stholmax;
    mapName = pathname;
    pathName = pathname;
    Logger::debug("Read %d reflections", refls.size());
    return refls.size();
}

bool EMapBase::FFTMap(int maptype, int gridlevel,
                      float resMin, float resMax)
{

    // swap if resmin/ resmax reversed
    if (resMin > resMax)
    {
        float t = resMax;
        resMax = resMin;
        resMin = t;
    }

    if (!CanDoMapType(maptype))
    {
        Logger::message("Before you can display the map, you need to use the\n"
                        "Calculate Structure Factors... command to\n"
                        "calculate phases from a model.");
        return false;
    }
    if (resMin != -1.0f)
    {
        mapheader->resmin = (float)std::max((double)resMin, (0.5/refls_stholmax));
    }
    if (resMax != -1.0f)
    {
        mapheader->resmax = (float)std::min((double)resMax, (0.5/refls_stholmin));
    }
    if (maptype >= 0)
    {
        mapheader->maptype = maptype;
    }
    float grid = 2.05F;
    if (gridlevel == 1 || gridlevel < 0)
    {
        grid = 2.9F;
    }
    else if (gridlevel == 2)
    {
        grid = 4.5F;
    }

    if (mapheader->nx == 0 || gridlevel >= 0)
        mapheader->nx = MIMapFactor((int)(grid*mapheader->a/mapheader->resmin), FFT_PRIME, EVEN, 2);
    if (mapheader->ny == 0 || gridlevel >= 0)
        mapheader->ny = MIMapFactor((int)(grid*mapheader->b/mapheader->resmin), FFT_PRIME, ODDOREVEN, 1);
    if (mapheader->nz == 0 || gridlevel >= 0)
        mapheader->nz = MIMapFactor((int)(grid*mapheader->c/mapheader->resmin), FFT_PRIME, ODDOREVEN, 1);
    mapheader->hmax = (int)(mapheader->a/mapheader->resmin);
    mapheader->kmax = (int)(mapheader->b/mapheader->resmin);
    mapheader->lmax = (int)(mapheader->c/mapheader->resmin);
    /* code in fftsub.c:  if ((nx <= 2*hmax) || (ny <= 2*kmax) || (nz <= 2*lmax)) */
    while (mapheader->nx <= 2*mapheader->hmax)
    {
        mapheader->nx = MIMapFactor(mapheader->nx+1, FFT_PRIME, EVEN, 2);
    }
    while (mapheader->ny <= 2*mapheader->kmax)
    {
        mapheader->ny = MIMapFactor(mapheader->ny+1, FFT_PRIME, ODDOREVEN, 1);
    }
    while (mapheader->nz <= 2*mapheader->lmax)
    {
        mapheader->nz = MIMapFactor(mapheader->nz+1, FFT_PRIME, ODDOREVEN, 1);
    }
    Logger::log("Resolution range:%0.2f - %0.2f\nMap grid %d %d %d\nUnit cell: %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f",
                mapheader->resmax, mapheader->resmin, mapheader->nx, mapheader->ny, mapheader->nz, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma);
    return FFTCalc();
}

bool EMapBase::FFTCalc()
{
    mapheader->hmax = (int)(mapheader->a/mapheader->resmin);
    mapheader->kmax = (int)(mapheader->b/mapheader->resmin);
    mapheader->lmax = (int)(mapheader->c/mapheader->resmin);
    /* prevent memory leaks - be sure to free emap  */
    float *emap;
    try
    {
        emap = fft3d(&refls[0], refls.size(), mapheader, 0);
    }
    catch (...)
    {
        Logger::message("Program Error in FFT'ing map");
        return false;
    }
    if (emap == NULL)
    {
        mapheader->nx = mapheader->ny = mapheader->nz = 0;
        return false;
    }
    try
    {
        map_points.resize(mapheader->nx*mapheader->ny*mapheader->nz);
        memcpy(&map_points[0], emap, mapheader->nx*mapheader->ny*mapheader->nz*sizeof(float));
    }
    catch (...)
    {
        map_points.clear();
        Logger::message("Program Error: Unable to allocate memory for map");
        return false;
    }
    free(emap);
    ScaleMap(CalcRMS());
    return true;
}

float EMapBase::avgrho(float fx, float fy, float fz)
{
    float dx0, dy0, dz0;
    int nx, ny, nz, x0, x1, y0, y1, z0, z1;

    nx = mapheader->nx;
    ny = mapheader->ny;
    nz = mapheader->nz;
    if (fx < 0.0)
    {
        fx += ((int)(-fx)) + 1.0F;
    }
    if (fy < 0.0)
    {
        fy += ((int)(-fy)) + 1.0F;
    }
    if (fz < 0.0)
    {
        fz += ((int)(-fz)) + 1.0F;
    }
    fx *= (float)nx;
    fy *= (float)ny;
    fz *= (float)nz;
    x0 = (int)fx;
    y0 = (int)fy;
    z0 = (int)fz;
    x1 = x0 +1;
    y1 = y0 +1;
    z1 = z0 +1;
    dx0 = fx - (float)x0;
    dy0 = fy - (float)y0;
    dz0 = fz - (float)z0;
    /* make sure indices are in bounds */
    x0 = (x0)%nx;
    y0 = (y0)%ny;
    z0 = (z0)%nz;
    x1 = (x1)%nx;
    y1 = (y1)%ny;
    z1 = (z1)%nz;
    return (pseudospline(map_points[nx*(ny*z0 + y0) + x0],
                         map_points[nx*(ny*z1 + y0) + x0],
                         map_points[nx*(ny*z0 + y1) + x0],
                         map_points[nx*(ny*z1 + y1) + x0],
                         map_points[nx*(ny*z0 + y0) + x1],
                         map_points[nx*(ny*z1 + y0) + x1],
                         map_points[nx*(ny*z0 + y1) + x1],
                         map_points[nx*(ny*z1 + y1) + x1],
                         dx0, dy0, dz0));
}

std::string EMapBase::MapID()
{
    std::string ID = mapName;
    if (HasPhases())
    {
        char buf[1024];
        sprintf(buf, " %s %.3f - %.3f", maptypes[mapheader->maptype],
                mapheader->resmax, mapheader->resmin);
        ID += std::string(buf);
    }
    return ID;
}

long EMapBase::LoadCNSMap(const char *pathname)
{
    FILE *fp = fopen(pathname, "r");
    if (!fp)
    {
        return 0;
    }
    Logger::log("Loading CNS format...");
    char buf[4096];
    std::string ubuf, mess;
    long nmap = 0;
    int im;
    int NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX;
    int ix, iy, iz;

    map_points.clear();
    // title (ignored)
    if (!next_line(fp, buf, sizeof(buf), ubuf))
    {
        fclose(fp);
        return 0;
    }

    //bounds and size
    if (!next_line(fp, buf, sizeof(buf), ubuf))
    {
        fclose(fp);
        return 0;
    }
    sscanf(buf, "%d%d%d%d%d%d%d%d%d", &NA, &AMIN, &AMAX, &NB, &BMIN, &BMAX, &NC, &CMIN, &CMAX);
    Logger::log("Map Size: %d %d %d\nMap Start:%d %d %d\nMap End: %d %d %d",
                NA, NB, NC, AMIN, BMIN, CMIN, AMAX, BMAX, CMAX);
    //bounds and size
    if (!next_line(fp, buf, sizeof(buf), ubuf))
    {
        fclose(fp);
        return 0;
    }
    sscanf(buf, "%f%f%f%f%f%f",
           &mapheader->a,
           &mapheader->b,
           &mapheader->c,
           &mapheader->alpha,
           &mapheader->beta,
           &mapheader->gamma);
    Logger::log("Unit cell: %f %f %f %f %f %f",
                mapheader->a,
                mapheader->b,
                mapheader->c,
                mapheader->alpha,
                mapheader->beta,
                mapheader->gamma);
    //sort order
    if (!next_line(fp, buf, sizeof(buf), ubuf))
    {
        fclose(fp);
        return 0;
    }

    nmap = (long)NA*(long)NB*(long)NC;
    mapheader->nx = NA;
    mapheader->ny = NB;
    mapheader->nz = NC;
    try
    {
        map_points.resize(nmap);
    }
    catch (...)
    {
        map_points.clear();
        Logger::log("Bummer: We ran out of memory in LoadCNSMap - aborting");
        fclose(fp);
        return 0;
    }

    int nsection = ((BMAX-BMIN)+1)*((AMAX-AMIN)+1);
    int i, ksect, ixt, iyt, izt, isym;
    float a, x, y, z, xp, yp, zp;

    for (iz = CMIN; iz <= CMAX; iz++)
    {
        if (fscanf(fp, "%d", &ksect) != 1)
        {
            Logger::log("Bummer: I/O error in LoadCNSMap - aborting");
            map_points.clear();
            fclose(fp);
            return 0;
        }
        ix = AMIN;
        iy = BMIN;
        for (i = 0; i < nsection; i++)
        {
            if (fscanf(fp, "%f", &a) != 1)
            {
                Logger::log("Bummer: I/O error in LoadCNSMap - aborting");
                map_points.clear();
                fclose(fp);
                return 0;
            }
            for (isym = 0; isym < mapheader->nsym; isym++)
            {
                x = (float)ix/(float)mapheader->nx;
                y = (float)iy/(float)mapheader->ny;
                z = (float)iz/(float)mapheader->nz;
                mapheader->symm_mh(x, y, z, &xp, &yp, &zp, isym);
                ixt = ROUND(xp*(float)mapheader->nx);
                iyt = ROUND(yp*(float)mapheader->ny);
                izt = ROUND(zp*(float)mapheader->nz);
                im = mapdex(ixt, iyt, izt);
                if (im < 0 || im > nmap)
                {
                    Logger::debug("Error: out of range %d; (%d %d %d)", im, ixt, iyt, izt);
                }
                else
                {
                    map_points[im] = a;
                }
            }
            ix++;
            if (ix > AMAX)
            {
                iy++;
                ix = AMIN;
            }
        }
    }
    float rave, rsigma;
    if (fscanf(fp, "%d", &ksect) != 1)
    {
        Logger::log("Bummer: I/O error in LoadCNSMap - aborting");
        map_points.clear();
        fclose(fp);
        return 0;
    }
    if (fscanf(fp, "%f%f", &rave, &rsigma) != 2)
    {
        Logger::log("Bummer: I/O error in LoadCNSMap - aborting");
        map_points.clear();
        fclose(fp);
        return 0;
    }
    for (unsigned int im = 0; im < map_points.size(); im++)
    {
        map_points[im] = map_points[im]*rsigma + rave;
    }

    mapName = pathname;
    pathName = pathname;
    fclose(fp);
    return nmap;
}

long EMapBase::LoadCCP4Map(const char *pathname, float *rms, float *min, float *max)
{
    FILE *fp = fopen(pathname, "rb");
    // controls swabbing - change to true for other endian machine ---
    char buf[4096];
    bool swab = true;
    unsigned int nc, nr, ns;
    int mode;
    int ncstart, nrstart, nsstart;
    int nx, ny, nz, ix = 0, iy = 0, iz = 0, isym, ixt, iyt, izt;
    int mapc, mapr, maps, ispg, nsymbt, lskflag;
    float amin, amax, amean;
    bool file_is_map = false;
    //std::vector<float> vmap;
    float vmap;
    double sum;
    unsigned int nmap, ic, ir, is;
    //int indx, indy, indz;
    int msize[3];
    float x, y, z, xp, yp, zp;

    if (!fp)
    {
        return 0;
    }
    else
    {
        pathName = pathname;
    }
    nc = 0;
    rewind(fp);
    fread(buf, 1, 4*52, fp);

    if (buf[0] == 0 && buf[4] == 0 && buf[8] == 0)
    {
        swab = true;
    }
    else
    {
        swab = false;
    }
    fread(buf, 1, 4, fp);
    if (buf[0] == 'M' && buf[1] == 'A' && buf[2] == 'P')
    {
        file_is_map = true;
    }

    if (!file_is_map)
    {
        Logger::log("File is not a CCP4 map file");
        fclose(fp);
        return 0;
    }
    if (swab)
    {
        Logger::log("unsigned char order will be swabbed");
    }

    rewind(fp);

    msize[0] = nc = map_in_int(fp, swab);
    msize[1] = nr = map_in_int(fp, swab);
    msize[2] = ns = map_in_int(fp, swab);
    Logger::log("Map size = %d %d %d", nc, nr, ns);

    mode = map_in_int(fp, swab);
    if (!(mode == 2 || mode == 0))
    {
        Logger::log("Unsupported map mode %d - send email to ccms-help@sdsc.edu with mode number and brief explanation\n", mode);
        fclose(fp);
        return 0;
    }
    Logger::log("Map mode = %d", mode);

    if (mode == 0)
    {
        Logger::log("Mask mode not implemented yet...");
        fclose(fp);
        return 0;
    }

    ncstart = -1 * map_in_int(fp, swab);
    nrstart = -1 * map_in_int(fp, swab);
    nsstart = -1 * map_in_int(fp, swab);
    nx = map_in_int(fp, swab);
    ny = map_in_int(fp, swab);
    nz = map_in_int(fp, swab);
    Logger::log("Start = %d %d %d; Extent = %d %d %d", ncstart, nrstart, nsstart, nx, ny, nz);
    mapheader->a = map_in_float(fp, swab);
    mapheader->b = map_in_float(fp, swab);
    mapheader->c = map_in_float(fp, swab);
    mapheader->alpha = map_in_float(fp, swab);
    mapheader->beta = map_in_float(fp, swab);
    mapheader->gamma = map_in_float(fp, swab);
    Logger::log("Unit cell: %f %f %f %f %f %f",
                mapheader->a,
                mapheader->b,
                mapheader->c,
                mapheader->alpha,
                mapheader->beta,
                mapheader->gamma);
    mapc = map_in_int(fp, swab);
    /*if(mapc == 1) indx = 0;
       if(mapc == 2) indy = 0;
       if(mapc == 3) indz = 0;
     */
    mapr = map_in_int(fp, swab);
    /*if(mapr == 1) indx = 1;
       if(mapr == 2) indy = 1;
       if(mapr == 3) indz = 1;
     */
    maps = map_in_int(fp, swab);
    /*if(maps == 1) indx = 2;
       if(maps == 2) indy = 2;
       if(maps == 3) indz = 2;
     */
    amin = map_in_float(fp, swab);
    *min = amin;
    amax = map_in_float(fp, swab);
    *max = amax;
    amean = map_in_float(fp, swab);
    ispg = map_in_int(fp, swab);
    mapheader->spgpno = ispg;
    mapheader->SetSymmOps();
    nsymbt = map_in_int(fp, swab);
    lskflag = map_in_int(fp, swab);
    if (lskflag)
    {
        Logger::log("Error: Map is skew map - sorry I can't read skew maps!");
        fclose(fp);
        return 0;
    }
    // Skip the skew matrices
    fread(buf, 1, 12*4, fp);
    // Skip reserved for future use
    fread(buf, 1, 15*4, fp);
    // Skip map type and machine type
    fread(buf, 1, 2*4, fp);

    float arms = map_in_float(fp, swab);
    *rms = arms;
    Logger::log("Density min = %f, max = %f", amin, amax);
    Logger::log("Density mean = %f, rms = %f", amean, arms);

    // rewind and skip the header
    rewind(fp);
    fread(buf, 1, 256*4, fp);

    // skip the symmetry card unsigned chars
    //fread(buf, 1, 80*4, fp);
    fread(buf, 1, nsymbt, fp);

    mapheader->nx = nx;
    mapheader->ny = ny;
    mapheader->nz = nz;

    nmap = (unsigned int)nx*(unsigned int)ny*(unsigned int)nz;
    try
    {
        map_points.resize(nmap);
        //vmap.resize(nc);
    }
    catch (...)
    {
        Logger::log("Error: Out of memory in LoadCCP4Map - aborting");
        fclose(fp);
        return 0;
    }
    if ( /*vmap.capacity() < nc ||*/ map_points.capacity() < nmap)
    {
        Logger::log("Error: Out of memory in LoadCCP4Map - aborting");
        fclose(fp);
        return 0;
    }

    // scarf up the map
    sum = 0.0;
    int im;
    for (is = 0; is < ns; is++)
    {
        if (maps == 1)
        {
            ix = is - nsstart;
        }
        else if (maps == 2)
        {
            iy = is - nsstart;
        }
        else
        {
            iz = is - nsstart;
        }
        for (ir = 0; ir < nr; ir++)
        {
            if (mapr == 1)
            {
                ix = ir - nrstart;
            }
            else if (mapr == 2)
            {
                iy = ir - nrstart;
            }
            else
            {
                iz = ir - nrstart;
            }
            for (ic = 0; ic < nc; ic++)
            {
                if (mapc == 1)
                {
                    ix = ic - ncstart;
                }
                else if (mapc == 2)
                {
                    iy = ic - ncstart;
                }
                else
                {
                    iz = ic - ncstart;
                }
                vmap = map_in_float(fp, swab);
                sum += vmap*vmap;
                // sort map into map_points
                for (isym = 0; isym < mapheader->nsym; isym++)
                {
                    x = (float)ix/(float)nx;
                    y = (float)iy/(float)ny;
                    z = (float)iz/(float)nz;
                    mapheader->symm_mh(x, y, z, &xp, &yp, &zp, isym);
                    ixt = ROUND(xp*(float)nx);
                    iyt = ROUND(yp*(float)ny);
                    izt = ROUND(zp*(float)nz);
                    im = mapdex(ixt, iyt, izt);
                    if (im < 0 || im > (int)nmap)
                    {
                        Logger::debug("Error: out of range %d; (%d %d %d)", im, ixt, iyt, izt);
                    }
                    else
                    {
                        map_points[im] = vmap;
                    }
                }
            }
        }
    }

    mapName = pathname;
    pathName = pathname;
    fclose(fp);
    return nmap;
}

bool EMapBase::LoadMapFile(const char *s)
{
    if (IsCCP4MapFile(s))
    {
        float rms;
        float min;
        float max;
        if (LoadCCP4Map(s, &rms, &min, &max) > 0)
        {
            mapmin = min;
            mapmax = max;
            // Although file's rms is used rather than calculated rms, it is calculated
            // here to set the predicted map type.
            CalcRMS();
            ScaleMap(rms);
            return true;
        }
        else
        {
            return false;
        }
    }
    else if (IsCNSMapFile(s))
    {
        if (LoadCNSMap(s) > 0)
        {
            ScaleMap(CalcRMS());
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        mapheader->SetSymmOps();
        if (LoadFSFOURMapFile(s) > 0)
        {
            ScaleMap(CalcRMS());
            return true;
        }
        else
        {
            return false;
        }
    }
    return false;
}

bool EMapBase::IsCCP4MapFile(const char *pathname)
{
    FILE *fp = fopen(pathname, "rb");
    if (!fp)
    {
        return false;
    }
    char buf[512];
    bool file_is_map = false;
    fread(buf, 1, 4*52, fp);
    fread(buf, 1, 4, fp);
    if (buf[0] == 'M' && buf[1] == 'A' && buf[2] == 'P')
    {
        file_is_map = true;
    }
    else if (buf[3] == 'M' && buf[2] == 'A' && buf[1] == 'P')
    {
        file_is_map = true;
    }
    fclose(fp);
    return (file_is_map);
}

bool EMapBase::IsCNSMapFile(const char *pathname)
{
    FILE *fp = fopen(pathname, "r");
    if (!fp)
    {
        return false;
    }
    char buf[512];
    std::string ubuf, last;
    bool cell = false;
    bool mode = false;
    int nlines = 0;
    while (fgets(buf, sizeof(buf), fp) != NULL && nlines < 10)
    {
        ubuf = buf;
        ubuf = MIToUpper(ubuf);
        MIStringTrim(ubuf, false); // trims space from left
        if (strncmp(ubuf.c_str(), "REMARK", 6) == 0 || ubuf.length() == 0)
        {
            continue;
        }
        if (strncmp(ubuf.c_str(), "ZYX", 3) == 0)
        {
            mode = true;
            // check to see if last line makes sense as a cell
            float a, b, c, alpha, beta, gamma;
            if (sscanf(last.c_str(), "%f%f%f%f%f%f", &a, &b, &c, &alpha, &beta, &gamma) == 6)
            {
                cell = a > 0 && b > 0 && c > 0 && alpha >= 30.0F && beta >= 30.0F && gamma >= 30.0F
                       && alpha < 180.0F && beta < 180.0F && gamma < 180.0F;
                break;
            }
        }
        last = ubuf;
        nlines++;
    }
    fclose(fp);
    return (mode && cell && nlines < 10);
}

long EMapBase::LoadFSFOURMapFile(const char *pathname)
{
    FILE *fp = fopen(pathname, "rb");
    if (!fp)
    {
        return 0;
    }
    else
    {
        pathName = pathname;
    }

    unsigned int nmap = 0;
    std::vector<int> rho;
    int swab;
    float rhoscale, rmsrho;
    char buf[1024];
    nmap = fssize_uni(fp, &swab);
    if (nmap <= 3)
    {
        Logger::log("Error:: Nonsensical map size - probable bad file, wrong file type or wrong unsigned char order - aborting");
        fclose(fp);
        return false;
    }
    Logger::log("Size of map %d", (int)nmap);
    rewind(fp);

    rho.resize(nmap);
    if (rho.capacity() < nmap)
    {
        Logger::log("Abort: Out of memory in Emap::LoadFSFOURMapFile");
        return 0;
    }
    fsread_uni(rho, &mapheader->nx, &mapheader->ny, &mapheader->nz, &rhoscale, &rmsrho, fp, buf, swab);
    if (rmsrho <= 0.00001F)
    {
        rmsrho = 1.0;
    }
    for (unsigned int i = 0; i < nmap; i++)
    {
        map_points.push_back((float)(rho[i]));
    }
    if (map_points.capacity() < nmap)
    {
        Logger::log("Abort: Out of memory in Emap::LoadFSFOURMapFile");
        return 0;
    }
    ScaleMap(CalcRMS());
    mapName = pathname;
    pathName = pathname;
    fclose(fp);
    return nmap;
}

bool EMapBase::SavePhases(const char *pathname, int type)
{
    FILE *fp = fopen(pathname, "w");
    bool result = false;
    if (fp)
    {
        if (type == EMapBase::mmCIF_phase)
        {
            result = SavemmCIF(fp) > 0;
        }
        else if (type == EMapBase::XtalView_phase)
        {
            result = SaveXtalViewPhase(fp) > 0;
        }
        else if (type == EMapBase::CCP4_phase)
        {
            result = SaveCCP4Phase(pathname) > 0;
        }
        else if (type == EMapBase::CNS_phase)
        {
            result = SaveCNSPhase(fp) > 0;
        }
        else if (type == EMapBase::Warp_phase)
        {
            result = SaveWarpPhase(fp) > 0;
        }
        fclose(fp);
    }
    return result;
}

long EMapBase::SavemmCIF(FILE *fp)
{
    int i;
    fprintf(fp, "#\n#  mmCIF phase file written by Mi-fit\n#\ndata_mifit_phase\n\n");
    fprintf(fp, "loop_\n _symmetry_equiv_pos_as_xyz\n");
    for (i = 0; i < mapheader->nsym; i++)
    {
        fprintf(fp, "\'%s\'\n", mapheader->SymopsString[i].c_str());
    }
    fprintf(fp, "\n_cell_length_a %10.4f\n_cell_length_b %10.4f\n_cell_length_c %10.4f\n_cell_angle_alpha %10.4f\n_cell_angle_beta %10.4f\n_cell_angle_gamma  %10.4f\n",
            mapheader->a,
            mapheader->b,
            mapheader->c,
            mapheader->alpha,
            mapheader->beta,
            mapheader->gamma);
    fprintf(fp, "\nloop_\n"
            " _refln.index_h\n"
            " _refln.index_k\n"
            " _refln.index_l\n"
            " _refln.F_meas\n"
            " _refln.F_sigma\n"
            " _refln.F_calc\n"
            " _refln.phase_calc\n");
    for (i = 0; (unsigned int)i < refls.size(); i++)
    {
        fprintf(fp, "%d %d %d %f %f %f %f\n",
                refls[i].ind[0], refls[i].ind[1], refls[i].ind[2],
                refls[i].fo, refls[i].sigma, refls[i].fc, refls[i].phi);
    }
    return refls.size();
}

long EMapBase::SaveXtalViewPhase(FILE *fp)
{
    unsigned int i, n = refls.size();
    for (i = 0; i < n; i++)
    {
        fprintf(fp, "%5d %4d %4d %8.2f %8.2f %8.2f\n",
                refls[i].ind[0], refls[i].ind[1], refls[i].ind[2], refls[i].fo, refls[i].fc, refls[i].phi);
    }
    return n;
}

long EMapBase::SaveWarpPhase(FILE *fp)
{
    unsigned int i, n = refls.size();
    for (i = 0; i < n; i++)
    {
        fprintf(fp, "%5d %4d %4d %8.2f %8.3f %8.2f %7.3f\n",
                refls[i].ind[0], refls[i].ind[1], refls[i].ind[2], refls[i].fo, refls[i].sigma, refls[i].coef, refls[i].phi);
    }
    return n;
}

long EMapBase::SaveCNSPhase(FILE*)
{
    //todo
    Logger::log("Not implemented yet");
    return 0;
}

void EMapBase::SetCrystal(const char *crystal)
{
    mapheader->LoadCrystal(crystal);
    Crystal = crystal;
}

bool EMapBase::IsCCP4MTZFile(const char *pathname)
{
    FILE *fp = fopen(pathname, "rb");
    if (!fp)
    {
        return false;
    }
    char buf[512];
    bool file_is_mtz = false;

    int n = fread(buf, 1, 4, fp);
    if (n != 4)
    {
        file_is_mtz = false;
    }
    else if (buf[0] == 'M' && buf[1] == 'T' && buf[2] == 'Z')
    {
        file_is_mtz = true;
    }
    fclose(fp);
    return (file_is_mtz);
}

bool EMapBase::SFCalc(Residue *res)
{
    float fft_scale;
    bool result = false;
    if (Monomer::isValid(res))
    {
        result = sfFFT(res, fft_scale) != 0;
    }
    return result;
}

int EMapBase::sfFFT(Residue *res, float &scale)
{
    micomplex *rho;
    float sc;
    int i, nx, ny, nz;
    CMapHeaderBase *mh = mapheader;

    /* save nx ny nz  res limits */
    ny = mapheader->ny;
    nx = mapheader->nx;
    nz = mapheader->nz;
    Logger::log("Building rho...");
    if ((rho = buildrho(res, mh)) == NULL)
    {
        return (0);
    }
    Logger::log("Inverting rho...");
    InvertMap(rho, mh, refls);
    free((void*)rho);
    Logger::log("Scaling Fc...");
    sc = ComputeScale(refls, mh);
    ApplyScale(refls, sc, mh);
    Logger::log("Fc Scale = %0.5f", sc);
    Logger::log("Phase Change = %0.3f degrees\n", RePhase(refls, mh));

    // these fields are set by RePhase
    PhicsValid = true;
    FcsValid = true;

    float fosum = 0.0, difsum = 0.0;
    for (i = 0; (unsigned int)i < refls.size(); i++)
    {
        refls[i].awhole = refls[i].acalc;
        refls[i].bwhole = refls[i].bcalc;
        difsum += (float)fabs(refls[i].fc - refls[i].fo);
        fosum += refls[i].fo;
    }
    RFactor = difsum/fosum;
    Logger::log("R-factor = %0.3f", RFactor);
    scale = sc;
    /* unsave nx ny nz */
    mapheader->ny = ny;
    mapheader->nx = nx;
    mapheader->nz = nz;
    return (1);
}

bool parseCoefficients(const char *str, int &mapType)
{
    const char *coeff = str;
    // Set default value
    mapType = MIMapType::DirectFFT;

    if (coeff == NULL)
    {
        return false;
    }

    // Find first character after first space
    while (*coeff != '\0' && !isspace(*coeff))
    {
        ++coeff;
    }
    while (*coeff != '\0' && isspace(*coeff))
    {
        ++coeff;
    }
    // Find last non-whitespace character
    const char *coeffEnd = str + strlen(str);
    while ((coeffEnd-1) > coeff && isspace(*(coeffEnd-1)))
    {
        --coeffEnd;
    }

    // Match map type string
    size_t len = coeffEnd - coeff;
    if (len > 0)
    {
        for (unsigned int i = 0; i < MAP_TYPE_COUNT; i++)
        {
            if (strncasecmp(coeff, maptypes[i], len) == 0 && len == strlen(maptypes[i]))
            {
                mapType = i;
                break;
            }
        }
    }
    return true;
}

bool EMapBase::ScriptCommand(const char *ibuf)
{
    static int c = 21; // Colors::MAP1;
    int d;
    float f;
    const char *ubuf;
    const char *buf;
    std::string s(ibuf);
    ubuf = s.c_str();
    std::string s2 = s;
    s2 = MIToLower(s2);
    buf = s2.c_str();
    if (strncmp(buf, "color", 5) == 0)
    {
        sscanf(buf, "%*s%d", &c);
        return true;
    }
    else if (strncmp(buf, "contourmap", 10) == 0)
    {
        if (sscanf(buf, "%*s%d", &d) == 1)
        {
            if (d == mapnumber)
            {
                changed = true;
            }
        }
        return true;
        //reading_commands = false;
    }
    else if (strncmp(buf, "maplinewidth", 12) == 0)
    {
        if (sscanf(buf, "%*s%f", &f) == 1)
        {
            settings->maplinewidth = f;
        }
        return true;
    }
    else if (strncmp(buf, "contourcolor", 12) == 0)
    {
        if (sscanf(buf, "%*s%d", &d) == 1)
        {
            settings->MapLevelColor[d-1] = c;
        }
        return true;
    }
    else if (strncmp(buf, "contourleveldefault", 19) == 0)
    {
        sscanf(buf, "%*s%f%f%f%f%f", &settings->MapLevel[0], &settings->MapLevel[1], &settings->MapLevel[2], &settings->MapLevel[3], &settings->MapLevel[4]);
        return true;
    }
    else if (strncmp(buf, "contourlevels", 12) == 0)
    {
        int d1, d2, d3, d4, d5;
        if (sscanf(buf, "%*s%d%d%d%d%d", &d1, &d2, &d3, &d4, &d5) == 5)
        {
            settings->MapLevelOn[0] = d1 != 0;
            settings->MapLevelOn[1] = d2 != 0;
            settings->MapLevelOn[2] = d3 != 0;
            settings->MapLevelOn[3] = d4 != 0;
            settings->MapLevelOn[4] = d5 != 0;
            return true;
        }
        else if (sscanf(buf, "%*s%d", &d) == 1)
        {
            settings->MapLevelOn[0] = (d&1) != 0;
            settings->MapLevelOn[1] = (d&2) != 0;
            settings->MapLevelOn[2] = (d&4) != 0;
            settings->MapLevelOn[3] = (d&8) != 0;
            settings->MapLevelOn[4] = (d&16) != 0;
            /* no make sense to do this me think - dem
               for(i=0; i<5; i++)
               if(MapLevelOn[i] > 0)
                MapLevel[i]=1;
             */
            return true;
        }
    }
    else if (strncmp(buf, "contourradius", 13) == 0)
    {
        if (sscanf(buf, "%*s%f", &f) == 1)
        {
            settings->Radius = f;
        }
        return true;
    }
    else if (strncmp(buf, "fftapply", 8) == 0)
    {
        return FFTMap(mapheader->maptype, -2, mapheader->resmin, mapheader->resmax) == true;
    }
    else if (strncmp(buf, "coefficient", 6) == 0)
    {
        if (parseCoefficients(buf, mapheader->maptype))
        {
            return true;
        }
    }
    else if (strncmp(buf, "resmin", 6) == 0)
    {
        if (sscanf(buf, "%*s%f", &f) == 1)
        {
            mapheader->resmin = f;
        }
        return true;
    }
    else if (strncmp(buf, "resmax", 6) == 0)
    {
        if (sscanf(buf, "%*s%f", &f) == 1)
        {
            mapheader->resmax = f;
        }
        return true;
    }
    else if (strncmp(buf, "unitcell", 6) == 0)
    {
        sscanf(buf, "%*s%f%f%f%f%f%f", &mapheader->a, &mapheader->b, &mapheader->c, &mapheader->alpha, &mapheader->beta, &mapheader->gamma);
        return true;
    }
    else if (strncmp(buf, "spacegroupno", 6) == 0)
    {
        if (sscanf(buf, "%*s%d", &d) == 1)
        {
            mapheader->spgpno = d;
            mapheader->SetSymmOps();
            return true;
        }
    }
    else if (strncmp(buf, "crystal", 7) == 0)
    {
        char name[2000];
        if (sscanf(ubuf, "%*s%s", name) == 1)
        {
            mapheader->LoadCrystal(name);
            return true;
        }
    }
    else if (strncmp(buf, "ctitle", 6) == 0 || strncmp(buf, "name", 4) == 0)
    {
        mapheader->title = ubuf;
        MIStringTrim(mapheader->title);
        if (strncmp(mapheader->title.c_str(), "ctitle ", 7) == 0)
        {
            mapheader->title = std::string(&mapheader->title[7]);
        }
        if (strncmp(mapheader->title.c_str(), "name ", 5) == 0)
        {
            mapheader->title = std::string(&mapheader->title[5]);
        }
        return true;
    }
    else if (strncmp(buf, "fftnx", 5) == 0)
    {
        if (sscanf(buf, "%*s%d", &d) == 1)
        {
            mapheader->nx = MIMapFactor(d, FFT_PRIME, EVEN, 2);
        }
        return true;
    }
    else if (strncmp(buf, "fftny", 5) == 0)
    {
        if (sscanf(buf, "%*s%d", &d) == 1)
        {
            mapheader->ny = MIMapFactor(d, FFT_PRIME, ODDOREVEN, 1);
        }
        return true;
    }
    else if (strncmp(buf, "fftnz", 5) == 0)
    {
        if (sscanf(buf, "%*s%d", &d) == 1)
        {
            mapheader->nz = MIMapFactor(d, FFT_PRIME, ODDOREVEN, 1);
        }
        return true;
    }
    return true;
}

int EMapBase::Read(const char *file)
{
    int n = 0;
    char ubuf[1000];
    const char *buf;
    bool reading_commands = false;
    //int c= MAP1;
    //int i;
    int d;
    //float f;
    FILE *fp = fopen(file, "r");
    if (!fp)
    {
        return 0;
    }

    while (fgets(ubuf, sizeof(ubuf), fp) != NULL)
    {
        std::string s(ubuf);
        s = MIToLower(s);
        buf = s.c_str();
        if (strncmp(buf, "maptocont", 9) == 0 || strncmp(buf, "loadmap", 7) == 0)
        {
            if (sscanf(buf, "%*s%d", &d) == 1)
            {
                if (d == mapnumber)
                {
                    reading_commands = true;
                }
                else
                {
                    reading_commands = false;
                }
            }
        }
        else if (reading_commands)
        {
            reading_commands = ScriptCommand(ubuf);
        }
    }
    fclose(fp);
    return n;
}

//DEL int EMapBase::ReadCrystal(const char *file)
//DEL {
//DEL /*  moved into ScriptCommand for greater code consistency
//DEL   int n=0;
//DEL   int d;
//DEL   float f;
//DEL   char ubuf[1000];
//DEL   char buf[1000];
//DEL   bool reading_commands=false;
//DEL   //int c= MAP1;
//DEL   int i;
//DEL
//DEL
//DEL   FILE * fp = fopen(file,"r");
//DEL   if(!fp) return 0;
//DEL
//DEL   while(fgets(ubuf, sizeof(ubuf), fp)!=NULL){
//DEL     std::string s(ubuf);
//DEL     s.LowerCase();
//DEL     strcpy(buf, s.c_str());
//DEL     if(strncmp(buf,"maptocont",9)==0 || strncmp(buf,"loadmap",7)==0){
//DEL       if(sscanf(buf,"%*s%d", &d)==1)
//DEL         if(d == mapnumber) reading_commands = true;
//DEL         else reading_commands = false;
//DEL     }
//DEL     if(reading_commands){
//DEL       if(strncmp(buf,"coefficients",6)==0){
//DEL         char coeff[200];
//DEL         if(sscanf(ubuf,"%*s%s", coeff)==1)
//DEL           for(i=0; i< MAP_TYPE_COUNT; i++)
//DEL             if(strncasecmp(coeff,maptypes[i],strlen(coeff))==0 && strlen(coeff)==strlen(maptypes[i]))
//DEL               mapheader->maptype = i;
//DEL       }else if(strncmp(buf,"resmin",6)==0){
//DEL         if(sscanf(buf,"%*s%f", &f)==1) mapheader->resmin = f;
//DEL       }else if(strncmp(buf,"resmax",6)==0){
//DEL         if(sscanf(buf,"%*s%f", &f)==1) mapheader->resmax = f;
//DEL       }else if(strncmp(buf,"unitcell",6)==0){
//DEL         sscanf(buf,"%*s%f%f%f%f%f%f", &mapheader->a,&mapheader->b,&mapheader->c,&mapheader->alpha,&mapheader->beta,&mapheader->gamma);
//DEL       }else if(strncmp(buf,"spacegroupno",6)==0){
//DEL         if(sscanf(buf,"%*s%d", &d)==1) {
//DEL           mapheader->spgpno = d;
//DEL           mapheader->SetSymmOps();
//DEL         }
//DEL       }else if(strncmp(buf,"crystal",7)==0){
//DEL         char name[2000];
//DEL         if(sscanf(ubuf,"%*s%s", name)==1) {
//DEL           mapheader->LoadCrystal(name);
//DEL         }
//DEL       }else if(strncmp(buf,"name",4)==0){
//DEL         mapheader->title = ubuf;
//DEL         mapheader->title.Replace("\n", "");
//DEL       }else if(strncmp(buf,"fftnx",5)==0){
//DEL         if(sscanf(buf,"%*s%d", &d)==1) mapheader->nx = MIMapFactor(d, FFT_PRIME, EVEN, 2);
//DEL       }else if(strncmp(buf,"fftny",5)==0){
//DEL         if(sscanf(buf,"%*s%d", &d)==1) mapheader->ny = MIMapFactor(d, FFT_PRIME, ODDOREVEN, 1);
//DEL       }else if(strncmp(buf,"fftnz",5)==0){
//DEL         if(sscanf(buf,"%*s%d", &d)==1) mapheader->nz = MIMapFactor(d, FFT_PRIME, ODDOREVEN, 1);
//DEL       }
//DEL     }
//DEL   }
//DEL   fclose(fp);
//DEL   return n;
//DEL
//DEL     return 0;
//DEL     */
//DEL   return 0;
//DEL }


void EMapBase::CheckCenter(float x, float y, float z)
{
    if (!visible)
    {
        return;
    }
    // check new center against old - if near the edge of the
    // map, recontour the map.

    float center[3] = {x, y, z};
    double dx = center[0] - map_center[0];
    double dy = center[1] - map_center[1];
    double dz = center[2] - map_center[2];
    double d = (dx*dx + dy*dy + dz*dz);
    float r = (0.05*settings->Radius);
    if (d > r*r)
    {
        Contour(center, CurrentAtoms);
    }
}

float EMapBase::Rho(float x, float y, float z)
{
    if (!HasDensity())
    {
        return 0.0;
    }
    float fx, fy, fz;
    fx = x;
    fy = y;
    fz = z;
    transform(mapheader->ctof, &fx, &fy, &fz);
    return avgrho(fx, fy, fz);
}

float EMapBase::RDensity(MIAtomList &atoms)
{
    float r = 0.0;
    float fx, fy, fz, rho, zweight;
    float rr = 0.0;
    if (atoms.size() <= 0)
    {
        return 0;
    }
    if (!HasDensity())
    {
        return 0.0;
    }
    int natoms = 0;
    for (unsigned int i = 0; i < atoms.size(); i++)
    {
        if (atoms[i]->name()[0] != 'H')
        {
            fx = atoms[i]->x();
            fy = atoms[i]->y();
            fz = atoms[i]->z();
            transform(mapheader->ctof, &fx, &fy, &fz);
            rho = avgrho(fx, fy, fz)-25.0F;
            zweight = ZByName(atoms[i]->name())/6.7F;
            if (mapheader->resmin > 2.5F && rho > 75.0F)
            {
                rho = 75.0f + 0.25f*(rho-75.0f);
            }
            else
            {
                if (rho > 150.0f)
                {
                    rho = 150.0f + 0.25f*(rho-150.0f);
                }
            }
            // penalize heavily breaks
            if (rho < 0)
            {
                rho *= 5.0f;
            }
            rr += rho/50.0f/zweight;
            natoms++;
        }
    }
    if (natoms > 0)
    {
        r = rr/(float)natoms;
    }
    else
    {
        r = 0;
    }
    return r;
}

float EMapBase::RCorrelation(MIAtomList &atoms)
{
    float r = 0.0;
    float fx, fy, fz, rho, zweight;
    float rr = 0.0, zz = 0.0, rz = 0.0;
    if (atoms.size() <= 0)
    {
        return 0;
    }
    if (!HasDensity())
    {
        return 0.0;
    }

    for (unsigned int i = 0; i < atoms.size(); i++)
    {
        fx = atoms[i]->x();
        fy = atoms[i]->y();
        fz = atoms[i]->z();
        transform(mapheader->ctof, &fx, &fy, &fz);
        rho = avgrho(fx, fy, fz);
        zweight = ZByName(atoms[i]->name());
        if (mapheader->resmin > 2.5f && rho > 75.0f)
        {
            rho = 75.0f + 0.25f*(rho-75.0f);
        }
        rr += rho*rho;
        zz += zweight*zweight;
        rz += zweight*rho;
    }
    r = (rz/(float)sqrt(rr*zz));
    return r;
}

#define N 4
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
float EMapBase::maxdens(MIAtom *CAprev, MIAtom *CA1, int l, float fx0, float fy0, float fz0)
{
    MIAtom p, p2;
    float fx, fy, fz, fx1, fy1, fz1;
    float phi, theta, r, r_sum, r_min, r_hi, ds, fx_inc, fy_inc, fz_inc;
    int rl[10], i, j, k, kcum;
    double ideal_dist = 3.80;
    r_hi = INT_MIN;
    for (i = 0; i <= 180; i += 5)
    {
        theta = (float)(i*M_PI/180.0);
        p.setPosition(p.x(), p.y(), (float)(CA1->z() + ideal_dist*cos(theta)));
        ds = (float)(ideal_dist*sin(theta));
        for (j = 0; j < 360; j += 5)
        {
            phi = (float)(j*M_PI/180.0);
            p.setPosition((float)(CA1->x() + ds*cos(phi)),
                          (float)(CA1->y() + ds*sin(phi)), p.z());
            /*
               Caution:  Do not change "!(CalcAtomAngle(CAprev,CA1,&p) >= 60.0)"
               to "CalcAtomAngle(CAprev,CA1,&p) < 60.0".  The former protects us
               against bad values returned by AtomAngle().
             */
            if (!(CalcAtomAngle(*CAprev, *CA1, p) >= 60.0))
            {
                goto nextpoint;
            }
            if (l >= 0)
            {
                for (k = 0; k < l; k++)
                {
                    p2.setPosition(PList[k].x, PList[k].y, PList[k].z);
                    if (!(CalcAtomAngle(p2, *CA1, p) >= 30.0))
                    {
                        goto nextpoint;
                    }
                }
            }
            fx1 = p.x()*mapheader->ctof[0][0]+p.y()*mapheader->ctof[0][1]+p.z()*mapheader->ctof[0][2];
            fy1 = p.x()*mapheader->ctof[1][0]+p.y()*mapheader->ctof[1][1]+p.z()*mapheader->ctof[1][2];
            fz1 = p.x()*mapheader->ctof[2][0]+p.y()*mapheader->ctof[2][1]+p.z()*mapheader->ctof[2][2];
            fx_inc = (fx1 - fx0) / N;
            fy_inc = (fy1 - fy0) / N;
            fz_inc = (fz1 - fz0) / N;
            r_sum = 0.0f;
            kcum = 0;
            r_min = (float)INT_MAX;
            kcum = 0;
            for (k = 1; k <= N; k++)
            {
                kcum += k;
                fx = fx0 + (k*fx_inc);
                fy = fy0 + (k*fy_inc);
                fz = fz0 + (k*fz_inc);
                r = avgrho(fx, fy, fz);
                rl[k] = (int)r;
                r_sum += r*k;
                if (r < r_min)
                {
                    r_min = r;
                }
            }
            rl[N+1] = INT_MIN;
            r = r_min;
            if (r > r_hi)
            {
                r_hi = r;
                if (l >= 0)
                {
                    PList[l].x = p.x();
                    PList[l].y = p.y();
                    PList[l].z = p.z();
                    PList[l].r = r;
                    for (k = 1; true; k++)
                    {
                        RList[l][k] = (float)rl[k];
                        if (rl[k] == INT_MIN)
                        {
                            break;
                        }
                    }
                }
            }
nextpoint:
            ;
        }
    }
    return r_hi;
#undef  N
}

static bool mapp_compare(const MAP_POINT &lt, const MAP_POINT &rt)
{
    return lt.r > rt.r;
}

std::vector<MAP_POINT>&EMapBase::PlaceCA(MIAtom *CAprev, MIAtom *CA1)
{
    MAP_POINT p;
    MIAtom a;
    float fx, fy, fz, fx0, fy0, fz0, r, r_hi, r2;
    int l, l_hi;
    PList.clear();

    if (HasDensity())
    {
        mapheader->getfx(CA1, &fx0, &fy0, &fz0);
        r_hi = INT_MIN;
        for (l = 0; l < 5; l++)
        {
            PList.push_back(p);
            r = maxdens(CAprev, CA1, l, fx0, fy0, fz0);
            mapheader->getfx(&(PList[l]), &fx, &fy, &fz);
            a.setPosition(PList[l].x, PList[l].y, PList[l].z);
            r2 = maxdens(CA1, &a, -1, fx, fy, fz);
            if (r2 > 200)
            {
                r2 = 200;
            }
            r += r2;
            if (r > r_hi)
            {
                r_hi = r;
                l_hi = l;
            }
        }
        sort(PList.begin(), PList.end(), mapp_compare);
    }
    return PList;
}

void EMapBase::SetCurrentAtoms(MIAtomList *current)
{
    CurrentAtoms = current;
}

//typedef int (__cdecl *) (*comparfunc)();
static int rcompare(const void *i, const void *j)
{
    return ((PEAK*)j)->rho - ((PEAK*)i)->rho;
}

/* the center in integer map coords */
static int cix, ciy, ciz;

static int dcompare(const void *ii, const void *jj)
{
    int djx, djy, djz, dix, diy, diz;
    PEAK *i = (PEAK*)ii;
    PEAK *j = (PEAK*)jj;
    djx = j->ix - cix;
    djy = j->iy - ciy;
    djz = j->iz - ciz;
    dix = i->ix - cix;
    diy = i->iy - ciy;
    diz = i->iz - ciz;
    return ((dix*dix+diy*diy+diz*diz)- (djx*djx+djy*djy+djz*djz));
}

static void
fixcell(MINATOM *a)
{
    a->tx = a->ty = a->tz = 0.0;
    while (a->fx < 0.0)
    {
        a->fx += 1.0;
        a->tx += 1.0;
    }
    while (a->fx >= 1.0)
    {
        a->fx -= 1.0;
        a->tx -= 1.0;
    }
    while (a->fy < 0.0)
    {
        a->fy += 1.0;
        a->ty += 1.0;
    }
    while (a->fy >= 1.0)
    {
        a->fy -= 1.0;
        a->ty -= 1.0;
    }
    while (a->fz < 0.0)
    {
        a->fz += 1.0;
        a->tz += 1.0;
    }
    while (a->fz >= 1.0)
    {
        a->fz -= 1.0;
        a->tz -= 1.0;
    }
}

static void PutNearProtein(float fx, float fy, float fz, MINATOM *atoms, int natoms, float *bx, float *by, float *bz, CMapHeaderBase *mh)
{
    int j, isymm;
    float dx, dy, dz;
    float sx, sy, sz;
    float cx, cy, cz;
    float tx, ty, tz;
    float dbest = FLT_MAX, dsquared;
    /* find nearest atom */
    for (isymm = 0; isymm < mh->nsym; isymm++)
    {
        symm_mh(fx, fy, fz, &sx, &sy, &sz, mh, isymm);
        for (j = 0; j < natoms; j++)
        {
            cx = atoms[j].cx;
            cy = atoms[j].cy;
            cz = atoms[j].cz;
            transform(mh->ctof, &cx, &cy, &cz);
            dx = sx-cx;
            dy = sy-cy;
            dz = sz-cz;
            while (dx < -0.5)
            {
                sx += 1.0;
                dx += 1.0;
            }
            while (dx > 0.5)
            {
                sx -= 1.0;
                dx -= 1.0;
            }
            while (dy < -0.5)
            {
                sy += 1.0;
                dy += 1.0;
            }
            while (dy > 0.5)
            {
                sy -= 1.0;
                dy -= 1.0;
            }
            while (dz < -0.5)
            {
                sz += 1.0;
                dz += 1.0;
            }
            while (dz > 0.5)
            {
                sz -= 1.0;
                dz -= 1.0;
            }
            transform(mh->ftoc, &dx, &dy, &dz);
            dsquared = dx*dx+dy*dy+dz*dz;
            if (dsquared < dbest)
            {
                tx = sx;
                ty = sy;
                tz = sz;
                transform(mh->ftoc, &tx, &ty, &tz);
                *bx = tx;
                *by = ty;
                *bz = tz;
                dbest = dsquared;
            }
        }
    }
}

static int FarCheck(float cx, float cy, float cz, MINATOM *atoms, int natoms, float dmax)
{
    /* returns true if coords closer than 5.0 A to another atom */
    int i;
    float dx, dy, dz, d;
    //float dmin =(float)atof((char*)xv_get(waterpop->maxdistance,PANEL_VALUE));;

    float dmin = dmax*dmax;
    for (i = 0; i < natoms; i++)
    {
        dx = cx - atoms[i].cx;
        dy = cy - atoms[i].cy;
        dz = cz - atoms[i].cz;
        d = dx*dx+dy*dy+dz*dz;
        if (d < dmin)
        {
            return 1;
        }
    }
    return 0;
}

static int NearCheck(float fx, float fy, float fz, MINATOM *atoms, int natoms, CMapHeaderBase *mh, int add_water, float dmin)
{
    float dC = 3.0*3.0;
    float dO = dmin*dmin;
    float dsquared;
    float dx, dy, dz;
    float fdx, fdy, fdz;
    int j;
    char type;

    if (add_water == 0)
    {
        dC = dO;
    }
    while (fx < 0.0)
    {
        fx += 1.0;
    }
    while (fx >= 1.0)
    {
        fx -= 1.0;
    }
    while (fy < 0.0)
    {
        fy += 1.0;
    }
    while (fy >= 1.0)
    {
        fy -= 1.0;
    }
    while (fz < 0.0)
    {
        fz += 1.0;
    }
    while (fz >= 1.0)
    {
        fz -= 1.0;
    }

    /* 3.3 Angstroms in fractional units */
    fdx = 3.3f/mh->a;
    fdy = 3.3f/mh->b;
    fdz = 3.3f/mh->c;
    /* look for an atom nearby */
    for (j = 0; j < natoms; j++)
    {
        dx = fx-atoms[j].fx;
        dy = fy-atoms[j].fy;
        dz = fz-atoms[j].fz;
        /* take care of cell-wrapping */
        while (dx < -0.5)
        {
            dx += 1.0;
        }
        while (dx > 0.5)
        {
            dx -= 1.0;
        }
        while (dy < -0.5)
        {
            dy += 1.0;
        }
        while (dy > 0.5)
        {
            dy -= 1.0;
        }
        while (dz < -0.5)
        {
            dz += 1.0;
        }
        while (dz > 0.5)
        {
            dz -= 1.0;
        }
        if (fabs(dx) < fdx && fabs(dy) < fdy && fabs(dz) < fdz)
        {
            /* within parallelopided - worth
             * further processing
             */
            transform(mh->ftoc, &dx, &dy, &dz);
            dsquared = dx*dx+dy*dy+dz*dz;
            type = atoms[j].type;
            if (type == 'C' || type == 'N' || type == 'O')
            {
                if (type == 'C' && dsquared < dC)
                {
                    return false;
                }
                if (type == 'N' && dsquared < dO)
                {
                    return false;
                }
                if (type == 'O' && dsquared < dO)
                {
                    return false;
                }
            }
            else
            {
                if (dsquared < dC)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

int EMapBase::HydrateMap(int minlevel, int maxadd, int add_water, MIMoleculeBase *model,
                         float dmin, float dmax, float xmin, float xmax,
                         float ymin, float ymax, float zmin, float zmax)
{
    /* add-water controls whether we add waters or fill the map with
     * psuedo-carbon atoms */
    int nadd = 0;
    int ix, iy, iz;
    float fx, fy, fz;
    float cx, cy, cz;
    Residue *res;
    Residue *residues = model->getResidues();
    CMapHeaderBase *mh = mapheader;
    //MINATOM * atoms=NULL;
    MINATOM atom, symatom;
    vector<MINATOM> atoms;
    vector<MINATOM> symatoms;
    PEAK peak;
    vector<PEAK> peaks;
    int i, j, kept = 0;
    int sx, sy, sz, ex, ey, ez;
    MIAtom *newatom;
    /*if(!add_water){
       dmin = 1.45F;
       dmax = 50.0F;
       }*/
    int imin;
    std::string buf;
    float dsquared;

    if (!HasDensity())
    {
        Logger::log("There is no map at position %d\n", mapnumber+1);
        return 0;
    }
    if (minlevel < 1)
    {
        minlevel = 1;
    }
    if (add_water)
    {
        Logger::log("HydrateMap: Minlevel = %d\n", minlevel);
    }
    else
    {
        Logger::log("FillMap: Minlevel = %d\n", minlevel);
    }
    /* convert the residue list into a list of minatoms including
     * all the symmetry mates
     */
    unsigned int num_atoms = 0;
    for (res = residues; res != NULL; res = res->next())
    {
        for (j = 0; j < res->atomCount(); j++)
        {
            // ignore Hydrogens
            if (res->atom(j)->name()[0] == 'H')
            {
                continue;
            }
            num_atoms++;
        }
    }
    if (num_atoms == 0)
    {
        Logger::log("There must be at least one non-H atom in model\n");
        return 0;
    }
    Logger::log("Protein has %d non-H atoms\n", num_atoms);
    unsigned int num_symatoms = (mh->nsym-1) * (num_atoms+maxadd);
    try
    {
        atoms.reserve(num_atoms+maxadd);
        symatoms.reserve(num_symatoms);
    }
    catch (bad_alloc)
    {
        Logger::log("Not enough memory in HydrateMap");
        return 0;
    }
    Logger::log("Protein has %d symatoms\n", num_symatoms-maxadd);
    Logger::footer("Protein has %d symatoms\n", num_symatoms-maxadd);
    i = 0;
    for (res = residues; res != NULL; res = res->next())
    {
        for (j = 0; j < res->atomCount(); j++)
        {
            // ignore Hydrogens
            if (res->atom(j)->name()[0] == 'H')
            {
                continue;
            }
            atom.fx = res->atom(j)->x();
            atom.fy = res->atom(j)->y();
            atom.fz = res->atom(j)->z();
            atom.cx = res->atom(j)->x();
            atom.cy = res->atom(j)->y();
            atom.cz = res->atom(j)->z();
            atom.symm = 0;
            atom.type = res->atom(j)->name()[0];
            transform(mh->ctof, &atom.fx, &atom.fy, &atom.fz);
            fixcell(&atom);
            i++;
            atoms.push_back(atom);
        }
    }
    for (j = 1; j < mh->nsym; j++)
    {
        for (i = 0; (unsigned int)i < atoms.size(); i++)
        {
            //jj = (j-1)*natoms+i;
            fx = atoms[i].cx;
            fy = atoms[i].cy;
            fz = atoms[i].cz;
            transform(mh->ctof, &fx, &fy, &fz);
            symm_mh(fx, fy, fz,
                    &symatom.fx, &symatom.fy, &symatom.fz,
                    mh, j);
            symatom.cx = symatom.fx;
            symatom.cy = symatom.fy;
            symatom.cz = symatom.fz;
            transform(mh->ftoc, &(symatom.cx), &(symatom.cy), &(symatom.cz));
            fixcell(&symatom);
            symatom.symm = j;
            symatom.type = atoms[i].type;
            symatoms.push_back(symatom);
        }
    }
    sx = ROUND(xmin * mh->nx);
    sy = ROUND(ymin * mh->ny);
    sz = ROUND(zmin * mh->nz);
    ex = ROUND(xmax * mh->nx);
    ey = ROUND(ymax * mh->ny);
    ez = ROUND(zmax * mh->nz);
    Logger::log("Search volume x=%0.2f,%0.2f y=%0.2f,%0.2f z=%0.2f,%0.2f\n", xmin, xmax, ymin, ymax, zmin, zmax);
    Logger::footer("Search volume x=%0.2f,%0.2f y=%0.2f,%0.2f z=%0.2f,%0.2f\n", xmin, xmax, ymin, ymax, zmin, zmax);
    for (ix = sx; ix < ex; ix++)
    {
        for (iy = sy; iy < ey; iy++)
        {
            for (iz = sz; iz < ez; iz++)
            {
                fx = (float)ix/(float)mh->nx;
                fy = (float)iy/(float)mh->ny;
                fz = (float)iz/(float)mh->nz;
                if (avgrho(fx, fy, fz) >= (float)minlevel)
                {
                    peak.ix = ix;
                    peak.iy = iy;
                    peak.iz = iz;
                    peak.rho = (int)avgrho(fx, fy, fz);
                    peaks.push_back(peak);
                }
            }
        }
    }
    Logger::log("Peaks before filtering = %d\n", (int)peaks.size());
    Logger::footer("Peaks before filtering = %d\n", (int)peaks.size());
    if (add_water)
    {
        /* sort the peaks by size */
        qsort(&peaks[0], peaks.size(), sizeof(PEAK), rcompare);
    }
    else
    {
        /* sort the peaks by size */
        qsort(&peaks[0], peaks.size(), sizeof(PEAK), dcompare);
    }

    /* work our way through the peak list */
    /* and find the unique peaks */
    dsquared = dmin * dmin;
    imin = ROUND(dmin/mh->a*(float)mh->nx);
    int psize;
    psize = peaks.size();
    for (i = 0; i < psize; i++)
    {
        if (peaks[i].rho > 0)
        {
            for (j = 0; j < psize; j++)
            {
                if (j != i && peaks[j].rho > 0)
                {
                    ix = abs(peaks[i].ix-peaks[j].ix);
                    if (ix <= imin)
                    {
                        iy = peaks[i].iy-peaks[j].iy;
                        iz = peaks[i].iz-peaks[j].iz;
                        fx = (float)ix/(float)mh->nx;
                        fy = (float)iy/(float)mh->ny;
                        fz = (float)iz/(float)mh->nz;
                        transform(mh->ftoc, &fx, &fy, &fz);
                        if (fx*fx+fy*fy+fz*fz <= dsquared)
                        {
                            peaks[j].rho = -peaks[j].rho;
                        }
                    }
                }
            }
        }
    }
    for (i = 0; (unsigned int)i < peaks.size(); i++)
    {
        if (peaks[i].rho > 0)
        {
            kept++;
        }
    }
    Logger::log("Kept %d peaks after sorting\n", kept);
    Logger::footer("Kept %d peaks after sorting\n", kept);

    /* evaluate peaks for wateryness */
    Residue *focusres;

    for (i = 0; (unsigned int)i < peaks.size(); i++)
    {
        if (peaks[i].rho > 0)
        {
            fx = (float)peaks[i].ix/(float)mh->nx;
            fy = (float)peaks[i].iy/(float)mh->ny;
            fz = (float)peaks[i].iz/(float)mh->nz;
            if (NearCheck(fx, fy, fz, &atoms[0], atoms.size(), mh, add_water, dmin))
            {
                if (NearCheck(fx, fy, fz, &symatoms[0], symatoms.size(), mh, add_water, dmin))
                {
                    PutNearProtein(fx, fy, fz, &atoms[0], atoms.size(), &cx, &cy, &cz, mh);
                    if (FarCheck(cx, cy, cz, &atoms[0], atoms.size(), dmax))
                    {
                        focusres = model->AddWater(cx, cy, cz, false);
                        if (focusres)
                        {
                            if (add_water)
                            {
                                Logger::log("Add water %s at %f %f %f\n", focusres->name().c_str(), cx, cy, cz);
                                Logger::footer("Add water %s at %f %f %f\n", focusres->name().c_str(), cx, cy, cz);
                            }
                            /* refine the water position */
                            //if(add_water) TranslationalSearchNotify(0,0);
                            /* add to the errorpop list */
                            newatom = focusres->atom(0);
                            res = focusres;
                            nadd++;
                            /* Must also add to atoms list
                             * to prevent adding two waters
                             * on top of each other *
                             */
                            atom.fx = newatom->x();
                            atom.fy = newatom->y();
                            atom.fz = newatom->z();
                            atom.cx = newatom->x();
                            atom.cy = newatom->y();
                            atom.cz = newatom->z();
                            atom.symm = 0;
                            atom.type = 'O';
                            transform(mh->ctof, &atom.fx, &atom.fy, &atom.fz);
                            fixcell(&atom);
                            for (j = 1; j < mh->nsym; j++)
                            {
                                symm_mh(atom.fx, atom.fy, atom.fz, &symatom.fx, &symatom.fy, &symatom.fz, mh, j);
                                symatom.cx = symatom.fx;
                                symatom.cy = symatom.fy;
                                symatom.cz = symatom.fz;
                                transform(mh->ftoc, &(symatom.cx), &(symatom.cy), &(symatom.cz));
                                fixcell(&symatom);
                                symatom.symm = j;
                                symatom.type = atom.type;
                                symatoms.push_back(symatom);
                            }
                            atoms.push_back(atom);
                            if (nadd >= maxadd)
                            {
                                goto cleanup;
                            }
                        }
                    }
                }
            }
        }
    }
cleanup:
    model->Build();
    if (add_water)
    {
        Logger::log("HydrateMap: Added %d waters", nadd);
        Logger::footer("HydrateMap: Added %d waters", nadd);
    }
    else
    {
        Logger::log("HydrateMap: Added %d atoms", nadd);
        Logger::footer("HydrateMap: Added %d atoms", nadd);
    }
    return nadd;
}

bool EMapBase::CalcSolventMask(float percentSolvent, float radius)
{
    if (SigmaMap())
    {
        if (SmoothMap(radius))
        {
            if (ScaleSolventPercent(percentSolvent))
            {
            }
        }
    }
    return false;
}

bool EMapBase::SigmaMap()
{
    int ix, iy, iz;
    int nx = mapheader->nx;
    int ny = mapheader->ny;
    int nz = mapheader->nz;
    long nmap = map_points.capacity();
    if (nmap < 9)
    {
        return false;
    }
    float rho, rhosq, sigma;
    std::vector<float> map2;
    try
    {
        map2.resize(nmap);
    }
    catch (...)
    {
        return false;
    }
    for (ix = 0; ix < nx; ix++)
    {
        for (iy = 0; iy < ny; iy++)
        {
            for (iz = 0; iz < nz; iz++)
            {
                rhosq = 0.0;
                rho =  map_points[mapdex(ix, iy, iz)];
                rhosq += rho*rho;
                rho =  map_points[mapdex(ix+1, iy, iz)];
                rhosq += rho*rho;
                rho =  map_points[mapdex(ix-1, iy, iz)];
                rhosq += rho*rho;
                rho =  map_points[mapdex(ix, iy+1, iz)];
                rhosq += rho*rho;
                rho =  map_points[mapdex(ix, iy-1, iz)];
                rhosq += rho*rho;
                rho =  map_points[mapdex(ix, iy, iz+1)];
                rhosq += rho*rho;
                rho =  map_points[mapdex(ix, iy, iz-1)];
                rhosq += rho*rho;

                sigma = sqrt((rhosq-rhosq/7.0f)/6.0f);

                map2[mapdex(ix, iy, iz)] = sigma;
            }
        }
    }
    copy(map2.begin(), map2.end(), map_points.begin());
    return true;
}

bool EMapBase::SmoothMap(float radius)
{
    int i, ix, iy, iz;
    float sum = 0.0;
    int nx = mapheader->nx;
    int ny = mapheader->ny;
    int nz = mapheader->nz;
    long nmap = map_points.capacity();
    if (nmap < 9)
    {
        return false;
    }
    std::vector<float> map2;
    try
    {
        map2.resize(nmap);
    }
    catch (...)
    {
        return false;
    }
    int ntimes = ROUND(radius/(mapheader->a/(float)mapheader->nx));
    for (i = 0; i < ntimes; i++)
    {
        for (ix = 0; ix < nx; ix++)
        {
            for (iy = 0; iy < ny; iy++)
            {
                for (iz = 0; iz < nz; iz++)
                {
                    sum =    2.0F*map_points[mapdex(ix-1, iy, iz)]
                          +2.0F*map_points[mapdex(ix+1, iy, iz)]
                          +2.0F*map_points[mapdex(ix, iy-1, iz)]
                          +2.0F*map_points[mapdex(ix, iy+1, iz)]
                          +2.0F*map_points[mapdex(ix, iy, iz-1)]
                          +2.0F*map_points[mapdex(ix, iy, iz+1)]
                          +map_points[mapdex(ix+1, iy+1, iz)]
                          +map_points[mapdex(ix-1, iy+1, iz)]
                          +map_points[mapdex(ix+1, iy-1, iz)]
                          +map_points[mapdex(ix-1, iy-1, iz)]
                          +map_points[mapdex(ix+1, iy, iz+1)]
                          +map_points[mapdex(ix-1, iy, iz+1)]
                          +map_points[mapdex(ix+1, iy, iz-1)]
                          +map_points[mapdex(ix-1, iy, iz-1)]
                          +map_points[mapdex(ix, iy+1, iz+1)]
                          +map_points[mapdex(ix, iy-1, iz+1)]
                          +map_points[mapdex(ix, iy+1, iz-1)]
                          +map_points[mapdex(ix, iy-1, iz-1)]
                          +6.0F*map_points[mapdex(ix, iy, iz)] ;
                    map2[mapdex(ix, iy, iz)] = sum/30.0f;
                }
            }
        }
        copy(map2.begin(), map2.end(), map_points.begin());
    }
    return true;
}

bool EMapBase::ScaleSolventPercent(float percentSolvent)
{
    if (percentSolvent < 0.0 || percentSolvent >= 1.0)
    {
        return false;
    }
    unsigned int nmap = map_points.capacity();
    std::vector<float> map2;
    try
    {
        map2.resize(nmap);
    }
    catch (...)
    {
        return false;
    }
    copy(map_points.begin(), map_points.end(), map2.begin());
    sort(map2.begin(), map2.end());
    float s0 = map2[(unsigned int)((float)nmap*percentSolvent)];
    float map_max = map2[nmap-1];
    float map_min = map2[0];
    float solv_add = 0 - s0;
    float solv_scale;
    if ((s0-map_min) != 0)
    {
        solv_scale = 100.0F/(s0-map_min);
    }
    else
    {
        solv_scale = 1.0F;
    }
    float prot_add = solv_add;
    float prot_scale;
    if ((map_max-s0) != 0)
    {
        prot_scale = 100.0F/(map_max-s0);
    }
    else
    {
        prot_scale = 1.0F;
    }
    mapmin = mapmax = map_points[0];
    for (unsigned int i = 0; i < nmap; i++)
    {
        if (map_points[i] < s0)
        {
            // solvent - put in range -100 to 0
            map_points[i] += solv_add;
            map_points[i] *= solv_scale;
        }
        else
        {
            // protein - put in range 0 to 100
            map_points[i] += prot_add;
            map_points[i] *= prot_scale;
        }
        if (map_points[i] > mapmax)
        {
            mapmax = map_points[i];
        }
        if (map_points[i] < mapmin)
        {
            mapmin = map_points[i];
        }
    }
    Logger::log("max = %f, min = %f", (float)mapmax, (float)mapmin);
    return true;
}

static const char *truth(bool t)
{
    if (t)
    {
        return "True";
    }
    else
    {
        return "False";
    }
}

std::string EMapBase::Info()
{
    std::string s;
    std::string ops;
    for (int i = 0; i < mapheader->nsym; i++)
    {
        ops += mapheader->SymopsString[i];
        ops += "\n";
    }

    char buf[1024];
    sprintf(buf, "ID: %s\n"
            "Pathname: %s\n"
            "Type: %s\n"
            "Has Phases: %s\n"
            "No. refls: %d\n"
            "Data Resolution Max-Min: %0.4f-%0.4f\n"
            "FFT Resolution Max-Min: %0.4f-%0.4f\n"
            "Has Density: %s\n"
            "Grid Size, nx, ny, nz: %d %d %d\n"
            "Map Center: %0.3f %0.3f %0.3f\n"
            "Radius: %0.3f\n"
            "Scale: %0.3f\n"
            "Visible: %s\n"
            "Crystal %s\n"
            "Unit Cell: %0.3f %0.3f %0.3f  %0.3f %0.3f %0.3f\n"
            "Spacegroup: %s\n"
            //"Spacegroup No.: %d\n"
            "Symmops: %s\n"
            "Reciprocal space extents: hmin %d hmax %d kmin %d kmax %d lmax %d lmin %d",
            MapID().c_str(), pathName.c_str(), maptypes[(mapheader->maptype)], truth(HasPhases()),
            (int)refls.size(), 0.5F/refls_stholmin, 0.5F/refls_stholmax, mapheader->resmax, mapheader->resmin,
            truth(HasDensity()),
            mapheader->nx, mapheader->ny, mapheader->nz, map_center[0], map_center[1], map_center[2],
            settings->Radius, scale, truth(visible != 0), Crystal.c_str(), mapheader->a, mapheader->b, mapheader->c,
            mapheader->alpha, mapheader->beta, mapheader->gamma,
            mapheader->SymInfoString.c_str(), ops.c_str(),
            mapheader->hmin, mapheader->hmax, mapheader->kmin, mapheader->kmax,
            mapheader->lmin, mapheader->lmax);
    return std::string(buf);
}

bool EMapBase::IsScaFile(const char *pathname)
{
    std::string ext;
    MISplitPath(pathname, 0, 0, &ext);
    return (strcasecmp(ext.c_str(), "sca") == 0);
}

bool EMapBase::IsPhsFile(const char *pathname)
{
    std::string ext;
    MISplitPath(pathname, 0, 0, &ext);
    return (strcasecmp(ext.c_str(), "phs") == 0);
}

bool EMapBase::IsFinFile(const char *pathname)
{
    std::string ext;
    MISplitPath(pathname, 0, 0, &ext);
    return (strcasecmp(ext.c_str(), "fin") == 0);
}

bool EMapBase::IsRefFile(const char *pathname)
{
    std::string ext;
    MISplitPath(pathname, 0, 0, &ext);
    return (strcasecmp(ext.c_str(), "ref") == 0);
}

void EMapBase::RecalcResolution()
{
    sthol(1, 1, 1, mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma, 1);
    refls_stholmax = -1.0F;
    refls_stholmin = 999999.0F;
    for (unsigned int i = 0; i < refls.size(); i++)
    {
        refls[i].sthol = sthol(refls[i].ind[0], refls[i].ind[1], refls[i].ind[2],
                               mapheader->a, mapheader->b, mapheader->c,
                               mapheader->alpha, mapheader->beta, mapheader->gamma, 0);
        if (refls[i].sthol > refls_stholmax)
        {
            refls_stholmax = refls[i].sthol;
        }
        if (refls[i].sthol < refls_stholmin)
        {
            refls_stholmin = refls[i].sthol;
        }
    }
    mapheader->resmax = 0.5F/refls_stholmin;
    mapheader->resmin = 0.5F/refls_stholmax;
    mapheader->nx = MIMapFactor((int)(3.0F*mapheader->a/mapheader->resmin), FFT_PRIME, EVEN, 2);
    mapheader->ny = MIMapFactor((int)(3.0F*mapheader->b/mapheader->resmin), FFT_PRIME, ODDOREVEN, 1);
    mapheader->nz = MIMapFactor((int)(3.0F*mapheader->c/mapheader->resmin), FFT_PRIME, ODDOREVEN, 1);
}

long EMapBase::AddFreeRFlag(int percent, bool use_shells)
{
    long nadded = 0;
    int range = ROUND(100.0F/(float)percent);
    if (range < 0)
    {
        range = 1;
    }
    if (use_shells)
    {
        Logger::message("Not implemented yet!");
    }
    else
    {
        int flag;
        for (size_t i = 0; i < refls.size(); i++)
        {
            flag = irand(range);
            refls[i].freeRflag = flag;
            if (flag == 0)
            {
                nadded++;
            }
        }
    }
    Logger::log("%0.1f percent data are now in R-Free set 0", (float)nadded/(float)refls.size()*100.0F);
    FreeRSet = true;
    SetModified(true);
    return nadded;
}

float EMapBase::CorrScore(MIAtomList atoms)
{
    micomplex *rho;
    //float sc;
    float corr;
    int i, nx, ny, nz;
    //char buf[2000];
    CMapHeaderBase *mh = mapheader;
    float fofcsum = 0.0;
    float fosum = 0.0, fcsum = 0.0;
    float fofosum = 0.0, fcfcsum = 0.0;
    //float Bscale;

    MIAtom *a;
    int n = 0;
    Residue res;
    res.setType("ALA");
    sfinit();

    Bscale = std::max(exp(mh->resmin+.1)/2.0, 0.0) + 2.0;

    /* first figure out the size of the map given the cell
     * and the resolution desired
     */
    nx = mapheader->nx;
    ny = mapheader->ny;
    nz = mapheader->nz;

    rho = (micomplex*)malloc(nx*ny*nz*sizeof(micomplex));

    for (i = 0; i < nx*ny*nz; i++)
    {
        rho[i].r = 0.0;
        rho[i].i = 0.0;
    }

    for (unsigned int i = 0; i < atoms.size(); i++)
    {
        a = atoms[i];
        addrho(a, rho, nx, ny, nz, mh, &res);
    }

    //if((rho = buildrho(res, mh, wait))==NULL)
    //  return(0);
    InvertMap(rho, mh, refls);
    free((void*)rho);
    //sc = ComputeScale( refls, mh);
    //ApplyScale(refls, sc, mh);
    n = 0;
    for (i = 0; (unsigned int)i < refls.size(); i++)
    {
        if (refls[i].sthol > 0.5/mh->resmin
            || refls[i].sthol < 0.5/mh->resmax)
        {
            continue;
        }
        refls[i].fc = (float)sqrt(refls[i].acalc*refls[i].acalc
                                  +refls[i].bcalc*refls[i].bcalc);
        fofcsum += refls[i].fc * refls[i].fo;
        fofosum += refls[i].fo * refls[i].fo;
        fcfcsum += refls[i].fc * refls[i].fc;
        fosum += refls[i].fo;
        fcsum += refls[i].fc;
        n++;
    }
    corr = ((float)n*fofcsum - fosum*fcsum)/ (float)sqrt(((double)n*fcfcsum-fcsum*fcsum)
                                                         *((double)n*fofosum-fosum*fosum));

    FcsValid = false; // yup, this function dorks them up.
    //corr = (ntot*sum22T - sum20T*sum02T)/
    //  sqrt((ntot*sum40T-sum20T*sum20T)*(ntot*sum04T-sum02T*sum02T));
    return (corr);
}

long EMapBase::Reindex(int index_mat[][3])
{
    long nindex = 0;

    if (refls.size() == 0)
    {
        return 0;
    }

    int h, k, l, hr, kr, lr;
    for (size_t i = 0; i < refls.size(); i++)
    {
        h = refls[i].ind[0];
        k = refls[i].ind[1];
        l = refls[i].ind[2];
        hr = h*index_mat[0][0] + k*index_mat[1][0] + l*index_mat[2][0];
        kr = h*index_mat[0][1] + k*index_mat[1][1] + l*index_mat[2][1];
        lr = h*index_mat[0][2] + k*index_mat[1][2] + l*index_mat[2][2];
        refls[i].ind[0] = hr;
        refls[i].ind[1] = kr;
        refls[i].ind[2] = lr;
        nindex++;
    }
    Logger::log("Re-indexed %d reflections by:", nindex);
    Logger::log("hnew = h*%2d + k*%2d + l*%2d\n", index_mat[0][0], index_mat[1][0], index_mat[2][0]);
    Logger::log("knew = h*%2d + k*%2d + l*%2d\n", index_mat[0][1], index_mat[1][1], index_mat[2][1]);
    Logger::log("lnew = h*%2d + k*%2d + l*%2d\n", index_mat[0][2], index_mat[1][2], index_mat[2][2]);

    SetModified(true);

    return nindex;
}

