#include <cstdio>
#include <cmath>
#include <cstring>
#ifdef __WXGTK__
#include <values.h>
#endif

#include <nongui/nonguilib.h>
#include "maplib.h"


#include "maptypes.h"
#include "fft.h"
#include "sfcalc.h"


#define c_mul(a, b, c) cmul(a.re, a.im, b.re, b.im, c.re, c.im)
#define cmul(a_x, a_y, b_x, b_y, c_x, c_y) { \
        register float _cmul1_, _cmul2_, _cmul3_, _cmul4_; \
        _cmul1_ = (a_x); \
        _cmul2_ = (a_y); \
        _cmul3_ = (b_x); \
        _cmul4_ = (b_y); \
        c_x = _cmul1_*_cmul3_ - _cmul2_*_cmul4_; \
        c_y = _cmul2_*_cmul3_ + _cmul1_*_cmul4_; \
}


static double raddeg;
static float g11, g12, g13, g22, g23, g33, dsqmin, dsqmax;
static float a, b, c, alpha, beta, gama, rmin, rmax, vol;
/*
   extern int r[3][4][96];
 */
static int rsym[3][4][MISymmop::MAXSYMMOPS];
static int nx, ny, nz, hmax, kmax, lmax;
static int xmin, ymin, zmin, xmax, ymax, zmax;
static int hin[4], hout[4];
static int nsymmops;
static int ix, iy, iz;
static int maptype = 0;
static char buf[2000];

static fcomplex *membuf;
/*
   extern nextf (), symops(), xpnd (), expansion_error ();
   extern void read_sf ();
   extern cycle (), permute_everything(), get_sectioning();
   extern get_unit_cell();
   extern wplane (), list_operator();
   extern report_error(), new_term(), initialize_term();
   extern cexp (), str_index(), my_index();
 */

float*
fft3d(CREFL *refl, int nrefl, CMapHeaderBase *mapheader, int usepsi)
{
    float scale, temp;
    int i, j, nyblk, nzblk, k, kl = 0, ku = 0;
    long int d[6];
    int offset, zl, zu, z1, z2, yl, yu, y_1, y_2;
    /*char title[81];*/
    long int jj;
    int iz1;


    a = mapheader->a;
    b = mapheader->b;
    c = mapheader->c;
    alpha = mapheader->alpha;
    beta = mapheader->beta;
    gama = mapheader->gamma;
    maptype = mapheader->maptype;
    get_unit_cell();

    for (k = 0; k < mapheader->nsym; k++)
    {
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                rsym[i][j][k] = (int)mapheader->symops[i][j][k];
            }
            /* translation is in 12ths */
            /* fixed 1-27-98 to round properly for spacegroups wqith 1/6 or 1/3 symmetry ops as suggested by K. Cowtan */
            rsym[i][3][k] = (int)(12.0*(mapheader->symops[i][3][k])+12.5)-12;
        }
    }
    nsymmops = mapheader->nsym -1; /* floatly the max index not number of ops */
#ifdef DEBUG
    sprintf(buf, "%d symmops =\n ", nsymmops+1);
    Logger::log(buf);
    for (k = 0; k <= nsymmops; k++)
    {
        sprintf(buf, "Symmetry operator:\n");
        Logger::log(buf);
        for (i = 0; i < 3; i++)
        {
            sprintf(buf, "  %3d %3d %3d %3d\n", rsym[i][0][k], rsym[i][1][k], rsym[i][2][k], rsym[i][3][k]);
            Logger::log(buf);
        }
    }
#endif

    nx = mapheader->nx;
    ny = mapheader->ny;
    nz = mapheader->nz;
    xmin = ymin = zmin = 0;
    xmax = nx -1;
    ymax = ny -1;
    zmax = nz -1;
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    {
        sprintf(buf, "Sampling must be at least 1 in each dimension.\n");
        Logger::log(buf);
        sprintf(buf, "Edges: %d %d %d\n", nx, ny, nz);
        Logger::log(buf);
        return (NULL);
    }
    if (nx%2 != 0)
    {
        Logger::log("NX must be even");
        return (NULL);
    }

    if (mapheader->resmin != 0.0)
    {
        dsqmax = 1.0F/(mapheader->resmin*mapheader->resmin);
        rmin = mapheader->resmin;
    }
    else
    {
        dsqmax = 1.0F/(2.0F*2.0F);
        rmin = 2.0F;
    }
    if (mapheader->resmax != 0.0)
    {
        dsqmin = 1.0F/(mapheader->resmax*mapheader->resmax);
        rmax = mapheader->resmax;
    }
    else
    {
        dsqmin = 0.0;
        rmax = 1000.0F;
    }
    hmax = mapheader->hmax;
    kmax = mapheader->kmax;
    lmax = mapheader->lmax;
#ifdef GLURP
    rmin = std::min(a/(float)hmax, b/(float)kmax);
    rmin = std::min(c/(float)lmax, rmin);
    dsqmax = 1.0F/rmin;
#endif
    if (nx == 1)
    {
        hmax = 0;
    }
    if (ny == 1)
    {
        kmax = 0;
    }
    if (nz == 1)
    {
        lmax = 0;
    }
    if ((nx <= 2*hmax) || (ny <= 2*kmax) || (nz <= 2*lmax))
    {
        sprintf(buf, "Grid too coarse for requested resolution\n");
        Logger::log(buf);
        sprintf(buf, "Resolution (d* max) = %f\n", 1.0/sqrt(dsqmax));
        Logger::log(buf);
        sprintf(buf, "Maximum indices = %d %d %d\n", hmax, kmax, lmax);
        Logger::log(buf);
        sprintf(buf, "Grid = %d %d %d\n", nx, ny, nz);
        Logger::log(buf);
        return (NULL);
    }

#ifdef DEBUG
    sprintf(buf, "Resolution: %5.2f Angstroms\n", rmin);
    Logger::log(buf);
    sprintf(buf, "%6.2f Angstroms inner limit\n", rmax);
    Logger::log(buf);
    sprintf(buf, " Index limits: h max = %d, k max = %d, l max = %d\n", hmax, kmax, lmax);
    Logger::log(buf);
    sprintf(buf, "Grid information:      X    Y    Z\n");
    Logger::log(buf);
    sprintf(buf, " Sampling          %5d%5d%5d\n", nx, ny, nz);
    Logger::log(buf);
    sprintf(buf, " Minimum boundary  %5d%5d%5d\n", xmin, ymin, zmin);
    Logger::log(buf);
    sprintf(buf, " Maximum boundary  %5d%5d%5d\n", xmax, ymax, zmax);
    Logger::log(buf);
#endif

    scale = 1.0;
    temp = 0.0;
    scale /= vol;

    /* if a SigmaA map we need to compute coefficients */
    if (maptype == (int)MIMapType::TwoFoFcSigmaA || maptype == (int)MIMapType::FoFcSigmaA)
    {
        if (SigmaA_C(refl, nrefl, mapheader) == 0)
        {
            return NULL;
        }
        /* SigmaA takes care of error message */
    }

    membuf = (fcomplex*)calloc((nx)*(ny)*(nz), sizeof(fcomplex));
    if (membuf == NULL)
    {
        sprintf(buf, "Error in fft3d: memory allocation...\n");
        Logger::log(buf);
        return (NULL);
    }

    if ((zmax - zmin + 1) >= nz)
    {
        nzblk = 1;
        zl = 0;
        zu = nz - 1;
    }
    else
    {
        nzblk = 1;
        if (zmax >= nz)
        {
            nzblk = 2;
        }
        zl = zmin;
        zu = zmax;
    }

    if ((ymax - ymin + 1) >= ny)
    {
        nyblk = 1;
        yl = 0;
        yu = ny - 1;
    }
    else
    {
        nyblk = 1;
        if (ymax >= ny)
        {
            nyblk = 2;
        }
        yl = ymin;
        yu = ymax;
    }

    read_sf((float*)membuf, nx/2, ny, nz, scale, temp, refl, nrefl, usepsi);

    /*      Transform on l -- to reduce paging, transform each k-value
     *      separately.
     */
    for (i = 0; i < 2; i++)
    {
        if (i == 0)
        {
            kl = 1;
            ku = kmax + 1;
        }
        if (i == 1)
        {
            kl = ny - kmax + 1;
            ku = ny;
        }
        d[0] = nx*ny*nz;
        d[1] = nx*ny;
        d[2] = d[0];
        d[3] = 2*hmax + 2;
        d[4] = 2;
        for (k = kl; k <= ku; k++)
        {
            long int nz_long = (long int) nz;
            offset = (k - 1)*nx;
            cmplft_((float*)membuf+offset, (float*)membuf+offset+1, &nz_long, d);
        }
    }

    /*      transform on k */

    z1 = zl;
    z2 = std::min(zu, nz-1);
    for (iz = 0; iz < nzblk; iz++)
    {
        d[0] = nx*ny;
        d[1] = nx;
        d[2] = nx*ny;
        d[3] = 2*hmax + 2;
        d[4] = 2;
        for (iz1 = z1; iz1 <= z2; iz1++)
        {
            long int ny_long = (long int) ny;
            offset = iz1*nx*ny;
            cmplft_((float*)membuf+offset,
                    (float*)membuf+offset+1, &ny_long, d);
        }
        z1 = 0;
        z2 = zu - nz;
    }

    /*      transform on h */

    z1 = zl;
    z2 = std::min(zu, nz-1);
    for (iz = 0; iz < nzblk; iz++)
    {
        y_1 = yl;
        y_2 = std::min(yu, ny-1);
        for (iy = 0; iy < nyblk; iy++)
        {
            d[0] = nx*ny*(z2 - z1 + 1);
            d[1] = 2;
            d[2] = nx*ny;
            d[3] = nx*(y_2 - y_1 + 1);
            d[4] = nx;
            offset = (z1*ny + y_1)*nx;
            jj = (long) nx/2;
            hermft((float*)membuf+offset, (float*)membuf+offset+1, &jj, d);
            y_1 = 0;
            y_2 = yu - ny;
        }
        z1 = 0;
        z2 = zu - nz;
    }


    /*
       wrtbuf = (float *)calloc((xmax - xmin + 1)*(ymax - ymin + 1),sizeof(float));

       z1 = zmin;
       while (z1 <= zmax) {
       z2 = z1%nz;
       if (z2 < 0)
          z2 += nz;
       fwrite(&z1,   sizeof(int),   1, mpfil);
       fwrite(&xmin, sizeof(int), 1, mpfil);
       fwrite(&xmax, sizeof(int), 1, mpfil);
       fwrite(&ymin, sizeof(int), 1, mpfil);
       fwrite(&ymax, sizeof(int), 1, mpfil);
       offset = z2*nx*ny;
       wplane ((float *)membuf+offset, nx, ny, wrtbuf,
          xmax-xmin+1, ymax-ymin+1, xmin, ymin, z1);
       z1++;
       }

       z1 = -1;
       fwrite(&z1,   sizeof(int),   1, mpfil);
       fwrite(&z1,   sizeof(int),   1, mpfil);
       fwrite(&z1,   sizeof(int),   1, mpfil);
       fwrite(&z1,   sizeof(int),   1, mpfil);
       fwrite(&z1,   sizeof(int),   1, mpfil);
     */
    Logger::log("Fourier Transform complete");
    return ((float*)membuf);
}

void get_unit_cell()
{
    double s, ca, cb, sa, sb, cg, sg, vf, astar, bstar, cstar;

    raddeg = 4.0*atan(1.0)/180.0;
    s = raddeg*alpha;
    ca = cos(s);
    sa = sin(s);
    s = raddeg*beta;
    cb = cos(s);
    sb = sin(s);
    s = raddeg*gama;
    cg = cos(s);
    sg = sin(s);

    vf = 1.0/sqrt(1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg);
    vol = (float)(a*b*c/vf);
    astar = sa*vf/a;
    bstar = sb*vf/b;
    cstar = sg*vf/c;

    g11 = (float)(astar*astar);
    g22 = (float)(bstar*bstar);
    g33 = (float)(cstar*cstar);
    g12 = (float)(2.0*astar*bstar*(ca*cb - cg)/(sa*sb));
    g13 = (float)(2.0*astar*cstar*(cg*ca - cb)/(sg*sa));
    g23 = (float)(2.0*bstar*cstar*(cb*cg - ca)/(sb*sg));
}

/* Cyclic permutation subroutine */

void cycle(int *x, int *y, int *z)
{
    int w;

    w = *x;
    *x = *y;
    *y = *z;
    *z = w;
}

/*
 * Inverse metric tensor elements -- for non-permuted indices!
 */

/* ----------------------------------------------------------------------- */
/*
        additions to read_sf for spline intrplnts begin here.
        this particular method will only work for orthorhombic
        cells because we scaled before the symmetry expansion.
 */

/*
   psi_(p, N, k). called by : read_sf.

   input:
        p = h, k, or l - the fourier index.
        N = nx, ny, or nz.
        k = order of the b-spline function
   output:
        psi, the scaling function.
 */

#define QUADRATIC 3
#define CUBIC 4

double psi_(int p, int N, int k)
{
    int i;
    double arg, psi;

#define USESIN
#ifdef USESIN
    arg = p*PI;
    arg /= (double) N;

    arg = (p != 0) ? sin(arg)/arg : 1.;

    psi = 1.;
    for (i = 0; i < k; ++i)
    {
        psi *= arg;
    }
#else
    arg = 2 * p * PI;
    arg /= (double)N;
    if (k == QUADRATIC)
    {
        psi = 0.75 + 0.25 * cos(arg);
    }
    else
    {
        psi = (2.0/3.0) + (1.0/3.0) * cos(arg);
    }
#endif // ifdef USESIN
    return psi;
}

void read_sf(float *x, int nx, int ny, int nz, float scale, float temp, CREFL *refl, int nrefl, int usepsi)
{
    int count, p, q, r, used;
    float dsq, fh, fk, fl, s, t;
    fcomplex f, /*ww,*/ fsym;
    int nsymop = 0;
    char buf[100];

    int i, k;         /* k = 3, quadratic splines */
    double psi;
    int /*n_i[3],*/ N[4];

    k = usepsi + 2;
    N[1] = 2*nx;
    N[2] = ny;
    N[3] = nz;

    /*      nextf is a logical function used to read structure factors.
     *      If nextf is true, the arguments have returned the indices and
     *      value of the next structure factor in the input list.
     *      If nextf is false, there are no more structure factors.
     *
     *      xpnd is a logical function used to apply symmetry relations
     *      to structure factors.  If xpnd is true, a symmetry operator
     *      has been applied and the indices and value of the related
     *      structure factor have been returned.  If xpnd is false, there
     *      are no more symmetry operators to apply for this hkl.
     */
    s = (float)scale;
    t = (float)(temp/4.0);
    used = count = 0;

    while (nextf(hin, &f, count, refl, nrefl))
    {
        count++;
        fh = (float)hin[0];
        fk = (float)hin[1];
        fl = (float)hin[2];
        dsq = fh*(fh*g11+fk*g12+fl*g13) + fk*(fk*g22+fl*g23) + g33*fl*fl;


        if ((dsq <= dsqmax) && (dsq >= dsqmin))
        {
            used++;
            f.re = (float)(s*exp(-t*dsq)*f.re);
            f.im = (float)(s*exp(-t*dsq)*f.im);
            if (usepsi)
            {
                psi = 1.;
                for (i = 0; i < 3; ++i)
                {
                    psi *= psi_(hin[i], N[i + 1], k);
                }
                f.re /= (float)psi;
                f.im /= (float)psi;
            }

            /*
             *                        the following is a complex exponent
             *                        it is different from the above which is correct
             *			  summer-1990
             *			cexp (&ww, -t*dsq);
             *			c_mul(f,ww,f);
             *			f.re *= s;
             *			f.im *= s;
             */
            nsymop = 0;
            while (xpnd(hin, hout, &f, &fsym, nsymop))
            {
                p = hout[0];
                if (p < 0)
                {
                    sprintf(buf, "Error in symmetry expansion\n");
                    Logger::log(buf);
                    sprintf(buf, "position %d hin= %d %d %d hout= %d %d %d\n",
                            nsymop+1, hin[0], hin[1], hin[2], hout[0], hout[1], hout[2]);
                    Logger::log(buf);
                    return;
                }
                q = hout[1];
                r = hout[2];
                if (p >= nx)
                {
                    if (nx != 1)
                    {
                        expansion_error(nsymop);
                    }
                    return;
                }
                else
                if (abs(q) >= ny)
                {
                    if (ny != 1)
                    {
                        expansion_error(nsymop);
                    }
                    return;
                }
                else
                if (abs(r) >= nz)
                {
                    if (nz != 1)
                    {
                        expansion_error(nsymop);
                    }
                    return;
                }
                else
                {
                    if (q < 0)
                    {
                        q += ny;
                    }
                    if (r < 0)
                    {
                        r += nz;
                    }
                    x[2*(r*nx*ny + q*nx + p)] = fsym.re;
                    x[2*(r*nx*ny + q*nx + p)+1] = fsym.im;
                    if (p <= 0)
                    {

                        /*      Fill in the 0,-k,-l reflections by conjugate symmetry */

                        q = -hout[1];
                        r = -hout[2];
                        if (q < 0)
                        {
                            q += ny;
                        }
                        if (r < 0)
                        {
                            r += nz;
                        }
                        x[2*(r*nx*ny + q*nx)] =  fsym.re;
                        x[2*(r*nx*ny + q*nx)+1] = -fsym.im;
                    }
                }
                nsymop++;
            }
        }
    }
    sprintf(buf, "  %d structure factors used out of %d in input.\n\n", used, count);
    Logger::log(buf);
}

void expansion_error(int nsymop)
{
    sprintf(buf, "Symmetry expansion error generated index out of bounds\n");
    Logger::log(buf);
    sprintf(buf, "Indices read: %d %d %d Indices calc: %d %d %d\n",
            hin[0], hin[1], hin[2], hout[0], hout[1], hout[2]);
    Logger::log(buf);
    sprintf(buf, "Position number %d\n", nsymop+1);
    Logger::log(buf);
    sprintf(buf, "Grid limits: %d %d %d\n", 2*nx, ny, nz);
    Logger::log(buf);
}

int xpnd(int h1[], int h2[], fcomplex *f1, fcomplex *f2, int n)
{

    /*      Symmetry operators are stored in a common block.  They include
     *      permutation operators to give proper sectioning.
     *
     *      This routine applies symmetry operation n to the structure
     *      factor given by (h1, f1) and returns the result in (h2, f2).
     *      If there is no nth symmetry operation the value of xpnd is
     *      .false. and (h2, f2) are not changed.
     *
     *      The symmetry operators are stored as a series of 3x4 matrices.
     *      The space group equivalent positions are specified in terms of
     *      these matrices r as
     *		x' = r(1,1)*x + r(1,2)*y + r(1,3)*z + r(1,4)/12
     *		y' = r(2,1)*x + r(2,2)*y + r(2,3)*z + r(2,4)/12
     *		z' = r(3,1)*x + r(3,2)*y + r(3,3)*z + r(3,4)/12
     *      because twelve is the lowest common denominator for all of
     *      the possible fractional unit cell translations.
     *
     *      The equivalent positions in reciprocal space are generated from
     *      r by the relationship
     *		F(h') = exp(2*pi*i*<h!s>)*F(h)
     *      where h' is the product of (r transpose) and h, and <h!s> is
     *      the dot product of h and the fourth column of r scaled by 1/12.
     *
     *      All equivalent positions are transposed by conjugate symmetry
     *      where necessary into the unique hemisphere defined by
     *      (h .ge. 0) without regard to the values of k and l.
     */

    int i, j, temp;
    float c4 = 0.0, s4 = 1.0, c6 = 0.5F, s6 = ((float)sqrt(3.0)/2.0F);

    if ((temp = ((n >= 0) && (n <= nsymmops))) != 0)
    {

        j = (h1[0]*rsym[0][3][n]+h1[1]*rsym[1][3][n]+h1[2]*rsym[2][3][n])%12;
        if (j < 0)
        {
            j += 12;
        }

        switch (j)
        {

        case 0:
            f2->re = f1->re;
            f2->im = f1->im;
            break;
        case 2:
            cmul(f1->re, f1->im, c6, -s6, f2->re, f2->im);
            break;
        case 3:
            cmul(f1->re, f1->im, c4, -s4, f2->re, f2->im);
            break;
        case 4:
            cmul(f1->re, f1->im, -c6, -s6, f2->re, f2->im);
            break;
        case 6:
            f2->re = -f1->re;
            f2->im = -f1->im;
            break;
        case 8:
            cmul(f1->re, f1->im, -c6, s6, f2->re, f2->im);
            break;
        case 9:
            cmul(f1->re, f1->im, c4, s4, f2->re, f2->im);
            break;
        case 10:
            cmul(f1->re, f1->im, c6, s6, f2->re, f2->im);
            break;
        default:
            printf("Invalid symmetry operator for position %d\n", n);
            printf("Input indices = %d %d %d\n", h1[0], h1[1], h1[2]);
            printf("Symmetry operator:\n");
            for (i = 0; i < 3; i++)
            {
                printf("%d %d %d %d\n", rsym[i][0][n], rsym[i][1][n], rsym[i][2][n], rsym[i][3][n]);
            }
            return (-1);
        }
    }

    for (i = 0; i < 3; i++)
    {
        h2[i] = 0;
        for (j = 0; j < 3; j++)
        {
            h2[i] += h1[j]*rsym[j][i][n];
        }
    }

    if (h2[0] < 0)
    {
        for (i = 0; i < 3; i++)
        {
            h2[i] = -h2[i];
        }
        f2->im = -f2->im;
    }
    return (temp);
}

void wplane(float *x, int nx, int ny, float *p, int mx, int my, int x0, int y_0, int level)
{
    int i, j;
    double r, rsq = 0.0;

    /*rmax = x[ny*x0+y_0];*/
    rmin = rmax = x[nx*y_0+x0];
    iy = y_0 + 1;
    for (j = 0; j < my; j++)
    {
        ix = x0 + 1;
        for (i = 0; i < mx; i++)
        {
            /*r = x[ny*ix+iy];*/
            r = x[nx*(iy-1)+(ix-1)];
            p[i + mx*j] = (float)r;
            rmin = (float)std::min(r, (double)rmin);
            rmax = (float)std::max(r, (double)rmax);
            rsq  += r*r;
            if (ix >= nx)
            {
                ix = 0;
            }
            ix++;
        }
        if (iy >= ny)
        {
            iy = 0;
        }
        iy++;
    }

    /*write (mpfil) p*/
    /*
       fwrite(p,sizeof(float),mx*my,mpfil);
     */

    r = sqrt(rsq/(mx*my));
    /*printf(" Level %d: Min = %f Max = %f rms = %f\n",level, rmin, rmax, r);*/
    printf(" Level %d: Min = %10.8e Max = %10.8e rms = %10.8e\n", level, rmin, rmax, r);
}

int nextf(int h[3], fcomplex *f, int i, CREFL *refl, int nrefl)
{
    static float phi, fobs;
    if (i >= nrefl)
    {
        return (false);
    }

    h[0] = refl[i].ind[0];
    h[1] = refl[i].ind[1];
    h[2] = refl[i].ind[2];
    switch (maptype)
    {
    case MIMapType::Fo:
    case MIMapType::DirectFFT:
        fobs = refl[i].fo;
        break;
    case MIMapType::Fc:
        fobs = refl[i].fc;
        break;
    case MIMapType::TwoFoFc:
        fobs = 2.0F*refl[i].fo - refl[i].fc;
        break;
    case MIMapType::FoFc:
        fobs = refl[i].fo- refl[i].fc;
        break;
    case MIMapType::FoFo:
        fobs = refl[i].fo* refl[i].fo;
        break;
    case MIMapType::Fofom:
        fobs = refl[i].fo*refl[i].fom;
        break;
    case MIMapType::ThreeFoTwoFc:
        fobs = 3.0F*refl[i].fo - 2.0F*refl[i].fc;
        break;
    case MIMapType::FiveFoThreeFc:
        fobs = 5.0F*refl[i].fo - 3.0F*refl[i].fc;
        break;
    case MIMapType::TwoFoFcSigmaA:
        fobs = refl[i].acalc;
        break;
    case MIMapType::FoFcSigmaA:
        fobs = refl[i].bcalc;
        break;
    default:
        fobs = refl[i].fo;
        break;
    }
    phi = refl[i].phi;
    phi *= (float)raddeg;
    f->re = fobs*(float)cos(phi);
    f->im = fobs*(float)sin(phi);
    refl[i].coef = (float)fabs(fobs);
    return (true);
}

void cexp(fcomplex *a, float b)
{
    a->re = (float)cos(b);
    a->im = (float)sin(b);
}

/* return true if str2 is in str1 */

int str_index(char *str1, char *str2)
{
    int i, len1, len2;

    len1 = strlen(str1);
    len2 = strlen(str2);
    for (i = 0; i < len1; i++)
    {
        if (!strncmp(str1+i, str2, len2))
        {
            return (1);
        }
    }
    return (0);
}

int my_index(char *str, char ch)
{
    int i = 0;

    while (*str)
    {
        if (*str++ == ch)
        {
            return (i+1);
        }
        i++;
    }
    return (0);
}

void put_80(char *cptr, FILE *file)
{
    char lineout[82];
    int len;

    strncpy(lineout, cptr, 80);
    len = strlen(lineout);
    printf(" header length is %d\n", len);
    while (((int)strlen(lineout)) < 80)
    {
        strncat(lineout, "               ", 80-strlen(lineout));
    }

    strncat(lineout, "\n", 1);
    fwrite(lineout, 81, 1, file);
}

/*
        below are additions to FFTApply to generate spline interpolants
        on the same grid points as a straigt fft would. the interpolated
    values of rho are returned by intrplnt_() in the array S[].
 */

static int SplineOrder;
static double HalfSplineOrder;

/* t_(j,k,n) computes the knot tj=[j-(k/2)]/n */

#define t_(j, N) (((double)j - HalfSplineOrder)/((double)N))

/* B_(x, j, k, N) computes the b-spline bj(x), at knot j, of orders k=4, and 3,    on the uniform knot sequence tj = [j - (k/2)]/N
 */

double B_(double x, int j, int N)
{
    double tj0, tj1, tj2, tj3, tj4;
    register double y;
    if (SplineOrder == QUADRATIC)
    {
        tj3 = t_(j + SplineOrder, N);
        if (x < tj3)
        {
            tj0 = t_(j, N);
            if (x > tj0)
            {
                tj1 = t_(j+1, N);
                if (x < tj1)
                {
                    y = N*(x - tj0);
                    return .5*y*y;
                }
                else
                {
                    tj2 = t_(j+2, N);
                    if (x < tj2)
                    {
                        y = N*(x - tj1);
                        return -y*y + y + .5;
                    }
                    else
                    {
                        y = N*(x - tj2);
                        return .5*y*y - y + .5;
                    }
                }
            }
        }
    }
    else
    {
        tj4 = t_(j + SplineOrder,  N);
        if (x < tj4)
        {
            tj0 = t_(j, N);
            if (x > tj0)
            {
                tj2 = t_(j+2, N);
                if (x < tj2)
                {
                    tj1 = t_(j+1, N);
                    if (x < tj1)
                    {
                        y = N*(x - tj0);
                        return y*y*y/6.;
                    }
                    else
                    {
                        y = N*(x - tj1);
                        return (1./6.)+(-y*y*y+y*y+y)/2.;
                    }
                }
                else
                {
                    tj3 = t_(j+3, N);
                    if (x < tj3)
                    {
                        y = N*(x - tj2);
                        return .5*y*y*y - y*y + (2./3.);
                    }
                    else
                    {
                        y = N*(x - tj3);
                        return (-y*y*y + 3.*y*y - 3.*y + 1)/6.;
                    }
                }
            }
        }
    }
    return (double)0.0;
}

#ifdef FUTURE
float
avgsplinerho(map, fx, fy, fz)
int map;
float fx, fy, fz;
{
    extern CMapHeaderBase MapHeader[];
    extern int *MapRho[];
    extern int SplineMap[];
    double bx[4], by[4];
    double x, y, z, fxf, fyf, fzf, r, rx, rxy;
    int px[4], py[4];
    int i, ix, iy, jx, jy, jz, jx0, jy0, jz0, jz1;
    int x0, y0, z0, nx, ny, nz, nxny;
    float *A, *az, *azy;

    SplineOrder = SplineMap[map] + 2;
    HalfSplineOrder = SplineOrder / 2.0;
    nx = MapHeader[map].nx;
    ny = MapHeader[map].ny;
    nz = MapHeader[map].nz;
    if (fx < 0.0)
    {
        fx += ((int)(-fx)) + 1.0;
    }
    if (fy < 0.0)
    {
        fy += ((int)(-fy)) + 1.0;
    }
    if (fz < 0.0)
    {
        fz += ((int)(-fz)) + 1.0;
    }
    fx *= (double)nx;
    fy *= (double)ny;
    fz *= (double)nz;
    x0 = (int)fx;
    y0 = (int)fy;
    z0 = (int)fz;
    fxf = fx - x0;
    fyf = fy - y0;
    fzf = fz - z0;
    /* make sure indices are in bounds */
    x0 = (x0)%nx;
    y0 = (y0)%ny;
    z0 = (z0)%nz;
    fx = x0 + fxf;
    fy = y0 + fyf;
    fz = z0 + fzf;
    x = fx/nx;
    y = fy/ny;
    z = fz/nz;

    A = (float*)(MapRho[map]);
    jx0 = 1 + (int)(HalfSplineOrder + fx) - SplineOrder;
    jy0 = 1 + (int)(HalfSplineOrder + fy) - SplineOrder;
    jz1 = 1 + (int)(HalfSplineOrder + fz);
    jz0 = jz1 - SplineOrder;
    nxny = nx * ny;
    for (i = 0; i < SplineOrder; i++)
    {
        jx = jx0 + i;
        bx[i] = B_(x, jx, nx);
        px[i] = (jx + nx) % nx;
        jy = jy0 + i;
        by[i] = B_(y, jy, ny);
        py[i] = ((jy + ny) % ny) * nx;
    }
    r = 0;
    for (jz = jz0; jz < jz1; jz++)
    {
        az = &(A[((jz + nz) % nz) * nxny]);
        rxy = 0;
        for (iy = 0; iy < SplineOrder; iy++)
        {
            azy = &(az[py[iy]]);
            rx = 0;
            for (ix = 0; ix < SplineOrder; ix++)
            {
                rx += azy[px[ix]] * bx[ix];
            }
            rxy += rx * by[iy];
        }
        r += rxy * B_(z, jz, nz);
    }
    r *= MapHeader[map].scale;
    return (float)r;
}

float
avgsplinerho2(map, fx, fy, fz)
int map;
float fx, fy, fz;
{
    extern CMapHeaderBase MapHeader[];
    extern int *MapRho[];
    static double bx[SplineOrder], by[SplineOrder], bz[SplineOrder];
    static int px[SplineOrder], py[SplineOrder], pz[SplineOrder];
    static float prevfx, prevfy, prevfz;
    static int prevmap = -2, nx, ny, nz, nxny;
    double x, y, z, fxf, fyf, fzf, r, rx, rxy;
    int i, ix, iy, jx, jy, jz, jx0, jy0, jz0, jz1;
    int x0, y0, z0;
    float *A, *az, *azy;

    if (map != prevmap)
    {
        prevmap = map;
        prevfz = prevfy = prevfx = MAXFLOAT;
        if (prevmap == -1)
        {
            return 0.0;
        }
        nx = MapHeader[map].nx;
        ny = MapHeader[map].ny;
        nz = MapHeader[map].nz;
        nxny = nx * ny;
    }
    A = (float*)(MapRho[map]);
    if (fz != prevfz)
    {
        prevfz = fz;
        if (fz < 0.0)
        {
            fz += ((int)(-fz)) + 1.0;
        }
        fz *= (double)nz;
        z0 = (int)fz;
        fzf = fz - z0;
        z0 = (z0)%nz;
        fz = z0 + fzf;
        z = fz/nz;
        jz0 = 1 + (int)(HalfSplineOrder + fz) - SplineOrder;
        for (i = 0; i < SplineOrder; i++)
        {
            jz = jz0 + i;
            bz[i] = B_(z, jz, nz);
            pz[i] = ((jz + nz) % nz) * nxny;
        }
    }
    if (fy != prevfy)
    {
        prevfy = fy;
        if (fy < 0.0)
        {
            fy += ((int)(-fy)) + 1.0;
        }
        fy *= (double)ny;
        y0 = (int)fy;
        fyf = fy - y0;
        y0 = (y0)%ny;
        fy = y0 + fyf;
        y = fy/ny;
        jy0 = 1 + (int)(HalfSplineOrder + fy) - SplineOrder;
        for (i = 0; i < SplineOrder; i++)
        {
            jy = jy0 + i;
            by[i] = B_(y, jy, ny);
            py[i] = ((jy + ny) % ny) * nx;
        }
    }
    if (fx != prevfx)
    {
        prevfx = fx;
        if (fx < 0.0)
        {
            fx += ((int)(-fx)) + 1.0;
        }
        fx *= (double)nx;
        x0 = (int)fx;
        fxf = fx - x0;
        x0 = (x0)%nx;
        fx = x0 + fxf;
        x = fx/nx;
        jx0 = 1 + (int)(HalfSplineOrder + fx) - SplineOrder;
        for (i = 0; i < SplineOrder; i++)
        {
            jx = jx0 + i;
            bx[i] = B_(x, jx, nx);
            px[i] = (jx + nx) % nx;
        }
    }
    r = 0;
    for (iz = 0; iz < SplineOrder; iz++)
    {
        az = &(A[pz[iz]]);
        rxy = 0;
        for (iy = 0; iy < SplineOrder; iy++)
        {
            azy = &(az[py[iy]]);
            rx = 0;
            for (ix = 0; ix < SplineOrder; ix++)
            {
                rx += azy[px[ix]] * bx[ix];
            }
            rxy += rx * by[iy];
        }
        r += rxy * bz[jz];
    }
    r *= MapHeader[map].scale;
    return (float)r;
}
#endif // ifdef FUTURE

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
#define MAXBINS 50
int scompare(const void *in, const void *jn)
{
    CREFL *i = (CREFL*)in;
    CREFL *j = (CREFL*)jn;
    /* use with qsort to sort reflections according to sthol */
    int k = 1;
    if (i->sthol == j->sthol)
    {
        return (0);
    }
    if (i->sthol < j->sthol)
    {
        return (-1);
    }
    return (k);
}

int SigmaA_C(CREFL *refl, int nrefl, CMapHeaderBase *mhin)
{
    float *epsilon;
    //extern int scompare();
    //char buf[200];
    unsigned char *centric;
    int i, j, k, ibin, nbin[MAXBINS];
    long int ih, ik, il;
    /* increased to 3 x 3 x 96 - dem 7/31/97 */
    long int iss[864];
    /* increased to 3 x 96 - dem 7/31/97 */
    long int its[288];
    float eps, wt;
    long int mult, mk, iflg, iflg2, nsym;
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
    if (ratio <= 0.0)
    {
        ratio = 1.0;
    }
    if (ratio > 3.0)
    {
        ratio = 1.0;
    }
#ifdef DEBUG
    sprintf(buf, "Ratio of reference to working set is %f\n", ratio);
    Logger::log(buf);
#endif

    if (nrefl < 100)
    {
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
    if (!epsilon || !centric)
    {
        Logger::log("Out of memory in SigmaA: sorry");
        return 0;
    }
    /* first get the centric flag and the epsilon */
    /* set up the symmetry information */
    for (k = 0; k < mhin->nsym; k++)
    {
        for (j = 0; j < 3; j++)
        {
            its[j+k*3] =  ROUND(mhin->symops[j][3][k]*24.0);
            for (i = 0; i < 3; i++)
            {
                iss[i+3*(j+3*k)] = (int)mhin->symops[i][j][k];
            }
        }
    }
    for (i = 0; i < nrefl; i++)
    {
        ih = refl[i].ind[0];
        ik = refl[i].ind[1];
        il = refl[i].ind[2];
        nsym = mhin->nsym;
        stdrefl_(&ih, &ik, &il, &mult, &eps, &mk, &iflg,
                 &iflg2, iss, its, &nsym);
        epsilon[i] = eps;
        centric[i] = (unsigned char)(mk-1);
        /* check to make sure that Fc column
         * does not contain figure-of-merits
         */
        if (refl[i].fc <= 1.0)
        {
            nfom++;
        }
    }
    if ((float)nfom/(float)nrefl > 0.80)
    {
        Logger::log("SigmaA: Fc's are figure-of-merits.");
        printf("SigmaA: Fc's are figure-of-merits.\n");
        for (i = 0; i < nrefl; i++)
        {
            refl[i].fom = refl[i].fc;
        }
        return 0;
    }

    /* zero the sums */
    for (i = 0; i < numbins; i++)
    {
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
    for (i = 0; i < nrefl; i++)
    {
        ibin = (i*numbins)/(nrefl);
        wt = 2.0;
        if (centric[i])
        {
            wt = 1.0;
        }
        sumwt[ibin] += wt;
        nbin[ibin]++;
        sigman[ibin] += refl[i].fo*refl[i].fo*wt/epsilon[i];
        sigmap[ibin] += refl[i].fc*refl[i].fc*wt/epsilon[i];
    }

    for (ibin = 0; ibin < numbins; ibin++)
    {
        if (sumwt[ibin])
        {
            sigman[ibin] /= sumwt[ibin];
            sigmap[ibin] /= sumwt[ibin];
        }
    }

    for (i = 0; i < nrefl; i++)
    {
        ibin = (i*numbins)/(nrefl);
        eo = refl[i].fo/(float)sqrt(sigman[ibin]*epsilon[i]);
        ec = refl[i].fc/(float)sqrt(sigmap[ibin]*epsilon[i]);
        eo *= eo;
        ec *= ec;
        sum22[ibin] += eo*ec;
        sum20[ibin] += eo;
        sum02[ibin] += ec;
        sum40[ibin] += eo*eo;
        sum04[ibin] += ec*ec;
    }

#ifdef DEBUG
    sprintf(buf, "SigmaA by resolution\n  Bin     Resolution   Number     SigmaA   Correlation\n");
    Logger::log(buf);
#endif
    for (ibin = 0; ibin < numbins; ibin++)
    {
        if (nbin[ibin] > 0)
        {
            ntot += nbin[ibin];
            sum22T += sum22[ibin];
            sum20T += sum20[ibin];
            sum02T += sum02[ibin];
            sum40T += sum40[ibin];
            sum04T += sum04[ibin];
            signum = nbin[ibin]*sum22[ibin] - sum20[ibin]*sum02[ibin];
            if (signum > 0.0)
            {
                sigden = (float)sqrt((nbin[ibin]*sum40[ibin]-sum20[ibin]*sum20[ibin])
                                     *(nbin[ibin]*sum04[ibin]-sum02[ibin]-sum02[ibin]));
                sigmaA[ibin] = (float)sqrt(signum/sigden)*ratio;
            }
            else
            {
#ifdef DEBUG
                sprintf(buf, "Warning: correlation between E**2's is non-positive for bin %d - sigmaA set to 0.05\n", ibin);
                Logger::log(buf);
#endif
                sigmaA[ibin] = 0.05F;
            }
            corr = (nbin[ibin]*sum22[ibin] - sum20[ibin]*sum02[ibin])
                   /(float)sqrt((nbin[ibin]*sum40[ibin]-sum20[ibin]*sum20[ibin])*(nbin[ibin]*sum04[ibin]-sum02[ibin]*sum02[ibin]));
#ifdef DEBUG
            sprintf(buf, "%5d  %6.2f-%6.2f   %6d    %6.5f   %6.5f\n", ibin+1,
                    (float)(0.5/refl[ibin*nrefl/numbins].sthol),
                    (float)(0.5/refl[std::min((ibin+1)*nrefl/numbins-1, nrefl-1)].sthol),
                    nbin[ibin], sigmaA[ibin], corr);
            Logger::log(buf);
#endif
        }
    }

    corr = (ntot*sum22T - sum20T*sum02T)
           /(float)sqrt((ntot*sum40T-sum20T*sum20T)*(ntot*sum04T-sum02T*sum02T));
#ifdef DEBUG
    sprintf(buf, "SigmaA: Overall correlation is %f for %d refl\n", corr, ntot);
    Logger::log(buf);
    Logger::log(buf);
#endif

    /* modify Fo and Fc according to sigmaA */
    for (i = 0; i < nrefl; i++)
    {
        ibin = i*numbins/nrefl;
        eo = refl[i].fo/(float)sqrt(sigman[ibin]*epsilon[i]);
        ec = refl[i].fc/(float)sqrt(sigman[ibin]*epsilon[i]);
        Xarg = sigmaA[ibin]*2.0F*eo*ec/(1.0F-sigmaA[ibin]*sigmaA[ibin]);
        if (!centric[i])
        {
            m = sim(Xarg);
            refl[i].acalc = (2.0F*m*eo-sigmaA[ibin]*ec)*(float)sqrt(sigman[ibin]*epsilon[i]);
            refl[i].bcalc = m*refl[i].fo-refl[i].fc;
        }
        else
        {
            m = (float)tanh(Xarg);
            refl[i].acalc = m*refl[i].fo;
            refl[i].bcalc = refl[i].fo-refl[i].fc;
        }
        refl[i].fom = m;
        mtot += m;
        /*
           if(i%100==0)printf("%3d %3d %3d %5.2f %5.2f %d %f\n",refl[i].ind[0],refl[i].ind[1],refl[i].ind[2],refl[i].fo, refl[i].fc, (int)centric[i], epsilon[i]);
         */
    }
#ifdef DEBUG
    sprintf(buf, "SigmaA: Average figure of merit is %0.3f\n", mtot/(float)ntot);
    Logger::log(buf);
#endif

    free(epsilon);
    free(centric);
    return nrefl;
}

float sim(float x)
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
    if (fabs(x) > 0.0001)
    {
        t = (float)(fabs(x)/3.75);
        if (t > 1.0F)
        {
            t = 1.0F/t;
            sim = ((((((((0.01787654F-t*0.00420059F)*t+ (-0.02895312F))*t
                        +0.02282967F)*t+ (-0.01031555F))*t+0.00163801F)*t
                     +(-0.00362018F))*t+ (-0.03988024F))*t+0.39894228F)
                  *SIGN(x)/ ((((((((-0.01647633F+t*0.00392377F)*t
                                   +0.02635537F)*t+ (-0.02057706F))*t+0.00916281F)*t
                                +(-0.00157565F))*t+0.00225319F)*t+0.01328592F)*t+0.39894228F);
        }
        else
        {
            t = t*t;
            sim = ((((((t*0.00032411F + 0.00301532F)*t+0.02658733F)*t
                      +0.15084934F)*t+0.51498869F)*t+0.87890594F)*t+0.5F)*x
                  /((((((t*0.0045813F+0.0360768F)*t+0.2659732F)*t
                       +1.2067492F)*t+3.0899424F)*t+3.5156229F)*t+1.0F);
        }
    }
    return sim;
}

/* Subroutine */
int matset(float *am, float *dv, long int *nm, long int *nv)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i, j, jk;

    /* Parameter adjustments */
    --dv;
    --am;

    /* Function Body */
    jk = *nm;
    i__1 = *nv;
    for (i = 1; i <= i__1; ++i)
    {
        i__2 = *nv;
        for (j = i; j <= i__2; ++j)
        {
            am[jk] += dv[i] * dv[j];
            /* L10: */
            --jk;
        }
        /* L20: */
    }
    return 0;
} /* matset_ */

/* Subroutine */
int solv(float *am, float *v, float *dv, float *diag, long int *nm, long int *nv, float *sig)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int i, j, ij, ijd;
    static float pdi;

    /* Parameter adjustments */
    --diag;
    --dv;
    --v;
    --am;

    /* Function Body */
    i__1 = *nv;
    for (i = 1; i <= i__1; ++i)
    {
        pdi = (float)0.;
        ij = *nm - i + 1;
        ijd = *nv - 1;
        i__2 = *nv;
        for (j = 1; j <= i__2; ++j)
        {
            pdi += am[ij] * v[j];
            if ((i__3 = j - i) < 0)
            {
                goto L4800;
            }
            else if (i__3 == 0)
            {
                goto L4850;
            }
            else
            {
                goto L4900;
            }
L4800:
            ij -= ijd;
            --ijd;
            goto L4950;
L4850:
            diag[i] = am[ij];
L4900:
            --ij;
L4950:
            ;
        }
        dv[i] = pdi;
        *sig -= pdi * v[i];
        /* L5000: */
    }
    return 0;
} /* solv_ */

int smi(float *am, long int*, long int *n, long int *nfail)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static int imax;
    static float suma, term;
    static int i, j, k, l, m;
    static float denom;
    static int ii, kdm, kli, kmi, iwn;

    /* Parameter adjustments */
    --am;

    /* Function Body */
    iwn = *n * (*n + 1) / 2 + 1;
    k = iwn - 1;
    if ((i__1 = *n - 1) < 0)
    {
        goto L1400;
    }
    else if (i__1 == 0)
    {
        goto L1050;
    }
    else
    {
        goto L1100;
    }
L1050:
    am[1] = (float)1. / am[1];
    goto L1900;
L1100:
    i__1 = *n;
    for (m = 1; m <= i__1; ++m)
    {
        imax = m - 1;
        i__2 = *n;
        for (l = m; l <= i__2; ++l)
        {
            suma = (float)0.;
            kli = iwn - l;
            kmi = iwn - m;
            if (imax <= 0)
            {
                goto L1250;
            }
            else
            {
                goto L1150;
            }
L1150:
            i__3 = imax;
            for (i = 1; i <= i__3; ++i)
            {
                suma += am[kli] * am[kmi];
                j = *n - i;
                kli -= j;
                /* L1200: */
                kmi -= j;
            }
L1250:
            term = am[k] - suma;
            if (l - m <= 0)
            {
                goto L1300;
            }
            else
            {
                goto L1450;
            }
L1300:
            if (term <= (float)0.)
            {
                goto L1400;
            }
            else
            {
                goto L1350;
            }
L1350:
            denom = (float)sqrt(term);
            am[k] = denom;
            goto L1500;
L1400:
            *nfail = iwn - k;
            goto L1950;
L1450:
            am[k] = term / denom;
L1500:
            --k;
        }
        /* L1550: */
    }
    am[iwn - 1] = (float)1. / am[iwn - 1];
    kdm = iwn - 1;
    i__1 = *n;
    for (l = 2; l <= i__1; ++l)
    {
        kdm = kdm - *n + l - 2;
        term = (float)1. / am[kdm];
        am[kdm] = term;
        kmi = iwn;
        kli = iwn - l;
        imax = l - 1;
        i__2 = imax;
        for (m = 1; m <= i__2; ++m)
        {
            k = kli;
            suma = (float)0.;
            i__3 = imax;
            for (i = m; i <= i__3; ++i)
            {
                ii = kmi - i;
                suma -= am[kli] * am[ii];
                /* L1600: */
                kli = kli - *n + i;
            }
            am[k] = suma * term;
            j = *n - m;
            kli = k - j;
            /* L1650: */
            kmi -= j;
        }
        /* L1700: */
    }
    k = iwn - 1;
    i__1 = *n;
    for (m = 1; m <= i__1; ++m)
    {
        kli = k;
        i__2 = *n;
        for (l = m; l <= i__2; ++l)
        {
            kmi = k;
            imax = *n - l + 1;
            suma = (float)0.;
            i__3 = imax;
            for (i = 1; i <= i__3; ++i)
            {
                suma += am[kli] * am[kmi];
                --kli;
                /* L1750: */
                --kmi;
            }
            am[k] = suma;
            /* L1800: */
            --k;
        }
        /* L1850: */
    }
L1900:
    *nfail = 0;
L1950:
    return 0;
} /* smi_ */

